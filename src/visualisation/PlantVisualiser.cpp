#include "PlantVisualiser.h"

#include "MappedOrganism.h"
#include "Quaternion.h"
#include "CatmullRomSpline.h"
#include "mymath.h"

#include <algorithm>

namespace CPlantBox {

  	// an iterator that can give me the vector mirrored around 0
	class MirrorIterator : public std::iterator<std::input_iterator_tag, std::pair<int,double> > {
	public:
		MirrorIterator(const std::vector<double>* v) : v(v) {
			//std::cout << "MirrorIterator was created " << v->size() << std::endl;
			// output all vector elements
      // std::copy(v->begin(), v->end(), std::ostream_iterator<double>(std::cout, " ")); std::cout << std::endl;
		}
		std::pair<int,double> operator*() { return std::make_pair(idx(), v->at(i)); }
		MirrorIterator& operator++() { inc(); return *this; }
		bool operator!=(const MirrorIterator& other) { return i != other.i && r != other.r; }
    bool operator==(const MirrorIterator& other) { return i == other.i && r == other.r; }
		std::size_t size() { return (v->size() > 0) ? (v->size() * 2 - (v->back() < std::numeric_limits<float>::epsilon() ? 1 : 0)) : 0; }
		double operator[](int i) {
			if(i < v->size())
				return v->at(i);
			else
			{
				return -v->at(v->size() - (i - v->size() + 1) - ((v->back() < std::numeric_limits<float>::epsilon()) ? 1 : 0));
			}
		 }
		// begin
		MirrorIterator begin() { return MirrorIterator(v); }
		// end
		MirrorIterator end() { return MirrorIterator(v, true, 0); }

    // computes the texture coordinate within the unit interval based on the index
    double texcoord(int i)
		{
		  double d = static_cast<double>(i) / static_cast<double>(size() - 1);
			return d;
		}
		void inc()
		{
			if(r)
			{
				--i;
			}
			else
			{
				++i;
				if(i == v->size())
				{
					r = true;
					i = v->size() - (v->back() < std::numeric_limits<float>::epsilon() ? 2 : 1);
				}
			}
		}
		// checks whether we are in the mirrored part for a given index
		bool isMirrored(int i)
		{
		  return i >= v->size();
    }
		int idx()
		{
			if(r)
				return v->size() + i - (v->back() < std::numeric_limits<float>::epsilon() ? 1 : 0);
			else
				return i;
		}
	private:
		MirrorIterator(const std::vector<double>* v, bool end, int i) : v(v), r(end), i(i) { }
		bool r = false;
		const std::vector<double>* v;
		int i = 0;
	};

  /** 
   * Hidden Template method to write a set of vectors into a buffer
   * @param buffer    the buffer
   * @param offset    the offset in the buffer
   * @param v...     the vectors
   */
  inline unsigned int vec2Buf(std::vector<double>& buffer, unsigned int offset, const Vector3d& v) {
    buffer[offset + 0] = v.x;
    buffer[offset + 1] = v.y;
    buffer[offset + 2] = v.z;
    return offset + 3;
  }
  template<typename... Args>
  inline unsigned int vec2Buf(std::vector<double>& buffer, unsigned int offset, const Vector3d& v, Args... args) {
    offset = vec2Buf(buffer, offset, v);
    return vec2Buf(buffer, offset, args...);
  }

  void ClampVectorBetweenLengths(Vector3d& v, double min, double max) {
    double l = v.length();
    if(l < 1e5)
    {
        return;
    }
    else if(l < min)
    {
      v.normalize();
      v = v.times(min);
    }
    else if(l > max)
    {
      v.normalize();
      v = v.times(max);
    }
  }

PlantVisualiser::PlantVisualiser() :
    geometry_resolution_(10), leaf_resolution_(10), include_midline_in_leaf_(false) {}

PlantVisualiser::PlantVisualiser(const PlantVisualiser& pv) = default;
PlantVisualiser::PlantVisualiser(std::shared_ptr<MappedPlant> plant) :
    geometry_resolution_(10), leaf_resolution_(10), include_midline_in_leaf_(false), plant_(plant) 
    {
      BuildAttachmentMap();
    }
PlantVisualiser::~PlantVisualiser() = default;

void PlantVisualiser::setPlant(const std::shared_ptr<MappedPlant>& plant) {
  plant_ = plant;
  BuildAttachmentMap();
}

void PlantVisualiser::ComputeGeometryForOrgan(int organId)
{
  auto result = plant_->getOrgans(-1, false);
  auto organ_it = std::find_if(result.begin(), result.end(), [organId](const auto& o) {
    return o->getId() == organId;
  });
  if(organ_it == result.end())
  {
    std::stringstream errMsg;
    errMsg << "MappedPlant::ComputeGeometryForOrgan: organ not found: " << organId;
    throw std::runtime_error(errMsg.str().c_str());
  }
  if(verbose_)
  {
    std::cout << "Computing geometry for organ: " << organId << " of type " << (*organ_it)->organType() << std::endl;
    std::cout << "Organ has " << (*organ_it)->getNumberOfNodes() << " nodes." << std::endl;
    std::cout << "Leaf resolution is " << leaf_resolution_ << std::endl;
  }
  auto organ = *organ_it;
  
  if(organ->organType() == 4)
  {
    // 4 SHOULD mean leaf, so we do not check for successful cast
		unsigned int point_space = 0, cell_space = 0;
		
    if(this->include_midline_in_leaf_)
    {
      GenerateStemGeometry(organ, point_space, cell_space);
      point_space += organ->getNumberOfNodes() * 3 * geometry_resolution_;
      point_space = geometry_.size();
      cell_space += (organ->getNumberOfNodes() - 1) * 6 * geometry_resolution_;
      cell_space = geometry_indices_.size();
    }
		auto leaf = std::dynamic_pointer_cast<Leaf>(organ);
		if(leaf->getLeafRandomParameter()->parametrisationType == 1)
		{
			GenerateRadialLeafGeometry(leaf, point_space, cell_space);
			point_space += leaf->getNumberOfNodes() * 6 * 3;
			cell_space += leaf->getNumberOfNodes() * 6;
		}
		else
		{
			int petiole_zone = 0;
			for(int i = 0; i < leaf->getNumberOfNodes(); i++)
			{
				if(leaf->nodeLeafVis(leaf->getLength(i)))
				{
					petiole_zone = i;
					break;
				}
			}
			if(petiole_zone + 1 < leaf->getNumberOfNodes())
			{
				//GenerateLeafGeometry(leaf, petiole_zone, point_space, cell_space);
        auto lrp = leaf->getLeafRandomParameter();
        lrp->createLeafGeometry(lrp->leafGeometryPhi,lrp->leafGeometryX,leaf_resolution_);

				GenerateRadialLeafGeometry(leaf, point_space, cell_space);
				//point_space += (organ->getNumberOfNodes() - petiole_zone) * 6 * 3;
        //point_space = geometry_.size();
				//cell_space += (organ->getNumberOfNodes() - petiole_zone) * 6;
        //cell_space = geometry_indices_.size();
			point_space += leaf->getNumberOfNodes() * 6 * 3;
			cell_space += leaf->getNumberOfNodes() * 6;
			}
		}
  }
  else
  {
    GenerateStemGeometry(organ);
  }
}

void PlantVisualiser::ComputeGeometryForOrganType(int organType, bool clearFirst)
{
  auto organ_list = plant_->getOrgans(-1, false);
		
  // First we check if we have enough memory to support the geometry_
  unsigned int point_space = 0;
  unsigned int cell_space = 0;
	unsigned int num_organs = 0;

  for(auto organ : organ_list)
  {
    // Check against object, because organType can be -1
		if(organ->organType() == organType || organType < 0)
		{
			if(organ->organType() == 4)
			{
				// 4 SHOULD mean leaf, so we do not check for successful cast
				if(include_midline_in_leaf_)
				{
				  point_space += organ->getNumberOfNodes() * 3 * geometry_resolution_;
				  cell_space += (organ->getNumberOfNodes() - 1) * 6 * geometry_resolution_;
				}
			  point_space += (organ->getNumberOfNodes()) * 4 * 3;
				cell_space += ((organ->getNumberOfNodes()) - 1) * 4 * 3 + 40;
			}
			else
			{
				point_space += organ->getNumberOfNodes() * 3 * geometry_resolution_;
				cell_space += (organ->getNumberOfNodes() - 1) * 6 * geometry_resolution_;
			}
		}
  }
  //std::cout << "Going to allocate " << point_space << " points and " << cell_space << " cells" << std::endl;
  if(clearFirst)
  {
    geometry_.clear();
    geometry_normals_.clear();
    geometry_colors_.clear();
    geometry_indices_.clear();
    geometry_texture_coordinates_.clear();
    geometry_node_ids_.clear();
  }
  geometry_.reserve(geometry_.size() + point_space);
  geometry_normals_.reserve(geometry_normals_.size() + point_space);
  geometry_colors_.reserve(geometry_colors_.size() + (point_space / 3 * 2));
  geometry_indices_.reserve(geometry_indices_.size() + cell_space);
  geometry_texture_coordinates_.reserve(geometry_texture_coordinates_.size() + point_space);
  geometry_node_ids_.reserve(geometry_node_ids_.size() + (point_space / 3 + 1));


  point_space = 0;
  cell_space = 0;
  unsigned int checked_organs = 0;
  for(auto organ : organ_list)
  {
    checked_organs++;
		//std::cout << "Going through organ " << organ->getId() << std::endl;

    if((organType >= 1 && organ->organType() != organType) || organ->getNumberOfNodes() <= 1)
    {
      continue;
    }
    if(organ->organType() == 4)
    {
			auto leaf = std::dynamic_pointer_cast<Leaf>(organ);
			// render petiole
      //std::cout << "Generating geometry_ for leaf " << organ->getId() << " with " << organ->getNumberOfNodes() << " nodes." << std::endl;
      //std::cout << "Stem part for petiole and rest" << std::endl;
			if(include_midline_in_leaf_)
			{
        GenerateStemGeometry(organ, point_space, cell_space);
			}
      //std::cout << "Updating buffer positions because the leaf is a two-parter" << std::endl;
      //point_space += organ->getNumberOfNodes() * 3 * geometry_resolution_;
      //cell_space += (organ->getNumberOfNodes() - 1) * 6 * geometry_resolution_;
      point_space = geometry_.size();
      cell_space = geometry_indices_.size();
			
			if(leaf->getLeafRandomParameter()->parametrisationType == 1)
			{
        if(verbose_)
				  std::cout << "Generating radial leaf geometry_" << std::endl;
				GenerateRadialLeafGeometry(leaf, point_space, cell_space);
				point_space = geometry_.size();
				cell_space = geometry_indices_.size();
			}
			else
			{
				int petiole_zone = 0;
				for(int i = 0; i < leaf->getNumberOfNodes(); i++)
				{
					if(!leaf->nodeLeafVis(leaf->getLength(i)))
					{
						petiole_zone = i;
					}
					else break;
				}
				if(petiole_zone + 1 < leaf->getNumberOfNodes())
				{
	        auto lrp = leaf->getLeafRandomParameter();
          lrp->createLeafGeometry(lrp->leafGeometryPhi,lrp->leafGeometryX,leaf_resolution_);

					GenerateRadialLeafGeometry(leaf, point_space, cell_space);
					point_space = geometry_.size();
					cell_space = geometry_indices_.size();
				}
			}
    }
    else
    {
			//std::cout << "Organ is a stem" << std::endl;
      //std::cout << "Generating geometry_ for stem " << organ->getId() << " with " << organ->getNumberOfNodes() << " nodes." << std::endl;
      auto prev_size = geometry_indices_.size();
      GenerateStemGeometry(organ, point_space, cell_space);
      //std::cout << "Organ " << organ->getId() << " pushed the size from " << prev_size << " to " << geometry_indices_.size() << std::endl;
      //point_space += organ->getNumberOfNodes() * 3 * geometry_resolution_;
      //cell_space += (organ->getNumberOfNodes() - 1) * 6 * geometry_resolution_;
      point_space = geometry_.size();
      cell_space = geometry_indices_.size();
    }
    //std::cout << "Done generating geometry_ for organ " << organ->getId() << std::endl;

  }
  assert(geometry_.size() % 3 == 0);
  assert(geometry_indices_.size() % 3 == 0);
}

void PlantVisualiser::ComputeGeometry()
{
  this->ComputeGeometryForOrganType(-1);
  if(verbose_)
  {
    std::cout << "Sanity Check for C++ " << std::endl;
    std::cout << "Geometry size: " << geometry_.size() << std::endl;
    std::cout << "Geometry Indices size: " << geometry_indices_.size() << std::endl;
    std::cout << "Geometry Normals size: " << geometry_normals_.size() << std::endl;
    std::cout << "Geometry Colors size: " << geometry_colors_.size() << std::endl;
    std::cout << "Geometry Texture Coordinates size: " << geometry_texture_coordinates_.size() << std::endl;
    std::cout << "Geometry Node Ids size: " << geometry_node_ids_.size() << std::endl;
    std::cout << "This would mean we have " << geometry_.size() / 3 << " points and " << geometry_indices_.size() / 3 << " cells." << std::endl;
  }
}

void PlantVisualiser::ResetGeometry()
{
  if(verbose_) std::cout << "Resetting the vis geometry" << std::endl;
  geometry_.clear();
  geometry_normals_.clear();
  geometry_colors_.clear();
  geometry_indices_.clear();
  geometry_texture_coordinates_.clear();
  geometry_node_ids_.clear();
}

void PlantVisualiser::MapPropertyToColors(std::vector<double> property, std::pair<double, double> minMax)
{
  // TODO: This is a very naive implementation. It should be improved or outsourced to Unreal.
  geometry_colors_.clear();
  geometry_colors_.resize(geometry_.size(), 0.0);
  for(int i = 0; i < geometry_colors_.size(); i++)
  {
    double value = property[i];
    double range = minMax.second - minMax.first;
    double normalized = (value - minMax.first) / range;
    // cast the double to the unsigned short
    geometry_colors_[i] = static_cast<unsigned short>(normalized * 65535.0);
  }
}

int PlantVisualiser::GetNumOrgans() const
{ return plant_->getOrgans(-1, false).size(); }

std::string PlantVisualiser::SelfCheck() const
{
  std::string result = "{";

  // check for NaNs in the geometry_
  result += "\"vertex_nan\":";
  bool has_nan = false;
  for(auto v : geometry_)
  {
    if(std::isnan(v))
    {
      has_nan = true;
      break;
    }
  }
  result += has_nan ? "true" : "false";

  // check for NaNs in the geometry_normals_
  result += ",\"normal_nan\":";
  has_nan = false;
  for(auto v : geometry_normals_)
  {
    if(std::isnan(v))
    {
      has_nan = true;
      break;
    }
  }
  result += has_nan ? "true" : "false";

  result += "}";
  return result;
}

void PlantVisualiser::BuildAttachmentMap()
{
  //std::cout << "Building attachment map" << std::endl;
  leaf_attachment_map_.clear();
  for(auto organ : plant_->getOrgans(Organism::OrganTypes::ot_leaf, true))
  {
    if(organ->getParent() != nullptr)
    {
      const Vector3d parent_position = organ->getParent()->getNode(organ->parentNI);
      const Vector3d organ_position = organ->getNode(0);
      leaf_attachment_map_[organ->getId()] = std::make_pair(organ->getParent()->getId(), vectorNormalized(organ_position - parent_position));
    }
  }
}

void PlantVisualiser::GenerateLeafGeometry(std::shared_ptr<Leaf> leaf, unsigned int petiole_zone, unsigned int p_o, unsigned int c_o)
{
  // std::vector::reserve should be idempotent.
  //std::cout << "Resizing geometry_ buffers for a leaf with n=" << leaf->getNumberOfNodes() << ", pet=" << petiole_zone << std::endl;
  if(verbose_)
  {
    std::cout << "Generating NORMAL geometry for leaf " << leaf->getId() << " with " << leaf->getNumberOfNodes() << " nodes and petiole zone " << petiole_zone << std::endl;
  }
	int first_surface_id = petiole_zone + 1;
	int total_points = leaf->getNumberOfNodes() - first_surface_id;
	geometry_.resize(std::max(static_cast<std::size_t>(p_o + total_points * 4 * 3), geometry_.size()), -1.0);
	geometry_normals_.resize(std::max(static_cast<std::size_t>(p_o + total_points * 4 * 3), geometry_normals_.size()), -1.0);
	geometry_indices_.resize(std::max(static_cast<std::size_t>(c_o + (total_points - 1) * 12), geometry_indices_.size()), static_cast<unsigned int>(-1));
	geometry_colors_.resize(std::max(static_cast<std::size_t>((p_o/3) + total_points * 4 * 3), geometry_colors_.size()), static_cast<unsigned short>(-1));
	geometry_texture_coordinates_.resize(std::max(static_cast<std::size_t>((p_o/3*2) + total_points * 4 * 2), geometry_texture_coordinates_.size()), -1.0);
	geometry_node_ids_.resize(std::max(static_cast<std::size_t>(p_o/3 + total_points * 4), geometry_node_ids_.size()), -1);

  if(verbose_) std::cout << "Orientation generation" << std::endl;
  unsigned int points_index = p_o;
  unsigned int cell_index = c_o;
  Quaternion rot = Quaternion::geodesicRotation({1,0,0},{0,0,1}).normalized();
  Vector3d lastPosition = leaf->getNode(first_surface_id);
  // TODO: check if getLength(true) is the correct length
  double totalLenght = leaf->getLength(true);
  double currentLength = leaf->getLength(first_surface_id);
  for(int i = first_surface_id; i < leaf->getNumberOfNodes(); ++i)
  {
    // we presume that the nodes follow the center line of the plant
    if(verbose_) std::cout << " going through node id " << i << " on this leaf." << std::endl;
    // and that the leaf is oriented along the x axis
    auto position = leaf->getNode(i);
    auto id = leaf->getNodeId(i);
    Vector3d dist;
    // This part is for the rotation of the leaf segment
    if(i + 1 < leaf->getNumberOfNodes())
    {
      dist = leaf->getNode(i + 1) - position;
    }
    else
    {
      dist = leaf->getNode(i - 1) - position;
    }
    // This, in contrast, is for the texture mapping
    currentLength += (position - lastPosition).length() / totalLenght;
    //rot *= Quaternion::geodesicRotation(rot.Forward(), dist);
    rot = Quaternion::FromForwardAndUp(rot.Forward(), rot.Up());
		if(verbose_) std::cout << "[Leaf] Rotating " << rot.toString() << " to get " << dist.toString() << std::endl;
    // TODO: Check with mona on what the Vector3d coordinates of
    // this function are, and if we need to change them
    auto vis = leaf->getLeafVis(i);
    // We don't normally split normals, but in this case we have a flat surface
    // and we want to have a smooth shading
    //std::cout << "Inserting some geometry_" << std::endl;
    if(verbose_) std::cout << "Vis has length " << vis.size() << std::endl;
    geometry_texture_coordinates_.insert(geometry_texture_coordinates_.begin() + (points_index/3*2),
                                            {currentLength, 0.0, currentLength, 1.0});
    //std::cout << "geometry_node_ids_[" << points_index << "] = " << id << std::endl;
    //std::cout << "Vector Data: Size=" << geometry_node_ids_.size() << ", capacity=" << geometry_node_ids_.capacity() << std::endl;
    //geometry_node_ids_[points_index] = id;
    //geometry_node_ids_[points_index + 1] = id;
    // TODO: it not obvious that points_index can be changed by the insert here
    geometry_node_ids_[points_index/3] = id;
    geometry_node_ids_[(points_index/3) + 1] = id;
    vec2Buf(geometry_, points_index, vis[0], vis[1]);
    points_index = vec2Buf(geometry_normals_, points_index, rot.Up(), rot.Up());
    vec2Buf(geometry_, points_index, vis[0], vis[1]);
    points_index = vec2Buf(geometry_normals_, points_index ,-rot.Up(), -rot.Up());
    // The triangles are defined clockwise for the front face and counter clockwise for the back face
		
		unsigned int point_index_offset = points_index / 3;
      //std::cout << "Inserting some indices: " << geometry_indices_.size() << " + 6 < " << geometry_indices_.capacity() << std::endl;
		if(i > first_surface_id)
		{
			geometry_indices_[cell_index +  0] = point_index_offset-7;
			geometry_indices_[cell_index +  1] = point_index_offset-8;
			geometry_indices_[cell_index +  2] = point_index_offset-4;
			geometry_indices_[cell_index +  3] = point_index_offset-7;
			geometry_indices_[cell_index +  4] = point_index_offset-4;
			geometry_indices_[cell_index +  5] = point_index_offset-3;
			geometry_indices_[cell_index +  6] = point_index_offset-5;
			geometry_indices_[cell_index +  7] = point_index_offset-6;
			geometry_indices_[cell_index +  8] = point_index_offset-2;
			geometry_indices_[cell_index +  9] = point_index_offset-5;
			geometry_indices_[cell_index + 10] = point_index_offset-2;
			geometry_indices_[cell_index + 11] = point_index_offset-1;
			cell_index += 12;
		}
		//std::cout << "Inserted" << std::endl;
  }
      //std::cout << "Done" << std::endl;
}

void PlantVisualiser::GenerateStemGeometry(std::shared_ptr<Organ> stem, unsigned int p_o, unsigned int c_o)
{
  //std::cout << "Generating Stem for " << stem->getId() << " and reserving buffers" << std::endl;
	geometry_.resize(std::max(static_cast<std::size_t>(p_o + (stem->getNumberOfNodes() * geometry_resolution_ * 3)), geometry_.size()),-1.0);
	geometry_normals_.resize(std::max(static_cast<std::size_t>(p_o + (stem->getNumberOfNodes() * geometry_resolution_ * 3)), geometry_normals_.size()),-1.0);
  //std::cout << "Wanting to resize " << geometry_indices_.size() << "/" << geometry_indices_.capacity() << " indices to " << static_cast<std::size_t>(c_o + (stem->getNumberOfNodes() -1) * geometry_resolution_ * 6) << std::endl;
  //std::cout << "c_o=" << c_o << ", stem->getNumberOfNodes()=" << stem->getNumberOfNodes() << ", geometry_resolution_=" << geometry_resolution_ << ", total is " << (stem->getNumberOfNodes() -1) * geometry_resolution_ * 6 << std::endl;
  geometry_indices_.resize(std::max(static_cast<std::size_t>(c_o + (stem->getNumberOfNodes() -1) * geometry_resolution_ * 6), geometry_indices_.size()),static_cast<unsigned int>(-1));
  geometry_colors_.resize(std::max(static_cast<std::size_t>(p_o + stem->getNumberOfNodes() * geometry_resolution_), geometry_colors_.size()),static_cast<unsigned short>(-1));
  geometry_texture_coordinates_.resize(std::max(static_cast<std::size_t>((p_o/3*2) + stem->getNumberOfNodes() * geometry_resolution_ * 2), geometry_texture_coordinates_.size()),-1.0);
  geometry_node_ids_.resize(std::max(static_cast<std::size_t>((p_o/3) + stem->getNumberOfNodes() * geometry_resolution_), geometry_node_ids_.size()),-1);
  Quaternion lastRotation = Quaternion::geodesicRotation({1,0,0},{0,0,1});

  for(auto i = 0; i < stem->getNumberOfNodes(); ++i)
  {
    double diameter = stem->getParameter("radius");
    const auto& node = stem->getNode(i);

		// if the current i is in the last 10% of the nodes
		if(static_cast<float>(i)/static_cast<float>(stem->getNumberOfNodes()) > 0.9f)
		{
		  // reduce the diameter to form a tip
			diameter *= (1.0 - ((static_cast<float>(i)/static_cast<float>(stem->getNumberOfNodes()) - 0.9f) * 10.0f));
		}
		if(i == stem->getNumberOfNodes() - 1)
		{
			diameter = 0.01;
		}

    Vector3d dist;
    if(i + 1 < stem->getNumberOfNodes())
    {
      dist = stem->getNode(i + 1) - node;
    }
    else
    {
      dist = node - stem->getNode(i - 1);
    }
    
		lastRotation = Quaternion::FromForwardAndUp(dist, lastRotation.Up());
    auto deltaphi = 2 * M_PI / geometry_resolution_;
		unsigned int point_index_offset = p_o / 3;
    for (auto j = 0; j < geometry_resolution_; ++j)
    {
      auto phi = j * deltaphi;
      Vector3d outer = {0.0, std::cos(phi) * diameter, std::sin(phi) * diameter};
      outer = lastRotation.Rotate(outer);
			//std::cout << "[Leaf] Rotating " << lastRotation.toString() << " to get " << outer.toString() << std::endl;
      // consequtive points are stored in the buffer
      // the plus 0 is not necessary, but it makes it easier to read
      geometry_[3 * (i * geometry_resolution_ + j) + 0 + p_o] = node.x + outer.x;
      geometry_[3 * (i * geometry_resolution_ + j) + 1 + p_o] = node.y + outer.y;
      geometry_[3 * (i * geometry_resolution_ + j) + 2 + p_o] = node.z + outer.z;
			// calculating the point index offset from the buffer offset
      // the indices are stored in the buffer, and are all front facing
			if (i > 0)
			{
				geometry_indices_[6 * ((i-1) * geometry_resolution_ + j) + 2 + c_o] = point_index_offset + ((i-1) + 1) * geometry_resolution_ + (j + 1) % geometry_resolution_;
				geometry_indices_[6 * ((i-1) * geometry_resolution_ + j) + 1 + c_o] = point_index_offset +  (i-1) * geometry_resolution_ + (j + 1) % geometry_resolution_;
				geometry_indices_[6 * ((i-1) * geometry_resolution_ + j) + 0 + c_o] = point_index_offset +  (i-1) * geometry_resolution_ + j;
				geometry_indices_[6 * ((i-1) * geometry_resolution_ + j) + 3 + c_o] = point_index_offset +  (i-1) * geometry_resolution_ + j;
				geometry_indices_[6 * ((i-1) * geometry_resolution_ + j) + 4 + c_o] = point_index_offset + ((i-1) + 1) * geometry_resolution_ + (j + 1) % geometry_resolution_;
				geometry_indices_[6 * ((i-1) * geometry_resolution_ + j) + 5 + c_o] = point_index_offset + ((i-1) + 1) * geometry_resolution_ + j;
			}
      // the normals are stored in the buffer
      geometry_normals_[3 * (i * geometry_resolution_ + j) + 0 + p_o] = outer.x;
      geometry_normals_[3 * (i * geometry_resolution_ + j) + 1 + p_o] = outer.y;
      geometry_normals_[3 * (i * geometry_resolution_ + j) + 2 + p_o] = outer.z;

      geometry_node_ids_[i * geometry_resolution_ + j + point_index_offset] = stem->getNodeId(i);

      geometry_texture_coordinates_[2 * (i * geometry_resolution_ + j) + 0 + point_index_offset*2] = i / (double)stem->getNumberOfNodes();
      geometry_texture_coordinates_[2 * (i * geometry_resolution_ + j) + 1 + point_index_offset*2] = phi / (2.0 * M_PI);
    }
  }
}

void PlantVisualiser::GenerateRadialLeafGeometry(std::shared_ptr<Leaf> leaf, unsigned int p_o, unsigned int c_o)
{
	// Fetch the phi array
  if(verbose_) std::cout << "Generating radial leaf geometry" << std::endl;
	double scaling_factor = leaf->getParameter("areaMax") * leaf->getLength(false) / leaf->getParameter("k");

	// resolution
	int resolution = leaf_resolution_;
	// Compute the mid vein of the leaf
	CatmullRomSplineManager midVein(leaf->getNodes(), leaf->getNodeIds());
  midVein.setAlpha(0.0, midVein.splineSize() - 1);
	// Compute the leaf length
	auto length = leaf->getLength(false);
	// save c_o
	unsigned int start_c_o = c_o;
  // std::cout << "Accessing leaf random parameter for leaf " << leaf->getId() << std::endl;
	// get leaf random parameter
	auto lrp = leaf->getLeafRandomParameter();
	auto stem = leaf->getParent();
  Vector3d origin = leaf->getOrigin();
  Vector3d origin_neighbor;
  if(leaf->parentNI < stem->getNumberOfNodes() - 1)
    origin_neighbor = stem->getNode(leaf->parentNI + 1);
  else 
    origin_neighbor = stem->getNode(leaf->parentNI - 1);
  Vector3d direction = (origin_neighbor - origin).normalized();
  if(verbose_)
  {
    bool has_petiole = false;
    for (int i = 0; i < leaf->getNumberOfNodes(); ++i)
    {
      if(!leaf->nodeLeafVis(leaf->getLength(i)))
      {
        has_petiole = true;
        break;
      }
    }
    if(has_petiole)
    {
      std::cout << "Leaf has a petiole" << std::endl;
    }
    else
    {
      std::cout << "Leaf has no petiole" << std::endl;
    }
  }

  // calculate the first right vector based on the angle between x axis and the tip of the leaf
  const Vector3d last_node = leaf->getNode(leaf->getNumberOfNodes() - 1);
  const double phi = std::atan2(last_node.y, last_node.x);
  Vector3d right = Vector3d(std::sin(phi), std::cos(phi), 0);
  // calculate the up vector
  Vector3d up = direction.cross(right).normalized();
  right = up.cross(direction).normalized();

  //Quaternion lastRotation = Quaternion::FromMatrix3d({direction, right, up});
	auto min_radius = (use_stem_radius_as_min_) ? stem->getParameter("radius") : this->leaf_minimum_width_;

  // std::cout << "Invoking create leaf radial geometry_ for leaf " << leaf->getId() << std::endl;

	// create leaf radial geometry_
	// greate points for the leaf outer
	// double check that this was not done in python
	if(lrp->leafGeometry.size() == 0)
	{
    if(verbose_) std::cout << "Needing to build leaf geometry nodes" << std::endl;
		lrp->createLeafRadialGeometry(lrp->leafGeometryPhi,lrp->leafGeometryX,resolution);
	}
	// retrieve the leaf geometry_
	auto& outer_geometry_points = lrp->leafGeometry;
  // set buffer sizes
  int point_buffer = 0;
  int index_buffer = 0;
  int last_amount = -1;
  int last_non_petiole = -1;
  double r_max = std::numeric_limits<float>::lowest();
  //std::cout << "Counting how much space we need for the leaf geometry_" << std::endl;
  for (auto i = 0; i < outer_geometry_points.size(); ++i)
  {
		MirrorIterator helper(&(outer_geometry_points[i]));
		auto current_amount = helper.size();
		r_max = std::max(r_max, *std::max_element(outer_geometry_points[i].begin(), outer_geometry_points[i].end()));
    if(current_amount < 2)
    {
      //std::cout << "Skipping petiole at " << i << " because it has size " << current_amount << std::endl;
			last_non_petiole = -1;
      continue;
    }
    //std::cout << "NOT Skipping petiole at " << i << " because it has size " << current_amount << std::endl;
    point_buffer += current_amount;
    if(last_non_petiole > -1 && i > last_non_petiole)
    {
				//std::cout << "Adding Ts : " << last_amount << "/" << current_amount << std::endl;
			if(last_amount != current_amount)
			{
				index_buffer += (std::min(last_amount, (int)current_amount) - 1) * 6;
				index_buffer += (std::abs(last_amount - (int)current_amount) - 1) * 3;
				point_buffer ++;
			}
			else
			{
				index_buffer += (current_amount - 1) * 6;
			}
    }
    if(current_amount > 1)
    {
      last_non_petiole = i;
			last_amount = current_amount;
    }
    if(i > last_non_petiole)
    {
      last_amount = current_amount;
    }
  }
  if(verbose_) std::cout << "Resizing geometry_ buffers by " << point_buffer << " points and " << index_buffer << " triangle values." << std::endl;
  // increase geometry_ buffers
  this->geometry_.resize(std::max(static_cast<std::size_t>(p_o + point_buffer * 3), this->geometry_.size()),-1.0);
	this->geometry_indices_.resize(std::max(static_cast<std::size_t>(c_o + index_buffer), this->geometry_indices_.size()),static_cast<unsigned int>(-1));
	this->geometry_normals_.resize(std::max(static_cast<std::size_t>(p_o + point_buffer * 3), this->geometry_normals_.size()),-1.0);
	this->geometry_texture_coordinates_.resize(std::max(static_cast<std::size_t>((p_o/3*2) + point_buffer * 2), this->geometry_texture_coordinates_.size()),-1.0);
	this->geometry_node_ids_.resize(std::max(static_cast<std::size_t>(p_o / 3 + point_buffer), this->geometry_node_ids_.size()),-1);
	// get the number of points
  // std::cout << "Iterating through the line intersections and generating the geometry_" << std::endl;
  last_amount = -1;
  last_non_petiole = -1;
	Quaternion last_orientation;
	Vector3d last_position;
	int last_index{-1};
  // create two random factors between 0 and 1
  float random_factor_1 = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
  float random_factor_2 = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);

  // computing the first quaternion from: stem is forward, right is angular direction and up is the cross product


	for(auto i = 0; i < outer_geometry_points.size(); ++i)
	{
    std::vector<double>& current = outer_geometry_points[i];

    // if we have a shape function, the current leaf geometry array should be
    // at maximum 1
    if (this->shape_function_.has_value())
    {
      auto max = *std::max_element(current.begin(), current.end());
      for (auto& c : current)
      {
        c = c /  max;
      }
    }

		MirrorIterator helper(&current);
		auto current_amount = helper.size();
    if(current_amount < 2)
    {
      if(verbose_ && i > 0)
      {
        std::cout << "Skipping petiole at " << i << "/" << outer_geometry_points.size()
        << " because it has size " << current_amount << std::endl;
      }
      last_non_petiole = -1;
      continue;
    }
		// compute the size of the current array, which is double its size unless one point is near zero which is only counted once
		int current_size = current_amount;
		// get the current point
    double t = static_cast<double>(i) / static_cast<double>(outer_geometry_points.size()-1);
		double l = t * length;
    auto midpoint = (i == 0) ? leaf->getNode(0) : midVein(t);
    auto current_nodeid = midVein.selectSpline(t).closestId(t);
    midpoint = (i == outer_geometry_points.size() - 1) ? leaf->getNode(leaf->getNumberOfNodes() - 1) : midpoint;
    // get the current point
		// get the best spline for the current point
		auto select_spline = midVein.selectSpline(t);
		
	  // input points, normaly, ids, texture coordinates
    // iterate through the points
		// const Vector3d derivative = select_spline.derivative(t);
		// Quaternion local_q = Quaternion::FromForwardAndRight(derivative, {0,0,1});
    // Vector3d local_right = local_q.Right();
		// auto up = local_q.Up();
    // find out whether right or up is more similar to {0,0,1}

    // local axis
    direction = midVein.derivative(t);
    right = up.cross(direction).normalized();
		right.z = right_penalty_ * right.z;
		right.normalize();
    up = direction.cross(right).normalized();
    direction.normalize();

    // iterate through the points
    //std::cout << "Iterating through the points of the current line intersection " << i << std::endl;

		float petiole_distance = leaf->getLeafRandomParameter()->lb;

		const Vector3d first_node = leaf->getNode(0);
		const Vector3d first_estimated = midVein(0);
		const auto first_distance = (first_node - first_estimated).length();
		//std::cout << "First distance is " << first_distance << std::endl;
    
    for(int p = 0; p < helper.size(); ++p)
    {
      //if(verbose_)
      //  std::cout << p_o << "/" << geometry_.size() << " ";

      auto r = helper[p];
      auto sav = r;
      //auto correction = leaf_minimum_width_ / leaf_width_scale_factor_;
      //r = (helper.isMirrored(p)) ? std::min(r, -correction) : std::max(r, correction);
      // get the point
			// get the wave effect which is a sine function along the length of the leaf

			float z_offset =  0.33 * (helper.isMirrored(p) ? -1.0 : 1.0);
			// make two different sine waves for each side
      if(!add_vertical_leaf_offset_)
      {
        z_offset = 0.0;
      }
			else if(helper.isMirrored(p))
			{
			  z_offset *= std::sin((2.0*random_factor_1 + 2.0) * l / M_PI);
      }
      else
      {
        z_offset *= std::sin((2.0 * random_factor_2 + 2.0) * l / M_PI + M_PI);
      }
			Vector3d base_direction = leaf_width_scale_factor_ * r * right * scaling_factor;
      if (this->shape_function_.has_value())
      {
        base_direction = shape_function_.value()(t) * right * leaf_width_scale_factor_ * r;
      }
			if(leaf_minimum_width_ > 0.0 && (p == 0 || p == helper.size() - 1) && base_direction.length() < leaf_minimum_width_)
      {
        base_direction = right * leaf_minimum_width_ * ((p == 0) ? 1.0 : -1.0);
      }
      else if (use_stem_influence_ && t < stem_influence_radius_)
      {
        base_direction = right * stem->getParameter("radius") * leaf_width_scale_factor_ * r;
      }

			
			//std::cout << base_direction.toString() << std::endl;
			//Vector3d updated_direction = local_q.Rotate(base_direction);
			//base_direction = (base_direction.length() > min_radius) ? base_direction : min_radius * vectorNormalized(base_direction);
			const Vector3d point = midpoint + base_direction +  up * z_offset;
      //ClampVectorBetweenLengths(base_direction, min_radius, 1000.f);
      if(verbose_ && i == outer_geometry_points.size() - 1)
      {
        std::cout << "Base direction: " << base_direction.toString() << std::endl;
        std::cout << "Right: " << right.toString() << std::endl;
        std::cout << "Up: " << up.toString() << std::endl;
        std::cout << "r: " << r << "/" << sav << ", leaf_width_scale_factor_: " << leaf_width_scale_factor_ << ", scaling_factor: " << scaling_factor << std::endl;
        if(this->shape_function_.has_value())
        {
          std::cout << "Shape function: " << this->shape_function_.value()(t) << std::endl;
        }
        std::cout << "Resulting point: " << point.toString() << std::endl;
      }
      //std::cout << "V: " << point.toString() << "; ";
      // set the point
      //std::cout << "p" << " ";
      geometry_[p_o + 0] = point.x;
      geometry_[p_o + 1] = point.y;
      geometry_[p_o + 2] = point.z;
      // set the normal
      //std::cout << "n" << " ";
      geometry_normals_[p_o + 0] = up.x;
      geometry_normals_[p_o + 1] = up.y;
      geometry_normals_[p_o + 2] = up.z;
      // set the texture coordinates
      //std::cout << "t" << " ";
      geometry_texture_coordinates_[(p_o/3*2)] = t;
      geometry_texture_coordinates_[(p_o/3*2) + 1] = helper.texcoord(p);
      if(p == 0) geometry_texture_coordinates_[(p_o/3*2) + 1] = 0.0;
      if(p == helper.size() - 1) geometry_texture_coordinates_[(p_o/3*2) + 1] = 1.0;
			// set the node id
      //std::cout << "i" << " ";
			geometry_node_ids_[p_o/3] = current_nodeid;
			// increase buffer
			p_o += 3;
    }
    //std::cout << std::endl;
    //std::cout << std::endl << "Generating the triangles for the current line intersection " << i << "(" << current.size() << ")" << std::endl;
    if(i > last_non_petiole && last_non_petiole >= 0)
    {
      // use the case distinction between number of intersections
      // for the triangulation between connected sections of the surface
      if(current_amount == last_amount)
      {
        // std::cout << " which is equal to the last one" << std::endl;
        // we construct pairwise triangles for the two sections
        for(auto j = 0; j < current_amount - 1; j += 2)
        {
					// first triangle
					geometry_indices_[c_o++] = (p_o/3) - current_amount + j - last_amount;
					geometry_indices_[c_o++] = (p_o/3) - current_amount + j + 1;
					geometry_indices_[c_o++] = (p_o/3) - current_amount + j;
					// second triangle
					geometry_indices_[c_o++] = (p_o/3) - current_amount + j - last_amount;
					geometry_indices_[c_o++] = (p_o/3) - current_amount + j - last_amount + 1;
					geometry_indices_[c_o++] = (p_o/3) - current_amount + j + 1;
        }
      }
      else
      {
        // set the normal
        geometry_normals_[p_o + 0] = up.x;
        geometry_normals_[p_o + 1] = up.y;
        geometry_normals_[p_o + 2] = up.z;
        // set the texture coordinates
        geometry_texture_coordinates_[(p_o/3*2)] = t;
        geometry_texture_coordinates_[(p_o/3*2) + 0] = 0.0;
        // set the node id
        geometry_node_ids_[p_o/3] = 1.0;
        if(current_amount > last_amount)
        {
					// since we have more points in one of the sections
					// we have to construct triangles with the midvein in mind
					// we construct pairwise triangles for the two sections
					geometry_[p_o + 0] = last_position.x;
					geometry_[p_o + 1] = last_position.y;
					geometry_[p_o + 2] = last_position.z;
          // std::cout << " which is larger to the last one" << std::endl;
          // set the triangles before we increase the buffer to keep the indices correct
          // outer triangles are the only ones connected to the last outline points
          // compute difference between the number of points
          auto diff = current_amount - last_amount;
          // iterate through the top points until we reach the midvein plus the difference
          // std::cout << "iterate through the top points until we reach the midvein plus the difference" << std::endl;
          for(auto j = 0; j < current_amount/2 -diff - 1; ++j)
          {
						// first triangle
						geometry_indices_[c_o++] = (p_o/3) - current_amount + j;
						geometry_indices_[c_o++] = (p_o/3) - current_amount + j + 1;
						geometry_indices_[c_o++] = (p_o/3) - current_amount + j - last_amount;
						// second triangle
						geometry_indices_[c_o++] = (p_o/3) - current_amount + j - last_amount;
						geometry_indices_[c_o++] = (p_o/3) - current_amount + j - last_amount + 1;
						geometry_indices_[c_o++] = (p_o/3) - current_amount + j + 1;
          }
          // iterate through the bottom points starting from the midvein plus the difference
          // std::cout << "iterate through the bottom points starting from the midvein plus the difference" << std::endl;
          for(auto j = current_amount/2+diff; j < current_amount - 1; ++j)
          {
						// first triangle
						geometry_indices_[c_o++] = (p_o/3) - current_amount + j;
						geometry_indices_[c_o++] = (p_o/3) - current_amount + j + 1;
						geometry_indices_[c_o++] = (p_o/3) - current_amount + j - last_amount - diff;
						// second triangle
						geometry_indices_[c_o++] = (p_o/3) - current_amount + j - last_amount - diff;
						geometry_indices_[c_o++] = (p_o/3) - current_amount + j - last_amount + 1 - diff;
						geometry_indices_[c_o++] = (p_o/3) - current_amount + j + 1;
          }
          // iterate through the midvein points
          // std::cout << "iterate through the midvein points" << std::endl;
          // note that this is a specific interpretation of what happens here and might not be correct
          // for a complete correct implementation of this, we would have to start with the phi array
          for(auto j = 0; j < diff - 1; ++j)
          {
            // we triangulate from one point against pairs of points
            // so we only need one triangle, but with the most recent point included
            geometry_indices_[c_o++] = (p_o/3) - current_amount/2 + j;
            geometry_indices_[c_o++] = (p_o/3) - current_amount/2 + j + 1;
            geometry_indices_[c_o++] = (p_o/3);
          }
          // increase buffer
          p_o += 3;
        }
        else
        {
					// since we have more points in one of the sections
					// we have to construct triangles with the midvein in mind
					// we construct pairwise triangles for the two sections
					geometry_[p_o + 0] = midpoint.x;
					geometry_[p_o + 1] = midpoint.y;
					geometry_[p_o + 2] = midpoint.z;
          //std::cout << " which is smaller to the last one" << std::endl;
          auto diff = last_amount - current_amount;
          for(auto j = 0; j < last_amount/2 -diff - 1; ++j)
          {
            geometry_indices_[c_o++] = (p_o/3) - current_amount - last_amount + j;
            geometry_indices_[c_o++] = (p_o/3) - current_amount + j; 
            geometry_indices_[c_o++] = (p_o/3) - current_amount - last_amount + j + 1;
            geometry_indices_[c_o++] = (p_o/3) - current_amount - last_amount + j + 1;
            geometry_indices_[c_o++] = (p_o/3) - current_amount + j;
            geometry_indices_[c_o++] = (p_o/3) - current_amount + j + 1;
          }
          for(auto j = last_amount/2+diff; j < last_amount - 1; ++j)
          {
            geometry_indices_[c_o++] = (p_o/3) - current_amount - last_amount + j;
            geometry_indices_[c_o++] = (p_o/3) - current_amount + j - diff;; 
            geometry_indices_[c_o++] = (p_o/3) - current_amount - last_amount + j + 1; 
            geometry_indices_[c_o++] = (p_o/3) - current_amount + j - diff;
            geometry_indices_[c_o++] = (p_o/3) - current_amount + j - diff + 1; 
            geometry_indices_[c_o++] = (p_o/3) - current_amount - last_amount + j + 1;
          }
          for(auto j = 0; j < diff - 1; ++j)
          {
            geometry_indices_[c_o++] = (p_o/3) - current_amount - last_amount/2 + j;
            geometry_indices_[c_o++] = (p_o/3); 
            geometry_indices_[c_o++] = (p_o/3) - current_amount - last_amount/2 + j + 1;
          }
          p_o += 3;
        }
      }
    }
    if(current_amount > 1)
    {
      last_non_petiole = i;
			last_amount = current_amount;
		}
		//last_orientation = local_q;
		last_position = midpoint;
	}
	//std::cout << "In the end I ended up adding " << c_o - start_c_o << " where I thought I'd add " << index_buffer << std::endl;
}

std::vector<double> PlantVisualiser::GetGeometry() const
{
  return geometry_;
}

std::vector<unsigned int> PlantVisualiser::GetGeometryIndices() const
{
  return geometry_indices_;
}

std::vector<double> PlantVisualiser::GetGeometryNormals() const
{
  return geometry_normals_;
}

std::vector<unsigned short> PlantVisualiser::GetGeometryColors() const
{
  return geometry_colors_;
}

std::vector<double> PlantVisualiser::GetGeometryTextureCoordinates() const
{
  return geometry_texture_coordinates_;
}

std::vector<int> PlantVisualiser::GetGeometryNodeIds() const
{
  return geometry_node_ids_;
}


} // namespace CPlantBox

