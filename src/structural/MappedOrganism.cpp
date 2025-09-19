// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "MappedOrganism.h"

#include "SegmentAnalyser.h"
#include "growth.h"
#include <algorithm>
#include <functional>
#include <cmath>

namespace CPlantBox {

/**
 * A static plant, as needed for flux computations, represented as
 *
 * @param nodes     	coordinates [cm]
 * @param nodeCTs   	node creation times [d]
 * @param segs      	describes a segment with two node indices segs.x, and segs.y [1]
 * @param radii     	segment radius [cm]
 * @param subTypes     	root type or order of the segment [1]
 * @param organTypes    organ type [1]
 */
MappedSegments::MappedSegments(std::vector<Vector3d> nodes, std::vector<double> nodeCTs, std::vector<Vector2i> segs,
		std::vector<double> radii, std::vector<int> subTypes, std::vector<int> organTypes) :
									nodes(nodes), nodeCTs(nodeCTs), segments(segs), radii(radii), subTypes(subTypes), organTypes(organTypes)
{
	assert((nodes.size()==nodeCTs.size()) && "MappedSegments::MappedSegments: Unequal vector sizes nodes and nodeCTs");
	assert((segments.size()==radii.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and radii");
	assert((segments.size()==subTypes.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and subTypes");
	assert((segments.size()==organTypes.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and organ types");
}

/**
 * A static root system, as needed for flux computations, represented as
 *
 * @param nodes     	coordinates [cm]
 * @param nodeCTs   	node creation times [d]
 * @param segs      	describes a segment with two node indices segs.x, and segs.y [1]
 * @param radii     	segment radius [cm]
 * @param subTypes     	root type or order of the segment [1]
 */
MappedSegments::MappedSegments(std::vector<Vector3d> nodes, std::vector<double> nodeCTs, std::vector<Vector2i> segs,
		std::vector<double> radii, std::vector<int> subTypes) :
									nodes(nodes), nodeCTs(nodeCTs), segments(segs), radii(radii), subTypes(subTypes)
{
	organTypes.resize(segments.size());
	std::fill(organTypes.begin(), organTypes.end(), Organism::ot_root);
	assert((nodes.size()==nodeCTs.size()) && "MappedSegments::MappedSegments: Unequal vector sizes nodes and nodeCTs");
	assert((segments.size()==radii.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and radii");
	assert((segments.size()==subTypes.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and subTypes");
	assert((segments.size()==organTypes.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and organ types");
}
/**
 *  A static root system, as needed for flux computations.
 *
 * @param nodes     [cm]
 * @param segs      describes a segment with two node indices segs.x, and segs.y [1]
 * @param radii     [cm] segment radius
 *
 * nodeCTs is set to 0. for all segments
 * subTypes are set to type 0 for all segments
 */

MappedSegments::MappedSegments(std::vector<Vector3d> nodes, std::vector<Vector2i> segs, std::vector<double> radii)
:nodes(nodes), segments(segs), radii(radii) {
	nodeCTs.resize(nodes.size());
	std::fill(nodeCTs.begin(), nodeCTs.end(), 0.);
	organTypes.resize(segments.size());
	std::fill(organTypes.begin(), organTypes.end(), Organism::ot_root);
	setSubTypes(0);
    assert((nodes.size()==nodeCTs.size()) && "MappedSegments::MappedSegments: Unequal vector sizes nodes and nodeCTs");
	assert((segments.size()==radii.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and radii");
	assert((segments.size()==subTypes.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and subTypes");
	assert((segments.size()==organTypes.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and organTypes");
}

/**
 * Sets the radius @param a [cm] for all segments.
 */
void MappedSegments::setRadius(double a) {
	radii.resize(segments.size());
	std::fill(radii.begin(), radii.end(), a);
}

/**
 * Sets the sub type of all segments to @param t.
 */
void MappedSegments::setSubTypes(int t) {
	subTypes.resize(segments.size());
	std::fill(subTypes.begin(), subTypes.end(), t);
}

/**
 * Sets the soil cell index call back function, resets and updates the mappers.
 *
 * @param s 		the callback function picks a cell with spatial coordinate [cm] and returns the index of the cell [1]
 */
void MappedSegments::setSoilGrid(const std::function<int(double,double,double)>& s, bool noChanges) {
	soil_index = s;
	seg2cell.clear(); // re-map all segments
	cell2seg.clear();
	constantLoc = noChanges;
	mapSegments(segments);
}

/**
 * Sets a rectangular grid, and cuts all segments along the grid cells
 *
 * @param min 		minimum of the soil domain [cm]
 * @param max		maximum of the soil domain [cm]
 * @param res	 	resolution, how many cells in each dimension [1]
 * @param cut 		determines if the segments are cut at the rectangular grid faces
 * @param noChanges	determines if the segments remain in the soil voxel they are created in
 */
void MappedSegments::setRectangularGrid(Vector3d min, Vector3d max, Vector3d res, bool cut, bool noChanges)
{
	minBound = min;
	maxBound = max;
	resolution = res;
	cutAtGrid = cut;
	constantLoc = noChanges;
	// std::cout << "setRectangularGrid: cutSegments \n" << std::flush;
	if (cutAtGrid) {
		cutSegments(); // re-add (for cutting)
	}
	// std::cout << "setRectangularGrid: sort \n" << std::flush;
	sort(); // todo should not be necessary, or only in case of cutting?
	seg2cell.clear(); // re-map all segments
	cell2seg.clear();
	// std::cout << "setRectangularGrid: map \n" << std::flush;
	mapSegments(segments);
}

/**
 * Update the mappers root2cell, which maps root segment index to soil cell index, and
 * cell2seg which maps soil cell index to multiple root segments.
 *
 * @param segs      the (new) segments that need to be mapped
 */
void MappedSegments::mapSegments(const std::vector<Vector2i>& segs) {
	for (auto& ns : segs) {
	    //std::cout<< "mapSegments():"<< ns.x <<", " << ns.y << "\n" << std::flush;
		Vector3d mid = (nodes[ns.x].plus(nodes[ns.y])).times(0.5);
		int cellIdx = soil_index(mid.x,mid.y,mid.z);
		int segIdx = ns.y-1; // this is unique in a tree like structured
		seg2cell[segIdx] = cellIdx;
		if (cell2seg.count(cellIdx)>0) {
			cell2seg[cellIdx].push_back(segIdx);
		} else {
			cell2seg[cellIdx] = std::vector<int>({segIdx});
		}
	}
}

/**
 * Cuts segments @param segs at a rectangular grid (@see MappedSegments::setSoilGrid)
 */
void MappedSegments::cutSegments() {
	assert(segments.size()==radii.size() && "MappedSegments::addSegments: number of segments and radii disagree!");
	assert(segments.size()==subTypes.size() && "MappedSegments::addSegments: number of segments and subTypes disagree!");
	assert(segments.size()==organTypes.size() && "MappedSegments::addSegments: number of segments and organTypes disagree!");
	int n = segments.size(); // segs.size() will change within the loop (recursive implementation)
	for (int i=0; i<n; i++ ) {
		addCutSegment(segments[i], radii[i], subTypes[i], organTypes[i],
		    leafBladeSurface[i], segVol[i], bladeLength[i], segO[i], i);
	}
}

/**
 * Adds and cuts a single segment at index @param ii. If the segment is cut, appends the remaining segments.
 * Used by cutSegments
 * This approach may run into problems if a segment is located exactly along a face.
 *
 * @param ns 		the segment to add and cut
 * @param r 		segment radius [cm]
 * @param st 		segment sub type
 * @param ot 		segment organ type
 * @param ii		index to insert the segment, -1 to append the segment
 */
void MappedSegments::addCutSegment(Vector2i ns, double r,  int st, int ot, double lbsurf, double vols, double blen, std::weak_ptr<Organ> segO, int ii) {
	Vector3d n1 = nodes[ns.x];
	Vector3d n2 = nodes[ns.y];
	Vector3d mid = (n1.plus(n2)).times(0.5);
	int im = soil_index(mid.x,mid.y,mid.z); // cell indices
	int in1 = soil_index(n1.x,n1.y,n1.z);
	int in2 = soil_index(n2.x,n2.y,n2.z);
	if ((im!=in1) || (im!=in2)) { // cut
		// build SDF
		auto width = maxBound.minus(minBound); // construct sdf
		Vector3d dx(width.x/resolution.x, width.y/resolution.y, width.z/resolution.z);
		auto mid0 = mid.minus(minBound);
		int x = std::floor(mid0.x/dx.x);
		int y = std::floor(mid0.y/dx.y);
		int z = std::floor(mid0.z/dx.z);
		SDF_Cuboid sdf; // create a signed distance function for cutting
		Vector3d minB(x*dx.x, y*dx.y, z*dx.z);
		minB = minB.plus(minBound);
		Vector3d maxB((x+1)*dx.x, (y+1)*dx.y, (z+1)*dx.z);
		maxB = maxB.plus(minBound);
		sdf.min = minB;
		sdf.max = maxB;
		// std::cout << minB.toString() << ", " << maxB.toString() << ", width " << width.toString() << ", " << sdf.min.toString() << ", " << sdf.max.toString() << "\n";

		im = sdf.getDist(mid)>0; // redo indices, since accuracy of pickking may differ
		in1 = sdf.getDist(n1)>0;
		in2 = sdf.getDist(n2)>0;
		if ((im!=in1) || (im!=in2)) {
			Vector3d cPoint;
			if (im==in1) { // is one node at mid (sort accordingly)
				// std::cout << "n1 " << sdf.getDist(n1) << " mid " << sdf.getDist(mid) << " n2 " << sdf.getDist(n2) << ",indices "<< in1 << ", " << im << ", " << in2 << "\n";
				if (sdf.getDist(n2)<0) {
					cPoint = SegmentAnalyser::cut(n2, mid, std::make_shared<SDF_Cuboid>(sdf), eps);
				} else {
					cPoint = SegmentAnalyser::cut(mid, n2, std::make_shared<SDF_Cuboid>(sdf), eps);
				}
				nodeCTs.push_back(nodeCTs[ns.y]); // todo: we might linearly interpolate
			} else if (im==in2) {
				// std::cout << "n2 " << sdf.getDist(n2) << " mid " << sdf.getDist(mid) << " n1 " << sdf.getDist(n1) << ", " << in2 << ", " << im << ", " << in1 << "\n";
				if (sdf.getDist(n1)<0) {
					cPoint = SegmentAnalyser::cut(n1, mid, std::make_shared<SDF_Cuboid>(sdf), eps);
				} else {
					cPoint = SegmentAnalyser::cut(mid, n1, std::make_shared<SDF_Cuboid>(sdf), eps);
				}
				nodeCTs.push_back(nodeCTs[ns.x]); // todo: we might linearly interpolate
			} else { // otherwise split in mid, use cutSegments on those
				cPoint = mid;
				nodeCTs.push_back(0.5*(nodeCTs[ns.x]+nodeCTs[ns.y]));
			}
			// std::cout << "[" << n1.toString() << n2.toString() << "] -> [" << nodes[ns.x].toString() << ", " << nodes.back().toString() << "], ["<< nodes.back().toString() << ", " << n2.toString() << "], " << "\n";
			nodes.push_back(cPoint);
			Vector2i s1(ns.x, nodes.size()-1);
			Vector2i s2(nodes.size()-1, ns.y);
			if ((length(s1)<eps) ||  (length(s2)<eps)) { // if the cut segments are too small, just give up
				addSegment(ns, r, st, ot, lbsurf, vols, blen, segO, ii);
				nodes.pop_back(); // remove cPoint
				nodeCTs.pop_back();
			} else {
				addCutSegment(s1, r, st , ot, lbsurf, vols, blen, segO, ii); // first segment replaces at index ii
				addCutSegment(s2, r, st , ot, lbsurf, vols, blen, segO, -1); // append second segment
			}
		} else { // im==in1==in2, dont't cut
			// std::cout << "ok " << ii <<": (" << ns.x <<", " << ns.y << ") [" << n1.toString() <<", "<< n2.toString() <<"\n";
			addSegment(ns, r, st, ot, lbsurf, vols, blen, segO, ii);
		}
	} else { // im==in1==in2, dont't cut
		// std::cout << "ok " << ii <<": (" << ns.x <<", " << ns.y << ") [" << n1.toString() <<", "<< n2.toString() <<"\n";
		addSegment(ns, r, st, ot, lbsurf, vols, blen, segO, ii);
	}
}

/**
 * Adds the segment at index i, or appends it, if i = -1
 *
 * @param s 		the segment to append or insert
 * @param r 		segment radius [cm]
 * @param st 		segment sub type
 * @param ot 		segment organ type
 * @param i			index to insert the segment, -1 to append the segment
 */
void MappedSegments::addSegment(Vector2i s, double a,  int st, int ot, double lbsurf, double vols, double blen, std::weak_ptr<Organ> o, int i) {
	if (i>=0) {
		segments.at(i) = s;
		radii.at(i) = a;
		subTypes.at(i) = st;
		organTypes.at(i) = ot;
        leafBladeSurface.at(i) = lbsurf;
        segVol.at(i) = vols;
        bladeLength.at(i) = blen;
        segO.at(i) = o;
	} else {
		segments.push_back(s);
		radii.push_back(a);
		subTypes.push_back(st);
		organTypes.push_back(ot);
		leafBladeSurface.push_back(lbsurf);
		segVol.push_back(vols);
		bladeLength.push_back(blen);
		segO.push_back(o);
	}
}

/**
 * Length of the segment @param s
 */
double MappedSegments::length(const Vector2i& s) const {
	return (nodes.at(s.y).minus(nodes.at(s.x))).length();
}

/**
 * Removes segments @param segs from the mappers
 */
void MappedSegments::unmapSegments(const std::vector<Vector2i>& segs) {
	for (auto& ns : segs) {
		int cellIdx = -1;
		int segIdx = ns.y-1;
		if (seg2cell.count(segIdx)>0) { // remove from seg2cell
			cellIdx = seg2cell[segIdx];
			auto it = seg2cell.find(segIdx);
			seg2cell.erase(it);
		} else {
			throw std::invalid_argument("MappedSegments::removeSegments: warning segment index "+ std::to_string(segIdx)+ " was not found in the seg2cell mapper");
		}
		if (cell2seg.count(cellIdx)>0) {
			auto& csegs= cell2seg[cellIdx];
			int c = 0;
			for (int i=0; i<csegs.size(); i++) {
				if (csegs[i] == segIdx) {
					csegs.erase(csegs.begin() + c, csegs.begin() + c + 1);
					break; // inner for
				}
				c++;
			}
		} else {
			throw std::invalid_argument("MappedSegments::removeSegments: warning cell index "+ std::to_string(cellIdx)+ " was not found in the cell2seg mapper");
		}
	}
}



/**
 * Maps a point into a cell and return the cells linear index (for a equidistant rectangular domain)
 */
int MappedSegments::soil_index_(double x, double y, double z) {
	Vector3d p(x,y,z);
	std::array<double,3>  r = { resolution.x, resolution.y, resolution.z};
	auto w = maxBound.minus(minBound);
	auto p0 = p.minus(minBound);
	std::array<double,3> i = { p0.x/w.x*r[0], p0.y/w.y*r[1], p0.z/w.z*r[2] };
	for (int k=0; k<3; k++) {
		if ((i[k] < 0) || (i[k] >= r[k])) {
			return -1; // point is out of domain
		}
	}
	return std::floor(i[2]) * r[0] * r[1] + std::floor(i[1]) * r[0] + std::floor(i[0]); // a linear index not periodic
}

/**
 * Sorts the segments, so that the segment index == second node index -1 (unique mapping in a tree)
 */
void MappedSegments::sort() {
	auto newSegs = segments;
	auto newRadii = radii;
	auto newSubTypes = subTypes;
	auto newTypesorgan = organTypes;
	for (int i=0; i<newSegs.size(); i++) {
		int ind = segments[i].y-1;
		newSegs[ind] = segments[i];
		newRadii[ind] = radii[i];
		newSubTypes[ind] = subTypes[i];
		newTypesorgan[ind] = organTypes[i];
	}
	segments = newSegs;
	radii = newRadii;
	subTypes = newSubTypes;
	organTypes = newTypesorgan;
}


/**
 * Sums segment fluxes over each cell
 *
 * @param segFluxes 	segment fluxes given per segment index [cm3/day]
 * @return hash map with cell indices as keys and fluxes as values [cm3/day]
 */
std::map<int,double> MappedSegments::sumSegFluxes(const std::vector<double>& segFluxes)
{
    std::map<int,double> fluxes;
    for (int si = 0; si<this->segments.size(); si++) {
        int j = this->segments[si].y;
        int segIdx = j-1;

        if (this->seg2cell.count(segIdx)>0)
		{
			int cellIdx = this->seg2cell[segIdx];
			if (cellIdx>=0)
			{
				if(this->organTypes[segIdx] == Organism::ot_root)//only divid the fluxes between the root segments
				{
					if (fluxes.count(cellIdx)==0) {
						fluxes[cellIdx] = segFluxes[segIdx];
					} else {
						fluxes[cellIdx] = fluxes[cellIdx] + segFluxes[segIdx]; // sum up fluxes per cell
					}
				}else{
					if(segFluxes[segIdx] != 0.)
					{
						std::stringstream errMsg;
						errMsg<<"MappedSegments::sumSegFluxes. ot:"<<this->organTypes[segIdx]<<" segIdx:"<<segIdx
						<<" cellIdx:"<<cellIdx<<" segFluxes[segIdx] :"
						<<segFluxes[segIdx]<<"=> shoot segment bellow ground ans exchanges water" <<std::endl;

						throw std::runtime_error(errMsg.str().c_str());
					}
				}
			}
        }
    }
    return fluxes;
}

/**
 * Splits soil fluxes per cell to the segments within the cell, so that the summed fluxes agree, @see sumSoilFluxes()
 *
 * @param soilFluxes 	cell fluxes per global index [cm3/day]
 * @param type 			split flux proportional to 0: segment volume, 1: segment surface, 2: segment length
 * @return fluxes for each segment [cm3/day]
 */
std::vector<double> MappedSegments::splitSoilFluxes(const std::vector<double>& soilFluxes, int type) const
{
    auto lengths =  this->segLength();
    std::vector<double> fluxes = std::vector<double>(this->segments.size());
    std::fill(fluxes.begin(), fluxes.end(), 0.);
    auto map = this->cell2seg;
	double fluxesTotTot =0;
    for(auto iter = map.begin(); iter != map.end(); ++iter) {
        int cellId =  iter->first;
        auto segs = map.at(cellId);
		if (cellId>=0) {
			double v = 0.;  // calculate sum over cell
			for (int i : segs) {

				if(this->organTypes[i] == Organism::ot_root)//only divid the fluxes between the root segments
				{
					if (type==0) { // volume
						v += M_PI*(this->radii[i]*this->radii[i])*lengths[i];
					} else if (type==1) { // surface
						v += 2*M_PI*this->radii[i]*lengths[i];
					} else if (type==2) { // length
						v += lengths[i];
					}
				}
			}
			double fluxesTot = 0;
			for (int i : segs) { // calculate outer radius

				if(this->organTypes[i] == Organism::ot_root)
				{
					double t =0.; // proportionality factor (must sum up to == 1 over cell)
					if (type==0) { // volume
						t = M_PI*(this->radii[i]*this->radii[i])*lengths[i]/v;
					} else if (type==1) { // surface
						t = 2*M_PI*this->radii[i]*lengths[i]/v;
					} else if (type==2) { // length
						t = lengths[i]/v;
					}
					if(fluxes[i] !=0){std::cout<<"fluxes "<<i<<" already set "<<std::endl;assert(false);}
					fluxes[i] = t*soilFluxes.at(cellId);
					fluxesTot +=  t*soilFluxes.at(cellId);
					fluxesTotTot += t*soilFluxes.at(cellId);
				}
			}
		}
    }
    return fluxes;
}


/**
 * Calculates segment lengths [cm]
 */
std::vector<double> MappedSegments::segLength() const {
	std::vector<double> lengths = std::vector<double>(segments.size());
	for(int i=0; i<lengths.size(); i++) {
		auto n1 = nodes[segments[i].x];
		auto n2 = nodes[segments[i].y];
		lengths[i] = (n2.minus(n1)).length();
	}
	return lengths;
}

/**
 * Returns soil matric potential per segment, for a given soil sx connected gy the mapper rs->seg2cell [TODO make more general perSegment(value_perCell)...]
 */
std::vector<double> MappedSegments::getHs(const std::vector<double> sx) const {
    double psi_air = -954378;
    std::vector<double> hs = std::vector<double>(this->segments.size());
    for (int si = 0; si<this->segments.size(); si++) {
        int cellIndex = this->seg2cell.at(si);
        if (cellIndex>=0) {
            if(sx.size()>1) {
                hs[si] = sx.at(cellIndex);
            } else {
                hs[si] = sx.at(0);
            }
        } else {
            hs[si] = psi_air;
        }
    }
    return hs;
}

/**
 * Calculates the z-coordinates of the segment
 */
std::vector<double> MappedSegments::getSegmentZ() const {
    std::vector<double> z = std::vector<double>(segments.size());
    for (int i=0; i<z.size(); i++) {
        z[i] = nodes[segments[i].y].z; // 0.5*(nodes[segments[i].x].z + nodes[segments[i].y].z);
    }
    return z;
}

/**
 * Calculates the total potential from the matric potential
 */
std::vector<double> MappedSegments::matric2total(std::vector<double> sx) const {
    std::vector<double> b = this->getSegmentZ();
    assert(sx.size() == b.size());
    std::transform(sx.begin( ), sx.end( ), b.begin( ), sx.begin( ),std::plus<double>( ));
    return sx;
}

/**
 * Calculates the matric potential from the tortal potential
 */
std::vector<double> MappedSegments::total2matric(std::vector<double> sx) const{
    std::vector<double> b = this->getSegmentZ();
    assert(sx.size() == b.size());
    std::transform(sx.begin( ), sx.end( ), b.begin( ), sx.begin( ),std::minus<double>( ));
    return sx;
}



/**
 * Returns seg2cell as vector
 */
std::vector<int> MappedSegments::getSegmentMapper() const {
    std::vector<int> mapper = std::vector<int>(segments.size());
    for (int i=0; i<mapper.size(); i++) {
        try {
            mapper[i] = seg2cell.at(i);
        } catch(...) {
            std::cout << "MappedSegments::getSegmentMapper(): Index "<< i << " not mapped\n" << std::flush;
            throw;
        }
    }
    return mapper;
}

/**
 * Calculates the minimum of node coordinates
 * (e.g. minimum corner of bounding box)
 * value not cached
 */
Vector3d MappedSegments::getMinBounds() {
    Vector3d min_ = Vector3d(nodes[0].x, nodes[0].y, nodes[0].z);
    for (const auto& n : nodes) {
        if (n.x < min_.x) {
            min_.x = n.x;
        }
        if (n.y < min_.y) {
            min_.y = n.y;
        }
        if (n.z < min_.z) {
            min_.z = n.z;
        }
    }
    return min_;
}

int MappedSegments::getSegment2leafId(int si_){
		throw std::runtime_error("MappedSegments::getsegment2leafId: tried to access leafId of");
		return -1;
}

std::vector<double> MappedSegments::getEffectiveRadii() {
	int n = radii.size();
	std::vector<double> radii_(n);
	for (int i = 0; i<n; i++) {
		radii_.at(i) = this->getEffectiveRadius(i);
	}
	return radii_;
}









/**
 * initialization of mappedplant
 * @param verbose 		indicates if status is written to the console (cout) (default = false)
 * @param LB		 	implement length-based waiting time before growth (true) of laterals or delay-based (false)? (default = true)
 *
 *
 * use setStochastic(stochastic); to turn on or off stochasticity
 */
void MappedPlant::initialize_(bool verbose, bool lengthBased) {

    reset(); // just in case (Plant::reset()) (careful, MappedPlant cannot reset, yet)
	auto stemP = getOrganRandomParameter(Organism::ot_stem);
	bool plantBox = stemP.size()>1; // prototype + a real parameter definition
	if ((extraNode == -1) && (plantBox)) {
		disableExtraNode(); // no meed for additional node to create the artificial stem
	} else {
	    enableExtraNode();
	}
	if (lengthBased){
	    if(verbose) {
	        std::cout << "MappedPlant::initializeLB \n" << std::flush;
	    }
	    Plant::initializeLB(verbose); // initializes plant
	} else {
        if(verbose) {
            std::cout << "MappedPlant::initializeDB \n" << std::flush;
        }
	    Plant::initializeDB(verbose); // initializes plant
	}
	if (extraNode==1) { // inserts a special seed segment (for root system only hydraulic simualtions)
        auto initial_nodes = this->getNodes();
        auto initial_ncts = this->getNodeCTs();
        auto n0 = initial_nodes.at(0);
        auto ct0 = initial_ncts.at(0);
        n0.x += 0.1; // 0.1 cm next to seed pos
        nodes.push_back(n0);
        for (auto n : initial_nodes) {
            nodes.push_back(n);
        }
        nodeCTs.push_back(ct0);
        for (auto ct : initial_ncts) {
            nodeCTs.push_back(ct);
        }
        addSegment(Vector2i(0,1), 0.1, 0, Organism::OrganTypes::ot_stem, 0, 0, 0, getSeed(), -1 ); // (Vector2i s, double a,  int st, int ot, double lbsurf, double vols, double blen, std::weak_ptr<Organ> o, int i)
	} else {
        nodes = this->getNodes();
        nodeCTs = this->getNodeCTs();
	}

	mapSegments(segments);
	mapSubTypes();
}

/**
 * creates a map to reorder sub types, so that
 * the N subtypes of one organ type go from 0 to N-1
 */
void MappedPlant::mapSubTypes(){
	for(int ot = 0; ot < organParam.size(); ot++) {
		//std::cout<<"MappedPlant::mapSubTypes for organtype "<<ot<<" with "<<organParam[ot].size()<<" subtypes "<<std::endl;
		int stNew = 0;
		for(int stOld_ = 1; stOld_ < organParam[ot].size();stOld_++) { //skipe stOld ==0, not realted to any organ st
			if (organParam[ot][stOld_] != NULL) {
				int stOld = organParam[ot][stOld_]->subType;
				st2newst[std::make_tuple(ot, stOld)] = stNew;
				//std::cout<<"old st: "<<stOld<<", new st: "<< stNew <<std::endl;
				stNew ++;
			} // else {std::cout<<"subType n#"<<stOld_<<" does not exist, skip "<<std::endl;}
		}
	}
}

/**
 * Simulates the development of the organism in a time span of @param dt days.
 *
 * @param dt        time step [day]
 * @param verbose   turns console output on or off
 */
void MappedPlant::simulate(double dt, bool verbose)
{
    size_t shift = (extraNode == 1);

	if (soil_index==nullptr) {
		throw std::invalid_argument("MappedPlant::simulate():soil was not set, use MappedPlant::simulate::setSoilGrid" );
	}

	Plant::simulate(dt,  verbose);

	auto uni = this->getUpdatedNodeIndices(); // move nodes
	for (int& i : uni) { // shift
		i += shift;
	}
	auto unodes = this->getUpdatedNodes();
	auto uncts = this->getUpdatedNodeCTs();
	assert(uni.size()==unodes.size() && "updated node indices and number of nodes must be equal");
	int c = 0;
	for (int i : uni) {
		nodes.at(i) = unodes[c];
		nodeCTs.at(i) = uncts[c];
		c++;
	}
	if (verbose) {
		std::cout << "nodes moved "<< uni.size() << "\n" << std::flush;
	}

	auto newnodes = this->getNewNodes(); // add nodes
	nodes.reserve(nodes.size()+newnodes.size());
	for (auto& nn : newnodes) {
		nodes.push_back(nn);
	}
	auto newnode_cts = this->getNewNodeCTs(); // add node cts
	nodeCTs.reserve(nodeCTs.size()+newnode_cts.size());
	for (auto& nct : newnode_cts) {
		nodeCTs.push_back(nct);
	}
	if (verbose) {
		std::cout << "new nodes added " << newnodes.size() << "\n" << std::flush;
	}

	auto newsegs = this->getSegments(); // add segments (TODO ALL) (TODO cutting, more nodes will lead to different shift???)
	segments.resize(newsegs.size()+shift);
	for (auto& ns : newsegs) { // shift
	    ns.x += shift;
	    ns.y += shift;
		segments[ns.y-1] = ns;
	}
	if (verbose) {
		std::cout << "segments added "<< newsegs.size() << "\n" << std::flush;
	}
	auto newsegO = this->getSegmentOrigins(); // (TODO ALL) to add radius and type (TODO cutting)
	radii.resize(newsegO.size()+shift);
	subTypes.resize(newsegO.size()+shift);
	organTypes.resize(newsegO.size()+shift);
	segVol.resize(newsegO.size()+shift);
	bladeLength.resize(newsegO.size()+shift);
	leafBladeSurface.resize(newsegO.size()+shift);
	this->segO.resize(newsegO.size()+shift);
	c = 0;
	if (verbose) {
		std::cout << "Number of segments " << radii.size() << ", including " << newsegO.size() << " new \n"<< std::flush;
	}
	for (auto& so : newsegO) {
		int segIdx = newsegs[c].y-1;

		radii.at(segIdx) = so->param()->a;
		organTypes.at(segIdx) = so->organType();
		subTypes.at(segIdx) = so->param()->subType; //  st2newst[std::make_tuple(organTypes[segIdx],so->param()->subType)];//new st
		this->segO.at(segIdx) = so; // useful when creating SegmentAnalyser from a mappedSegment

		if (organTypes.at(segIdx) == Organism::ot_leaf) //leaves can be cylinder, cuboid or characterized by user-defined 2D shape
		{
			int index;
			auto nodeIds = so->getNodeIds();
		    for (int& i : nodeIds) { // shift
		        i += shift;
		    }
			auto it = find(nodeIds.begin(), nodeIds.end(), newsegs[c].y);
			if (it != nodeIds.end()){ index = it - nodeIds.begin() -1;
			} else {
				throw std::runtime_error("MappedPlant::simulate: global segment index not found in organ");
			}
			int localSegId = index;
			bool realized = true;
			bool withPetiole = false;
			segVol.at(segIdx) = -1;
			bladeLength.at(segIdx) = std::static_pointer_cast<Leaf>(so)->leafLengthAtSeg(localSegId, withPetiole);
			leafBladeSurface.at(segIdx) =  std::static_pointer_cast<Leaf>(so)->leafAreaAtSeg(localSegId,realized, withPetiole);
			withPetiole = true;
			segVol.at(segIdx) = std::static_pointer_cast<Leaf>(so)->leafVolAtSeg(localSegId, realized, withPetiole);//* thickness;
			if(segVol.at(segIdx) < 0) {
				std::stringstream errMsg;
				errMsg <<"MappedPlant::simulate: computation of leaf volume failed "<<segVol.at(segIdx)<<"\n";
				throw std::runtime_error(errMsg.str().c_str());
			}

		} else { //stems and roots are cylinder
			auto s = segments.at(segIdx);
			double length_seg = (nodes.at(s.x).minus(nodes.at(s.y))).length();
			segVol.at(segIdx) = radii.at(segIdx) * radii.at(segIdx) * M_PI * length_seg;
			bladeLength.at(segIdx) = 0;
			leafBladeSurface.at(segIdx) = 0;
		}
		c++;
	}

	// map new segments
	newsegs = this->getNewSegments();
    for (auto& ns : newsegs) { // shift
        ns.x += shift;
        ns.y += shift;
    }
	this->mapSegments(newsegs);

	// update segments of moved nodes
	std::vector<Vector2i> rSegs;
	if(!constantLoc)//for 1d-3d coupling need to have segments remain in the same voxel
	{//also, if soil_index is in parallel, this blocks the program as plant only grows on
		// one thread
		for (int i : uni) {
			int segIdx = i -1;
			int cellIdx = seg2cell[segIdx];
			auto s = segments[segIdx];
			Vector3d mid = (nodes[s.x].plus(nodes[s.y])).times(0.5);
			int newCellIdx = soil_index(mid.x,mid.y,mid.z);
			// 1. check if mid is still in same cell (otherwise, remove, and add again)
			// 2. if cut is on, check if end point is in same cell than mid point (otherwise remove and add again)
			bool remove = false;
			if (cellIdx==newCellIdx) {
				if (cutAtGrid) {
					auto endPoint = nodes[s.y];
					newCellIdx = soil_index(endPoint.x,endPoint.y,endPoint.z);
					remove = (newCellIdx!=cellIdx);
				}
			} else {
				if (!constantLoc) {
					remove = true;
				}
			}
			if (remove) {
				rSegs.push_back(s);
			}
		}
	}
	MappedSegments::unmapSegments(rSegs);
	MappedSegments::mapSegments(rSegs);

	if (kr_length > 0. || rootHairs) {
	    calcExchangeZoneCoefs();
	}

	getSegment2leafIds();

}

/**
 * computes coeficients for kr
 * when root kr > 0 up to kr_length cm from the root tip
 * see @XylemFlux::kr_RootExchangeZonePerType()
 **/
void MappedPlant::calcExchangeZoneCoefs() { //
    size_t shift = (extraNode==1);
	exchangeZoneCoefs.resize(segments.size(), -1.0);
	distanceTip.resize(segments.size(), -1.0);
	distanceTip.at(0) = 0.;
	exchangeZoneCoefs.at(0) =1.;
	auto orgs = getOrgans(-1);
	for(auto org: orgs) {
		for(int localIdx = 1; localIdx < org->getNumberOfNodes();localIdx++) {
			int globalIdx_x = org->getNodeId(localIdx -1 ) + shift; // add shift
			int globalIdx_y = org->getNodeId(localIdx) + shift;
			if(org->organType() != Organism::ot_root){
				exchangeZoneCoefs.at(globalIdx_y-1) = 1;
				distanceTip.at(globalIdx_y-1) = -1;
			} else {
				auto n1 = nodes.at(globalIdx_x);
				auto n2 = nodes.at(globalIdx_y);
				auto v = n2.minus(n1);
				double l = v.length();
				double distance2Tip_y = org->getLength(true) - org->getLength(localIdx);
				double length_in_exchangeZone = std::min(l,std::max(kr_length - std::max(distance2Tip_y,0.),0.));
				exchangeZoneCoefs.at(globalIdx_y-1) = length_in_exchangeZone/l;
				distanceTip.at(globalIdx_y-1) = distance2Tip_y + l/2; //distance of the segment center to the tip
			}
		}
	}
	const int notFound = std::count(exchangeZoneCoefs.cbegin(), exchangeZoneCoefs.cend(), -1.0);
	if (notFound != 0) {
		std::stringstream errMsg;
		errMsg <<"MappedPlant::calcExchangeZoneCoefs(): "<<notFound<<" elements not initalized";
		throw std::runtime_error(errMsg.str().c_str());
		std::cout<<"notFound "<<notFound<<std::endl;
	}
	// std::cout << "distanceTip " << distanceTip.size()  << "\n" << std::flush;
}


/**
 *Gives an overview of the mappedplant object (for debugging)
 *
 **/
void MappedPlant::printNodes() {
	std::cout << "\n MappedPlant::printnodes \n"<< std::flush;
	std::cout << "\n nodes \n"<< std::flush;
	nodes = this->getNodes();
	for (auto nd : nodes) {
		std::cout <<nd.toString()<< std::flush;
	}
	std::cout << "\n nodes size \n" <<  nodes.size() << std::flush;
	std::cout << "\n organ types \n" << std::flush;
	for (auto ot : organTypes) {
		std::cout << ot << std::flush;
	}
	std::cout << "\n organTypes size\n"<<  organTypes.size() << std::flush;
	std::cout << "\n subtypes \n"<< std::flush;
	for (auto st : subTypes) {
		std::cout << st << std::flush;
	}
	std::cout << "\n subtypes size \n"<< subTypes.size() << std::flush;
	std::cout << "\n cts \n"<< std::flush;
	for (auto to : nodeCTs) {
		std::cout << to << std::flush;
	}
	std::cout << "\n cts size \n"<< nodeCTs.size() << std::flush;
	std::cout << "\n segments \n"<< std::flush;
	for (auto to : segments) {
		std::cout << to.toString() << std::flush;
	}
	std::cout << "\n segments size \n"<< segments.size() << std::flush;
}

/**
 *	index of node of organtype ot
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          Id of each segment
 */
std::vector<int> MappedPlant::getSegmentIds(int ot) const
{
	std::vector<int> segId;// = std::vector<int>(segments.size());
    for (int i=0; i<segments.size(); i++) {
		if((ot == -1)||(organTypes[segments[i].y-1]== ot)){
			segId.push_back(segments[i].y-1);
		}
    }
	return segId;
}

/**
 *	index of node of organtype ot
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          Id of each segment
 */
void MappedPlant::getSegment2leafIds()
{
	segment2leafIds = std::vector<int>(organTypes.size(),-1);
	int leafId = 0;
    for (int i=0; i<organTypes.size(); i++) {
		if(organTypes.at(i) == Organism::ot_leaf){
			segment2leafIds.at(i) = leafId;
			leafId +=1;
		}
    }
}

/**
 * index of segment of organ type ot
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          Id of each segment
 */
std::vector<int> MappedPlant::getNodeIds(int ot) const {
	std::vector<int> nodeId;// = std::vector<int>(segments.size());
    for (int i=0; i<segments.size(); i++) {
		if((ot == -1)||(organTypes[segments[i].y-1]== ot)){
			nodeId.push_back(segments[i].y);
		}
    }
	return nodeId;
}

/**
 *  returns the plant radius plus root hair length
 */
double MappedPlant::getEffectiveRadius(int si) {
    int ot = organTypes.at(si);
    if (ot == Organism::ot_root) {
        int st = subTypes.at(si);
        double l = distanceTip.at(si);
        auto rrp = std::static_pointer_cast<RootRandomParameter>(this->getOrganRandomParameter(ot, st));
        double zl = rrp->hairsZone;
        double el = rrp->hairsElongation;
        if (l<(zl+el) && l>el) { // add effective root length
            double hl = rrp->hairsLength;
            return this->radii.at(si)+hl;
        } else {
            return this->radii.at(si);
        }
    } else {
        return this->radii.at(si);
    }
}

/**
 * index of segment of organ type ot
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          Id of each segment
 */
double MappedPlant::getPerimeter(int si_, double l_) {
	if (organTypes.at(si_) == Organism::ot_leaf) {
        //perimeter of the leaf blade
        // "*2" => C3 plant has stomatas on both sides.
        //later make it as option to have C4, i.e., stomatas on one side
        //int leafId = getSegment2leafId(si_);
	    return leafBladeSurface.at(si_) / l_ *2;
    } else {
    	return 2 * M_PI * this->getEffectiveRadius(si_);
    }
}

int MappedPlant::getSegment2leafId(int si_) {
		int leafSi = segment2leafIds.at(si_);
		if(leafSi == -1) {
			throw std::runtime_error("MappedSegments::getsegment2leafId: tried to access leafId of non-leaf segment");
		}
		return leafSi;
	}

} // namespace
