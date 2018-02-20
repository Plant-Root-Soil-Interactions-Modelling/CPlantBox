// *** ADDED BY HEADER FIXUP ***
#include <istream>
// *** END ***
#include "analysis.h"
#include <iostream>
#include <iomanip>

/**
 * Copies the segments of the roots system into the analysis class
 *
 * @param rs      the root system that is analysed
 */
SegmentAnalyser::SegmentAnalyser(const Plant& rs)
{
	nodes = rs.getNodes();
	segments = rs.getSegments();
	ctimes = rs.getNETimes();
	segO = rs.getSegmentsOrigin();
	assert(segments.size()==ctimes.size());
	assert(segments.size()==segO.size());
    this->rs = &rs; //needed for dgf writer, only
}

/**
 * Adds all segmentes from root system @param rs to the analysis.
 */
void SegmentAnalyser::addSegments(const Plant& rs)
{
	addSegments(SegmentAnalyser(rs));
}

/**
 * Adds all segments from the analyser @param a to this analysis.
 */
void SegmentAnalyser::addSegments(const SegmentAnalyser& a)
{
	int offset = nodes.size();
	nodes.insert(nodes.end(),a.nodes.begin(),a.nodes.end()); // copy nodes
	auto ns = a.segments;
	for (auto& s : ns) { // shift indices
		s.x += offset;
		s.y += offset;
	}
	segments.insert(segments.end(),ns.begin(),ns.end()); // copy segments
	ctimes.insert(ctimes.end(),a.ctimes.begin(),a.ctimes.end()); // copy times
	segO.insert(segO.end(),a.segO.begin(),a.segO.end());// copy origins
	assert(segments.size()==ctimes.size());
	assert(segments.size()==segO.size());
}

/**
 * Returns a specific parameter per root segment
 *
 * @param st    parameter type @see RootSystem::ScalarType per segment
 * \return      vector containing parameter value per segment
 */
std::vector<double> SegmentAnalyser::getScalar(std::string name) const
{
	std::vector<double> data(segO.size());

	if (name=="st_time") {
		data = ctimes;
		return data;
	}
//	if (st==Plant::st_userdata1) {
//		data = userData.at(0);
//		return data;
//	}
//	if (st==Plant::st_userdata2) {
//		data = userData.at(1);
//		return data;
//	}
//	if (st==Plant::st_userdata3) {
//		data = userData.at(2);
//		return data;
//	}
//
	for (size_t i=0; i<segO.size(); i++) {
		data.at(i) = segO.at(i)->getScalar(name);
	}
	return data;
}

/**
 * Returns the length of a segment
 *
 * @param i 	index of the segment
 * \return 		the length of segment i
 */
double SegmentAnalyser::getSegmentLength(int i) const
{
	Vector2i s = segments.at(i);
	Vector3d x = nodes.at(s.x);
	Vector3d y = nodes.at(s.y);
	return (x.minus(y)).length();
}

/**
 * Crops the segments with some geometry
 *
 * @param geometry      signed distance function of the geometry
 */
void SegmentAnalyser::crop(SignedDistanceFunction* geometry)
{
	//std::cout << "cropping " << segments.size() << " segments...";
	std::vector<Vector2i> seg;
	std::vector<Organ*> sO;
	std::vector<double> ntimes;
	for (size_t i=0; i<segments.size(); i++) {
		auto s = segments.at(i);
		Vector3d x = nodes.at(s.x);
		Vector3d y = nodes.at(s.y);
		bool x_ = geometry->getDist(x)<=0; // in?
		bool y_ = geometry->getDist(y)<=0; // in?
		if ((x_==true) && (y_==true)) { //segment is inside
			seg.push_back(s);
			sO.push_back(segO.at(i));
			ntimes.push_back(ctimes.at(i));
		} else if ((x_==false) && (y_==false)) { // segment is outside

		} else { // one node is inside, one outside
			// sort
			Vector3d in;
			Vector3d out;
			int ini;
			if (x_==true) {
				in = x;
				ini = s.x;
				out = y;
			} else {
				in = y;
				ini = s.y;
				out = x;
			}
			// cut
			Vector3d newnode = cut(in, out, geometry);
			// add new segment
			nodes.push_back(newnode);
			Vector2i newseg(ini,nodes.size()-1);
			seg.push_back(newseg);
			sO.push_back(segO.at(i));
			ntimes.push_back(ctimes.at(i));
		}

	}
	segments = seg;
	segO  = sO;
	ctimes = ntimes;
	//std::cout << " cropped to " << segments.size() << " segments " << "\n";
}

/**
 * Filters the segments to the ones, where data is within [min,max], @see AnalysisSDF::getData,
 * i.e. all other segments are deleted.
 *
 * @param st    parameter type @see RootSystem::ScalarType
 * @param min   minimal value
 * @param max   maximal value
 */
void SegmentAnalyser::filter(std::string pname, double min, double max)
{
	std::vector<double> data = getScalar(pname);
	std::vector<Vector2i> seg;
	std::vector<Organ*> sO;
	std::vector<double> ntimes;
	for (size_t i=0; i<segments.size(); i++) {
		if ((data.at(i)>=min) && (data.at(i)<=max)) {
			seg.push_back(segments.at(i));
			sO.push_back(segO.at(i));
			ntimes.push_back(ctimes.at(i));
		}
	}
	segments = seg;
	segO  = sO;
	ctimes = ntimes;
}

/**
 * Filters the segments to the ones, where data equals value, @see AnalysisSDF::getData,
 * i.e. all other segments are deleted.
 *
 * @param st        parameter type @see RootSystem::ScalarType
 * @param value     parameter value of the segments that are kept
 */
void SegmentAnalyser::filter(std::string pname, double value)
{
	std::vector<double> data = getScalar(pname);
	std::vector<Vector2i> seg;
	std::vector<Organ*> sO;
	std::vector<double> ntimes;
	for (size_t i=0; i<segments.size(); i++) {
		if (data.at(i)==value) {
			seg.push_back(segments.at(i));
			sO.push_back(segO.at(i));
			ntimes.push_back(ctimes.at(i));
		}
	}
	segments = seg;
	segO  = sO;
	ctimes = ntimes;
}

/**
 * Sorts nodes and deletes unused nodes.
 * This can save a lot of memory, since AnalysisSDF::crop and AnalysisSDF::filter only delete segments, not unused nodes
 */
void SegmentAnalyser::pack() {
	std::vector<double> ni(nodes.size());
	std::fill(ni.begin(),ni.end(), 0.);
	std::vector<Vector3d> newnodes;
	for (auto& s:segments) {
		if (ni.at(s.x)==0) { // the node is new
			newnodes.push_back(nodes.at(s.x));
			ni.at(s.x) = newnodes.size()-1; // set index of the new node
		}
		s.x = ni.at(s.x);
		if (ni.at(s.y)==0) { // the node is new
			newnodes.push_back(nodes.at(s.y));
			ni.at(s.y) = newnodes.size()-1; // set index of the new node
		}
		s.y = ni.at(s.y);
	}
	// std::cout << "pack(): nodes: " << nodes.size() << " -> " << newnodes.size() << ", " << double(newnodes.size())/double(nodes.size()) << " \n";
	nodes = newnodes; // kabum!
}

/**
 *  Numerically computes the intersection point
 *
 * @param in       the node within the domain
 * @param out      the node outside of the domain
 * @param geometry signed distance function of the geometry
 * \return         the intersection point
 */
Vector3d SegmentAnalyser::cut(Vector3d in, Vector3d out, SignedDistanceFunction* geometry)
{
	assert(geometry->getDist(in)<=0);
	assert(geometry->getDist(out)>=0);
	if (std::abs(geometry->getDist(out))>1e-6) {
		Vector3d c =  in.plus(out).times(0.5); // mid
		if (geometry->getDist(c)<0) { // in
			return cut(c,out,geometry);
		} else { // out
			return cut(in,c,geometry);
		}
	} else {
		return out;
	}
}

/**
 * \return The summed parameter of type @param st (@see RootSystem::ScalarType)
 */
double SegmentAnalyser::getSummed(std::string pname) const {
	std::vector<double> v_ = getScalar(pname);
	return std::accumulate(v_.begin(), v_.end(), 0.0);
}

/**
 * \return The summed parameter of type @param st (@see RootSystem::ScalarType),
 * that is within geometry @param g based on the segment mid point (i.e. not exact).
 * To sum exactly, first crop to the geometry, then run SegmentAnalyser::getSummed(st).
 */
double SegmentAnalyser::getSummed(std::string pname, SignedDistanceFunction* g) const {
	std::vector<double> data = getScalar(pname);
	double v = 0;
	for (size_t i=0; i<segments.size(); i++) {
		double d = data.at(i);
		Vector2i s = segments.at(i);
		Vector3d n1 = nodes.at(s.x);
		Vector3d n2 = nodes.at(s.y);
		Vector3d mid = n1.plus(n2).times(0.5);
		if (g->getDist(mid)<0) {
			v += d;
		}
	}
	return v;
}



/**
 * Return the unique origins of the segments
 */
std::vector<Organ*> SegmentAnalyser::getRoots() const
{
	std::set<Organ*> rootset;  // praise the stl
	for (Organ* r : segO) {
		rootset.insert(r);
	}
	return std::vector<Organ*>(rootset.begin(), rootset.end());
}

/**
 * \return The number of roots
 */
int SegmentAnalyser::getNumberOfRoots() const
{
	const auto& rootset = getRoots();
	return rootset.size();
}



/**
 * Projects the segments to an image plane (todo verify this code)
 *
 * @param pos       position of camera
 * @param ons       orthonormal system, row 1 is orthogonal to the image plane given by [row 2,row 3]
 * @param fl        focal length, alpha = 2*arctan(d/(2*fl)), were alpha is the angle of field, and d the image diagonal
 *
 * \return The image segments in the x-y plane (z=0)
 */
SegmentAnalyser SegmentAnalyser::foto(const Vector3d& pos, const Matrix3d& ons, double fl) const
{
	SegmentAnalyser f(*this); // copy
	for (auto& n : f.nodes) { // translate
		n = n.minus(pos);
	}
	Matrix3d m = ons.inverse(); // rotate
	for (auto& n : f.nodes) {
		n = m.times(n);
	}
	//	// crop to objects in front of the camera
	Vector3d o(0.,0.,0.);
	Vector3d plane(0.,0., -1);
	SDF_HalfPlane sdf = SDF_HalfPlane(o,plane);
	f.crop(&sdf);
	f.pack();
	// project
	for (auto& a : f.nodes) {
		a = a.times(fl/(-plane.times(a)));
		a.z = 0;
	}
	// final image crop
	SDF_PlantBox box(1.,1.,2.); // TODO --> d = sqrt(2) cm
	f.crop(&box);
	f.pack();
	return f;
}

/**
 * Keeps the segments that intersect with a plane
 *
 * @param plane 	half plane
 */
SegmentAnalyser SegmentAnalyser::cut(const SDF_HalfPlane& plane) const
{
	SegmentAnalyser f;
	f.nodes = nodes; // copy all nodes
	for (size_t i=0; i<segments.size(); i++) {
		Vector2i s = segments.at(i);
		Vector3d n1 = nodes.at(s.x);
		Vector3d n2 = nodes.at(s.y);
		double d = plane.getDist(n1)*plane.getDist(n2);
		if (d<=0) { // one is inside, one is outside
			f.segments.push_back(s);
			f.ctimes.push_back(ctimes.at(i));
			f.segO.push_back(segO.at(i));
		}
	}
	f.pack(); // delete unused nodes
	return f;
}

/**
 *  Creates a vertical distribution of the parameter of type @param st (@see RootSystem::ScalarType)
 *
 * @param st        parameter type @see RootSystem::ScalarType
 * @param top       vertical top position (cm)
 * @param bot       vertical bot position (cm)
 * @param n         number of layers (each with a height of (bot-top)/n )
 * @param exact     calculates the intersection with the layer boundaries (true), only based on segment midpoints (false)
 * \return Vector of size @param n containing the summed parameter in this layer
 */
std::vector<double> SegmentAnalyser::distribution(std::string pname, double top, double bot, int n, bool exact) const
{
	std::vector<double> d(n);
	double dz = (bot-top)/double(n);
	SDF_PlantBox* layer = new SDF_PlantBox(1e100,1e100,dz);
	for (int i=0; i<n; i++) {
		Vector3d t(0,0,top-i*dz);
		SDF_RotateTranslate g(layer,t);
		if (exact) {
			SegmentAnalyser a(*this); // copy everything
			a.crop(&g); // crop exactly
			d.at(i) = a.getSummed(pname);
		} else {
			d.at(i) = this->getSummed(pname, &g);
		}
	}
	delete layer;
	return d;
}

/**
 *  Creates a vertical distribution
 *
 * @param top       vertical top position (cm)
 * @param bot       vertical bot position (cm)
 * @param n         number of layers (each with a height of (bot-top)/n )
 * \return Vector of size @param n containing an Analysis object of the layers (cropped exactly)
 */
std::vector<SegmentAnalyser> SegmentAnalyser::distribution(double top, double bot, int n) const
{
	std::vector<SegmentAnalyser> d(n);
	double dz = (bot-top)/double(n);
	SDF_PlantBox* layer = new SDF_PlantBox(1e100,1e100,dz);
	for (int i=0; i<n; i++) {
		Vector3d t(0,0,top-i*dz);
		SDF_RotateTranslate g(layer,t);
		SegmentAnalyser a = SegmentAnalyser(*this); // copy everything
		a.crop(&g); // crop exactly
		d.at(i) = a;
	}
	delete layer;
	return d;
}

/**
 *  Creates a two-dimensional distribution of the parameter of type @param st (@see RootSystem::ScalarType)
 *
 * @param st        parameter type @see RootSystem::ScalarType
 * @param top       vertical top position (cm)
 * @param bot       vertical bot position (cm)
 * @param left      left along x-axis (cm)
 * @param right     right along x-axis (cm)
 * @param n         number of vertical grid elements (each with height of (bot-top)/n )
 * @param m 		number of horizontal grid elements (each with length of (right-left)/m)
 * @param exact     calculates the intersection with the layer boundaries (true), only based on segment midpoints (false)
 * \return Vector of size @param n containing the summed parameter in this layer
 */
std::vector<std::vector<double>> SegmentAnalyser::distribution2(std::string pname, double top, double bot, double left, double right, int n, int m, bool exact) const
{
	std::vector<std::vector<double>> d(n);
	double dz = (bot-top)/double(n);
	double dx = (right-left)/double(m);
	SDF_PlantBox* layer = new SDF_PlantBox(dx,1e9,dz);

	for (int i=0; i<n; i++) {

		std::vector<double> row(m); // m columns
		for (int j=0; j<m; j++) {

			Vector3d t(left+(j+0.5)*dx,0.,top-i*dz); // box is [-x/2,-y/2,0] - [x/2,y/2,-z]
			SDF_RotateTranslate g(layer,t);

			if (exact) {
				SegmentAnalyser a(*this); // copy everything
				a.crop(&g); // crop exactly
				row.at(j) = a.getSummed(pname);
			} else {
				row.at(j) = this->getSummed(pname, &g);
			}

		}

		d.at(i)=row; // store the row (n rows)
	}
	delete layer;
	return d;
}

/**
 *  Creates a vertical distribution
 *
 * @param top       vertical top position (cm)
 * @param bot       vertical bot position (cm)
 * @param left      left along x-axis (cm)
 * @param right     right along x-axis (cm)
 * @param n         number of vertical grid elements (each with height of (bot-top)/n )
 * @param m 		number of horizontal grid elements (each with length of (right-left)/m)
 * \return Vector of size @param n containing the summed parameter in this layer
 */
std::vector<std::vector<SegmentAnalyser>> SegmentAnalyser::distribution2(double top, double bot, double left, double right, int n, int m) const
{
	std::vector<std::vector<SegmentAnalyser>> d(n);
	double dz = (bot-top)/double(n);
	double dx = (right-left)/double(m);
	SDF_PlantBox* layer = new SDF_PlantBox(dx,1e4,dz);
	// std::cout << "dx " << dx  <<", dz "<< dz << "\n";
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			Vector3d t(left+(j+0.5)*dx,0.,top-i*dz); // box is [-x/2,-y/2,0] - [x/2,y/2,-z]
			SDF_RotateTranslate g(layer,t);
			SegmentAnalyser a(*this); // copy everything
			a.crop(&g); // crop exactly
			d.at(i).push_back(a);

			if (i==0) {
				std::stringstream ss;
				ss << "(" << i << ":" << j << ").py";
				std::string name = "test"+ss.str();
				std::ofstream fos;
				fos.open(name.c_str());
				fos << "from paraview.simple import *\n";
				fos << "paraview.simple._DisableFirstRenderCameraReset()\n";
				fos << "renderView1 = GetActiveViewOrCreate('RenderView')\n\n";
				g.writePVPScript(fos);
				fos.close();
			}
		}
	}
	delete layer;
	return d;
}

/**
 * Exports the simulation results with the type from the extension in name
 * (that must be lower case)
 *
 * @param name      file name e.g. output.vtp
 */
void SegmentAnalyser::write(std::string name)
{
	this->pack(); // a good idea before writing any file
	std::ofstream fos;
	fos.open(name.c_str());
	std::string ext = name.substr(name.size()-3,name.size()); // pick the right writer
	if (ext.compare("vtp")==0) {
		std::cout << "writing VTP: " << name << "\n";
		this->writeVTP(fos,{ /* todo */ }); // depricated Plant::st_radius, Plant::st_subtype, Plant::st_time, Plant::st_otype
	} else if (ext.compare("txt")==0)  {
		std::cout << "writing text file for Matlab import: "<< name << "\n";
		this->writeRBSegments(fos);
	} else if (ext.compare("dgf")==0)  {
		std::cout << "writing dgf file: "<< name << "\n";
		this->writeDGF(fos);
	} else {
		throw std::invalid_argument("SegmentAnalyser::write: Unknown file type");
	}
	fos.close();
}

/**
 * Writes a VTP file with @param types data per segment.
 *
 * @param os        typically a file out stream
 * @param types     multiple parameter types (@see RootSystem::ScalarType) that are saved in the VTP file,
 * 					additionally, all userdata is saved per default
 */
void SegmentAnalyser::writeVTP(std::ostream & os, std::vector<std::string> typeNames) const
{
	assert(segments.size() == segO.size());
	assert(segments.size() == ctimes.size());
	os << "<?xml version=\"1.0\"?>";
	os << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	os << "<PolyData>\n";
	os << "<Piece NumberOfLines=\""<< segments.size() << "\" NumberOfPoints=\""<< nodes.size()<< "\">\n";
	// data (CellData)
	os << "<CellData Scalars=\" CellData\">\n";
	for (auto n : typeNames) {
		std::vector<double> data = getScalar(n);
		os << "<DataArray type=\"Float32\" Name=\"" << n << "\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
		for (auto const& t : data) {
			os << t << " ";
		}
		os << "\n</DataArray>\n";
	}
	// write user data
	for (size_t i=0; i<userData.size(); i++) {
		const auto& data = userData.at(i);
		std::string name = userDataNames.at(i);
		os << "<DataArray type=\"Float32\" Name=\"" << name << "\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
		for (auto const& t : data) {
			os << t << " ";
		}
		os << "\n</DataArray>\n";
	}
	os << "\n</CellData>\n";
	// nodes (Points)
	os << "<Points>\n"<<"<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\" >\n";
	for (auto const& n:nodes) {
		os << n.x << " "<< n.y <<" "<< n.z<< " ";
	}
	os << "\n</DataArray>\n"<< "</Points>\n";
	// segments (Lines)
	os << "<Lines>\n"<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	for (auto const& s:segments) {
		os << s.x << " " << s.y << " ";
	}
	os << "\n</DataArray>\n"<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	for (size_t i=0; i<segments.size(); i++) {
		os << 2*i+2 << " ";
	}
	os << "\n</DataArray>\n";
	os << "\n</Lines>\n";
	//
	os << "</Piece>\n";
	os << "</PolyData>\n" << "</VTKFile>\n";
}

/**
 * Writes the (line)segments of the root system, and
 * mimics the Matlab script getSegments() of RootBox
 *
 * @param os      typically a file out stream
 */
void SegmentAnalyser::writeRBSegments(std::ostream & os) const
{
         double radius ;
                double red ;
                double green;
                double blue ;
                double time ;
                double type ;
                double organ;

	os << "x1 y1 z1 x2 y2 z2 radius R G B time type organ \n";
	for (size_t i=0; i<segments.size(); i++) {
		Vector2i s = segments.at(i);
		Vector3d n1 = nodes.at(s.x);
		Vector3d n2 = nodes.at(s.y);
		Organ* r = segO.at(i);
		organ=r->organType();
		// std::cout << "#" << i <<": organ " << r->organType() << " (" << s.x<< ", "<< s.y << ")  \n";
		switch (r->organType()){

            case Organ::ot_root :
                 radius = ((Root*)r)->tParam()->a;
                red = ((Root*)r)->tParam()->colorR;
                green = ((Root*)r)->tParam()->colorG;
                blue = ((Root*)r)->tParam()->colorB;
                time = ctimes.at(i);
                type = r->param->subType;

                os << std::fixed << std::setprecision(4)<< n1.x << " " << n1.y << " " << n1.z << " " << n2.x << " " << n2.y << " " << n2.z << " " <<
                        radius << " " << red << " " << green << " " << blue << " " << time<< " " << type << " " << organ <<" \n";
break;
            case Organ::ot_stem :
                radius = ((Stem*)r)->stParam()->a;
                red = ((Stem*)r)->stParam()->colorR;
                green = ((Stem*)r)->stParam()->colorG;
                blue = ((Stem*)r)->stParam()->colorB;
                time = ctimes.at(i);
                type = r->param->subType;
                os << std::fixed << std::setprecision(4)<< n1.x << " " << n1.y << " " << n1.z << " " << n2.x << " " << n2.y << " " << n2.z << " " <<
                        radius << " " << red << " " << green << " " << blue << " " << time<< " " << type << " " << organ << " \n";
                        break;
                        //
                        case Organ::ot_leafe :
                radius = ((Leaf*)r)->ltParam()->a;
                red = ((Leaf*)r)->ltParam()->colorR;
                green = ((Leaf*)r)->ltParam()->colorG;
                blue = ((Leaf*)r)->ltParam()->colorB;
                time = ctimes.at(i);
                type = r->param->subType;
                os << std::fixed << std::setprecision(4)<< n1.x << " " << n1.y << " " << n1.z << " " << n2.x << " " << n2.y << " " << n2.z << " " <<
                        radius << " " << red << " " << green << " " << blue << " " << time<< " " << type << " " << organ << " \n";
                        break;
		}
	}
	//std::cout << "fin\n";
}

/**
 * Writes the (line)segments of the root system in dgf format used by DuMux
 *
 * @param os      typically a file out stream
 */

void SegmentAnalyser::writeDGF(std::ostream & os) const
{
	os << "DGF \n";
	os << "Vertex \n";
	for (auto& n : nodes) {
		os << n.x/100 << " " << n.y/100 << " " << n.z/100 << " \n";
	}

	os << "# \n";
	os << "SIMPLEX \n";
	os << "parameters 10 \n";
	// node1ID, node2ID, type, branchID, surfaceIdx, length, radiusIdx, massIdx, axialPermIdx, radialPermIdx, creationTimeId
	for (size_t i=0; i<segments.size(); i++) {
		Vector2i s = segments.at(i);
		Vector3d n1 = nodes.at(s.x);
		Vector3d n2 = nodes.at(s.y);
		Organ* r = segO.at(i);
		int branchnumber = r->id;
		double radius = 0; // r->param.a; // TODO
		double length = sqrt((n1.x-n2.x)*(n1.x-n2.x)+(n1.y-n2.y)*(n1.y-n2.y)+(n1.z-n2.z)*(n1.z-n2.z));
		double surface = 2*radius*M_PI*length;
		double time = ctimes.at(i);
		double type = 0; // r->param.type; // TODO
		os << s.x << " " << s.y << " " << type << " " << branchnumber << " " << surface/10000 << " " << length/100 <<" " << radius/100 << " " << "0.00" << " " << "0.0001" << " "<< "0.00001" << " " << time*3600*24 << " \n";
	}

	os << "# \n";
	os << "BOUNDARYDOMAIN \n";
	os << "default 1 \n";
	os << "# \n";
}

