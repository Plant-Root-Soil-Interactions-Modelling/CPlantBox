#include "Plant.h"
#include <memory>
#include <iostream>

namespace CPlantBox {


unsigned int Plant::noParamFile[5] = {0, 0, 1, 1, 1}; // check if there are parameter files TODO to make it simpler used in Plant::readParameter

Plant::Plant()
{
	initOTP();
	setParameter(new SeedTypeParameter());
	seed = new Seed(this);
}

Plant::~Plant()
{
	for (auto& otp_:organParam) {
		for (auto& otp:otp_) {
			delete otp;
		}
	}
	delete seed;
	delete geometry;
}

/**
 *	Deletes the old parameter, sets the new one
 */
void Plant::setParameter(OrganTypeParameter*  otp)
{
	unsigned int ot = otp->organType;
	unsigned int i = ot2index(ot);
	delete organParam.at(i).at(otp->subType);
	organParam.at(i).at(otp->subType) = otp;
}

OrganTypeParameter* Plant::getParameter(int otype, int subtype) const
{
	return organParam.at(ot2index(otype)).at(subtype);
}

unsigned int Plant::ot2index(unsigned int ot) {
	switch (ot) { //check the type of the organ
	case Organ::ot_seed: return 0;
	case Organ::ot_root: return 1;
	case Organ::ot_stem: return 2;
	case Organ::ot_leafe: return 3;
	default:
		throw std::invalid_argument("Plant::setOrganTypeParameter: pure organ type expected");
	}
}

/**
 * Deletes the old geometry, sets the new one, passes it to the tropisms
 */
void Plant::setGeometry(SignedDistanceFunction* geom)
{
	delete geometry;
	geometry = geom;
	for (int i=0; i<maxtypes; i++) {
		RootTypeParameter* rtp = (RootTypeParameter*) getParameter(Organ::ot_root,i);
		if (rtp->subType!=-1) { // defined
			delete rtp->tropism;
			rtp->createTropism(geom);
		}
	}
}

/**
 * Resets the root system: deletes all roots, does not change the parameters
 */
void Plant::reset()
{

	delete seed; // TODO??????
	seed = new Seed(this); // make_shared<Seed>(this), seed of type shared_ptr<Seed>
	simtime=0;
	rid = -1;
	nid = -1;
}

/**
 * Puts default values into the root type parameters vector
 */
void Plant::initOTP()
{
	organParam = std::vector<std::vector<OrganTypeParameter*> >(maxorgans);
	for (auto& otp:organParam) {
		otp = std::vector<OrganTypeParameter*>(maxtypes);
		for (size_t i=0; i<otp.size(); i++) {
			otp.at(i) = new OrganTypeParameter();
		}
	}
}

/**
 * Reads the root parameter from a file. Opens plant parameters with the same filename if available,
 * othterwise assumes a tap root system at position (0,0,-3).
 *
 * @param name          filename without file extension
 * @param subdir        directory ("modelparameter/" by default)
 */
void Plant::openXML(std::string name, std::string subdir) //The first run will convert the rparam file to the XML, and then use the XML later on.
{

	std::ifstream fis;
	std::string XMLname = subdir;
	XMLname.append(name);
	XMLname.append(".xml");
	fis.open(XMLname.c_str());
	if (fis.good()){

		tinyxml2::XMLDocument xmlParamFile;
		xmlParamFile.LoadFile(XMLname.c_str());
		int c = 0;
		for (const tinyxml2::XMLElement* organ_param = xmlParamFile.FirstChildElement( "Plant" )->FirstChildElement( "organ" ); organ_param != 0 ; organ_param = organ_param->NextSiblingElement("organ") )
		{
			if (organ_param->Attribute("type", "seed")) {
				SeedTypeParameter* stp = (SeedTypeParameter*)getParameter(Organ::ot_seed,0);
				stp->readXML(organ_param);
				c++;
//				std::cout << " Read from XML " << c << " seed type parameters \n";
			}


			if (organ_param->Attribute("type", "root")) {
				RootTypeParameter* p  = new RootTypeParameter();
				p->readXML(organ_param);
				setParameter(p);
				//                root_element = root_element->NextSiblingElement("organ") ;
				c++;
				Plant::noParamFile[1] = 0;
				std::cout << " Read from XML " << c << " root type parameters \n";
			}

			if (organ_param->Attribute("type", "stem")) {

				StemTypeParameter* stem_p  = new StemTypeParameter();
				stem_p->readXML(organ_param);
				setParameter(stem_p);
				c++;
//				std::cout << " Read from XML " << c << " stem type parameters \n";
				Plant::noParamFile[2] = 0;
			} else //{Plant::noParamFile[2] = 1;}
			if (organ_param->Attribute("type", "leaf")) {

				LeafTypeParameter* leaf_p  = new LeafTypeParameter();
				leaf_p->readXML(organ_param);
				setParameter(leaf_p);
				c++;
//				std::cout << " Read from XML " << c << " leaf type parameters \n";
				Plant::noParamFile[3] = 0;
			}
		}

		//        xmlParamFile.FirstChildElement( "Plant" )->FirstChildElement( "organ" )->FirstChildElement( "parameter" )->QueryDoubleAttribute("location_z", &seedz);

		//            if (root_p->Attribute("type" , "root" ) && root_p->Attribute("name", "taproot") )
		//            {
		//              const  tinyxml2::XMLElement* taproot_p  = root_p->FirstChildElement("parameter");
		//            }
		//            const char* aa = "xx";
		//            double cc = 0;
		//            double lbs = 0;
		//           cc = xmlParamFile.FirstChildElement( "Plant" )->FirstChildElement( "organ" )->NextSiblingElement("organ")->FirstChildElement("parameter")->DoubleAttribute("value");
		//            taproot_p->QueryStringAttribute("name", &aa);
		//

		writeAlltoXML(name.append("_new"));

	} else{ openFile(name, subdir);
	writeAlltoXML(name);
	}
}




void Plant::openFile(std::string name, std::string subdir)
{
	std::ifstream fis;
	std::string rp_name = subdir;
	rp_name.append(name);
	rp_name.append(".rparam");
	fis.open(rp_name.c_str());
	int c = 0;
	if (fis.good()) { // did it work?
		c = readRootParameters(fis);
		fis.close();
		std::cout << "Read " << c << " root type parameters \n";
	} else {
		std::cout << "No root system parameter file found \n";
		Plant::noParamFile [1] = 1;

	}
	// debug

	// open seed parameter
	SeedTypeParameter* stp = (SeedTypeParameter*)getParameter(Organ::ot_seed,0);
	std::string pp_name = subdir;
	pp_name.append(name);
	pp_name.append(".pparam");
	fis.open(pp_name.c_str());
	if (fis.good()) { // did it work?
		stp->read(fis);
		fis.close();
	} else { // create a tap root system
		std::cout << "No seed system parameter file found, using default tap root system \n";
		delete stp;
		setParameter(new SeedTypeParameter());
	}


	// open stem parameter

	std::string stp_name = subdir;
	stp_name.append(name);
	stp_name.append(".stparam");
	fis.open(stp_name.c_str());
	int stem_c = 0;
	if (fis.good()) { // did it work?
		stem_c = readStemParameters(fis);
		fis.close();
	} else {
		std::cout << "No stem parameter file found, skipping leaf  \n";
		Plant::noParamFile [2] = 1;
	}
	std::cout << "Read " << stem_c << " stem type parameters \n"; // debug

	std::string lp_name = subdir;
	lp_name.append(name);
	lp_name.append(".leparam");
	fis.open(lp_name.c_str());
	int leaf_c = 0;
	if (fis.good()) { // did it work?
		leaf_c = readLeafParameters(fis);
		fis.close();
	} else {
		std::cout << "No leaf parameters found, skipping leaf  \n";
		Plant::noParamFile [3] = 1;


	}
	std::cout << "Read " << leaf_c << " leaf type parameters \n"; // debug
    writeAlltoXML(name);
}

/**
 * Reads parameter from input stream (there is a Matlab script exporting these, @see writeParams.m) TODO
 *
 * @param cin  in stream
 */
int Plant::readRootParameters(std::istream& cin)
{
	// initOTP();
	int c = 0;
	while (cin.good()) {
		RootTypeParameter* p  = new RootTypeParameter();
		p->read(cin);
		setParameter(p); // sets the param to the index (p.type-1) TODO
		c++;
	}
	return c;
}

int Plant::readStemParameters(std::istream& cin)
{
	// initOTP();
	int stem_c = 0;
	while (cin.good()) {

		StemTypeParameter* stem_p  = new StemTypeParameter();///added copypaste
		stem_p->read(cin);

		//		setOrganTypeParameter(p); // sets the param to the index (p.type-1) TODO
		setParameter(stem_p);

		stem_c++;
	}
	return stem_c;
}


int Plant::readLeafParameters(std::istream& cin)
{
	// initOTP();
	int leaf_c = 0;
	while (cin.good()) {
		LeafTypeParameter* leaf_p  = new LeafTypeParameter();///added copypaste
		leaf_p->read(cin);
		setParameter(leaf_p);
		leaf_c++;
	}
	return leaf_c;

}
/**
 * Writes root parameters (for debugging) TODO
 *
 * @param os  out stream
 */
void Plant::writeParameters(std::ostream& os) const
{
	const char* declaration ="<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>";
	tinyxml2::XMLDocument doc;
	tinyxml2::XMLDocument dxml;

	for (auto const& otp_:organParam) {
		unsigned int t = 0;
		for (auto const& otp : otp_) {
			if ((otp->organType>=0) && (otp->subType>=0) && ((otp->subType)==t)) {
								assert((otp->subType)==t); // check if index really equals subType-1
				os << otp->writeXML(0);
			}
			t++;
		}
	}



}

void Plant::writeAlltoXML(std::string name, std::string subdir){
	std::string xmlname = subdir;
	xmlname.append(name);
	xmlname.append(".xml");
	std::ofstream xmloutput;
	xmloutput.open( xmlname.c_str());
	xmloutput<<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
	xmloutput<<"\n<Plant name=\""<<name<<"\" filetype=\"parameters\">\n";
	for (auto const& otp_:organParam) {
		unsigned int t = 0;
		for (auto const& otp : otp_) {
			if ((otp->organType>=1) && (otp->subType>=0) && ((otp->subType)==t) ) {
								assert((otp->subType)==t); // check if index really equals subType-1
				xmloutput<<otp->writeXML(0);
				xmloutput<<std::endl;
			} else {
			}

			t++;

		}
	}
	xmloutput<<"\n</Plant>\n";
	xmloutput.close();
}
/**
 * Sets up the plant according to the given parameters
 */
void Plant::initialize()
{

	reset(); // deletes the old seed, makes a new one
	seed->initialize(seed->initializeparam());// randomness in the python binding

	/* the following code will be moved to the shoot */
	//  // Shoot borne roots
	//  if ((rs.nC>0) && (rs.delaySB<maxT)) { // if the type is not defined, copy basal root
	//      //		if (getRootTypeParameter(shootbornetype)->type<1) {
	//      //			std::cout << "Shootborne root type #" << shootbornetype << " was not defined, using tap root parameters instead\n";
	//      //			RootTypeParameter srtp = RootTypeParameter(*getRootTypeParameter(1));
	//      //			srtp.type = shootbornetype;
	//      //			setRootTypeParameter(srtp);
	//      //		}
	//      Vector3d sbpos = rs.seedPos;
	//      sbpos.z=sbpos.z/2.; // half way up the mesocotyl
	//      int maxSB = ceil((maxT-rs.firstSB)/rs.delayRC); // maximal number of root crowns
	//      double delay = rs.firstSB;
	//      for (int i=0; i<maxSB; i++) {
	//          //			for (int j=0; j<rs.nC; j++) {
	//          //				Organ* shootborne = new Organ(this, shootbornetype, iheading ,delay, nullptr, 0, 0);
	//          //				// TODO fix the initial radial heading
	//          //				shootborne->addNode(sbpos,delay);
	//          //				baseRoots.push_back(shootborne);
	//          //				delay += rs.delaySB;
	//          //			}
	//          sbpos.z+=rs.nz;  // move up, for next root crown
	//          delay = rs.firstSB + i*rs.delayRC; // reset age
	//      }
	//   }

}

/**
 * Simulates root system growth for time span dt
 *
 * @param dt    	time step [days]
 * @param silence 	indicates if status is written to the console (cout) (default = false)
 */
void Plant::simulate(double dt, bool silence)
{
	if (!silence) {
//		std::cout << "Plant.simulate(dt) from "<< simtime << " to " << simtime+dt << " days \n";
	}
	old_non = getNumberOfNodes();
	old_nor = organs.size();
	seed->simulate(dt);
	simtime+=dt;
	organs.clear(); // empty buffer
	organs_type = -1;
}

/**
 * Simulates root system growth for the time span defined in the parameter set TODO
 */
void Plant::simulate()
{
	// this->simulate(seedParam->simtime); TODO
}


int Plant::getNumberOfSegments() const
{
	getOrgans(Organ::ot_organ);
	int c = 0;
	for (const auto& o : organs) {
		c += (o->r_nodes.size()-1);
	}
	return c;
}

/**
 * Returns a reference to the sequential list of organs,
 * and updates the list if necessary.
 *
 * @param otype 	organ type
 */
std::vector<Organ*>& Plant::getOrgans(unsigned int otype) const
{
	if (organs_type!=otype) { // create buffer
		organs = seed->getOrgans(otype);
		organs_type = otype;
		return organs;
	} else { // return buffer
		return organs;
	}
}

/**
 * Copies the nodes of the root systems into a sequential vector,
 * nodes are unique (default). See also RootSystem::getSegments
 */
std::vector<Vector3d> Plant::getNodes() const
{
	this->getOrgans(Organ::ot_organ); // update roots (if necessary)
	int non = getNumberOfNodes();
	std::vector<Vector3d> nv = std::vector<Vector3d>(non); // reserve big enough vector
	for (auto const& r: organs) {
		for (size_t i=0; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
			nv.at(r->getNodeID(i)) = r->getNode(i); // pray that ids are correct
		}
	}
	return nv;
}

/**
 * Returns the root system as polylines, i.e. each root is represented by its nodes
 */
std::vector<std::vector<Vector3d>> Plant::getPolylines(unsigned int otype) const
{
	this->getOrgans(otype); // update roots (if necessary)
	std::vector<std::vector<Vector3d>> nodes = std::vector<std::vector<Vector3d>>(organs.size()); // reserve big enough vector
	for (size_t j=0; j<organs.size(); j++) {
		std::vector<Vector3d>  rn = std::vector<Vector3d>(organs[j]->getNumberOfNodes());
		for (size_t i=0; i<organs[j]->getNumberOfNodes(); i++) { // loop over all nodes of all roots
			rn.at(i) = organs[j]->getNode(i);
		}
		nodes[j] = rn;
	}
	return nodes;
}

/**
 * Return the segments of the root system at the current simulation time
 */
std::vector<Vector2i> Plant::getSegments(unsigned int otype) const
{
	this->getOrgans(otype); // update roots (if necessary)
	int nos=getNumberOfSegments();
	std::vector<Vector2i> s(nos);
	int c=0;
	for (auto const& r:organs) {
		for (size_t i=0; i<r->getNumberOfNodes()-1; i++) {
			Vector2i v(r->getNodeID(i),r->getNodeID(i+1));
			s.at(c) = v;
			c++;
		}
	}
	return s;
}


std::vector<int> Plant::getNodesOrganType() const
{
	this->getOrgans(15);
	int nos = getNumberOfSegments();
	std::vector<int> s(nos);
	int c = 0;
	for (auto const& r : organs) {
		for (size_t i = 0; i<r->getNumberOfNodes() - 1; i++) {
			int v(r->getScalar("organtype"));
			s.at(c) = v;
			c++;
		}
	}
	return s;
}
/**
 * Returns pointers to the organs corresponding to each segment
 */
std::vector<Organ*> Plant::getSegmentsOrigin(unsigned int otype) const
{
	this->getOrgans(otype); // update (if necessary)
	int nos=getNumberOfSegments();
	std::vector<Organ*> s(nos);
	int c=0;
	for (auto const& o:organs) {
		for (size_t i=0; i<o->getNumberOfNodes()-1; i++) {
			s.at(c) = o;
			c++;
		}
	}
	return s;
}

/**
 * Copies the node emergence times of the root systems into a sequential vector,
 * see RootSystem::getNodes()
 */
std::vector<double> Plant::getNETimes() const
{
	this->getOrgans(Organ::ot_organ); // update roots (if necessary)
	int nos=getNumberOfSegments();
	std::vector<double> netv = std::vector<double>(nos); // reserve big enough vector
	int c=0;
	for (auto const& r: organs) {
		for (size_t i=1; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
			netv.at(c) = r->getNodeCT(i); // pray that ids are correct
			c++;
		}
	}
	return netv;
}

/**
 *  Returns the node emergence times to the corresponding polylines, see also RootSystem::getPolylines
 */
std::vector<std::vector<double>> Plant::getPolylinesNET(unsigned int otype) const
{
	this->getOrgans(otype); // update roots (if necessary)
	std::vector<std::vector<double>> times = std::vector<std::vector<double>>(organs.size()); // reserve big enough vector
	for (size_t j=0; j<organs.size(); j++) {
		std::vector<double>  rt = std::vector<double>(organs[j]->getNumberOfNodes());
		for (size_t i=0; i<organs[j]->getNumberOfNodes(); i++) {
			rt[i] = organs[j]->getNodeCT(i);
		}
		times[j] = rt;
	}
	return times;
}

/**
 * Copies a scalar that is constant per organ to a sequential vector (one scalar per organ).
 *
 * @param stype     a scalar type (@see RootSystem::ScalarTypes). st_time is the emergence time of the root
 */
std::vector<double> Plant::getScalar(unsigned int otype, std::string name) const
{
	this->getOrgans(otype); // update roots (if necessary)
	std::vector<double> scalars(organs.size());
	for (size_t i=0; i<organs.size(); i++) {
		scalars.at(i) = organs[i]->Organ::getScalar(name);
	}
	return scalars;
}



/**
 * todo most important debug informations
 */
std::string Plant::toString() const
{
	return "todo";
}

/**
 * Exports the simulation results with the type from the extension in name
 * (that must be lower case)
 *
 * @param name      file name e.g. output.vtp
 */
void Plant::write(std::string name, int otype) const
{
	std::ofstream fos;
	fos.open(name.c_str());
	std::string ext = name.substr(name.size()-3,name.size()); // pick the right writer
	if (ext.compare("sml")==0) {
		std::cout << "writing RSML... "<< name.c_str() <<"\n";
		writeRSML(fos);
	} else if (ext.compare("vtp")==0) {
		std::cout << "writing VTP... "<< name.c_str() <<"\n";
		TiXMLwriteVTP(otype, fos);   ///< Write .VTP file by using TinyXML2 performance slowed by 0.5 seconds but precision increased.
	} else if (ext.compare(".py")==0)  {
		std::cout << "writing Geometry ... "<< name.c_str() <<"\n";
		writeGeometry(fos);
	} else {
		throw std::invalid_argument("RootSystem::write(): Unkwown file type");
	}
	fos.close();
}

/**
 * Creates an RSML file
 *
 * @param os      typically a file out stream
 */
void Plant::writeRSML(std::ostream & os) const
{
	os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"; // i am not using utf-8, but not sure if ISO-8859-1 is correct
	os << "<rsml>\n";
	writeRSMLMeta(os);
	os<< "<scene>\n";
	writeRSMLPlant(os);
	os << "</scene>\n";
	os << "</rsml>\n";
}

/**
 * Writes RML meta data tag
 *
 * @param os      typically a file out stream
 */
void Plant::writeRSMLMeta(std::ostream & os) const
{
	os << "<metadata>\n";
	os << "\t<version>" << 1 << "</version>\n";
	os << "\t<unit>" << "cm" << "</unit>\n";
	os << "\t<resolution>" << 1 << "</resolution>\n";
	// fetch time
	//    os << "<last-modified>";
	//    auto t = std::time(nullptr);
	//    auto tm = *std::localtime(&t);
	//    os << std::put_time(&tm, "%d-%m-%Y"); // %H-%M-%S" would do the job for gcc 5.0
	//    os << "</last-modified>\n";
	os << "\t<software>CPlantBox</software>\n";
	os << "</metadata>\n";
}

/**
 * Writes RSML plant tag
 *
 * @param os      typically a file out stream
 */
void Plant::writeRSMLPlant(std::ostream & os) const
{
	os << "<plant>\n";
	seed->writeRSML(os,"");
	os << "</plant>\n";
}

/**
 * Writes current simulation results as VTP (VTK polydata file),
 * where each root is represented by a polyline.
 *
 * Use SegmentAnalyser::writeVTP() for a representation based on segments,
 * e.g. for creating a movie (and run the animate.py script), or mapping values to segments
 *
 * @param os      typically a file out stream
 */
void Plant::writeVTP(int otype, std::ostream & os) const // currently not using, if TiXMLwriteVtP has no bug, than we could delete this member function
{
	this->getOrgans(otype); // update roots (if necessary)
	const auto& nodes = getPolylines(otype);
	const auto& times = getPolylinesNET(otype);

	os << "<?xml version=\"1.0\"?>";
	os << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	os << "<PolyData>\n";
	int non = 0; // number of nodes
	for (auto const& r : organs) {
		non += r->getNumberOfNodes();
	}
	int nol=organs.size(); // number of lines
	os << "<Piece NumberOfLines=\""<< nol << "\" NumberOfPoints=\""<<non<<"\">\n";

	// POINTDATA
	os << "<PointData Scalars=\" PointData\">\n" << "<DataArray type=\"Float32\" Name=\"time\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	for (const auto& r: times) {
		for (const auto& t : r) {
			os << t << " ";
		}
	}
	os << "\n</DataArray>\n" << "\n</PointData>\n";

	// CELLDATA (live on the polylines)
	os << "<CellData Scalars=\" CellData\">\n";
	std::vector<std::string> sTypeNames = { "organtype", "id", "subtype"}; //  , "order", "radius"
	for (size_t i=0; i<sTypeNames.size(); i++) {
		os << "<DataArray type=\"Float32\" Name=\"" << sTypeNames[i] <<"\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
		std::vector<double> scalars = getScalar(otype, sTypeNames[i]);
		for (double s : scalars) {
			os << s << " ";
		}
		os << "\n</DataArray>\n";
	}
	os << "\n</CellData>\n";

	// POINTS (=nodes)
	os << "<Points>\n"<<"<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\" >\n";
	for (auto const& r:nodes) {
		for (auto const& n : r) {
			os << n.x << " "<< n.y <<" "<< n.z<< " ";
		}
	}
	os << "\n</DataArray>\n"<< "</Points>\n";

	// LINES (polylines)
	os << "<Lines>\n"<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	int c=0;
	for (auto const& r:organs) {
		for (size_t i=0; i<r->getNumberOfNodes(); i++) {
			os << c << " ";
			c++;
		}
	}
	os << "\n</DataArray>\n"<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	c = 0;
	for (auto const& r:organs) {
		c += r->getNumberOfNodes();
		os << c << " ";
	}
	os << "\n</DataArray>\n";
	os << "\n</Lines>\n";

	os << "</Piece>\n";
	os << "</PolyData>\n" << "</VTKFile>\n";
}

/**
write VTP using tinyXML
 **/
void Plant::TiXMLwriteVTP(int otype, std::ostream & os) const // Write .VTP file by using TinyXML2 performance slowed by 0.5 seconds but precision increased
{
	tinyxml2::XMLPrinter printer( 0, false, 0 );

	this->getOrgans(otype); // update roots (if necessary)
	auto nodes = getPolylines(otype);
	auto times = getPolylinesNET(otype);

	os << "<?xml version=\"1.0\"?>";
	printer.OpenElement("VTKFile"); printer.PushAttribute("type", "PolyData"); printer.PushAttribute("version", "0.1"); printer.PushAttribute("byte_order", "LittleEndian");
	printer.OpenElement("PolyData");
	int non = 0; // number of nodes
	for (auto const& r : organs) {
		non += r->getNumberOfNodes();
	}
	int nol=organs.size(); // number of lines
	printer.OpenElement("Piece"); printer.PushAttribute("NumberOfLines",  nol); printer.PushAttribute("NumberOfPoints", non);

	// POINTDATA
	printer.OpenElement("PointData"); printer.PushAttribute("Scalars", "Pointdata");
	printer.OpenElement("DataArray"); printer.PushAttribute("type", "Float32");  printer.PushAttribute("Name", "time"); printer.PushAttribute("NumberOfComponents", "1"); printer.PushAttribute("format", "ascii" );
	for (std::vector<double> r: times) {
		for (double t : r) {
			printer.PushText(t); printer.PushText(" ");
		}
	}
	printer.CloseElement();
	printer.CloseElement();

	// CELLDATA (live on the polylines)
	printer.OpenElement("CellData"); printer.PushAttribute("Scalars", "CellData" );
	std::vector<std::string> sTypeNames = { "organtype", "id",  "emergencetime", "creationtime", "age", "subtype"}; //  , "order", "radius", "subtype" ,
	for (size_t i=0; i<sTypeNames.size(); i++) {
		std::string sType = sTypeNames[i];
		char const *schar = sType.c_str();
		printer.OpenElement("DataArray"); printer.PushAttribute("type", "Float32");  printer.PushAttribute("Name", schar); printer.PushAttribute("NumberOfComponents", "1"); printer.PushAttribute("format", "ascii" );
		std::vector<double> scalars = getScalar(otype, sTypeNames[i]);
		for (double s : scalars) {
			printer.PushText(s); printer.PushText(" ");
		}
		printer.CloseElement();
	}
	printer.CloseElement();

	// POINTS (=nodes)
	printer.OpenElement("Points");
	printer.OpenElement("DataArray"); printer.PushAttribute("type", "Float32");  printer.PushAttribute("Name", "Coordinates"); printer.PushAttribute("NumberOfComponents", "3"); printer.PushAttribute("format", "ascii" );
	for (auto const& r:nodes) {
		for (auto const& n : r) {
			printer.PushText(n.x); printer.PushText(" "); printer.PushText(n.y); printer.PushText(" "); printer.PushText(n.z); printer.PushText(" ");
		}
	}
	printer.CloseElement();
	printer.CloseElement();

	// LINES (polylines)
	printer.OpenElement("Lines");
	printer.OpenElement("DataArray"); printer.PushAttribute("type", "Float32");  printer.PushAttribute("Name", "connectivity"); printer.PushAttribute("NumberOfComponents", "1"); printer.PushAttribute("format", "ascii" );
	int c=0;
	for (auto const& r:organs) {
		for (size_t i=0; i<r->getNumberOfNodes(); i++) {
			printer.PushText(c); printer.PushText(" ");
			c++;
		}
	}
	printer.CloseElement();

	printer.OpenElement("DataArray"); printer.PushAttribute("type", "Float32");  printer.PushAttribute("Name", "offsets"); printer.PushAttribute("NumberOfComponents", "1"); printer.PushAttribute("format", "ascii" );
	c = 0;
	for (auto const& r:organs) {
		c += r->getNumberOfNodes();
		printer.PushText(c); printer.PushText(" ");
	}
	printer.CloseElement();
	printer.CloseElement();
	printer.CloseElement();
	printer.CloseElement();
	printer.CloseElement();
	os << std::string(printer.CStr());
}




/**
 * Writes the current confining geometry (e.g. a plant container) as Paraview Python script
 * Just adds the initial lines, before calling the method of the SDF.
 *
 * @param os      typically a file out stream
 */
void Plant::writeGeometry(std::ostream & os) const
{
	os << "from paraview.simple import *\n";
	os << "paraview.simple._DisableFirstRenderCameraReset()\n";
	os << "renderView1 = GetActiveViewOrCreate('RenderView')\n\n";
	geometry->writePVPScript(os);
}


} // namespace CPlantBox
