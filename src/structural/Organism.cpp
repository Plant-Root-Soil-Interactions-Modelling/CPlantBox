// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Organism.h"

#include "Organ.h"
#include "Seed.h"
#include "organparameter.h"

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <ctime>
#include <numeric>

namespace CPlantBox {

std::vector<std::string> Organism::organTypeNames = { "organ", "seed", "root", "stem", "leaf" };

int Organism::instances = 0; // number of instances

/**
 * Constructs organism, initializes random number generator
 * @param seednum    option to set seed (for creation of random number) default = 0.
 */
Organism::Organism(unsigned int seednum)
{
    instances++;
    if(seednum >0){
        seed_val = seednum;
    } else {
        seed_val = std::chrono::system_clock::now().time_since_epoch().count()+instances;
    }
    gen = std::mt19937(seed_val);
    plantId = instances; // for debugging
    // std::cout << "Created Organism: " << plantId << "\n" << std::flush;
};


/**
 * @return the organ type number of an organ type name @param name
 */
int Organism::organTypeNumber(std::string name)
{
    for (int ot=0; ot < organTypeNames.size(); ot++ ) {
        if (name.compare(organTypeNames.at(ot))==0) {
            return ot;
        }
    }
    std::cout << "Organism::organTypeNumber: unknown organ type name " << name << "\n" << std::flush;
    return -1;
}

/**
 * @return the organ type name of an organ type number @param ot
 */
std::string Organism::organTypeName(int ot)
{
    try {
        return organTypeNames.at(ot);
    } catch (const std::exception& e) {
        throw std::invalid_argument("Organism::organTypeName: unknown organ type number "+ std::to_string(ot));
    }
}

/**
 * Deep copies the organism
 */
std::shared_ptr<Organism> Organism::copy()
{
    std::cout << "Warning Organism (copy) should not be called directly: " << plantId << "\n" << std::flush;
    auto no = std::make_shared<Organism>(*this); // copy constructor
    for (int i = 0; i < baseOrgans.size(); i++) {
        no->baseOrgans[i] = baseOrgans[i]->copy(no);
    }
    for (int ot = 0; ot < numberOfOrganTypes; ot++) { // copy organ type parameters
        for (auto& otp : no->organParam[ot]) {
            // no->setOrganRandomParameter(otp.second->copy(no));
        	std::shared_ptr<OrganRandomParameter> new_params= otp.second->copy(no);
        	// std::cout << "Plant::copy() " << new_params->plant.lock()->plantId << std::flush <<  "\n";
			no->setOrganRandomParameter(new_params);
        }
    }
    return no;
}

/**
 * Copies the organ random parameters of all sub types of one specific organ type into a vector
 *
 * @param ot    the organ type
 * @return      all random parameters as a list (list index != subType)
 */
std::vector<std::shared_ptr<OrganRandomParameter>> Organism::getOrganRandomParameter(int ot) const
    {
    auto  otps = std::vector<std::shared_ptr<OrganRandomParameter>>(0);
    for (auto& otp : organParam[ot]) {
        otps.push_back(otp.second);
    }
    return otps;
    }

/**
 * Returns an organ random parameter of a specific organ type and sub type
 *
 * @param ot       the organ type (e.g. ot_root)
 * @param subType  the sub type (e.g. root type)
 * @return         the respective random parameter
 */
std::shared_ptr<OrganRandomParameter> Organism::getOrganRandomParameter(int ot, int subType) const
{
    try {
        //            std::cout << "reading organ type " << ot << " sub type " << subType <<": ";
        //           for (auto& p : organParam.at(ot)) {
        //              std::cout << p.first;
        //          }
        //          std::cout << "\n" << std::flush;
        return organParam.at(ot).at(subType);
    } catch(...) { // const std::out_of_range& oor
        throw std::invalid_argument("Organism::getOrganTypeParameter: OrganRandomParameter for "+ Organism::organTypeName(ot) + ", of sub type " +
            std::to_string(subType) + " was not set");
    }
}

/**
 *  Sets the random parameter.
 *  subType and organType are defined within @param p.
 *
 *  @param p    the organ random parameter
 */
void Organism::setOrganRandomParameter(std::shared_ptr<OrganRandomParameter> p)
{
    assert(p->plant.lock().get()==this && "OrganTypeParameter::plant should be this organism");
    organParam[p->organType][p->subType] = p;
    // std::cout << "setting organ type " << p->organType << ", sub type " << p->subType << ", name "<< p->name << "\n";
}

/**
 * @return Get the organ sub type by its name
 *
 * @param organtype 		the organ type
 * @param str 				the sub type name
 */

int Organism::getParameterSubType(int organtype, std::string str) const
{
    auto orp = this->getOrganRandomParameter(organtype);
    for (auto& o :orp) {
        if (o->name == str) {
            return o->subType;
        }
    }
    return -1;
}

/**
 * Overwrite if there is the need for additional initializations,
 * before the simulation starts.
 *
 * e.g. for initialization of GrowthFunctions, TropismFunctions, set up base Organs
 */
void Organism::initialize(bool verbose)
{ }

/**
 * Simulates the development of the organism in a time span of @param dt days.
 *
 * @param dt        time step [day]
 * @param verbose   turns console output on or off
 */
void Organism::simulate(double dt, bool verbose)
{
    if (verbose) {
        std::cout << "Organism::simulate: from "<< simtime << " to " << simtime+dt << " days" << std::endl;
    }
    this->dt = dt;
    oldNumberOfNodes = getNumberOfNodes();
    oldNumberOfOrgans = getNumberOfOrgans();
    for (const auto& r : baseOrgans) {
        r->simulate(dt, verbose);
    }
    simtime+=dt;
}

/**
 * Creates a sequential list of organs. Considers only organs with more than 1 node.
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default).
 * @param all       get also the organs with only one node? default: false.  Sometimes true for carbon-limited growth
 * @return Sequential list of organs. For all == false, if there is less than one node,
 * or another organ type is expected, an empty vector is returned.
 */
std::vector<std::shared_ptr<Organ>> Organism::getOrgans(int ot, bool all) const
{
    auto organs = std::vector<std::shared_ptr<Organ>>(0);
    organs.reserve(getNumberOfOrgans()); // just for speed up
    for (const auto& o : this->baseOrgans) {
        o->getOrgans(ot, organs, all);
    }
    return organs;
}

/**
 * Returns a single scalar parameter for each organ as sequential list,
 * corresponding to the sequential organ list, see Organism::getOrgans.
 *
 * This method is mostly for post processing, since it is flexible but slow.
 *
 * @param name      name of the parameter
 * @param ot        the expected organ type, where -1 denotes all organ types (default).
 * @param organs    optionally, a predefined sequential organ list can be used (@param ot is ignored in this case)
 * @return A vector of one parameter values per each organ, if unknown NaN
 */
std::vector<double> Organism::getParameter(std::string name, int ot, std::vector<std::shared_ptr<Organ>> organs) const
{
    if (organs.empty()) {
        organs = getOrgans(ot);
    }
    std::vector<double> p = std::vector<double>(organs.size());
    for (int i=0; i<organs.size(); i++) {
        p[i] = organs[i]->getParameter(name);
    }
    return p;
}

/**
 * Returns the summed parameter, obtained by Organism::getParameters (e.g. getSummed("length"))
 *
 * @param name      name of the parameter
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          the summed up value
 */
double Organism::getSummed(std::string name, int ot) const
{
    auto v = getParameter(name, ot);
    return std::accumulate(v.begin(), v.end(), 0.0);
}

/**
 * Returns the organisms number of segments of a specific organ type
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          total number of segments in the organism of type ot
 */
int Organism::getNumberOfSegments(int ot) const
{
    int s=0;
    auto organs = getOrgans(ot);
    for (const auto& o : organs) {
        s += o->getNumberOfSegments();
    }
    return s;
}

/**
 * Represents the organ's nodes as polyline
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          for each organ a vector of nodes
 */
std::vector<std::vector<Vector3d>> Organism::getPolylines(int ot) const
    {
    auto organs = getOrgans(ot);
    std::vector<std::vector<Vector3d>> nodes = std::vector<std::vector<Vector3d>>(organs.size());
    for (size_t j=0; j<organs.size(); j++) {
        std::vector<Vector3d>  n = std::vector<Vector3d>(organs[j]->getNumberOfNodes());
        for (size_t i=0; i<organs[j]->getNumberOfNodes(); i++) { // loop over all nodes of all organs
            n.at(i) = organs[j]->getNode(i);
        }
        nodes[j] = n;
    }
    return nodes;
    }

/**
 * The corresponding node creation times to the organ polyline representation (@see Organism::getPolylines)
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          for each organ a vector of nodes
 */
std::vector<std::vector<double>> Organism::getPolylineCTs(int ot) const
    {
    auto organs = getOrgans(ot);
    std::vector<std::vector<double>> nodes = std::vector<std::vector<double>>(organs.size());
    for (size_t j=0; j<organs.size(); j++) {
        std::vector<double>  nct = std::vector<double>(organs[j]->getNumberOfNodes());
        for (size_t i=0; i<organs[j]->getNumberOfNodes(); i++) { // loop over all nodes of all organs
            nct.at(i) = organs[j]->getNodeCT(i);
        }
        nodes[j] = nct;
    }
    return nodes;
    }

/**
 * All nodes of emerged organs are ordered by their node index,
 * initial nodes of base organs are copied, even if not emerged.
 *
 * @return          a vector of nodes
 */
std::vector<Vector3d> Organism::getNodes() const
{
    auto organs = getOrgans();
    std::vector<Vector3d> nv = std::vector<Vector3d>(getNumberOfNodes()); // reserve big enough vector
    for (const auto& o : baseOrgans) { // copy initial nodes (even if organs have not developed)
        nv.at(o->getNodeId(0)) = o->getNode(0);
    }
    for (const auto& o : organs) { // copy all organ nodes
        for (size_t i = 1; i<o->getNumberOfNodes(); i++) { // since all base nodes are stored twice as base and along root
            nv.at(o->getNodeId(i)) = o->getNode(i);
        }
    }
    return nv;
}

/**
 * The organim's node creation times of a specific organ type corresponding to Organism::getNodes.
 *
 * At a branching point, there are two creation times attached to one node,
 * the creation time of the base node, and the emergence time of the lateral root.
 * This method copies the node creation time of time of the base root.
 *
 * @return          a vector of node creation times
 */
std::vector<double> Organism::getNodeCTs() const
{
    auto organs = getOrgans();
    std::vector<double> cts = std::vector<double>(getNumberOfNodes()); // reserve big enough vector
    for (const auto& o : baseOrgans) { // copy initial nodes (even if organs have not developed)
        cts.at(o->getNodeId(0)) = o->getNodeCT(0);
    }
    for (const auto& o : organs) { // copy all organ creation times
        for (size_t i=1; i<o->getNumberOfNodes(); i++) {
            // we start at 1 because we don't copy the emergence time of the lateral root
            cts.at(o->getNodeId(i)) = o->getNodeCT(i);
        }
    }
    return cts;
}

/**
 * The line segments of the organism, each segment consisting of two node indices,
 * corresponding to the node list Organism::getNode(-1),
 * or the node indices from Organism::getNodeId(ot).
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          line segments
 */
std::vector<Vector2i> Organism::getSegments(int ot) const
{
    auto organs = getOrgans(ot);
    std::vector<Vector2i> segs = std::vector<Vector2i>(0);
    segs.reserve(this->getNumberOfSegments(ot)); // for speed up
    for (const auto& o : organs) {
        auto s = o->getSegments();
        segs.insert(segs.end(), s.begin(), s.end()); // append s; todo check if it works..
    }
    return segs;
}

/**
 * The creation time per segment, corresponding to Organism::getSegments(ot).
 * The segment creation time equals the node creation time of the second node.
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          creation times of each segment
 */
std::vector<double> Organism::getSegmentCTs(int ot) const
{
    auto nodeCT = getNodeCTs(); // of all nodes (otherwise indices are tricky)
    auto segs = getSegments(ot);
    std::vector<double> cts = std::vector<double>(segs.size());
    for (int i=0; i<cts.size(); i++) {
        cts[i] = nodeCT[segs[i].y]; // segment creation time is the node creation time of the second node
    }
    return cts;
}
/**
 * The segment indexes, corresponding to Organism::getSegments(ot).
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          Id of each segment
 */
std::vector<int> Organism::getSegmentIds(int ot) const
{
    auto segments = getSegments(ot); // of all nodes (otherwise indices are tricky)
    std::vector<int> segId = std::vector<int>(segments.size());
    for (int i=0; i<segments.size(); i++) {
        segId[i] = segments[i].y-1;
    }
    return segId;
}


/**
 * A vector of pointers to the organs containing each segment, corresponding to Organism::getSegments.
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          creation times of each segment
 */
std::vector<std::shared_ptr<Organ>> Organism::getSegmentOrigins(int ot) const
    {
    auto organs = getOrgans(ot);
    auto segs = std::vector<std::shared_ptr<Organ>>(0);
    segs.reserve(this->getNumberOfSegments(ot)); // for speed up
    for (const auto& o : organs) {
        auto s = o->getSegments();
        for (int i=0; i<s.size(); i++) {
            segs.push_back(o);
        }
    }
    return segs;
    }

/**
 * @return the indices of the nodes that were moved during the last time step,
 * update the node coordinates using Organism::getUpdatedNodes(),
 * and creation times using Organism::getUpdatedNodeCTs()
 * for stem and leaves, assume that all the old nodes need to be updated
 */
std::vector<int> Organism::getUpdatedNodeIndices() const
{
    auto organs = this->getOrgans();
    std::vector<int> ni = std::vector<int>(0);
    for (const auto& o : organs) {
        if (o->hasMoved()&&(o->getOldNumberOfNodes()>1)) {
            if((o->organType() >2)){//is stem or leaf
                int onon =  o->getOldNumberOfNodes();//because of tropism and internodal growth, all nodes can move
                for(int i = 1; i < onon; i++){
                    ni.push_back(o->getNodeId(i));
                }
            }else{ni.push_back(o->getNodeId(o->getOldNumberOfNodes()-1));}
        }
    }
    return ni;
}

/**
 * @return the new coordinates of nodes that were updated during the last time step,
 * corresponding to Organism::getUpdatedNodeIndices
 * for stem and leaves, assume that all the old nodes need to be updated
 */
std::vector<Vector3d> Organism::getUpdatedNodes() const
{
    auto organs = this->getOrgans();
    std::vector<Vector3d> nv = std::vector<Vector3d>(0);
    for (const auto& o : organs) {
        if (o->hasMoved()&&(o->getOldNumberOfNodes()>1)) {
            if((o->organType() > 2)){
                int onon =  o->getOldNumberOfNodes();//in case of internodal growth, not just last node needs to be updated
                for(int i = 1; i < onon; i++){
                    nv.push_back(o->getNode(i));
                }
            }else{nv.push_back(o->getNode(o->getOldNumberOfNodes()-1));}
        }
    }
    return nv;
}

/**
 * @return the new creation times of nodes that were updated during the last time step,
 * corresponding to Organism::getUpdatedNodeIndices
 */
std::vector<double> Organism::getUpdatedNodeCTs() const
{
    auto organs = this->getOrgans();
    std::vector<double> nv = std::vector<double>(0);
    for (const auto& o : organs) {
        if (o->hasMoved()&&(o->getOldNumberOfNodes()>1)) {
            if((o->organType() > 2)){
                int onon =  o->getOldNumberOfNodes();//in case of internodal growth, not just last node needs to be updated
                for(int i = 1; i < onon; i++){
                    nv.push_back(o->getNodeCT(i));
                }
            }else{nv.push_back(o->getNodeCT(o->getOldNumberOfNodes()-1));}
        }
    }
    return nv;
}

/**
 * @return a vector of all nodes created during the last time step
 * to dynamically add to the old node list, see also RootSystem::getNodes()
 */
std::vector<Vector3d> Organism::getNewNodes() const
{
    auto organs = this->getOrgans();
    std::vector<Vector3d> nv(this->getNumberOfNewNodes());
    for (const auto& o : organs) {
        int onon = o->getOldNumberOfNodes();
        for (size_t i=onon; i<o->getNumberOfNodes(); i++) { // loop over all new nodes
            nv.at(o->getNodeId(i)-this->oldNumberOfNodes) = o->getNode(i);
        }
    }
    return nv;
}

/**
 * @return a vector of the creation times of the nodes created in the last time step,
 * initially start with a Organism::getNodeCTs() call
 *
 * At a branching point we copy the creation time of the base root @see Organism::getNodeCTs()
 */
std::vector<double> Organism::getNewNodeCTs() const
{
    auto organs = this->getOrgans();
    std::vector<double> nv(this->getNumberOfNewNodes());
    for (const auto& o : organs) {
        int onon = o->getOldNumberOfNodes();
        if (onon==0) { // make sure to never copy the root emergence time
            onon = 1;
        }
        for (size_t i=onon; i<o->getNumberOfNodes(); i++) { // loop over all new nodes
            nv.at(o->getNodeId(i)-this->oldNumberOfNodes) = o->getNodeCT(i);
        }
    }
    return nv;
}

/**
 * Creates a vector of new created segments that were created during the last time step
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          a vector of newly created segments
 */
std::vector<Vector2i> Organism::getNewSegments(int ot) const
{
    auto organs = this->getOrgans(ot);
    std::vector<Vector2i> si = std::vector<Vector2i>(0);
    si.reserve(this->getNumberOfNewNodes());
    for (const auto& o :organs) {
        int onon = o->getOldNumberOfNodes();
        for (size_t i=onon-1; i<o->getNumberOfNodes()-1; i++) { // loop over new segments
            Vector2i v(o->getNodeId(i),o->getNodeId(i+1));
            si.push_back(v);
        }
    }
    return si;
}

/**
 * Creates a vector of pointers to the organ class containing the segments,
 * for each newly created segment of organ type @param ot,
 * corresponding to the segments obtained by Organism::getNewSegments(ot)
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          a vector of pointers to organs
 */
std::vector<std::shared_ptr<Organ>> Organism::getNewSegmentOrigins(int ot) const
    {
    auto organs = this->getOrgans(ot);
    std::vector<std::shared_ptr<Organ>> so;
    for (auto& o :organs) {
        int onon = o->getOldNumberOfNodes();
        for (size_t i=onon-1; i<o->getNumberOfNodes()-1; i++) { // loop over new segments
            so.push_back(o);
        }
    }
    return so;
    }


int Organism::getDelayDefinition(int ot_lat)
{
	auto srp = std::static_pointer_cast<SeedRandomParameter>(this->getOrganRandomParameter(Organism::ot_seed,0 ));
	if((ot_lat ==  Organism::ot_stem)||(ot_lat ==  Organism::ot_leaf)){return srp->delayDefinitionShoot;
	}else{return srp->delayDefinition;}
}

/**
 * @return Quick info about the object for debugging
 */
std::string Organism::toString() const
{
    std::stringstream str;
    str << "Organism with "<< baseOrgans.size() <<" base organs, " << getNumberOfNodes()
                                                		        << " nodes, and a total of " << getNumberOfOrgans() << " organs, after " << getSimTime() << " days";
    return str.str();
}

/**
 * Polymorphic XML parameter file reader:
 * adds all organ parameter types of a XML file to the organism's parameters
 *
 * Each OrganTypeParameter must have already one prototype defined
 * in the organ type parameters Organism::organParam, since
 * the organ type parameters are created by OrganTypeParameter::copy.
 *
 * @param name      file name or data as string
 * @param basetag   name of the base tag (e.g. "organism", or "plant")
 * @param fromFile  @param name == file name (true) or data as string (false)
 */
void Organism::readParameters(std::string name, std::string basetag, bool fromFile, bool verbose)
{
    tinyxml2::XMLDocument doc;
    if(fromFile){doc.LoadFile(name.c_str()); //open xml file and read data
    }else{doc.Parse((const char*)name.c_str());} //get data directly from string
    if(doc.ErrorID() == 0) {
        tinyxml2::XMLElement* base = doc.FirstChildElement(basetag.c_str());
        if(base != nullptr){
            auto p = base->FirstChildElement();
            while((p!=nullptr) && (p->Name()!=nullptr)) {

                std::string tagname = p->Name();
                int ot = Organism::organTypeNumber(tagname);
                std::shared_ptr<OrganRandomParameter> prototype;
                if (organParam[ot].count(0)) { // is the prototype defined?
                    prototype = organParam[ot][0];
                } else {
                    prototype = nullptr;
                }

                if ((ot==0) && (prototype==nullptr)) { // read depricated xml format, in case organ is not used
                    tagname = p->Attribute("type");
                    ot = Organism::organTypeNumber(tagname);
                    if (ot>0) { // tagname known?
                        if (organParam[ot].count(0)) { // is the prototype defined?
                            prototype = organParam[ot][0];
                        } else {
                            prototype = nullptr;
                        }
                    } else {
                        prototype = nullptr;
                    }
                }

                if (prototype!=nullptr) { // read prototype
                    auto otp = prototype->copy(shared_from_this());
                    otp->readXML(p, verbose);
                    otp->organType = ot; // in depricated case, readXML will a give wrong value
                    setOrganRandomParameter(otp);
                } else { // skip prototype
                    if(verbose){
                        std::cout << "Organism::readParameters: warning, skipping " << tagname <<
                            ", no random parameter class defined, use initializeReader()\n" << std::flush;
                    }
                }
                p = p->NextSiblingElement();
            } // while
        } else {
            if (basetag.compare("plant") == 0) { // try old spelling
                if(verbose){
                    std::cout << "Organism::readParameters: plant tag was not found in xml file, retrying with Plant " << std::endl;
                }
                readParameters(name, "Plant", fromFile, verbose); // rerun
                return;
            }
            throw std::invalid_argument ("Organism::readParameters: " + std::string(basetag.c_str()) + " tag was not found in xml file");
        }
    } else {
        std::cout << "Organism::readParameters: could not open file " << name << "\n" << std::flush;
    }
}

/**
 * XML parameter file writer
 * writes all organ parameter types into a XML File
 *
 * @param name      file name
 * @param basetag   name of the base tag (e.g. "organism", or "plant") *
 * @param comments  write parameter descriptions
 */
void Organism::writeParameters(std::string name, std::string basetag, bool comments) const
{
    std::setlocale(LC_NUMERIC, "en_US.UTF-8");
    tinyxml2::XMLDocument xmlDoc;
    tinyxml2:: XMLElement* xmlParams = xmlDoc.NewElement(basetag.c_str()); // RSML
    for (int ot = 0; ot < numberOfOrganTypes; ot++) {
        for (auto& otp : organParam[ot]) {
            if ((ot!=ot_seed) && (otp.second->subType > 0)) { // subType 0 is unused, and used as a prototype for reading the xml tags
                xmlParams->InsertEndChild(otp.second->writeXML(xmlDoc, comments));
            }
            if (ot==ot_seed) {
                xmlParams->InsertEndChild(otp.second->writeXML(xmlDoc, comments));
            }
        }
    }
    xmlDoc.InsertEndChild(xmlParams);
    xmlDoc.SaveFile(name.c_str());
}

/**
 * Exports the simulation results with the type from the extension in name
 * (that must be lower case)
 *
 * @param name      file name e.g. output.vtp
 */
void Organism::write(std::string name) const
{
    std::string ext = name.substr(name.size()-3,name.size()); // pick the right write}r
    if (ext.compare("sml")==0) {
        std::cout << "writing RSML... "<< name.c_str() <<"\n";
        writeRSML(name); // use base class writer
    } else if (ext.compare("vtp")==0) {
        std::cout << "writing VTP... "<< name.c_str() <<"\n";
        std::ofstream fos;
        fos.open(name.c_str());
        writeVTP(-1, fos);
        fos.close();
    } else if (ext.compare(".py")==0)  {
        std::cout << "writing Geometry ... "<< name.c_str() <<"\n";
        std::ofstream fos;
        fos.open(name.c_str());
        writeGeometry(fos);
        fos.close();
    } else {
        throw std::invalid_argument("Plant::write(): Unkwown file type");
    }
}

/**
write VTP using tinyXML
 **/
void Organism::writeVTP(int otype, std::ostream & os) const // Write .VTP file by using TinyXML2 performance slowed by 0.5 seconds but precision increased
{
    tinyxml2::XMLPrinter printer( 0, false, 0 );

    auto organs = this->getOrgans(otype); // update roots (if necessary)
    auto nodes = getPolylines(otype);
    auto times = getPolylineCTs(otype);

    os << "<?xml version=\"1.0\"?>";
    printer.OpenElement("VTKFile"); printer.PushAttribute("type", "PolyData"); printer.PushAttribute("version", "0.1"); printer.PushAttribute("byte_order", "LittleEndian");
    printer.OpenElement("PolyData");
    int non = 0; // number of nodes
    for (const auto& r : organs) {
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
    std::vector<std::string> sTypeNames = { "organType", "id", "creationTime", "age", "subType", "order", "radius"}; //  , "order", "radius", "subtype" ,
    for (size_t i=0; i<sTypeNames.size(); i++) {
        std::string sType = sTypeNames[i];
        const char *schar = sType.c_str();
        printer.OpenElement("DataArray"); printer.PushAttribute("type", "Float32");  printer.PushAttribute("Name", schar); printer.PushAttribute("NumberOfComponents", "1"); printer.PushAttribute("format", "ascii" );
        std::vector<double> scalars = getParameter(sTypeNames[i], otype);
        for (double s : scalars) {
            printer.PushText(s); printer.PushText(" ");
        }
        printer.CloseElement();
    }
    printer.CloseElement();

    // POINTS (=nodes)
    printer.OpenElement("Points");
    printer.OpenElement("DataArray"); printer.PushAttribute("type", "Float32");  printer.PushAttribute("Name", "Coordinates"); printer.PushAttribute("NumberOfComponents", "3"); printer.PushAttribute("format", "ascii" );
    for (const auto& r : nodes) {
        for (const auto& n : r) {
            printer.PushText(n.x); printer.PushText(" "); printer.PushText(n.y); printer.PushText(" "); printer.PushText(n.z); printer.PushText(" ");
        }
    }
    printer.CloseElement();
    printer.CloseElement();

    // LINES (polylines)
    printer.OpenElement("Lines");
    printer.OpenElement("DataArray"); printer.PushAttribute("type", "Int32");  printer.PushAttribute("Name", "connectivity"); printer.PushAttribute("NumberOfComponents", "1"); printer.PushAttribute("format", "ascii" );
    int c=0;
    for (const auto& r : organs) {
        for (size_t i=0; i<r->getNumberOfNodes(); i++) {
            printer.PushText(c); printer.PushText(" ");
            c++;
        }
    }
    printer.CloseElement();

    printer.OpenElement("DataArray"); printer.PushAttribute("type", "Int32");  printer.PushAttribute("Name", "offsets"); printer.PushAttribute("NumberOfComponents", "1"); printer.PushAttribute("format", "ascii" );
    c = 0;
    for (const auto& r : organs) {
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
 * Writes the current confining geometry (e.g. a plant container) as paraview python script
 * Just adds the initial lines, before calling the method of the sdf.
 *
 * todo move to Organism (including geometry)
 *
 * @param os      typically a file out stream
 */
void Organism::writeGeometry(std::ostream & os) const
{
    os << "from paraview.simple import *\n";
    os << "paraview.simple._DisableFirstRenderCameraReset()\n";
    os << "renderView1 = GetActiveViewOrCreate('RenderView')\n\n";
    geometry->writePVPScript(os);
}


/**
 * Creates a rsml file with filename @param name.
 *
 * @param name      name of the rsml file
 */
void Organism::writeRSML(std::string name) const
{
    std::setlocale(LC_NUMERIC, "en_US.UTF-8");
    tinyxml2::XMLDocument xmlDoc;
    tinyxml2:: XMLElement* rsml = xmlDoc.NewElement("rsml"); // RSML
    tinyxml2:: XMLElement* meta = getRSMLMetadata(xmlDoc);
    tinyxml2:: XMLElement* scene = getRSMLScene(xmlDoc);
    rsml->InsertEndChild(meta);
    rsml->InsertEndChild(scene);
    xmlDoc.InsertEndChild(rsml);
    xmlDoc.SaveFile(name.c_str());
}

/**
 * @return the meta tag of the rsml file
 */
tinyxml2:: XMLElement* Organism::getRSMLMetadata(tinyxml2::XMLDocument& xmlDoc) const
{
    tinyxml2:: XMLElement* metadata = xmlDoc.NewElement("metadata"); // META
    tinyxml2:: XMLElement* version = xmlDoc.NewElement("version");
    version->SetText(1);
    tinyxml2:: XMLElement* unit = xmlDoc.NewElement("unit");
    unit->SetText("cm");
    tinyxml2:: XMLElement* resolution = xmlDoc.NewElement("resolution");
    resolution->SetText(1);
    tinyxml2:: XMLElement* last = xmlDoc.NewElement("last-modified");
    std::time_t t = std::time(0);
    std::tm* now = std::localtime(&t);
    std::string s = std::to_string(now->tm_mday)+"-"+std::to_string(now->tm_mon+1)+"-"+std::to_string(now->tm_year + 1900);
    last->SetText(s.c_str());
    tinyxml2:: XMLElement* software = xmlDoc.NewElement("software");
    software->SetText("OrganicBox");
    // todo no image tag (?)
    // todo property-definitions
    // todo time sequence (?)
    metadata->InsertEndChild(version);
    metadata->InsertEndChild(unit);
    metadata->InsertEndChild(resolution);
    metadata->InsertEndChild(last);
    metadata->InsertEndChild(software);
    // todo insert remaining tags
    return metadata;
}

/**
 * @return the scene tag of the RSML document, calls base organs to write their tags
 */
tinyxml2:: XMLElement* Organism::getRSMLScene(tinyxml2::XMLDocument& xmlDoc) const
{
    tinyxml2:: XMLElement* scene = xmlDoc.NewElement("scene");
    tinyxml2:: XMLElement* plant = xmlDoc.NewElement("plant");

    //	auto p = plant.lock();
    //	auto stemP = p->getOrganRandomParameter(Organism::ot_stem);
    //	bool plantBox = stemP.size()>0;
    //	if (verbose) {
    //		if (plantBox) {
    //			std::cout << "Seed::initialize: Plant \n";
    //		} else {
    //			std::cout << "Seed::initialize: RootSystem \n";
    //		}
    //	}

    scene->InsertEndChild(plant);
    for (auto& o: baseOrgans) {
        o->writeRSML(xmlDoc, plant);
    }
    return scene;
}

/**
 * Sets the seed of the organisms random number generator.
 * In order to obtain two exact same organisms call before Organism::initialize().
 *
 * @param seed      the random number generator seed
 */
void Organism::setSeed(unsigned int seed)
{
    this->gen = std::mt19937(seed);
}

/**
 * Returns the seed of the plant
 */
std::shared_ptr<Seed> Organism::getSeed()
{
    return std::static_pointer_cast<Seed>(baseOrgans.at(0));
}

} // namespace
