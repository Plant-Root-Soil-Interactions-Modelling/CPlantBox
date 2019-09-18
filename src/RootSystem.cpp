// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "RootSystem.h"

#include "organparameter.h"
#include "Organism.h"
#include "Seed.h"

namespace CRootBox {

/**
 * Sets up the xml reader Organism::readParameters
 */
RootSystem::RootSystem(): Organism()
{
   auto rrp = new RootRandomParameter(this);
   rrp->subType = 0;
   setOrganRandomParameter(rrp);
   auto srp = new SeedRandomParameter(this);
   srp->subType = 0;
   setOrganRandomParameter(srp);
}

/**
 * Copy Constructor
 *
 * deep copies the root system
 * does not deep copy geometry, elongation functions, and soil (all not owned by rootsystem)
 * empties buffer
 *
 * @param rs        root system that is copied
 */
RootSystem::RootSystem(const RootSystem& rs): Organism(rs), geometry(rs.geometry), soil(rs.soil)
{
    std::cout << "Copying root system with "<<rs.baseOrgans.size()<< " base roots \n";
    roots.clear();
}

/**
 * @return the i-th root parameter of sub type @param type.
 */
RootRandomParameter* RootSystem::getRootTypeParameter(int type) const
{
    return (RootRandomParameter*) getOrganRandomParameter(Organism::ot_root, type);
}

/**
 * @return all root type parameters in a vector
 */
std::vector<RootRandomParameter*> RootSystem::getRootTypeParameter() const
{
    std::vector<RootRandomParameter*>  otps = std::vector<RootRandomParameter*>(0);
    for (auto& otp : organParam[Organism::ot_root]) {
        otps.push_back((RootRandomParameter*)otp.second);
    }
    return otps;
}

/**
 * Sets the root system parameters @param rsp
 */
void RootSystem::setRootSystemParameter(SeedRandomParameter& rsp)
{
    assert(rsp.subType==0 && "RootSystem::setRootSystemParameter: In CRootBox must have subType 0");
    this->setOrganRandomParameter(&rsp);
}

/**
 * @return the root system parameter
 */
SeedRandomParameter* RootSystem::getRootSystemParameter()
{
    return (SeedRandomParameter*)this->getOrganRandomParameter(ot_seed,0 );
}

/**
 * Resets the root system: deletes all roots, sets simulation time to 0.
 */
void RootSystem::reset()
{
    for(auto b :baseOrgans) {
        delete b;
    }
    baseOrgans.clear();
    simtime = 0;
    organId = -1;
    nodeId = -1;
}

/**
 * Reads the root parameter from a file. Opens plant parameters with the same filename if available,
 * otherwise assumes a tap root system at position (0,0,-3).
 *
 * @param name          filename without file extension
 * @param subdir        directory ("modelparameter/" by default)
 */
void RootSystem::openFile(std::string name, std::string subdir)
{
	std::cout << "RootSystem::openFile is deprecated, use readParameter instead \n";
    std::ifstream fis;
    // open root parameter
    std::string rp_name = subdir;
    rp_name.append(name);
    rp_name.append(".rparam");
    fis.open(rp_name.c_str());
    if (fis.good()) { // did it work?
        readParameters(fis);
        fis.close();
    } else {
        std::string s = "RootSystem::openFile() could not open root parameter file ";
        throw std::invalid_argument(s.append(rp_name));
    }
    // std::cout << "Read " << c << " root type parameters \n"; // debug

    // open plant parameter
    std::string pp_name = subdir;
    pp_name.append(name);
    pp_name.append(".pparam");
    fis.open(pp_name.c_str());
    auto randomSeed = new SeedRandomParameter(this);
    randomSeed->subType = 0;
    if (fis.good()) { // did it work?
        randomSeed->read(fis); // reads the random parameters
        fis.close();
    } else { // create a tap root system
        std::cout << "No root system parameters found, using default tap root system \n";
    }
    this->setRootSystemParameter(*randomSeed);
}

/**
 * Reads root type parameter from an input stream @param is
 * (there is a Matlab script exporting these, @see writeParams.m)
 *
 * @param cin  in stream
 */
int RootSystem::readParameters(std::istream& is)
{
	std::cout << "RootSystem::readParameters(std::istream) is deprecated, use readParameters(std::string filename) instead \n";
    int c = 0;
    while (is.good()) {
        RootRandomParameter* p = new RootRandomParameter(this);
        p->read(is);
        p->organType = Organism::ot_root;
        setOrganRandomParameter(p);
        c++;
    }
    return c;
}

/**
 * Writes root type parameters to an output stream @param os
 *
 * @param os  out stream
 */
void RootSystem::writeParameters(std::ostream& os) const
{
	std::cout << "RootSystem::writeParameters is deprecated, use writeParameters(std::string filename) instead \n";
    for (auto& otp :organParam[Organism::ot_root]) {
        ((RootRandomParameter*)otp.second)->write(os);
    }
}
/**
 * Sets up the base roots according to the plant parameters,
 * a confining geometry, the tropism functions, and the growth functions.
 *
 * Call this method before simulation and after setting geometry, plant and root parameters
 *
 * @parm basaltype      the type of the basal roots (default = 4)
 * @parm basaltype      the type of the shootborne roots (default = 5)
 */
void RootSystem::initialize(int basaltype, int shootbornetype)
{
    reset(); // just in case

    // introduce an extra node at nodes[0]
    getNodeIndex(); // increase node index

    // create seed
    Seed seed = Seed(this);
    seed.initialize();
    seedParam = SeedSpecificParameter(*seed.param()); // copy the specific parameters
    // std::cout << "RootSystem::initialize:\n" <<  seedParam.toString() ;
    baseOrgans = seed.copyBaseOrgans();

    oldNumberOfNodes = baseOrgans.size();
    initCallbacks();
}

/**
 * Called by RootSystem::initialize.
 * Sets up tropism and growth functions call backs using
 * RootSystem::createTropismFunction and RootSystem::createGrowthFunction
 */
void RootSystem::initCallbacks()
{
    // Create tropisms and growth functions per root type
    for (auto& p_otp :organParam[Organism::ot_root]) {
        RootRandomParameter* rtp = (RootRandomParameter*)p_otp.second;
        Tropism* tropism = this->createTropismFunction(rtp->tropismT, rtp->tropismN, rtp->tropismS);
        tropism->setGeometry(geometry);
        delete rtp->f_tf; // delete old tropism
        rtp->f_tf = tropism; // set new one
        GrowthFunction* gf_ = this->createGrowthFunction(rtp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        delete rtp->f_gf;
        rtp->f_gf  = gf_;
    }
}

/**
 * Manually sets a tropism function for a specific or for all root types.
 * Must be called after RootSystem::initialize(), otherwise its overwritten.
 *
 * @param tf_           a tropism
 * @param rt            root type, if rt = -1 all types are set to this tropism (default).
 */
void RootSystem::setTropism(Tropism* tf_, int rt)
{
    if (rt>-1) { // set for a specific root type
        getRootTypeParameter(rt)->f_tf=tf_;
    } else { // set for all root types (default)
        for (auto& p_otp :organParam[Organism::ot_root]) {
            RootRandomParameter* rtp = (RootRandomParameter*)p_otp.second;
            rtp->f_tf = tf_;
        }
    }
}

/**
 * Simulates root system growth for time span dt
 *
 * @param dt    	time step [day]
 * @param verbose 	indicates if status is written to the console (cout) (default = false)
 */
void RootSystem::simulate(double dt, bool verbose)
{
    Organism::simulate(dt,verbose);
    roots.clear(); // empty buffer
}

/**
 * Simulates root system growth for the time span defined in the root system parameters
 */
void RootSystem::simulate()
{
    this->simulate(seedParam.simtime);
}

/**
 * Simulates root system growth for a time span, elongates a maximum of @param maxinc total length [cm/day]
 * using the proportional elongation @param se to impede overall growth.
 *
 * @param dt        time step [day]
 * @param maxinc_   maximal total length [cm/day] the root system is allowed to grow in this time step
 * @param se        The class ProportionalElongation is used to scale overall root growth
 * @param verbose   indicates if status is written to the console (cout) (default = false)
 */
void RootSystem::simulate(double dt, double maxinc_, ProportionalElongation* se, bool verbose)
{
    const double accuracy = 1.e-3;
    const int maxiter = 20;
    double maxinc = dt*maxinc_; // [cm]
    double ol = getSummed("length");
    int i = 0;

    push();
    se->setScale(1.);
    simulate(dt, verbose);
    double l = getSummed("length");
    double inc_ = l - ol;
    if (verbose) {
        std::cout << "expected increase is " << inc_ << " maximum is " << maxinc
            << "\n";
    }
    pop();

    if ((inc_>maxinc) && (std::abs(inc_-maxinc)>accuracy)) { // check if we have to perform a binary search

        double sl = 0.; // left
        double sr = 1.; // right

        while ( ((std::abs(inc_-maxinc)) > accuracy) && (i<maxiter) )  { // binary search

            double m = (sl+sr)/2.; // mid
            push();
            se->setScale(m);
            simulate(dt, verbose);
            l = getSummed("length");
            inc_ = l - ol;
            pop();
            if (verbose) {
                std::cout << "\t(sl, mid, sr) = (" << sl << ", " <<  m << ", " <<  sr << "), inc " <<  inc_ << ", err: " << std::abs(inc_-maxinc) << " > " << accuracy << "\n";
            }
            if (inc_>maxinc) { // concatenate
                sr = m;
            } else {
                sl = m;
            }
            i++;

        }
    }
    this->simulate(dt, verbose);
}

/**
 * Creates a specific tropism from the tropism type index.
 * the function must be extended or overwritten to add more tropisms.
 *
 * @param tt        the tropism type index, given in the root type parameters
 * @param N         tropism parameter (passsed to the tropism class)
 * @param sigma     tropism parameter (passsed to the tropism class)
 * @return          the tropism class containing with the callback functions
 */
Tropism* RootSystem::createTropismFunction(int tt, double N, double sigma) {
    // std::cout << "Creating (" << tt << ", " << N << ", " << sigma <<")\n";
    switch (tt) {
    case tt_plagio: return new Plagiotropism(this,N,sigma);
    case tt_gravi: return new Gravitropism(this,N,sigma);
    case tt_exo: return new Exotropism(this,N,sigma);
    case tt_hydro: {
        Tropism* gt =  new Gravitropism(this,N,sigma);
        Tropism* ht= new Hydrotropism(this,N,sigma,soil);
        Tropism* cht = new CombinedTropism(this,N,sigma,ht,10.,gt,1.); // does only use the objective functions from gravitropism and hydrotropism
        return cht;
    }
    default: throw std::invalid_argument( "RootSystem::createTropismFunction() tropism type not implemented" );
    }
}

/**
 * Creates a growth functions from the growth function index.
 * This function must bee extended or overwritten to add more growth functions
 *
 * @param gft       the growth function index, given in the root type parameters
 * @return          the growth function class containing with the callback functions
 */
GrowthFunction* RootSystem::createGrowthFunction(int gft) {
    switch (gft) {
    case gft_negexp: return new ExponentialGrowth();
    case gft_linear: return new LinearGrowth();
    default: throw std::invalid_argument( "RootSystem::createGrowthFunction() growth function type not implemented" );
    }
}

/**
 * Represents the root system as sequential vector of roots, copies the root only, if it has more than 1 node.
 * buffers the result, until next call of simulate(dt), a bit redundant to @see Organsim::getOrgans().
 *
 * \return sequential vector of roots with more than 1 node
 */
std::vector<Root*> RootSystem::getRoots() const
{
    if (roots.empty()) { // create buffer
        std::vector<Organ*> organs;
        organs.reserve(getNumberOfOrgans()); // for speed up
        for (const auto& br : this->baseOrgans) {
            ((Root*)br)->getOrgans(ot_root, organs);
        }
        for (auto& o :organs) {
            roots.push_back((Root*)o);
        }
        return roots;
    } else { // return buffer
        return roots;
    }
}

/**
 * @copydoc Organism::getNodes
 *
 * adds an artificial shoot node
 */
std::vector<Vector3d> RootSystem::getNodes() const
{
    auto v = Organism::getNodes();
    v.at(0) = Vector3d(0.,0.,0.); // add artificial shoot
    return v;
}


/**
 * @return the node indices of the root tips
 */
std::vector<int> RootSystem::getRootTips() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<int> tips;
    for (auto& r : roots) {
        tips.push_back(r->getNodeId(r->getNumberOfNodes()-1));
    }
    return tips;
}

/**
 * @return the node indices of the root bases
 */
std::vector<int> RootSystem::getRootBases() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<int> bases;
    for (auto& r : roots) {
        bases.push_back(r->getNodeId(0));
    }
    return bases;
}

/**
 * Return the segments connecting tap root, basal roots, and shoot borne roots.
 *
 * The upper node represents the oldest emerged shoot-borne root, or if none, the node where
 * all basal roots emerge.
 */
std::vector<Vector2i> RootSystem::getShootSegments() const
{
    std::vector<Vector2i> seg = std::vector<Vector2i>(0);
    seg.push_back(Vector2i(0,baseOrgans.at(0)->getNodeId(0))); // connect  basal roots node (seed) to artificial shoot
    int n1=0, n2=0;
    for (int i=0; i<numberOfCrowns-1; i++) { // connecting root crowns
        int brn = baseOrgans.size()-1;
        n1 = baseOrgans.at(brn-i*seedParam.nC)->getNodeId(0);
        n2 = baseOrgans.at(brn-(i+1)*seedParam.nC)->getNodeId(0);
        seg.push_back(Vector2i(n1,n2));
    }
    if (numberOfCrowns>0) { // connect to basal roots node (seed)
        int ti = baseOrgans.at(0)->getNodeId(0);
        seg.push_back(Vector2i(n2,ti));
    }
    return seg;
}

/**
 * Pushes current root system state to the stack
 */
void RootSystem::push()
{
    stateStack.push(RootSystemState(*this));
}

/**
 * Retrieves previous root system state from the stack
 */
void RootSystem::pop()
{
    RootSystemState& rss = stateStack.top();
    rss.restore(*this);
    stateStack.pop();
}

/**
 * @return quick info about the root system for debugging
 */
std::string RootSystem::toString() const
{
    std::stringstream str;
    str << "Rootsystem with "<< baseOrgans.size() <<" base roots, " << getNumberOfNodes()
                            << " nodes, and a total of " << getNumberOfOrgans() << " roots, after " << getSimTime() << " days";
    return str.str();
}

/**
 * Exports the simulation results with the type from the extension in name
 * (that must be lower case)
 *
 * @param name      file name e.g. output.vtp
 */
void RootSystem::write(std::string name) const
{
    std::ofstream fos;
    fos.open(name.c_str());
    std::string ext = name.substr(name.size()-3,name.size()); // pick the right writer
    if (ext.compare("sml")==0) {
        std::cout << "writing RSML... "<< name.c_str() <<"\n";
        //writeRSML(fos);
    } else if (ext.compare("vtp")==0) {
        std::cout << "writing VTP... "<< name.c_str() <<"\n";
        writeVTP(fos);
    } else if (ext.compare(".py")==0)  {
        std::cout << "writing Geometry ... "<< name.c_str() <<"\n";
        writeGeometry(fos);
    } else {
        throw std::invalid_argument("RootSystem::write(): Unkwown file type");
    }
    fos.close();
}

/**
 * Writes current simulation results as VTP (VTK polydata file),
 * where each root is represented by a polyline.
 *
 * Use SegmentAnalyser::writeVTP() for a representation based on segments,
 * e.g. for creating a movie (and run the animate.py script), or mapping values to segments
 *
 * todo use tinyxml2, move to Organism
 *
 * @param os      typically a file out stream
 */
void RootSystem::writeVTP(std::ostream & os) const
{
    this->getRoots(); // update roots (if necessary)
    const auto& nodes = getPolylines();
    const auto& times = getPolylineCTs();

    os << "<?xml version=\"1.0\"?>";
    os << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    os << "<PolyData>\n";
    int non = 0; // number of nodes
    for (const auto& r : roots) {
        non += r->getNumberOfNodes();
    }
    int nol=roots.size(); // number of lines
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
    const size_t N = 3; // SCALARS
    std::string scalarTypeNames[N] = {"type", "order", "radius" };
    for (size_t i=0; i<N; i++) {
        os << "<DataArray type=\"Float32\" Name=\"" << scalarTypeNames[i] <<"\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
        auto scalars = getParameter(scalarTypeNames[i]);
        for (auto s : scalars) {
            os << s<< " ";
        }
        os << "\n</DataArray>\n";
    }
    os << "\n</CellData>\n";
    // POINTS (=nodes)
    os << "<Points>\n"<<"<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\" >\n";
    for (const auto& r : nodes) {
        for (const auto& n : r) {
            os << n.x << " "<< n.y <<" "<< n.z<< " ";
        }
    }
    os << "\n</DataArray>\n"<< "</Points>\n";
    // LINES (polylines)
    os << "<Lines>\n"<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
    int c=0;
    for (const auto& r : roots) {
        for (size_t i=0; i<r->getNumberOfNodes(); i++) {
            os << c << " ";
            c++;
        }
    }
    os << "\n</DataArray>\n"<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
    c = 0;
    for (const auto& r : roots) {
        c += r->getNumberOfNodes();
        os << c << " ";
    }
    os << "\n</DataArray>\n";
    os << "\n</Lines>\n";

    os << "</Piece>\n";
    os << "</PolyData>\n" << "</VTKFile>\n";
}

/**
 * Writes the current confining geometry (e.g. a plant container) as paraview python script
 * Just adds the initial lines, before calling the method of the sdf.
 *
 * @param os      typically a file out stream
 */
void RootSystem::writeGeometry(std::ostream & os) const
{
    os << "from paraview.simple import *\n";
    os << "paraview.simple._DisableFirstRenderCameraReset()\n";
    os << "renderView1 = GetActiveViewOrCreate('RenderView')\n\n";
    geometry->writePVPScript(os);
}



/**
 * Create a root system state object from a rootsystem, use RootSystemState::restore to go back to that state.
 *
 * @param rs        the root system to be stored
 */
RootSystemState::RootSystemState(const RootSystem& rs) : simtime(rs.simtime), rid(rs.organId), old_non(rs.oldNumberOfNodes), old_nor(rs.oldNumberOfOrgans),
    numberOfCrowns(rs.numberOfCrowns), gen(rs.gen), UD(rs.UD), ND(rs.ND)
{
    baseRoots = std::vector<RootState>(rs.baseOrgans.size()); // store base roots
    for (size_t i=0; i<baseRoots.size(); i++) {
        baseRoots[i] = RootState(*((Root*)rs.baseOrgans[i]));
    }
}

/**
 * Restore evolved rootsystem back to its previous state
 *
 * @param rs    the root system to be restored
 */
void RootSystemState::restore(RootSystem& rs)
{
    rs.roots.clear(); // clear buffer
    rs.simtime = simtime; // copy back everything
    rs.organId = rid;
    rs.nodeId = nid;
    rs.oldNumberOfNodes = old_non;
    rs.oldNumberOfOrgans = old_nor;
    rs.numberOfCrowns = numberOfCrowns;
    rs.gen = gen;
    rs.UD = UD;
    rs.ND = ND;
    for (size_t i=0; i<baseRoots.size(); i++) { // restore base roots
        baseRoots[i].restore(*((Root*)rs.baseOrgans[i]));
    }
}

/**
 * Create a root state object from a root, use RootState::restore to go back to that state.
 *
 * @param r        the root to be stored
 */
RootState::RootState(const Root& r): alive(r.alive), active(r.active), age(r.age), length(r.length), old_non(r.oldNumberOfNodes)
{
    lNode = r.nodes.back();
    lNodeId = r.nodeIds.back();
    lneTime = r.nodeCTs.back();
    non = r.nodes.size();
    laterals = std::vector<RootState>(r.children.size());
    for (size_t i=0; i<laterals.size(); i++) {
        laterals[i] = RootState(*((Root*)r.children[i]));
    }
}

/**
 * Restore evolved root back to its previous state
 *
 * @param r    the root to be restored
 */
void RootState::restore(Root& r)
{
    r.alive = alive; // copy things that changed
    r.active = active;
    r.age = age;
    r.length = length;
    r.oldNumberOfNodes = old_non;
    r.nodes.resize(non); // shrink vectors
    r.nodeIds.resize(non);
    r.nodeCTs.resize(non);
    r.nodes.back() = lNode; // restore last value
    r.nodeIds.back() = lNodeId;
    r.nodeCTs.back() = lneTime;
    for (size_t i = laterals.size(); i<r.children.size(); i++) { // delete roots that have not been created
        delete r.children[i];
    }
    r.children.resize(laterals.size()); // shrink and restore laterals
    for (size_t i=0; i<laterals.size(); i++) {
        laterals[i].restore(*((Root*)r.children[i]));
    }
}

} // end namespace CRootBox
