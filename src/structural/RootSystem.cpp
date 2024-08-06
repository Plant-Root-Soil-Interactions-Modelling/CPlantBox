// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "RootSystem.h"

#include "organparameter.h"
#include "Organism.h"
#include "RootDelay.h"

namespace CPlantBox {

/**
 * Creates a root system
 */
RootSystem::RootSystem(unsigned int seednum): Organism(seednum)
{ }

/**
 * Deep copies the organism
 */
std::shared_ptr<Organism> RootSystem::copy()
{
    roots.clear(); // clear buffer
    auto nrs = std::make_shared<RootSystem>(*this); // copy constructor
    nrs->seed = std::static_pointer_cast<Seed>(seed->copy(nrs));
    baseOrgans = nrs->seed->copyBaseOrgans();
    for (int ot = 0; ot < numberOfOrganTypes; ot++) { // copy organ type parameters
        for (auto& otp : nrs->organParam[ot]) {
            otp.second = otp.second->copy(nrs);
			nrs->setOrganRandomParameter(otp.second);
        }
    }
    return nrs;
}

/**
 * @return the i-th root parameter of sub type @param type.
 */
std::shared_ptr<RootRandomParameter> RootSystem::getRootRandomParameter(int type) const
{
    return std::static_pointer_cast<RootRandomParameter>(getOrganRandomParameter(Organism::ot_root, type));
}

/**
 * @return all root type parameters in a vector
 */
std::vector<std::shared_ptr<RootRandomParameter>> RootSystem::getRootRandomParameter() const
{
    auto otps = std::vector<std::shared_ptr<RootRandomParameter>>(0);
    for (auto& otp : organParam[Organism::ot_root]) {
        otps.push_back(std::static_pointer_cast<RootRandomParameter>(otp.second));
    }
    return otps;
}

/**
 * Sets the root system parameters @param rsp
 */
void RootSystem::setRootSystemParameter(std::shared_ptr<SeedRandomParameter> sp)
{
    assert(sp->subType==0 && "RootSystem::setRootSystemParameter: In CPlantBox must have subType 0");
    this->setOrganRandomParameter(sp);
}

/**
 * @return the root system parameter
 */
std::shared_ptr<SeedRandomParameter> RootSystem::getRootSystemParameter()
{
    return std::static_pointer_cast<SeedRandomParameter>(this->getOrganRandomParameter(ot_seed,0 ));
}

/**
 * Resets the root system: deletes all roots, sets simulation time to 0.
 */
void RootSystem::reset()
{
    roots.clear(); // clear buffer
    baseOrgans.clear();
    simtime = 0;
    organId = -1;
    nodeId = -1;
}

/**
 * Initializes the xml parameter reader, specifying prototypes for roots (RootRandomParameter), and seed (SeedRandomParameter)
 *
 * Overwrite this function if you want to use other parameter classes.
 */
void RootSystem::initializeReader()
{
    auto rrp = std::make_shared<RootRandomParameter>(shared_from_this());
    rrp->subType = 0;
    setOrganRandomParameter(rrp);
    auto srp = std::make_shared<SeedRandomParameter>(shared_from_this());
    srp->subType = 0;
    setOrganRandomParameter(srp);
}

/**
 * DEPRICATED Reads the root parameter from a file. Opens plant parameters with the same filename if available,
 * otherwise assumes a tap root system at position (0,0,-3).
 *
 * @param name          filename without file extension
 * @param subdir        directory ("modelparameter/" by default)
 */
void RootSystem::openFile(std::string name, std::string subdir)
{
    std::cout << "RootSystem::openFile is deprecated, use readParameters instead \n";
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
    auto randomSeed = std::make_shared<SeedRandomParameter>(shared_from_this());
    randomSeed->subType = 0;
    if (fis.good()) { // did it work?
        randomSeed->read(fis); // reads the random parameters
        fis.close();
    } else { // create a tap root system
        std::cout << "No root system parameters found, using default tap root system \n";
    }
    this->setRootSystemParameter(randomSeed);
}

/**
 * DEPRICATED Reads root type parameter from an input stream @param is
 * (there is a Matlab script exporting these, @see writeParams.m)
 *
 * @param cin  in stream
 */
int RootSystem::readParameters(std::istream& is)
{
    std::cout << "RootSystem::readParameters(std::istream) is deprecated, use readParameters(std::string filename) instead \n";
    int c = 0;
    while (is.good()) {
        auto p = std::make_shared<RootRandomParameter>(shared_from_this());
        p->read(is);
        p->organType = Organism::ot_root;
        setOrganRandomParameter(p);
        c++;
    }
    return c;
}

/**
 * (DEPRICATED) Writes root type parameters to an output stream @param os
 *
 * @param os  out stream
 */
void RootSystem::writeParameters(std::ostream& os) const
{
    std::cout << "RootSystem::writeParameters is deprecated, use writeParameters(std::string filename) instead \n";
    for (auto& otp :organParam[Organism::ot_root]) {
        std::static_pointer_cast<RootRandomParameter>(otp.second)->write(os);
    }
}

/**
 * Sets up the base roots according to the plant parameters,
 * a confining geometry, the tropism functions, and the growth functions.
 *
 * LB, Length based: Delay for lateral root is calculated from the apical length (classical RootBox approach)
 *
 * Call this method before simulation and after setting geometry, plant and root parameters
 *
 * @parm basal  	    the type of the basal roots (default = 4)
 * @parm shootborne     the type of the shootborne roots (default = 5)
 * @param verbose 	 	chatty with the std::couts
 */
void RootSystem::initializeLB(int basal, int shootborne, bool verbose)
{
	reset(); // just in case
    getNodeIndex(); // introduce an extra node at nodes[0]
    seed = std::make_shared<Seed>(shared_from_this()); // introduce a 2nd node =>  2 nodes to make a seed segment
	initialize_(basal, shootborne, verbose);
}

/**
 * Sets up the base roots according to the plant parameters,
 * a confining geometry, the tropism functions, and the growth functions.
 *
 * DB, Delay based: Delay for lateral root is predefined, apical length therefore not constant
 *
 * Call this method before simulation and after setting geometry, plant and root parameters
 *
 * @parm basal      	the type of the basal roots (default = 4)
 * @parm shootborne     the type of the shootborne roots (default = 5)
 * @param verbose 	 	chatty with the std::couts
 */
void RootSystem::initializeDB(int basal, int shootborne, bool verbose)
{
	reset(); // just in case
    getNodeIndex(); // introduce an extra node used with node created by seed to have a seed segment.

    class SeedDB :public Seed { // make the seed use the RootDelay class
    	using Seed::Seed;
    	std::shared_ptr<Organ> createRoot(std::shared_ptr<Organism> plant, int type, double delay) override {
    		return std::make_shared<RootDelay>(plant, type,  delay, shared_from_this(), 0);
    	};
    };

    seed = std::make_shared<SeedDB>(shared_from_this());
    initialize_(basal, shootborne, verbose);
}

/**
 * Initializes the seed (@see initialize, initializeDB, initializeLB)
 */
void RootSystem::initialize_(int basal, int shootborne, bool verbose)
{
    seed->basalType = basal;
    seed->shootborneType = shootborne;
    seed->initialize(verbose);
    seedParam = SeedSpecificParameter(*seed->param()); // copy the specific parameters
    baseOrgans = seed->copyBaseOrgans();
    numberOfCrowns = seed->getNumberOfRootCrowns(); // a bit redundant...
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
        auto rp = std::static_pointer_cast<RootRandomParameter>(p_otp.second);
        auto tropism = this->createTropismFunction(rp->tropismT, rp->tropismN, rp->tropismS);
        tropism->setGeometry(geometry);
        rp->f_tf = tropism; // set new one
        auto gf_ = this->createGrowthFunction(rp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        rp->f_gf  = gf_;
    }
}

/**
 * Manually sets a tropism function for a specific or for all root types.
 * Must be called after RootSystem::initialize(), otherwise its overwritten.
 *
 * @param tf_           a tropism
 * @param rt            root type, if rt = -1 all types are set to this tropism (default).
 */
void RootSystem::setTropism(std::shared_ptr<Tropism> tf_, int rt)
{
    if (rt>-1) { // set for a specific root type
        getRootRandomParameter(rt)->f_tf=tf_;
    } else { // set for all root types (default)
        for (auto& p_otp :organParam[Organism::ot_root]) {
            auto rtp = std::static_pointer_cast<RootRandomParameter>(p_otp.second);
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
std::shared_ptr<Tropism>  RootSystem::createTropismFunction(int tt, double N, double sigma) {
    // std::cout << "Creating (" << tt << ", " << N << ", " << sigma <<")\n";
    switch (tt) {
    case tt_plagio: return std::make_shared<Plagiotropism>(shared_from_this(),N,sigma);
    case tt_gravi: return std::make_shared<Gravitropism>(shared_from_this(),N,sigma);
    case tt_exo: return std::make_shared<Exotropism>(shared_from_this(),N,sigma);
    case tt_hydro: { // uses weighted objective functions from gravitropism and hydrotropism
        auto gt =  std::make_shared<Gravitropism>(shared_from_this(),N,sigma);
        auto ht= std::make_shared<Hydrotropism>(shared_from_this(),N, sigma, soil);
        return std::make_shared<CombinedTropism>(shared_from_this(),N,sigma,ht,10.,gt,1.);
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
std::shared_ptr<GrowthFunction>  RootSystem::createGrowthFunction(int gft) {
    switch (gft) {
    case gft_negexp: return std::make_shared<ExponentialGrowth>();
    case gft_linear: return std::make_shared<LinearGrowth>();
    default: throw std::invalid_argument( "RootSystem::createGrowthFunction() growth function type not implemented" );
    }
}

/**
 * Represents the root system as sequential vector of roots, copies the root only, if it has more than 1 node.
 * buffers the result, until next call of simulate(dt), a bit redundant to @see Organsim::getOrgans().
 *
 * \return sequential vector of roots with more than 1 node
 */
std::vector<std::shared_ptr<Root>> RootSystem::getRoots() const
{
    if (roots.empty()) { // create buffer
        std::vector<std::shared_ptr<Organ>> organs;
        organs.reserve(getNumberOfOrgans()); // for speed up
        for (const auto& br : this->baseOrgans) {
            br->getOrgans(ot_root, organs);
        }
        for (auto& o :organs) {
            roots.push_back(std::static_pointer_cast<Root>(o));
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
    v.at(0) = Vector3d(0.,0.,0.); // the artifical node is created by initialize()
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
    int n1=0, n2=0;
    if (numberOfCrowns>0) { // connect to basal roots node (seed)
        seg.push_back(Vector2i(0,baseOrgans.back()->getNodeId(0)));
    }
    for (int i=0; i<numberOfCrowns-1; i++) { // connecting root crowns
        int brn = baseOrgans.size()-1;
        n1 = baseOrgans.at(brn-i*seedParam.nC)->getNodeId(0);
        n2 = baseOrgans.at(brn-(i+1)*seedParam.nC)->getNodeId(0);
        seg.push_back(Vector2i(n1,n2));
    }
    seg.push_back(Vector2i(n2,baseOrgans.at(0)->getNodeId(0))); // connect  basal roots node (seed) to artificial shoot
	std::cout<<"RootSystem::getShootSegments "<<numberOfCrowns<<" "<<seg.size()<<std::endl;
    return seg;
}

/**
 * Pushes current root system state to the stack
 */
void RootSystem::push()
{
    stateStack.push_back(RootSystemState(*this));
}

/**
 * Retrieves previous root system state from the stack
 */
void RootSystem::pop()
{
    RootSystemState& rss = stateStack.back();
    rss.restore(*this);
    stateStack.pop_back();
}

/**
 * @return quick info about the root system for debugging
 */
std::string RootSystem::toString() const
{
    std::stringstream str;
    str << "Rootsystem with "<< baseOrgans.size() <<" base roots, " << getNumberOfNodes()
                                << " nodes, and a total of " << getNumberOfOrgans() << " organs, after " << getSimTime() << " days";
    return str.str();
}

/**
 * Create a root system state object from a rootsystem, use RootSystemState::restore to go back to that state.
 *
 * @param rs        the root system to be stored
 */
RootSystemState::RootSystemState(const RootSystem& rs) : simtime(rs.simtime), dt(rs.dt), organId(rs.organId), nodeId(rs.nodeId),
    oldNumberOfOrgans(rs.oldNumberOfOrgans), numberOfCrowns(rs.numberOfCrowns), gen(rs.gen), UD(rs.UD), ND(rs.ND)
{
    baseRoots = std::vector<RootState>(rs.baseOrgans.size()); // store base roots
    for (size_t i=0; i<baseRoots.size(); i++) {
        baseRoots[i] = RootState(*(std::static_pointer_cast<Root>(rs.baseOrgans[i])));
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
    rs.dt = dt;
    rs.organId = organId;
    rs.nodeId = nodeId;
    rs.oldNumberOfNodes = oldNumberOfNodes;
    rs.oldNumberOfOrgans = oldNumberOfOrgans;
    rs.numberOfCrowns = numberOfCrowns;

    rs.gen = gen;
    rs.UD = UD;
    rs.ND = ND;
    for (size_t i=0; i<baseRoots.size(); i++) { // restore base roots
        baseRoots[i].restore(*(std::static_pointer_cast<Root>(rs.baseOrgans[i])));
    }
}

/**
 * Create a root state object from a root, use RootState::restore to go back to that state.
 *
 * @param r        the root to be stored
 */
RootState::RootState(const Root& r): alive(r.alive), active(r.active), age(r.age), length(r.getLength(true)),
    epsilonDx(r.epsilonDx), moved(r.moved), oldNumberOfNodes(r.oldNumberOfNodes), firstCall(r.firstCall)
{
    lNode = r.nodes.back();
    lNodeId = r.nodeIds.back();
    lneTime = r.nodeCTs.back();
    non = r.nodes.size();
    laterals = std::vector<RootState>(r.children.size());
    for (size_t i=0; i<laterals.size(); i++) {
        laterals[i] = RootState(*(std::static_pointer_cast<Root>(r.children[i])));
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
    r.epsilonDx = epsilonDx;
    r.moved = moved;
    r.oldNumberOfNodes = oldNumberOfNodes;
    r.firstCall = firstCall; //

    r.nodes.resize(non); // shrink vectors
    r.nodeIds.resize(non);
    r.nodeCTs.resize(non);
    r.nodes.back() = lNode; // restore last value
    r.nodeIds.back() = lNodeId;
    r.nodeCTs.back() = lneTime;
    r.children.resize(laterals.size()); // shrink and restore laterals
    for (size_t i=0; i<laterals.size(); i++) {
        laterals[i].restore(*(std::static_pointer_cast<Root>(r.children[i])));
    }
}

} // end namespace CPlantBox
