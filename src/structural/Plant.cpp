#include "Plant.h"
#include "RootDelay.h"

#include <memory>
#include <iostream>
#include <sys/stat.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <string>
#include <fstream>

namespace CPlantBox {

/*
 * Constructs plant, initializes random number generator
 * @param seednum    option to set seed (for creation of random number) default = 0.
 */
Plant::Plant(unsigned int seednum): Organism(seednum)
{ }

/**
 * Deep copies the plant
 */
std::shared_ptr<Organism> Plant::copy()
{
    auto no = std::make_shared<Plant>(*this); // default copy constructor
    // std::cout << "instances before increase (copy) " << instances << "\n";
    instances++;
    no->plantId = instances; // add code for instance counting
    // std::cout << "Created Organism (copy): " << no->plantId << " from " << plantId << "\n" << std::flush;
    for (int i=0; i<baseOrgans.size(); i++) {
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
 * todo docme , this could be made unique? and probably should be protected
 */
void Plant::initializeReader()
{
    auto rrp = std::make_shared<RootRandomParameter>(shared_from_this());
    rrp->subType = 0;
    setOrganRandomParameter(rrp);
    auto srp = std::make_shared<SeedRandomParameter>(shared_from_this());
    srp->subType = 0;
    setOrganRandomParameter(srp);
    auto strp = std::make_shared<StemRandomParameter>(shared_from_this());
    strp->subType = 0;
    setOrganRandomParameter(strp);
    // auto strp1 = std::make_shared<StemRandomParameter>(shared_from_this()); // Dummy stem, in case there is no stem defined
    // strp1->subType = 1;
    // setOrganRandomParameter(strp1);
    auto lrp = std::make_shared<LeafRandomParameter>(shared_from_this());
    lrp->subType = 0;
    setOrganRandomParameter(lrp);
}

/**
 * Resets the root system: deletes all roots, sets simulation time to 0.
 */
void Plant::reset()
{
    baseOrgans.clear();
    simtime = 0;
    organId = -1;
    nodeId = -1;
}

/**
 * Sets up the plant according to the plant parameters,
 * a confining geometry, the tropism functions, and the growth functions.
 *
 * If not used for test file: Call this method before simulation and after setting geometry,
 * plant and root parameters
 * @param verbose       print information
 */
void Plant::initialize_(bool verbose)
{
    oldNumberOfNodes = getNumberOfNodes(); // todo check what this does

    // further initializations
	initCallbacks();
}

/**
 * Sets up the plant according to the plant parameters,
 * a confining geometry, the tropism functions, and the growth functions.
 *
 * LB, Length based: Delay for lateral root is calculated from the apical length (classical RootBox approach)
 *
 * Call this method before simulation and after setting geometry,
 * plant and root parameters
 * @param verbose       print information
 */
void Plant::initializeLB(bool verbose)
{
    reset(); // just in case
    auto seed = std::make_shared<Seed>(shared_from_this());
    this->addOrgan(seed);
	seed->initialize(verbose);
    initialize_(verbose);
}

/**
 * Sets up the plant according to the plant parameters,
 * a confining geometry, the tropism functions, and the growth functions.
 *
 * DB, Delay based: Delay for lateral root is predefined, apical length therefore not constant
 *
 * Call this method before simulation and after setting geometry,
 * plant and root parameters
 * @param verbose       print information
 */
void Plant::initializeDB(bool verbose)
{
	reset(); // just in case

    class SeedDB :public Seed { // make the seed use the RootDelay class
    	using Seed::Seed;
    	std::shared_ptr<Organ> createRoot(std::shared_ptr<Organism> plant, int type,  double delay) override {
    		return std::make_shared<RootDelay>(plant, type, delay, shared_from_this(), 0);
    	};
    };

    auto seed = std::make_shared<SeedDB>(shared_from_this());
    this->addOrgan(seed);
	seed->initialize(verbose);
    initialize_(verbose);
}

/**
 * Called by Plant::initialize.
 * Sets up tropism and growth functions call backs using
 * Plant::createTropismFunction and Plant::createGrowthFunction
 */
void Plant::initCallbacks()
{
    // Create tropisms and growth functions per random root parameter
    for (auto& p_otp :organParam[Organism::ot_root]) {
		auto rp = std::static_pointer_cast<RootRandomParameter>(p_otp.second);
        auto tropism = this->createTropismFunction(rp->tropismT, rp->tropismN, rp->tropismS);
        tropism->setGeometry(geometry);
        rp->f_tf = tropism; // set new one
        auto gf_ = this->createGrowthFunction(rp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        rp->f_gf  = gf_;
    }
    // Create tropisms and growth functions per random leaf parameter
    for (auto& p_otp :organParam[Organism::ot_leaf]) {
		auto rp = std::static_pointer_cast<LeafRandomParameter>(p_otp.second);
		double Tage =  rp->tropismAge +  rp->tropismAges * randn();
        auto tropism = this->createTropismFunction(rp->tropismT, rp->tropismN, rp->tropismS, Tage);
        tropism->setGeometry(geometry);
        rp->f_tf = tropism; // set new one
        auto gf_ = this->createGrowthFunction(rp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        rp->f_gf  = gf_;
    }
    // Create tropisms and growth functions per random stem parameter
    for (auto& p_otp :organParam[Organism::ot_stem]) {
		auto rp = std::static_pointer_cast<StemRandomParameter>(p_otp.second);
		double Tage =  rp->tropismAge +  rp->tropismAges * randn();
        auto tropism = this->createTropismFunction(rp->tropismT, rp->tropismN, rp->tropismS, Tage);
        tropism->setGeometry(geometry);
        rp->f_tf = tropism; // set new one
        auto gf_ = this->createGrowthFunction(rp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        rp->f_gf  = gf_;
    }

}

/**
 * Manually sets a tropism function for a specific or for all root types. TODO
 * Must be called after RootSystem::initialize(), otherwise its overwritten.
 *
 * @param tf_           a tropism
 * @param organType     organ type, if organType = -1 all types are set to this tropism (default).
 * @param subtype       organ subtype
 */
void Plant::setTropism(std::shared_ptr<Tropism> tf, int organType, int subType) // todo
{
	if (organType==-1) { // call for all relevant
        setTropism(tf, Organism::ot_root, subType);
        setTropism(tf, Organism::ot_stem, subType);
        setTropism(tf, Organism::ot_leaf, subType);
        return;
    }
    // put it into a vector
    std::vector<std::shared_ptr<OrganRandomParameter>> orp;
    if (subType>-1) { // set for a specific root type
        orp.push_back(getOrganRandomParameter(organType, subType));
    } else { // set for all root types (default)
        for (auto& orp_ :organParam[organType]) {
            orp.push_back(orp_.second);
        }
    }
    // static cast since tropisms are not in the base class
    for (auto& orp_ : orp)
    switch(organType) {
    case Organism::ot_root: {
        auto rtp = std::static_pointer_cast<RootRandomParameter>(orp_);
        rtp->f_tf = tf;
    } break;
    case Organism::ot_stem: {
        auto rtp = std::static_pointer_cast<StemRandomParameter>(orp_);
        rtp->f_tf = tf;
    } break;
    case Organism::ot_leaf: {
        auto rtp = std::static_pointer_cast<LeafRandomParameter>(orp_);
        rtp->f_tf = tf;
    } break;
    default: throw std::invalid_argument( "Plant::setTropism() cannot set tropism for organ type "+ Organism::organTypeName(organType) );
    }
}

/**
 * Simulates plant growth
 * @param dt		duration of the simulation
 * @param verbose	whether to print information
 */
void Plant::simulate(double dt, bool verbose)
{
	abs2rel();
    Organism::simulate(dt, verbose);
	rel2abs();
}

/**
 * Simulates plant growth for the time span defined in the root system parameters
 */
void Plant::simulate()
{
    auto srp = std::static_pointer_cast<SeedRandomParameter>(organParam[Organism::ot_seed][0]);
    Plant::simulate(srp->simtime);
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
void Plant::simulate(double dt, double maxinc_, std::shared_ptr<ProportionalElongation> se, bool verbose)
{
    this->simulateLimited(dt, maxinc_, "lengthTh", {1.,1.,1.,1.,1.}, se, verbose);
}

/**
 * Simulates root system growth for a time span, elongates a maximum of @param maxinc total length [cm/day]
 * using the proportional elongation @param se to impede overall growth.
 *
 * @param dt            time step [day]
 * @param maxinc_       maximal paramName [(paramName units)/day] the root system is allowed to grow in this time step
 * @param paramName     e.g. length or volume
 * @param scales        per organ type { ot_organ = 0, ot_seed = 1, ot_root = 2, ot_stem = 3, ot_leaf = 4 };
 * @param se            The class ProportionalElongation is used to scale overall root growth
 * @param verbose       indicates if status is written to the console (cout) (default = false)
 */
void Plant::simulateLimited(double dt, double max_, std::string paramName, std::vector<double> scales,
    std::shared_ptr<ProportionalElongation> se, bool verbose)
{
    const double accuracy = 1.e-3;
    const int maxiter = 20;
    double maxinc = dt*max_; // [cm]

    double ol = this->weightedSum(paramName, scales);
    int i = 0;

    // test run with scale == 1 (on copy)
    std::shared_ptr<Plant> rs = std::static_pointer_cast<Plant>(this->copy());
    se->setScale(1.);
    rs->simulate(dt, verbose);
    double l = rs->weightedSum(paramName, scales);
    double inc_ = l - ol;
    if (verbose) {
        std::cout << "expected increase is " << inc_ << " maximum is " << maxinc << "\n";
    }

    if ((inc_>maxinc) && (std::abs(inc_-maxinc)>accuracy)) { // if necessary, perform a binary search

        double sl = 0.; // left
        double sr = 1.; // right

        while ( ((std::abs(inc_-maxinc)) > accuracy) && (i<maxiter) )  { // binary search

            double m = (sl+sr)/2.; // mid

            // test run (on copy) with scale m
            std::shared_ptr<Plant> rs = std::static_pointer_cast<Plant>(this->copy()); // reset to old #############
            se->setScale(m);
            rs->simulate(dt, false);
            l = rs->weightedSum(paramName, scales);
            inc_ = l - ol;

            if (verbose) {
                std::cout <<  i << "\t(sl, mid, sr) = (" << sl << ", " <<  m << ", " <<  sr << "), inc " <<  inc_ << ", err: " << std::abs(inc_-maxinc) << "<>" << accuracy << "\n";
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

double Plant::weightedSum(std::string paramName, std::vector<double> scales) const {
    double sum=0.;
    auto organs = this->getOrgans();
    for (const auto& o : organs) {
        int i  = o->organType();
        sum += scales.at(i)*o->getParameter(paramName);
    }
    return sum;
}


/**
 * go from absolute to relative coordinates for aboveground organs
 */
void Plant::abs2rel()
{
	auto s = getSeed();
	for (int i = 0; i< s->getNumberOfChildren();i++) {
		auto child = s->getChild(i);
		//if(child->organType() >2){ //if aboveground-organ
			child->abs2rel(); //apply to all organs
		//}

    }

}

/**
 * go from relative to absolute coordinates for aboveground organs
 */
void Plant::rel2abs()
{
	auto s = getSeed();
	for (int i = 0; i< s->getNumberOfChildren();i++) {
		auto child = s->getChild(i);
		//if(child->organType() >2){ //if aboveground-organ
			child->rel2abs();//apply to all organs
		//}

    }
}

/**
 * Creates a specific tropism from the tropism type index.
 * the function must be extended or overwritten to add more tropisms.
 *
 * @param tt        the tropism type index, given in the root type parameters
 * @param N         tropism parameter (passsed to the tropism class)
 * @param sigma     tropism parameter (passsed to the tropism class)
 * @param ageSwitch age at which new tropism funciton is implemented
 * @return          the tropism class containing with the callback functions
 */
std::shared_ptr<Tropism> Plant::createTropismFunction(int tt, double N, double sigma, double ageSwitch) {
    switch (tt) {
    case tt_plagio: return std::make_shared<Plagiotropism>(shared_from_this(),N,sigma);
    case tt_gravi: return std::make_shared<Gravitropism>(shared_from_this(),N,sigma);
    case tt_exo: return std::make_shared<Exotropism>(shared_from_this(),N,sigma);
    case tt_hydro: { // uses weighted objective functions from gravitropism and hydrotropism
        auto gt =  std::make_shared<Gravitropism>(shared_from_this(),N,sigma);
        auto ht= std::make_shared<Hydrotropism>(shared_from_this(),N, sigma, soil);
        return std::make_shared<CombinedTropism>(shared_from_this(),N,sigma,ht,10.,gt,1.);
    }
    case tt_twist:  return std::make_shared<TwistTropism>(shared_from_this(),N,sigma);
    case tt_antigravi: return std::make_shared<AntiGravitropism>(shared_from_this(),N,sigma);
	case tt_antigravi2gravi: return std::make_shared<AntiGravi2Gravitropism>(shared_from_this(),N,sigma, ageSwitch);
    default: throw std::invalid_argument( "Plant::createTropismFunction() tropism type not implemented" );
    }
}

/**
 * Creates a growth functions from the growth function index.
 * This function must bee extended or overwritten to add more growth functions
 *
 * @param gft       the growth function index, given in the root type parameters
 * @return          the growth function class containing with the callback functions
 */
std::shared_ptr<GrowthFunction>Plant::createGrowthFunction(int gft) {
    switch (gft) {
    case gft_negexp: return std::make_shared<ExponentialGrowth>();
    case gft_linear: return std::make_shared<LinearGrowth>();
    case gft_CWLim: return std::make_shared<CWLimitedGrowth>();
    default: throw std::invalid_argument( "Plant::createGrowthFunction() growth function type not implemented" );
    }
}

/**
 * todo most important debug informations
 */
std::string Plant::toString() const
{
    return "Plant " + Organism::toString();
}

} // namespace CPlantBox
