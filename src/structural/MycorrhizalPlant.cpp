#include "Plant.h"
#include "Organ.h"
#include "Organism.h"
#include "MycorrhizalPlant.h"
#include "MycorrhizalRoot.h"
#include "mycorrhizalrootparameter.h"
#include "hyphaeparameter.h"
#include "sdf.h"
#include "soil.h"



namespace CPlantBox {
MycorrhizalPlant::MycorrhizalPlant(unsigned int seednum): Plant(seednum) {}

void MycorrhizalPlant::initializeReader()
{
    //new Parameters start here
    auto mycrrp = std::make_shared<MycorrhizalRootRandomParameter>(shared_from_this());
    mycrrp -> subType = 0;
    setOrganRandomParameter(mycrrp);
    //new Parameters end here
    auto srp = std::make_shared<SeedRandomParameter>(shared_from_this());
    srp->subType = 0;
    setOrganRandomParameter(srp);
    auto strp = std::make_shared<StemRandomParameter>(shared_from_this());
    strp->subType = 0;
    setOrganRandomParameter(strp);
    auto strp1 = std::make_shared<StemRandomParameter>(shared_from_this()); // Dummy stem, in case there is no stem defined
    strp1->subType = 1;
    setOrganRandomParameter(strp1);
    auto lrp = std::make_shared<LeafRandomParameter>(shared_from_this());
    lrp->subType = 0;
    setOrganRandomParameter(lrp);
}

void MycorrhizalPlant::initializeLB(bool verbose)
{
    reset(); // just in case
    class MycorrhizalSeed :public Seed{
        using Seed::Seed;
        std::shared_ptr<Organ> createRoot(std::shared_ptr<Organism> plant, int type, double delay) override {
            return std::make_shared<MycorrhizalRoot>(plant, type, delay, shared_from_this(), 0);
        };
    };
    auto seed = std::make_shared<MycorrhizalSeed>(shared_from_this());
    baseOrgans.push_back(seed);
    seed->initialize(verbose);
    initialize_(verbose);
}


std::vector<int> MycorrhizalPlant::getNodeInfections(int ot) const {
    auto organs = this -> getOrgans(ot);
    std::vector<int> infs = std::vector<int>(getNumberOfNodes());
    for (const auto& o : baseOrgans)
    {
        if(o->organType() == Organism::ot_root)
        {

            infs.at(o->getNodeId(0)) = std::dynamic_pointer_cast<MycorrhizalRoot> (o) -> getNodeInfection(0);
        }
    }

    for (const auto & o : organs)
    {
        for (size_t i = 1; i < o ->getNumberOfNodes()-1; i++)
        {
            infs.at(o->getNodeId(i)) = std::dynamic_pointer_cast<MycorrhizalRoot> (o) -> getNodeInfection(i);
        }

    }
    return infs;
}

std::vector<double> MycorrhizalPlant::getNodeInfectionTime(int ot) const {
    auto organs = this -> getOrgans(ot);
    std::vector<double> infTime = std::vector<double>(getNumberOfNodes());
    for (const auto& o : baseOrgans)
    {
        if(o->organType() == Organism::ot_root)
        {

            infTime.at(o->getNodeId(0)) = std::dynamic_pointer_cast<MycorrhizalRoot> (o) -> getNodeInfectionTime(0);
        }
    }

    for (const auto & o : organs)
    {
        for (size_t i = 1; i < o ->getNumberOfNodes()-1; i++)
        {
            infTime.at(o->getNodeId(i)) = std::dynamic_pointer_cast<MycorrhizalRoot> (o) -> getNodeInfectionTime(i);
        }

    }
    return infTime;
}
/**
 * Simulates plant growth
 * @param dt		duration of the simulation
 * @param verbose	whether to print information
 */
void MycorrhizalPlant::simulate(double dt, bool verbose)
{
    auto organs = getOrgans();
    abs2rel();
    Organism::simulate(dt, verbose);
    rel2abs();
}

void MycorrhizalPlant::simulatePrimaryInfection(double dt, bool verbose) {
    for (const auto& r : baseOrgans) {
        std::dynamic_pointer_cast<MycorrhizalRoot>(r)->simulatePrimaryInfection(dt);
    }
}
void MycorrhizalPlant::simulateSecondaryInfection(double dt, bool verbose) {
    for (const auto& r : baseOrgans) {
        std::dynamic_pointer_cast<MycorrhizalRoot>(r)->simulateSecondaryInfection(dt);
    }
}

void MycorrhizalPlant::simulateHyphalGrowth(double dt)
{
    oldNumberOfOrgans = getNumberOfOrgans();
    auto organs = getOrgans(); // TODO  getOrgans(Organism::ot_root) empty list????
    for (const auto & o : organs) {
        if (o->organType() == Organism::ot_root) {
            std::dynamic_pointer_cast<MycorrhizalRoot>(o) -> simulateHyphalGrowth();
        }
    }
};

void MycorrhizalPlant::initCallbacks() {

    // std::cout << "MycorrhizalPlant::initCallbacks()\n";

    Plant::initCallbacks(); 

    for (auto& p_otp :organParam[Organism::ot_root]) {
        auto rp = std::static_pointer_cast<MycorrhizalRootRandomParameter>(p_otp.second);
        auto bigbox = std::make_shared<SDF_PlantBox>(1.e100,1.e100,1.e100); // TODO Fix this
        auto inf_ = std::make_shared<SoilLookUpSDF>();
        auto sdf = SignedDistanceFunction();
        // auto sdf = SoilLookUp();
        sdf = *bigbox;
        inf_->sdf = &sdf;
        inf_->fmax = rp->p;
        inf_->fmin = 0;
        inf_->slope = 0;
        // std::cout<< inf_->sdf->toString() << std::endl;
        // rp->f_inf  = sdf; // set new one
    }

    // Create tropisms and growth functions per random hyphae parameter
    for (auto& p_otp :organParam[Organism::ot_hyphae]) {
        auto rp = std::static_pointer_cast<HyphaeRandomParameter>(p_otp.second);
        auto tropism = this->createTropismFunction(rp->tropismT, rp->tropismN, rp->tropismS);
        tropism->setGeometry(geometry);
        rp->f_tf = tropism; // set new one
        // growth function is set to LinearGrowth in constructor of HyphaeRandomParameter
    }

}


// void MycorrhizalPlant::addTree() {
//     // TODO adapt from exudation model or push to hyphae???
//     // right now just straight copy does not make sense
//     dx3 = (length/nx)*(width/ny)*(depth/nz); // for integration of eqn 13
//         roots = rs->getRoots();

//         for (const auto& r : roots) {
//             if (r->getNumberOfNodes()>1) { // started growing
//                 // time when the root stopped growing
//                 double sTime = r->getNodeCT(r->getNumberOfNodes()-1);
//                 if (r->isActive()) {
//                     stopTime.push_back(0);
//                 } else {
//                     stopTime.push_back(sTime);
//                 }
//                 // root tip
//                 Vector3d t = r->getNode(r->getNumberOfNodes()-1);
//                 tip.push_back(t);
//                 // direction towards root base
//                 Vector3d base = r->getNode(0);
//                 double a = r->getNodeCT(r->getNumberOfNodes()-1) - r->getNodeCT(0);
//                 v.push_back(base.minus(t).times(1./a));
//             }

//             sdfs.push_back(SDF_RootSystem(*r, observationRadius));

//         }
// }

}
