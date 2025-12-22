#include "Plant.h"
#include "Organ.h"
#include "Organism.h"
#include "MycorrhizalPlant.h"
#include "MycorrhizalRoot.h"
#include "mycorrhizalrootparameter.h"
#include "hyphaeparameter.h"
// #include "aabbcc/AABB.h"
// #include "sdf.h"
#include "sdf_rs.h"
// #include "soil.h"

#include <functional>

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

std::vector<Vector3d> MycorrhizalPlant::getAnastomosisPoints(int ot) const {
    auto organs = this -> getOrgans(ot);
    std::vector<Vector3d> anaPoints;
    for (const auto& o : baseOrgans)
    {
        if(o->organType() == Organism::ot_hyphae)
        {
            auto h = std::dynamic_pointer_cast<Hyphae> (o);
            if (h->mergePointID != -1) {
                anaPoints.push_back(h->getMergePoint(h->mergePointID));
            }
        }
    }

    for (const auto & o : organs)
    {
        auto h = std::dynamic_pointer_cast<Hyphae> (o);
            if (h->mergePointID != -1) {
                anaPoints.push_back(h->getMergePoint(h->mergePointID));
            }
    }
    return anaPoints;
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
    sdf = std::make_shared<SDF_RootSystem>(*this);
    sdf->selectedOrganType = Organism::ot_hyphae;
    simulateAnastomosis();
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

void MycorrhizalPlant::simulateAnastomosis() {
    auto hyphae = this->getOrgans(Organism::ot_hyphae);
    double dist = 1000;

    for (const auto & h : hyphae) {
        sdf->excludeTreeId = h->getParameter("hyphalTreeIndex");
        // if (std::dynamic_pointer_cast<Hyphae>(h)->mergedHyphae != nullptr) {
        //     std::cout<< std::dynamic_pointer_cast<Hyphae>(h)->mergedHyphae <<std::endl;// already merged
        // }
        // std::cout << "Does this hypha already have merged hyphae? " << (std::dynamic_pointer_cast<Hyphae>(h)->mergedHyphae != nullptr) << "\n";
        if (h->isActive()) {
            auto tip = h->getNode(h->getNumberOfNodes()-1);

            dist = sdf->getDist(tip);

            // std::cout<<"Distance to nearest hyphae from tip " << tip.toString() << " is " << dist << " cm." << std::endl;
            if (fabs(dist) < h->getParameter("distTH")) 
            {
                auto lastIndex = sdf->distIndex; 
                // std::cout<< "Anastomosis occurred at distance: " << dist << " cm.\n";
                // std::cout <<"Anastomosis at tip: " << tip.toString() <<" with distance id: " << lastIndex << std::endl;
                // std::cout << "Hyphal tree index " << h->getParameter("hyphalTreeIndex") << "\n";
                auto connected_to_hyphae = std::dynamic_pointer_cast<Hyphae>(sdf->segO.at(lastIndex).lock());                
                // std::cout << "connected to " << connected_to_hyphae->hyphalTreeIndex << "\n";

                //std::cout<< "OrganID: " << h->getId() << " SDF" << sdf->treeIds_.at(distID)<< std::endl;
                h->setActive(false); // deactivate hyphae after anastomosis
                std::dynamic_pointer_cast<Hyphae>(h)->setMergePointID(lastIndex+1); // set node ID where anastomosis happened
                std::dynamic_pointer_cast<Hyphae>(h)->setMergedHyphae(connected_to_hyphae); // set merged hyphae
            }
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
        inf_->sdf = std::make_shared<SignedDistanceFunction>(sdf);
        inf_->fmax = rp->lmbd;
        inf_->fmin = 0;
        inf_->slope = 0;
        // std::cout<< inf_->sdf->toString() << std::endl;
        // rp->f_inf = inf_; // set new one
    }

    // Create tropisms and growth functions per random hyphae parameter
    for (auto& p_otp :organParam[Organism::ot_hyphae]) {
        auto rp = std::static_pointer_cast<HyphaeRandomParameter>(p_otp.second);
        auto tropism = this->createTropismFunction(rp->tropismT, rp->tropismN, rp->tropismS);
        tropism->setGeometry(geometry);
        rp->f_tf = tropism; // set new one
        // growth function is set to LinearGrowth in constructor of HyphaeRandomParameter
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

};



// void MycorrhizalPlant::updateAnastomosisTree(double dt) {
//     // tree = aabb::Tree();
//     double currentSimTime = this->getSimTime();
//     double previousSimTime = currentSimTime - dt;
//     Vector3d pos;
//     std::vector<double> posVec;
//     double radius;
//     double creationTime;

//     auto hyphae = this->getOrgans(Organism::ot_hyphae);
//     for (const auto& h : hyphae) {
//         for (unsigned int i=0; i < h->getNumberOfNodes(); i++) {
//             pos = h->getNode(i);
//             posVec = {pos.x, pos.y, pos.z};
//             radius = h->getParameter("a");
//             creationTime = h->getNodeCT(i);
//             if (creationTime <= currentSimTime && creationTime > previousSimTime) {
//                 // node is new if creation time is within current sim time and sim time - dt
//                 tree.insertParticle(h->getNodeId(i), posVec, radius);
//             } else {
//                 tree.updateParticle(h->getNodeId(i), posVec, radius);
//             }
//         }
//     }
//     tree.rebuild();
// }




}
