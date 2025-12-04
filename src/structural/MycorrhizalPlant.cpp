#include "Plant.h"
#include "Organ.h"
#include "Organism.h"
#include "MycorrhizalPlant.h"
//#include "MycorrhizalRoot.h"
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

    // buildAnastomosisTree();
    simulateAnastomosis();
    // sdfs = {};
    // buildAnastomosisTree();
    // updateAnastomosisTree(dt);
    // simulateAnastomosisTree();
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
    // auto numberofHyphae = sdfs.size();
    // auto numberofHyphae2 = hyphae.size();
    // std::cout<< numberofHyphae << " "<< numberofHyphae2 <<std::endl;
    double dist = 1000;
    // Vector3d closestNode;
    

    for (const auto & h : hyphae) {
        sdf->excludeTreeId = h->getParameter("hyphalTreeIndex");
        if (h->isActive()) {
            auto tip = h->getNode(h->getNumberOfNodes()-1);

            dist = sdf->getDist(tip);
            // std::cout<<h->getParameter("distTH") << std::endl;
            // std::cout<<"Distance to nearest hyphae from tip " << tip.toString() << " is " << dist << " cm." << std::endl;

            if (dist < h->getParameter("distTH"))
            {
                auto distID = sdf->distIndex;
                std::cout <<"Anastomosis at tip: " << tip.toString() <<" with distance id: " << distID << std::endl;
                std::cout<< "OrganID: " << h->getId() << " SDF" << sdf->treeIds_.at(distID)<< std::endl;
                // std::cout <<"Node for Anastomosis: " << closestNode.toString() << std::endl;

                // for (const auto & hh : hyphae)
                // {

                //     auto nodes = hh->getNodes();
                //     for (const auto & n : nodes)
                //     {
                //         if (n == closestNode)
                //         {
                //             std::cout<< "Found the closes node!" << std::endl;
                //         }
                        
                //     }
                    
                // }
                
            }
        }

        
        
    }
    
};

void MycorrhizalPlant::simulateAnastomosisTree() {
    auto hyphae = this->getOrgans(Organism::ot_hyphae);
    unsigned int dist;
    Vector3d closestNode;
    for (const auto & h : hyphae) {
        auto tipID = h->getNodeId(h->getNumberOfNodes()-1);
        auto tip = h->getNode(h->getNumberOfNodes()-1);
        auto radius = h->getParameter("a"); // hyphal diameter plus threshold
        dist = getDistTree(tipID, tip, radius*2 + h->getParameter("distTH"));

        std::cout<<"Anastomosis between: " << tipID << " and " << dist << " ." << std::endl;
        // std::vector<double> posVec = {tip.x, tip.y, tip.z};
        // std::vector<aabb::Tree::ParticleInfo> results;
        // tree.queryPoint(posVec, results);
        // for (const auto& r : results) {
        //     // get position of the particle
        //     // double px = r.position[0];
        //     // double py = r.position[1];
        //     // double pz = r.position[2];
        //     Vector3d partPos = Vector3d(r.position[0], r.position[1], r.position[2]);
        //     double distfromsdf = partPos.minus(tip).length();
        //     closestNode = partPos;
        //     if ( distfromsdf > 0 && distfromsdf < dist) dist = distfromsdf;
        // }
        // if (dist < h->getParameter("distTH"))
        // {
        //     std::cout <<"Anastomosis at tip: " << tip.toString() <<" with distance: " << dist << std::endl;
        //     std::cout <<"Node for Anastomosis: " << closestNode.toString() << std::endl;
        // }
        
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

}


// void MycorrhizalPlant::addTree() {
//     // TODO adapt from exudation model or push to hyphae???
//         auto hyphae = this->getOrgans(Organism::ot_hyphae);

//         for (const auto& h : hyphae) {
//             if (h->getNumberOfNodes()>1) { // started growing
//                 // time when the root stopped growing
//                 // root tip
//                 Vector3d t = h->getNode(h->getNumberOfNodes()-1);
//                 // tips.push_back(t);
//                 // // direction towards root base
//                 // Vector3d base = h->getNode(0);
//                 // double a = h->getNodeCT(h->getNumberOfNodes()-1) - h->getNodeCT(0);
//                 // v.push_back(base.minus(t).times(1./a));
//             }
//             // sdfs.push_back(SDF_RootSystem(*(std::dynamic_pointer_cast<Hyphae>(h)), h->getParameter("dx")));

//         }
// }

// void MycorrhizalPlant::buildAnastomosisTree() {
//     auto hyphae = this->getOrgans(Organism::ot_hyphae);
//     std::vector<Vector3d> nodes;
//     std::vector<Vector2i> segments;
//     std::vector<double> radii;
//     double dx_;
//     for (const auto& h : hyphae) {
//         for (unsigned int i=0; i < h->getNumberOfNodes(); i++) {
//             dx_ = h->getParameter("dx");
//             Vector3d pos = h->getNode(i);
//             nodes.push_back(pos);
//             segments.push_back(Vector2i(i,i+1));
//             radii.push_back(h->getParameter("a"));
            
//             // tree.insertParticle(h->getNodeId(i), posVec, radius);
//             localNodes.push_back(i);
//             localHyphae.push_back(std::dynamic_pointer_cast<Hyphae>(h));
//         }
//     }
//     sdf = SDF_RootSystem(nodes,segments,radii,dx_);
// }   

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

unsigned int MycorrhizalPlant::getDistTree(unsigned int p, Vector3d  tip, double dist) {
    
    // get the needed distance
    double dx_ = dist; 

    // get the tree minus the particle we are trying to find the distance for
    // tree.removeParticle(p);


    // make the box around the point
    std::vector<double> a = { tip.x-dx_, tip.y-dx_, tip.z-dx_ };
    std::vector<double> b = { tip.x+dx_, tip.y+dx_, tip.z+dx_ };
    aabb::AABB box = aabb::AABB(a,b);
    double mdist = 1e100; // far far away

    // get the indices of all particles in the box
    auto indices = tree.query(box);
    std::cout << indices.size() << " segments in range\n";
    if (indices.size() == 1) {
        return p; // only itself
    }
    else if (indices.size()==2)
    {
        if (indices.at(0) == p) {
            return indices.at(1);
        }
        else {
            return indices.at(0);
        }
    }
    else {
        std::cout<< "more than two hyphal nodes nearby"<< std::endl;
        return p;
        // iterate over all hyphae to find the ones with the accurate global node id
        // safe the location and ids of the nodes
        // check the distance from the tip to these nodes
        // only if there are more than two
    }
}



}
