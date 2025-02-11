#include "Plant.h"
#include "Organ.h"
#include "Organism.h"
#include "MycorrhizalPlant.h"
#include "Mycorrhizalrootparameter.h"



namespace CPlantBox {
    MycorrhizalPlant::MycorrhizalPlant(unsigned int seednum): Plant(seednum) { std::cout << "MycorrhizalPlant::MycorrhizalPlant called" << std::endl; }

    void MycorrhizalPlant::initializeReader()
{
    std::cout << "MycorrhizalPlant::initializeReader called" << std::endl;
    //new Parameters start here
    auto mycrrp = std::make_shared<MycorrhizalRootRandomParameter>(shared_from_this());
    mycrrp -> subType = 0;
    setOrganRandomParameter(mycrrp);
    //new Parameters end here
    // auto rrp = std::make_shared<RootRandomParameter>(shared_from_this());
    // rrp -> subType = 0;
    // setOrganRandomParameter(rrp);
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

std::shared_ptr<Organ> MycorrhizalPlant::createRoot(std::shared_ptr<Organism> plant, int type, double delay) {
    class MycorrhizalSeed : public Seed{
        using Seed::Seed;
        std::shared_ptr<Organ> createRoot(std::shared_ptr<Organism> plant, int type, double delay) override {
            return std::make_shared<MycorrhizalRoot>(plant, type, delay, shared_from_this(), 0);
        };
    };
    auto root = std::make_shared<MycorrhizalSeed>(shared_from_this());
    return root;
}
// std::vector<int> MycorrhizalPlant::getNodeInfections(int ot) const {
//     auto organs = this -> getOrgans(ot);
//     std::vector<int> infs = std::vector<int>(getNumberOfNodes());
//     std::cout << organs.size() << std::flush;
//     for (const auto& o : baseOrgans)
//     {
//         infs.at(o->getNodeId(0)) = std::dynamic_pointer_cast<MycorrhizalRoot> (o) -> getNodeInfection(0);
        
//     }
    
//     for (const auto & o : organs)
//     {
//         for (size_t i = 1; i < o ->getNumberOfNodes(); i++)
//         {
//             infs.at(o->getNodeId(i)) = std::dynamic_pointer_cast<MycorrhizalRoot> (o) -> getNodeInfection(i);
//         }
        
//     }
//     return infs;
// }

// std::vector<int> MycorrhizalPlant::getSegmentInfections(int ot) const {
//     auto nodeInfection = getNodeInfections(ot);
//     auto segments = getSegments(ot);
//     std::vector<int> Infections = std::vector<int>(segments.size());
//     for (int i = 0; i < Infections.size(); i++)
//     {
//         Infections[i] = nodeInfection[segments[i].y];
//     }
//     return Infections;
// }

}