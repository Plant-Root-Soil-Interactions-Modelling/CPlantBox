#include "Plant.h"
#include "Organ.h"
#include "Organism.h"
#include "MycorrhizalPlant.h"
#include "Mycorrhizalrootparameter.h"



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