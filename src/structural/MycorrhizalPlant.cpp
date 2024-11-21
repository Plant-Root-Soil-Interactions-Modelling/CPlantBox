#include "Plant.h"
#include "MycorrhizalPlant.h"


namespace CPlantBox {
    MycorrhizalPlant::MycorrhizalPlant(unsigned int seednum): Plant(seednum) {}

    void Plant::initializeReader()
{
    auto rrp = std::make_shared<RootRandomParameter>(shared_from_this());
    rrp->subType = 0;
    setOrganRandomParameter(rrp);
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
}