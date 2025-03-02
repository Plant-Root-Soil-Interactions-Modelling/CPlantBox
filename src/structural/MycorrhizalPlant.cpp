#include "Plant.h"
#include "Organ.h"
#include "Organism.h"
#include "MycorrhizalPlant.h"
#include "MycorrhizalRoot.h"
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

std::vector<double> MycorrhizalPlant::getNodeIT(int ot) const {
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
    for (const auto & o : organs)
    {
        if (o->organType() == Organism::ot_root)
        {
            std::dynamic_pointer_cast<MycorrhizalRoot>(o) -> simulateInfection(dt,verbose);
        }
    }
	rel2abs();
}
}