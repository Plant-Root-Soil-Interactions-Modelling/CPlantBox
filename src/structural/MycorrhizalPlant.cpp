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
    auto hyphrp = std::make_shared<HyphaeRandomParameter>(shared_from_this());
    hyphrp -> subType = 0;
    setOrganRandomParameter(hyphrp);
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
    for (const auto& o : baseOrgans) {
        if(o->organType() == Organism::ot_root) {
            infs.at(o->getNodeId(0)) = std::dynamic_pointer_cast<MycorrhizalRoot> (o) -> getNodeInfection(0);
        }
    }

    for (const auto & o : organs) {
        for (size_t i = 1; i < o ->getNumberOfNodes()-1; i++) { // TOOO check if here we need to check for organ type!!!
            infs.at(o->getNodeId(i)) = std::dynamic_pointer_cast<MycorrhizalRoot> (o) -> getNodeInfection(i);
        }
    }
    return infs;
}

std::vector<double> MycorrhizalPlant::getNodeInfectionTime(int ot) const {
    auto organs = this -> getOrgans(ot);
    std::vector<double> infTime = std::vector<double>(getNumberOfNodes());
    for (const auto& o : baseOrgans) {
        if(o->organType() == Organism::ot_root){
            infTime.at(o->getNodeId(0)) = std::dynamic_pointer_cast<MycorrhizalRoot> (o) -> getNodeInfectionTime(0);
        }
    }

    for (const auto & o : organs){
        for (size_t i = 1; i < o ->getNumberOfNodes()-1; i++){ // TOOO check if here we need to check for organ type!!!
            infTime.at(o->getNodeId(i)) = std::dynamic_pointer_cast<MycorrhizalRoot> (o) -> getNodeInfectionTime(i);
        }
    }
    return infTime;
}

std::vector<int> MycorrhizalPlant::getAnastomosisPoints(int ot) const {
    auto organs = this -> getOrgans(ot);
    std::vector<int> anaPoints = std::vector<int>(getNumberOfNodes());
    for (const auto& o : baseOrgans) {
        if(o->organType() == Organism::ot_hyphae) {
            auto h = std::dynamic_pointer_cast<Hyphae>(o);
            if (h->mergePointID != -1) {
                anaPoints.at(o->getNodeId(h->getNumberOfNodes()-1)) = 1;
            }
            else {
                anaPoints.at(o->getNodeId(h->getNumberOfNodes()-1)) = 0;
            }
        }
    }

    for (const auto & o : organs) {
        auto h = std::dynamic_pointer_cast<Hyphae>(o);
            if (o->organType() == Organism::ot_hyphae) {
                if (h->mergePointID != -1) {
                    anaPoints.at(o->getNodeId(h->getNumberOfNodes()-1)) = 1;
                }
                else {
                    anaPoints.at(o->getNodeId(h->getNumberOfNodes()-1)) = 0;
                }
            }
    }
    return anaPoints;
}
/**
 * Simulates mycorrhizal plant growth
 * @param dt		duration of the simulation
 * @param verbose	whether to print information
 */
void MycorrhizalPlant::simulate(double dt, bool verbose)
{
    abs2rel();
    Organism::simulate(dt, verbose);
    sdf = std::make_shared<SDF_RootSystem>(*this);
    sdf->selectedOrganType = Organism::ot_hyphae;
    simulateAnastomosis(dt, verbose);
    rel2abs();
}

void MycorrhizalPlant::simulateHyphae(double dt, bool verbose) {
    auto hyphae = getOrgans(5);
    this->dt = dt;
    oldNumberOfNodes = getNumberOfNodes();
    oldNumberOfOrgans = getNumberOfOrgans();
    for (const auto & h : hyphae) {
        std::dynamic_pointer_cast<Hyphae>(h) -> simulate(dt,verbose);
    }
    simtime += dt;
    sdf = std::make_shared<SDF_RootSystem>(*this);
    sdf->selectedOrganType = Organism::ot_hyphae;
    simulateAnastomosis(dt, verbose);
}
/*
 * Simulates primary infection for all mycorrhizal roots
 * @param dt		duration of the simulation
 * @param verbose	whether to print information
*/
void MycorrhizalPlant::simulatePrimaryInfection(double dt, bool verbose) {
    for (const auto& r : baseOrgans) {
        if (r->organType() == Organism::ot_root) {
            std::dynamic_pointer_cast<MycorrhizalRoot>(r)->simulatePrimaryInfection(dt);
        }
    }
    auto organs = getOrgans();
    for (const auto& r : organs) {
        if (r->organType() == Organism::ot_root) {
            std::dynamic_pointer_cast<MycorrhizalRoot>(r)->simulatePrimaryInfection(dt);
        }
    }
}

/*
 * Simulates secondary infection for all mycorrhizal roots
 * @param dt		duration of the simulation
 * @param verbose	whether to print information
 */
void MycorrhizalPlant::simulateSecondaryInfection(double dt, bool verbose) {
    for (const auto& r : baseOrgans) {
        if (r->organType() == Organism::ot_root) {
            std::dynamic_pointer_cast<MycorrhizalRoot>(r)->simulateSecondaryInfection(dt);
        }
    }

    auto organs = getOrgans();
    for (const auto& r : organs) {
        if (r->organType() == Organism::ot_root) {
            std::dynamic_pointer_cast<MycorrhizalRoot>(r)->simulateSecondaryInfection(dt);
        }
    }
}

void MycorrhizalPlant::simulateInfection(double dt, bool verbose) {
    for (const auto& r : baseOrgans) {
        if (r->organType() == Organism::ot_root) {
            std::dynamic_pointer_cast<MycorrhizalRoot>(r)->simulateInfection(dt, verbose);
        }
    }
    
    auto organs = getOrgans();
    for (const auto& r : organs) {
        if (r->organType() == Organism::ot_root) {
            std::dynamic_pointer_cast<MycorrhizalRoot>(r)->simulateInfection(dt, verbose);
        }
    }
}

/*
 * Simulates hyphal growth
 * @param dt		duration of the simulation
 * @param verbose	whether to print information
 */
void MycorrhizalPlant::simulateHyphalGrowth(double dt, bool verbose)
{
    oldNumberOfOrgans = getNumberOfOrgans();
    auto organs = getOrgans();
    for (const auto & o : organs) {
        if (o->organType() == Organism::ot_root) {
            std::dynamic_pointer_cast<MycorrhizalRoot>(o) -> simulateHyphalGrowth(dt,verbose);
        }
    }
};

/*
 * Simulates anastomosis for all hyphae
 * @param dt		duration of the simulation
 * @param verbose	whether to print information
 */
void MycorrhizalPlant::simulateAnastomosis(double dt, bool verbose) {
    auto hyphae = this->getOrgans(Organism::ot_hyphae);
    double dist = 1000;

    for (const auto & h : hyphae) {
        sdf->excludeTreeId = h->getParameter("hyphalTreeIndex");

        if (h->isActive() && h->getParameter("subType") == 1) { // only check for anastomosis for active hyphae tips
            auto tip = h->getNode(h->getNumberOfNodes()-1);
            dist = sdf->getDist(tip);
            if (fabs(dist) < h->getParameter("distTH") && rand() < h->getParameter("ana")) 
            {
                auto lastIndex = sdf->distIndex; 
                auto connected_to_hyphae = std::dynamic_pointer_cast<Hyphae>(sdf->lastOrgan.lock());
                int locallastIndex = -1;

                for (size_t i = 0; i < connected_to_hyphae->getNumberOfNodes(); i++)
                {
                    if (connected_to_hyphae->getNodeId(i) == lastIndex)
                    {
                        locallastIndex = i;
                        break;
                    }
                }
                
                if (verbose) {
                    std::cout <<"Anastomosis at tip: " << tip.toString() <<" to node: " << connected_to_hyphae->getNode(locallastIndex).toString() << std::endl;
                    std::cout << "Anastomosis occurred at distance: " << dist << " cm.\n";
                    std::cout << "Hyphal tree index " << h->getParameter("hyphalTreeIndex") << "\n";
                }
                
                h->setActive(false); // deactivate hyphae after anastomosis
                std::dynamic_pointer_cast<Hyphae>(h)->setMergePointID(lastIndex); // set node ID where anastomosis happened
                std::dynamic_pointer_cast<Hyphae>(h)->setMergedHyphae(connected_to_hyphae); // set merged hyphae
                std::dynamic_pointer_cast<Hyphae>(h)->addNode(connected_to_hyphae->getNode(locallastIndex), getNodeIndex(), h->getNodeCT(h->getNumberOfNodes()-1), h->getNumberOfNodes(),false); // add anastomosis point as new node to the hyphae
            }
        }

    }

};

void MycorrhizalPlant::initCallbacks() {
    Plant::initCallbacks();

    // Create tropisms and growth functions per random hyphae parameter
    for (auto& p_otp :organParam[Organism::ot_hyphae]) {
        auto rp = std::static_pointer_cast<HyphaeRandomParameter>(p_otp.second);
        auto tropism = this->createTropismFunction(rp->tropismT, rp->tropismN, rp->tropismS);
        // std::cout << "Created tropism function for hyphae with type " << rp->tropismT << " and parameters n=" << rp->tropismN << " and sigma=" << rp->tropismS << std::endl;
        tropism->setGeometry(geometry);
        rp->f_tf = tropism; // set new one
    }

};

void MycorrhizalPlant::changeGeometry(int ot,std::shared_ptr<SignedDistanceFunction> geom) {
    for (auto& p_otp :organParam[ot]) {
        auto rp = std::static_pointer_cast<HyphaeRandomParameter>(p_otp.second);
        auto tropism = this->createTropismFunction(rp->tropismT, rp->tropismN, rp->tropismS);
        tropism->setGeometry(geom);
        rp->f_tf = tropism; // set new one
    }
};

std::vector<int> MycorrhizalPlant::getNodeTips( int ot) const {
    auto organs = this -> getOrgans(ot);
    std::vector<int> tips(getNumberOfNodes(), 0);
    for (const auto& o : baseOrgans) {
        if (o->organType() == ot && o->isActive()) {
            tips.at(o->getNodeId(o->getNumberOfNodes()-1)) = 1;
        }    
    }

    for (const auto & o : organs) {
        if (o->isActive()) {
            tips.at(o->getNodeId(o->getNumberOfNodes()-1)) = 1;
        }
    }
    return tips;
};
}
