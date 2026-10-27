#include "mycorrhizalrootparameter.h"
#include "hyphaeparameter.h"
#include "MycorrhizalRoot.h"
#include "Hyphae.h"
#include "Root.h"
#include "Stem.h"
#include "Leaf.h"
#include "Organ.h"
#include "Organism.h"
#include "math.h"
#include <cmath>

namespace CPlantBox {

MycorrhizalRoot::MycorrhizalRoot(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
    Vector3d partialIHeading_, int pni, bool moved, int oldNON)
         :Root(id, param, alive, active, age, length,
             partialIHeading_,pni, moved,  oldNON ) {}

MycorrhizalRoot::MycorrhizalRoot(std::shared_ptr<Organism> rs, int type,  double delay, std::shared_ptr<Organ> parent, int pni)
:Root(rs,type, delay,parent, pni) {
    if(!(parent->organType()==Organism::ot_seed))
    {
        if (!parent->hasRelCoord())  // the first node of the base roots must be created in RootSystem::initialize()
        {
            colonized.push_back(0);
            emergedHyphae.push_back(0);
            hyphalEmergenceDone.push_back(false);
            colonizationTime.push_back(-1);

        }else{
            if ((parent->organType()==Organism::ot_stem)&&(parent->getNumberOfChildren()>0)) {
            }
            colonized.push_back(0);
            emergedHyphae.push_back(0);
            hyphalEmergenceDone.push_back(false);
            colonizationTime.push_back(-1);
        }
    }
}
/**
 * Adds a node to the root.
 *
 * For simplicity nodes can not be deleted, roots can only become deactivated or die
 *
 * @param n        new node
 * @param id       global node index
 * @param t        exact creation time of the node
 * @param index	   position were new node is to be added (between the nodes index-1 and index, if shift)
 * @param shift	   do we need to shift the nodes? (i.e., is the new node inserted between existing nodes to refine a colonized segment?)
 *
 * An inserted node takes the colonization state and time of the node at index (the node behind it), as in Stem::addNode
 * the global node indices are not kept in order along the root
 */
void MycorrhizalRoot::addNode(Vector3d n, int id, double t, size_t index, bool shift) {
    if (!shift)
    {
        Organ::addNode(n, id,  t,  index, shift);
        colonized.push_back(0);
        emergedHyphae.push_back(0);
        hyphalEmergenceDone.push_back(false);
        colonizationTime.push_back(-1);
    }
    else {
        assert((index > 0) && (index < nodes.size()) && "MycorrhizalRoot::addNode: nodes can only be inserted between existing nodes");
        int col = colonized.at(index); // copies, the vectors are modified below
        double colTime = colonizationTime.at(index);
		nodes.insert(nodes.begin() + index, n);//add the node at index
		nodeIds.push_back(id);
		nodeCTs.insert(nodeCTs.begin() + index, t);
        colonized.insert(colonized.begin() + index, col);
        emergedHyphae.insert(emergedHyphae.begin() + index, 0);
        hyphalEmergenceDone.insert(hyphalEmergenceDone.begin() + index, false);
        colonizationTime.insert(colonizationTime.begin() + index, std::max(colTime, t));

        for(auto kid : children){//if carries children after the added node, update their "parent node index"
			if(kid->parentNI >= index){
				kid->moveOrigin(kid->parentNI + 1);
			}
		}
    }
}

/**
 * Organ::calcCreationTime returns nodeCTs[0] + age for the newest node, which is summed up separately from the
 * simulation time and can exceed it by rounding errors. Hyphae and colonization use the end of the current time step
 * (simtime + dt, simtime is increased after the organs are simulated), so node creation times are limited to it.
 */
double MycorrhizalRoot::calcCreationTime(double length, double dt)
{
    auto plant_ = plant.lock();
    return std::min(Organ::calcCreationTime(length, dt), plant_->getSimTime() + plant_->getDt());
}

std::shared_ptr<Organ> MycorrhizalRoot::copy(std::shared_ptr<Organism> rs)
{
    auto r = std::make_shared<MycorrhizalRoot>(*this); // shallow copy
    r->parent = std::weak_ptr<Organ>();
    r->plant = rs;
    r->param_ = std::make_shared<MycorrhizalRootSpecificParameter>(*param()); // copy parameters
    for (size_t i=0; i< children.size(); i++) {
        r->children[i] = children[i]->copy(rs); // copy laterals
        r->children[i]->setParent(r);
    }
    return r;
}

void MycorrhizalRoot::primaryColonization(double dt, bool silence){
    double lmbd;
    double highres = getRootRandomParameter()->highresolution;
    // absolute time at the end of the current time step (age is the age of this root, not the simulation time;
    // plant->getSimTime() is only increased after all organs are simulated), node creation times are limited to it (see calcCreationTime)
    double t = plant.lock()->getSimTime() + plant.lock()->getDt();
    // Determine the appropriate colonization rate for each node based on soil properties and age
    for (size_t i = 1; i < nodes.size(); i++){
        if (std::dynamic_pointer_cast<SoilLookUpSDF>(getRootRandomParameter()->f_inf))
        {
            lmbd = getRootRandomParameter()->f_inf->getValue(nodes.at(i), shared_from_this());
        }
        else {
            lmbd = getRootRandomParameter()->lmbd;
        }
        if (t - nodeCTs.at(i) < getRootRandomParameter() ->minAge) {lmbd = 0;}//account for minimal age (of the node) in rate
        lmbd = (1 - (t - nodeCTs.at(i))/getRootRandomParameter()->maxAge) * lmbd; // account for maximal age (of the node) in rate
        // Determine the probability for colonization of the current node and infect if successful
        // Also refine the root if the segment is too long and a new node is inserted, which inherits the colonization status based on parent node
        double cursegLength = (nodes.at(i).minus(nodes.at(i-1))).length();
        if (colonized.at(i) == 0 && plant.lock()->rand() < lmbd*cursegLength*dt)
        {
            // insert node here if segment too long and set all nodes to be colonized
            setColonization(i,1,t);
            if (highres >= 1. && cursegLength > getRootRandomParameter() ->dx_inf) {
                int newNodesNumber = std::max( int(cursegLength / getRootRandomParameter() ->dx_inf) - 1, 0);
                Vector3d fromNode = nodes.at(i-1);   // freeze endpoints before any insertion
                Vector3d toNode   = nodes.at(i);
                double toCT = nodeCTs.at(i);
                for (size_t j = 0; j < newNodesNumber; j++)
                {
                    double newx = fromNode.x + (toNode.x - fromNode.x) *(j+1)/(newNodesNumber +1) ;
                    double newy = fromNode.y + (toNode.y - fromNode.y) *(j+1)/(newNodesNumber +1);
                    double newz = fromNode.z + (toNode.z - fromNode.z) *(j+1)/(newNodesNumber +1);
                    Vector3d newNode = Vector3d(newx,newy,newz);
                    addNode(newNode,plant.lock()->getNodeIndex(), toCT, i, true);
                    // insertion happens between i-1 and i, shifting i (the colonized endpoint) up by 1
                    i++;
                }
            }
        }
    }
}

// TODO: this functino is a bit hard to read, might be good to break it up in smaller parts.
void MycorrhizalRoot::secondaryColonization(bool silence, double dt){
    double max_length_colonization = age*getRootRandomParameter()->vi; // age used as duration here
    double t = plant.lock()->getSimTime() + plant.lock()->getDt(); // absolute time at the end of the current time step, see primaryColonization
    double infTime;
    double highres = getRootRandomParameter()->highresolution;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        if (colonized.at(i) == 1 || colonized.at(i)== 3)
        {
            int oldNode = i;
            double colonizationLength = 0;
            if (i>=1) {  // secondary colonization in basal direction can only occur if there is another node in basal direction in this root
                int basalnode = i-1;
                double cursegLength;
                while(basalnode > 0) {
                    cursegLength = abs(nodes.at(oldNode).minus(nodes.at(basalnode)).length());
                    colonizationLength += cursegLength;

                    if (colonizationLength > max_length_colonization) {break;}

                    if (colonized.at(basalnode) != 0) {
                        // already colonized (this call or earlier) — valid anchor, keep walking
                        oldNode = basalnode;
                        basalnode--;
                        continue;
                    }

                    infTime = colonizationTime.at(oldNode) + cursegLength/getRootRandomParameter()->vi;
                    if (infTime > t) {break;}  // front hasn't reached here yet; farther nodes are worse

                    setColonization(basalnode,2,std::max(infTime, nodeCTs.at(basalnode))); // node can only be colonized after it is created

                    if (highres >= 1. && cursegLength > getRootRandomParameter() ->dx_inf) {
                        int newNodesNumber = std::max( int(cursegLength / getRootRandomParameter() ->dx_inf) - 1, 0);
                        Vector3d fromNode = nodes.at(oldNode);
                        Vector3d toNode   = nodes.at(basalnode);
                        double toCT = nodeCTs.at(basalnode);
                        for (size_t j = 0; j < newNodesNumber; j++)
                        {
                            double newx = fromNode.x + (toNode.x - fromNode.x) *(j+1)/(newNodesNumber +1);
                            double newy = fromNode.y + (toNode.y - fromNode.y) *(j+1)/(newNodesNumber +1);
                            double newz = fromNode.z + (toNode.z - fromNode.z) *(j+1)/(newNodesNumber +1);
                            Vector3d newNode = Vector3d(newx,newy,newz);
                            // insertion happens between basalnode and basalnode+1 (each new node closer to basalnode than the one before),
                            // the new node is colonized like basalnode, not like the node behind it
                            addNode(newNode,plant.lock()->getNodeIndex(), toCT, basalnode+1, true);
                            setColonization(basalnode+1, colonized.at(basalnode), colonizationTime.at(basalnode));
                            i++; // the start node i is shifted by the insertion
                        }
                    }

                    // TODO: Y is this hier? we will never get basalnode==0  in that loop
                    if(basalnode==0 && std::dynamic_pointer_cast<MycorrhizalRoot>(getParent()))
                    {
                        std::dynamic_pointer_cast<MycorrhizalRoot>(getParent())->setColonization(parentNI,3,infTime);
                        std::dynamic_pointer_cast<MycorrhizalRoot>(getParent())->simulateColonization(dt,silence);
                    }

                    oldNode = basalnode;
                    basalnode--;
                }
            }

            auto apicalnode = i+1;
            oldNode = i;
            colonizationLength = 0;
            double cursegLength;
            while (apicalnode < nodes.size())
            {
                cursegLength = abs(nodes.at(oldNode).minus(nodes.at(apicalnode)).length());
                colonizationLength += cursegLength;

                if (colonizationLength > max_length_colonization) {break;}

                if (colonized.at(apicalnode) != 0) {
                    oldNode = apicalnode;
                    apicalnode++;
                    continue;
                }

                infTime = colonizationTime.at(oldNode) + cursegLength/getRootRandomParameter()->vi;
                if (infTime > t) {break;}

                setColonization(apicalnode,2,std::max(infTime, nodeCTs.at(apicalnode))); // node can only be colonized after it is created (the tip node gets a new creation time when it is moved)
                if (highres >= 1. && cursegLength > getRootRandomParameter() ->dx_inf) {
                    int newNodesNumber = std::max( int(cursegLength / getRootRandomParameter() ->dx_inf) - 1, 0);
                    Vector3d fromNode = nodes.at(oldNode);
                    Vector3d toNode   = nodes.at(apicalnode);
                    double toCT = nodeCTs.at(apicalnode);
                    for (size_t j = 0; j < newNodesNumber; j++)
                    {
                        double newx = fromNode.x + (toNode.x - fromNode.x) *(j+1)/(newNodesNumber +1);
                        double newy = fromNode.y + (toNode.y - fromNode.y) *(j+1)/(newNodesNumber +1);
                        double newz = fromNode.z + (toNode.z - fromNode.z) *(j+1)/(newNodesNumber +1);
                        Vector3d newNode = Vector3d(newx,newy,newz);
                        addNode(newNode,plant.lock()->getNodeIndex(), toCT, apicalnode, true); // insertion between apicalnode-1 and apicalnode
                        apicalnode++;
                    }
                }
                oldNode = apicalnode;
                apicalnode++;
            }
        }
    }
}

void MycorrhizalRoot::simulateSecondaryColonization(double dt) {
    secondaryColonization(false,dt);
}

void MycorrhizalRoot::simulatePrimaryColonization(double dt) {
    primaryColonization(dt,false);
}

void MycorrhizalRoot::simulateHyphalGrowth(double dt, bool verbose) {
    if (getRootRandomParameter()->highresolution >= 1) { // Version where at every node there is one hypha created
        for (size_t i = 1; i < nodes.size(); i++) {
            if (colonized.at(i) > 0 && emergedHyphae.at(i) == 0){ // if the current node is colonized and the number of hyphae to be created is reached
                createHyphae(i); // increases emergedHyphae.at(i)
            }
        }
    } else { // Version where the number of hyphae of a node is based on hyphal emergence density and the length of the segment ending at that node
        double hed = getRootRandomParameter()->hyphalEmergenceDensity;
        double s = 0.; // length along the root from its base to node i
        for (size_t i = 1; i < nodes.size(); i++) {
            double s_prev = s;
            s += nodes.at(i).minus(nodes.at(i-1)).length();
            if (colonized.at(i) > 0 && !hyphalEmergenceDone.at(i)) { // hyphae emerge in the time step the node is colonized
                // one hypha every 1/hed cm along the root: the segment ending at node i gets the multiples of 1/hed within (s_prev, s],
                // so the hyphae are evenly distributed independent of the order of colonization (tolerance for rounding)
                int noh = int(std::floor(hed * s + 1.e-9) - std::floor(hed * s_prev + 1.e-9));
                for (int k = 0; k < noh; k++) {
                    createHyphae(i);
                }
                hyphalEmergenceDone.at(i) = true;
            }
        }
    }
}



void MycorrhizalRoot::simulateColonization(double dt, bool verbose) {

    if (this->nodes.size()>1) {

        //Primary Colonization
        primaryColonization(dt,verbose);

        // Secondary Colonization
        secondaryColonization(verbose,dt);

        // Transfer to the laterals. The laterals are not colonized here, they colonize themselves in their own simulate()
        // (called by Root::simulate() before), MycorrhizalPlant::simulateColonization() calls every root
        for (auto l : children)
        {
            if (l->organType()==Organism::ot_root) {
                if (colonized.at(l->parentNI) != 0) { // the base of root l is colonized
                    if (l->getNumberOfNodes() > 1 && std::dynamic_pointer_cast<MycorrhizalRoot>(l) -> getNodeColonization(1) == 0) {
                        std::dynamic_pointer_cast<MycorrhizalRoot>(l) ->setColonization(0, 3,
								std::max( l->getNodeCT(0), colonizationTime.at(l->parentNI))); //root can only be colonized after it is born (and need to account for growth delay)
                    }
                }
            }
        }
    }
}


void MycorrhizalRoot::simulate(double dt, bool verbose)
{
    Root::simulate(dt,verbose);
    simulateColonization(dt,verbose);
    simulateHyphalGrowth(dt,verbose);
}

std::shared_ptr<const MycorrhizalRootSpecificParameter> MycorrhizalRoot::param() const
{
    return std::static_pointer_cast<const MycorrhizalRootSpecificParameter>(param_);
}

std::shared_ptr<MycorrhizalRootRandomParameter> MycorrhizalRoot::getRootRandomParameter() const
{
    return std::static_pointer_cast<MycorrhizalRootRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_root, param_->subType));
}

double MycorrhizalRoot::getParameter(std::string name, std::map<std::string, double> addParams) const {
    if (name == "primaryColonization")
    {
        double primaryColonizedLength = 0;
        for (size_t i = 1; i < nodes.size(); i++)
        {
            if (colonized.at(i)==1){primaryColonizedLength += nodes.at(i).minus(nodes.at(i-1)).length();}
        }
        return primaryColonizedLength;
    }
    if (name == "secondaryColonization")
    {
        double secondaryColonizedLength = 0;
        for (size_t i = 1; i < nodes.size(); i++)
        {
            if (colonized.at(i)>1) {secondaryColonizedLength += nodes.at(i).minus(nodes.at(i-1)).length();}
        }
        return secondaryColonizedLength;
    }
    if (name == "colonizationLength")
    {
        double colonizedLength = 0;
        for (size_t i = 1; i < nodes.size(); i++)
        {
            if (colonized.at(i)>0) {colonizedLength += nodes.at(i).minus(nodes.at(i-1)).length();}
        }
        return colonizedLength;
    }
    if (name == "hyphalTreeIndex") {return hyphalTreeIndex;}
    return Root::getParameter(name, addParams);
}

void MycorrhizalRoot::setColonization(int i, int colonization, double t)
{
    colonized.at(i) = colonization;
    colonizationTime.at(i) = t;
	assert(colonizationTime.at(i) >= 0 && "MycorrhizalRoot::setColonization colonizationTime.at(i) < 0");
	assert(colonizationTime.at(i) >= nodeCTs.at(i) && "MycorrhizalRoot::setColonization colonizationTime.at(i) < nodeCTs.at(i)");
}

void MycorrhizalRoot::createHyphae(int pni)
{
    double delay = getRootRandomParameter()->hyphalDelay;
    double emergenceTime = colonizationTime.at(pni) + delay;
    double parentCT = getNodeCT(pni);
	double delay_for_creation = emergenceTime - parentCT; //difference between creation of parent node and that of the hyphae
    if (parentCT + delay_for_creation > emergenceTime) { // the Hyphae constructor uses parentCT + delay, which must not round above emergenceTime
        delay_for_creation = std::nextafter(delay_for_creation, 0.);
    }
    assert(delay_for_creation >= 0 && "MycorrhizalRoot::createHyphae delay_for_creation < 0");
    // as in Organ::createLateral: the constructor sets age = -delay_for_creation, so the hypha has to be simulated for the time
    // since the creation of the parent node (not since emergence), then it has age 0 at emergence (colonizationTime + delay).
    // End of the current time step, plant->getSimTime() is only increased after all organs are simulated
    double dt_ = plant.lock()->getSimTime() + plant.lock()->getDt() - parentCT;
    int subType = 1;
    auto hyphae = std::make_shared<Hyphae>(plant.lock(), subType,  delay_for_creation, shared_from_this(), pni);
    children.push_back(hyphae);
    emergedHyphae.at(pni) += 1;
    hyphae->setHyphalTreeIndex(-1);
    hyphae->simulate(dt_);
}

std::string MycorrhizalRoot::toString() const
{
    // TODO add additional stuff
    std::stringstream newstring;
    newstring << "; number of colonized Nodes " << getNumberofColonizedNodes() << "; length of colonized root segments "<< getParameter("colonizationLength")<< ".";
    return  Root::toString()+newstring.str();
}

int MycorrhizalRoot::getNumberofColonizedNodes() const
{
    int numberColonizedNodes =0;
    for (size_t i = 0; i < getNumberOfNodes()-1; i++)
    {
        if (colonized.at(i)!=0)
        {
            numberColonizedNodes++;
        }
    }
    return numberColonizedNodes;
}

}
