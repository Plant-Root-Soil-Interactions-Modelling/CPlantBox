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
            colonizationTime.push_back(-1);

        }else{
            if ((parent->organType()==Organism::ot_stem)&&(parent->getNumberOfChildren()>0)) {
            }
            colonized.push_back(0);
            emergedHyphae.push_back(0);
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
 * @param index	   position were new node is to be added
 * @param shift	   do we need to shift the nodes? (i.e., is the new node inserted between existing nodes because of internodal growth?)
 */
void MycorrhizalRoot::addNode(Vector3d n, int id, double t, size_t index, bool shift) {
    if (!shift)
    {
        Organ::addNode(n, id,  t,  index, shift);
        colonized.push_back(0);
        emergedHyphae.push_back(0);
        colonizationTime.push_back(-1);
    }
    else {
		nodes.insert(nodes.begin() + index-1, n);//add the node at index
		nodeIds.push_back(id);
		nodeCTs.insert(nodeCTs.begin() + index-1, t);
        colonized.insert(colonized.begin()+index-1, colonized.at(index-1));
        emergedHyphae.insert(emergedHyphae.begin()+index-1, 0);
        colonizationTime.insert(colonizationTime.begin()+index-1, std::max(colonizationTime.at(index-1), t));

        for(auto kid : children){//if carries children after the added node, update their "parent node index"
			if((kid->parentNI >= index-1 )&&(kid->parentNI > 0)){
				kid->moveOrigin(kid->parentNI + 1);
			}
		}
    }
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
    // Determine the appropriate colonization rate for each node based on soil properties and age  
    for (size_t i = 1; i < nodes.size(); i++){
        if (std::dynamic_pointer_cast<SoilLookUpSDF>(getRootRandomParameter()->f_inf))
        {
            lmbd = getRootRandomParameter()->f_inf->getValue(nodes.at(i), shared_from_this());
        }
        else {
            lmbd = getRootRandomParameter()->lmbd;
        }
        if (age - nodeCTs.at(i) < getRootRandomParameter() ->minAge) {lmbd = 0;}//account for minimal age in rate
        lmbd = (1 - (age- nodeCTs.at(i))/getRootRandomParameter()->maxAge) * lmbd; // account for maximal age in rate
        // Determine the probability for colonization of the current node and colonize if successful
        // Also refine the root if the segment is too long and a new node is inserted, which inherits the colonization status based on parent node
        double cursegLength = (nodes.at(i).minus(nodes.at(i-1))).length();
        if (colonized.at(i) == 0 && plant.lock()->rand() < lmbd*cursegLength*dt)
        {
            // insert node here if segment too long and set all nodes to be colonized
            setColonization(i,1,age);
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
                    // insertion happens at i-1, shifting i (the colonized endpoint) up by 1
                    i++;
                }
            }
        }
    }
}

bool MycorrhizalRoot::isSecondaryColonizationAnchor(size_t node) const
{
    return colonized.at(node) == 1 || colonized.at(node) == 3;
}

double MycorrhizalRoot::segmentLength(size_t a, size_t b) const
{
    return abs(nodes.at(a).minus(nodes.at(b)).length());
}

double MycorrhizalRoot::colonizationArrivalTime(size_t fromNode, double cursegLength) const
{
    return colonizationTime.at(fromNode) + cursegLength / getRootRandomParameter()->vi;
}

void MycorrhizalRoot::refineSegment(size_t fromNode, size_t toNode, double toCT, int &oldNode, int &currentNode, int step)
{
    int newNodesNumber = std::max(int(segmentLength(fromNode, toNode) / getRootRandomParameter()->dx_inf) - 1, 0);
    if (newNodesNumber == 0) {
        return;
    }

    Vector3d fromNodePos = nodes.at(fromNode);
    Vector3d toNodePos = nodes.at(toNode);

    for (int j = 0; j < newNodesNumber; j++)
    {
        double alpha = double(j + 1) / (newNodesNumber + 1);
        Vector3d newNode(
            fromNodePos.x + (toNodePos.x - fromNodePos.x) * alpha,
            fromNodePos.y + (toNodePos.y - fromNodePos.y) * alpha,
            fromNodePos.z + (toNodePos.z - fromNodePos.z) * alpha
        );
        addNode(newNode, plant.lock()->getNodeIndex(), toCT, currentNode, true);

        if (step < 0) {
            currentNode++;
            oldNode++;
        } else {
            oldNode++;
            currentNode++;
        }
    }
}

/*
*/
void MycorrhizalRoot::propagateSecondaryColonization(size_t anchorNode, int step, double maxLength, bool silence, double dt)
{
    int currentNode = static_cast<int>(anchorNode) + step;
    int oldNode = static_cast<int>(anchorNode);
    double colonizationLength = 0;
    double highres = getRootRandomParameter()->highresolution;
    double dxInf = getRootRandomParameter()->dx_inf;

    while (currentNode >= 0 && currentNode < static_cast<int>(nodes.size()))
    {
        double cursegLength = segmentLength(oldNode, currentNode);
        colonizationLength += cursegLength;
        if (colonizationLength > maxLength) {
            break;
        }

        if (colonized.at(currentNode) != 0)
        {
            oldNode = currentNode;
            currentNode += step;
            continue;
        }

        double infTime = colonizationArrivalTime(oldNode, cursegLength);
        if (infTime > age) {
            break;
        }

        setColonization(currentNode, 2, infTime);
        if (highres >= 1. && cursegLength > dxInf)
        {
            refineSegment(oldNode, currentNode, nodeCTs.at(currentNode), oldNode, currentNode, step);
        }

        if (step == -1 && currentNode == 0 && std::dynamic_pointer_cast<MycorrhizalRoot>(getParent()))
        {
            std::dynamic_pointer_cast<MycorrhizalRoot>(getParent())->setColonization(parentNI, 3, infTime);
            std::dynamic_pointer_cast<MycorrhizalRoot>(getParent())->simulateColonization(dt, silence);
        }

        oldNode = currentNode;
        currentNode += step;
    }
}

void MycorrhizalRoot::secondaryColonization(bool silence, double dt)
{
    double max_length_colonization = age * getRootRandomParameter()->vi;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        if (isSecondaryColonizationAnchor(i)) {
            propagateSecondaryColonization(i, -1, max_length_colonization, silence, dt);
            propagateSecondaryColonization(i, 1, max_length_colonization, silence, dt);
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
                createHyphae(i);
                emergedHyphae.at(i) += 1; // increase the number of hyphae at this node
            }
        }
    } else { // Version where a set number of hyphae are created based on hyphal emergence density and root segment length
        auto rrp = getRootRandomParameter(); // param()
        double hed = rrp->hyphalEmergenceDensity;

        // double cumLength = 0.; // cumulative colonized length
        double numberOfHyphae= 0;

        for (size_t i = 0; i < nodes.size(); i++) {
            numberOfHyphae += emergedHyphae.at(i);
        }

        int new_noh = int(hed * getParameter("colonizationLength") - numberOfHyphae); // account for rounding errors
        if (hed * getParameter("colonizationLength") - numberOfHyphae - new_noh > 0.5) {
            new_noh += 1; // round up if the difference is larger than 0.5
        }
        double new_total_noh = numberOfHyphae + new_noh;

        int currentNode = 1;
        while (new_noh > 0 && numberOfHyphae < new_total_noh)
        {
            if (colonized.at(currentNode) > 0){ // if the current node is colonized and the number of hyphae to be created is reached
                createHyphae(currentNode);
                numberOfHyphae += 1;
                new_noh -= 1;
            }
            currentNode++;
            if (currentNode >= nodes.size() && new_noh > 0) {
                currentNode = 1; // reset to the first node if the end of the nodes vector is reached
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

        for (auto l : children)
        {
            if (l->organType()==Organism::ot_root) {
                if (colonized.at(l->parentNI) != 0) { // the base of root l is colonized
                    if (l->getNumberOfNodes() > 1 && std::dynamic_pointer_cast<MycorrhizalRoot>(l) -> getNodeColonization(1) == 0) {
                        std::dynamic_pointer_cast<MycorrhizalRoot>(l) ->setColonization(0, 3, 
								std::max( l->getNodeCT(0), colonizationTime.at(l->parentNI))); //root can only be colonized after it is born (and need to account for growth delay)
                    }
                }
                std::dynamic_pointer_cast<MycorrhizalRoot>(l) -> simulateColonization(dt, verbose);
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

double MycorrhizalRoot::getParameter(std::string name) const {
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
    return Root::getParameter(name);
}

void MycorrhizalRoot::setColonization(int i, int colonization, double t)
{
    colonized.at(i) = colonization;
    colonizationTime.at(i) = t;
	assert(colonizationTime.at(i) >= 0 && "MycorrhizalRoot::setColonization colonizationTime.at(i) < 0");
	assert(colonizationTime.at(i) >= nodeCTs.at(i) && "MycorrhizalRoot::setColonization colonizationTime.at(i) < nodeCTs.at(i)");
}

void MycorrhizalRoot::createLateral(double dt, bool verbose)
{
    auto rp = getOrganRandomParameter(); // rename

    for(int i = 0; i < rp->successorST.size(); i++){//go through each successor rule
        //found id
        bool applyHere = getApplyHere(i);

        if(applyHere)
        {
            int numlats = 1;//how many laterals? default = 1
            if(rp->successorNo.size()>i){numlats =  rp->successorNo.at(i);}
            for(int nn = 0; nn < numlats; nn++)
            {
                const Vector3d& pos = Vector3d();
                int p_id = rp->getLateralType(pos, i);//if probabilistic branching

                if(p_id >=0)
                {
                    int ot;

                    if((rp->successorOT.size()>i)&&(rp->successorOT.at(i).size()>p_id)){
                        ot = rp->successorOT.at(i).at(p_id);
                    }else{ot = getParameter("organType");}//default

                    int st = rp->successorST.at(i).at(p_id);

                    double delay = getLatGrowthDelay(ot, st, dt);// forDelay*multiplyDelay
                    double growth_dt = getLatInitialGrowth(dt);

                    switch(ot){
                    case Organism::ot_root:{
                        // std::cout << "Marco!" << std::endl;
                        auto lateral = std::make_shared<MycorrhizalRoot>(plant.lock(), st,  delay, shared_from_this(),  nodes.size() - 1);
                        // std::cout<< "Polo!"<< std::endl;
                        children.push_back(lateral);
                        lateral->simulate(growth_dt,verbose);
                        break;}
                    case Organism::ot_stem:{
                        auto lateral = std::make_shared<Stem>(plant.lock(), st, delay, shared_from_this(),  nodes.size() - 1);
                        children.push_back(lateral);
                        lateral->simulate(growth_dt,verbose);
                        break;}
                    case Organism::ot_leaf:{
                        auto lateral = std::make_shared<Leaf>(plant.lock(), st,  delay, shared_from_this(),  nodes.size() - 1);
                        children.push_back(lateral);
                        lateral->simulate(growth_dt,verbose);//age-ageLN,verbose);
                        break;}
                    }
                }
            }
        }
    }
    created_linking_node ++;
    storeLinkingNodeLocalId(created_linking_node,verbose);//needed (currently) only for stems when doing nodal growth
}

void MycorrhizalRoot::createHyphae(int pni)
{
    double delay = getRootRandomParameter()->hyphalDelay;
    double dt_ = plant.lock()->getSimTime() - colonizationTime.at(pni) - delay; // time the hyphae should have grown
	double delay_for_creation = colonizationTime.at(pni) + delay - shared_from_this()->getNodeCT(pni);//difference between creation of parent node and that of the hyphae
    assert(delay_for_creation >= 0 && "MycorrhizalRoot::createHyphae delay_for_creation < 0");
    int subType = 1;
    auto hyphae = std::make_shared<Hyphae>(plant.lock(), subType,  delay_for_creation, shared_from_this(), pni); // delay - dt_
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
