#include "MycorrhizalRoot.h"
#include "Root.h"
#include "Stem.h"
#include "Leaf.h"
#include "Organ.h"
#include "Organism.h"

#include "Hyphae.h"
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
            infected.push_back(0);
            emergedHyphae.push_back(0);
            infectionTime.push_back(-1);

        }else{
            if ((parent->organType()==Organism::ot_stem)&&(parent->getNumberOfChildren()>0)) {
            }
            infected.push_back(0);
            emergedHyphae.push_back(0);
            infectionTime.push_back(-1);
        }
    }
    // std::cout<< "infected size " << infected.size() << std::endl;
    // std::cout<< "nodes size " << nodes.size() << std::endl;
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
        infected.push_back(0);
        emergedHyphae.push_back(0);
        infectionTime.push_back(-1);
    } 
    else {
        //Organ::addNode(n, id,  t,  index, shift);
		nodes.insert(nodes.begin() + index-1, n);//add the node at index
		nodeIds.push_back(id);
		nodeCTs.insert(nodeCTs.begin() + index-1, t);
        infected.insert(infected.begin()+index-1, infected.at(index-1));
        emergedHyphae.insert(emergedHyphae.begin()+index-1, 0);
        infectionTime.insert(infectionTime.begin()+index-1, infectionTime.at(index-1));

        for(auto kid : children){//if carries children after the added node, update their "parent node index"

			if((kid->parentNI >= index-1 )&&(kid->parentNI > 0)){
				kid->moveOrigin(kid->parentNI + 1);
			}
		}
    }
    // std::cout<< "infected size " << infected.size() << std::endl;
    // std::cout<< "nodes size " << nodes.size() << std::endl;
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

void MycorrhizalRoot::primaryInfection(double dt, bool silence){
    double lmbd;
    double highres = getRootRandomParameter()->highresolution;
    for (size_t i = 1; i < nodes.size(); i++){
        if (getRootRandomParameter()->f_inf->getValue(nodes.at(i), shared_from_this()) != 1.)
        {
            lmbd = getRootRandomParameter()->f_inf->getValue(nodes.at(i), shared_from_this());
            // std::cout << "MycorrhizalRoot::primaryInfection(): Infection rate at node " << i << ": " << lmbd << std::endl;
        } 
        else {
            lmbd = getRootRandomParameter()->lmbd;
        }
        if (age - nodeCTs.at(i) < getRootRandomParameter() ->minAge) {lmbd = 0;}//account for minimal age in rate
        lmbd = (1 - (age- nodeCTs.at(i))/getRootRandomParameter()->maxAge) * lmbd; // account for maximal age in rate
        double cursegLength = (nodes.at(i).minus(nodes.at(i-1))).length();
        if (infected.at(i) == 0 && plant.lock()->rand() < lmbd*cursegLength*dt)
        {
            // insert node here if segment too long and set all nodes to be infected
            setInfection(i,1,age);
            if (highres >= 1. && cursegLength > getRootRandomParameter() ->dx_inf) {
                int newNodesNumber = std::max( int(cursegLength / getRootRandomParameter() ->dx_inf) - 1, 0);
                for (size_t j = 0; j < newNodesNumber; j++)
                {
                    double newx = nodes.at(i-1).x + (nodes.at(i).x - nodes.at(i-1).x) *(j+1)/(newNodesNumber +1) ;
                    double newy = nodes.at(i-1).y + (nodes.at(i).y - nodes.at(i-1).y) *(j+1)/(newNodesNumber +1);
                    double newz = nodes.at(i-1).z + (nodes.at(i).z - nodes.at(i-1).z) *(j+1)/(newNodesNumber +1);
                    Vector3d newNode = Vector3d(newx,newy,newz);
                    //std::cout << "inserting node at index " << i << " at position " << newNode.toString() << "\n" << "Current node size: " << nodes.size() << std::endl;
                    addNode(newNode,plant.lock()->getNodeIndex(), nodeCTs.at(i), i, true);
                }   
            }
        }
    }
}

void MycorrhizalRoot::secondaryInfection(bool silence, double dt){

    double max_length_infection = age*getRootRandomParameter()->vi;

    double infTime;

    double highres = getRootRandomParameter()->highresolution;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        if (infected.at(i) == 1 || infected.at(i)== 3)
        {
            int oldNode = i;
            double infectionLength = 0;
            if (i>=1) {  // secondary infection in basal direction can only occur if there is another node in basal direction in this root
                int basalnode = i-1;
                double cursegLength;
                // std::cout << infectionLength << std::endl;
                while(basalnode > 0) {
                    // std::cout << "basalnode " << basalnode << std::endl;
                    cursegLength = abs(nodes.at(oldNode).minus(nodes.at(basalnode)).length());
                    infectionLength += cursegLength;
                    infTime = infectionTime.at(oldNode) + cursegLength/getRootRandomParameter()->vi;

                    if (infectionLength > max_length_infection) {break;}

                    if (infected.at(basalnode) == 0 && infTime <= age)
                    {
                        // insert node here if segment too long and set all nodes to be infected
                        setInfection(basalnode,2,infTime);
                        if (highres >= 1. && cursegLength > getRootRandomParameter() ->dx_inf) {
                            int newNodesNumber = std::max( int(cursegLength / getRootRandomParameter() ->dx_inf) - 1, 0);
                            for (size_t j = 0; j < newNodesNumber; j++)
                            {
                                double newx = nodes.at(oldNode).x + (nodes.at(basalnode).x - nodes.at(oldNode).x) *(j+1)/(newNodesNumber +1);
                                double newy = nodes.at(oldNode).y + (nodes.at(basalnode).y - nodes.at(oldNode).y) *(j+1)/(newNodesNumber +1);
                                double newz = nodes.at(oldNode).z + (nodes.at(basalnode).z - nodes.at(oldNode).z) *(j+1)/(newNodesNumber +1);
                                Vector3d newNode = Vector3d(newx,newy,newz);
                                // infTime = infectionTime.at(oldNode) + abs(nodes.at(oldNode).minus(newNode).length())/getRootRandomParameter()->vi;
                                // std::cout << "inserting node at index " << basalnode << " at position " << newNode.toString() << "\n" << "Current node size: " << nodes.size() << std::endl;
                                addNode(newNode,plant.lock()->getNodeIndex(), nodeCTs.at(basalnode), basalnode, true);
                                // infectionTime.at(basalnode) = infTime;
                            }   
                        }
                        // std::cout<< "secondary infection from " << i << " to " << basalnode << std::endl;
                        if(basalnode==0 && std::dynamic_pointer_cast<MycorrhizalRoot>(getParent()))
                        {
                            // std::cout << "basalnode is 0" << std::endl;
                            std::dynamic_pointer_cast<MycorrhizalRoot>(getParent())->setInfection(parentNI,3,infTime);
                            std::dynamic_pointer_cast<MycorrhizalRoot>(getParent())->simulateInfection(dt,silence);
                        }
                    }
                    
                    oldNode = basalnode;
                    basalnode--;
                }
            }

            auto apicalnode = i+1;
            oldNode = i;
            infectionLength = 0;
            double cursegLength;
            while (apicalnode < nodes.size())
            {
                cursegLength = abs(nodes.at(oldNode).minus(nodes.at(apicalnode)).length());
                infectionLength += cursegLength;
                infTime = infectionTime.at(oldNode) + cursegLength/getRootRandomParameter()->vi;
                if (infectionLength > max_length_infection) {break;}
                if (infected.at(apicalnode) == 0 && infTime <= age)
                {
                    // insert node here if segment too long and set all nodes to be infected
                    setInfection(apicalnode,2,infTime);
                    if (highres >= 1. && cursegLength > getRootRandomParameter() ->dx_inf) {
                            int newNodesNumber = std::max( int(cursegLength / getRootRandomParameter() ->dx_inf) - 1, 0);
                            for (size_t j = 0; j < newNodesNumber; j++)
                            {
                                double newx = nodes.at(oldNode).x + (nodes.at(apicalnode).x - nodes.at(oldNode).x) *(j+1)/(newNodesNumber +1);
                                double newy = nodes.at(oldNode).y + (nodes.at(apicalnode).y - nodes.at(oldNode).y) *(j+1)/(newNodesNumber +1);
                                double newz = nodes.at(oldNode).z + (nodes.at(apicalnode).z - nodes.at(oldNode).z) *(j+1)/(newNodesNumber +1);
                                Vector3d newNode = Vector3d(newx,newy,newz);
                                infTime = infectionTime.at(oldNode) + abs(nodes.at(oldNode).minus(newNode).length())/getRootRandomParameter()->vi;
                                //std::cout << "inserting node at index " << i << " at position " << newNode.toString() << "\n" << "Current node size: " << nodes.size() << std::endl;
                                addNode(newNode,plant.lock()->getNodeIndex(), nodeCTs.at(apicalnode), apicalnode, true);
                            }   
                    }
                }    
                oldNode = apicalnode;
                apicalnode++;
            }
        }
    }
}

void MycorrhizalRoot::simulateSecondaryInfection(double dt) {
    secondaryInfection(false,dt);
}

void MycorrhizalRoot::simulatePrimaryInfection(double dt) {
    primaryInfection(dt,false);
}

void MycorrhizalRoot::simulateHyphalGrowth() {
    if (getRootRandomParameter()->highresolution >= 1) { // Version where at every node there is one hypha created
        for (size_t i = 1; i < nodes.size(); i++) {
            if (infected.at(i) > 0 && emergedHyphae.at(i) == 0){ // if the current node is infected and the number of hyphae to be created is reached
                createHyphae(i);
                emergedHyphae.at(i) += 1; // increase the number of hyphae at this node
            }
        }
    } else { // Version where a set number of hyphae are created based on hyphal emergence density and root segment length
        auto rrp = getRootRandomParameter(); // param()
        double hed = rrp->hyphalEmergenceDensity;

        // double cumLength = 0.; // cumulative infected length
        double numberOfHyphae= 0;

        for (size_t i = 0; i < nodes.size(); i++) {
            numberOfHyphae += emergedHyphae.at(i);
        }

        int new_noh = int(hed * getParameter("infectionLength") - numberOfHyphae); // account for rounding errors
        if (hed * getParameter("infectionLength") - numberOfHyphae - new_noh > 0.5) {
            new_noh += 1; // round up if the difference is larger than 0.5
        }
        double new_total_noh = numberOfHyphae + new_noh;
        // std::cout << "MycorrhizalRoot::simulateHyphalGrowth(): " << "Hyphal Emergence density " << hed << ", infectionLength:" << getParameter("infectionLength") << ", noh " << numberOfHyphae <<  ", new noh " << new_noh << std::endl;

        int currentNode = 1;
        while (new_noh > 0 && numberOfHyphae < new_total_noh) 
        {
            if (infected.at(currentNode) > 0){ // if the current node is infected and the number of hyphae to be created is reached
                createHyphae(currentNode);
                numberOfHyphae += 1;
                new_noh -= 1;
                // lastEmergedNode = currentNode; // update the last emerged node
            }
            currentNode++;
            if (currentNode >= nodes.size() && new_noh > 0) {
                currentNode = 1; // reset to the first node if the end of the nodes vector is reached
            }
        }
    }
    
}



void MycorrhizalRoot::simulateInfection(double dt, bool verbose) {

    if (this->nodes.size()>1) {

        //Primary Infection
        primaryInfection(dt,verbose);

        // Secondary Infection
        secondaryInfection(verbose,dt);

        for (auto l : children)
        {
            if (l->organType()==Organism::ot_root) {

                if (infected.at(l->parentNI) != 0) { // the base of root l is infected

                    if (l->getNumberOfNodes() > 1 && std::dynamic_pointer_cast<MycorrhizalRoot>(l) -> getNodeInfection(1) == 0) {

                        std::dynamic_pointer_cast<MycorrhizalRoot>(l) ->setInfection(0, 3, infectionTime.at(l->parentNI));
                    }
                }
                std::dynamic_pointer_cast<MycorrhizalRoot>(l) -> simulateInfection(dt, verbose);
            }
        }
    }
}


void MycorrhizalRoot::simulate(double dt, bool verbose)
{
    // std::cout << "\nstart " << getId() <<  std::flush;
    Root::simulate(dt,verbose);
    simulateInfection(dt,verbose);
    simulateHyphalGrowth();
    // std::cout << getRootRandomParameter()->la << std::endl;
    // std::cout << "\nend " << getId() <<  std::flush;
}

std::shared_ptr<const MycorrhizalRootSpecificParameter> MycorrhizalRoot::param() const
{
    // std::cout << "MycorrhizalRoot::param called" << std::endl;
    return std::static_pointer_cast<const MycorrhizalRootSpecificParameter>(param_);
}

std::shared_ptr<MycorrhizalRootRandomParameter> MycorrhizalRoot::getRootRandomParameter() const
{
    // std::cout << "MycorrhizalRoot::getRootRandomParameter called" << std::endl;
    return std::static_pointer_cast<MycorrhizalRootRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_root, param_->subType));
}

double MycorrhizalRoot::getParameter(std::string name) const {
    // std::cout << "MycorrhizalRoot::getParameter called" << std::endl;
    if (name == "primaryInfection") 
    {
        double primaryInfectedLength = 0;
        for (size_t i = 1; i < nodes.size(); i++)
        {
            if (infected.at(i)==1)
            {
                primaryInfectedLength += nodes.at(i).minus(nodes.at(i-1)).length();
            }
        }
        return primaryInfectedLength;
    }
    if (name == "secondaryInfection")
    {
        double secondaryInfectedLength = 0;
        for (size_t i = 1; i < nodes.size(); i++)
        {
            if (infected.at(i)>1)
            {
                secondaryInfectedLength += nodes.at(i).minus(nodes.at(i-1)).length();
            }
        }
        return secondaryInfectedLength;
    }
    if (name == "infectionLength") 
    {
        double infectedLength = 0;
        for (size_t i = 1; i < nodes.size(); i++)
        {
            if (infected.at(i)>0)
            {
                infectedLength += nodes.at(i).minus(nodes.at(i-1)).length();
            }
        }
        return infectedLength;
    }
    
    return Root::getParameter(name);
}

void MycorrhizalRoot::setInfection(int i, int infection, double t)
{
    // std::cout << i << " " << infection << " " << t << std::endl;
    infected.at(i) = infection;
    infectionTime.at(i) = t;
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
    double dt_ = plant.lock()->getSimTime() - infectionTime.at(pni); // time the hyphae should have grown
    double delay = getRootRandomParameter()->hyphalDelay; // todo specific (with std)
    int subType = 1;
    auto hyphae = std::make_shared<Hyphae>(plant.lock(), subType,  delay, shared_from_this(), pni); // delay - dt_
    children.push_back(hyphae);
    emergedHyphae.at(pni) += 1;
    // std::cout << "********* simulate "  << ", "<< plant.lock()->getSimTime() <<", " << dt_ << "\n";
    hyphae->simulate(dt_);
}

std::string MycorrhizalRoot::toString() const
{
    // TODO this does not actually return the number of infected nodes fix this and add additional stuff
    std::stringstream newstring;
    newstring << "; number of infected Nodes " << getNumberofInfectedNodes() << "; length of infected root segments "<< getParameter("infectionLength")<< ".";
    return  Root::toString()+newstring.str();
}

int MycorrhizalRoot::getNumberofInfectedNodes() const
{
    int numberInfectedNodes =0;
    for (size_t i = 0; i < getNumberOfNodes()-1; i++)
    {
        if (infected.at(i)!=0)
        {
            numberInfectedNodes++;
        }
    }
    return numberInfectedNodes;
}


double MycorrhizalRoot::prob(double t, double segLength, double p)
{
    return 1 - pow(1-p,t*segLength);
}

}
