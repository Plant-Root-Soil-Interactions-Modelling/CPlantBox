#include "MycorrhizalRoot.h"
#include "Root.h"
#include "Stem.h"
#include "Leaf.h"
#include "Organ.h"
#include "Organism.h"

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
            infectionTime.push_back(-1);

		}else{
			if ((parent->organType()==Organism::ot_stem)&&(parent->getNumberOfChildren()>0)) {
		}
            infected.push_back(0);
            infectionTime.push_back(-1);
		}
	}
    
    // std::cout<< "infected size " << infected.size() << std::endl;
    // std::cout<< "nodes size " << nodes.size() << std::endl;
}

void MycorrhizalRoot::addNode(Vector3d n, int id, double t, size_t index, bool shift) {
    Organ::addNode(n, id,  t,  index, shift);
    infected.push_back(0);
    infectionTime.push_back(-1);
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

void MycorrhizalRoot::simulateInfection(double dt, bool verbose){
    if (this->nodes.size()>1) {
        // std::cout << "\nstart Infection " << getId() <<  std::flush;
		//Primary Infection
        if (infected.size()!= nodes.size()) {std::cout<<"danger infection size not like node size!!!"<<std::endl;}
        double infTime;
        // get the soil from plant here or set it to some parameter
        // i.e. set a grid?? and then proceed as in the original code but with a if getDist <0 then localized infection
        if (getRootRandomParameter()->infradius > 0) // check if localized infection should be applied
        {
            Vector3d startPos = Vector3d(getRootRandomParameter()->posX, getRootRandomParameter()->posY, getRootRandomParameter()->posZ); // save the start position
            double infrad = getRootRandomParameter()->infradius; // easier access to the radius
            
            for (size_t i = 0; i < nodes.size()-1; i++)
            {
                if (startPos.minus(nodes.at(i)).length() < infrad && infected.at(i) == 0) // if within radius from start position then 100% gets infected
                {
                    infTime = plant.lock() ->rand()*dt;
                    setInfection(i,1,age + infTime);
                }
                // else if (plant.lock()->rand() < 0.01)
                // {
                //     infTime = plant.lock() ->rand()*dt;
                //     setInfection(i,1,age + infTime);
                // }
            }
        } else { //if this is not a loclized infection use equal probability everywhere
            for (size_t i=1; i<nodes.size(); i++) {
                double cursegLength = (nodes.at(i).minus(nodes.at(i-1))).length();
                if ((plant.lock()->rand() < 1 - pow(1-getRootRandomParameter()->p,dt*cursegLength) && (infected.at(i-1) == 0)))
                {
                    infTime = plant.lock() ->rand()*dt;
                    setInfection(i-1,1,age+infTime); 
                }
            }
        }
        // Secondary Infection
        auto max_length_infection = (age+dt)*getRootRandomParameter()->vi;
        
        for (size_t i = 0; i < nodes.size()-1; i++)
        {   
            if (infected.at(i) == 1 || infected.at(i)== 3)
            {
                
                auto max_length_basal = nodes.at(i).length() - max_length_infection;
                auto basalnode = i-1;
                double infTime;
                
                while (basalnode > 1 && basalnode< nodes.size()-1 && nodes.at(basalnode).length()> max_length_basal && infected.at(basalnode) == 0)
                {
                    
                    infTime = infectionTime.at(i) + nodes.at(i).minus(nodes.at(basalnode)).length()/getRootRandomParameter()->vi; 
                    if (infTime > nodeCTs.at(basalnode)){ // TODO make sure that new infected segments are not infected again!! also in apical direction
                        setInfection(basalnode,2,infTime);
                        if(basalnode==0) {
                            std::dynamic_pointer_cast<MycorrhizalRoot>(getParent()) ->setInfection(parentNI,3,infTime);
                            std::dynamic_pointer_cast<MycorrhizalRoot>(getParent()) ->simulateInfection(dt,verbose);
                        }
                    }
                    basalnode--;
                }

                auto max_length_apical = nodes.at(i).length() + max_length_infection;
                auto apicalnode = i+1;
                
                while (apicalnode < nodes.size()-1 && nodes.at(apicalnode).length() < max_length_apical && infected.at(apicalnode) == 0)
                {
                    infTime = infectionTime.at(i) + nodes.at(apicalnode).minus(nodes.at(i)).length()/getRootRandomParameter()->vi;
                    if (infTime > nodeCTs.at(apicalnode))
                    {
                        setInfection(apicalnode,2,infTime);
                    }
                    apicalnode++;
                }
            }
        }
        for (auto l : children)
        {
            // std::cout << "third infection" << std::endl;
            if (infected.at(l->parentNI) != 0)
            {
                if (l->getNumberOfNodes() > 1 && std::dynamic_pointer_cast<MycorrhizalRoot>(l) -> getNodeInfection(1) == 0)
                {
                    if(std::dynamic_pointer_cast<MycorrhizalRoot>(l)){
                        std::dynamic_pointer_cast<MycorrhizalRoot>(l) ->setInfection(0, 3, infectionTime.at(l->parentNI));
                    }
                    else std::cout<< "dynamic_cast failed!" <<std::endl;                
                }
            }
            std::dynamic_pointer_cast<MycorrhizalRoot>(l) -> simulateInfection(dt, verbose);
        }
    }
}

void MycorrhizalRoot::simulate(double dt, bool verbose)
{   
    // std::cout << "\nstart " << getId() <<  std::flush;
    Root::simulate(dt,verbose);
    simulateInfection(dt,verbose);
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
    // if (name == "infected") {return param() -> infected;}
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

std::string MycorrhizalRoot::toString() const
{
    // TODO this does not actually return the number of infected nodes fix this and add additional stuff
    std::stringstream newstring;
    newstring << "; infected Nodes " << getNumberofInfectedNodes() << ".";
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

}