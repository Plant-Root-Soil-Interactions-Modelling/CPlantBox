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
:Root(rs,type, delay,parent, pni) {}

void MycorrhizalRoot::addNode(Vector3d n, int id, double t, size_t index, bool shift) {
    // std::cout << "MycorrhizalRoot::addNode called" << std::endl;
    Organ::addNode(n, id,  t,  index, shift);
    infected.push_back(0);
    infectionTime.push_back(-1);
}


std::shared_ptr<Organ> MycorrhizalRoot::copy(std::shared_ptr<Organism> rs)
{
    // std::cout<< "MycorrhizalRoot::copy called" << std::endl;
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

void MycorrhizalRoot::simulate(double dt, bool verbose)
{   
    // std::cout << "MycorrhizalRoot::simulate" << std::endl;
    Root::simulate(dt,verbose);
 

    if (this->nodes.size()>1) {
		//Primary Infection
        if (getRootRandomParameter()->radius > 0) // check if localized infection should be applied
        {
            Vector3d startPos = Vector3d(getRootRandomParameter()->posX, getRootRandomParameter()->posY, getRootRandomParameter()->posZ); // save the start position
            double radius = getRootRandomParameter()->radius;
            double infectionage;
            for (size_t i = 0; i < nodes.size()-1; i++)
            {
                if (startPos.minus(nodes.at(i)).length() < radius && infected.at(i) == 0) // if within radius from start position then 100% gets infected
                {
                    setInfection(i,1,age + dt); // TODO this time stamp is not right yet
                }
                else if (plant.lock()->rand() < startPos.minus(nodes.at(i)).length()/radius) // TODO if not within radius probability decreases need to see how excactly
                {
                    infectionage= startPos.minus(nodes.at(i)).length()-radius;// TODO infection age not right right now
                    setInfection(i,1,age + infectionage);
                }
            }
        } else { //if this is not a loclized infection use equalprobability everywhere
            for (size_t i=1; i<nodes.size(); i++) {
                double cursegLength = (nodes.at(i).minus(nodes.at(i-1))).length();
                if ((plant.lock()->rand() < 1 - pow(1-getRootRandomParameter()->p,dt*cursegLength) && (infected.at(i-1) == 0)))
                {
                    setInfection(i-1,1,age + dt); // TODO age + dt does not make a whole lot of sense
                }
            }
        }
        

        // Secondary Infection
        auto max_length_infection = dt*getRootRandomParameter()->vi;
        for (size_t i = 0; i < nodes.size()-1; i++)
        {   
            if (infected.at(i) == 1)
            {
                auto max_length_basal = nodes.at(i).length() - max_length_infection;
                auto basalnode = i-1;
                double infectionage;
                
                while (basalnode > 1 && basalnode< nodes.size()-1 && nodes.at(basalnode).length()> max_length_basal && infected.at(basalnode) == 0)
                {
                    infectionage = age + nodes.at(i).minus(nodes.at(basalnode)).length()/getRootRandomParameter()->vi; 
                    setInfection(basalnode,2,infectionage); 
                    basalnode--;
                }

                auto max_length_apical = nodes.at(i).length() + max_length_infection;
                auto apicalnode = i+1;
                
                while (apicalnode < nodes.size()-1 && nodes.at(apicalnode).length() < max_length_apical && infected.at(apicalnode) == 0)
                {
                    infectionage = age + nodes.at(apicalnode).minus(nodes.at(i)).length()/getRootRandomParameter()->vi;
                    setInfection(apicalnode,2,infectionage);
                    apicalnode++;
                }
            }
        }
        
            for (auto l : children)
            {
                if (infected.at(l->parentNI) == 2)
                {
                    auto mnodes = std::dynamic_pointer_cast<MycorrhizalRoot>(l) -> getNodes();
                    if (mnodes.size() > 1 && std::dynamic_pointer_cast<MycorrhizalRoot>(l) -> getNodeInfection(0) == 0)
                    {
                        std::dynamic_pointer_cast<MycorrhizalRoot>(l) ->setInfection(0, 3, infectionTime.at(l->parentNI));
                    }
                }
            }
        
    }

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
    // std::cout << "MycorrhizalRoot::setInfection called" << std::endl;
    infected.at(i) = infection;
    infectionTime.at(i) = t;
}

void MycorrhizalRoot::createLateral(double dt, bool verbose)
{
    // std::cout << "MycorrhizalRoot::createLateral called" << std::endl;
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
							auto lateral = std::make_shared<MycorrhizalRoot>(plant.lock(), st,  delay, shared_from_this(),  nodes.size() - 1);
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
    newstring << "; infected Nodes " << infected.size() << ".";
    return  Root::toString()+newstring.str();
}

}