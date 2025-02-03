#include "MycorrhizalRoot.h"
#include "Root.h"
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
    Organ::addNode(n, id,  t,  index, shift);
    infected.push_back(0);
    infectionTime.push_back(-1);
}



void MycorrhizalRoot::simulate(double dt, bool verbose)
{
    Root::simulate(dt,verbose);
 

    if (this->nodes.size()>1) {
		//Primary Infection
		for (size_t i=1; i<nodes.size(); i++) {
			
            double cursegLength = (nodes.at(i).minus(nodes.at(i-1))).length();
            if ((plant.lock()->rand() < (getRootRandomParameter()->p*dt*cursegLength)) && (infected.at(i-1) == 0))
            {
                infected.at(i-1) = 1;
                infectionTime.at(i-1) = age + dt;
                
		    }
	    }

        //Secondary Infection
        auto max_length_infection = dt*getRootRandomParameter()->vi;
        for (size_t i = 1; i < nodes.size(); i++)
        {
            if (infected.at(i-1) == 1)
            {
                auto max_length_basal = nodes.at(i-1).length() - max_length_infection;
                auto basalnode = i-2;
                
                while (basalnode > 1 && basalnode< nodes.size() && nodes.at(basalnode).length()> max_length_basal && infected.at(basalnode) == 0)
                {
                    infected.at(basalnode) = 2;
                    infectionTime.at(basalnode) = age + nodes.at(i-1).minus(nodes.at(basalnode)).length()/getRootRandomParameter()->vi; 
                    basalnode--;
                }

                auto max_length_apical = nodes.at(i-1).length() + max_length_infection;
                auto apicalnode = i;
                
                while (apicalnode < nodes.size()-1 && nodes.at(apicalnode).length() < max_length_apical && infected.at(apicalnode) == 0)
                {
                    infected.at(apicalnode) = 2;
                    infectionTime.at(apicalnode) = age + nodes.at(apicalnode).minus(nodes.at(i-1)).length()/getRootRandomParameter()->vi; 
                    apicalnode++;
                }
            }
        }
        
        for (auto l:children)
        {
            for (size_t i = 1; i < nodes.size(); i++)
            {
                if (getNodeId(i-1) == getNodeId(l->parentNI) && infected.at(i-1) == 2)
                {
                     double infectionage = age + nodes.at(i).minus(nodes.at(i-1)).length()/getRootRandomParameter()->vi;
                     std::dynamic_pointer_cast<MycorrhizalRoot>(l) -> setInfection(0,3,infectionage);
                }
            }
            
            
        }
        
    }

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
    // if (name == "infected") {return param() -> infected;}
    return Root::getParameter(name);
}

void MycorrhizalRoot::setInfection(int i, int infection, double t)
{
    infected.at(i) = infection;
    infectionTime.at(i) = t;
}


}