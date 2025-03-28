#include "MycorrhizalRoot.h"
#include "Root.h"
#include "Stem.h"
#include "Leaf.h"
#include "Organ.h"
#include "Organism.h"
#include "Hyphae.h"

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

void MycorrhizalRoot::primaryInfection(double dt, bool silence){
    double infTime;
    
    if (getRootRandomParameter()->infradius > 0) // check if localized infection should be applied
    {
        if (getRootRandomParameter()->f_inf)
        {
            for (size_t i = 1; i <nodes.size(); i++)
            {
                double p = getRootRandomParameter()->f_inf->getValue(nodes.at(i-1), shared_from_this());
                p = (1 - (age- nodeCTs.at(i-1))/getRootRandomParameter()->maxAge) * p;
                if (age - nodeCTs.at(i-1) < getRootRandomParameter() ->minAge) {p = 0;}
                double cursegLength = (nodes.at(i).minus(nodes.at(i-1))).length();
                infTime = plant.lock()->rand()*age + nodeCTs.at(i-1);
                if (infected.at(i-1) == 0 && plant.lock()->rand() < prob(infTime,cursegLength,p))
                {
                    setInfection(i-1,1,infTime); 
                }
            }
        }
    } else { //if this is not a loclized infection use equal probability everywhere
        for (size_t i=1; i<nodes.size(); i++) {
            double cursegLength = (nodes.at(i).minus(nodes.at(i-1))).length();
            infTime = plant.lock()->rand()*age + nodeCTs.at(i-1);
            double p = getRootRandomParameter()->p;
            p = (1 - (age- nodeCTs.at(i-1))/getRootRandomParameter()->maxAge) * p;
            if (age - nodeCTs.at(i-1) < getRootRandomParameter() ->minAge) {p = 0;}
            if ((plant.lock()->rand() <  prob(infTime,cursegLength,p) && (infected.at(i-1) == 0)))
            {
                setInfection(i-1,1,infTime); 
            }
        }
    }    
}

void MycorrhizalRoot::secondaryInfection(double maxLength, bool silence, double dt){
    for (size_t i = 0; i < nodes.size()-1; i++)
    {   
        if (infected.at(i) == 1 || infected.at(i)== 3)
        {
            int oldNode = i;
            double infTime;
            if (i>=1) {  // secondary infection in basal direction can only occur if there is another node in basal direction in this root
                int basalnode = i-1;
                while (basalnode >= 0 && basalnode < nodes.size()-1 && abs(nodes.at(i).minus(nodes.at(basalnode)).length()) < maxLength)
                {   
                    infTime = infectionTime.at(oldNode) + abs(nodes.at(basalnode).minus(nodes.at(oldNode)).length())/getRootRandomParameter()->vi; 
                    if (infected.at(basalnode) == 0 && infTime < age)
                    {
                        setInfection(basalnode,2,std::max(infTime,nodeCTs.at(basalnode)));
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

            while (apicalnode < nodes.size()-1 && abs(nodes.at(i).minus(nodes.at(apicalnode)).length()) < maxLength)
            {
                infTime = infectionTime.at(oldNode) + abs(nodes.at(oldNode).minus(nodes.at(apicalnode)).length())/getRootRandomParameter()->vi;
                if (infected.at(apicalnode) == 0 && infTime < age)
                {
                    setInfection(apicalnode,2,std::max(infTime,nodeCTs.at(apicalnode)));
                }
                oldNode = apicalnode;
                apicalnode++;
            }
            const MycorrhizalRootSpecificParameter& p = *param(); // rename
            if (alive) { // dead roots wont grow

        // // probabilistic branching model
        // if ((age>0) && (age-dt<=0)) { // the root emerges in this time step
        //     double P = getRootRandomParameter()->f_sbp->getValue(nodes.back(),shared_from_this());
        //     if (P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
        //         double p = 1.-std::pow((1.-P), dt); //probability of emergence in this time step
        //         if (plant.lock()->rand()>p) { // not rand()<p
        //             age -= dt; // the root does not emerge in this time step
        //         }
        //     }
        // }

        if (age>0) { // unborn  roots have no children

            // children first (lateral roots grow even if base root is inactive)
            for (auto l:children) {
                l->simulate(dt,verbose);
            }


            if (active) {

                // // length increment
                // double age_ = calcAge(length); // root age as if grown unimpeded (lower than real age)
                // double dt_; // time step
                // if (age<dt) { // the root emerged in this time step, adjust time step
                //     dt_= age;
                // } else {
                //     dt_=dt;
                // }

                // double targetlength = calcLength(age_+dt_)+ this->epsilonDx;

                // double e = targetlength-length; // unimpeded elongation in time step dt
                // double scale = getRootRandomParameter()->f_se->getValue(nodes.back(), shared_from_this());
                // double dl = std::max(scale*e, 0.);//  length increment = calculated length + increment from last time step too small to be added
                length = getLength();
                // this->epsilonDx = 0.; // now it is "spent" on targetlength (no need for -this->epsilonDx in the following)

                // create geometry
                if (length / getRootRandomParameter()->nEntryP > getRootRandomParameter()->d_e) { // root has children
                    double maxEntries = length / getRootRandomParameter()->d_e;
                    for (size_t i = 0; i < maxEntries; i++)
                    {
                        if (i == created_linking_node) {
                            createHyphae(dt,silence);
                        }
                    }
                
                        if ((p.ln.size()==created_linking_node)&& (getLength(true)-s>-1e-9)){
                            createLateral(dt_, verbose);
                        }
                    }

                    /* apical zone */
                    if (dl>0) {
                        createSegments(dl,dt_,verbose);
                        length+=dl; // - this->epsilonDx;
                    }
                } else { // no laterals

                    if (dl>0) {
                        createSegments(dl,dt_,verbose);
                        length+=dl; //- this->epsilonDx;
                    }
                } // if lateralgetLengths
            } // if active
            active = getLength(false)<=(p.getK()*(1 - 1e-11)); // become inactive, if final length is nearly reached
        }
    }
        }
    }

}

void MycorrhizalRoot::simulateInfection(double dt, bool verbose){
    if (this->nodes.size()>1) {
		//Primary Infection
        primaryInfection(dt,verbose);

        // Secondary Infection
        auto max_length_infection = age*getRootRandomParameter()->vi;

        secondaryInfection(max_length_infection,verbose,dt);
    
        for (auto l : children)
        {
            if (infected.at(l->parentNI) != 0)
            {
                if (l->getNumberOfNodes() > 1 && std::dynamic_pointer_cast<MycorrhizalRoot>(l) -> getNodeInfection(1) == 0)
                {
                    if(std::dynamic_pointer_cast<MycorrhizalRoot>(l)){
                        // std::cout<< "dynamic_cast successful!" <<std::endl;
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

void MycorrhizalRoot::createHyphae(double dt, bool silence)
{
    auto rp = getOrganRandomParameter(); // rename

	for(int i = 0; i < rp->successorST.size(); i++){//go through each successor rule
		//found id
		bool applyHere = getApplyHere(i);

		if(applyHere)
		{
			int numlats = 1;//how many laterals? default = 1
            int nEntryP = getRootRandomParameter()->nEntryP;
			// if(rp->successorNo.size()>i){numlats =  rp->successorNo.at(i);}
			for(int nn = 0; nn < numlats; nn++)
			{

				const Vector3d& pos = Vector3d();
				// int p_id = rp->getLateralType(pos, i);//if probabilistic branching

				// if(p_id >=0)
				// {
					int ot = getParameter("organType");

					// if((rp->successorOT.size()>i)&&(rp->successorOT.at(i).size()>p_id)){
					// 	ot = rp->successorOT.at(i).at(p_id);
					// }else{ot = getParameter("organType");}//default

					// int st = rp->successorST.at(i).at(p_id);
                    int type = 0;

					// double delay = getLatGrowthDelay(ot, st, dt);// forDelay*multiplyDelay
                    double delay = getRootRandomParameter()->hyphal_delay; // TODO make this parameter
					double growth_dt = getHyphalInitialGrowth(dt); // TODO add this function

					if (ot == Organism::ot_root) {
                        auto hyphae = std::make_shared<Hyphae>(plant.lock(), type,  delay, shared_from_this(),  nodes.size() - 1);
                        children.push_back(hyphae);
                        hyphae->simulate(growth_dt, silence);
                        nEntryP++;
                        getRootRandomParameter()->nEntryP = nEntryP;
                        break;
                    }
				// }
			}
		}

	}
	created_linking_node ++;
	storeLinkingNodeLocalId(created_linking_node,silence);//needed (currently) only for stems when doing nodal growth
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
double MycorrhizalRoot::prob(double  t, double segLength, double p)
{
    return 1 - pow(1-p,t*segLength);
}

}