#include "Stem.h"

#include "Leaf.h"
#include "Root.h"
#include "Plant.h"
#include "Seed.h"
#include <algorithm>

namespace CPlantBox {

std::vector<int> Stem::phytomerId = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/**
 * Constructs a root from given data.
 * The organ tree must be created, @see Organ::setPlant, Organ::setParent, Organ::addChild
 * Organ geometry must be created, @see Organ::addNode, ensure that this->getNodeId(0) == parent->getNodeId(pni)
 *
 * @param id        		the organ's unique id (@see Organ::getId)
 * @param param     		the organs parameters set, ownership transfers to the organ
 * @param alive     		indicates if the organ is alive (@see Organ::isAlive)
 * @param active    		indicates if the organ is active (@see Organ::isActive)
 * @param age       		the current age of the organ (@see Organ::getAge)
 * @param length    		the current length of the organ (@see Organ::getLength)
 * @param partialIHeading 	the initial partial heading of this root
 * @param pbl       		base length of the parent root, where this root emerges
 * @param pni       		local node index, where this root emerges
 * @deprecated moved		indicates if nodes were moved in the previous time step (default = false)
 * @param oldNON    		the number of nodes of the previous time step (default = 0)
 */
Stem::Stem(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
		Vector3d partialIHeading_, int pni, bool moved, int oldNON)
:Organ(id, param, alive, active, age, length, partialIHeading_,  pni, moved,  oldNON)
{}

/**
 * Constructor
 * This is a Copy Paste of the Root.cpp but it works independently, it has its own parameter file (in .stparam file) tropism, growth function, txt and vtp writing system.
 * All of those can be modified to fit the real growth of the Plant.
 *
 * Typically called by the Plant::Plant(), or Stem::createNewStem().
 * For stem the initial node and node emergence time (netime) must be set from outside
 *
 * @param plant 		points to the plant
 * @param parent 		points to the parent organ
 * @param subtype		sub type of the stem
 * @param delay 		delay after which the organ starts to develop (days)
 * @param rheading		relative heading (within parent organ)
 * @param pni			parent node index
 * @param pbl			parent base length
 */
Stem::Stem(std::shared_ptr<Organism> plant, int type, double delay,  std::shared_ptr<Organ> parent, int pni)
:Organ(plant, parent, Organism::ot_stem, type, delay, pni)
{
	assert(parent!=nullptr && "Stem::Stem parent must be set");
	auto p = this->param();
	addPhytomerId(p->subType);
	double beta = getphytomerId(p->subType)*M_PI*getStemRandomParameter()->rotBeta +
			M_PI*plant->rand()*getStemRandomParameter()->betaDev;
	beta = beta + getStemRandomParameter()->initBeta*M_PI;
	if (getStemRandomParameter()->initBeta >0 && getphytomerId(p->subType)==0 ){
		beta = beta + getStemRandomParameter()->initBeta*M_PI;
	}
	double theta = p->theta;//M_PI*p->theta;
	if (parent->organType()!=Organism::ot_seed) { // scale if not a base organ, to delete?
		double scale = getStemRandomParameter()->f_sa->getValue(parent->getNode(pni), parent);
		theta *= scale;
	}
	//used when computing actual heading, @see Stem::getIHeading
	this->partialIHeading = Vector3d::rotAB(theta,beta);
	if (parent->organType()!=Organism::ot_seed) { // initial node
		//if lateral of stem, initial creation time: 
		//time when stem reached end of basal zone (==CT of parent node of first lateral) + delay
		// @see stem::createLateral
		double creationTime;
		if (parent->getNumberOfChildren() == 0){creationTime = parent->getNodeCT(pni)+delay;
		}else{creationTime = parent->getChild(0)->getParameter("creationTime") + delay;}
		
		Organ::addNode(Vector3d(0.,0.,0.), parent->getNodeId(pni), creationTime);//do not know why, but i have to add "Organ::" now
	}
}

/**
 * Deep copies the organ into the new plant @param rs.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
std::shared_ptr<Organ> Stem::copy(std::shared_ptr<Organism> p)
{
	auto s = std::make_shared<Stem>(*this); // shallow copy
	s->parent = std::weak_ptr<Organ>();
	s->plant = p;
	s->param_ = std::make_shared<StemSpecificParameter>(*param()); // copy parameters
	for (size_t i=0; i< children.size(); i++) {
		s->children[i] = children[i]->copy(p); // copy laterals
		s->children[i]->setParent(s);
	}
	return s;
}

/**
 * Simulates growth of this stem for a time span dt
 *
 * @param dt       time step [day]
 * @param verbose  indicates if status messages are written to the console (cout) (default = false)
 */
void Stem::simulate(double dt, bool verbose)
{
	if(!hasRelCoord()){
		throw std::runtime_error("organism no set in rel coord");
	}
	const StemSpecificParameter& p = *param(); // rename
	firstCall = true;
	oldNumberOfNodes = nodes.size();
	auto p_all = plant.lock();
	auto p_stem = p_all->getOrganRandomParameter(Organism::ot_stem);
	

	
	if (alive) { // dead roots wont grow

		// increase age
		if (age+dt>p.rlt) { // root life time
			dt=p.rlt-age; // remaining life span
			alive = false; // this root is dead
		}
		age+=dt;

		// probabilistic branching model (todo test)
		if ((age>0) && (age-dt<=0)) { // the root emerges in this time step
			//use relative coordinates for this function. Delete as it s not a root?
			double P = getStemRandomParameter()->f_sbp->getValue(nodes.back(),shared_from_this());
			if (P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
				double p = 1.-std::pow((1.-P), dt); //probability of emergence in this time step
				if (plant.lock()->rand()>p) { // not rand()<p
					age -= dt; // the root does not emerge in this time step
				}
			}
		}

		if (age>0) { // unborn  roots have no children

			// children first (lateral roots grow even if base root is inactive)
			for (auto l:children) {
				l->simulate(dt,verbose);
			}

			if (active) {
                double dt_; // time step
                if (age<dt) { // the root emerged in this time step, adjust time step
                    dt_= age;
                } else {
                    dt_=dt;
                }

				// length increment
				double age__ = age;
				if(age > p.delayNGStart){//simulation ends after start of growth pause
					if(age < p.delayNGEnd){age__ =p.delayNGStart;//during growth pause
					}else{
						age__ = age - (p.delayNGEnd - p.delayNGStart);//simulation ends after end of growth pause
					}
				}//delay to apply 
				/*as we currently do not implement impeded growth for stem and leaves
				*we can use directly the organ's age to cumpute the target length
				*/
				double targetlength = calcLength(age__)+ this->epsilonDx;
				double e = targetlength-length; // store value of elongation to add
				//can be negative
				double dl = e;//length increment = calculated length + increment from last time step too small to be added
				length = getLength(true);
				this->epsilonDx = 0.; // now it is "spent" on targetlength (no need for -this->epsilonDx in the following)
				// create geometry
				if (p.laterals) { // stem has laterals
					/* basal zone */
					if ((dl>0)&&(length<p.lb)) { // length is the current length of the root
						if (length+dl<=p.lb) {
							createSegments(dl,dt_,verbose);
							length+=dl;
							dl=0;
						} else {
							double ddx = p.lb-length;
							createSegments(ddx,dt_,verbose);
							dl-=ddx; // ddx already has been created
							length=p.lb;
						}
					}
					/* branching zone */
					//go into branching zone if organ has laterals and has reached 
					//the end of the basal zone
					if (((created_linking_node)<(p.ln.size()+1))&&(length>=p.lb)) 
					{
						for (size_t i=0; (i<p.ln.size()); i++) {
							createLateral(dt_, verbose);
							if(p.ln.at(created_linking_node-1)>0){
								createSegments(this->dxMin(),dt_,verbose);
								dl-=this->dxMin();
								length+=this->dxMin();
							}
						}
						createLateral(dt_, verbose);
					}
					//we can have (p.ln.size()+1)>(created_linking_node) if one ln == 0cm
					if((length>=p.lb)&&((p.ln.size()+1)<(created_linking_node))){
						std::stringstream errMsg;
						errMsg <<"Stem::simulate(): higher number of realized linking nodes ("<<created_linking_node<<
						") than of max laterals ("<<p.ln.size()+1<<")";
						throw std::runtime_error(errMsg.str().c_str());
					}
					//internodal elongation, if the basal zone of the stem is created and still has to grow
					double maxInternodeDistance = p.getK()-p.la - p.lb;//maximum length of branching zone
					if((dl>0)&&(length>=p.lb)&&(maxInternodeDistance>0)){
							int nn = localId_linking_nodes.back(); //node carrying the last lateral == end of branching zone
							double currentInternodeDistance = getLength(nn) - p.lb; //actual length of branching zone
							double ddx = std::min(maxInternodeDistance-currentInternodeDistance, dl);//length to add to branching zone 

							if(ddx > 0){
								internodalGrowth(ddx,dt_, verbose);
								dl -= ddx;
							length += ddx;
								
							}
						}
					/* apical zone */
					//only grows once the basal and branching nodes are developped
					if ((dl>0)&&(length-(maxInternodeDistance + p.lb)>-1e-9)) {
						createSegments(dl,dt_,verbose);
						length+=dl;
					} 
				} else { // no laterals
					if (dl>0) {
						createSegments(dl,dt_,verbose);
						length+=dl;
						
						}
				} // if lateralgetLengths
			if(dl <0){ //to keep in memory that realised length is too long, as created nodes to carry children
									
				this->epsilonDx = dl;//targetlength + e - length;
				length += this->epsilonDx;//go back to having length = theoratical length
			}
			} // if active
			//set limit below 1e-10, as the test files see if correct length 
			//once rounded at the 10th decimal
			//@see test/test_stem_ng.py
			active = getLength(false)<=(p.getK()*(1 - 1e-11)); // become inactive, if final length is nearly reached
		}
	} // if alive
}


/**
 *  @see Organ::createLateral
 *  @param dt       time step recieved by parent organ [day]	
 *  @return growth period to send to lateral after creation
 */
double Stem::getLatInitialGrowth(double dt)
{
	double ageLN = this->calcAge(param()->lb); // MINIMUM age of root when lateral node is created
    ageLN = std::max(ageLN, age-dt);
	return age-ageLN;
}


/**
 *  @see Organ::createLateral
 *  @param ot_lat       organType of lateral to create	
 *  @param st_lat       subType of lateral to create	
 *  @param dt       time step recieved by parent organ [day]	
 *  @return emergence delay to send to lateral after creation
 */
double Stem::getLatGrowthDelay(int ot_lat, int st_lat, double dt) const //override for stems
{
	
	bool verbose = false;
	auto rp = getOrganRandomParameter(); // rename
	double forDelay; //store necessary variables to define lateral growth delay
	int delayDefinition = std::static_pointer_cast<const SeedRandomParameter>(getOrganism()->getOrganRandomParameter(Organism::ot_seed,0))->delayDefinition;


	assert(delayDefinition >= 0);

			if(verbose){std::cout<<"create lat, delay def "<<delayDefinition<<" "
			<<getId()<<" "<< (nodes.size() - 1)<<" "<<age
			<<" "<<getNodeId(nodes.size() - 1)<<" "<<getNodeId(0)<<std::endl;
			}
	if(verbose){std::cout<<"create lat, delay def "<<delayDefinition<<std::endl;}
	//count the number of laterals of subtype st already created on this organ std::function<double(int, int, std::shared_ptr<Organ>)> 
	auto correctST = [ot_lat, st_lat](std::shared_ptr<Organ> org) -> double
		{
			return double((org->getParameter("organType") == ot_lat)&&(org->getParameter("subType")==st_lat));
		};//return 1. if organ of correct type and subtype, 0. otherwise
	
	double multiplyDelay = double(std::count_if(children.begin(), children.end(),
									 correctST));

	switch(delayDefinition){
		case Organism::dd_distance:
		{
			double meanLn = getParameter("lnMean"); // mean inter-lateral distance
			double effectiveLa = std::max(getParameter("la")-meanLn/2, 0.); // effective apical distance, observed apical distance is in [la-ln/2, la+ln/2]
			if(verbose)
			{
				std::cout<<"case Organism::dd_distance "<<organType()<<" "<<getParameter("subType")<<" "<<getLength(true)
				<<" "<<effectiveLa<<" "<<getParameter("la")<<" "<<meanLn<<std::endl;
			}
			double ageLN = this->calcAge(param()->lb); // age of root when lateral node is created
			ageLN = std::max(ageLN, age-dt);
			double ageLG = this->calcAge(param()->lb+effectiveLa); // age of the root, when the lateral starts growing (i.e when the apical zone is developed)
			forDelay = ageLG-ageLN; // time the lateral has to wait
			multiplyDelay = 1;//in this case, even for stems, it does not matter how many laterals there were before.
			break;
		}
		case Organism::dd_time_lat:
		{
			// time the lateral has to wait
			forDelay = std::max(rp->ldelay + plant.lock()->randn()*rp->ldelays, 0.);
			if(verbose){std::cout<<"Organism::dd_time_lat "<<rp->ldelay <<" "<<rp->ldelays<<" "<<forDelay<<std::endl;}
			break;
		}
		case Organism::dd_time_self:
		{
			
			//get delay per lateral
			auto latRp = plant.lock()->getOrganRandomParameter(ot_lat, st_lat); // random parameter of lateral to create
			forDelay = std::max(latRp->ldelay + plant.lock()->randn()*latRp->ldelays, 0.);
			if(verbose){
				std::cout<<"create lat, delay output "<<forDelay<<std::endl;
				std::cout<<"						 "<<ot_lat<<", "<<st_lat <<" "<<latRp->ldelay<<" "
				<< latRp->ldelays<<" "<<forDelay<<" "<<nodes.size()<<std::endl;
			}
			break;
		}
		default:
		{
			std::cout<<"delayDefinition "<<delayDefinition<<" "<<Organism::dd_distance<<" ";
			std::cout<< Organism::dd_time_lat<<" "<< Organism::dd_time_self<<std::endl<<std::flush;
			std::cout<<"				"<<(delayDefinition==Organism::dd_distance)<<" ";
			std::cout<<(delayDefinition== Organism::dd_time_lat)<<" "<< (delayDefinition==Organism::dd_time_self)<<std::endl<<std::flush;
			throw std::runtime_error("Delay definition type (delayDefinition) not recognised");
		}
	}
	if(verbose){std::cout<<"create lat, delay defEND "<<forDelay<<" "<<multiplyDelay<<std::endl;}
	return forDelay*multiplyDelay;
}
/**
 * Simulates internodal growth of dl for this stem
 * divid total stem growth between the phytomeres
 * currently two option:
 * growth devided equally between the phytomeres or 
 * the phytomere grow sequentially
 *
 * @param 	dl			total length of the segments that are created [cm]
 * @param	verbose		print information
 */
void Stem::internodalGrowth(double dl,double dt, bool verbose)
{
	const StemSpecificParameter& p = *param(); // rename
	std::vector<double> toGrow(p.ln.size());
	double dl_;
	const int ln_0 = std::count(p.ln.cbegin(), p.ln.cend(), 0);//number of laterals wich grow on smae branching point as the one before
	if(p.nodalGrowth==0){//sequentiall growth
		toGrow[0] = dl;
		std::fill(toGrow.begin()+1,toGrow.end(),0) ;
	}
	if(p.nodalGrowth ==1)
	{//equal growth
		std::fill(toGrow.begin(),toGrow.end(),dl/(p.ln.size()-ln_0)) ; 
	}
	int loopId = 0;
	size_t phytomerId = 0;
	while( (dl >0)&&(loopId<2) ) {//do the loop at most twice over the children
		//if the phytomere can do a growth superior to the mean phytomere growth, we add the value of "missing" 
		//(i.e., length left to grow to get the predefined total growth of the branching zone)
		int nn1 = localId_linking_nodes.at(phytomerId); //node at the beginning of phytomere		
		int nn2 = localId_linking_nodes.at(phytomerId+1); //node at end of phytomere (if nn1 != nn2)
		
		double length1 = getLength(nn1);
		double availableForGrowth = p.ln.at(phytomerId) -( getLength(nn2) - length1 ) ;//difference between maximum and current length of the phytomer
		if(availableForGrowth<-1e-3)
		{
			std::stringstream errMsg;
			errMsg <<"Stem::internodalGrowth phytomere "<<phytomerId<<" is too long: "<<availableForGrowth<<" "<<
			p.ln.at(phytomerId)<<" "<<getLength(nn2)<<" "<<length1<<std::endl;
			throw std::runtime_error(errMsg.str().c_str());
		}
		dl_ = std::max(0.,std::min(std::min(toGrow[phytomerId],availableForGrowth), dl));
		if(dl_ > 0)
		{
			createSegments(dl_,dt,verbose, nn2 ); dl -= dl_;
		}	
		if((phytomerId+1)< p.ln.size()){
			toGrow.at(phytomerId+1) +=  toGrow.at(phytomerId) - dl_ ;
			phytomerId ++;
		}else{
			toGrow.at(0) +=  toGrow.at(phytomerId) - dl_ ;
			loopId++; phytomerId = 0;
		}	//loop twice other the children
		
	}
	if(std::abs(dl)> 1e-6){//this sould not happen as computed dl to be <= sum(availableForGrowth)
		std::stringstream errMsg;
		errMsg <<"Stem::internodalGrowth length left to grow: "<<dl;
		throw std::runtime_error(errMsg.str().c_str());
	}
}
/**
 * Returns a parameter per organ
 *
 * @param name 		parameter name (returns nan if not available)
 *
 */
double Stem::getParameter(std::string name) const
{
	if (name=="lb") { return param()->lb; } // basal zone [cm]
	if (name=="delayNGStart") { return param()->delayNGStart; } // delay for nodal growth [day]
	if (name=="delayNGEnd") { return param()->delayNGEnd; } // delay for nodal growth [day]
	if (name=="la") { return param()->la; } // apical zone [cm]
	if (name=="nob") { return param()->nob(); } // number of branching points
	if (name=="r"){ return param()->r; }  // initial growth rate [cm day-1]
	if (name=="radius") { return param()->a; } // root radius [cm]
	if (name=="a") { return param()->a; } // root radius [cm]
	if (name=="theta") { return param()->theta; } // angle between root and parent root [rad]
	if (name=="rlt") { return param()->rlt; } // root life time [day]
	if (name=="k") { return param()->getK(); }; // maximal root length [cm]
	if (name=="lnMean") { // mean lateral distance [cm]
        auto& v =param()->ln;
		if(v.size()>0){
			return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
		}else{
			return 0;
		}
	}
	if (name=="lnDev") { // standard deviation of lateral distance [cm]
		auto& v =param()->ln;
		double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
		double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
		return std::sqrt(sq_sum / v.size() - mean * mean);
	}
	if (name=="volume") { return param()->a*param()->a*M_PI*getLength(true); } // // root volume [cm^3]
	if (name=="surface") { return 2*param()->a*M_PI*getLength(true); }
	if (name=="type") { return this->param_->subType; }  // delete to avoid confusion?
	if (name=="subType") { return this->param_->subType; }  // organ sub-type [-]
	if (name=="parentNI") { return parentNI; } // local parent node index where the lateral emerges
	return Organ::getParameter(name);
}




/**
 * Analytical length of the stem at a given age
 *
 * @param age          age of the stem [day]
 */
double Stem::calcLength(double age)
{
	assert(age>=0 && "Stem::calcLength() negative root age");
	return getStemRandomParameter()->f_gf->getLength(age,getStemRandomParameter()->r,param()->getK(),shared_from_this());
}

/**
 * Analytical age of the stem at a given length
 * no scaling of organ growth , so can return age directly
 * otherwise cannot compute exact age between delayNGStart and delayNGEnd
 * @param length   length of the stem [cm]
 */
double Stem::calcAge(double length) const
{
	assert(length>=0 && "Stem::calcAge() negative root age");
	double age__ = getStemRandomParameter()->f_gf->getAge(length,getStemRandomParameter()->r,param()->getK(),shared_from_this());
	if(age__ >param()->delayNGStart ){age__ += (param()->delayNGEnd - param()->delayNGStart);}
	return age__;
}


/**
 * stores the local id of the linking node. used by @see Stem::internodalGrowth()
 */
void Stem::storeLinkingNodeLocalId(int numCreatedLN, bool verbose)
{
	localId_linking_nodes.push_back(nodes.size()-1);
	if(numCreatedLN!=localId_linking_nodes.size())
	{
		throw std::runtime_error("wrong number of linking nodes in stem: "+std::to_string(numCreatedLN)
		+" against "+std::to_string(localId_linking_nodes.size()));
	}
	if(verbose)
	{
		std::cout<<"Stem::storeLinkingNodeLocalId "<<numCreatedLN<<" "<<(nodes.size()-1)<<" "<<localId_linking_nodes.size()<<std::endl;
	}
}

/**
 * Adds a node to the organ.
 *
 * For simplicity nodes can not be deleted, organs can only become deactivated or die
 *
 * @param n        new node
 * @param id       global node index
 * @param t        exact creation time of the node
 * @param index	   position were new node is to be added
 * @param shift	   do we need to shift the nodes? (i.e., is the new node inserted between existing nodes because of internodal growth?)
 */
void Stem::addNode(Vector3d n, int id, double t, size_t index, bool shift)
{
	bool verbose = false;
	if(verbose)
	{
		std::cout<<"Organ::addNode "<<id<<" "<<getId()<<" "<<organType()<<" "<<getParameter("subType")<<std::endl;
		std::cout<<"Organ::addNode "<<n.toString()<<" "<<t<<" "<<index<<" "<<shift<<std::endl;
		
	}
	if(!shift){//node added at the end of organ
		nodes.push_back(n); // node
		nodeIds.push_back(id); //unique id
		nodeCTs.push_back(t); // exact creation time
	}
	else{//could be quite slow  to insert, but we won t have that many (node-)tillers (?)
		nodes.insert(nodes.begin() + index-1, n);//add the node at index
		//add a global index.
		//no need for the nodes to keep the same global index and makes the update of the nodes position for MappedPlant object more simple)
		//if(verbose){
			//			std::cout<<"Organ::addNode "<<organType()<<" "<<id<<" "<<index<<std::endl<<std::flush;
		//}
		nodeIds.push_back(id);
		nodeCTs.insert(nodeCTs.begin() + index-1, t);
		for(auto kid : children){//if carries children after the added node, update their "parent node index"
		
			if((kid->parentNI >= index-1 )&&(kid->parentNI > 0)){
				kid->moveOrigin(kid->parentNI + 1);
				}

		}
		for(int numnode = 0; numnode < localId_linking_nodes.size();numnode++){//update the local ids of the linking nodes
			if((localId_linking_nodes.at(numnode) >= index-1 )&&(localId_linking_nodes.at(numnode) > 0))
			{
				localId_linking_nodes.at(numnode) += 1;
			}
		}

	}
}


/**
 * @return The StemTypeParameter from the plant
 */
std::shared_ptr<StemRandomParameter> Stem::getStemRandomParameter() const
{
	return std::static_pointer_cast<StemRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_stem, param_->subType));
}

/**
 * @return Parameters of the specific root
 */
std::shared_ptr<const StemSpecificParameter> Stem::param() const
{
	return std::static_pointer_cast<const StemSpecificParameter>(param_);
}

/*
 * Quick info about the object for debugging
 * additionally, use param()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Stem::toString() const
{
	std::stringstream newstring;
	newstring << "; initial heading: " << getiHeading0().toString() << ", parent node index" << parentNI << ".";
	return Organ::toString()+newstring.str();
}

} // namespace CPlantBox
