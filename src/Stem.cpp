#include "Stem.h"

#include "Leaf.h"
#include "Root.h"
#include "Plant.h"

namespace CPlantBox {

std::vector<int> Stem::phytomerId = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/**
 * Constructs a organ from given data.
 * The organ tree must be created, @see Organ::setPlant, Organ::setParent, Organ::addChild
 * Organ geometry must be created, @see Organ::addNode, ensure that this->getNodeId(0) == parent->getNodeId(pni)
 *
 * @param id        	the organ's unique id (@see Organ::getId)
 * @param param     	the organs parameters set, ownership transfers to the organ
 * @param alive     	indicates if the organ is alive (@see Organ::isAlive)
 * @param active    	indicates if the organ is active (@see Organ::isActive)
 * @param age       	the current age of the organ (@see Organ::getAge)
 * @param length    	the current length of the organ (@see Organ::getLength)
 * @param iheading  	the initial heading of this organ
 * @param pbl       	base length of the parent organ, where this organ emerges
 * @param pni       	local node index, where this organ emerges
 * @deprecated moved	indicates if nodes were moved in the previous time step (default = false)
 * @param oldNON    	the number of nodes of the previous time step (default = 0)
 */
Stem::Stem(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
		Matrix3d iHeading, int pni, bool moved, int oldNON)
:Organ(id, param, alive, active, age, length, iHeading, pni, moved,  oldNON)
{}

/**
 * Constructor
 * This is a Copy Paste of the organ.cpp but it works independently, it has its own parameter file (in .stparam file) tropism, growth function, txt and vtp writing system.
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
Stem::Stem(std::shared_ptr<Organism> plant, int type, Matrix3d iHeading, double delay,  std::shared_ptr<Organ> parent, int pni)
:Organ(plant, parent, Organism::ot_stem, type, delay, iHeading, pni)
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
		
		addNode(Vector3d(0.,0.,0.), parent->getNodeId(pni), creationTime);
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
	const StemSpecificParameter& p = *param(); // rename
	firstCall = true;
	oldNumberOfNodes = nodes.size();
	auto p_all = plant.lock();
	auto p_stem = p_all->getOrganRandomParameter(Organism::ot_stem);

	int nC = getPlant()->getSeed()->param()->nC; //number of the shoot born root
	double nZ = getPlant()->getSeed()->param()->nz; // distance between shoot born root and the seed
	double res = nZ -floor(nZ / dx())*dx();
	if(res < dxMin() && res != 0){
		if(res <= dxMin()/2){ nZ -= res;
		}else{nZ =  floor(nZ / dx())*dx() + dxMin();}
		if(verbose){std::cout<<"\nStem::simulate: nZ changed to "<<nZ<<" for compatibility with dx and dxMin"<<std::endl;}
	}			//make nZ compatible with dx() and dxMin()


	int additional_childern;
	if (p.subType == 1)
	{
		additional_childern= (int)round(nC); //if it is the main stem, the children should include the shoot borne root
	} else {
		additional_childern = 0;
	}

	if (alive) { // dead organs wont grow

		// increase age
		if (age+dt>p.rlt) { // organ life time
			dt=p.rlt-age; // remaining life span
			alive = false; // this organ is dead
		}
		age+=dt;

		// probabilistic branching model (todo test)
		if ((age>0) && (age-dt<=0)) { // the organ emerges in this time step
			//use relative coordinates for this function. Delete as it s not a root?
			double P = getStemRandomParameter()->f_sbp->getValue(nodes.back(),shared_from_this());
			if (P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
				double p = 1.-std::pow((1.-P), dt); //probability of emergence in this time step
				if (plant.lock()->rand()>p) { // not rand()<p
					age -= dt; // the organ does not emerge in this time step
				}
			}
		}

		if (age>0) { // unborn  organs have no children

			// children first (lateral organs grow even if base organ is inactive)
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
				double delayNG_ = 0.; //delay between organ creation and expansion
				double ageAtLb = calcAge(p.lb);//age at which organ reached lb
				if((age > ageAtLb)){
					delayNG_ = std::min(p.delayNG, age - ageAtLb );
				}//delay to apply 
				/*as we currently do not implement impeded growth for stem and leaves
				*we can use directly the organ's age to cumpute the target length
				*/
				double targetlength = calcLength(age - delayNG_)+ this->epsilonDx;
				double e = targetlength-length; // store value of elongation to add
				//can be negative
				double dl = e;//length increment = calculated length + increment from last time step too small to be added
				length = getLength();
				this->epsilonDx = 0.; // now it is "spent" on targetlength (no need for -this->epsilonDx in the following)
				// create geometry
				if (p.laterals) { // stem has laterals
					//std::cout<<"sim seed nC is"<< nC<<"\n";
					//std::cout<<"sim seed nZ is"<< nZ<<"\n";
					/*
                    shoot born root
					 */
					if ((dl>0)&&(length< nZ)) {
						if (length+dl <= nZ) {
							createSegments(dl,dt_,verbose);
							length+=dl ;
							dl=0;
						} else {
							double ddx = nZ - length;
							createSegments(ddx,dt_,verbose);//should it not be ddx here?

							dl-=ddx;
							shootBorneRootGrow(verbose);
							length = nZ;
						}


					}

					/* basal zone */
					if ((dl>0)&&(length<p.lb)) { // length is the current length of the organ
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
					if (((children.size()-additional_childern)<(p.ln.size()+1))&&(length>=p.lb)) 
					{
						for (size_t i=0; (i<p.ln.size()); i++) {
							createLateral(verbose);
							if (getStemRandomParameter()->getLateralType(getNode(nodes.size()-1))==2)
							{
								leafGrow(verbose);
							}
							createSegments(this->dxMin(),dt_,verbose);
							dl-=this->dxMin();
							length+=this->dxMin();}
						createLateral(verbose);
						if (getStemRandomParameter()->getLateralType(getNode(nodes.size()-1))==2){
										leafGrow(verbose);
						}
					}
					if((length>=p.lb)&&((p.ln.size()+1)!=(children.size()))){
						std::stringstream errMsg;
						errMsg <<"Stem::simulate(): different number of realized laterals ("<<children.size()<<
						") and max laterals ("<<p.ln.size()+1<<")";
						throw std::runtime_error(errMsg.str().c_str());
					}
					//internodal elongation, if the basal zone of the stem is created and still has to grow
					double maxInternodeDistance = p.getK()-p.la - p.lb;//maximum length of branching zone
					if((dl>0)&&(length>=p.lb)){
							int nn = children.at(p.ln.size())->parentNI; //node carrying the last lateral == end of branching zone
							double currentInternodeDistance = getLength(nn) - p.lb; //actual length of branching zone
							double ddx = std::min(maxInternodeDistance-currentInternodeDistance, dl);//length to add to branching zone 

							if(ddx > 0){
								internodalGrowth(ddx, dt_,verbose);
								dl -= ddx;
							length += ddx;
								
							}
						}
					/* apical zone */
					//only grows once the basal and branching nodes are developped
					if ((dl>0)&&(length>=(maxInternodeDistance + p.lb))) {
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
 * Simulates internodal growth of dl for this stem
 * divid total stem growth between the phytomeres
 * currently two option:
 * growth devided equally between the phytomeres or 
 * the phytomere grow sequentially
 *
 * @param 	dl			total length of the segments that are created [cm]
 * @param	verbose		print information
 */
void Stem::internodalGrowth(double dl,double dt_, bool verbose)
{
	
	const StemSpecificParameter& p = *param(); // rename
	const StemRandomParameter& Prandom = *getStemRandomParameter();
	std::vector<double> toGrow(p.ln.size());
	double dl_;
	double missing = 0.;
	if(Prandom.nodalGrowth==0){
		toGrow[0] = dl;
		std::fill(toGrow.begin()+1,toGrow.end(),0) ;
	}
	if(Prandom.nodalGrowth ==1){std::fill(toGrow.begin(),toGrow.end(),dl/p.ln.size()) ; }
	size_t i=1;
	while( (dl >0)&&(i<=p.ln.size()) ) {
		//if the phytomere can do a growth superior to the mean phytomere growth, we add the value of "missing" 
		//(i.e., length left to grow to get the predefined total growth of the branching zone)
		int nn1 = children.at(i-1)->parentNI; //node at the beginning of phytomere
		double length1 = getLength(nn1);
		int nn2 = children.at(i)->parentNI; //node at end of phytomere
		double availableForGrowth = p.ln.at(i-1) -( getLength(nn2) - length1 ) ;//difference between maximum and current length of the phytomere
		dl_ = std::min(toGrow[i-1],availableForGrowth);
		if(i< p.ln.size()){
			toGrow[i] +=  toGrow[i-1] - dl_ ;
		}else{missing = toGrow[i-1] - dl_;}
		if(dl_ > 0){createSegments(dl_,dt_,verbose, i ); dl -= dl_;}
		i++;
		
	}
	if(missing > 1e-6){//this sould not happen as computed dl to be <= sum(availableForGrowth)
		std::stringstream errMsg;
		errMsg <<"Stem::internodalGrowth length left to grow: "<<missing;
		throw std::runtime_error(errMsg.str().c_str());
	}
}


/**
 * Creates a new lateral by calling Stem::createNewstem().
 *
 * Overwrite this method to implement more specialized stem classes.
 */
void Stem::createLateral(bool silence)
{
	auto sp = param(); // rename
	int lt = getStemRandomParameter()->getLateralType(getNode(nodes.size()-1));//if lt ==2, don't add lateral as leaf is added instead
	double ageLN = this->calcAge(sp->lb); // age of stem when first lateral node is created
	double delay = sp->delayLat * children.size();	//time the lateral has to wait before growing				  
	Matrix3d h = Matrix3d(); //not needed anymore
	int lnf = getStemRandomParameter()->lnf;
	if (lnf == 2&& lt !=2) {
		auto lateral = std::make_shared<Stem>(plant.lock(), lt, h, delay/2, shared_from_this(),  nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		auto lateral2 = std::make_shared<Stem>(plant.lock(), lt, h, delay/2, shared_from_this(),  nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	} else if (lnf==3&& lt !=2) {
		auto lateral = std::make_shared<Stem>(plant.lock(), lt, h, delay/2, shared_from_this(), nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		auto lateral2 = std::make_shared<Stem>(plant.lock(), lt, h, delay/2, shared_from_this(), nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	} else if (lnf==4 && lt !=2) {
		auto lateral = std::make_shared<Stem>(plant.lock(), lt, h, delay, shared_from_this(),nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	} else if (lnf==5 && lt>0) {
		auto lateral = std::make_shared<Stem>(plant.lock(), lt, h, delay/2,  shared_from_this(), nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

		auto lateral2 = std::make_shared<Stem>(plant.lock(), lt, h, delay/2,  shared_from_this(),  nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	} else if (lt !=2) {
		auto lateral = std::make_shared<Stem>(plant.lock(), lt, h, delay, shared_from_this(),  nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}

}
/**
 * @return Current absolute heading of the organ at node n, based on initial heading, or segment before
 */
Vector3d Stem::getiHeading()  const
{	
	Matrix3d iHeading;
	if (getParent()->organType()==Organism::ot_seed) { // from seed?
		iHeading = Matrix3d(Vector3d(0, 0, 1), Vector3d(0, 1, 0), Vector3d(1, 0, 0));
	}else{
		Vector3d vIHeading = getParent()->heading(parentNI);
		iHeading = Matrix3d::ons(vIHeading);};
	auto heading = iHeading.column(0);
	Vector3d new_heading = Matrix3d::ons(heading).times(this->partialIHeading);
	return Matrix3d::ons(new_heading).column(0);
}

/*
 *
 */
void Stem::leafGrow(bool silence)
{
	auto sp = param(); // rename
	double ageLN = this->calcAge(sp->lb); // age of stem when lateral node is created
	double delay = sp->delayLat * children.size();	//time the lateral has to wait before growing	
	Matrix3d h = Matrix3d(); // current heading in absolute coordinates TODO (revise??)
	int lt =2;
	int lnf = getStemRandomParameter()->lnf;
	if (lnf==2) {
		auto lateral = std::make_shared<Leaf>(plant.lock(), lt,  h, delay/2, shared_from_this(), nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		auto lateral2 = std::make_shared<Leaf>(plant.lock(), lt,  h, delay/2 ,shared_from_this(), nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	} else if (lnf==3) {
		auto lateral = std::make_shared<Leaf>(plant.lock(), lt, h, delay/2,  shared_from_this(), nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		auto lateral2 = std::make_shared<Leaf>(plant.lock(), lt, h,  delay/2, shared_from_this(),nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	} else if (lnf==4) {
		auto lateral = std::make_shared<Leaf>(plant.lock(), lt, h, delay/2,  shared_from_this(),  nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		auto lateral2 = std::make_shared<Leaf>(plant.lock(), lt,  h, delay/2, shared_from_this(), nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass
	} else if (lnf==5) {
		auto lateral = std::make_shared<Leaf>(plant.lock(), lt, h, delay/2, shared_from_this(),  nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		//std::cout <<"leaf heading is "<<h.toString()<< "\n";
		auto lateral2 = std::make_shared<Leaf>(plant.lock(), lt, h, delay/2, shared_from_this(), nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	} else { // TODO error message or warning?
		auto lateral = std::make_shared<Leaf>(plant.lock(), lt,  h, delay, shared_from_this(),  nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
}

/*
 *
 */
void Stem::shootBorneRootGrow(bool verbose)
{
	//const double maxT = 365.; // maximal simulation time
	auto sp = this->param(); // rename
	auto stem_p = this->param();
	auto p = plant.lock();
	auto p_seed = p->getOrganRandomParameter(Organism::ot_seed,0);
	int st = p->getParameterSubType(Organism::ot_root, "shootborne");
	if (st>0) {
		shootborneType = st;
	} // otherwise stick with default
	try {
		p->getOrganRandomParameter(Organism::ot_root, shootborneType); // if the type is not defined an exception is thrown
	} catch (...) {
		if (verbose) {
			std::cout << "Seed::initialize:Shootborne root type #" << shootborneType << " was not defined, using tap root parameters instead\n";
		}
		shootborneType =1;
	}
	//        std::cout <<"subtype ="<<stem_p->subType <<"stem getPhytomerId =" <<getphytomerId(stem_p->subType)<< "\n";
	addPhytomerId(stem_p->subType);

	int nC = getPlant()->getSeed()->param()->nC;
	double nZ = getPlant()->getSeed()->param()->nz;
	double res = nZ-floor(nZ / dx())*dx();
	if(res < dxMin() && res != 0){
		if(res <= dxMin()/2){ nZ -= res;
		}else{nZ =  floor(nZ / dx())*dx() + dxMin();}
		if(verbose){std::cout<<"\nStem::shootBorneRootGrow: nZ changed to "<<nZ<<" for compatibility with dx and dxMin"<<std::endl;}
	}
	if ( nC>0 ) { // only if there are any shootborne roots && (p_seed->firstSB+p_seed->delaySB<maxT)
		std::cout<<"seed nC is "<< nC <<"\n";
		std::cout<<"seed nZ is "<< nZ <<"\n";
		double ageLN = this->calcAge(getLength(true)); // age of stem when lateral node is created
		double ageLG = this->calcAge(getLength(true)+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		for (int i=0; i< nC; i++) {
			double  beta = i*M_PI*getStemRandomParameter()->rotBeta;
			Vector3d newHeading = iHeading.times(Vector3d::rotAB(0,beta));
			auto shootBorneRoot = std::make_shared<Root>(plant.lock() , shootborneType, newHeading, delay,
					shared_from_this(), nodes.size() - 1);
			children.push_back(shootBorneRoot);
			shootBorneRoot->simulate(age-ageLN,verbose);
			std::cout<<"root grow number "<<i<<"\n";
		}
	}

}




/**
 * @return The organTypeParameter from the plant
 */
std::shared_ptr<StemRandomParameter> Stem::getStemRandomParameter() const
{
	return std::static_pointer_cast<StemRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_stem, param_->subType));
}

/**
 * @return Parameters of the specific organ
 */
std::shared_ptr<const StemSpecificParameter> Stem::param() const
{
	return std::static_pointer_cast<const StemSpecificParameter>(param_);
}

/*
 * Quick info about the object for debugging
 * additionally, use getParam()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Stem::toString() const
{
	std::stringstream newstring;
	newstring << "; initial heading: " << iHeading.toString() << ", parent node index" << parentNI << ".";
	return Organ::toString()+newstring.str();
}

} // namespace CPlantBox
