#include "Stem.h"

#include "Leaf.h"
#include "Root.h"
#include "Plant.h"

namespace CPlantBox {

std::vector<int> Stem::phytomerId = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/**
 * Constructs a root from given data.
 * The organ tree must be created, @see Organ::setPlant, Organ::setParent, Organ::addChild
 * Organ geometry must be created, @see Organ::addNode, ensure that this->getNodeId(0) == parent->getNodeId(pni)
 *
 * @param id        the organ's unique id (@see Organ::getId)
 * @param param     the organs parameters set, ownership transfers to the organ
 * @param alive     indicates if the organ is alive (@see Organ::isAlive)
 * @param active    indicates if the organ is active (@see Organ::isActive)
 * @param age       the current age of the organ (@see Organ::getAge)
 * @param length    the current length of the organ (@see Organ::getLength)
 * @param iheading  the initial heading of this root
 * @param pbl       base length of the parent root, where this root emerges
 * @param pni       local node index, where this root emerges
 * @param moved     indicates if nodes were moved in the previous time step (default = false)
 * @param oldNON    the number of nodes of the previous time step (default = 0)
 */
Stem::Stem(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
		Vector3d iheading, double pbl, int pni, bool moved, int oldNON)
:Organ(id, param, alive, active, age, length, iheading, pbl, pni, moved,  oldNON)
{ }

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
Stem::Stem(std::shared_ptr<Organism> plant, int type, Vector3d iheading, double delay,  std::shared_ptr<Organ> parent, double pbl, int pni)
:Organ(plant, parent, Organism::ot_stem, type, delay, iheading, pbl, pni)
{
	/*the relative heading is maulfunctioning
	so it is disabled and rerolled to old heading*/
	assert(parent!=nullptr && "Stem::Stem parent must be set");

	auto stem_p = this->param();

	//        std::cout <<"subtype ="<<stem_p->subType <<"stem getPhytomerId =" <<getphytomerId(stem_p->subType)<< "\n";
	addPhytomerId(stem_p->subType);
	double beta = getphytomerId(stem_p->subType)*M_PI*getStemRandomParameter()->rotBeta +
			M_PI*plant->rand()*getStemRandomParameter()->betaDev ;  //+ ; //2 * M_PI*plant->rand(); // initial rotation
	Matrix3d ons = Matrix3d::ons(iHeading);
		beta = beta + getStemRandomParameter()->initBeta;

	if (getStemRandomParameter()->initBeta >0 && getphytomerId(stem_p->subType)==0 ){
		beta = beta + getStemRandomParameter()->initBeta;
	}
	//ons.times(Matrix3d::rotX(beta));
	double theta = M_PI*stem_p->theta;
	if (parent->organType()!=Organism::ot_seed) { // scale if not a base root
		double scale = getStemRandomParameter()->f_sa->getValue(parent->getNode(pni), parent);
		theta *= scale;
	}
	//ons.times(Matrix3d::rotZ(theta));
	this->iHeading = ons.times(Vector3d::rotAB(theta,beta)); // new initial heading
	if (parent->organType()!=Organism::ot_seed) { // initial node
		// assert(pni+1 == parent->getNumberOfNodes() && "at object creation always at last node");
		addNode(parent->getNode(pni), parent->getNodeId(pni), parent->getNodeCT(pni)+delay);
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
	moved = false;
	oldNumberOfNodes = nodes.size();
	auto p_all = plant.lock();
	auto p_stem = p_all->getOrganRandomParameter(Organism::ot_stem);

	int nC = getPlant()->getSeed()->param()->nC; //number of the shoot born root
	double nZ = getPlant()->getSeed()->param()->nz; // distance between shoot born root and the seed

	int additional_childern;
	if (p.subType == 1)
	{
		additional_childern= (int)round(nC); //if it is the main stem, the children should include the shoot borne root
	} else {
		additional_childern = 0;
	}

	if (alive) { // dead roots wont grow

		// increase age
		if (age+dt>p.rlt) { // root life time
			dt=p.rlt-age; // remaining life span
			alive = false; // this root is dead
		}
		age+=dt;

		// probabilistic branching model (todo test)
		if ((age>0) && (age-dt<=0)) { // the root emerges in this time step
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

				// length increment
				double age_ = calcAge(length); // root age as if grown unimpeded (lower than real age)
				double dt_; // time step
				if (age<dt) { // the root emerged in this time step, adjust time step
					dt_= age;
				} else {
					dt_=dt;
				}

				double targetlength = calcLength(age_+dt_);
				double e = targetlength-length; // unimpeded elongation in time step dt
				double scale = getStemRandomParameter()->f_sa->getValue(nodes.back(),shared_from_this());
				double dl = std::max(scale*e, 0.); // length increment

				// create geometry
				if (p.ln.size()>0) { // stem has laterals
					//std::cout<<"sim seed nC is"<< nC<<"\n";
					//std::cout<<"sim seed nZ is"<< nZ<<"\n";
					/*
                    shoot born root
					 */
					if ((dl>0)&&(length< nZ)) {
						if (length+dl <= nZ) {
							createSegments(dl,verbose);
							length+=dl;
							dl=0;
						} else {
							double ddx = nZ - length;
							createSegments(dl,verbose);

							dl-=ddx; // ddx already has been created
							shootBorneRootGrow(verbose);
							length = nZ;
						}


					}

					/* basal zone */
					if ((dl>0)&&(length<p.lb)) { // length is the current length of the root
						if (length+dl<=p.lb) {
							createSegments(dl,verbose);
							length+=dl;
							dl=0;
						} else {
							double ddx = p.lb-length;
							createSegments(ddx,verbose);
							dl-=ddx; // ddx already has been created
							length=p.lb;
						}
					}
					/* branching zone */
					if ((dl>0)&&(length>=p.lb)) {
						double s = p.lb; // summed length
						for (size_t i=0; ((i<p.ln.size()) && (dl>0)); i++) {
							s+=p.ln.at(i);
							if (length<s) {
								if (i==children.size()-additional_childern) { // new lateral
									//this decide which successor grow the leaf (TODO) adding it to parameterfile

									createLateral(verbose);
									if (getStemRandomParameter()->getLateralType(getNode(nodes.size()-1))==2){
										leafGrow(verbose, nodes.back());
									}else
									{
										// terminates the loop if the stem is not organtype 2
										//break;
									}

								}
								if (length+dl<=s) { // finish within inter-lateral distance i
									createSegments(dl,verbose);
									length+=dl;
									dl=0;
								} else { // grow over inter-lateral distance i
									double ddx = s-length;
									createSegments(ddx,verbose);
									dl-=ddx;
									length=s;
								}
							}
						}
						if (p.ln.size()==children.size()-additional_childern) { // new lateral (the last one)
							createLateral(verbose);
						}
					}
					/* apical zone */
					if (dl>0) {
						createSegments(dl,verbose);
						length+=dl;
					} else {
						//if (p.subType == 3) {

						//}
					}
				} else { // no laterals
					if (dl>0) {
						createSegments(dl,verbose);
						length+=dl;
						leafGrow(verbose, nodes.back());
					}
				} // if lateralgetLengths
			} // if active
			active = length<(p.getK()-dx()/10); // become inactive, if final length is nearly reached
		}
	} // if alive
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
	if (name=="la") { return param()->la; } // apical zone [cm]
	if (name=="nob") { return param()->nob(); } // number of branches
	if (name=="r"){ return param()->r; }  // initial growth rate [cm day-1]
	if (name=="radius") { return param()->a; } // root radius [cm]
	if (name=="a") { return param()->a; } // root radius [cm]
	if (name=="theta") { return param()->theta; } // angle between root and parent root [rad]
	if (name=="rlt") { return param()->rlt; } // root life time [day]
	if (name=="k") { return param()->getK(); }; // maximal root length [cm]
	if (name=="lnMean") { // mean lateral distance [cm]
		auto& v =param()->ln;
		return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
	}
	if (name=="lnDev") { // standard deviation of lateral distance [cm]
		auto& v =param()->ln;
		double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
		double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
		return std::sqrt(sq_sum / v.size() - mean * mean);
	}
	if (name=="volume") { return param()->a*param()->a*M_PI*getLength(); } // // root volume [cm^3]
	if (name=="surface") { return 2*param()->a*M_PI*getLength(); }
	if (name=="type") { return this->param_->subType; }  // in CPlantBox the subType is often called just type
	if (name=="iHeadingX") { return iHeading.x; } // root initial heading x - coordinate [cm]
	if (name=="iHeadingY") { return iHeading.y; } // root initial heading y - coordinate [cm]
	if (name=="iHeadingZ") { return iHeading.z; } // root initial heading z - coordinate [cm]
	if (name=="parentBaseLength") { return parentBaseLength; } // length of parent root where the lateral emerges [cm]
	if (name=="parentNI") { return parentNI; } // local parent node index where the lateral emerges
	return Organ::getParameter(name);
}


/**
 * Analytical creation (=emergence) time of a node at a length along the stem
 *
 * @param length   length of the stem [cm]
 */
double Stem::calcCreationTime(double length)
{
	assert(length >= 0 && "Stem::getCreationTime() negative length");
	double stemage = calcAge(length);
	stemage = std::min(stemage, age);
	assert(stemage >= 0 && "Stem::getCreationTime() negative stem age");
	return stemage+nodeCTs[0];
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
 *
 * @param length   length of the stem [cm]
 */
double Stem::calcAge(double length)
{
	assert(length>=0 && "Stem::calcAge() negative root age");
	return getStemRandomParameter()->f_gf->getAge(length,getStemRandomParameter()->r,param()->getK(),shared_from_this());
}

/**
 * Creates a new lateral by calling Stem::createNewstem().
 *
 * Overwrite this method to implement more specialized stem classes.
 */
void Stem::createLateral(bool silence)
{
	//std::cout<<"lnf is = "<<getStemRandomParameter()->lnf<<"\n";
	auto sp = param(); // rename
	int lt = getStemRandomParameter()->getLateralType(getNode(nodes.size()-1));
        double ageLN = this->calcAge(length); // age of stem when lateral node is created
		double ageLG = this->calcAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
	if (getStemRandomParameter()->lnf == 2&& lt>0) {


		auto lateral = std::make_shared<Stem>(plant.lock(), lt, h, delay, shared_from_this(), length, nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		auto lateral2 = std::make_shared<Stem>(plant.lock(), lt, h, delay, shared_from_this(), length, nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
	else if (getStemRandomParameter()->lnf==3&& lt>0) {

		auto lateral = std::make_shared<Stem>(plant.lock(), lt, h, delay, shared_from_this(), length, nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		auto lateral2 = std::make_shared<Stem>(plant.lock(), lt, h, delay, shared_from_this(), length, nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
	else if (getStemRandomParameter()->lnf==4 && lt>0) {

		auto lateral = std::make_shared<Stem>(plant.lock(), lt, h, delay, shared_from_this(), length, nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

	}
	else if (getStemRandomParameter()->lnf==5 && lt>0) {

		auto lateral = std::make_shared<Stem>(plant.lock(), lt, h, delay,  shared_from_this(), length, nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

		auto lateral2 = std::make_shared<Stem>(plant.lock(), lt, h, delay,  shared_from_this(), length, nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

	} else if (lt>0) {

		auto lateral = std::make_shared<Stem>(plant.lock(), lt, h, delay, shared_from_this(), length, nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}

}


/*
 *
 */
void Stem::leafGrow(bool silence, Vector3d bud)
{
	auto sp = param(); // rename
	//int lt = getStemRandomParameter()->getLateralType(getNode(nodes.size()-1));
	//  	std::cout << "LeafGrow createLateral()\n";
	//  	std::cout << "LeafGrow type " << lt << "\n";
	double ageLN = this->calcAge(length); // age of stem when lateral node is created
		double ageLG = this->calcAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
	int lt =2;
	if (getStemRandomParameter()->lnf==2) {

		auto lateral = std::make_shared<Leaf>(plant.lock(), lt,  h, delay, shared_from_this(), length, nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		auto lateral2 = std::make_shared<Leaf>(plant.lock(), lt,  h, delay, shared_from_this(), length, nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
	else if (getStemRandomParameter()->lnf==3) {

		auto lateral = std::make_shared<Leaf>(plant.lock(), lt, h, delay,  shared_from_this(), length, nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		auto lateral2 = std::make_shared<Leaf>(plant.lock(), lt, h,  delay, shared_from_this(), length, nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
	else if (getStemRandomParameter()->lnf==4 ) {

		auto lateral = std::make_shared<Leaf>(plant.lock(), lt, h, delay,  shared_from_this(), length, nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		auto lateral2 = std::make_shared<Leaf>(plant.lock(), lt,  h, delay, shared_from_this(), length, nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass
	}else if (getStemRandomParameter()->lnf==5) {


		auto lateral = std::make_shared<Leaf>(plant.lock(), lt, h, delay, shared_from_this(), length, nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);

		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

		//std::cout <<"leaf heading is "<<h.toString()<< "\n";

		auto lateral2 = std::make_shared<Leaf>(plant.lock(), lt, h, delay, shared_from_this(), length, nodes.size() - 1);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);

		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

	}else{
		auto lateral = std::make_shared<Leaf>(plant.lock(), lt,  h, delay, shared_from_this(), length, nodes.size() - 1);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
	//}
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
	if ( nC>0 ) { // only if there are any shootborne roots && (p_seed->firstSB+p_seed->delaySB<maxT)
		std::cout<<"seed nC is"<< nC <<"\n";
		std::cout<<"seed nZ is"<< nZ <<"\n";
		double ageLN = this->calcAge(length); // age of stem when lateral node is created
		double ageLG = this->calcAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		for (int i=0; i< nC; i++) {
			double  beta = i*M_PI*getStemRandomParameter()->rotBeta;
			Matrix3d ons = Matrix3d::ons(iHeading);
			auto shootBorneRoot = std::make_shared<Root>(plant.lock() , shootborneType, ons.times(Vector3d::rotAB(0,beta)), delay ,shared_from_this(), length, nodes.size() - 1);
			children.push_back(shootBorneRoot);
			shootBorneRoot->simulate(age-ageLN,verbose);
			std::cout<<"root grow number"<<i<<"\n";
		}
	}

	//    auto sp = param(); // rename
	//    int lt = getStemRandomParameter()->getLateralType(getNode(nodes.size()-1));
	//    //    std::cout << "ShootBorneRootGrow createLateral()\n";
	//    //    std::cout << "ShootBorneRootGrow lateral type " << lt << "\n";
	//
	//    if ( lt > 0 ) {
	//        double ageLN = this->calcAge(length); // age of stem when lateral node is created
	//        double ageLG = this->calcAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
	//        double delay = ageLG-ageLN; // time the lateral has to wait
	//        int nodeToGrowShotBorneRoot = 2;
	//        Vector3d sbrheading(0,0,-1); //just a test heading
	//        auto shootBorneRootGrow = std::make_shared<Root>(plant.lock() , 5, sbrheading, delay ,shared_from_this(), length, nodeToGrowShotBorneRoot);
	//        if (nodes.size() > nodeToGrowShotBorneRoot ) {
	//            //                                ShootBorneRootGrow->addNode(getNode(NodeToGrowShotBorneRoot), length);
	//            children.push_back(shootBorneRootGrow);
	//            shootBorneRootGrow->simulate(age-ageLN,silence);// pass time overhead (age we want to achieve minus current age)
	//        }
	//    }
}


/**
 * @return Current root heading
 */
Vector3d Stem::heading() const
{
	if (nodes.size()>1) {
		auto h = nodes.back().minus(nodes.at(nodes.size()-2)); // a->b = b-a
		h.normalize();
		return h;
	} else {
		return iHeading;
	}
}

/**
 * Returns the increment of the next segments
 *
 *  @param p       coordinates of previous node
 *  @param sdx     length of next segment [cm]
 *  @return        the vector representing the increment
 */
Vector3d Stem::getIncrement(const Vector3d& p, double sdx)
{
	Vector3d h = heading();
	Matrix3d ons = Matrix3d::ons(h);
	Vector2d ab = getStemRandomParameter()->f_tf->getHeading(p, ons, sdx, shared_from_this());
	Vector3d sv = ons.times(Vector3d::rotAB(ab.x,ab.y));
	return sv.times(sdx);
}

/**
 *  Creates nodes and node emergence times for a length l
 *
 *  Checks that each new segments length is <= dx but >= smallDx
 *
 *  @param l        total length of the segments that are created [cm]
 *  @param verbose  turns console output on or off
 */
void Stem::createSegments(double l, bool verbose)
{
	if (l==0) {
		std::cout << "Stem::createSegments: zero length encountered \n";
		return;
	}
	if (l<0) {
		std::cout << "Stem::createSegments: negative length encountered \n";
	}

	// shift first node to axial resolution
	double shiftl = 0; // length produced by shift
	int nn = nodes.size();
	if (firstCall) { // first call of createSegments (in Root::simulate)
		firstCall = false;

		int pni = -1;
		if (!children.empty()) {
			auto o = children.back();
			pni = o->parentNI;
		}

		if ((nn>1) && (children.empty() || (nn-1 != pni) )) { // don't move a child base node
			Vector3d n2 = nodes[nn-2];
			Vector3d n1 = nodes[nn-1];
			double olddx = n1.minus(n2).length(); // length of last segment
			if (olddx<dx()*0.99) { // shift node instead of creating a new node
				shiftl = std::min(dx()-olddx, l);
				double sdx = olddx + shiftl; // length of new segment
				Vector3d newdxv = getIncrement(n2, sdx);
				nodes[nn-1] = Vector3d(n2.plus(newdxv));
				double et = this->calcCreationTime(length+shiftl);
				nodeCTs[nn-1] = et; // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
				moved = true;
				l -= shiftl;
				if (l<=0) { // ==0 should be enough
					return;
				}
			} else {
				moved = false;
			}
		} else {
			moved = false;
		}
	}
	// create n+1 new nodes
	double sl = 0; // summed length of created segment
	int n = floor(l/dx());
	for (int i = 0; i < n + 1; i++) {

		double sdx; // segment length (<=dx)
		if (i<n) {  // normal case
			sdx = dx();
		} else { // last segment
			sdx = l-n*dx();
			if (sdx<smallDx) { // quit if l is too small
				if (verbose) {
					std::cout << "skipped small segment ("<< sdx <<" < "<< smallDx << ") \n";
				}
				return;
			}
		}
		sl += sdx;
		Vector3d newdx = getIncrement(nodes.back(), sdx);
		Vector3d newnode = Vector3d(nodes.back().plus(newdx));
		double et = this->calcCreationTime(length+shiftl+sl);
		// in case of impeded growth the node emergence time is not exact anymore,
		// but might break down to temporal resolution
		addNode(newnode, et);
	}
}

/**
 *
 */
std::shared_ptr<Plant> Stem::getPlant() {
	return std::static_pointer_cast<Plant>(plant.lock());
}


/**
 *
 */
double Stem::dx() const
{
	return getStemRandomParameter()->dx;
}

/**
 * @return The RootTypeParameter from the plant
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
 * additionally, use getParam()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Stem::toString() const
{
	std::string str = Organ::toString();
	str.replace(0, 5, "Stem");
	std::stringstream newstring;
	newstring << "; initial heading: " << iHeading.toString() << ", parent base length " << parentBaseLength << ", parent node index" << parentNI << ".";
	return str+newstring.str();
}




} // namespace CPlantBox
