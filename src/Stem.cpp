#include "Stem.h"
#include "Leaf.h"
#include "Root.h"

#include <memory>

namespace CPlantBox {

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
//Stem::Stem(Organism* plant, int type, Vector3d pheading, double delay, Organ* parent, double pbl, int pni) :Organ(plant,parent,Organism::ot_stem,type,delay), pni(pni), pbl(pbl)
//{
//	/*the relative heading is maulfunctioning
//	so it is disabled and rerolled to old heading*/
//
//    initialStemHeading = pheading;
//	//  std::cout << "stem pni = "<< pni<< std::endl;
//	//  std::cout << "Stem constructor \n";
//	StemRandomParameter* sttp = (StemRandomParameter*) plant->getOrganRandomParameter(Organism::ot_stem, type);
//	param_ = sttp->realize(); // throw the dice
//	StemSpecificParameter* stem_p = (StemSpecificParameter*) param_;
//
////        std::cout <<"subtype ="<<stem_p->subType <<"stem getPhytomerId =" <<getphytomerId(stem_p->subType)<< "\n";
//		addPhytomerId(stem_p->subType);
//	double beta = getphytomerId(stem_p->subType)*M_PI*sttp->RotBeta + M_PI*plant->rand()*sttp->BetaDev ;  //+ ; //2 * M_PI*plant->rand(); // initial rotation
//	Matrix3d ons = Matrix3d::ons(initialStemHeading);
//	if (sttp->InitBeta >0 && getphytomerId(stem_p->subType)==0 ){
//		beta = beta + sttp->InitBeta;
//	}
//
//	//ons.times(Matrix3d::rotX(beta));
//
//	double theta = M_PI*stem_p->theta;
//	if (parent->organType() != Organism::ot_seed) { // scale if not a base root
//		double scale = sttp->f_sa->getValue(parent->getNode(pni), this);
//		theta *= scale;
//	}
//	//ons.times(Matrix3d::rotZ(theta));
//	this->initialStemHeading = ons.times(Vector3d::rotAB(theta,beta)); // new initial heading
//    age = -delay; // the root starts growing when age>0
//	alive = 1; // alive per default
//
//	this->parent = parent;
//	parent_base_length = pbl;
//	parent_ni = pni;
//	length = 0;
//
//	// initial node
//	if (parent->organType()!=Organism::ot_seed) { // the first node of the base stems must be created in Seed::initialize()
//		// otherwise, don't use addNode for the first node of the stem,
//		// since this node exists already and does not need a new identifier
//		nodes.push_back(Vector3d());
//		nodeIds.push_back(parent->getNodeId(pni));
//		nodeCTs.push_back(parent->getNodeCT(pni)+delay);
//	}
//
//}

/**
 * Simulates growth of this stem for a time span dt
 *
 * @param dt       time step [day]
 * @param silence  indicates if status messages are written to the console (cout) (default = false)
 */
void Stem::simulate(double dt, bool silence)
{
    moved = false;
    oldNumberOfNodes = nodes.size();

	const StemSpecificParameter* sp = param(); // rename
	const StemRandomParameter* sttp = getStemRandomParameter();

	// increase age
	if (age+dt>sp->rlt) { // stem life time
		dt=sp->rlt-age; // remaining life span
		alive = false; // this stem is dead
	}
	age+=dt;

	if (alive) { // dead stem wont grow
		// probabilistic branching model (todo test)
		if ((age>0) && (age-dt<=0)) { // the stem emerges in this time step
			double stem_P = sttp->f_sbp->getValue(nodes.back(),this);
			if (stem_P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
				double stem_p = 1.-std::pow((1.-stem_P), dt); //probability of emergence in this time step
				//        std::cout <<stem_P<<", "<<stem_p<< "\n";
				if (plant->rand()>stem_p) { // not rand()<p
					age -= dt; // the stem does not emerge in this time step
				}
			}
		}

		if (age>0) {

			// children first (lateral stems grow even if base stem is inactive)
			for (auto c:children) {
				c->simulate(dt,silence);
			}

			if (active) {

				//      std::cout << "type\t" << sp->subType << "\n"  << "lb\t"<< sp->lb <<"\n" << "la\t"<< sp->la <<"\n" << "ln\t";
				//	for (size_t i=0; i<sp->ln.size(); i++) {
				//		std::cout << sp->ln[i] << "\t";
				//	}

				// length increment
				double length_ = StemGetLength(std::max(age-dt,0.)); // length of the stem for unimpeded growth (i.e. length_==length for unimpeded growth)
				double targetlength = StemGetLength(age);
				double e = targetlength-length_; //elongation in time step dt
				double scale = sttp->f_se->getValue(nodes.back(),this); // hope some of this is optimized out if not set
				double dl = std::max(scale*e, double(0)); // length increment, dt is not used anymore

				// create geometry
				if (sp->ln.size()>0) { // stem has laterals
					// basal zone
					if ((dl>0)&&(length<sp->lb)) { // length is the current length of the stem
						if (length+dl<=sp->lb) {
							createSegments(dl,silence);
							length+=dl;
							dl=0;
						} else {
							double ddx = sp->lb-length;
							createSegments(ddx,silence);
							dl-=ddx; // ddx already has been created
							length=sp->lb;
						}
					}
					// branching zone Condition will be changed later
					if ((dl>0)&&(length>=sp->lb)) {
						double s = sp->lb; // summed length
						for (size_t i=0; ((i<sp->ln.size()) && (dl>0)); i++) {
							s+=sp->ln.at(i);
							if (length<s) {
								if (i==children.size()) { // new internode leaf and shootBorneRoot
									//if (sp->subType==3)
									{ //this decide which successor grow the leaf (TODO) adding it to parameterfile
										leafGrow(silence, nodes.back());
										//ShootBorneRootGrow(silence);
									//} else {
										createLateral(silence);
										//LeafGrow(silence);
									}
									//                        ShootBorneRootGrow(silence);
								}
								if (length+dl<=s) { // finish within inter-lateral distance i
									createSegments(dl,silence);
									length+=dl;
									dl=0;
								} else { // grow over inter-lateral distance i
									double ddx = s-length;
									createSegments(ddx,silence);
									dl-=ddx;
									length=s;
								}
							}
						}
						if (dl>0) {
							if (sp->ln.size()== children.size()) { // new lateral (the last one)
									createLateral(silence);
							}
						}
					}
					// apical zone
					if (dl>0) {
						createSegments(dl,silence);
						length+=dl;
					}
					else {
						if (sp->subType == 3)
						{
							leafGrow(silence, nodes.back());
						}
					}
				} else { // no laterals
					if (dl>0) {
						createSegments(dl,silence);
						length+=dl;
					}
				} // if laterals
			} // if active
			active = StemGetLength(std::max(age,0.))<(sp->getK()-dx()/10); // become inactive, if final length is nearly reached
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
    if (name=="nob") { return param()->nob; } // number of branches
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
double Stem::getCreationTime(double length)
{
	assert(length>=0);
	double age = StemGetAge(length);
	assert(age>=0);
	if (parent->organType()!=Organism::ot_seed) {
		if (parent->organType()==Organism::ot_stem) {
			double pl = pbl+((Stem*)parent)->getStemRandomParameter()->la; // parent length, when this stem was created
			double pAge=((Stem*)parent)->getCreationTime(pl);
			return age+pAge;
		} else { // organ type is seed
			return age;
		}
	} else {
		return age;
	}
}

/**
 * Analytical length of the stem at a given age
 *
 * @param age          age of the stem [day]
 */
double Stem::StemGetLength(double age)
{
	assert(age>=0);
	return getStemRandomParameter()->f_gf->getLength(age,getStemRandomParameter()->r,getStemRandomParameter()->getK(),this);
}

/**
 * Analytical age of the stem at a given length
 *
 * @param length   length of the stem [cm]
 */
double Stem::StemGetAge(double length)
{
	assert(length>=0);
	return getStemRandomParameter()->f_gf->getAge(length,getStemRandomParameter()->r,getStemRandomParameter()->getK(),this);
}

/**
 *
 */
double Stem::dx() const
{
	return getStemRandomParameter()->dx;
}

/**
 * Creates a new lateral by calling Stem::createNewstem().
 *
 * Overwrite this method to implement more specialized stem classes.
 */
void Stem::createLateral(bool silence)
{
	const StemSpecificParameter* sp = param(); // rename
	int lt = getStemRandomParameter()->getLateralType(getNode(nodes.size()-1));

	if (getStemRandomParameter()->lnf==2&& lt>0) {
		double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading

		Stem* lateral = new Stem(plant, lt, h, delay, this, nodes.size() - 1, length);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		Stem* lateral2 = new Stem(plant, lt, h, delay, this, nodes.size() - 1, length);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
	else if (getStemRandomParameter()->lnf==3&& lt>0) {
		double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Stem* lateral = new Stem(plant, lt, h, delay, this, nodes.size() - 1, length);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		Stem* lateral2 = new Stem(plant, lt, h, delay, this, nodes.size() - 1, length);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
	else if (getStemRandomParameter()->lnf==4 && lt>0) {
		double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Stem* lateral = new Stem(plant, lt, h, delay, this, nodes.size() - 1, length);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

    }
	else if (getStemRandomParameter()->lnf==5 && lt>0) {
		double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Stem* lateral = new Stem(plant, lt, h, delay,  this, nodes.size() - 1, length);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

		Stem* lateral2 = new Stem(plant, lt, h, delay,  this, nodes.size() - 1, length);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

	} else if (lt>0) {
		double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Stem* lateral = new Stem(plant, lt, h, delay, this, nodes.size() - 1, length);
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
	const StemSpecificParameter* sp = param(); // rename
	//int lt = getStemRandomParameter()->getLateralType(getNode(nodes.size()-1));
	//  	std::cout << "LeafGrow createLateral()\n";
	//  	std::cout << "LeafGrow type " << lt << "\n";

	//if (lt>0) {
	//int lt = getStemRandomParameter()->getLateralType(getNode(nodes.size()-1));
	int lt =2;
	if (getStemRandomParameter()->lnf==2&& lt>0) {
		double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Leaf* lateral = new Leaf(plant, lt,  h, delay, this, nodes.size() - 1, length);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		Leaf* lateral2 = new Leaf(plant, lt,  h, delay, this, nodes.size() - 1, length);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
	else if (getStemRandomParameter()->lnf==3&& lt>0) {
		double ageLN = this->StemGetAge(length); // athis, ge of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Leaf* lateral = new Leaf(plant, lt, h, delay,  this, nodes.size() - 1, length);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		Leaf* lateral2 = new Leaf(plant, lt, h,  delay, this, nodes.size() - 1, length);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
	else if (getStemRandomParameter()->lnf==4 && lt>0) {
		double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Leaf* lateral = new Leaf(plant, lt, h, delay,  this, nodes.size() - 1, length);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		Leaf* lateral2 = new Leaf(plant, lt,  h, delay, this, nodes.size() - 1, length);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass
	}else if (getStemRandomParameter()->lnf==5&& lt>0) {
				double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading

		Leaf* lateral = new Leaf(plant, lt, h, delay, this, nodes.size() - 1, length);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);

		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

        //std::cout <<"leaf heading is "<<h.toString()<< "\n";

		Leaf* lateral2 = new Leaf(plant, lt, h, delay, this, nodes.size() - 1, length);
		//lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);

		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

	}else if (lt>0) {
		double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Leaf* lateral = new Leaf(plant, lt,  h, delay,this, nodes.size() - 1, length);
		//lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		}
	//}


}


/*
 *
 */
void Stem::shootBorneRootGrow(bool silence)
{
	StemSpecificParameter* sp = param(); // rename
	int lt = getStemRandomParameter()->getLateralType(getNode(nodes.size()-1));
	//    std::cout << "ShootBorneRootGrow createLateral()\n";
	//    std::cout << "ShootBorneRootGrow lateral type " << lt << "\n";
	if ( lt > 0 ) {
		double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		int nodeToGrowShotBorneRoot = 2;
		Vector3d sbrheading(0,0,-1); //just a test heading
		Root* shootBorneRootGrow = new Root(plant , 5, sbrheading, delay ,this, length, nodeToGrowShotBorneRoot);
		if (nodes.size() > nodeToGrowShotBorneRoot ) {
			//                                ShootBorneRootGrow->addNode(getNode(NodeToGrowShotBorneRoot), length);
			children.push_back(shootBorneRootGrow);
			shootBorneRootGrow->simulate(age-ageLN,silence);// pass time overhead (age we want to achieve minus current age)
		}
	}
}


/**
 *  Creates nodes and node emergence times for length l,
 *  and updates the stem heading
 *
 *  Cecks that each new segments length is <= dx but >= ddx
 *
 *  @param l       length the stem growth [cm]
 */
void Stem::createSegments(double l, bool silence)
{
	//   std::cout << "create Stem Segments("<< l << ")\n";
	assert(l>0);
	double sl=0; // summed length of created segment

	// shift first node to axial resolution
	int nn = nodes.size();
	if (oldNumberOfNodes==0) { // first call of createSegments (in stem::simulate)
		if (nn>1) {
			auto n2 = nodes.at(nn-2);
			auto n1 = nodes.at(nn-1);
			double olddx = n1.minus(n2).length();
			if (olddx<dx()*0.99) { // shift node instead of creating a new node

				double newdx = std::min(dx()-olddx, l);
                double sdx = olddx + newdx; // length of new segment


                Vector3d newdxv = getIncrement(n2, sdx);
				nodes[nn - 1] = Vector3d(n2.plus(newdxv));

				double ct = this->getCreationTime(length + newdx);
				nodeCTs[nn - 1] = std::max(ct, plant->getSimTime()); // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
				oldNumberOfNodes = nn;
				l -= newdx;

				if (l<=0) { // ==0 should be enough

					//if (tipleaf) {
						leafGrow(silence, nodes.back()); //the stem tip will grow leafs
						std::cout << "leaf grow\n";
					//	tipleaf = false;

					//}
									return;
				}
			}
		}


		oldNumberOfNodes = nn; // CHECK
	}

	if (l<smallDx) {
		if (!silence) {
			//      std::cout << "skipped small segment l<Dx (<"<< smallDx << ") \n";
		}
		return;
	}

	int n = floor(l/dx());
	// create n+1 new nodes
	for (int i=0; i<n+1; i++) {

		double sdx; // segment length (<=dx)
		if (i<n) {  // normal case
			sdx = dx();
		} else { // last segment
			sdx = l-n*dx();
			if (sdx<smallDx) {
				if (!silence) {
					//          std::cout << "skipped small segment i<n (<"<< smallDx << ") \n";
				}
				return;
			}
		}
		sl+=sdx;
		Vector3d newdx = getIncrement(nodes.back(), sdx);

			Vector3d newnode = Vector3d(nodes.back().plus(newdx));
		double ct = this->getCreationTime(length+sl);
		ct = std::max(ct,plant->getSimTime()); // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
		//     std::cout<<"add node "<<newnode.toString()<<"\n";
		//if (tipleaf) {
			addNode(newnode, ct);
		//}

	} // for

}

Vector3d Stem::getIncrement(const Vector3d& p, double sdx) {
    Vector3d h = heading(); // current heading
    Matrix3d ons = Matrix3d::ons(h);
    Vector2d ab = getStemRandomParameter()->f_tf->getHeading(p, ons, sdx, this);
    Vector3d sv = ons.times(Vector3d::rotAB(ab.x,ab.y));
    return sv.times(sdx);
}

/**
 * @return The RootTypeParameter from the plant
 */
StemRandomParameter* Stem::getStemRandomParameter() const
{
    return (StemRandomParameter*)plant->getOrganRandomParameter(Organism::ot_root, param_->subType);
}

/**
 * @return Parameters of the specific root
 */
StemSpecificParameter* Stem::param() const
{
    return (StemSpecificParameter*)param_;
}


Vector3d Stem::heading() const {
	Vector3d h;
	if (nodes.size()>1) {
		h = nodes.back().minus(nodes.at(nodes.size() - 2)); // getHeading(b-a)
        h.normalize();
	} else {
		h = initialStemHeading;
	}
	return h;
}

/*
 * Quick info about the object for debugging
 */
std::string Stem::toString() const
{
	std::stringstream str;
	str << "Stem #"<< id <<": type "<<param()->subType << ", length: "<< length << ", age: " <<age<<" with "<< children.size() << " laterals\n";
	return str.str();
}


} // namespace CPlantBox
