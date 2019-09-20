#include "Leaf.h"

namespace CPlantBox {


/**
 * Constructor
 * This is a Copy Paste of the Root.cpp but it works independently, it has its own parameter file (in .lParam file) tropism, f_gf function, txt and vtp writing syleaf.
 * All of those can be modified to fit the real f_gf of the Plant.
 *
 * Typically called by the Plant::Plant(), or Leaf::createNewLeaf().
 * For leaf the initial node and node emergence time (netime) must be set from outside
 *
 * @param plant 		points to the plant
 * @param parent 		points to the parent organ
 * @param subtype		sub type of the leaf
 * @param delay 		delay after which the organ starts to develop (days)
 * @param rheading		relative heading (within parent organ)
 * @param pni			parent node index
 * @param pbl			parent base length
 */
Leaf::Leaf(Organism* plant, int type, Vector3d iheading, double delay, Organ* parent, double pbl, int pni) :Organ(plant,parent,Organism::ot_leaf, type,delay), pni(pni), pbl(pbl)
{
	  initialLeafHeading=iheading;
	//  std::cout << "Leaf pni = "<< pni<< std::endl;
	//  std::cout << "Organism* plant ="<< plant <<" "<< parent<<std::endl;
		param_ = getLeafRandomParameter()->realize(); // throw the dice, specific

	//  std::cout <<", "<<(LeafParameter*) param<< "\n";


        std::cout <<"subtype ="<<param()->subType <<"getleafphytomerID =" <<getleafphytomerID(param()->subType)<< "\n";
		addleafphytomerID(param()->subType);
	double beta = getleafphytomerID(param()->subType)*M_PI*getLeafRandomParameter()->RotBeta + M_PI*plant->rand()*getLeafRandomParameter()->BetaDev ;  //+ ; //2 * M_PI*plant->rand(); // initial rotation
	Matrix3d ons = Matrix3d::ons(initialLeafHeading);
//	if (getLeafRandomParameter()->InitBeta >0 && getleafphytomerID(param()->subType)==0 ){
		beta = beta + getLeafRandomParameter()->InitBeta;
//	}

		if (getLeafRandomParameter()->InitBeta >0 && getLeafRandomParameter()->subType==2 && getLeafRandomParameter()->lnf==5 && getleafphytomerID(2)%4==2 )
		{beta = beta + getLeafRandomParameter()->InitBeta*M_PI;}
		else if (getLeafRandomParameter()->InitBeta >0 && getLeafRandomParameter()->subType==2 && getLeafRandomParameter()->lnf==5 && getleafphytomerID(2)%4==3 )
        {beta = beta + getLeafRandomParameter()->InitBeta*M_PI + M_PI;}
	//ons.times(Matrix3d::rotX(beta));

	double theta = M_PI*param()->theta;
	if (parent->organType() != Organism::ot_seed) { // scale if not a base root
		double scale = getLeafRandomParameter()->f_sa->getValue(parent->getNode(pni), this);
		theta *= scale;
	}
	//ons.times(Matrix3d::rotZ(theta));
	this->initialLeafHeading = ons.times(Vector3d::rotAB(theta,beta)); // new initial heading
    age = -delay; // the root starts growing when age>0
	alive = 1; // alive per default

	this->parent = parent;
	parent_base_length = pbl;
	parent_ni = pni;
	length = 0;
	// initial node
	 if (parent->organType()!=Organism::ot_seed ) { // the first node of the base leafs must be created in Seed::initialize()
	// otherwise, don't use addNode for the first node of the leaf,
	// since this node exists already and does not need a new identifier
	nodes.push_back(Vector3d());
	nodeIds.push_back(parent->getNodeId(pni));
	nodeCTs.push_back(parent->getNodeCT(pni)+delay);
	  }
}

/**
 * Deep copies the organ into the new plant @param rs.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
Organ* Leaf::copy(Organism* plant)
{
    Leaf* l = new Leaf(*this); // shallow copy
    l->parent = nullptr;
    l->plant = plant;
    l->param_ = new LeafSpecificParameter(*param()); // copy parameters
    for (size_t i=0; i< children.size(); i++) {
        l->children[i] = children[i]->copy(plant); // copy laterals
        l->children[i]->setParent(this);
    }
    return l;
}

/**
 * Simulates f_gf of this leaf for a time span dt
 *
 * @param dt       time step [day]
 * @param silence  indicates if status messages are written to the console (cout) (default = false)
 */
void Leaf::simulate(double dt, bool silence)
{
	old_non = 0; // is set in Leaf:createSegments, the zero indicates the first call to createSegments

	// increase age
	if (age+dt>param()->rlt) { // leaf life time
		dt=param()->rlt-age; // remaining life span
		alive = false; // this leaf is dead
	}
	age+=dt;

	if (alive) { // dead leaf wont grow
		// probabilistic branching model (todo test)
		if ((age>0) && (age-dt<=0)) { // the leaf emerges in this time step
			double f_soil = getLeafRandomParameter()->f_sbp->getValue(nodes.back(),this);
			if (f_soil<1.) { // P==1 means the lateral emerges with probability 1 (default case)
				f_soil = 1.-std::pow((1.-f_soil), dt); //probability of emergence in this time step
				//        std::cout <<param()-><<", "<<param()-><< "\n";
				if (plant->rand()>f_soil) { // not rand()<p
					age -= dt; // the leaf does not emerge in this time step
				}
			}
		}

		if (age>0) {

			// children first (lateral leafs grow even if base leaf is inactive)
			for (auto c:children) {
				c->simulate(dt,silence);
			}

			if (active) {

				// length increment
				double length_ = LeafGetLength(std::max(age-dt,0.)); // length of the leaf for unimpeded f_gf (i.e. length_==length for unimpeded f_gf)
				double targetlength = LeafGetLength(age);
				double e = targetlength-length_; //elongation in time step dt
				double scale = getLeafRandomParameter()->f_se->getValue(nodes.back(),this); // hope some of this is optimized out if not set
				double dl = std::max(scale*e, double(0)); // length increment, dt is not used anymore

				// create geometry
				if (param()->ln.size()>0) { // leaf has laterals
					// basal zone
					if ((dl>0)&&(length<param()->lb)) { // length is the current length of the leaf
						if (length+dl<=param()->lb) {
							createSegments(dl,silence);
							length+=dl;
							dl=0;
						} else {
							double ddx = param()->lb-length;
							createSegments(ddx,silence);
							dl-=ddx; // ddx already has been created
							length=param()->lb;
						}
					}
					// branching zone Condition will be changed later
					if ((dl>0)&&(length>=param()->lb)) {
						double s = param()->lb; // summed length
						for (size_t i=0; ((i<param()->ln.size()) && (dl>0)); i++) {
							s+=param()->ln.at(i);
							if (length<s) {
								if (i==children.size()) { // new lateral


									createLateral(silence);

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
							if (param()->ln.size()==children.size()) { // new lateral (the last one)
								createLateral(silence);

							}
						}
					}
					// apical zone
					if (dl>0) {
						createSegments(dl,silence);
						length+=dl;
					}
				} else { // no laterals
					if (dl>0) {
						createSegments(dl,silence);
						length+=dl;
					}
				} // if laterals
			} // if active
			active = LeafGetLength(std::max(age,0.))<(param()->getK()-dx()/10); // become inactive, if final length is nearly reached
		}
	} // if alive



}

/**
 *
 */
double Leaf::getParameter(std::string name) const {
double r = Organ::getParameter(name);
	if (name=="basal zone") { r = param()->lb; }
	if (name=="apical zone") { r = param()->la; }
	if (name=="initial f_gf rate") { r = param()->r; }
	if (name=="radius") { r = param()->a; }
	if (name=="insertion angle") { r = param()->theta; }
	if (name=="root life time") { r = param()->rlt; }
//	if (name=="mean internodal distance") {
//		r = std::accumulate(param()->ln.begin(), param()->ln.end(), 0.0) / param()->ln.size();
//	}
//	if (name=="sd internodal distance") {
//		const std::vector<double>& v_ = param()->ln;
//		double mean = std::accumulate(v_.begin(), v_.end(), 0.0) / v_.size();
//		double sq_sum = std::inner_product(v_.begin(), v_.end(), v_.begin(), 0.0);
//		r = std::sqrt(sq_sum / v_.size() - mean * mean);
//	}
//	if (name=="surface") { r = param()->a*param()->a*M_PI*length; }
//	if (name=="number of branches") { r = param()->ln.size() +1; }
	return r;

}

/**
 * Analytical creation (=emergence) time of a node at a length along the leaf
 *
 * @param length   length of the leaf [cm]
 */
double Leaf::getCreationTime(double length)
{
	assert(length>=0);
	double age = LeafGetAge(length);
	assert(age>=0);
	if (parent->organType()!=Organism::ot_stem) {
		if (parent->organType()==Organism::ot_leaf) {
			double pl = pbl+((Leaf*)parent)->getLeafRandomParameter()->la; // parent length, when this leaf was created
			double pAge=((Leaf*)parent)->getCreationTime(pl);
			return age+pAge;
		} else { // organ type is seed
			return age+nodeCTs[0];
		}
	} else {
		return age+nodeCTs[0];
	}
}

/**
 * Analytical length of the leaf at a given age
 *
 * @param age          age of the leaf [day]
 */
double Leaf::LeafGetLength(double age)
{
//    std::string organ_name = std::string(getLeafRandomParameter()->organName);
//    //std::cout<<"organName is "<<name()<<"\n";
//    if (name()  == "maize1eaf"){
        assert(age>=0);
        std::cout<<"organName is ? \n";
        //return getLeafRandomParameter()->f_gf->LeafgetLength(age,getLeafRandomParameter()->r,getLeafRandomParameter()->getK(),this);
        return getLeafRandomParameter()->f_gf->getLength(age,getLeafRandomParameter()->r,getleafphytomerID(getLeafRandomParameter()->subType)*3,this);
//	}else {
//assert(age>=0);
//	    return getLeafRandomParameter()->f_gf->LeafgetLength(age,getLeafRandomParameter()->r,getLeafRandomParameter()->getK(),this);
//	    }
}

/**
 * Analytical age of the leaf at a given length
 *
 * @param length   length of the leaf [cm]
 */
double Leaf::LeafGetAge(double length)
{
//    std::string organ_name = std::string(getLeafRandomParameter()->organName);
//    //std::cout<<getLeafRandomParameter()->name<< "\n";
//     if ( name() == "maize1eaf"){
        assert(length>=0);
        //std::cout<<"length subtype is"<<getLeafRandomParameter()->subType<<"\n";
        return getLeafRandomParameter()->f_gf->getAge(length,getLeafRandomParameter()->r,getLeafRandomParameter()->getK(),this);
//        return getLeafRandomParameter()->f_gf->LeafgetAge(length,getLeafRandomParameter()->r,getleafphytomerID(getLeafRandomParameter()->subType)*3,this);
//        }else {
//assert(age>=0);
//	    return getLeafRandomParameter()->f_gf->LeafgetAge(length,getLeafRandomParameter()->r,getLeafRandomParameter()->getK(),this);
//	    }
}

/**
 *
 */
//LeafRandomParameter* Leaf::getLeafRandomParameter() const {
//	return (LeafRandomParameter*)getLeafRandomParameter();
//}

/**
 *
 */
double Leaf::dx() const
{
	return getLeafRandomParameter()->dx;
}

//std::string Leaf::name() const
//{
//	return ((LeafRandomParameter*)getOrganRandomOrganParameter())->name;
//}

/**
 * Creates a new lateral by calling Leaf::createNewleaf().
 *
 * Overwrite this method to implement more specialized leaf classes.
 */
void Leaf::createLateral(bool silence)
{

    int lt = getLeafRandomParameter()->getLateralType(getNode(nodes.size()-1));

	if (lt>0) {

		if (getLeafRandomParameter()->lnf==2&& lt>0) {
		double ageLN = this->LeafGetAge(length); // age of Leaf when lateral node is created
		double ageLG = this->LeafGetAge(length+param()->la); // age of the Leaf, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Leaf* lateral = new Leaf(plant, lt, h, delay, this, nodes.size() - 1, length);
//		lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		Leaf* lateral2 = new Leaf(plant, lt, h, delay, this, nodes.size() - 1, length);
//		lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
	else if (getLeafRandomParameter()->lnf==3&& lt>0) { //ln equal and both side leaf
		double ageLN = this->LeafGetAge(length); // age of Leaf when lateral node is created
		double ageLG = this->LeafGetAge(length+param()->la); // age of the Leaf, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Leaf* lateral = new Leaf(plant,  lt, h, delay, this, nodes.size() - 1, length);
//		lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		Leaf* lateral2 = new Leaf(plant,  lt, h, delay, this, nodes.size() - 1, length);
//		lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
	else if (getLeafRandomParameter()->lnf==4 && lt>0) {//ln exponential decreasing and one side leaf
		double ageLN = this->LeafGetAge(length); // age of Leaf when lateral node is created
		double ageLG = this->LeafGetAge(length+param()->la); // age of the Leaf, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Leaf* lateral = new Leaf(plant,  lt, h, delay, this, nodes.size() - 1, length);
//		lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

	} else if (getLeafRandomParameter()->lnf==5&& lt>0) { //ln exponential decreasing and both side leaf
		double ageLN = this->LeafGetAge(length); // age of Leaf when lateral node is created
		double ageLG = this->LeafGetAge(length+param()->la); // age of the Leaf, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Leaf* lateral = new Leaf(plant, lt,  h, delay,  this, nodes.size() - 1, length);
//		lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)

		addleafphytomerID(getLeafRandomParameter()->subType);

		Leaf* lateral2 = new Leaf(plant, lt, h, delay,  this, nodes.size() - 1, length);
//		lateral2->setRelativeOrigin(nodes.back());
		children.push_back(lateral2);
		lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	} else if (lt>0) {
		double ageLN = this->LeafGetAge(length); // age of Leaf when lateral node is created
		double ageLG = this->LeafGetAge(length+param()->la); // age of the Leaf, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading

		Leaf* lateral = new Leaf(plant, lt, h, delay, this, nodes.size() - 1, length);
//		lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		}else {


		double ageLN = this->LeafGetAge(length); // age of leaf when lateral node is created
		double ageLG = this->LeafGetAge(length+param()->la); // age of the leaf, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Leaf* lateral = new Leaf(plant, lt, h, delay, this, nodes.size()-1, length);
//		lateral->setRelativeOrigin(nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
}
}



/**
 *  Creates nodes and node emergence times for length l,
 *  and updates the leaf heading
 *
 *  Cecks that each new segments length is <= dx but >= ddx
 *
 *  @param l       length the leaf f_gf [cm]
 */
void Leaf::createSegments(double l, bool silence)
{
	//   std::cout << "create Leaf Segments("<< l << ")\n";
	assert(l>0);
	double sl=0; // summed length of created segment

	// shift first node to axial resolution
	int nn = nodes.size();
	if (old_non==0) { // first call of createSegments (in leaf::simulate)
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
				nodeCTs[nn - 1] = std::max(ct, plant->getSimTime()); // in case of impeded f_gf the node emergence time is not exact anymore, but might break down to temporal resolution
				old_non = nn;
				l -= newdx;
				if (l<=0) { // ==0 should be enough
					return;
				}
			}
		}
		old_non = nn; // CHECK
	}

	if (l<smallDx) {
		if (!silence) {
			//      std::cout << "skipped small segment (<"<< smallDx << ") \n";
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
					//          std::cout << "skipped small segment (<"<< i << ") \n";
				}
				return;
			}
		}
		sl+=sdx;

            Vector3d newdx = getIncrement(nodes.back(), sdx);

			Vector3d newnode = Vector3d(nodes.back().plus(newdx));
		double ct = this->getCreationTime(length+sl);
		ct = std::max(ct,plant->getSimTime()); // in case of impeded f_gf the node emergence time is not exact anymore, but might break down to temporal resolution
		//     std::cout<<"add node "<<newnode.toString()<<"\n";
		addNode(newnode,ct);

	} // for

}

Vector3d Leaf::getIncrement(const Vector3d& p, double sdx) {
    Vector3d h = heading(); // current heading
    Matrix3d ons = Matrix3d::ons(h);
    Vector2d ab = getLeafRandomParameter()->f_tf->getHeading(p, ons, sdx, this);
    Vector3d sv = ons.times(Vector3d::rotAB(ab.x,ab.y));
//    if (rootsystem->poreGeometry==nullptr) { // no pores defined
//        return sv.times(sdx);
//    } else {
//        if (rootsystem->poreGeometry->getDist(p)<0) { // inside the pore
//            auto sv1 = rootsystem->applyPoreConductivities(sv);
//            // std::cout << "Length before " << sv.length() << ", length after " << sv1.length() << "\n";
//            sv1.normalize();
//            return sv1.times(sdx);
//        } else {
            return sv.times(sdx);
//        }
//    }
}


/**
 * Relative heading of organ tip
 */
//Vector3d Leaf::relHeading() const {
//	Vector3d h;
//	if (nodes.size()>1) {
//		h = nodes.back().minus(nodes.at(nodes.size()-2)); // getHeading(b-a)
//		h.normalize();
//	} else {
//		h = getRelativeHeading().column(2);
//	}
//	return h;
//}

/**
 * Absolute heading of organ tip
 */
//Vector3d Leaf::absHeading() const {
//	return  getHeading().column(2);;
//}
Vector3d Leaf::heading() const {
	Vector3d h;
	if (nodes.size()>1) {
		h = nodes.back().minus(nodes.at(nodes.size() - 2)); // getHeading(b-a)
	}
	else {
		h = initialLeafHeading;
	}
	return h;
}

/**
 * @return The RootTypeParameter from the plant
 */
LeafRandomParameter* Leaf::getLeafRandomParameter() const
{
    return (LeafRandomParameter*)plant->getOrganRandomParameter(Organism::ot_root, param_->subType);
}

/**
 * @return Parameters of the specific root
 */
LeafSpecificParameter* Leaf::param() const
{
    return (LeafSpecificParameter*)param_;
}


/**
 * Quick info about the object for debugging
 */
std::string Leaf::toString() const
{
	std::stringstream str;
	str << "Leaf #"<< id <<": type "<<param()->subType << ", length: "<< length << ", age: " <<age<<" with "<< children.size() << " laterals\n";
	return str.str();
}





} // namespace CPlantBox
