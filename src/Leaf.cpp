#include "Stem.h"
#include "Leaf.h"
#include "Root.h"
#include "Plant.h"

namespace CPlantBox {

/**
 * Constructs a leaf from given data.
 * The organ tree must be created, @see Organ::setPlant, Organ::setParent, Organ::addChild
 * Organ geometry must be created, @see Organ::addNode, ensure that this->getNodeId(0) == parent->getNodeId(pni)
 *
 * @param id        	the organ's unique id (@see Organ::getId)
 * @param param     	the organs parameters set, ownership transfers to the organ
 * @param alive     	indicates if the organ is alive (@see Organ::isAlive)
 * @param active    	indicates if the organ is active (@see Organ::isActive)
 * @param age       	the current age of the organ (@see Organ::getAge)
 * @param length    	the current length of the organ (@see Organ::getLength)
 * @param iheading  	the initial heading of this leaf
 * @param pbl       	base length of the parent leaf, where this leaf emerges
 * @param pni       	local node index, where this leaf emerges
 * @deprecated moved	as long as stem is active, nodes are assumed to have moved (@see Organism::getUpdatedNodes)
 * @param oldNON    	the number of nodes of the previous time step (default = 0)
 */
Leaf::Leaf(int id, const std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
		Matrix3d iHeading,int pni, bool moved, int oldNON)
		:Organ(id, param, alive, active, age, length, iHeading, pni, moved,  oldNON )
{}

/**
 * Constructor
 * Typically called by the Plant::Plant(), or Leaf::createNewLeaf().
 * For leaf the initial node (in nodes) and node emergence time (in nodeCTs) must be set from outside
 *
 * @param plant 		points to the plant
 * @param parent 		points to the parent organ
 * @param subtype		sub type of the leaf
 * @param delay 		delay after which the organ starts to develop (days)
 * @param rheading		relative heading (within parent organ)
 * @param pni			parent node index
 * @param pbl			parent base length
 */
Leaf::Leaf(std::shared_ptr<Organism> plant, int type, Matrix3d iHeading, double delay,  std::shared_ptr<Organ> parent, int pni)
:Organ(plant, parent, Organism::ot_leaf, type, delay, iHeading, pni)
{
	assert(parent!=nullptr && "Leaf::Leaf parent must be set");
	addleafphytomerID(param()->subType);
	ageDependentTropism = getLeafRandomParameter()->f_tf->ageSwitch > 0;
	// Calculate the rotation of the leaves. The code begins here needs to be rewritten, because another following project will work on the leaves. The code here is just temporally used to get some nice visualizations. When someone rewrites the code, please take "gimbal lock" into consideration.  
	//Rewritten Begin: 															 
	double beta = getleafphytomerID(param()->subType)*M_PI*getLeafRandomParameter()->rotBeta
			+ M_PI*plant->rand()*getLeafRandomParameter()->betaDev ;  //+ ; //2 * M_PI*plant->rand(); // initial rotation
	beta = beta + getLeafRandomParameter()->initBeta*M_PI;
	if (getLeafRandomParameter()->initBeta >0 && getLeafRandomParameter()->subType==2 && getLeafRandomParameter()->lnf==5 && getleafphytomerID(2)%4==2) {
		beta = beta + getLeafRandomParameter()->initBeta*M_PI;
	} else if (getLeafRandomParameter()->initBeta >0 && getLeafRandomParameter()->subType==2 && getLeafRandomParameter()->lnf==5 && getleafphytomerID(2)%4==3) {
		beta = beta + getLeafRandomParameter()->initBeta*M_PI + M_PI;
	}
	double theta = param()->theta;
	if (parent->organType()!=Organism::ot_seed) { // scale if not a base leaf
		double scale = getLeafRandomParameter()->f_sa->getValue(parent->getNode(pni), parent);
		theta *= scale;
	}
	//used when computing actual heading, @see LEaf::getIHeading
	this->partialIHeading = Vector3d::rotAB(theta,beta);
	// Rewritten ends 
	if (parent->organType()!=Organism::ot_seed) { // if not base organ
	
		double creationTime;
		if (parent->organType()==Organism::ot_stem) {
			//if lateral of stem, initial creation time: 
			//time when stem reached end of basal zone (==CT of parent node of first lateral) + delay
			// @see stem::leafGrow
			if (parent->getNumberOfChildren() == 0){creationTime = parent->getNodeCT(pni)+delay;
			}else{creationTime = parent->getChild(0)->getParameter("creationTime") + delay;}
		}else{
			creationTime = parent->getNodeCT(pni)+delay;
		}
		addNode(Vector3d(0.,0.,0.), parent->getNodeId(pni), creationTime);//create first node. relative coordinate = (0,0,0)
}}

/**
 * Deep copies the organ into the new plant @param rs.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
std::shared_ptr<Organ> Leaf::copy(std::shared_ptr<Organism> p)
{
	auto l = std::make_shared<Leaf>(*this); // shallow copy
	l->parent = std::weak_ptr<Organ>();
	l->plant = p;
	l->param_ = std::make_shared<LeafSpecificParameter>(*param()); // copy parameters
	for (size_t i=0; i< children.size(); i++) {
		l->children[i] = children[i]->copy(p); // copy laterals
		l->children[i]->setParent(l);
	}
	return l;
}

/**
 * Simulates f_gf of this leaf for a time span dt
 *
 * @param dt       time step [day]
 * @param verbose  indicates if status messages are written to the console (cout) (default = false)
 */
void Leaf::simulate(double dt, bool verbose)
{
	firstCall = true;
	oldNumberOfNodes = nodes.size();

	const LeafSpecificParameter& p = *param(); // rename

	if (alive) { // dead leafs wont grow

		// increase age
		if (age+dt>p.rlt) { // leaf life time
			dt=p.rlt-age; // remaining life span
			alive = false; // this leaf is dead
		}
		age+=dt;

		// probabilistic branching model (todo test)
		if ((age>0) && (age-dt<=0)) { // the leaf emerges in this time step
			//currently, does not use absolute coordinates for these function. 
			double P = getLeafRandomParameter()->f_sbp->getValue(nodes.back(),shared_from_this());
			if (P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
				double p = 1.-std::pow((1.-P), dt); //probability of emergence in this time step
				if (plant.lock()->rand()>p) { // not rand()<p
					age -= dt; // the leaf does not emerge in this time step
				}
			}
		}

		if (age>0) { // unborn leafs have no children

			// children first (lateral leafs grow even if base leaf is inactive)
			for (auto l:children) {
				l->simulate(dt,verbose);
			}

			if (active) {

				// length increment
				double age_ = calcAge(length); // leaf age as if grown unimpeded (lower than real age)
				double dt_; // time step
				if (age<dt) { // the leaf emerged in this time step, adjust time step
					dt_= age;
				} else {
					dt_=dt;
				}

				double targetlength = calcLength(age_+dt_)+ this->epsilonDx;
				double e = targetlength-length; // unimpeded elongation in time step dt
				double dl = std::max(e, 0.);// length increment = calculated length + increment from last time step too small to be added
				length = getLength();
				this->epsilonDx = 0.; // now it is "spent" on targetlength (no need for -this->epsilonDx in the following)
				// create geometry
				if (p.laterals) { // leaf has laterals
					/* basal zone */
					if ((dl>0)&&(length<p.lb)) { // length is the current length of the leaf
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
					double s = p.lb; // summed length
					/* branching zone */
					if ((dl>0)&&(length>=p.lb)) {
						for (size_t i=0; ((i<p.ln.size()) && (dl>0)); i++) {
							s+=p.ln.at(i);
							if (length<s) {
								if (i==children.size()) { // new lateral
									createLateral(verbose);
								}
								if (length+dl<=s) { // finish within inter-lateral distance i
									createSegments(dl,verbose);
									length+=dl;//- this->epsilonDx;
									dl=0;
								} else { // grow over inter-lateral distance i
									double ddx = s-length;
									createSegments(ddx,verbose);
									dl-=ddx;
									length=s;
								}
							}
						}
						if (p.ln.size()==children.size()&& (getLength(true)>=s)) { // new lateral (the last one)
							createLateral(verbose);
						}
					}
					/* apical zone */
					if (dl>0) {
						createSegments(dl,verbose);//y not with dt_?
						length+=dl;//- this->epsilonDx;
					}
				} else { // no laterals
					if (dl>0) {
						createSegments(dl,verbose);
						length+=dl;//- this->epsilonDx;
					}
				} // if lateralgetLengths
			} // if active
			//level of precision = 1e-10 to not create an error in the test files
			active = getLength(false)<=(p.getK()*(1 - 1e-11)); // become inactive, if final length is nearly reached
		}
	} // if alive

}

/**
 *
 */
double Leaf::getParameter(std::string name) const {
	if (name=="lb") { return param()->lb; } // basal zone [cm]
	if (name=="la") { return param()->la; } // apical zone [cm]
	//if (name=="nob") { return param()->nob; } // number of branches
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
	if (name=="volume") { return param()->a*param()->a*M_PI*getLength(true); } // // realized root volume [cm^3]
	if (name=="surface") { return 2*param()->a*M_PI*getLength(true); } // // realized root surface [cm^2]
	if (name=="type") { return this->param_->subType; }  // in CPlantBox the subType is often called just type
	if (name=="parentNI") { return parentNI; } // local parent node index where the lateral emerges
	return Organ::getParameter(name);
}


/**
 * in case there are no lateral leafs return leaf surface area [cm2]
 */
double Leaf::leafArea()
{
	if (param()->laterals) {
		return 0.;
	} else {
		return param()->areaMax * (leafLength()/param()->leafLength());
	}
};

/**
 * indicates if the node is in the leaf surface are and should be viusalized as polygon
 *
 * leaf base (false), branched leaf (false), or leaf surface area (true)
 */
bool Leaf::nodeLeafVis(double l)
{
	if (param()->laterals) {
		return false;
	} else {
		return l >= param()->lb; // true if not in basal zone
	}
}

/**
 * Parameterization x value, at position l along the leaf axis
 */
std::vector<double> Leaf::getLeafVisX_(double l) {
	auto& lg = getLeafRandomParameter()->leafGeometry;
	int n = lg.size();
	int ind = int( ((l - param()->lb) /leafLength())*(n-1) + 0.5); // index within precomputed normalized geometry
	auto x_ = lg.at(ind); // could be more than one point for non-convex geometries
	return x_;
}

/**
 * for Python binding
 */
std::vector<double> Leaf::getLeafVisX(int i) {
	return getLeafVisX_(getLength(i));
}

/**
 * Scales unit leaf shape to the specific leaf,
 * and returns leaf shape coordinates per node (normally 2 points, for convex domain it could be more points)
 * see used by vtk_plot.py to create a polygon representation of the leaf area
 */
std::vector<Vector3d> Leaf::getLeafVis(int i)
{
	double l = getLength(i);
	if (nodeLeafVis(l)) {
		auto& lg = getLeafRandomParameter()->leafGeometry;
		int n = lg.size();
		if (n>0) {
			std::vector<Vector3d> coords;
			auto x_ = getLeafVisX_(l);
			Vector3d x1= iHeading.column(0);
			x1.normalize();
			Vector3d y1 = Vector3d(0,0,-1).cross(x1); // todo angle between leaf - halfs
			y1.normalize();
			double a  = leafArea() / leafLength(); // scale radius
			for (double x :x_) {
				coords.push_back(getNode(i).plus(y1.times(x*a)));
			}
			for (double x :x_) {
				coords.push_back(getNode(i).minus(y1.times(x*a)));
			}
			return coords;
		} else {
			std::cout << "Leaf::getLeafVis: WARNING leaf geometry was not set \n";
			return std::vector<Vector3d>();
		}
	} else { // no need for polygonal visualisation
		return std::vector<Vector3d>();
	}
}

/**
 * Analytical creation (=emergence) time of a node at a length along the leaf
 *
 * @param length   length of the leaf [cm]
 */
double Leaf::calcCreationTime(double length)
{
	assert(length >= 0 && "Leaf::getCreationTime() negative length");
	double leafage = calcAge(length);
	leafage = std::min(leafage, age);
	assert(leafage >= 0 && "Leaf::getCreationTime() negative leaf age");
	return leafage+nodeCTs[0];
}

/**
 * Analytical length of the leaf at a given age
 *
 * @param age          age of the leaf [day]
 */
double Leaf::calcLength(double age)
{
	assert(age>=0  && "Leaf::calcLength() negative root age");
	return getLeafRandomParameter()->f_gf->getLength(age,getLeafRandomParameter()->r,param()->getK(),shared_from_this());
}

/**
 * Analytical age of the leaf at a given length
 *
 * @param length   length of the leaf [cm]
 */
double Leaf::calcAge(double length)
{
	assert(length>=0 && "Leaf::calcAge() negative root length");
	return getLeafRandomParameter()->f_gf->getAge(length,getLeafRandomParameter()->r,param()->getK(),shared_from_this());
}

/**
 *
 */
void Leaf::minusPhytomerId(int subtype)
{
	getPlant()->leafphytomerID[subtype]--;
}

/**
 *
 */
int Leaf::getleafphytomerID(int subtype)
{
	return getPlant()->leafphytomerID[subtype];
}

/**
 *
 */
void Leaf::addleafphytomerID(int subtype)
{
	getPlant()->leafphytomerID.at(subtype)++;
}

/**
 * Creates a new lateral by calling Leaf::createNewleaf().
 *
 * Overwrite this method to implement more specialized leaf classes.
 * 
 * This method was done in a rush, because there will be a following project specifically
 * focus on the leaves. This should be rewritten instead of kept. 
 * 
 */
void Leaf::createLateral(bool silence)
{

	int lt = getLeafRandomParameter()->getLateralType(getNode(nodes.size()-1));

	if (lt>0) {

		int lnf = getLeafRandomParameter()->lnf;
		double ageLN = this->calcAge(getLength(true)); // age of Leaf when lateral node is created
		double meanLn = getLeafRandomParameter()->ln; // mean inter-lateral distance
		double effectiveLa = std::max(param()->la-meanLn/2, 0.); // effective apical distance, observed apical distance is in [la-ln/2, la+ln/2]
		double ageLG = this->calcAge(getLength(true)+effectiveLa); // age of the Leaf, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Matrix3d h = Matrix3d(); //heading not need anymore
		if (lnf==2&& lt>0) {
			auto lateral = std::make_shared<Leaf>(plant.lock(), lt, h, delay, shared_from_this(),nodes.size() - 1);
			children.push_back(lateral);
			lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
			auto lateral2 = std::make_shared<Leaf>(plant.lock(), lt, h, delay, shared_from_this(),  nodes.size() - 1);
			children.push_back(lateral2);
			lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		} else if (lnf==3&& lt>0) { //ln equal and both side leaf
			auto lateral = std::make_shared<Leaf>(plant.lock(),  lt, h, delay, shared_from_this(),  nodes.size() - 1);
			children.push_back(lateral);
			lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
			auto lateral2 = std::make_shared<Leaf>(plant.lock(),  lt, h, delay, shared_from_this(),  nodes.size() - 1);
			children.push_back(lateral2);
			lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		} else if (lnf==4 && lt>0) {//ln exponential decreasing and one side leaf
			auto lateral = std::make_shared<Leaf>(plant.lock(),  lt, h, delay, shared_from_this(), nodes.size() - 1);
			children.push_back(lateral);
			lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		} else if (lnf==5&& lt>0) { //ln exponential decreasing and both side leaf
			auto lateral = std::make_shared<Leaf>(plant.lock(), lt,  h, delay,  shared_from_this(), nodes.size() - 1);
			children.push_back(lateral);
			lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
			addleafphytomerID(getLeafRandomParameter()->subType);
			auto lateral2 = std::make_shared<Leaf>(plant.lock(), lt, h, delay,  shared_from_this(), nodes.size() - 1);
			children.push_back(lateral2);
			lateral2->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		} else if (lt>0) {
			auto lateral = std::make_shared<Leaf>(plant.lock(), lt, h, delay, shared_from_this(), nodes.size() - 1);
			children.push_back(lateral);
			lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		} else {
			auto lateral = std::make_shared<Leaf>(plant.lock(), lt, h, delay, shared_from_this(), nodes.size() - 1);
			children.push_back(lateral);
			lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		}
	}
}

/**
 * computes absolute coordinates from relative coordinates
 * when this function is called, the parent organ has already
 * its absolute coordinates
 * called by @see Plant::rel2abs
 */
void Leaf::rel2abs() 
{
	double ageSwitch = getLeafRandomParameter()->f_tf->ageSwitch; //rename
	bool tropismChange = (age > ageSwitch);
	nodes[0] = getOrigin();//get absolute coordinates of first node via coordinates of parent
	for(size_t i=1; i<nodes.size(); i++){
		Vector3d newdx = nodes[i];
		//if new node or has an age-dependent tropism + reached age at which tropism changes. Might need to update the conditions if do new tropism functions
		//i.e., gradual change according to age
		if((i>= oldNumberOfNodes )|| (this->ageDependentTropism&& tropismChange)){
			double sdx = nodes[i].length();
			newdx = getIncrement(nodes[i-1], sdx, i-1);
		}
		nodes[i] = nodes[i-1].plus(newdx);
		
	}
	if(this->ageDependentTropism && tropismChange){
		this->ageDependentTropism = false; //switch done
	}
	for(size_t i=0; i<children.size(); i++){
		(children[i])->rel2abs();
	}//if carries children, update their coordinates from relative to absolute
	
	
}

/**
 * computes relative coordinates from absolute coordinates
 * when this function is called, the parent organ has already
 * its relative coordinates
 * called by @see Plant::abs2rel
 */
void Leaf::abs2rel()
{
	for (int j = nodes.size(); j>1; j--) {
		nodes[j-1] = nodes.at(j-1).minus(nodes.at(j-2));
		}
	nodes[0] = Vector3d(0.,0.,0.);
	for(size_t i=0; i<children.size(); i++){
		(children[i])->abs2rel();
	}//if carry children, update their pos
	
}

/**
 * Returns the increment of the next segments
 *
 *  @param p       coordinates of previous node (in absolute coordinates)
 *  @param sdx     length of next segment [cm]
 *  @param n       index at which new node is to be inserted
 *  @return        the vector representing the increment
 */
Vector3d Leaf::getIncrement(const Vector3d& p, double sdx, int n)
{
	Vector3d h = heading(n);
	Matrix3d ons = Matrix3d::ons(h);
	//use dx() rather rhan sdx to compute heading
	//to make tropism independante from growth rate
	Vector2d ab = getLeafRandomParameter()->f_tf->getHeading(p, ons, dx(),shared_from_this());
	Vector3d sv = ons.times(Vector3d::rotAB(ab.x,ab.y));
	return sv.times(sdx);
}


/**
 * @return Current absolute heading of the organ at node n, based on initial heading, or segment before
 */
Vector3d Leaf::heading(int n) const
{
	if(n<0){n=nodes.size()-1 ;}
	if ((nodes.size()>1)&&(n>0)) {
		n = std::min(int(nodes.size()),n);
		Vector3d h = getNode(n).minus(getNode(n-1));
		h.normalize();
		return h;
	} else {
		return getiHeading();
	}
}


/**
 * @return Current absolute heading of the organ at node n, based on initial heading, or segment before
 */
Vector3d Leaf::getiHeading()  const
{
	Vector3d vIHeading = getParent()->heading(parentNI);
	Matrix3d iHeading = Matrix3d::ons(vIHeading);
	auto heading = iHeading.column(0);
	Vector3d new_heading = Matrix3d::ons(heading).times(this->partialIHeading);
	return Matrix3d::ons(new_heading).column(0);
}


/**
 *  Creates nodes and node emergence times for a length l
 *
 *  Checks that each new segments length is <= dx but >= smallDx
 *
 *  @param l        total length of the segments that are created [cm]
 *  @param verbose  turns console output on or off
 */
void Leaf::createSegments(double l, bool verbose)
{
	if (l==0) {
		std::cout << "Leaf::createSegments: zero length encountered \n";
		return;
	}
	if (l<0) {
		std::cout << "Leaf::createSegments: negative length encountered \n";
	}

	// shift first node to axial resolution
	double shiftl = 0; // length produced by shift
	int nn = nodes.size();
	if (firstCall) { // first call of createSegments (in Leaf::simulate)
		firstCall = false;
		if (nn>1){ // don't move first node of organ
			Vector3d h = nodes[nn-1];
			double olddx = h.length(); // length of last segment
			if (olddx<dx()*0.99) { // shift node instead of creating a new node
				shiftl = std::min(dx()-olddx, l);
				double sdx = olddx + shiftl; // length of new segment
				h.normalize();  
				nodes[nn-1] = h.times(sdx);
				double et = this->calcCreationTime(getLength(true)+shiftl);
				nodeCTs[nn-1] = et; // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
				l -= shiftl;
				if (l<=0) { // ==0 should be enough
					return;
				}
			}
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
			if (sdx<dxMin()*0.99) { // quit if l is too small
				if (verbose&& sdx != 0) {
					std::cout <<  "Leaf::createSegments(): length increment below dxMin threshold ("<< sdx <<" < "<< dxMin() << ") and kept in memory\n";
				}
				this->epsilonDx = sdx;
				return;
			}
			this->epsilonDx = 0; //no residual
		}
		sl += sdx;
		Vector3d newnode = Vector3d(sdx, 0., 0.);//set relative position (@see rel2abs for addition of tropism)
		double et = this->calcCreationTime(getLength(true)+shiftl+sl);
		// in case of impeded growth the node emergence time is not exact anymore,
		// but might break down to temporal resolution
		addNode(newnode, et);
	}
}


/**
 * @return The LeafTypeParameter from the plant
 */
std::shared_ptr<LeafRandomParameter> Leaf::getLeafRandomParameter() const
{
	return std::static_pointer_cast<LeafRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_leaf, param_->subType));
}

/**
 * @return Parameters of the specific leaf
 */
std::shared_ptr<const LeafSpecificParameter>  Leaf::param() const
{
	return std::static_pointer_cast<const LeafSpecificParameter>(param_);
}

/**
 * Quick info about the object for debugging
 * additionally, use getParam()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Leaf::toString() const
{
	std::stringstream newstring;
	newstring << "; initial heading: " << iHeading.column(0).toString()  << ", parent node index " << parentNI << ".";
	return  Organ::toString()+newstring.str();
}
} // namespace CPlantBox
