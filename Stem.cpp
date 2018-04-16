#include "Stem.h"
#include "Leaf.h"
#include "Root.h"
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
Stem::Stem(Plant* plant, Organ* parent, int subtype, double delay, Vector3d rheading ,int pni, double pbl) :Organ(plant,parent,subtype,delay), pni(pni), pbl(pbl)
{

	//  initialStemHeading = isheading;
	//  std::cout << "stem pni = "<< pni<< std::endl;
	//  std::cout << "Stem constructor \n";
	StemTypeParameter* sttp = (StemTypeParameter*) plant->getParameter(Organ::ot_stem, subtype);
	param = sttp->realize(); // throw the dice
	StemParameter* stem_p = (StemParameter*) param;
	//  std::cout <<", "<<(StemParameter*) param<< "\n";
          std::cout<<"theta "<<stem_p->theta<<"\n";
	Matrix3d heading = Matrix3d::ons(rheading); // isheading is the z direction, i.e. column 2 in the matrix

	double beta = M_PI*plant->getSTPIndex()*0.5 ;//0.25*M_PI;//  +  initial rotation M_PI*plant->getSTPIndex()  +
	Matrix3d rotX = Matrix3d::rotX(beta);
	double theta = M_PI*stem_p->theta;
	Matrix3d rotZ = Matrix3d::rotY(theta);
	
	
	heading.times(rotX);
	//parent->setRelativeHeading(rotX); // now the parent is rotating, so the beta is working as before
	heading.times(rotZ);

	setRelativeHeading(heading);

	// initial node
	if (parent->organType()!=Organ::ot_seed) { // the first node of the base stems must be created in Seed::initialize()
		// otherwise, don't use addNode for the first node of the stem,
		// since this node exists already and does not need a new identifier
		r_nodes.push_back(Vector3d());
		nodeIDs.push_back(parent->getNodeID(pni));
		nctimes.push_back(parent->getNodeCT(pni)+delay);
	}

}

/**
 * Simulates growth of this stem for a time span dt
 *
 * @param dt       time step [day]
 * @param silence  indicates if status messages are written to the console (cout) (default = false)
 */
void Stem::simulate(double dt, bool silence)
{
	//  vector
	old_non = 0; // is set in Stem:createSegments, the zero indicates the first call to createSegments

	const StemParameter* sp = sParam(); // rename
	const StemTypeParameter* sttp = stParam();

	// increase age
	if (age+dt>sp->rlt) { // stem life time
		dt=sp->rlt-age; // remaining life span
		alive = false; // this stem is dead
	}
	age+=dt;

	if (alive) { // dead stem wont grow
		// probabilistic branching model (todo test)
		if ((age>0) && (age-dt<=0)) { // the stem emerges in this time step
			double stem_P = sttp->sbp->getValue(r_nodes.back(),this);
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
				double scale = sttp->se->getValue(r_nodes.back(),this); // hope some of this is optimized out if not set
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
								if (i==Organ::getChildren(2).size()) { // new internode leaf and shootBorneRoot
									//if (sp->subType==3) { //this decide which successor grow the leaf (TODO) adding it to parameterfile
									//	//LeafGrow(silence, newnode);
									//	//ShootBorneRootGrow(silence);
									//} else {
										createLateral(silence);
										//LeafGrow(silence);
									//}
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
							if (sp->ln.size()== Organ::getChildren(2).size()) { // new lateral (the last one)
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
							LeafGrow(silence, r_nodes.back());
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
double Stem::getScalar(std::string name) const {
double r = Organ::getScalar(name);
	if (name=="basal zone") { r = stParam()->lb; }
	if (name=="apical zone") { r = stParam()->la; }
	if (name=="initial growth rate") { r = stParam()->r; }
	if (name=="radius") { r = stParam()->a; }
	if (name=="insertion angle") { r = stParam()->theta; }
	if (name=="root life time") { r = stParam()->rlt; }
//	if (name=="mean internodal distance") {
//		r = std::accumulate(stParam()->ln.begin(), stParam()->ln.end(), 0.0) / stParam()->ln.size();
//	}
//	if (name=="sd internodal distance") {
//		const std::vector<double>& v_ = stParam()->ln;
//		double mean = std::accumulate(v_.begin(), v_.end(), 0.0) / v_.size();
//		double sq_sum = std::inner_product(v_.begin(), v_.end(), v_.begin(), 0.0);
//		r = std::sqrt(sq_sum / v_.size() - mean * mean);
//	}
//	if (name=="surface") { r = stParam()->a*stParam()->a*M_PI*length; }
//	if (name=="number of branches") { r = stParam()->ln.size() +1; }
	return r;
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
	if (parent->organType()!=Organ::ot_seed) {
		if (parent->organType()==Organ::ot_stem) {
			double pl = pbl+((Stem*)parent)->stParam()->la; // parent length, when this stem was created
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
	return stParam()->growth->StemgetLength(age,stParam()->r,stParam()->getK(),this);
}

/**
 * Analytical age of the stem at a given length
 *
 * @param length   length of the stem [cm]
 */
double Stem::StemGetAge(double length)
{
	assert(length>=0);
	return stParam()->growth->StemgetAge(length,stParam()->r,stParam()->getK(),this);
}

/**
 *
 */
StemTypeParameter* Stem::stParam() const {
	return (StemTypeParameter*)getOrganTypeParameter();
}

/**
 *
 */
double Stem::dx() const
{
	return ((StemTypeParameter*)getOrganTypeParameter())->dx;
}

/**
 * Creates a new lateral by calling Stem::createNewstem().
 *
 * Overwrite this method to implement more specialized stem classes.
 */
void Stem::createLateral(bool silence)
{
	const StemParameter* sp = sParam(); // rename
	int lt = stParam()->getLateralType(getNode(r_nodes.size()-1));

	if (lt>0) {
		double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = absHeading(); // current heading
		Stem* lateral = new Stem(plant, this, lt, delay, h,  r_nodes.size()-1, length);
		lateral->setRelativeOrigin(r_nodes.back());
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
}



void Stem::LeafGrow(bool silence, Vector3d bud)
{

	const StemParameter* sp = sParam(); // rename
	//int lt = stParam()->getLateralType(getNode(r_nodes.size()-1));
	//  	std::cout << "LeafGrow createLateral()\n";
	//  	std::cout << "LeafGrow type " << lt << "\n";

	//if (lt>0) {
		double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = absHeading(); // current heading
		Leaf* LeafGrow = new Leaf(plant, this , 3, 1, h, r_nodes.size()-1, length);
		//    LeafGrow->addNode(getNode(r_nodes.size()-1), length);
		LeafGrow->setRelativeOrigin(bud);
		this->children.push_back(LeafGrow);
		LeafGrow->simulate(age-ageLN, silence);// pass time overhead (age we want to achieve minus current
	//}


}



void Stem::ShootBorneRootGrow(bool silence)
{

	const StemParameter* sp = sParam(); // rename
	int lt = stParam()->getLateralType(getNode(r_nodes.size()-1));
	//    std::cout << "ShootBorneRootGrow createLateral()\n";
	//    std::cout << "ShootBorneRootGrow lateral type " << lt << "\n";

	if (lt>0) {
		double ageLN = this->StemGetAge(length); // age of stem when lateral node is created
		double ageLG = this->StemGetAge(length+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		int NodeToGrowShotBorneRoot = 2 ;
		Vector3d sbrheading(0,0,-1); //just a test heading
		Root* ShootBorneRootGrow = new Root(plant, this , 5, 0., sbrheading ,NodeToGrowShotBorneRoot, length);
		if (r_nodes.size() > NodeToGrowShotBorneRoot ) {
			//                                ShootBorneRootGrow->addNode(getNode(NodeToGrowShotBorneRoot), length);
			children.push_back(ShootBorneRootGrow);
			ShootBorneRootGrow->simulate(age-ageLN,silence);// pass time overhead (age we want to achieve minus current age)
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
	int nn = r_nodes.size();
	if (old_non==0) { // first call of createSegments (in stem::simulate)
		if (nn>1) {
			auto n2 = r_nodes.at(nn-2);
			auto n1 = r_nodes.at(nn-1);
			double olddx = n1.minus(n2).length();
			if (olddx<dx()*0.99) { // shift node instead of creating a new node

				double sdx = std::min(dx()-olddx,l);

				Vector2d ab = stParam()->tropism->getHeading(getNode(nn-2),getHeading(),olddx+sdx,this);

				Vector3d h = absHeading();
				Matrix3d heading = Matrix3d::ons(h);
				heading.times(Matrix3d::rotX(ab.y));
				heading.times(Matrix3d::rotZ(ab.x));
				Vector3d newdx = Vector3d(heading.column(0).times(olddx+sdx));

				Vector3d newnode = Vector3d(r_nodes.at(nn-2).plus(newdx));
				sl = sdx;
				double ct = this->getCreationTime(length+sl);
				r_nodes[nn-1] = newnode;
				nctimes[nn-1] = std::max(ct,plant->getSimTime()); // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
				old_non = nn-1;
				l -= sdx;

				if (l<=0) { // ==0 should be enough
									return;
				}
			}
		}
		

		old_non = nn; // CHECK
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

		Vector2d ab = stParam()->tropism->getHeading(getNode(nn-1),getHeading(),sdx,this);
		Vector3d h= absHeading();
		Matrix3d heading = Matrix3d::ons(h);
		heading.times(Matrix3d::rotX(ab.y));
		heading.times(Matrix3d::rotZ(ab.x));
		Vector3d newdx = Vector3d(heading.column(0).times(sdx));
		Vector3d newnode = Vector3d(r_nodes.back().plus(newdx));
		double ct = this->getCreationTime(length+sl);
		ct = std::max(ct,plant->getSimTime()); // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
		//     std::cout<<"add node "<<newnode.toString()<<"\n";
		addNode(newnode,ct);
			
	} // for

}

/**
 * Relative heading of organ tip
 */
Vector3d Stem::relHeading() const {
	Vector3d h;
	if (r_nodes.size()>1) {
		h = r_nodes.back().minus(r_nodes.at(r_nodes.size()-2)); // getHeading(b-a)
		h.normalize();
	} else {
		h = getRelativeHeading().column(2);
	}
	return h;
}

/**
 * Absolute heading of organ tip
 */
Vector3d Stem::absHeading() const {
	return  getHeading().column(2);;
}

/**
 * Adds the next node to the stem.
 *
 * Add nodes only with this function! For simplicity nodes can not be deleted, and stems can only become deactivated by dying
 *
 * @param n        the new node
 * @param t        exact creation time of the node
 */
void Stem::addNode(Vector3d n, double t)
{
	assert(t>=0.);
	r_nodes.push_back(n); // node
	nodeIDs.push_back(plant->getNodeIndex()); // new unique id
	nctimes.push_back(t); // exact creation time
}

/**
 * writes RSML stem tag
 *
 * @param cout      typically a file out stream
 * @param indent    we care for looks
 */
void Stem::writeRSML(std::ostream & cout, std::string indent) const
{
	if (this->r_nodes.size()>1) {
		cout << indent << "<stem id=\"" <<  id << "\">\n";  // open stem

		/* geometry tag */
		cout << indent << "\t<geometry>\n"; // open geometry
		cout << indent << "\t\t<polyline>\n"; // open polyline
		// polyline nodes
		cout << indent << "\t\t\t" << "<point ";
		Vector3d v = r_nodes.at(0);
		cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
		int n = 5; // this->plant->rsmlReduction;
		for (size_t i = 1; i<r_nodes.size()-1; i+=n) {
			cout << indent << "\t\t\t" << "<point ";
			Vector3d v = r_nodes.at(i);
			cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
		}
		cout << indent << "\t\t\t" << "<point ";
		v = r_nodes.at(r_nodes.size()-1);
		cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
		cout << indent << "\t\t</polyline>\n"; // close polyline
		cout << indent << "\t</geometry>\n"; // close geometry

		/* properties */
		cout << indent <<"\t<properties>\n"; // open properties
		// TODO
		cout << indent << "\t</properties>\n"; // close properties

		cout << indent << "\t<functions>\n"; // open functions
		cout << indent << "\t\t<function name='emergence_time' domain='polyline'>\n"; // open functions
		cout << indent << "\t\t\t" << "<sample>" << nctimes.at(0) << "</sample>\n";
		for (size_t i = 1; i<nctimes.size()-1; i+=n) {
			cout << indent << "\t\t\t" << "<sample>" << nctimes.at(i) << "</sample>\n";

		}
		cout << indent << "\t\t\t" << "<sample>" << nctimes.at(nctimes.size()-1) << "</sample>\n";

		cout << indent << "\t\t</function>\n"; // close functions
		cout << indent << "\t</functions>\n"; // close functions

		/* laterals stems */
		for (size_t i = 0; i<children.size(); i++) {
			children[i]->writeRSML(cout,indent+"\t");
		}

		cout << indent << "</root>\n"; // close stem
	}
}

/**
 * Quick info about the object for debugging
 */
std::string Stem::toString() const
{
	std::stringstream str;
	str << "Root #"<< id <<": type "<<param->subType << ", length: "<< length << ", age: " <<age<<" with "<< children.size() << " laterals\n";
	return str.str();
}


