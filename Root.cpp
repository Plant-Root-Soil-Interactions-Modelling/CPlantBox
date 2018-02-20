#include "Root.h"

/**
 * Constructor
 *
 * Typically called by the RootSystem::RootSystem(), or Root::createNewRoot().
 * For base roots the initial node and node emergence time (netime) must be set from outside
 *
 * @param rs 			points to RootSystem
 * @param type 		    type of root that is created
 * @param pheading		heading of parent root at emergence
 * @param delay 		to give apical zone of parent time to develop
 * @param parent		parent root
 * @param pni			parent node index
 * @param pbl			parent base length
 */
Root::Root(Plant* plant, Organ* parent, int type, double delay, Vector3d iheading ,int pni, double pbl) :Organ(plant,parent,type,delay), pni(pni), pbl(pbl)
{
	initialHeading = iheading;
//	std::cout << "Root pni = "<< pni<< std::endl;
	//std::cout << "Root constructor \n";

	RootTypeParameter* rtp = (RootTypeParameter*) plant->getParameter(Organ::ot_root, type);
	param = rtp->realize(); // throw the dice
	RootParameter* root_p = (RootParameter*) param;
	double beta = 2*M_PI*plant->rand(); // initial rotation
	Matrix3d ons = Matrix3d::ons(iheading);
	ons.times(Matrix3d::rotX(beta));
	double theta = root_p->theta;
	if (parent->organType()!=Organ::ot_seed) { // scale if not a base root
		double scale = rtp->sa->getValue(parent->getNode(pni),this);
		theta*=scale;
	}
	ons.times(Matrix3d::rotZ(theta));
	this->initialHeading= ons.column(0);
	// initial node
	if (parent->organType()!=Organ::ot_seed) { // the first node of the base roots must be created in Seed::initialize()
		// otherwise, don't use addNode for the first node of the root,
		// since this node exists already and does not need a new identifier
		r_nodes.push_back(parent->getNode(pni));
		nodeIDs.push_back(parent->getNodeID(pni));
		nctimes.push_back(parent->getNodeCT(pni)+delay);
	}
}


/**
 * Simulates growth of this root for a time span dt
 *
 * @param dt       time step [day]
 * @param silence  indicates if status messages are written to the console (cout) (default = false)
 */
void Root::simulate(double dt, bool silence)
{


	old_non = 0; // is set in Root:createSegments, the zero indicates the first call to createSegments

	const RootParameter* rp = rParam(); // rename
	const RootTypeParameter* tp = tParam();

	// increase age
	if (age+dt>rp->rlt) { // root life time
		dt=rp->rlt-age; // remaining life span
		alive = false; // this root is dead
	}
	age+=dt;

	if (alive) { // dead roots wont grow

		// probabilistic branching model (todo test)
		if ((age>0) && (age-dt<=0)) { // the root emerges in this time step
			double root_P = tp->sbp->getValue(r_nodes.back(),this);
			if (root_P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
				double root_p = 1.-std::pow((1.-root_P), dt); //probability of emergence in this time step
				std::cout <<root_P<<", "<<root_p<< "\n";
				if (plant->rand()>root_p) { // not rand()<p
					age -= dt; // the root does not emerge in this time step
				}
			}
		}

		if (age>0) {

			// children first (lateral roots grow even if base root is inactive)
			for (auto c:children) {
				c->simulate(dt,silence);
			}

			if (active) {

				// length increment
				double length_ = getLength(std::max(age-dt,0.)); // length of the root for unimpeded growth (i.e. length_==length for unimpeded growth)
				double targetlength = getLength(age);
				double e = targetlength-length_; //elongation in time step dt
				double scale = tp->se->getValue(r_nodes.back(),this); // hope some of this is optimized out if not set
				double dl = std::max(scale*e, double(0)); // length increment, dt is not used anymore

				// create geometry
				if (rp->ln.size()>0) { // root has laterals
					// basal zone
					if ((dl>0)&&(length<rp->lb)) { // length is the current length of the root
						if (length+dl<=rp->lb) {
							createSegments(dl,silence);
							length+=dl;
							dl=0;
						} else {
							double ddx = rp->lb-length;
							createSegments(ddx,silence);
							dl-=ddx; // ddx already has been created
							length=rp->lb;
						}
					}
					// branching zone
					if ((dl>0)&&(length>=rp->lb)) {
						double s = rp->lb; // summed length
						for (size_t i=0; ((i<rp->ln.size()) && (dl>0)); i++) {
							s+=rp->ln.at(i);
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
							if (rp->ln.size()==children.size()) { // new lateral (the last one)
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
			active = getLength(std::max(age,0.))<(rp->getK()-dx()/10); // become inactive, if final length is nearly reached
		}
	} // if alive


}

/**
 *
 */
double Root::getScalar(std::string name) const {
	double r = Organ::getScalar(name);
	if (name=="basalzone") { r = rParam()->lb; } // todo and so on
//	switch(stype) {
//	// st_rlt, st_meanln, st_stdln , st_nob, st_surface, , // root level
//	case Plant::st_lb:
//		return rParam()->lb;
//	case Plant::st_la:
//		return rParam()->la;
//	case Plant::st_r:
//		return rParam()->r;
//	case Plant::st_radius:
//		return rParam()->a;
//	case Plant::st_theta:
//		return rParam()->theta;
//	case Plant::st_rlt:
//		return rParam()->rlt;
//	case Plant::st_meanln:
//		return std::accumulate(rParam()->ln.begin(), rParam()->ln.end(), 0.0) / rParam()->ln.size();
//	case Plant::st_stdln: {
//		const std::vector<double>& v_ = rParam()->ln;
//		double mean = std::accumulate(v_.begin(), v_.end(), 0.0) / v_.size();
//		double sq_sum = std::inner_product(v_.begin(), v_.end(), v_.begin(), 0.0);
//		return std::sqrt(sq_sum / v_.size() - mean * mean);
//	}
//	case Plant::st_surface:
//		return rParam()->a*rParam()->a*M_PI*length;
//	case Plant::st_nob:
//		return rParam()->ln.size();
//	default:
//		return  Organ::getScalar(stype);
//	}
	return r;
}

/**
 * Analytical creation (=emergence) time of a node at a length along the root
 *
 * @param length   length of the root [cm]
 */
double Root::getCreationTime(double length)
{
	assert(length>=0);
	double age = getAge(length);
	assert(age>=0);
	if (parent->organType()!=Organ::ot_seed) {
		if (parent->organType()==Organ::ot_root) {
			double pl = pbl+((Root*)parent)->rParam()->la; // parent length, when this root was created
			double pAge=((Root*)parent)->getCreationTime(pl);
			return age+pAge;
		} else { // organ type is seed
			return age;
		}
	} else {
		return age;
	}
}

/**
 * Analytical length of the root at a given age
 *
 * @param age          age of the root [day]
 */
double Root::getLength(double age)
{
	assert(age>=0);
	return tParam()->growth->getLength(age,rParam()->r,rParam()->getK(),this);
}

/**
 * Analytical age of the root at a given length
 *
 * @param length   length of the root [cm]
 */
double Root::getAge(double length)
{
	assert(length>=0);
	return tParam()->growth->getAge(length,rParam()->r,rParam()->getK(),this);
}

/**
 *
 */
RootTypeParameter* Root::tParam() const {
	return (RootTypeParameter*)getOrganTypeParameter();
}

/**
 *
 */
double Root::dx() const
{
	return ((RootTypeParameter*)getOrganTypeParameter())->dx;
}

/**
 * Creates a new lateral by calling RootSystem::createNewRoot().
 *
 * Overwrite this method to implement more sezialized root classes.
 */
void Root::createLateral(bool silence)
{
	const RootParameter* rp = rParam(); // rename
	int lt = tParam()->getLateralType(r_nodes.back());
	//	std::cout << "createLateral()\n";
	//	std::cout << "lateral type " << lt << "\n";

	if (lt>0) {
		double ageLN = this->getAge(length); // age of root when lateral node is created
		double ageLG = this->getAge(length+rp->la); // age of the root, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		Vector3d h = heading(); // current heading
		Root* lateral = new Root(plant, this, lt, delay, h,  r_nodes.size()-1, length);
		children.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
	}
}

/**
 *  Creates nodes and node emergence times for length l,
 *  and updates the root heading
 *
 *  Cecks that each new segments length is <= dx but >= ddx
 *
 *  @param l       length the root growth [cm]
 */
void Root::createSegments(double l, bool silence)
{
	// std::cout << "createSegments("<< l << ")\n";
	assert(l>0);
	double sl=0; // summed length of created segment

	// shift first node to axial resolution
	int nn = r_nodes.size();
	if (old_non==0) { // first call of createSegments (in Root::simulate)
		if (nn>1) {
			auto n2 = r_nodes.at(nn-2);
			auto n1 = r_nodes.at(nn-1);
			double olddx = n1.minus(n2).length();
			if (olddx<dx()*0.99) { // shift node instead of creating a new node

				Vector3d h; // current heading
				if (nn>2) {
					h = n2.minus(r_nodes.at(nn-3));
					h.normalize();
				} else {
					h = initialHeading;
				}
				double sdx = std::min(dx()-olddx,l);

				Matrix3d ons = Matrix3d::ons(h);
				Vector2d ab = tParam()->tropism->getHeading(r_nodes.at(nn-2),ons,olddx+sdx,this);
				ons.times(Matrix3d::rotX(ab.y));
				ons.times(Matrix3d::rotZ(ab.x));
				Vector3d newdx = Vector3d(ons.column(0).times(olddx+sdx));

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
//			std::cout << "skipped small segment (<"<< smallDx << ") \n";
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
//					std::cout << "skipped small segment (<"<< smallDx << ") \n";
				}
				return;
			}
		}
		sl+=sdx;

		Vector3d h= heading();
		Matrix3d ons = Matrix3d::ons(h);
		Vector2d ab = tParam()->tropism->getHeading(r_nodes.back(),ons,sdx,this);
		ons.times(Matrix3d::rotX(ab.y));
		ons.times(Matrix3d::rotZ(ab.x));
		Vector3d newdx = Vector3d(ons.column(0).times(sdx));
		Vector3d newnode = Vector3d(r_nodes.back().plus(newdx));
		double ct = this->getCreationTime(length+sl);
		ct = std::max(ct,plant->getSimTime()); // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
		// std::cout<<"add node "<<newnode.toString()<<"\n";
		addNode(newnode,ct);

	} // for

}

Vector3d Root::heading() const {
	Vector3d h;
	if (r_nodes.size()>1) {
		h = r_nodes.back().minus(r_nodes.at(r_nodes.size()-2)); // getHeading(b-a)
	} else {
		h = initialHeading;
	}
	return h;
}

/**
 * Adds the next node to the root.
 *
 * Add nodes only with this function! For simplicity nodes can not be deleted, and roots can only become deactivated by dying
 *
 * @param n        the new node
 * @param t        exact creation time of the node
 */
void Root::addNode(Vector3d n, double t)
{
	assert(t>=0.);
	r_nodes.push_back(n); // node
	nodeIDs.push_back(plant->getNodeIndex()); // new unique id
	nctimes.push_back(t); // exact creation time
}

/**
 * writes RSML root tag
 *
 * @param cout      typically a file out stream
 * @param indent    we care for looks
 */
void Root::writeRSML(std::ostream & cout, std::string indent) const
{
	if (this->r_nodes.size()>1) {
		cout << indent << "<root id=\"" <<  id << "\">\n";  // open root

		/* geometry tag */
		cout << indent << "\t<geometry>\n"; // open geometry
		cout << indent << "\t\t<polyline>\n"; // open polyline
		// polyline nodes
		cout << indent << "\t\t\t" << "<point ";
		Vector3d v = r_nodes.at(0);
		cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
		int n = 5; //this->plant->rsmlReduction;
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

		/* laterals roots */
		for (size_t i = 0; i<children.size(); i++) {
			children[i]->writeRSML(cout,indent+"\t");
		}

		cout << indent << "</root>\n"; // close root
	}
}

/**
 * Quick info about the object for debugging
 */
std::string Root::toString() const
{
	std::stringstream str;
	str << "Root #"<< id <<": type "<<param->subType << ", length: "<< length << ", age: " <<age<<" with "<< children.size() << " laterals\n";
	return str.str();
}


