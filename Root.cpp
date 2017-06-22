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
 * @param pbl			parent base length
 * @param pni			parent node index
 */
Root::Root(Plant* rs, int type, Vector3d pheading, double delay,  Organ* parent, double pbl, int pni)
{
	//std::cout << "Root constructor \n";
	rootsystem=rs; // remember
	param = rs->getRootTypeParameter(type)->realize(); // throw the dice
	double beta = 2*M_PI*rs->rand(); // initial rotation
	Matrix3d ons = Matrix3d::ons(pheading);
	ons.times(Matrix3d::rotX(beta));
	double theta = param.theta;
	if (parent!=nullptr) { // scale if not a baseRoot
		double scale = rs->getRootTypeParameter(type)->sa->getValue(parent->getNode(pni),this);
		theta*=scale;
	}
	ons.times(Matrix3d::rotZ(theta));
	this->iheading = ons.column(0);  // new initial heading
	//
	age = -delay; // the root starts growing when age>0
	alive = 1; // alive per default
	id = rs->getRootIndex(); // root id
	this->parent = parent;
	parent_base_length=pbl;
	parent_ni=pni;
	length = 0;
	// initial node
	if (parent!=nullptr) { // the first node of the base roots must be created in RootSystem::initialize()
		// otherwise, don't use addNode for the first node of the root,
		// since this node exists already and does not need a new identifier
		nodes.push_back(parent->getNode(pni));
		nodeIds.push_back(parent->getNodeId(pni));
		netimes.push_back(parent->getNodeETime(pni)+delay);
	}
}

/**
 * Destructor, spread the word
 */
Root::~Root()
{
	for(auto l : laterals) {
		delete l;
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

	const RootParameter &p = param; // rename

	// increase age
	if (age+dt>p.rlt) { // root life time
		dt=p.rlt-age; // remaining life span
		alive = false; // this root is dead
	}
	age+=dt;

	if (alive) { // dead roots wont grow

		// probabilistic branching model (todo test)
		if ((age>0) && (age-dt<=0)) { // the root emerges in this time step
			double P = rootsystem->getRootTypeParameter(param.type)->sbp->getValue(nodes.back(),this);
			if (P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
				double p = 1.-std::pow((1.-P), dt); //probability of emergence in this time step
				std::cout <<P<<", "<<p<< "\n";
				if (rootsystem->rand()>p) { // not rand()<p
					age -= dt; // the root does not emerge in this time step
				}
			}
		}

		if (age>0) {

			// children first (lateral roots grow even if base root is inactive)
			for (auto l:laterals) {
				l->simulate(dt,silence);
			}

			if (active) {

				// length increment
				double length_ = getLength(std::max(age-dt,0.)); // length of the root for unimpeded growth (i.e. length_==length for unimpeded growth)
				double targetlength = getLength(age);
				double e = targetlength-length_; //elongation in time step dt
				double scale = rootsystem->getRootTypeParameter(param.type)->se->getValue(nodes.back(),this); // hope some of this is optimized out if not set
				double dl = std::max(scale*e, double(0)); // length increment, dt is not used anymore

				// create geometry
				if (p.nob>0) { // root has laterals
					// basal zone
					if ((dl>0)&&(length<p.lb)) { // length is the current length of the root
						if (length+dl<=p.lb) {
							createSegments(dl,silence);
							length+=dl;
							dl=0;
						} else {
							double ddx = p.lb-length;
							createSegments(ddx,silence);
							dl-=ddx; // ddx already has been created
							length=p.lb;
						}
					}
					// branching zone
					if ((dl>0)&&(length>=p.lb)) {
						double s = p.lb; // summed length
						for (size_t i=0; ((i<p.ln.size()) && (dl>0)); i++) {
							s+=p.ln.at(i);
							if (length<s) {
								if (i==laterals.size()) { // new lateral
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
							if (p.ln.size()==laterals.size()) { // new lateral (the last one)
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
			active = getLength(std::max(age,0.))<(p.getK()-dx()/10); // become inactive, if final length is nearly reached
		}
	} // if alive
}

/**
 * Analytical creation (=emergence) time of a node at a length along the root
 *
 * @param length   length of the root [cm]
 */
double Root::getCreationTime(double length)
{
	assert(length>=0);
	double rootage = getAge(length);
	if (rootage<0) {
		std::cout << "Root::getCreationTime() negative root age "<<rootage<<" at length "<< length;
		std::cout.flush();
		throw std::invalid_argument( "bugbugbug" );
	}
	if (parent!=nullptr) {
		double pl = parent_base_length+parent->param.la; // parent length, when this root was created
		double page=parent->getCreationTime(pl);
		assert(page>=0);
		return rootage+page;
	} else {
		return rootage;
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
	return rootsystem->gf.at(param.type-1)->getLength(age,param.r,param.getK(),this);
}

/**
 * Analytical age of the root at a given length
 *
 * @param length   length of the root [cm]
 */
double Root::getAge(double length)
{
	assert(length>=0);
	return rootsystem->gf.at(param.type-1)->getAge(length,param.r,param.getK(),this);
}

RootTypeParameter* Root::getRootTypeParameter() const
{
	return rootsystem->getRootTypeParameter(param.type);
}

/**
 * Creates a new lateral by calling RootSystem::createNewRoot().
 *
 * Overwrite this method to implement more sezialized root classes.
 */
void Root::createLateral(bool silence)
{
	// std::cout << "createLateral()\n";
	const RootParameter &p = param; // rename

	int lt = rootsystem->getRootTypeParameter(p.type)->getLateralType(nodes.back());
	//std::cout << "lateral type " << lt << "\n";

	if (lt>0) {

		Vector3d h; // old heading
		if (nodes.size()>1) {
			h = nodes.back().minus(nodes.at(nodes.size()-2)); // getHeading(b-a)
			// std::cout << "Heading " << h.toString() << "\n";
		} else {
			h= iheading;
		}

		double ageLN = this->getAge(length); // age of root when lateral node is created
		double ageLG = this->getAge(length+p.la); // age of the root, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait

		Organ* lateral = rootsystem->createRoot(lt,  h, delay,  this, length, nodes.size()-1);
		laterals.push_back(lateral);
		lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
		//cout << "time overhead " << age-ageLN << "\n";
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
	int nn = nodes.size();
	if (old_non==0) { // first call of createSegments (in Root::simulate)
		if (nn>1) {
			auto n2 = nodes.at(nn-2);
			auto n1 = nodes.at(nn-1);
			double olddx = n1.minus(n2).length();
			if (olddx<dx()*0.99) { // shift node instead of creating a new node

				Vector3d h; // current heading
				if (nn>2) {
					h = n2.minus(nodes.at(nn-3));
					h.normalize();
				} else {
					h = iheading;
				}
				double sdx = std::min(dx()-olddx,l);

				Matrix3d ons = Matrix3d::ons(h);
				Vector2d ab = rootsystem->tf.at(param.type-1)->getHeading(nodes.back(),ons,olddx+sdx,this);
				ons.times(Matrix3d::rotX(ab.y));
				ons.times(Matrix3d::rotZ(ab.x));
				Vector3d newdx = Vector3d(ons.column(0).times(sdx));

				Vector3d newnode = Vector3d(nodes.back().plus(newdx));
				sl = sdx;
				double et = this->getCreationTime(length+sl);
				nodes[nn-1] = newnode;
				netimes[nn-1] = std::max(et,rootsystem->getSimTime()); // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
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
			std::cout << "skipped small segment (<"<< smallDx << ") \n";
		}
		return;
	}

	int n = floor(l/dx());
	// create n+1 new nodes
	for (int i=0; i<n+1; i++) {

		Vector3d h; // current heading
		if (nodes.size()>1) {
			h = nodes.back().minus(nodes.at(nodes.size()-2));
			h.normalize();
		} else {
			h = iheading;
		}

		double sdx; // segment length (<=dx)
		if (i<n) {  // normal case
			sdx = dx();
		} else { // last segment
			sdx = l-n*dx();
			if (sdx<smallDx) {
				if (!silence) {
					std::cout << "skipped small segment (<"<< smallDx << ") \n";
				}
				return;
			}
		}
		sl+=sdx;

		Matrix3d ons = Matrix3d::ons(h);
		Vector2d ab = rootsystem->tf.at(param.type-1)->getHeading(nodes.back(),ons,sdx,this);
		ons.times(Matrix3d::rotX(ab.y));
		ons.times(Matrix3d::rotZ(ab.x));
		Vector3d newdx = Vector3d(ons.column(0).times(sdx));
		Vector3d newnode = Vector3d(nodes.back().plus(newdx));
		double et = this->getCreationTime(length+sl);
		et = std::max(et,rootsystem->getSimTime()); // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
		addNode(newnode,et);

	} // for

}

/**
 * Returns the root system as sequential list,
 * copies only roots with more than 1 node.
 *
 * \return sequential list of roots
 */
std::vector<Organ*> Root::getRoots()
{
	std::vector<Organ*> v = std::vector<Organ*>();
	getRoots(v);
	return v;
}

/**
 * Returns the root system as sequential list,
 * copies only roots with more than 1 node.
 *
 * @param v     adds the subrootsystem to this vector
 */
void Root::getRoots(std::vector<Organ*>& v)
{
	if (this->nodes.size()>1) {
		v.push_back(this);
	}
	for (auto const& l:this->laterals) {
		l->getRoots(v);
	}
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
	nodes.push_back(n); // node
	nodeIds.push_back(rootsystem->getNodeIndex()); // new unique id
	netimes.push_back(t); // exact creation time
}

/**
 * writes RSML root tag
 *
 * @param cout      typically a file out stream
 * @param indent    we care for looks
 */
void Root::writeRSML(std::ostream & cout, std::string indent) const
{
	if (this->nodes.size()>1) {
		cout << indent << "<root id=\"" <<  id << "\">\n";  // open root

		/* geometry tag */
		cout << indent << "\t<geometry>\n"; // open geometry
		cout << indent << "\t\t<polyline>\n"; // open polyline
		// polyline nodes
		cout << indent << "\t\t\t" << "<point ";
		Vector3d v = nodes.at(0);
		cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
		int n = this->rootsystem->rsmlReduction;
		for (size_t i = 1; i<nodes.size()-1; i+=n) {
			cout << indent << "\t\t\t" << "<point ";
			Vector3d v = nodes.at(i);
			cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
		}
		cout << indent << "\t\t\t" << "<point ";
		v = nodes.at(nodes.size()-1);
		cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
		cout << indent << "\t\t</polyline>\n"; // close polyline
		cout << indent << "\t</geometry>\n"; // close geometry

		/* properties */
		cout << indent <<"\t<properties>\n"; // open properties
		// TODO
		cout << indent << "\t</properties>\n"; // close properties

		cout << indent << "\t<functions>\n"; // open functions
		cout << indent << "\t\t<function name='emergence_time' domain='polyline'>\n"; // open functions
		cout << indent << "\t\t\t" << "<sample>" << netimes.at(0) << "</sample>\n";
		for (size_t i = 1; i<netimes.size()-1; i+=n) {
			cout << indent << "\t\t\t" << "<sample>" << netimes.at(i) << "</sample>\n";

		}
		cout << indent << "\t\t\t" << "<sample>" << netimes.at(netimes.size()-1) << "</sample>\n";

		cout << indent << "\t\t</function>\n"; // close functions
		cout << indent << "\t</functions>\n"; // close functions

		/* laterals roots */
		for (size_t i = 0; i<laterals.size(); i++) {
			laterals[i]->writeRSML(cout,indent+"\t");
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
	str << "Root #"<< id <<": type "<<param.type << ", length: "<< length << ", age: " <<age<<" with "<< laterals.size() << " laterals\n";
	return str.str();
}


