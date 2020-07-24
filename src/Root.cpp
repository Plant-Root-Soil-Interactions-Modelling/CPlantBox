// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Root.h"

namespace CPlantBox {

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
Root::Root(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
    Vector3d iheading, double pbl, int pni, bool moved, int oldNON)
 :Organ(id, param, alive, active, age, length, iheading, pbl, pni, moved,  oldNON )
{
    insertionAngle = this->param()->theta;
}


/**
 * Constructor: Should be only called during simulation by Root::createLateral().
 * For base roots the initial node and node creation time must be set from outside
 *
 * @param rs 			points to RootSystem
 * @param type 		    type of root that is created
 * @param heading		heading of parent root at emergence
 * @param delay 		to give apical zone of parent time to develop
 * @param parent		parent root
 * @param pbl			parent base length
 * @param pni			parent node index
 */
Root::Root(std::shared_ptr<Organism> rs, int type, Vector3d heading, double delay, std::shared_ptr<Organ> parent, double pbl, int pni)
    :Organ(rs, parent, Organism::ot_root, type, delay, heading, pbl, pni) // <- OrganRandomParameter::realize() is called here
{
    assert(parent!=nullptr && "Root::Root parent must be set");
    double beta = 2*M_PI*plant.lock()->rand(); // initial rotation
    double theta = param()->theta;
    if (parent->organType()!=Organism::ot_seed) { // scale if not a baseRoot
        double scale = getRootRandomParameter()->f_sa->getValue(parent->getNode(pni), parent);
        theta*=scale;
    }
    insertionAngle = theta;
    iHeading = Matrix3d::ons(heading).times(Vector3d::rotAB(theta,beta)); // new initial heading
    if (parent->organType()!=Organism::ot_seed) { // the first node of the base roots must be created in RootSystem::initialize()
        // assert(pni+1 == parent->getNumberOfNodes() && "Root::Root: at object creation always at last node");
        addNode(parent->getNode(pni), parent->getNodeId(pni), parent->getNodeCT(pni)+delay);
    }
}

/**
 * Deep copies the organ into the new plant @param rs.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
std::shared_ptr<Organ> Root::copy(std::shared_ptr<Organism> rs)
{
    auto r = std::make_shared<Root>(*this); // shallow copy
    r->parent = std::weak_ptr<Organ>();
    r->plant = rs;
    r->param_ = std::make_shared<RootSpecificParameter>(*param()); // copy parameters
    for (size_t i=0; i< children.size(); i++) {
        r->children[i] = children[i]->copy(rs); // copy laterals
        r->children[i]->setParent(r);
    }
    return r;
}

/**
 * Simulates the development of the organ in a time span of @param dt days.
 *
 * @param dt        time step [day]
 * @param verbose   turns console output on or off
 */
void Root::simulate(double dt, bool verbose)
{
    firstCall = true;
    moved = false;
    oldNumberOfNodes = nodes.size();

    const RootSpecificParameter& p = *param(); // rename

    if (alive) { // dead roots wont grow

        // increase age
        if (age+dt>p.rlt) { // root life time
            dt=p.rlt-age; // remaining life span
            alive = false; // this root is dead
        }
        age+=dt;

        // probabilistic branching model
        if ((age>0) && (age-dt<=0)) { // the root emerges in this time step
            double P = getRootRandomParameter()->f_sbp->getValue(nodes.back(),shared_from_this());
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
                double scale = getRootRandomParameter()->f_se->getValue(nodes.back(), shared_from_this());
                double dl = std::max(scale*e, 0.); // length increment
                // std::cout << "Root::simulate: " << dt_ << ", " << e <<  ", " << scale <<"\n"; // for debugging

                // create geometry
                if (p.ln.size()>0) { // root has children
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
                    if ((dl>0)&&(length>=p.lb)) {
                        double s = p.lb; // summed length
                        for (size_t i=0; ((i<p.ln.size()) && (dl>0)); i++) {
                            s+=p.ln.at(i);
                            if (length<s) {
                                if (i==children.size()) { // new lateral
                                    createLateral(dt_, verbose);
                                }
                                if (length+dl<=s) { // finish within inter-lateral distance i
                                    createSegments(dl,dt_,verbose);
                                    length+=dl;
                                    dl=0;
                                } else { // grow over inter-lateral distance i
                                    double ddx = s-length;
                                    createSegments(ddx,dt_,verbose);
                                    dl-=ddx;
                                    length=s;
                                }
                            }
                        }
                        if (p.ln.size()==children.size()) { // new lateral (the last one)
                            createLateral(dt_, verbose);
                        }
                    }
                    /* apical zone */
                    if (dl>0) {
                        createSegments(dl,dt_,verbose);
                        length+=dl;
                    }
                } else { // no laterals
                    if (dl>0) {
                        createSegments(dl,dt_,verbose);
                        length+=dl;
                    }
                } // if lateralgetLengths
            } // if active
            active = length<(p.getK()-dx()/10); // become inactive, if final length is nearly reached
        }
    } // if alive

}

/**
 * Analytical creation (=emergence) time of a point along the already grown root
 *
 * @param length   length along the root, where the point is located [cm]
 * @param dt 	   current time step [day]
 * @return         the analytic time when this point was reached by the growing root [day],
 * 				   if growth is impeded, the value is not exact, but approximated dependent on the temporal resolution.
 */
double Root::calcCreationTime(double length, double dt)
{
    assert(length >= 0 && "Root::getCreationTime() negative length");
    double age_ = calcAge(length); // root age as if grown unimpeded (lower than real age)
    double a = std::max(age_, age-dt);
    a = std::min(a, age); // a in [age-dt, age]
    return a+nodeCTs[0];
}

/**
 * Analytical length of the single root at a given age
 *
 * @param age          age of the root [day]
 * @return             root length [cm]
 */
double Root::calcLength(double age)
{
    assert(age >= 0 && "Root::calcLength() negative root age");
    return getRootRandomParameter()->f_gf->getLength(age,param()->r,param()->getK(), shared_from_this());
}

/**
 * Analytical age of the single root at a given length
 *
 * @param length   length of the root [cm]
 * @return local age [day]
 */
double Root::calcAge(double length)
{
    assert(length >= 0 && "Root::calcAge() negative root length");
    return getRootRandomParameter()->f_gf->getAge(length,param()->r,param()->getK(), shared_from_this());
}

/**
 * @return The RootTypeParameter from the plant
 */
std::shared_ptr<RootRandomParameter> Root::getRootRandomParameter() const
{
    return std::static_pointer_cast<RootRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_root, param_->subType));
}

/**
 * @return Parameters of the specific root
 */
 std::shared_ptr<const RootSpecificParameter> Root::param() const
{
    return std::static_pointer_cast<const RootSpecificParameter>(param_);
}

/**
 * Creates a new lateral root and passes time overhead
 *
 * Overwrite this method to implement more spezialized root classes.
 *
 * @param verbose   turns console output on or off
 */
void Root::createLateral(double dt, bool verbose)
{
    int lt = getRootRandomParameter()->getLateralType(nodes.back());
    if (lt>0) {
        double ageLN = this->calcAge(length); // age of root when lateral node is created
        ageLN = std::max(ageLN, age-dt);
        double meanLn = getRootRandomParameter()->ln; // mean inter-lateral distance
        double effectiveLa = std::max(param()->la-meanLn/2, 0.); // effective apical distance, observed apical distance is in [la-ln/2, la+ln/2]
        double ageLG = this->calcAge(length+effectiveLa); // age of the root, when the lateral starts growing (i.e when the apical zone is developed)
        double delay = ageLG-ageLN; // time the lateral has to wait
        auto lateral = std::make_shared<Root>(plant.lock(), lt,  heading(), delay,  shared_from_this(), length, nodes.size()-1);
        children.push_back(lateral);
        lateral->simulate(age-ageLN,verbose); // pass time overhead (age we want to achieve minus current age)
    }
}

/**
 * @return Current root heading
 */
Vector3d Root::heading() const
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
Vector3d Root::getIncrement(const Vector3d& p, double sdx)
{
    Vector3d h = heading();
    Matrix3d ons = Matrix3d::ons(h);
    Vector2d ab = getRootRandomParameter()->f_tf->getHeading(p, ons, dx(), shared_from_this());
    Vector3d sv = ons.times(Vector3d::rotAB(ab.x,ab.y));
    return sv.times(sdx);
}

/**
 *  Creates nodes and node emergence times for a length l
 *
 *  Checks that each new segments length is <= dx but >= parent->minDx
 *
 *  @param l        total length of the segments that are created [cm]
 *  @param dt       time step [day]
 *  @param verbose  turns console output on or off
 */
void Root::createSegments(double l, double dt, bool verbose)
{
    if (l==0) {
        std::cout << "Root::createSegments: zero length encountered \n";
        return;
    }
    if (l<0) {
        std::cout << "Root::createSegments: negative length encountered \n";
    }

    // shift first node to axial resolution
    double shiftl = 0; // length produced by shift
    int nn = nodes.size();
    if (firstCall) { // first call of createSegments (in Root::simulate)
        firstCall = false;
        if ((nn>1) && (children.empty() || (nn-1 != std::static_pointer_cast<Root>(children.back())->parentNI)) ) { // don't move a child base node
            Vector3d n2 = nodes[nn-2];
            Vector3d n1 = nodes[nn-1];
            Vector3d h = n1.minus(n2);
            double olddx = h.length(); // length of last segment
            if (olddx<dx()*0.99) { // shift node instead of creating a new node
                shiftl = std::min(dx()-olddx, l);
                double sdx = olddx + shiftl; // length of new segment
                // Vector3d newdxv = getIncrement(n2, sdx);
                h.normalize();
                nodes[nn-1] = Vector3d(n2.plus(h.times(sdx))); // n2.plus(newdxv)
                double et = this->calcCreationTime(length+shiftl, dt);
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
            if (sdx<plant.lock()->getMinDx()) { // quit if l is too small
                if (verbose) {
                    std::cout << "skipped small segment ("<< sdx <<" < "<< plant.lock()->getMinDx() << ") \n";
                }
                return;
            }
        }
        sl += sdx;
        Vector3d newdx = getIncrement(nodes.back(), sdx);
        Vector3d newnode = Vector3d(nodes.back().plus(newdx));
        double et = this->calcCreationTime(length+shiftl+sl, dt);
        // in case of impeded growth the node emergence time is not exact anymore,
        // but might break down to temporal resolution
        addNode(newnode, et);
    }
}

/**
 * @copydoc Organ::getParameter
 *
 * Note:
 * lnMean, and lnDev denotes the mean and standard deviation of the inter-lateral distance of this organ
 * ln_mean, and ln_dev is the mean and standard deviation from the RootRandomParmaeters
 */
double Root::getParameter(std::string name) const
{
    // specific parameters
    if (name=="type") { return this->param_->subType; }  // in CPlantBox the subType is often called just type
    if (name=="lb") { return param()->lb; } // basal zone [cm]
    if (name=="la") { return param()->la; } // apical zone [cm]
    if (name=="r"){ return param()->r; }  // initial growth rate [cm day-1]
    if (name=="theta") { return insertionAngle; } // angle between root and parent root [rad]
    if (name=="rlt") { return param()->rlt; } // root life time [day]
    // specific parameters member functions
    if (name=="nob") { return param()->nob(); } // number of lateral emergence nodes
    if (name=="k") { return param()->getK(); }; // maximal root length [cm]
    if (name=="lmax") { return param()->getK(); }; // maximal root length [cm]
    // member functions
    if (name=="numberOfLaterals") { return getNumberOfLaterals(); }
    // further
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
    if (name=="rootLength") { return getLength(); } // root length [cm], same as length, but a SegmentAnalyser::getParameter call would give the segment length
    if (name=="volume") { return param()->a*param()->a*M_PI*getLength(); } // root volume [cm^3]
    if (name=="surface") { return 2*param()->a*M_PI*getLength(); } // root surface [cm^2]
    return Organ::getParameter(name); // pass to base class
}

/**
 * @return The number of emerged lateral roots (i.e. number of children with age>0)
 * @see Organ::getNumberOfChildren
 */
int Root::getNumberOfLaterals() const {
	int nol = 0;
	for (auto& c : children)  {
		if (c->getAge()>0) { // born
			nol ++;
		}
	}
	return nol;
}

/**
 * @return Quick info about the object for debugging
 * additionally, use getParam()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Root::toString() const
{
    std::string str = Organ::toString();
    str.replace(0, 5, "Root");
    std::stringstream newstring;
    newstring << "; initial heading: " << iHeading.toString() << ", parent base length " << parentBaseLength << ", parent node index" << parentNI << ".";
    return str+newstring.str();
}

} // end namespace CPlantBox
