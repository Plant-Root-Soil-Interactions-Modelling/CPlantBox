// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Root.h"
#include <numeric>

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
    Vector3d partialIHeading_, int pni, bool moved, int oldNON)
     :Organ(id, param, alive, active, age, length, partialIHeading_,pni, moved,  oldNON )
{
    insertionAngle = this->param()->theta;
}


/**
 * Constructor: Should be only called during simulation by Root::createLateral().
 * For base roots the initial node and node creation time must be set from outside
 *
 * @param plant_      	points to RootSystem
 * @param subType 		subType of root that is created
 * @param heading		heading of parent root at emergence
 * @param delay 		to give apical zone of parent time to develop
 * @param parent		parent root
 * @param pbl			parent base length
 * @param pni			parent node index
 */
Root::Root(std::shared_ptr<Organism> plant_, int subType,  double delay, std::shared_ptr<Organ> parent, int pni)
:Organ(plant_, parent, Organism::ot_root, subType, delay,  pni) // <- OrganRandomParameter::realize() is called here
{
    assert(parent!=nullptr && "Root::Root parent must be set");
    double beta = 2*M_PI*plant.lock()->rand(); // initial rotation
    double theta = param()->theta;
    if (parent->organType()!=Organism::ot_seed) { // scale if not a baseRoot
        double scale = getRootRandomParameter()->f_sa->getValue(parent->getNode(pni), parent);
        theta*=scale;
    }
    insertionAngle = theta;
	this->partialIHeading = Vector3d::rotAB(theta,beta);

	if(!(parent->organType()==Organism::ot_seed))
	{
			double creationTime= parent->getNodeCT(pni)+delay;//default
		if (!parent->hasRelCoord())  // the first node of the base roots must be created in RootSystem::initialize()
		{
			addNode(parent->getNode(pni), parent->getNodeId(pni), creationTime);

		}else{
			if ((parent->organType()==Organism::ot_stem)&&(parent->getNumberOfChildren()>0)) {
			//if lateral of stem, initial creation time:
			//time when stem reached end of basal zone (==CT of parent node of first lateral) + delay
			// @see stem::leafGrow
			creationTime = parent->getChild(0)->getParameter("creationTime") + delay;
		}
			addNode(Vector3d(0.,0.,0.), parent->getNodeId(pni), creationTime);
		}
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
                double p = 1.-(1.-P*dt); //probability of emergence in this time step
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

                double targetlength = calcLength(age_+dt_)+ this->epsilonDx;

                double e = targetlength-length; // unimpeded elongation in time step dt
                double scale = getRootRandomParameter()->f_se->getValue(nodes.back(), shared_from_this());
                double dl = std::max(scale*e, 0.);//  length increment = calculated length + increment from last time step too small to be added
                length = getLength();
                this->epsilonDx = 0.; // now it is "spent" on targetlength (no need for -this->epsilonDx in the following)
                // create geometry
                if (p.laterals ) { // root has children
                    /* basal zone */
                    if ((dl>0)&&(length<p.lb)) { // length is the current length of the root
                        if (length+dl<=p.lb) {
                            createSegments(dl,dt_,verbose);
                            length+=dl; // - this->epsilonDx;
                            dl=0;
                        } else {
                            double ddx = p.lb-length;
                            createSegments(ddx,dt_,verbose);
                            dl-=ddx; // ddx already has been created
                            length=p.lb;
                            //							if(this->epsilonDx != 0){//this sould not happen as p.lb was redefined in rootparameter::realize to avoid this
                            //								throw std::runtime_error("Root::simulate: p.lb - length < dxMin");
                            //							} // this could happen, if the tip ends in this section
                        }
                    }
                    /* branching zone */
                    if ((dl>0)&&(length>=p.lb)) {
                        double s = p.lb; // summed length
                        for (size_t i=0; ((i<p.ln.size()) && (dl > 0)); i++) {
                            s+=p.ln.at(i);
                            if (length<=s) {//need "<=" instead of "<" => in some cases ln.at(i) == 0 when adapting ln to dxMin (@see rootrandomparameter::realize())

                                if (i==created_linking_node) { // new lateral
                                    createLateral(dt_, verbose);
                                }

                                if(length < s)//because with former check we have (length<=s)
                                {
                                    if (length+dl<=s) { // finish within inter-lateral distance i
                                        createSegments(dl,dt_,verbose);
                                        length+=dl; //- this->epsilonDx;
                                        dl=0;
                                    } else { // grow over inter-lateral distance i
                                        double ddx = s-length;
                                        createSegments(ddx,dt_,verbose);
                                        dl-=ddx;
                                        length=s;
                                        //									if(this->epsilonDx != 0){//this sould not happen as p.lb was redefined in rootparameter::realize to avoid this
                                        //										throw std::runtime_error( "Root::simulate: p.ln.at(i) - length < dxMin");
                                        //									} // this could happen, if the tip ends in this section
                                    }

                                }
                            }
                        }

                        if ((p.ln.size()==created_linking_node)&& (getLength(true)-s>-1e-9)){
                            createLateral(dt_, verbose);
                        }
                    }
                    /* apical zone */
                    if (dl>0) {
                        createSegments(dl,dt_,verbose);
                        length+=dl; // - this->epsilonDx;
                    }
                } else { // no laterals

                    if (dl>0) {
                        createSegments(dl,dt_,verbose);
                        length+=dl; //- this->epsilonDx;
                    }
                } // if lateralgetLengths
            } // if active
            active = getLength(false)<=(p.getK()*(1 - 1e-11)); // become inactive, if final length is nearly reached
        }
    } // if alive
    // std::cout << "end" << getId() << "\n" << std::flush;
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
double Root::calcAge(double length) const
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
 * @copydoc Organ::getParameter
 *
 * Note:
 * lnMean, and lnDev denotes the mean and standard deviation of the inter-lateral distance of this organ
 * ln_mean, and ln_dev is the mean and standard deviation from the RootRandomParmaeters
 */
double Root::getParameter(std::string name) const
{
    // specific parameters
    if (name=="type") { return this->param_->subType; }  // delete to avoid confusion?
    if (name=="subType") { return this->param_->subType; }  // organ sub-type [-]
    if (name=="lb") { return param()->lb; } // basal zone [cm]
    if (name=="la") { return param()->la; } // apical zone [cm]
    if (name=="r"){ return param()->r; }  // initial growth rate [cm day-1]
    if (name=="theta") { return insertionAngle; } // angle between root and parent root [rad]
    if (name=="theta_deg") { return insertionAngle/M_PI*180.; } // angle between root and parent root [°]

    if (name=="rlt") { return param()->rlt; } // root life time [day]
    // specific parameters member functions
    if (name=="nob") { return param()->nob(); } // number of lateral emergence nodes/branching points
    if (name=="k") { return param()->getK(); }; // maximal root length [cm]
    if (name=="lmax") { return param()->getK(); }; // maximal root length [cm]
    // further
    if (name=="lnMean") { // mean lateral distance [cm]
        auto& v =param()->ln;
        if(v.size()>0){
            return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        } else {
            return 0;
        }
    }
    if (name=="lnDev") { // standard deviation of lateral distance [cm]
        auto& v =param()->ln;
        double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
        return std::sqrt(sq_sum / v.size() - mean * mean);
    }
    if (name=="rootLength") { return getLength(true); } // root length [cm], same as length, but a SegmentAnalyser::getParameter call would give the segment length
    if (name=="volume") { return param()->a*param()->a*M_PI*getLength(true); } // root volume [cm^3]
    if (name=="surface") { return 2*param()->a*M_PI*getLength(true); } // root surface [cm^2]
    return Organ::getParameter(name); // pass to base class
}


/**
 * @return Quick info about the object for debugging
 * additionally, use param()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Root::toString() const
{
    std::stringstream newstring;
    newstring << "; initial heading: " << getiHeading0().toString() << ", parent node index" << parentNI << ".";
    return  Organ::toString()+newstring.str();
}

/**
 * Static root
 */

StaticRoot::StaticRoot(int id, std::shared_ptr<const OrganSpecificParameter> param, double length, int pni)
        :Root(id, param, true, false, 0., length, Vector3d(0.,0.,-1.), pni, false,  0 )
{ }

void StaticRoot::initializeLaterals() {
    assert(lateralNodeIndices.size()==lateralTypes.size() && lateralNodeIndices.size()==lateralDelays.size() && "Root::Root parent must be set");
    // std::cout << lateralNodeIndices.size() << "\n";
    for (int i=0; i<lateralNodeIndices.size(); i++) {
        int lni = lateralNodeIndices.at(i); // rename
        std::shared_ptr<Organ> lateral = std::make_shared<Root>(plant.lock(), lateralTypes.at(i) , lateralDelays.at(i), shared_from_this(), lni);
        // lateral->addNode(getNode(lni), getNodeId(lni), lateralDelays.at(i));
        this->addChild(lateral);
    }
}

} // end namespace CPlantBox
