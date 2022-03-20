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
    Vector3d iHeading, int pni, bool moved, int oldNON)
 :Organ(id, param, alive, active, age, length, Matrix3d::ons(iHeading), pni, moved,  oldNON )
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
Root::Root(std::shared_ptr<Organism> rs, int type, Vector3d iHeading_, double delay, std::shared_ptr<Organ> parent, int pni)
    :Organ(rs, parent, Organism::ot_root, type, delay, Matrix3d::ons(iHeading_), pni) // <- OrganRandomParameter::realize() is called here
{
    assert(parent!=nullptr && "Root::Root parent must be set");
    double beta = 2*M_PI*plant.lock()->rand(); // initial rotation
    double theta = param()->theta;
    if (parent->organType()!=Organism::ot_seed) { // scale if not a baseRoot
        double scale = getRootRandomParameter()->f_sa->getValue(parent->getNode(pni), parent);
        theta*=scale;
    }
    insertionAngle = theta;
    Vector3d newHeading = iHeading.times(Vector3d::rotAB(theta,beta));
    iHeading = Matrix3d::ons(newHeading); // new initial heading
    if (parent->organType()!=Organism::ot_seed) { // the first node of the base roots must be created in RootSystem::initialize()
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
                            if (length<s) {
                                if (i==children.size()) { // new lateral
                                    createLateral(dt_, verbose);
                                }
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
                        if (p.ln.size()==children.size()&& (getLength(true)>=s)){
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
        double ageLN = this->calcAge(getLength(true)); // age of root when lateral node is created
        ageLN = std::max(ageLN, age-dt);
        double meanLn = getRootRandomParameter()->ln; // mean inter-lateral distance
        double effectiveLa = std::max(param()->la-meanLn/2, 0.); // effective apical distance, observed apical distance is in [la-ln/2, la+ln/2]
        double ageLG = this->calcAge(getLength(true)+effectiveLa); // age of the root, when the lateral starts growing (i.e when the apical zone is developed)
        double delay = ageLG-ageLN; // time the lateral has to wait
        auto lateral = std::make_shared<Root>(plant.lock(), lt,  heading(), delay,  shared_from_this(), nodes.size()-1);
        children.push_back(lateral);
        lateral->simulate(age-ageLN,verbose); // pass time overhead (age we want to achieve minus current age)
    }
}






/**
 * @return Quick info about the object for debugging
 * additionally, use getParam()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Root::toString() const
{
    std::stringstream newstring;
    newstring << "; initial heading: " << iHeading.toString() << ", parent node index" << parentNI << ".";
    return  Organ::toString()+newstring.str();
}

} // end namespace CPlantBox
