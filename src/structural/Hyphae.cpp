// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Hyphae.h"
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
Hyphae::Hyphae(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
    Vector3d partialIHeading_, int pni, bool moved, int oldNON)
                     :Organ(id, param, alive, active, age, length, partialIHeading_,pni, moved,  oldNON )
                      {
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
Hyphae::Hyphae(std::shared_ptr<Organism> rs, int type,  double delay, std::shared_ptr<Organ> parent, int pni)
:Organ(rs, parent, Organism::ot_hyphae, type, delay,  pni) // <- OrganRandomParameter::realize() is called here
{
    // std::cout << "create Hyphae\n" << std::flush;
    assert(parent!=nullptr && "Hyphae::Hyphae parent must be set");
    double beta = 2*M_PI*plant.lock()->rand(); // initial rotation
    double theta = M_PI/2.;
    this->partialIHeading = Vector3d::rotAB(theta,beta);
    double creationTime= parent->getNodeCT(pni)+delay;//default
    addNode(parent->getNode(pni), parent->getNodeId(pni), creationTime);
}

/**
 * Deep copies the organ into the new plant @param rs.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
std::shared_ptr<Organ> Hyphae::copy(std::shared_ptr<Organism> rs)
{
    auto r = std::make_shared<Hyphae>(*this); // shallow copy
    r->parent = std::weak_ptr<Organ>();
    r->plant = rs;
    r->param_ = std::make_shared<HyphaeSpecificParameter>(*param()); // copy parameters
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
void Hyphae::simulate(double dt, bool verbose)
{

//    firstCall = true;
//    moved = false;
    oldNumberOfNodes = nodes.size();

    const HyphaeSpecificParameter& p = *param(); // rename

    if (alive) { // dead hypaes won't grow

        // increase age
        if (age+dt>p.hlt) { // root life time
            dt=p.hlt-age; // remaining life span
            alive = false; // this root is dead
        }
        age+=dt;

        if (age>0) { // unborn  roots have no children

            if (children.size() == 0) { // ELONGATE

                if (active) {

                    double targetlength = p.v*(age+dt);
                    double dl = targetlength-length; // unimpeded elongation in time step dt
                    // double scale = 1.; //getHyphaeRandomParameter()->f_se->getValue(nodes.back(), shared_from_this());
                    // double dl = std::max(scale*e, 0.);//  length increment = calculated length + increment from last time step too small to be added
                    createSegments(dl,dt,verbose);
                    // std::cout << "*";
                    length+=dl;
                }

                // active = getLength(false)<=(p.getK()*(1 - 1e-11)); // become inactive, if final length is nearly reached

            } else { // NOT ACTIVE (children grow)

                // children first (lateral roots grow even if base root is inactive)
                for (auto l:children) {
                    l->simulate(dt,verbose);
                }
            }

        } // age>0
    } // if alive

}

///**
// * Analytical length of the single root at a given age
// *
// * @param age          age of the root [day]
// * @return             root length [cm]
// */
double Hyphae::calcLength(double age)
{
   assert(age >= 0 && "Hyphae::calcLength() negative hyphae age");
   return getHyphaeRandomParameter()->f_gf->getLength(age,param()->v,0., shared_from_this());
}

///**
// * Analytical age of the single root at a given length
// *
// * @param length   length of the root [cm]
// * @return local age [day]
// */

double Hyphae::calcAge(double length) const
{
    assert(length >= 0 && "Hyphae::calcAge() negative hyphae length");
    return getHyphaeRandomParameter()->f_gf->getAge(length,param()->v, 0., shared_from_this());
}

/**
 * @return The RootTypeParameter from the plant
 */
std::shared_ptr<HyphaeRandomParameter> Hyphae::getHyphaeRandomParameter() const
{
    return std::static_pointer_cast<HyphaeRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_hyphae, param_->subType));
}

/**
 * @return Parameters of the specific root
 */
std::shared_ptr<const HyphaeSpecificParameter> Hyphae::param() const
{
    return std::static_pointer_cast<const HyphaeSpecificParameter>(param_);
}

/**
 * @copydoc Organ::getParameter
 * @param name  name of the parameter
 */
double Hyphae::getParameter(std::string name) const 
{
    // specific parameters
        if (name=="type") { return this->param_->subType; }  // delete to avoid confusion?
        if (name=="subType") { return this->param_->subType; }  // organ sub-type [-]
        if (name=="v") { return param()->v; } // Tip elongation rate [cm day-1]
        if (name=="b") { return param()->b; } // Branching rate [1 day-1]
        if (name=="hlt") { return param()->hlt; } // Hyphal life time [day]
        if (name=="theta") { return param()->theta; } // Branching angle [rad]
    return Organ::getParameter(name); // pass to base class
}


/**
 * @return Quick info about the object for debugging
 * additionally, use param()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Hyphae::toString() const
{
    std::stringstream newstring;
    newstring << "."; // TODO
    return  Organ::toString()+newstring.str();
}

} // end namespace CPlantBox
