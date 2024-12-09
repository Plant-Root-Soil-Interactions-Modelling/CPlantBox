#include "MycorrhizalRoot.h"
#include "Root.h"
#include "Organ.h"

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
MycorrhizalRoot::MycorrhizalRoot(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
    Vector3d partialIHeading_, int pni, bool moved, int oldNON)
     :Root(id, param, alive, active, age, length,
	 partialIHeading_,pni, moved,  oldNON )
      {}
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
MycorrhizalRoot::MycorrhizalRoot(std::shared_ptr<Organism> rs, int type,  double delay, std::shared_ptr<Organ> parent, int pni)
:Root(rs,type,delay, parent,pni) // <- OrganRandomParameter::realize() is called here
{}

std::shared_ptr<Organ> MycorrhizalRoot::copy(std::shared_ptr<Organism> rs) {
    auto r = std::make_shared<MycorrhizalRoot>(*this);
    r->parent = std::weak_ptr<Organ>();
    r->plant = rs;
    r->param_ = std::make_shared<MycorrhizalRootSpecificParameter>(*param());
    for (size_t i=0; i< children.size(); i++) {
    r->children[i] = children[i]->copy(rs); // copy laterals
    r->children[i]->setParent(r);
    }
    return r;
}
    //simulate
    // organType ?
    // simualte
    // double MycorrhizalRoot::getParameter(std::string name) const {
    //     //if (name=="vi") {return param()->vi;}
    //     return Root::getParameter(name);
    // }
    //toString
    //getMycRootRandomParameter
    //Mycparam
}