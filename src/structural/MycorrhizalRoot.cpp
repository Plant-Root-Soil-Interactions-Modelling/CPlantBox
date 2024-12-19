#include "MycorrhizalRoot.h"
#include "Root.h"
#include "Organ.h"

namespace CPlantBox {

MycorrhizalRoot::MycorrhizalRoot(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
    Vector3d partialIHeading_, int pni, bool moved, int oldNON)
     :Root(id, param, alive, active, age, length,
	 partialIHeading_,pni, moved,  oldNON ) {}

MycorrhizalRoot::MycorrhizalRoot(std::shared_ptr<Organism> rs, int type,  double delay, std::shared_ptr<Organ> parent, int pni)
:Root(rs,type, delay,parent, pni) {}

void MycorrhizalRoot::addNode(Vector3d n, int id, double t, size_t index, bool shift) {
    Organ::addNode(n, id,  t,  index, shift);
    infected.push_back(0);
}



void MycorrhizalRoot::simulate(double dt, bool verbose)
{
    Root::simulate(dt,verbose);
    // branching points updaten
    // TODO how to update branching points?
    // indices entlang polyline
    // TODO find out what is meant by this

    // TODO how to iterate over all nodes? just get Number of Nodes? (minus the first one)

    // TODO separate primary infection? from secondary? - yes have to
    
    int n = getNumberOfNodes();

    for (size_t i = 1 ; i < n; i++)
    {
        if (rand()< getRootRandomParameter()->p && infected.at(i-1)== 0)
        {
            infected.at(i-1) = 1;
        }
        
    }
    
    // über alle knoten iterieren (-1)
    // "spontane" infektion mit wahrscheinlihckeit p
    // 
    // TODO "länge" der infektion

    // TODO how to check if neighbors are infected
    // immer von der urpsrünglichen infektion "länge" der infektion
    // schauen ob beide nachbarn infiziert sind dann nix
    // falls nur einer 
    // ausrechnen wie viele infiziert werden
    // dann diese anzahl als infiziert setzen

    // in beide richtungen wachsen lassen bei branching points mit gleicher geschwindigkeit (erstmal annehmen)

    // bleibt in der eigenen polyline, ruft dann laterals als children auf

}
    
std::shared_ptr<const MycorrhizalRootSpecificParameter> MycorrhizalRoot::param() const
{
    return std::static_pointer_cast<const MycorrhizalRootSpecificParameter>(param_);
}

std::shared_ptr<MycorrhizalRootRandomParameter> MycorrhizalRoot::getRootRandomParameter() const
{
    return std::static_pointer_cast<MycorrhizalRootRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_root, param_->subType));
}

double MycorrhizalRoot::getParameter(std::string name) const {
    if (name == "infected") {return param() -> infected;}
    return Root::getParameter(name);
}

}