#include "MycorrhizalRoot.h"
#include "Root.h"
#include "Organ.h"
#include "Organism.h"

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
    infectionTime.push_back(NAN);
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
    

    // Primary Infection
    int m = getNumberOfSegments();
    auto seg = getSegments();	
    // "spontane" infektion mit wahrscheinlihckeit p
    // über alle knoten iterieren (-1)
    for (size_t i = 0; i < m; i++)
    {
        auto seglength = getLength(seg[i].y) - getLength(seg[i].x);
        if (plant.lock()->rand() < (getRootRandomParameter()->p)*dt*seglength && infected.at(seg[i].y) == 0)
        {
            infected.at(seg[i].y) = 1;
            infectionTime.at(seg[i].y) = age + dt;// TODO check if + dt necessary or not
        }
    }
    
    // Secondary Infection
    auto max_length_infection = dt*getRootRandomParameter()->vi;
    for (size_t i = 0; i < m; i++)
    {
        if (infected.at(seg[i].y) == 1)
        { // if length at x is more than length at y + max_length_infection
            auto seg_length = getLength(seg[i].y) - getLength(seg[i].x);
            if (infected.at(seg[i].x) == 0 && dt*seg_length < max_length_infection)
            {
                infected.at(seg[i].x) = 2;
            }
            
        }
    }
    
    
 
   

    // // length always measured from start so max distance "infection" can travel is
    // // length of first node - max_length_infection
    // for (size_t i = 0; i < segments.size() ; i++)
    // {
    //     if (getNodeInfection(segments[i].y) == 1)
    //     {
    //         int length = getLength(segments[i].y);
    //         int min_length = length - max_length_infection;
    //         int max_length = length + max_length_infection;
    //         if (getNodeInfection(segments[i].x) == 0)
    //         {
    //             if (getLength(segments[i].x) < max_length && getLength(segments[i].x) > min_length)
    //             {
    //                 infected.at(segments[i].x) = 2;
    //             } 
    //         }
            
    //     }
    // }   
    double infection_length = dt*getRootRandomParameter()->vi;

    auto segs = getSegments();
        std::vector<int> segId = std::vector<int>(segs.size());
        for (int i=0; i<segs.size(); i++) {
            segId[i] = segs[i].y-1;
        }

        for (int i = 0; i < segs.size(); i++)
        {
            if (infected.at(segs[i].y) == 1)
            {
                int length_origin = getLength(segs[i].y);
                int newlength = getLength(segs[i].x);
                while (newlength  > length_origin +  infection_length) 
                    // TODO darf nur =0
                    infected.at(segs[i].x) = 2;
                    int j;
                    for (int k=0; k<segs.size(); k++) {
                        if (segs[i].x == segId[k] + 1) {
                            j = k;
                        }
                    }
                    newlength = getLength(segs[j].x);
                }
            }
        
    
    // IDEE: durch segmente durchiterieren und dann bei einem infizierten knoten die nachbarn anschauen
    // wenn beide infiziert sind dann nix
    // wenn nur einer infiziert ist dann ausrechnen schauen wie weit die infektion wandert oder ob sie durch eine andere infizierte branch gestoppt wird

    
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