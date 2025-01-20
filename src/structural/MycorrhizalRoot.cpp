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
    infectionTime.push_back(-1);
}



void MycorrhizalRoot::simulate(double dt, bool verbose)
{
    Root::simulate(dt,verbose);
    // branching points updaten
    // TODO how to update branching points?
    // indices entlang polyline
    // TODO find out what is meant by this


    if (this->nodes.size()>1) {
		//Primary Infection
		for (size_t i=1; i<nodes.size(); i++) {
			
            // double cursegLength = nodes.at(i).length() - nodes.at(i-1).length();
            if (plant.lock()->rand() < (getRootRandomParameter()->p)*dt && infected.at(i-1) == 0)
            {
                infected.at(i-1) = 1;
                infectionTime.at(i-1) = age + dt;
                
		    }
	    }

        //Secondary Infection
        auto max_length_infection = dt*getRootRandomParameter()->vi;
        for (size_t i = 1; i < nodes.size(); i++)
        {
            if (infected.at(i-1) == 1)
            {
                auto max_length_basal = nodes.at(i).length() - max_length_infection;
                auto currentbasalnode = i-1;
                while (nodes.at(currentbasalnode).length()> max_length_basal && infected.at(currentbasalnode) == 0)
                {
                    infected.at(currentbasalnode) = 2;
                    infectionTime.at(currentbasalnode) = age + dt; // FIXME make it the proper time and not just the time step
                    currentbasalnode--;
                }

                auto max_length_apical = nodes.at(i).length() + max_length_infection;
                auto currentapicalnode = i-1;
                while (nodes.at(currentapicalnode).length() < max_length_apical && infected.at(currentapicalnode) == 0)
                {
                    infected.at(currentapicalnode) = 2;
                    infectionTime.at(currentapicalnode) = age + dt; // FIXME make it the proper time and not just the time step
                    currentapicalnode++;
                }
            }
            
        }
        

    }
    
    // Secondary Infection
    // auto max_length_infection = dt*getRootRandomParameter()->vi;
    // std::vector<int> segId = std::vector<int>(seg.size());
    //     for (int i=0; i<seg.size(); i++) {
    //         segId[i] = seg[i].y-1;
    //     }
    
    // for (size_t i = 0; i < m; i++)
    // {
    //     if (infected.at(seg[i].y) == 1)
    //     { 
    //         auto max_length_basal = getLength(seg[i].y) - max_length_infection;
    //         auto currentbasalnode = seg[i].x;
    //         while (getLength(currentbasalnode)> max_length_basal && infected.at(currentbasalnode) == 0)
    //         {
    //             infected.at(currentbasalnode) = 2;
    //             infectionTime.at(currentbasalnode) = age + dt;
    //             int newnode;
    //             for (int k=0; k<segId.size(); k++) {
    //                 if (segId[k] == currentbasalnode-1) {
    //                     newnode = k;
    //                 }
    //             }
    //             currentbasalnode = seg[newnode].x;
    //         }

    //         auto max_length_apical = getLength(seg[i].y) + max_length_infection;
    //         auto currentapicalnode = seg[i].y;
    //         bool found = false;
    //         int j = 0;
    //         while (found == false)
    //         {
    //             if (seg[j].x == currentapicalnode)
    //             {
    //                 found = true;
    //                 currentapicalnode = seg[j].y;
    //             }
    //             else
    //             {
    //                 j++;
    //             }
                
    //         }
            
    //         while (getLength(currentapicalnode) < max_length_apical && infected.at(currentapicalnode) == 0)
    //         {
    //             infected.at(currentapicalnode) = 2;
    //             infectionTime.at(currentapicalnode) = age + dt;
    //             int newnode;
    //             for (int k=0; k<seg.size(); k++) {
    //                 if (seg[k].x == currentapicalnode) {
    //                     newnode = k;
    //                 }
    //             }
    //             currentapicalnode = seg[newnode].y;
    //         }
            
    //     }
    // }
        
    
    // IDEE: durch segmente durchiterieren und dann bei einem infizierten knoten die nachbarn anschauen
    // wenn beide infiziert sind dann nix
    // wenn nur einer infiziert ist dann ausrechnen schauen wie weit die infektion wandert oder ob sie durch eine andere infizierte branch gestoppt wird

    
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
    // if (name == "infected") {return param() -> infected;}
    return Root::getParameter(name);
}

}