#ifndef MYCORRHIZALROOT_H_
#define MYCORRHIZALROOT_H_

#include "Mycorrhizalrootparameter.h"
#include "Root.h"
#include "Organ.h"
#include "Organism.h"


namespace CPlantBox {
    	
    class MycorrhizalRoot :public Root {
        MycorrhizalRoot(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
        Vector3d partialIHeading_, int pni, bool moved= false, int oldNON = 0);

        MycorrhizalRoot(std::shared_ptr<Organism> rs, int type, double delay, std::shared_ptr<Organ> parent, int pni);

        virtual ~MycorrhizalRoot() { };

        // copy
        // organType ?
        // simualte
        //get PArameter
        //toString
        //getMycRootRandomParameter
        //Mycparam


    };

}


#endif