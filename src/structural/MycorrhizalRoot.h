#ifndef MYCORRHIZALROOT_H_
#define MYCORRHIZALROOT_H_

#include "Mycorrhizalrootparameter.h"
#include "Root.h"
#include "Organ.h"
#include "Organism.h"


namespace CPlantBox {
    	
    class MycorrhizalRoot :public Root {
        public:
        MycorrhizalRoot(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
        Vector3d partialIHeading_, int pni, bool moved= false, int oldNON = 0);

        MycorrhizalRoot(std::shared_ptr<Organism> rs, int type, double delay, std::shared_ptr<Organ> parent, int pni);

        virtual ~MycorrhizalRoot() { };

        //std::shared_ptr<Organ> copy(std::shared_ptr<Organism> rs) override;
        // organType ?
        // simualte
        //double getParameter(std::string name) const override;
        //toString
        //std::shared_ptr<MycorrhizalRootRandomParameter> getMycorrhizalRootRandomParameter() const;
        //std::shared_ptr<const MycorrhizalRootSpecificParameter> param() const;
        //Mycparam


    };

}


#endif