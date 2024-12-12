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
        Vector3d partialIHeading_, int pni, bool infected = false, bool moved= false, int oldNON = 0): 
        Root(id,param, alive,active,age,length,partialIHeading_,pni, moved,oldNON), infected(infected) {};

        MycorrhizalRoot(std::shared_ptr<Organism> rs, int type, double delay, std::shared_ptr<Organ> parent, bool infected, int pni):
        Root(rs,type, delay, parent,pni), infected(infected) {};

        bool infected = false;

        virtual ~MycorrhizalRoot() { };

        // std::shared_ptr<Organ> copy(std::shared_ptr<Organism> rs) override;
        // organType ?
        void simulate(double dt, bool silence = false) override; ///< root growth for a time span of @param dt
        double getParameter(std::string name) const override;
        // toString
        std::shared_ptr<MycorrhizalRootRandomParameter> getRootRandomParameter() const;
        std::shared_ptr<const MycorrhizalRootSpecificParameter> param() const;


    };

}


#endif