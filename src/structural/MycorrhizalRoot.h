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

        std::vector<int> infected;
        std::vector<double> infectionTime;

        virtual ~MycorrhizalRoot() { };

        std::shared_ptr<Organ> copy(std::shared_ptr<Organism> rs) override;
        int organType() const override { return Organism::ot_root; };

        void simulate(double dt, bool silence = false) override; ///< root growth for a time span of @param dt
        
        double getParameter(std::string name) const override;

        void addNode(Vector3d n, int id, double t, size_t index, bool shift) override;
        void createLateral(double ageLN, bool silence) override;
        // toString
        std::shared_ptr<MycorrhizalRootRandomParameter> getMycorrhizalRootRandomParameter() const;
        std::shared_ptr<const MycorrhizalRootSpecificParameter> param() const;

        int getNodeInfection(int i) const {std::cout<< "MycorrhizalRoot::getNodeInfection called"<<std::endl; return infected.at(i);}
        void setInfection(int i, int inf, double t);


    };

}


#endif