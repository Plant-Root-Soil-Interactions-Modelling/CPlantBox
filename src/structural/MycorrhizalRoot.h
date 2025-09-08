#ifndef MYCORRHIZALROOT_H_
#define MYCORRHIZALROOT_H_

#include "mycorrhizalrootparameter.h"
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
    std::vector<int> emergedHyphae;
    std::vector<double> infectionTime;

    virtual ~MycorrhizalRoot() { };

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> rs) override;

    void simulate(double dt, bool silence = false) override; ///< root growth for a time span of @param dt
    void simulatePrimaryInfection(double dt);
    void simulateSecondaryInfection(double dt);
    void simulateHyphalGrowth();


    void simulateInfection(double dt, bool silence = false);

    double getParameter(std::string name) const override;

    void addNode(Vector3d n, int id, double t, size_t index, bool shift) override;
    void createLateral(double ageLN, bool silence) override;


    std::string toString() const override;

    std::shared_ptr<MycorrhizalRootRandomParameter> getRootRandomParameter() const;
    std::shared_ptr<const MycorrhizalRootSpecificParameter> param() const;

    int getNodeInfection(int i) const {return infected.at(i);}
    int getNumberofInfectedNodes() const;
    double getNodeInfectionTime(int i) const {return infectionTime.at(i);}
    void setInfection(int i, int inf, double t);

protected:
    void createHyphae(int pni);

    double prob(double  t, double segLength, double p);
    void primaryInfection(double dt, bool silence);
    void secondaryInfection(bool silence, double dt);

    void insertInfectedNode(int i);



};

}


#endif
