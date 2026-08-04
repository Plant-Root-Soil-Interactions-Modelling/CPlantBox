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
    std::vector<double> colonizationTime;

    virtual ~MycorrhizalRoot() { };

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> rs) override;

    void simulate(double dt, bool silence = false) override; ///< root growth for a time span of @param dt
    void simulatePrimaryColonization(double dt);
    void simulateSecondaryColonization(double dt);
    void simulateHyphalGrowth(double dt, bool verbose);
    void simulateColonization(double dt, bool silence = false);

    double getParameter(std::string name) const override;

    void addNode(Vector3d n, int id, double t, size_t index, bool shift) override;
    void createLateral(double ageLN, bool silence) override;

    std::string toString() const override;

    std::shared_ptr<MycorrhizalRootRandomParameter> getRootRandomParameter() const;
    std::shared_ptr<const MycorrhizalRootSpecificParameter> param() const;

    int getNodeColonization(int i) const {return infected.at(i);}
    int getNumberofInfectedNodes() const;
    double getNodeColonizationTime(int i) const {return colonizationTime.at(i);}
    void setColonization(int i, int inf, double t);

private:
    bool isSecondaryColonizationAnchor(size_t node) const;
    double segmentLength(size_t a, size_t b) const;
    double colonizationArrivalTime(size_t fromNode, double cursegLength) const;
    void refineSegment(size_t fromNode, size_t toNode, double toCT, int &oldNode, int &currentNode, int step);
    void propagateSecondaryColonization(size_t anchorNode, int step, double maxLength, bool silence, double dt);

    double primaryColonizationRate(size_t node) const;
    int refinePrimarySegment(size_t fromNode, size_t toNode, double toCT, size_t insertIndex);

protected:
    void createHyphae(int pni);
    int hyphalTreeIndex = -1;
    void primaryColonization(double dt, bool silence);
    void secondaryColonization(bool silence, double dt);
};

}


#endif
