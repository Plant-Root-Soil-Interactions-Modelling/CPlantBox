#ifndef STEM_H_
#define STEM_H_

#include "Organ.h"
#include "Organism.h"
#include "stemparameter.h"

#include <iostream>
#include <assert.h>

namespace CPlantBox {

class Plant;

/**
 * Stem
 *
 * Describes a single stem, by a vector of nodes representing the stem.
 * The method simulate() creates new nodes of this stem, and lateral stems in the stem's branching zone.
 *
 */
class Stem : public Organ
{
public:

    static std::vector<int> phytomerId;

    Stem(int id,  std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
    		Vector3d partialIHeading_, int pni, bool moved = true, int oldNON = 0);
    Stem(std::shared_ptr<Organism> plant, int type, Matrix3d iHeading, double delay, std::shared_ptr<Organ> parent, int pni); ///< used within simulation
    virtual ~Stem() { };

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> plant) override;   ///< deep copies the root tree

    int organType() const override { return Organism::ot_stem; } ///< returns the organs type

    void simulate(double dt, bool silence = false) override; ///< stem growth for a time span of \param dt
	void internodalGrowth(double dl, bool silence = false); ///< internodal growth of \param dl [cm]
    Vector3d getNode(int i) const override { return nodes.at(i); } ///< i-th node of the organ

    double getParameter(std::string name) const override; ///< returns an organ parameter
	double getLength(bool realized = true) const; ///< length of the organ (realized => dependent on dx() and dxMin())
    double getLength(int i) const override;
    std::string toString() const override;

    /* exact from analytical equations */
    double calcCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
    double calcLength(double age); ///< analytical length of the stem
    double calcAge(double length); ///< analytical age of the stem

    /* abbreviations */
    std::shared_ptr<StemRandomParameter> getStemRandomParameter() const;  ///< root type parameter of this root
    std::shared_ptr<const StemSpecificParameter> param() const; ///< root parameter

    int shootborneType = 5;

	/* orientation */
	Vector3d heading(int n)  const override; ///< current (absolute) heading of the organs at node n
    Vector3d heading() const override {return heading( -1 ); } 
	Vector3d getiHeading0()  const override;
	
    void rel2abs() override;
	void abs2rel() override;
	bool hasMoved() const override { return true; }; ///< have any nodes moved during the last simulate call
																										 
protected:
	int getLeafSubType();
	Vector3d partialIHeading;
    void minusPhytomerId(int subtype) { phytomerId[subtype]--;  }
    int getphytomerId(int subtype) { return phytomerId[subtype]; }
    void addPhytomerId(int subtype) { phytomerId[subtype]++;  }

    void createLateral(bool silence); ///< creates a new lateral, called by Stem::simulate()
    void leafGrow(bool silence);
    void shootBorneRootGrow(bool silence);

    Vector3d getIncrement(const Vector3d& p, double sdx, int n = -1); ///< called by createSegments, to determine growth direction
    void createSegments(double l, bool silence, int PhytoIdx = -1); ///< creates segments of length l, called by stem::simulate()

    bool firstCall = true;

};

} // namespace CPlantBox

#endif
