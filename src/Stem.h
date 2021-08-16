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
    		Matrix3d iHeading, int pni, bool moved = false, int oldNON = 0);
    Stem(std::shared_ptr<Organism> plant, int type, Matrix3d iHeading, double delay, std::shared_ptr<Organ> parent, int pni); ///< used within simulation
    virtual ~Stem() { };

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> plant) override;   ///< deep copies the root tree

    int organType() const override { return Organism::ot_stem; } ///< returns the organs type

    void simulate(double dt, bool silence = false) override; ///< stem growth for a time span of \param dt

    Vector3d getNode(int i) const override { return rel2abs(nodes.at(i)); } ///< i-th node of the organ

    double getParameter(std::string name) const override; ///< returns an organ parameter

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
	Matrix3d inv; // inverse matrix of M

protected:

    Vector3d rel2abs(const Vector3d& n) const;
	Vector3d abs2rel(const Vector3d& n) const;

    void minusPhytomerId(int subtype) { phytomerId[subtype]--;  }
    int getphytomerId(int subtype) { return phytomerId[subtype]; }
    void addPhytomerId(int subtype) { phytomerId[subtype]++;  }

    void createLateral(bool silence); ///< creates a new lateral, called by Stem::simulate()
    void leafGrow(bool silence);
    void shootBorneRootGrow(bool silence);

    virtual Vector3d getIncrement(const Vector3d& p, double sdx); ///< called by createSegments, to determine growth direction
    void createSegments(double l, bool silence); ///< creates segments of length l, called by stem::simulate()

    bool firstCall = true;

};

} // namespace CPlantBox

#endif
