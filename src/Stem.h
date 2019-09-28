#ifndef STEM_H_
#define STEM_H_

#include <iostream>
#include <assert.h>

#include "Organ.h"
#include "Organism.h"
#include "stemparameter.h"

namespace CPlantBox {

static int phytomerId[10]= {0};

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

    Stem(int id,  std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
        Vector3d iheading, double pbl, int pni, bool moved= false, int oldNON = 0);
    Stem(std::shared_ptr<Organism> plant, int type, Vector3d pheading, double delay, std::shared_ptr<Organ> parent, double pbl, int pni); ///< used within simulation
    virtual ~Stem() { }; // base class constructor is called automatically in c++

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> plant) override;   ///< deep copies the root tree

    int organType() const override { return Organism::ot_stem; } ///< returns the organs type

    /* simulation */
    void simulate(double dt, bool silence = false) override; ///< stem growth for a time span of \param dt

    /* get results */
    double getParameter(std::string name) const override; ///< returns an organ parameter

    /* exact from analytical equations */
    double calcCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
    double calcLength(double age); ///< analytical length of the stem
    double calcAge(double length); ///< analytical age of the stem

    /* abbreviations */
    std::shared_ptr<StemRandomParameter> getStemRandomParameter() const;  ///< root type parameter of this root
    std::shared_ptr<const StemSpecificParameter> param() const; ///< root parameter

    double dx() const; ///< returns the axial resolution

    std::string toString() const;

    Vector3d iHeading; ///< the initial heading of the root, when it was created
    double parentBaseLength; ///< length [cm]
    int parentNI; ///< parent node index

protected:


    void minusPhytomerId(int subtype) { phytomerId[subtype]--;  }
    int getphytomerId(int subtype) { return phytomerId[subtype]; }
    void addPhytomerId(int subtype) { phytomerId[subtype]++;  }

    void createLateral(bool silence); ///< creates a new lateral, called by Stem::simulate()
    void leafGrow(bool silence, Vector3d bud);
    void shootBorneRootGrow(bool silence);

    Vector3d heading() const; /// current heading of the root tip
    virtual Vector3d getIncrement(const Vector3d& p, double sdx); ///< called by createSegments, to determine growth direction
    void createSegments(double l, bool silence); ///< creates segments of length l, called by stem::simulate()

    bool firstCall = true;
    const double smallDx = 1e-6; ///< threshold value, smaller segments will be skipped (otherwise stem tip direction can become NaN)
};

} // namespace CPlantBox

#endif
