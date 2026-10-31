#ifndef STEM_H_
#define STEM_H_

#include "Organ.h"
#include "Organism.h"
#include "stemparameter.h"

#include <assert.h>
#include <iostream>

namespace CPlantBox {

class Plant;

/**
 * Stem
 *
 * Describes a single stem, by a vector of nodes representing the stem.
 * The method simulate() creates new nodes of this stem, and lateral stems in the stem's branching zone.
 *
 */
class Stem : public Organ {

  public:
    Stem(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length, Vector3d partialIHeading_, int pni,
         bool moved = true, int oldNON = 0);
    Stem(std::shared_ptr<Organism> plant, int type, double delay, std::shared_ptr<Organ> parent, int pni); ///< used within simulation
    virtual ~Stem() {};

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> plant) override; ///< deep copies the root tree

    int organType() const override { return Organism::ot_stem; } ///< returns the organs type

    void simulate(double dt, bool silence = false) override;           ///< stem growth for a time span of \param dt
    void internodalGrowth(double dl, double dt, bool silence = false); ///< internodal growth of \param dl [cm]
    double getLatInitialGrowth(double dt) override;
    double getLatGrowthDelay(int ot_lat, int st_lat, double dt) const override;

    double getParameter(std::string name, std::map<std::string, double> addParams = {}) const override; ///< returns an organ parameter
    std::string toString() const override;

    /* exact from analytical equations */
    double calcLength(double age) override;       ///< analytical length of the stem
    double calcAge(double length) const override; ///< analytical age of the stem

    /* abbreviations */
    std::shared_ptr<StemRandomParameter> getStemRandomParameter() const; ///< root type parameter of this root
    std::shared_ptr<const StemSpecificParameter> param() const;          ///< root parameter

    OrganState getState() const override; ///< Adds the Stem members to the state of Organ::getState()
    void setState(const OrganState &state, const std::map<int, std::shared_ptr<Organ>> &organs = {}) override; ///< Restores state, see Organ::setState()

    std::vector<int> stemphytomerId = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int shootborneType = 5;

    Vector3d getNode(int i) const override { return nodes.at(i); }                 ///< i-th node of the organ
    void addNode(Vector3d n, int id, double t, size_t index, bool shift) override; //< adds a node to the root

  protected:
    void createLateral(double dt, bool verbose) override;         ///< creates the laterals and stores the local id of the new linking node
    void storeLinkingNodeLocalId(int numCreatedLN, bool verbose); ///< stores the local id of the newly created linking node
    std::vector<int> localId_linking_nodes;
    void minusPhytomerId(int subtype) { stemphytomerId[subtype]--; }
    int getphytomerId(int subtype) { return stemphytomerId[subtype]; }
    void addPhytomerId(int subtype) { stemphytomerId[subtype]++; }
};

} // namespace CPlantBox

#endif
