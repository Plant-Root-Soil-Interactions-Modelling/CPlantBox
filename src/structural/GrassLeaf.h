// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef GRASSLEAF_H_
#define GRASSLEAF_H_

#include "Organ.h"
#include "Organism.h"
#include "grassleafparameter.h"
#include "mymath.h"

#include <memory>
#include <string>

namespace CPlantBox {

class Plant;

/**
 * @brief A single grass leaf organ whose geometry is described by a Meristem polyline.
 *
 * The polyline is built by prepending nodes with addNodeFront():
 *
 *  - The first segment(s) added describe the **sheath** (the part wrapped around
 *    the stem); they grow upward relative to the stem.
 *  - Once the sheath is complete, subsequent segments describe the **midrib** (blade),
 *    which emerges at the ligule and bends outward with the blade angle.
 *
 * During simulate() the leaf advances through three phases:
 *  1. Sheath elongation  – sheath nodes are prepended until sheathLength is reached.
 *  2. Blade delay        – a waiting period before blade emergence.
 *  3. Blade elongation   – blade nodes are prepended until bladeLength is reached.
 *
 * getNode(i) delegates to the internal Meristem, returning Cartesian coordinates
 * by replaying turtle commands from the anchor.
 */
class GrassLeaf : public Organ
{
public:

    // ------------------------------------------------------------------ //
    //  Constructors
    // ------------------------------------------------------------------ //

    /// Restoration constructor (used when loading from file / copy).
    GrassLeaf(int id,
              std::shared_ptr<const OrganSpecificParameter> param,
              bool alive, bool active,
              double age, double length,
              Vector3d partialIHeading_,
              int pni,
              bool moved = false,
              int oldNON = 0);

    /// Simulation constructor – called by the plant during growth.
    GrassLeaf(std::shared_ptr<Organism> plant,
              int subtype,
              double delay,
              std::shared_ptr<Organ> parent,
              int pni);

    virtual ~GrassLeaf() {}

    // ------------------------------------------------------------------ //
    //  Organ interface
    // ------------------------------------------------------------------ //

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> plant) override;

    int organType() const override { return Organism::ot_leaf; }

    void simulate(double dt, bool verbose = false) override;

    /// Returns Cartesian node i from the internal Meristem.
    Vector3d getNode(int i) const override { return meristem.getNode(i); }

    int getNumberOfNodes() const { return meristem.size(); }

    double getParameter(std::string name) const override;

    std::string toString() const override;

    double calcLength(double age) override;
    double calcAge(double length) const override;

    // ------------------------------------------------------------------ //
    //  Grass-leaf specific accessors
    // ------------------------------------------------------------------ //

    /// Returns the current sheath length grown so far [cm].
    double getSheathLength() const { return sheathLengthGrown; }

    /// Returns the current blade length grown so far [cm].
    double getBladeLengthGrown() const { return bladeLengthGrown; }

    /// True once the sheath phase is complete.
    bool isSheathComplete() const { return sheathLengthGrown >= param()->sheathLength; }

    /// True once the blade delay has elapsed.
    bool isBladeEmerged() const { return isSheathComplete() && (age >= param()->bladeDelay + calcAgeAtSheathComplete()); }

    /// Read-only access to the underlying meristem.
    const Meristem& getMeristem() const { return meristem; }

    // ------------------------------------------------------------------ //
    //  Parameter convenience
    // ------------------------------------------------------------------ //

    std::shared_ptr<GrassLeafRandomParameter> getGrassLeafRandomParameter() const;
    std::shared_ptr<const GrassLeafSpecificParameter> param() const;

private:

    // ------------------------------------------------------------------ //
    //  Internal helpers
    // ------------------------------------------------------------------ //

    /// Age at which sheath elongation completes (approximate, from growth function).
    double calcAgeAtSheathComplete() const;

    /// Grows the sheath by dl [cm] by prepending nodes into the Meristem.
    void growSheath(double dl);

    /// Grows the blade by dl [cm] by prepending nodes into the Meristem.
    void growBlade(double dl);

    // ------------------------------------------------------------------ //
    //  State
    // ------------------------------------------------------------------ //

    Meristem meristem;             ///< Polyline of the entire leaf in turtle coordinates
    double sheathLengthGrown = 0.; ///< Accumulated sheath length [cm]
    double bladeLengthGrown  = 0.; ///< Accumulated blade length [cm]
};

} // end namespace CPlantBox

#endif
