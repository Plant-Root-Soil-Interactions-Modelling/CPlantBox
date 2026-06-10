// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef GRASSLEAF_H_
#define GRASSLEAF_H_

#include "Organ.h"
#include "Organism.h"
#include "grassleafparameter.h"
#include "mymath.h"
#include "Plant.h"

#include <memory>
#include <string>

namespace CPlantBox {

class Plant;

/**
 * @brief A grass leaf organ whose geometry is stored as a TurtlePolyline.
 *
 * The leaf is divided into two zones based on the ratio sheathLength/totalLength:
 *  - Sheath: straight segments (pitch = 0) wrapping around the stem.
 *  - Blade: gently curving segments (pitch proportional to segment length) emerging above the leaf collar.
 *
 * Growth is driven by a single logistic/linear growth function (f_gf) parameterised by
 * leafGrowthDuration and the total leaf length (sheathLength + bladeLength).  New segments
 * are prepended at the base each time step (intercalary meristem model).  After each growth
 * step the pitch of every segment is updated by fixPitch() to reflect the current
 * sheath/blade partition.
 *
 * Node 0 is the attachment point on the parent organ; turtle nodes are indexed from 1.
 */
class GrassLeaf : public Organ
{
public:

    GrassLeaf(int id,
              std::shared_ptr<const OrganSpecificParameter> param,
              bool alive, bool active,
              double age, double length,
              Vector3d partialIHeading_,
              int pni,
              bool moved = false,
              int oldNON = 0); ///< Restoration constructor; tree links and geometry must be set by the caller.

    GrassLeaf(std::shared_ptr<Organism> plant,
              int subtype,
              double delay,
              std::shared_ptr<Organ> parent,
              int pni); ///< Simulation constructor – called by the plant during growth.

    virtual ~GrassLeaf() {}

  
    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> plant) override;

    int organType() const override { return Organism::ot_leaf; }

    void simulate(double dt, bool verbose = false) override;

    void abs2rel() override {} ///< No-op: GrassLeaf geometry is meristem-managed, not nodes-based.
    void rel2abs() override {} ///< No-op: GrassLeaf geometry is meristem-managed, not nodes-based.

    int getNumberOfNodes() const override { return turtle.size()+1; } ///< +1 because the anchor point is not stored in the meristem
    Vector3d getNode(int i) const override; ///< Returns Cartesian node i from the internal Meristem, anchored at the parent attachment point.

    double getParameter(std::string name) const override;

    std::string toString() const override;
   
    double calcLength(double age) override; 
    double calcAge(double length) const override; 

    double getSheathLength() const { return sheathLength; }      ///< Returns the current sheath length grown so far [cm].
    double getBladeLength() const { return bladeLength; }   ///< Returns the current blade length grown so far [cm].

    std::shared_ptr<GrassLeafRandomParameter> getGrassLeafRandomParameter() const;
    std::shared_ptr<const GrassLeafSpecificParameter> param() const;

    std::vector<Vector3d> getLeafVis(int i);    ///< 3D edge points at node i; returns the two leaf edge points for both blade and sheath nodes.

private:

    void growLeaf(double dl, double dt); ///< Prepends new base segments totalling dl [cm] to the TurtlePolyline.
    void fixPitch(); ///< Assigns pitch angles to every turtle node based on its position along the leaf (sheath vs blade).

    mutable TurtlePolyline turtle; ///< Turtle-coordinate polyline of the entire leaf; mutable so getNode() can set the anchor.
    double sheathLength = 0.;        ///< Current sheath length = length * sheathLength/(sheathLength+bladeLength) [cm]
    double bladeLength  = 0.;        ///< Current blade length  = length * bladeLength/(sheathLength+bladeLength) [cm]

};

} // end namespace CPlantBox

#endif
