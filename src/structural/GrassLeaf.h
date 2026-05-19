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
 * @brief A grass leaf organ whose geometry is described by a polyline with Turtle coordinates.
 *
 * Growth phases in simulate():
 *  1. Elongation by creation of small segments (mimics cell division in the meristem)
 *  2. Seperation into sheath and blade segments (angular change at leaf collar)
 *  3. Blade and sheath grow by elongation of the existing segments
 * 
 * In reality cell divisions and elongation happen simultaneously, but we separate them here for simplicity. 
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

    int getNumberOfNodes() const override { return meristem.size()+1; } ///< +1 because the anchor point is not stored in the meristem
    Vector3d getNode(int i) const override; ///< Returns Cartesian node i from the internal Meristem, anchored at the parent attachment point.

    double getParameter(std::string name) const override;

    std::string toString() const override;
   
    double calcLength(double age) override; 
    double calcAge(double length) const override; 

    double getSheathLengthGrown() const { return sheathLengthGrown; }      ///< Returns the current sheath length grown so far [cm].
    double getBladeLengthGrown() const { return bladeLengthGrown; }   ///< Returns the current blade length grown so far [cm].

    std::shared_ptr<GrassLeafRandomParameter> getGrassLeafRandomParameter() const;
    std::shared_ptr<const GrassLeafSpecificParameter> param() const;

    std::vector<Vector3d> getLeafVis(int i);    ///< 3D edge points at node i: 2 for blade, 0 for sheath

private:

    int initialLength;  // legnth until all segments are created

    void growLeaf(double dl, double dt); ///< Grows the sheath by dl [cm] by prepending nodes into the Meristem.
    void elongateLeaf(double dl,double dt); ///< Elongates the leaf by dl [cm] by prepending nodes into the Meristem.

    mutable TurtlePolyline meristem;     ///< Polyline of the entire leaf in turtle coordinates
    double sheathLengthGrown = 0.; ///< Accumulated sheath length [cm]
    double bladeLengthGrown  = 0.; ///< Accumulated blade length [cm]
};

} // end namespace CPlantBox

#endif
