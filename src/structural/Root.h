// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ROOT_H_
#define ROOT_H_

#include "mymath.h"
#include "Organ.h"
#include "Organism.h"
#include "rootparameter.h"

namespace CPlantBox {

class RootState;

/**
 * Root
 *
 * Describes a single root, by a vector of nodes representing the root.
 * The method simulate() creates new nodes of this root, and lateral roots in the root's branching zone.
 *
 */
class Root :public Organ
{

    friend RootState;

public:

    Root(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
        Vector3d partialIHeading_, int pni, bool moved= false, int oldNON = 0); // ///< creates everything from scratch
    Root(std::shared_ptr<Organism> rs, int type, double delay, std::shared_ptr<Organ> parent, int pni); ///< used within simulation
    virtual ~Root() { }; ///< no need to do anything, children are deleted in ~Organ()

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> rs) override;  ///< deep copies the root tree

    int organType() const override { return Organism::ot_root; }; ///< returns the organs type

    void simulate(double dt, bool silence = false) override; ///< root growth for a time span of @param dt

    double getParameter(std::string name) const override; ///< returns an organ pa:vector<CPlantBox::Vector3d>::size_type)’
    std::string toString() const override;

    /* From analytical equations */
    double calcLength(double age); ///< analytical length of the root
    double calcAge(double length) const; ///< analytical age of the root

    /* Abbreviations */
    std::shared_ptr<RootRandomParameter> getRootRandomParameter() const;  ///< root type parameter of this root
    std::shared_ptr<const RootSpecificParameter> param() const; ///< root parameter

    double insertionAngle=0.; ///< differs to (const) theta, if angle is scaled by soil properties with RootRandomParameter::f_sa TODO some better idea?

protected:


};

} // end namespace CPlantBox

#endif /* ROOT_H_ */
