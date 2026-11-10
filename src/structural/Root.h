// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ROOT_H_
#define ROOT_H_

#include "Organ.h"
#include "Organism.h"
#include "mymath.h"
#include "rootparameter.h"

namespace CPlantBox {

/**
 * Root
 *
 * Describes a single root, by a vector of nodes representing the root.
 *
 * The method simulate() creates new nodes of this root, and lateral roots in the root's branching zone.
 *
 * getNode(), getNodes(), heading(), getLength() are overwritten to use absolute coordinates,
 * for shootborne they are translated relative to the parent node.
 */
class Root : public Organ {

  public:
    Root(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length, Vector3d partialIHeading_, int pni, bool moved = false,
         int oldNON = 0);                                                                                      ///< creates everything from scratch
    Root(std::shared_ptr<Organism> plant_, int subType, double delay, std::shared_ptr<Organ> parent, int pni); ///< used within simulation
    virtual ~Root() {};                                                                                        ///< no need to do anything, children are deleted in ~Organ()

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> rs) override; ///< deep copies the root tree

    int organType() const override { return Organism::ot_root; }; ///< returns the organs type

    /* use absolute coordinates for root nodes, translate for shootborne roots */
    Vector3d getNode(int i) const override;          ///< returns the absolute coordinate of node @p i
    std::vector<Vector3d> getNodes() const override; ///< returns the absolute coordinates of all nodes
    Vector3d heading(int n) const override;          ///< returns the heading of node @p n
    double getLength(int i) const override;          ///< length from the first node up to node @p i (absolute nodes, or relative if shootborne)
    using Organ::getLength;                          ///< bring all other overloads of getLength from Organ into scope

    /* development */
    void simulate(double dt, bool silence = false) override; ///< root growth for a time span of @param dt

    /* last time step */
    bool hasMoved() const override { return moved || shootborne; } ///< nodes of shoot-borne roots move with the parent

    /* parameters */
    bool isShootborne() const { return shootborne; }                                                    ///< true if the root emerges from a shoot organ (stem or leaf)
    double getParameter(std::string name, std::map<std::string, double> addParams = {}) const override; ///< returns an organ pa:vector<CPlantBox::Vector3d>::size_type)'
    std::string toString() const override;

    /* from analytical equations */
    double calcLength(double age);                ///< analytical length of the root
    double calcAge(double length) const override; ///< analytical age of the root

    /* abbreviations */
    std::shared_ptr<RootRandomParameter> getRootRandomParameter() const; ///< root type parameter of this root
    std::shared_ptr<const RootSpecificParameter> param() const;          ///< root parameter

    /* state (de)serialisation */
    OrganState getState() const override;                                                                      ///< Adds the Root members to the state of Organ::getState()
    void setState(const OrganState &state, const std::map<int, std::shared_ptr<Organ>> &organs = {}) override; ///< Restores state, see Organ::setState()

    double insertionAngle = 0.; ///< differs to (const) theta, if angle is scaled by soil properties with RootRandomParameter::f_sa

  protected:
    bool shootborne = false; ///< true if the root emerges from a shoot organ (stem or leaf)

    void createSegments(double l, double dt, bool silence);           ///< Creates segments totalling length @p l; called from simulate()
    Vector3d getIncrement(const Vector3d &p, double sdx, int n = -1); ///< Returns the growth increment vector of length @p sdx from point @p p
};

/**
 * StaticRoot
 *
 * Describes a single root, by a vector of nodes representing the root.
 * The root is static (i.e. does not change over time), all laterals and corresponding emergence times are predefined
 *
 */
class StaticRoot : public Root {

  public:
    StaticRoot(int id, std::shared_ptr<const OrganSpecificParameter> param, double length, int pni);

    virtual ~StaticRoot() {}; ///< no need to do anything, children are deleted in ~Organ()

    void simulate(double dt, bool silence = false) override { Organ::simulate(dt, silence); }

    OrganState getState() const override;                                                                      ///< Adds the StaticRoot members to the state of Organ::getState()
    void setState(const OrganState &state, const std::map<int, std::shared_ptr<Organ>> &organs = {}) override; ///< Restores state, see Organ::setState()

    void initializeLaterals();

    void addLateral(int pni, int lateralType, double lateralDelay) {
        lateralNodeIndices.push_back(pni);
        lateralTypes.push_back(lateralType);
        lateralDelays.push_back(lateralDelay);
    }

  protected:
    std::vector<int> lateralNodeIndices;
    std::vector<int> lateralTypes;
    std::vector<double> lateralDelays;
};

} // end namespace CPlantBox

#endif /* ROOT_H_ */
