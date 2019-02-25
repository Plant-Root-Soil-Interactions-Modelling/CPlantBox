// *** ADDED BY HEADER FIXUP ***
#include <cassert>
#include <istream>
// *** END ***
#ifndef ROOT_H_
#define ROOT_H_

#include <iostream>
#include <assert.h>

#include "Organ.h"
#include "mymath.h"
#include "sdf.h"
#include "RootTropism.h"
#include "RootGrowth.h"
#include "ModelParameter.h"

namespace CPlantBox {


class Plant;
class RootParameter;
class RootTypeParameter;

/**
 * Root
 *
 * Describes a single root, by a vector of nodes representing the root.
 * The method simulate() creates new nodes of this root, and lateral roots in the root's branching zone.
 *
 * Since roots in soil cannot move the roots nodes are represented by absolute coordinates, in contrast to
 * other organs like Leaf or Stem.
 *
 */
class Root : public Organ
{

public:

    Root(Plant* plant, Organ* parent, int subtype, double delay, Vector3d iheading, int pni, double pbl); ///< typically called by constructor of RootSystem, or Root::createLaterals()
	Root(const Organ& r, Plant& plant);
	virtual ~Root() { }; // base class constructor is called automatically in c++

    virtual int organType() const override { return Organ::ot_root; };

    /* simulation */
    virtual void simulate(double dt, bool silence = false) override; ///< root growth for a time span of \param dt

    /* get results */
    virtual double getScalar(std::string name) const override; ///< returns an organ parameter of Plant::ScalarType

    /* exact from analytical equations */
    double getCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
    double getLength(double age); ///< analytical length of the root
    double getAge(double length); ///< analytical age of the root

    /* abbreviations */
    RootParameter* rParam() const { return (RootParameter*)param;  } ///< type cast
    RootTypeParameter* tParam() const; // type cast
    double dx() const; ///< returns the axial resolution
    Vector3d heading() const; /// current heading of the root tip

    /* IO */
    void writeRSML(std::ostream & cout, std::string indent) const; ///< writes a RSML root tag
    std::string toString() const;

    /* nodes */
    void addNode(Vector3d n, double t); //< adds a node to the root

    /* parameters that are given per root that are constant*/
	double parent_base_length; ///< length [cm]
	int parent_ni; ///< parent node index
	int pni;
	double pbl;


    const double smallDx = 1e-6; ///< threshold value, smaller segments will be skipped (otherwise root tip direction can become NaN)

    Vector3d initialHeading;

protected:

    void createSegments(double l, bool silence); ///< creates segments of length l, called by Root::simulate()
    void createLateral(bool silence); ///< creates a new lateral, called by Root::simulate()

    int old_non = 0;

};





} // namespace CPlantBox

#endif /* ROOT_H_ */
