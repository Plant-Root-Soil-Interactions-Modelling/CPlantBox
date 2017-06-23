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

class Plant;

/**
 * Root
 *
 * Describes a single root, by a vector of nodes representing the root.
 * The method simulate() creates new nodes of this root, and lateral roots in the root's branching zone.
 *
 */
class Root : public Organ
{

public:

    Root(Plant* plant, Organ* parent, int type, double delay, Matrix3d heading, int pni, double pbl); ///< typically called by constructor of RootSystem, or Root::createLaterals()
    virtual ~Root();

    void simulate(double dt, bool silence = false); ///< root growth for a time span of \param dt

    /* exact from analytical equations */
    double getCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
    double getLength(double age); ///< analytical length of the root
    double getAge(double length); ///< analytical age of the root

    double dx() const { return ((RootTypeParameter*)getOrganTypeParameter())->dx; } ///< returns the axial resolution

    std::vector<Organ*> getRoots(); ///< return the root including laterals as sequential vector
    void getRoots(std::vector<Organ*>& v); ///< return the root system as sequential vector

    /* Nodes of the root */
    void addNode(Vector3d n,double t); //< adds a node to the root

    /* IO */
    void writeRSML(std::ostream & cout, std::string indent) const; ///< writes a RSML root tag
    std::string toString() const;

    /* parameters that are given per root that are constant*/
    double parent_base_length; ///< length [cm]
    int parent_ni; ///< parent node index

    const double smallDx = 1e-6; ///< threshold value, smaller segments will be skipped (otherwise root tip direction can become NaN)

protected:

    void createSegments(double l, bool silence); ///< creates segments of length l, called by Root::simulate()
    void createLateral(bool silence); ///< creates a new lateral, called by Root::simulate()

};

#endif /* ROOT_H_ */
