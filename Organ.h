#ifndef ORGAN_H_
#define ORGAN_H_

#include <cassert>
#include <istream>
#include <iostream>
#include <assert.h>
#include <stdexcept>

#include "mymath.h"
#include "sdf.h"
#include "ModelParameter.h"

class Plant;
class OrganTypeParameter;
class OrganParameter;

/**
 * Organ
 *
 * Base class of seed, root, shoot and leaf
 *
 */
class Organ
{

public:

    enum OrganTypes { ot_seed = 1, ot_root = 2, ot_stem = 4, ot_leafe = 8, ot_shoot = ot_stem | ot_leafe, ot_organ = ot_seed | ot_root | ot_stem | ot_leafe}; ///< organ types bit wise

    Organ(Plant* plant, Organ* parent, int subtype, double delay);

    virtual ~Organ();

    virtual int organType() const { return Organ::ot_organ; }  ///< returns the organs type, overwrite for each organ

    /* scene graph for upper plant parts */

    virtual Vector3d getRelativeOrigin() const { return Vector3d(); };
    ///< the relative position within the parent organ
    virtual void setRelativeOrigin(const Vector3d& o) { throw std::invalid_argument("Organ::setRelativeOrigin not implemented"); };
    ///< the relative position within the parent organ
    virtual Matrix3d getRelativeHeading() const { return Matrix3d(); };
    ///< the heading in the parent organ
    virtual void setRelativeHeading(const Matrix3d& m) { throw std::invalid_argument("Organ::getRelativeHeading not implemented"); };
    ///< the heading in the parent organ

    Vector3d getOrigin() const; ///< absolute coordinates of the organ origin
    Matrix3d getHeading() const; ///< absolute heading of the organ

    Vector3d getNode(int i) const;  ///< i-th node of the organ in absolute coordinates

    /* parameters */
    OrganTypeParameter* getOrganTypeParameter() const;  ///< organ type parameter

    /* simulation */
    virtual void simulate(double dt, bool silence = false); ///< growth for a time span of \param dt

    /* post processing */
    std::vector<Organ*> getOrgans(unsigned int otype); ///< the organ including successors in a sequential vector
    void getOrgans(unsigned int otype, std::vector<Organ*>& v); ///< the organ including successors in a sequential vector

    virtual double getScalar(int stype) const; ///< returns an organ parameter of Plant::ScalarType

    /* IO */
    virtual std::string toString() const;
    virtual void writeRSML(std::ostream & cout, std::string indent) const; ///< writes a RSML root tag

    size_t getNumberOfNodes() const { return r_nodes.size(); } ///< number of nodes of the organ
    int getNodeID(int i) const { return nodeIDs.at(i); } ///< unique identifier of i-th node
    double getNodeCT(int i) const { return nctimes.at(i); } ///< creation time of i-th node

    /* up and down the organ tree */
    Plant* plant; ///< the plant of which this organ is part of
    Organ* parent; ///< pointer to the parent organ (equals nullptr if it has no parent)
    std::vector<Organ*> children; ///< the successive organs

    /* Parameters that are constant*/
    int id; ///< unique organ id, (not used so far)
    OrganParameter* param = nullptr; ///< the parameters of this root

    /* Parameters that may change with time */
    bool alive = 1; ///< true: alive, false: dead
    bool active = 1; ///< true: active, false: root stopped growing
    double age = 0; ///< current age [days]
    double length = 0; ///< actual length [cm] of the root. might differ from getLength(age) in case of impeded root growth

    /* node data */
    std::vector<Vector3d> r_nodes; ///< relative nodes of the root
    std::vector<int> nodeIDs; ///< unique node identifier
    std::vector<double> nctimes; ///< node creation times [days]

};

#include "Plant.h" // why?

#endif /* ORGAN_H_ */
