#ifndef ORGAN_H_
#define ORGAN_H_

#include <iostream>
#include <assert.h>
#include <stdexcept>

#include "mymath.h"
#include "sdf.h"
#include "ModelParameter.h"

class Plant;

/**
 * Organ
 *
 * Base class of seed, root, shoot and leaf
 *
 */
class Organ
{

public:

    Organ(Plant* plant, Organ* parent, int type, double delay);
    virtual ~Organ();

    virtual int organType() const;///< overwrite for each organ

    Vector3d getRelativeOrigin(Vector3d o) const { return r_origin; };
    void setRelativeOrigin(Vector3d o) { r_origin = o; };
    virtual Matrix3d getRelativeInitialHeading() const { throw std::invalid_argument("Organ::getRelativeInitialHeading not implemented"); };
    virtual void setRelativeInitialHeading(Matrix3d m) { throw std::invalid_argument("Organ::setRelativeInitialHeading not implemented");  };
    Vector3d getAbsoluteOrigin() const;
    Matrix3d getAbsoluteInitialHeading() const;

    virtual void initialize() { }; ///< create call backs from the parameters set TODO
    virtual void simulate(double dt, bool silence = false) { age+=dt; }; ///< growth for a time span of \param dt
    virtual double getScalar(int stype);

    OrganTypeParameter* getOrganTypeParameter() const;  ///< returns the root type parameter of the root
    std::vector<Organ*> getOrgans(int otype); ///< return the organ including successors
    void getOrgans(int otype, std::vector<Organ*>& v); ///< returns the plant as sequential vector

    /* IO */
    virtual std::string toString() const;
    virtual void writeRSML(std::ostream & cout, std::string indent) const; ///< writes a RSML root tag

    /* up and down the organ tree */
    Plant* plant; ///< the plant of which this organ is part of
    Organ* parent; ///< pointer to the parent organ (equals nullptr if it has no parent)
    std::vector<Organ*> children; ///< the successive organs

    /* getter for the node data */
    size_t getNumberOfNodes() const {return r_nodes.size(); }  ///< return the number of the nodes of the root
    std::vector<Vector3d> getRelativeNodes() { return r_nodes; }
    std::vector<int> getNodeIDs() { return nodeIDs; }
    std::vector<double> getNodeETs() { return nctimes; }
    Vector3d getNode(int i) const { return r_nodes.at(i); } ///< i-th node of the root
    int getNodeID(int i) const {return nodeIDs.at(i); } ///< unique identifier of i-th node
    double getNodeCT(int i) const { return nctimes.at(i); } ///< emergence time of i-th node

    /* Parameters that are constant*/
    int id; ///< unique organ id, (not used so far)
    OrganParameter* param = nullptr; ///< the parameters of this root

    /* Parameters that may change with time */
    bool alive = 1; ///< true: alive, false: dead
    bool active = 1; ///< true: active, false: root stopped growing
    double age = 0; ///< current age [days]
    double length = 0; ///< actual length [cm] of the root. might differ from getLength(age) in case of impeded root growth

    /* node data (todo private?) */
    Vector3d r_origin =Vector3d(); // relative origin

    std::vector<Vector3d> r_nodes; ///< relative nodes of the root
    std::vector<int> nodeIDs; ///< unique node identifier
    std::vector<double> nctimes; ///< node creation times [days]

};

#include "Plant.h" // why?

#endif /* ORGAN_H_ */
