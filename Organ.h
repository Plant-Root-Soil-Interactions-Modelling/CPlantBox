// *** ADDED BY HEADER FIXUP ***
#include <cassert>
#include <istream>
// *** END ***
#ifndef ORGAN_H_
#define ORGAN_H_

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
    virtual Vector3d getRelativeOrigin() const
    {
        throw std::invalid_argument("Organ::getRelativeInitialHeading not implemented");
    };
    virtual void setRelativeOrigin(Vector3d o)
    {
        throw std::invalid_argument("Organ::getRelativeInitialHeading not implemented");
    };
    virtual Matrix3d getRelativeInitialHeading() const
    {
        throw std::invalid_argument("Organ::getRelativeInitialHeading not implemented");
    };
    virtual void setRelativeInitialHeading(Matrix3d m)
    {
        throw std::invalid_argument("Organ::setRelativeInitialHeading not implemented");
    };
    Vector3d getAbsoluteOrigin() const; ///< the absolute origin coordinates
    Matrix3d getAbsoluteInitialHeading() const; ///<  the absolute heading coordinates

    /* parameters */
    OrganTypeParameter* getOrganTypeParameter() const;  ///< organ type parameter

    /* simulation */
    virtual void simulate(double dt, bool silence = false)
    {
        age+=dt;
    }; ///< growth for a time span of \param dt

    /* post processing */
    std::vector<Organ*> getOrgans(unsigned int otype); ///< the organ including successors in a sequential vector
    void getOrgans(unsigned int otype, std::vector<Organ*>& v); ///< the organ including successors in a sequential vector

    virtual double getScalar(int stype) const; ///< returns an organ parameter of Plant::ScalarType

    /* IO */
    virtual std::string toString() const;
    virtual void writeRSML(std::ostream & cout, std::string indent) const; ///< writes a RSML root tag

    /* getter for the node data */
    size_t getNumberOfNodes() const
    {
        return r_nodes.size();    ///< return the number of the nodes of the root
    }
    std::vector<Vector3d> getRelativeNodes()
    {
        return r_nodes;
    }
    std::vector<int> getNodeIDs()
    {
        return nodeIDs;
    }
    std::vector<double> getNodeETs()
    {
        return nctimes;
    }
    Vector3d getNode(int i) const
    {
        return r_nodes.at(i);    ///< i-th node of the root TODO
    }
    int getNodeID(int i) const
    {
        return nodeIDs.at(i);    ///< unique identifier of i-th node
    }
    double getNodeCT(int i) const
    {
        return nctimes.at(i);    ///< emergence time of i-th node
    }

    /* up and down the organ tree */
    Plant* plant; ///< the plant of which this organ is part of
    Organ* parent; ///< pointer to the parent organ (equals nullptr if it has no parent)
    std::vector<Organ*> children; ///< the successive organs

    /* Parameters that are constant*/
    int id; ///< unique organ id, (not used so far)
    OrganParameter* root_param = nullptr; ///< the parameters of this root
    OrganParameter* stem_param = nullptr; ///< the parameters of this stem
    OrganParameter* leaf_param = nullptr;///< the parameters of this leaf


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
