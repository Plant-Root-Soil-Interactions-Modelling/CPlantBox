#ifndef ROOT_H_
#define ROOT_H_

#include <iostream>
#include <assert.h>

#include "mymath.h"
#include "sdf.h"
#include "RootTropism.h"
#include "RootGrowth.h"
#include "ModelParameter.h"
#include "Plant.h"

class Plant;

/**
 * Root
 *
 * Describes a single root, by a vector of nodes representing the root.
 * The method simulate() creates new nodes of this root, and lateral roots in the root's branching zone.
 *
 */
class Organ
{

public:

    Organ(Plant* rs, int type, Vector3d pheading, double delay, Organ* parent, double pbl, int pni); ///< typically called by constructor of RootSystem, or Root::createLaterals()
    virtual ~Organ();

    void simulate(double dt, bool silence = false); ///< root growth for a time span of \param dt

    /* exact from analytical equations */
    double getCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
    double getLength(double age); ///< analytical length of the root
    double getAge(double length); ///< analytical age of the root

    RootTypeParameter* getRootTypeParameter() const;  ///< returns the root type parameter of the root
    double dx() const { return getRootTypeParameter()->dx; } ///< returns the axial resolution

    std::vector<Organ*> getRoots(); ///< return the root including laterals as sequential vector
    void getRoots(std::vector<Organ*>& v); ///< return the root system as sequential vector

    /* Nodes of the root */
    Vector3d getNode(int i) const { return nodes.at(i); } ///< i-th node of the root
    double getNodeETime(int i) const { return netimes.at(i); } ///< creation time of i-th node
    int getNodeId(int i) const {return nodeIds.at(i); } ///< unique identifier of i-th node
    size_t getNumberOfNodes() const {return nodes.size(); }  ///< return the number of the nodes of the root
    void addNode(Vector3d n,double t); //< adds a node to the root

    /* IO */
    void writeRSML(std::ostream & cout, std::string indent) const; ///< writes a RSML root tag
    std::string toString() const;

    Plant* rootsystem; ///< the root system this root is part of

    /* parameters that are given per root that are constant*/
    RootParameter param; ///< the parameters of this root
    Vector3d iheading; ///< the initial heading of the root, when it was created
    int id; ///< unique root id, (not used so far)
    double parent_base_length; ///< length [cm]
    int parent_ni; ///< parent node index

    /* parameters that are given per root that may change with time */
    bool alive = 1; ///< true: alive, false: dead
    bool active = 1; ///< true: active, false: root stopped growing
    double age = 0; ///< current age [days]
    double length = 0; ///< actual length [cm] of the root. might differ from getLength(age) in case of impeded root growth
    int old_non = 0; ///< index of the node that was update last time step (==0 if no update was performed)

    /* up and down */
    Organ* parent; ///< pointer to the parent root (equals nullptr if it is a base root)
    std::vector<Organ*> laterals; ///< the lateral roots of this root

    const double smallDx = 1e-6; ///< threshold value, smaller segments will be skipped (otherwise root tip direction can become NaN)

protected:

    void createSegments(double l, bool silence); ///< creates segments of length l, called by Root::simulate()
    void createLateral(bool silence); ///< creates a new lateral, called by Root::simulate()

    /* parameters that are given per node */
    std::vector<Vector3d> nodes; ///< nodes of the root
    std::vector<int> nodeIds; ///< unique node identifier
    std::vector<double> netimes; ///< node emergence times [days]

};

#endif /* ROOT_H_ */
