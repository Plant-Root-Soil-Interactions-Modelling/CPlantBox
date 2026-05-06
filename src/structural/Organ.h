// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ORGAN_H_
#define ORGAN_H_

/**
 * @file Organ.h
 * @brief Defines the Organ base class for all plant organs (seed, root, stem, leaf).
 *
 * An Organ stores its identity, geometry (nodes, node IDs, creation times), growth state,
 * and parent/child tree links.  Derived classes specialise simulate(), organType(), and
 * geometry helpers.  Post-processing utilities (getOrgans, getParameter, writeRSML) are
 * provided here so they work uniformly across all organ types.
 */

// project headers
#include "mymath.h"
#include "tinyxml2.h"

// standard library
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <vector>

namespace CPlantBox {

class OrganSpecificParameter;
class OrganRandomParameter;
class Organism;
class Plant;

/**
 * @brief Base class for all plant organs (seed, root, stem, leaf).
 *
 * Stores organ identity, geometry (nodes, node IDs, creation times), growth state,
 * and parent/child tree links.  Simulation, organ-type classification, and geometry
 * helpers are virtual so derived classes can specialise them.
 *
 * Post-processing utilities (getOrgans(), getParameter(), writeRSML()) are implemented
 * here and work uniformly across all organ types.
 *
 * @see Root, Stem, Leaf, Seed, Organism
 */
class Organ : public std::enable_shared_from_this<Organ> {

  public:
    Organ(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length, Vector3d partialIHeading_, int pni, bool moved = false, int oldNON = 0); ///< Constructs an organ from explicit state; tree links must be set by the caller
    Organ(std::shared_ptr<Organism> plant, std::shared_ptr<Organ> parent, int organtype, int subtype, double delay, int pni); ///< Simulation constructor: draws parameters via OrganRandomParameter::realize() and starts with age = -delay
    virtual ~Organ() {}

    virtual std::shared_ptr<Organ> copy(std::shared_ptr<Organism> plant); ///< Deep-copies the organ subtree into @p plant; caller owns the returned root

    virtual int organType() const; ///< Returns the coarse organ type (ot_organ=0, ot_seed=1, ot_root=2, ot_stem=3, ot_leaf=4); override in each derived class

    /* development */
    virtual void simulate(double dt, bool verbose = false); ///< Advances organ development by @p dt days; override in derived classes

    /* tree */
    void setOrganism(std::shared_ptr<Organism> p) { plant = p; }           ///< Sets the organism that owns this organ
    std::shared_ptr<Organism> getOrganism() const { return plant.lock(); } ///< Returns the owning organism
    void setParent(std::shared_ptr<Organ> p) { parent = p; }               ///< Sets the parent organ
    std::shared_ptr<Plant> getPlant() const;                               ///< Returns the owning organism downcast to Plant (nullptr if not a Plant)
    std::shared_ptr<Organ> getParent() const { return parent.lock(); }     ///< Returns the parent organ (nullptr if none)
    void addChild(std::shared_ptr<Organ> c);                               ///< Appends @p c as a child and sets its parent pointer
    int getNumberOfChildren() { return children.size(); }                  ///< Returns the number of direct child organs
    std::shared_ptr<Organ> getChild(int i) { return children.at(i); }      ///< Returns the child at local index @p i

    /* parameters */
    int getId() const { return id; }                                               ///< Returns the globally unique organ id
    std::shared_ptr<const OrganSpecificParameter> param() const { return param_; } ///< Returns the realised (instance-level) organ parameters
    std::shared_ptr<OrganRandomParameter> getOrganRandomParameter() const;         ///< Returns the stochastic parameter set for this organ's sub-type
    bool isAlive() const { return alive; }                                         ///< Returns true if the organ is alive
    bool isActive() const { return active; }                                       ///< Returns true if the organ is still growing
    double getAge() const { return age; }                                          ///< Returns the current age [days]
    double getLength(bool realized = true) const;                                  ///< Returns organ length [cm]; realized=true subtracts unresolved residual growth epsilonDx
    double getLength(int i) const;                  ///< Returns length from the first node up to node index @p i [cm]
    double getEpsilon() const { return epsilonDx; } ///< Returns residual growth not yet converted to a new segment [cm]
    int getParentNI() const { return parentNI; }    ///< Returns the local node index in the parent organ where this organ is attached
    virtual double calcAge(double length) const { throw std::runtime_error("calcAge() not implemented"); }  ///< Maps length -> age; must be overridden; needed by getOrgans()
    virtual double calcLength(double age) { throw std::runtime_error("calcLength() not implemented"); }     ///< Maps age -> length; must be overridden

    /* geometry */
    int getNumberOfNodes() const { return nodes.size(); }                  ///< Returns the number of nodes
    int getNumberOfSegments() const { return nodes.size() - 1; }           ///< Returns the number of polyline segments (nodes - 1)
    Vector3d getOrigin() const { return getParent()->getNode(parentNI); }; ///< Returns the absolute coordinate of the organ's attachment point on the parent
    virtual Vector3d getNode(int i) const { return nodes.at(i); }          ///< Returns the absolute (or relative) coordinate of node @p i
    int getNodeId(int i) const { return nodeIds.at(i); }                   ///< Returns the global node index of local node @p i
    std::vector<int> getNodeIds() const { return nodeIds; }                ///< Returns all global node indices
    const std::vector<Vector3d> &getNodes() const { return nodes; }        ///< Returns all node coordinates
    double getNodeCT(int i) const { return nodeCTs.at(i); }                ///< Returns the creation time of node @p i [days]
    void addNode(Vector3d n, double t, size_t index, bool shift);          ///< Adds a node using the next global index from the organism
    virtual void addNode(Vector3d n, int id, double t, size_t index, bool shift); ///< Adds a node with an explicit global id; override in Stem for phytomere handling
    void addNode(Vector3d n, int id, double t) { addNode(n, id, t, size_t(0), false); } ///< pybind convenience overload (no index or shift)
    void addNode(Vector3d n, double t) { addNode(n, t, size_t(0), false); };            ///< pybind convenience overload (no id, index, or shift)
    std::vector<Vector2i> getSegments() const;                                          ///< Returns all polyline segments as pairs of global node indices

    double dx() const;                                 ///< Returns the maximal axial segment length from the organ random parameter [cm]
    double dxMin() const;                              ///< Returns the minimal axial segment length from the organ random parameter [cm]
    void rel2abs();                                    ///< Converts node coordinates from relative to absolute (recursive over children)
    void abs2rel();                                    ///< Converts node coordinates from absolute to relative (recursive over children)
    void moveOrigin(int idx);                          ///< Updates the parent attachment node index (used during internodal growth)
    double calcCreationTime(double length, double dt); ///< Returns the analytical creation time of a node at @p length along the organ [days]

    /* last time step */
    virtual bool hasMoved() const { return moved; };             ///< Returns true if any node moved during the last simulate() call
    int getOldNumberOfNodes() const { return oldNumberOfNodes; } ///< Returns the node count before the last simulate() call

    /* for post processing */
    std::vector<std::shared_ptr<Organ>> getOrgans(int ot = -1, bool all = false);        ///< Returns this organ and all descendants with >1 node; ot=-1 means all types
    void getOrgans(int otype, std::vector<std::shared_ptr<Organ>> &v, bool all = false); ///< Appends this organ and descendants to @p v (in-place variant)
    virtual double getParameter(std::string name) const;                                 ///< Returns a named scalar parameter; delegates to OrganRandomParameter if not found locally
    int getNumberOfLaterals() const; ///< Returns the number of children whose age > 0

    /* IO */
    virtual std::string toString() const;                                           ///< Returns a compact human-readable summary for debugging
    void writeRSML(tinyxml2::XMLDocument &doc, tinyxml2::XMLElement *parent) const; ///< Appends an RSML @c <root> element for this organ to @p parent; not exposed to Python

    /* useful */
    virtual Vector3d heading(int n) const; ///< Returns the absolute heading at node @p n (direction of segment [n-1, n], or initial heading if n==0)
    Vector3d getiHeading0() const;         ///< Returns the initial absolute heading when the organ was created
    bool hasRelCoord() const;              ///< Returns true if node coordinates are currently stored in relative form
    bool has_rel_coord = false;            ///< Flag: true when nodes are in relative coordinates

    int parentNI; ///< Local node index in the parent organ where this organ is attached

    /* for carbon-limited growth (know future (or past) volume (or length))*/
    virtual double orgVolume(double length_ = -1., bool realized = false) const; ///< Returns organ volume [cm³] for the current or a specified length; cylinder approximation by default
    virtual double orgVolume2Length(double volume_) {
        return volume_ / (M_PI * getParameter("radius") * getParameter("radius"));
    } ///< Inverse of orgVolume(): returns the organ length [cm] that corresponds to @p volume_

  protected:
    mutable Vector3d partialIHeading; ///< Initial local heading used to compute getiHeading0(); mutable to allow derived-class caching

    virtual void createLateral(double ageLN, bool silence); ///< Creates a new lateral organ at the current tip; called from simulate(); override in RootDelay
    virtual void storeLinkingNodeLocalId(int numCreatedLN, bool silence) { ; }; ///< Stores the local id of the newly created linking node; override in Stem for nodal growth
    virtual Vector3d getIncrement(const Vector3d &p, double sdx, int n = -1); ///< Returns the growth increment vector of length @p sdx from point @p p; override in Leaf
    void createSegments(double l, double dt, bool silence, int phytoIdx = -1); ///< Creates segments totalling length @p l; called from simulate(); phytoIdx>=0 for internodal growth
    virtual double getLatInitialGrowth(double dt);                             ///< Returns the initial growth period passed to a newly created lateral [days]
    virtual double getLatGrowthDelay(int ot_lat, int st_lat, double dt) const; ///< Returns the emergence delay for a new lateral of type (@p ot_lat, @p st_lat) [days]
    bool getApplyHere(int i) const;                                            ///< Returns true if successor rule @p i applies at the current linking node
    std::vector<int> addIncrement;                                             ///< Scratch indices tracking inserted nodes during coordinate-system switching

    /* Up and down the organ tree */
    std::weak_ptr<Organism> plant;                ///< Owning organism (weak to avoid reference cycles)
    std::weak_ptr<Organ> parent;                  ///< Parent organ (nullptr if this is a base organ)
    std::vector<std::shared_ptr<Organ>> children; ///< Direct child organs (laterals)

    /* Parameters that are constant over the organ life time */
    const int id;                                         ///< Globally unique organ id assigned at construction
    std::shared_ptr<const OrganSpecificParameter> param_; ///< Realised parameter set for this organ instance

    /* Parameters changing over time */
    bool alive = true;               ///< true while the organ is alive
    bool active = true;              ///< true while the organ is growing
    double age = 0;                  ///< Current age [days]
    double length = 0;               ///< Theoretical organ length [cm] (may exceed realised length by epsilonDx)
    double epsilonDx = 0;            ///< Residual growth too small for a new segment; carried forward to the next time step [cm]
    size_t created_linking_node = 0; ///< Number of linking nodes (lateral attachment points) created so far

    /* Node data */
    std::vector<Vector3d> nodes; ///< Node coordinates [cm] (absolute or relative, see has_rel_coord)
    std::vector<int> nodeIds;    ///< Global node indices corresponding to each node
    std::vector<double> nodeCTs; ///< Node creation times [days]

    /* Last time step */
    bool moved = false;       ///< True if any node moved during the last simulate() call
    int oldNumberOfNodes = 0; ///< Node count at the start of the last simulate() call
    bool firstCall = true;    ///< True until createSegments() is called for the first time
};

} // namespace CPlantBox

#endif
