// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ORGANISM_H_
#define ORGANISM_H_

#include "mymath.h"
#include "sdf.h"

#include <array>
#include <iosfwd>
#include <map>
#include <memory>
#include <random>

namespace tinyxml2 {
class XMLDocument;
class XMLElement;
} // namespace tinyxml2

namespace CPlantBox {

class Organ;
class OrganRandomParameter;
class Seed;
class Root;
class Stem;
class Leaf;

/**
 * @brief Base container for plant organs and global simulation state.
 *
 * Manages organ parameter prototypes, simulation time stepping, geometry extraction,
 * and export utilities (RSML, VTP, geometry script).
 */
class Organism : public std::enable_shared_from_this<Organism> {
  public:
    enum DelayDefinition {
        dd_distance = 0,
        dd_time_lat = 1,
        dd_time_self = 2
    }; ///< definition of the growth delay. 0: distance based, 1: delay-based carried by the parent for all lateral, 2: delay-based carried by each lateral type
    enum OrganTypes { ot_organ = 0, ot_seed = 1, ot_root = 2, ot_stem = 3, ot_leaf = 4 }; ///< coarse organ classification
    static std::vector<std::string> organTypeNames;                                       ///< names of the organ types
    static int instances;                                                                 ///< the number of instances of this or derived classes

    static int organTypeNumber(std::string name); ///< Returns organ type enum from name.
    static std::string organTypeName(int ot);     ///< Returns organ type name from enum value.

    Organism(unsigned int seednum = 0); ///< Constructs an organism and initializes RNG state.
    virtual ~Organism() {};             ///< destructor

    virtual std::shared_ptr<Organism> copy(); ///< Deep-copies the organism including organ tree and parameters.

    /* organs */
    std::shared_ptr<Seed> getSeed(); ///< Returns the seed organ (expected at base index 0).

    /* organ parameter management */
    std::shared_ptr<OrganRandomParameter> getOrganRandomParameter(int ot, int subType) const; ///< Returns parameter set for organ type and subtype.
    std::vector<std::shared_ptr<OrganRandomParameter>> getOrganRandomParameter(int ot) const; ///< Returns all subtype parameter sets for organ type.
    void setOrganRandomParameter(std::shared_ptr<OrganRandomParameter> p);                    ///< Registers/overwrites one organ parameter set.
    int getParameterSubType(int organtype, std::string str) const;                            ///< Finds subtype index by subtype name.

    /* initialization and simulation */
    void setGeometry(std::shared_ptr<SignedDistanceFunction> geom) { geometry = geom; } ///< Sets optional confining geometry.
    virtual void initialize(bool verbose = true, std::string mode = "");                ///< Model-specific setup hook before simulation.
    virtual void simulate(double dt, bool verbose = false);                             ///< Advances all base organs by dt.
    double getSimTime() const { return simtime; }                                       ///< returns the current simulation time
    double getDt() const { return dt; }                                                 ///< returns the simulation time step of the last simulate() call
    void addOrgan(std::shared_ptr<Organ> o) { baseOrgans.push_back(o); }                ///< adds an organ, takes ownership

    /* as sequential list */
    std::vector<std::shared_ptr<Organ>> getOrgans(int ot = -1, bool all = false) const; ///< sequential list of organs
    virtual std::vector<double>
    getParameter(std::string name, int ot = -1,
                 std::vector<std::shared_ptr<Organ>> organs = std::vector<std::shared_ptr<Organ>>(0)) const; ///< parameter value per organ
    double getSummed(std::string name, int ot = -1) const;                                                   ///< summed up parameters

    /* geometry */
    int getNumberOfOrgans() const { return organId + 1; }               ///< Number of created organs.
    int getNumberOfNodes() const { return nodeId + 1; }                 ///< Number of created nodes.
    virtual int getNumberOfSegments(int ot = -1) const;                 ///< number of segments of the organism
    std::vector<std::vector<Vector3d>> getPolylines(int ot = -1) const; ///< nodes per organ
    std::vector<std::vector<double>> getPolylineCTs(int ot = -1) const; ///< node creation times per organ
    virtual std::vector<Vector3d> getNodes() const;                     ///< nodes of the organ
    virtual std::vector<double> getNodeCTs() const;                     ///< node creation times, corresponding to Organism::getNodes
    virtual std::vector<Vector2i> getSegments(int ot = -1) const;       ///< line segment containing two node indices, corresponding to Organism::getNodes
    virtual std::vector<double> getSegmentCTs(int ot = -1) const;       ///< line creation times, corresponding to Organism::getSegments
    virtual std::vector<int> getSegmentIds(int ot = -1) const;          ///< line segment indices, corresponding to Organism::getSegments
    virtual std::vector<std::shared_ptr<Organ>>
    getSegmentOrigins(int ot = -1) const; ///< points to the organ which contains the segment, corresponding to Organism::getSegments

    /* last time step */
    int getNumberOfNewNodes() const {
        return getNumberOfNodes() - oldNumberOfNodes;
    } ///< The number of new nodes created in the previous time step (same number as new segments)
    int getNumberOfNewOrgans() const { return getNumberOfOrgans() - oldNumberOfOrgans; } ///< The number of new roots created in the previous time step
    std::vector<int> getUpdatedNodeIndices() const;                                      ///< Indices of nodes that were updated in the previous time step
    std::vector<Vector3d> getUpdatedNodes() const;                                       ///< new coordinates of the updated nodes
    std::vector<double> getUpdatedNodeCTs() const;                                       ///< Updated creation times of moved nodes.
    std::vector<Vector3d> getNewNodes() const;                                           ///< nodes created in the previous time step
    std::vector<double> getNewNodeCTs() const;                                           ///< nodes created in the previous time step
    std::vector<Vector2i> getNewSegments(int ot = -1) const;                             ///< Segments created in the previous time step
    std::vector<std::shared_ptr<Organ>> getNewSegmentOrigins(int ot = -1) const;         ///< Copies a pointer to the root containing the new segments

    /* io */
    virtual std::string toString() const; ///< quick info for debugging
    virtual void initializeReader() {}    ///< initializes parameter reader
    virtual void readParameters(std::string name, std::string basetag = "plant", bool fromFile = true, bool verbose = false); ///< reads all organ type parameters from a xml file
    virtual std::string writeParameters(std::string name, std::string basetag = "plant", bool intoFile = true, bool comments = true) const;             ///< write all organ type parameters into a xml file
    virtual std::string write(std::string name, bool intoFile = true) const;     /// Writes output based on file extension.
    virtual void writeVTP(int otype, std::ostream &os) const;                    ///< Writes VTP polyline output.
    virtual void writeGeometry(std::ostream &os) const;                          ///< Writes confining geometry as ParaView Python script.
    virtual std::string writeRSML(std::string name, bool intoFile = true) const; ///< writes a RSML file
    int getRSMLSkip() const { return rsmlSkip; }                                 ///< skips points in the RSML output (default = 0)
    void setRSMLSkip(int skip) { assert(skip >= 0 && "Organism::setRSMLSkip(): skip must be >= 0"); rsmlSkip = skip; } ///< skips points in the RSML output (default = 0)
    std::vector<std::string> &getRSMLProperties() { return rsmlProperties; } ///< Returns mutable list of RSML property names.

    /* id management */
    int getOrganIndex() { organId++; return organId; } ///< returns next unique organ id, only organ constructors should call this
    int getNodeIndex() { nodeId++; return nodeId; } ///< returns next unique node id, only organ constructors should call this

    /* discretisation*/
    void setMinDx(double dx) { minDx = dx; } ///< minimum segment size, smaller segments will be skipped
    double getMinDx() { return minDx; }      ///< minimum segment size, smaller segments will be skipped

    /* random numbers */
    void setSeed(unsigned int seed); ///< Sets RNG seed.
    void setRandomSeed(unsigned int seed) { setSeed(seed); } ///< renamed setSeed to avoid confusion with plant seed, but keep setSeed for backward compatibility
    unsigned int getSeedVal() { return seed_val; }
    double rand() { if (stochastic) { return UD(gen); } else { return 0.5; } } ///< uniformly distributed random number [0, 1)
    double randn() { if (stochastic) { return ND(gen); } else { return 0.0; } } ///< normally distributed random number
    void setStochastic(bool stochastic_) { stochastic = stochastic_; } ///< Enables/disables stochastic sampling.
    bool getStochastic() { return stochastic; }                        ///< Returns stochastic sampling flag.

    std::vector<std::shared_ptr<Organ>> baseOrgans; ///< base organs of the organism (e.g. seed, or initial plant)

    virtual bool hasRelCoord() { return false; } ///< Overridden by derived classes with relative coordinates.
    int getDelayDefinition(int ot_lat);          ///< Returns delay-definition mode for a lateral organ type.

    int plantId; ///< Unique organism id (mainly for debugging/copy tracing).

    virtual std::shared_ptr<Seed> createSeed();                                                                    ///< Factory method for creating a seed organ.
    virtual std::shared_ptr<Root> createRoot(int subType, double delay, std::shared_ptr<Organ> parent, int pni);  ///< Factory method for creating roots, overridden by Plant to return Root instances.
    virtual std::shared_ptr<Stem> createStem(int subType, double delay, std::shared_ptr<Organ> parent, int pni);  ///< Factory method for creating stems, overridden by Plant to return Stem instances.
    virtual std::shared_ptr<Leaf> createLeaf(int subType, double delay, std::shared_ptr<Organ> parent, int pni);  ///< Factory method for creating leaves, overridden by Plant to return Leaf instances.

  protected:
    virtual tinyxml2::XMLElement *getRSMLMetadata(tinyxml2::XMLDocument &doc) const; ///< Builds RSML metadata node.
    virtual tinyxml2::XMLElement *getRSMLScene(tinyxml2::XMLDocument &doc) const;    ///< Builds RSML scene node.

    static const int numberOfOrganTypes = 5;
    std::array<std::map<int, std::shared_ptr<OrganRandomParameter>>, numberOfOrganTypes> organParam;

    std::shared_ptr<SignedDistanceFunction> geometry = std::make_shared<SignedDistanceFunction>(); ///< Confining geometry (unconfined by default)

    double simtime = 0;
    double dt = 0;
    int organId = -1;
    int nodeId = -1;
    int oldNumberOfNodes = 0;
    int oldNumberOfOrgans = 0;

    std::vector<std::string> rsmlProperties = {"organType", "subType", "length", "age", "parent-node", "diameter"};
    int rsmlSkip = 0;     ///< Point decimation factor for RSML output.
    double minDx = 1.e-6; ///< threshold value, smaller segments will be skipped, otherwise root tip direction can become NaN

    unsigned int seed_val; ///< value to use as seed, keep in memory to send to tropism
    std::mt19937 gen;
    std::uniform_real_distribution<double> UD = std::uniform_real_distribution<double>(0.0, 1.0);
    std::normal_distribution<double> ND = std::normal_distribution<double>(0.0, 1.0);
    bool stochastic = true; ///<  Whether to implement stochasticity
};

} // namespace CPlantBox

#endif
