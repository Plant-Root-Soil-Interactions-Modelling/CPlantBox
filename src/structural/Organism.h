// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ORGANISM_H_
#define ORGANISM_H_

#include "mymath.h"
#include "sdf.h"

#include "tinyxml2.h"

#include <chrono>
#include <random>
#include <map>
#include <array>
#include <memory>
#include <iostream>

namespace CPlantBox {

class Organ;
class OrganRandomParameter;
class Seed;

/**
 * Organism
 *
 *
 * Manages the OrganRandomParameters
 * Offers an interface for the simulation loop (initialize, simulate, ...)
 * Collects node and line segment geometry from the organ tree
 * Collect parameters from the organs
 * Can collect information about the last time step
 * Supports RSML
 * Holds global node index and organ index counter
 * Holds random numbers generator for the organ classes
 */
class Organism : public std::enable_shared_from_this<Organism> {
public:

	//0: distance based, 1: delay-based carried by the parent for all lateral, 2: delay-based carried by each lateral type
	enum DelayDefinition { dd_distance = 0, dd_time_lat = 1, dd_time_self = 2}; ///< definition of the growth delay
    enum OrganTypes { ot_organ = 0, ot_seed = 1, ot_root = 2, ot_stem = 3, ot_leaf = 4 }; ///< coarse organ classification
    static std::vector<std::string> organTypeNames; ///< names of the organ types
    static int instances; ///< the number of instances of this or derived classes

    static int organTypeNumber(std::string name); ///< organ type number from a string
    static std::string organTypeName(int ot); ///< organ type name from an organ type number

    Organism(unsigned int seednum  = 0); ///< constructor
    virtual ~Organism() { }; ///< destructor

    virtual std::shared_ptr<Organism> copy(); ///< deep copies the organism

    /* organs */
	std::shared_ptr<Seed> getSeed(); ///< the plant seed

    /* organ parameter management */
    std::shared_ptr<OrganRandomParameter> getOrganRandomParameter(int ot, int subType) const; ///< returns the respective the type parameter
    std::vector<std::shared_ptr<OrganRandomParameter>> getOrganRandomParameter(int ot) const; ///< returns all type parameters of an organ type (e.g. root)
    void setOrganRandomParameter(std::shared_ptr<OrganRandomParameter> p); ///< sets an organ type parameter, subType and organType defined within p
    int getParameterSubType(int organtype, std::string str) const; ///< returns the parameter sub type index of name @param str

    /* initialization and simulation */
    void setGeometry(std::shared_ptr<SignedDistanceFunction> geom) { geometry = geom; } ///< optionally, sets a confining geometry (call before RootSystem::initialize())
    void addOrgan(std::shared_ptr<Organ> o) { baseOrgans.push_back(o); } ///< adds an organ, takes ownership
    virtual void initialize(bool verbose = true); ///< overwrite for initialization jobs
    virtual void simulate(double dt, bool verbose = false); ///< calls the base organs simulate methods
    double getSimTime() const { return simtime; } ///< returns the current simulation time
    double getDt() const { return dt; } ///< returns the current simulation duration/time step

    /* as sequential list */
    std::vector<std::shared_ptr<Organ>> getOrgans(int ot=-1, bool all = false) const; ///< sequential list of organs
    virtual std::vector<double> getParameter(std::string name, int ot = -1, std::vector<std::shared_ptr<Organ>> organs = std::vector<std::shared_ptr<Organ>>(0)) const; ///< parameter value per organ
    double getSummed(std::string name, int ot = -1) const; ///< summed up parameters
    // std::shared_ptr<Organ> pickOrgan(int nodeId); // TODO

    /* geometry */
    int getNumberOfOrgans() const { return organId+1; } ///< number of nodes of the organism
    int getNumberOfNodes() const { return nodeId+1; } ///< number of nodes of the organism
    virtual int getNumberOfSegments(int ot=-1) const; ///< number of segments of the organism
    std::vector<std::vector<Vector3d>> getPolylines(int ot=-1) const; ///< nodes per organ
    std::vector<std::vector<double>> getPolylineCTs(int ot=-1) const; ///< node creation times per organ
    virtual std::vector<Vector3d> getNodes() const; ///< nodes of the organ
    virtual std::vector<double> getNodeCTs() const; ///< node creation times, corresponding to Organism::getNodes
    virtual std::vector<Vector2i> getSegments(int ot=-1) const; ///< line segment containing two node indices, corresponding to Organism::getNodes
    virtual std::vector<double> getSegmentCTs(int ot=-1) const; ///< line creation times, corresponding to Organism::getSegments
	virtual std::vector<int> getSegmentIds(int ot=-1) const; ///< line segment indices, corresponding to Organism::getSegments
    virtual std::vector<std::shared_ptr<Organ>> getSegmentOrigins(int ot=-1) const; ///< points to the organ which contains the segment, corresponding to Organism::getSegments

    /* last time step */
    int getNumberOfNewNodes() const { return getNumberOfNodes()- oldNumberOfNodes; } ///< The number of new nodes created in the previous time step (same number as new segments)
    int getNumberOfNewOrgans() const { return getNumberOfOrgans() - oldNumberOfOrgans; }  ///< The number of new roots created in the previous time step
    std::vector<int> getUpdatedNodeIndices() const; ///< Indices of nodes that were updated in the previous time step
    std::vector<Vector3d> getUpdatedNodes() const; ///< new coordinates of the updated nodes
    std::vector<double> getUpdatedNodeCTs() const; ///< new coordinates of the updated nodes
    std::vector<Vector3d> getNewNodes() const; ///< nodes created in the previous time step
    std::vector<double> getNewNodeCTs() const; ///< nodes created in the previous time step
    std::vector<Vector2i> getNewSegments(int ot=-1) const; ///< Segments created in the previous time step
    std::vector<std::shared_ptr<Organ>> getNewSegmentOrigins(int ot=-1) const; ///< Copies a pointer to the root containing the new segments

    /* io */
    virtual std::string toString() const; ///< quick info for debugging
    virtual void initializeReader() { } ///< initializes parameter reader
    virtual void readParameters(std::string name, std::string  basetag = "plant", bool fromFile = true, bool verbose = false); ///< reads all organ type parameters from a xml file
    virtual void writeParameters(std::string name, std::string basetag = "plant", bool comments = true) const; ///< write all organ type parameters into a xml file
    virtual void write(std::string name) const; /// writes simulation results (type is determined from file extension in name)
    virtual void writeVTP(int otype, std::ostream & os) const;
    virtual void writeGeometry(std::ostream & os) const;
    virtual void writeRSML(std::string name) const; ///< writes a RSML file
    int getRSMLSkip() const { return rsmlSkip; } ///< skips points in the RSML output (default = 0)
    void setRSMLSkip(int skip) { assert(rsmlSkip>=0 && "rsmlSkip must be >= 0" ); rsmlSkip = skip;  } ///< skips points in the RSML output (default = 0)
    std::vector<std::string>& getRSMLProperties() { return rsmlProperties; } ///< reference to the vector<string> of RSML property names, default is { "organType", "subType","length", "age"  }

    /* id management */
    int getOrganIndex() { organId++; return organId; } ///< returns next unique organ id, only organ constructors should call this
    int getNodeIndex() { nodeId++; return nodeId; } ///< returns next unique node id, only organ constructors should call this

    /* discretisation*/
    void setMinDx(double dx) { minDx = dx; } ///< Minimum segment size, smaller segments will be skipped
    double getMinDx() { return minDx; } ///< Minimum segment size, smaller segments will be skipped

    virtual void setSeed(unsigned int seed); ///< sets the seed of the organisms random number generator
    unsigned int getSeedVal(){return seed_val;}

	virtual double rand() { if (stochastic) { return UD(gen); } else { return 0.5; } }  ///< uniformly distributed random number [0, 1[
    virtual double randn() { if (stochastic) { return ND(gen); } else { return 0.0; } }  ///< normally distributed random number [-3, 3] in 99.73% of cases
	void setStochastic(bool stochastic_) { stochastic = stochastic_; }
	bool getStochastic(){ return stochastic; }
	std::vector<std::shared_ptr<Organ>> baseOrgans;  ///< base organs of the orgnism
	virtual bool hasRelCoord(){ return false; } ///< overriden by @Plant::hasRelCoord()
	int getDelayDefinition(int ot_lat);

    int plantId; // unique plant id (for debugging copy)

protected:

    virtual tinyxml2:: XMLElement* getRSMLMetadata(tinyxml2::XMLDocument& doc) const;
    virtual tinyxml2:: XMLElement* getRSMLScene(tinyxml2::XMLDocument& doc) const;

    static const int numberOfOrganTypes = 5;
    std::array<std::map<int, std::shared_ptr<OrganRandomParameter>>, numberOfOrganTypes> organParam;

    std::shared_ptr<SignedDistanceFunction> geometry = std::make_shared<SignedDistanceFunction>();  ///< Confining geometry (unconfined by default)

    double simtime = 0;
	double dt = 0;
    int organId = -1;
    int nodeId = -1;
    int oldNumberOfNodes = 0;
    int oldNumberOfOrgans = 0;

    std::vector<std::string> rsmlProperties = { "organType", "subType", "length", "age", "parent-node", "diameter" };
    int rsmlSkip = 0; // skips points
    double minDx = 1.e-6; ///< threshold value, smaller segments will be skipped, otherwise root tip direction can become NaN

	unsigned int seed_val;///<value to use as seed, keep in memory to send to tropism
    std::mt19937 gen;
    std::uniform_real_distribution<double> UD;
    std::normal_distribution<double> ND;
	bool stochastic = true;///<  Whether to implement stochasticity

};

} // namespace

#endif
