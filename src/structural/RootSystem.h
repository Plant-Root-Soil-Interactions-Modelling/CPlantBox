// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ROOTSYSTEM_H_
#define ROOTSYSTEM_H_

#include <fstream>

#include "soil.h"
#include "tropism.h"
#include "Organism.h"
#include "Root.h"
#include "seedparameter.h"
#include "Seed.h"

namespace CPlantBox {

class RootState;
class RootSystemState;

/**
 * RootSystem
 *
 * This class manages model parameter, the simulation,
 * stores the base roots of the root system,
 * and offers utility functions for post processing.
 * More post processing functions can be found in the class SegmentAnalyser
 */
class RootSystem :public Organism
{

    friend RootSystemState;

public:

    enum TropismTypes { tt_plagio = 0, tt_gravi = 1, tt_exo = 2, tt_hydro = 3 };  ///< root tropism types
    enum GrowthFunctionTypes { gft_negexp = 1, gft_linear = 2 }; // root growth function

    RootSystem(unsigned int seednum  = 0.); ///< empty root system
    virtual ~RootSystem() { };

    std::shared_ptr<Organism> copy() override; ///< deep copies the organism

    /* Parameter input output */
    std::shared_ptr<RootRandomParameter> getRootRandomParameter(int type) const;///< returns the i-th root parameter set (i=1..n)
    std::vector<std::shared_ptr<RootRandomParameter>> getRootRandomParameter() const; ///< all root type parameters as a vector
    void setRootSystemParameter(std::shared_ptr<SeedRandomParameter> rsp); ///< sets the root system parameters
    std::shared_ptr<SeedRandomParameter> getRootSystemParameter(); ///< gets the root system parameters

    void initializeReader() override; ///< initializes XML reader
    void readParameters(std::string name, std::string  basetag = "plant", bool fromFile = true, bool verbose = true) override { this->initializeReader(); Organism::readParameters(name, basetag, fromFile, verbose); };
    void openFile(std::string filename, std::string subdir="modelparameter/"); ///< DEPRICATED reads root parameter and plant parameter
    int readParameters(std::istream & cin); ///< DEPRICATED reads root parameters from an input stream
    void writeParameters(std::ostream & os) const; ///< DEPRICATED writes root parameters

    /* Simulation */
    void setSoil(std::shared_ptr<SoilLookUp> soil_) { soil = soil_; } ///< optionally sets a soil for hydro tropism (call before RootSystem::initialize())
    void reset(); ///< resets the root class, keeps the root type parameters
    virtual void initialize(bool verbose = true) override { initializeLB(4,5, verbose); };
    ///< creates the base roots, call before simulation and after setting the plant and root parameters
    virtual void initializeLB(int basal = 4, int shootborne = 5, bool verbose = true);
    ///< creates the base roots (length based lateral emergence times), call before simulation and after setting plant and root parameters
    virtual void initializeDB(int basal = 4, int shootborne = 5, bool verbose = true);
    ///< creates the base roots (delay based lateral emergence times), call before simulation and after setting plant and root parameters
    void setTropism(std::shared_ptr<Tropism> tf, int rt = -1); ///< sets a tropism function for a single root type or all root types (defaut)
    void simulate(double dt, bool verbose = false) override; ///< simulates root system growth for time span dt
    void simulate(); ///< simulates root system growth for the time defined in the root system parameters
    void simulate(double dt, double maxinc, ProportionalElongation* se, bool silence = false);
    ///< simulates the root system with a maximal overall elongation

    /* sequential */
    std::vector<std::shared_ptr<Root>> getRoots() const; ///< represents the root system as sequential vector of roots and buffers the result

    /* call back function creation */
    void initCallbacks(); ///< sets up callback functions for tropisms and growth functions, called by initialize()
    virtual std::shared_ptr<Tropism> createTropismFunction(int tt, double N, double sigma);
    ///< creates the tropisms, overwrite or change this method to add more tropisms
    virtual std::shared_ptr<GrowthFunction> createGrowthFunction(int gft);
    ///< creates the growth function per root type, overwrite or change this method to add more tropisms

    /* Analysis of simulation results */
    int getNumberOfSegments(int ot = -1) const override { return nodeId-numberOfCrowns-1; }
    ///< number of segments of the root system ((nid+1)-1) - numberOfCrowns - 1 (artificial shoot)
    int getNumberOfRoots(bool all = false) const { if (all) return organId+1; else return getRoots().size(); }
    std::vector<Vector3d> getNodes() const override;
    std::vector<std::shared_ptr<Organ>> getBaseRoots() const { return baseOrgans; } ///< base roots are tap root, basal roots, and shoot borne roots TODO
    std::vector<Vector2i> getShootSegments() const; ///< copies the segments connecting tap, basal root, shootborne roots
    std::vector<int> getRootTips() const; ///< node indices of the root tips
    std::vector<int> getRootBases() const; ///< node indices of the root bases

    /* dynamics */
    void push(); ///< push current state to a stack
    void pop(); ///< retrieve previous state from stack

    std::string toString() const override; ///< infos about current root system state (for debugging)

protected:

    int numberOfCrowns = 0;

private:

    void initialize_(int basal = 4, int shootborne = 5, bool verbose = true);

    std::shared_ptr<Seed> seed = nullptr;
    SeedSpecificParameter seedParam;

    std::shared_ptr<SoilLookUp> soil; ///< callback for hydro, or chemo tropism (needs to set before initialize()) TODO should be a part of tf, or rtparam

    mutable std::vector<std::shared_ptr<Root>> roots = std::vector<std::shared_ptr<Root>>(); // buffer for getRoots()

    std::vector<RootSystemState> stateStack;
};



/**
 * Sores a state of the RootSystem,
 * i.e. all data that changes over time (*), i.e. excluding node data that cannot change
 *
 * (*) excluding changes regarding RootSystemParameter, any RootTypeParameter, confining geometry, and soil
 */
class RootSystemState
{

    friend RootSystem;

public:

    RootSystemState(const RootSystem& rs); ///< create root system state from a rootsystem

    void restore(RootSystem& rs); ///< restore evolved rootsystem back to its previous state

private:

    std::vector<RootState> baseRoots;  ///< Base roots of the root system

    double simtime = 0;
    double dt = 0; //< last time step that was used
    int organId = -1;
    int nodeId = -1;
    int oldNumberOfNodes = 0;
    int oldNumberOfOrgans = 0;
    int numberOfCrowns = 0; ///< old number of root crowns

    mutable std::mt19937 gen; ///< random generator state
    mutable std::uniform_real_distribution<double> UD;  ///< random generator state
    mutable std::normal_distribution<double> ND; ///< random generator state

};



/**
 * Stores a state of the root that can be restored at a later point
 * (for RootSystem::push and RootSystem::pop)
 */
class RootState {

public:

    RootState() { };

    RootState(const Root& r); ///< create the root state from a root

    void restore(Root& r); ///< restore evolved root back to its previous state

private:

    /* parameters that are given per root that may change with time */
    bool alive = 1;
    bool active = 1;
    double age = 0;
    double length = 0;
    double epsilonDx = 0;
    bool moved = false;
    int oldNumberOfNodes = 0;
    bool firstCall = true;

    /* down the root branch*/
    std::vector<RootState> laterals = std::vector<RootState>(0); ///< the lateral roots of this root

    /* last node */
    Vector3d lNode = Vector3d(0.,0.,0.); ///< last node
    int lNodeId = 0; ///< last node id
    double lneTime = 0.;  ///< last creation time
    size_t non = 0; ///< number of nodes

};

} // end namespace CPlantBox

#endif /* ROOTSYSTEM_H_ */
