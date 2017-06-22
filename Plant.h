#ifndef ROOTSYSTEM_H_
#define ROOTSYSTEM_H_

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <chrono>
#include <random>
#include <numeric>

#include "ModelParameter.h"
#include "Organ.h"
#include "soil.h"

class Organ;
class TropismFunction;



/**
 * RootSystem
 *
 * This class manages all model parameter and the simulation,
 * stores the base roots of the root system,
 * and offers utility functions for post processing
 */
class Plant
{

  friend Organ;  // obviously :-)

public:

  enum TropismTypes { tt_plagio = 0, tt_gravi = 1, tt_exo = 2, tt_hydro = 3 };  ///< root tropism
  enum GrowthFunctionTypes { gft_negexp = 1, gft_linear = 2 }; // root growth function
  enum ScalarTypes { st_type = 0, st_radius = 1, st_order = 2, st_time = 3, st_length = 4, st_surface = 5, st_one = 6,
    st_userdata1 = 7, st_userdata2 = 8, st_userdata3 = 9, st_parenttype = 10,
    st_lb = 11, st_la = 12, st_nob = 13, st_r = 14, st_theta = 15, st_rlt = 16,
    st_meanln = 17, st_stdln = 18}; ///< @see RootSystem::getScalar
  static const std::vector<std::string> scalarTypeNames; ///< the corresponding names

  Plant() { initRTP(); };
  virtual ~Plant();

  // Parameter input output
  void setRootTypeParameter(RootTypeParameter p) { rtparam.at(p.type-1) = p; } ///< set the root type parameter to the index type-1
  RootTypeParameter* getRootTypeParameter(int type) { return &rtparam.at(type-1); } ///< Returns the i-th root parameter set (i=1..n)
  void setRootSystemParameter(const RootSystemParameter& rsp) { rsparam = rsp; } ///< sets the root system parameters
  RootSystemParameter* getRootSystemParameter() { return &rsparam; } ///< gets the root system parameters

  void openFile(std::string filename, std::string subdir="modelparameter/"); ///< Reads root paramter and plant parameter
  int readParameters(std::istream & cin); ///< Reads root parameters from an input stream
  void writeParameters(std::ostream & os) const; ///< Writes root parameters

  // Simulation
  void setGeometry(SignedDistanceFunction* geom) { geometry = geom; } ///< optionally, sets a confining geometry (call before RootSystem::initialize())
  void setSoil(SoilProperty* soil_) { soil = soil_; } ///< optionally sets a soil for hydro tropism (call before RootSystem::initialize())
  void reset(); ///< resets the root class, keeps the root type parameters
  void initialize(int basal=4, int shootborne=5); ///< creates the base roots, call before simulation and after setting the plant and root parameters
  void simulate(double dt, bool silence = false); ///< simulates root system growth for time span dt
  void simulate(); ///< simulates root system growth for the time defined in the root system parameters
  double getSimTime() const { return simtime; } ///< returns the current simulation time

  // call back functions (todo simplify)
  virtual Organ* createRoot(int lt, Vector3d  h, double delay, Organ* parent, double pbl, int pni);
  ///< Creates a new lateral root, overwrite or change this method to use more specialized root classes
  virtual TropismFunction* createTropismFunction(int tt, double N, double sigma);
  ///< Creates the tropisms, overwrite or change this method to add more tropisms TODO a vector<tropism*> might be easier to use
  virtual GrowthFunction* createGrowthFunction(int gft);
  ///< Creates the growth function per root type, overwrite or change this method to add more tropisms

  // Analysis of simulation results
  int getNumberOfNodes() const { return nid+1; } ///< Number of nodes of the root system
  int getNumberOfSegments() const { return nid+1-baseRoots.size(); } ///< Number of segments of the root system (the number of nodes-1 for tap root systems)
  std::vector<Organ*> getRoots() const; ///< Represents the root system as sequential vector of roots and buffers the result
  std::vector<Organ*> getBaseRoots() const { return baseRoots; } ///< Base roots are tap root, basal roots, and shoot borne roots
  std::vector<Vector3d> getNodes() const; ///< Copies all root system nodes into a vector
  std::vector<std::vector<Vector3d>> getPolylines() const; ///< Copies the nodes of each root into a vector return all resulting vectors
  std::vector<Vector2i> getSegments() const; ///< Copies all segments indices into a vector
  std::vector<Organ*> getSegmentsOrigin() const; ///< Copies a pointer to the root containing the segment
  std::vector<double> getNETimes() const; ///< Copies all node emergence times into a vector
  std::vector<std::vector<double>> getPolylinesNET() const; ///< Copies the node emergence times of each root into a vector and returns all resulting vectors
  std::vector<double> getScalar(int stype=Plant::st_length) const; ///< Copies a scalar root parameter that is constant per root to a vector
  std::vector<int> getRootTips() const; ///< Node indices of the root tips
  std::vector<int> getRootBases() const; ///< Node indices of the root bases

  // Dynamic information what happened last time step
  int getNumberOfNewNodes() { return getNumberOfNodes()-old_non; } ///< returns the number of new nodes, which is exactly the same number as new segments
  int getNumberOfNewRoots() { return getRoots().size() -old_nor; }  ///< returns the number of new roots
  std::vector<int> getNodeUpdateIndices(); // todo test and comment
  std::vector<Vector3d> getUpdatedNodes(); // to replace to the old node vector
  std::vector<Vector3d> getNewNodes(); // to dynamically add to the old node vector
  std::vector<Vector2i> getNewSegments(); // to dynamically add to the list of segments
  // restore(); ///< restore old time step

  // Output Simulation results
  void write(std::string name) const; /// writes simulation results (type is determined from file extension in name)
  void writeRSML(std::ostream & os) const; ///< writes current simulation results as RSML
  void writeVTP(std::ostream & os) const; ///< writes current simulation results as VTP (VTK polydata file)
  void writeGeometry(std::ostream & os) const; ///< writes the current confining geometry (e.g. a plant container) as paraview python script

  std::string toString() const; ///< infos about current root system state (for debugging)

  // random stuff
  void setSeed(double seed); ///< help fate (sets the seed of all random generators)
  double rand() { return UD(gen); } ///< Uniformly distributed random number (0,1)
  double randn() { return ND(gen); } ///< Normally distributed random number (0,1)

  int rsmlReduction = 5; ///< only each n-th node is written to the rsml file (to coarsely adjust axial resolution for output)

private:

  RootSystemParameter rsparam; ///< Plant parameter
  std::vector<RootTypeParameter> rtparam; ///< Parameter set for each root type
  std::vector<Organ*> baseRoots;  ///< Base roots of the root system
  std::vector<GrowthFunction*> gf; ///< Growth function per root type
  std::vector<TropismFunction*> tf;  ///< Tropism per root type
  SignedDistanceFunction* geometry = new SignedDistanceFunction(); ///< Confining geometry (unconfined by default)
  SoilProperty* soil = nullptr; ///< callback for hydro, or chemo tropism (needs to set before initialize()) TODO should be a part of tf, or rtparam

  double simtime = 0;
  int rid = -1; // unique root id counter
  int nid = -1; // unique root id counter

  int old_non=0;
  int old_nor=0;
  mutable std::vector<Organ*> roots = std::vector<Organ*>(); // buffer for getRoots()

  const int maxtypes = 100;
  void initRTP(); // default values for rtparam vector

  void writeRSMLMeta(std::ostream & os) const;
  void writeRSMLPlant(std::ostream & os) const;

  int getRootIndex() { rid++; return rid; } ///< returns next unique root id, called by the constructor of Root
  int getNodeIndex() { nid++; return nid; } ///< returns next unique node id, called by Root::addNode()

  std::mt19937 gen = std::mt19937(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> UD = std::uniform_real_distribution<double>(0,1);
  std::normal_distribution<double> ND = std::normal_distribution<double>(0,1);

};

#endif /* ROOTSYSTEM_H_ */
