// *** ADDED BY HEADER FIXUP ***
#include <istream>
// *** END ***
#ifndef PLANT_H_
#define PLANT_H_

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <chrono>
#include <random>
#include <numeric>

#include "ModelParameter.h"
#include "Organ.h"
#include "Root.h"
#include "Seed.h"
#include "Stem.h"

#include "soil.h"
#include "RootTropism.h"
#include "RootGrowth.h"
#include "StemGrowth.h"
#include "StemTropism.h"

class Seed;
class SeedParameter;

/**
 * Plant
 *
 * This class manages all model parameter and the simulation,
 * stores the seed of the plant,
 * and offers utility functions for post processing
 *
 */
class Plant
{

public:

  enum ScalarTypes { st_id, st_otype, st_subtype, st_alive, st_active, st_age, st_length, st_one, st_order, st_parenttype, st_time, // organ level
	 st_lb, st_la, st_r, st_radius, st_theta, st_rlt, st_meanln, st_stdln , st_nob, st_surface, // root level
	 st_userdata1, st_userdata2, st_userdata3 // analyser
  }; ///< @see RootSystem::getScalar
  static const std::vector<std::string> scalarTypeNames; ///< the corresponding names

  Plant();
  virtual ~Plant();

  /* Parameter input output */
  void setOrganTypeParameter(OrganTypeParameter*  p );///< set the organ type parameter TODO change p
  OrganTypeParameter* getOrganTypeParameter(int otype, int subtype) const;
  void openFile(std::string filename, std::string subdir="modelparameter/"); ///< Reads root paramter and plant parameter
  int readRootParameters(std::istream & cin); ///< Reads root parameters from an input stream
  int readStemParameters(std::istream & cin); ///< Reads stem parameters from an input stream
  void writeParameters(std::ostream & os) const; ///< Writes root parameters

  /* Simulation */
  void setGeometry(SignedDistanceFunction* geom); ///< optionally, sets a confining geometry (call before RootSystem::initialize())
  void reset(); ///< resets the root class, keeps the root type parameters
  void initialize(int basal=4, int shootborne=5); ///< creates the base roots, call before simulation and after setting the plant and root parameters
  void simulate(double dt, bool silence = false); ///< simulates root system growth for time span dt
  void simulate(); ///< simulates root system growth for the time defined in the root system parameters
  double getSimTime() const { return simtime; } ///< returns the current simulation time

  /* Analysis of simulation results */
  // Organs
  int getNumberOfNodes() const { return nid+1; } ///< Number of nodes of the root system
  int getNumberOfSegments() const {
	  getOrgans(Organ::ot_organ);
	  int c = 0;
	  for (const auto& o : organs) {
		  c += (o->r_nodes.size()-1);
	  }
	  return c;
  } ///< todo -baseRoots.size() Number of segments of the root system (the number of nodes-1 for tap root systems)
  std::vector<Organ*> getOrgans(unsigned int otype) const; ///< Represents the root system as sequential vector of roots and buffers the result
  std::vector<Vector3d> getNodes() const; ///< Copies all root system nodes into a vector
  std::vector<std::vector<Vector3d> > getPolylines(unsigned int otype=Organ::ot_organ) const; ///< Copies the nodes of each root into a vector return all resulting vectors
  std::vector<Vector2i> getSegments(unsigned int otype=Organ::ot_organ) const; ///< Copies all segments indices into a vector
  std::vector<Organ*> getSegmentsOrigin(unsigned int otype=Organ::ot_organ) const; ///< Copies a pointer to the root containing the segment
  std::vector<double> getNETimes() const; ///< Copies all node emergence times into a vector
  std::vector<std::vector<double> > getPolylinesNET(unsigned int otype=Organ::ot_organ) const; ///< Copies the node emergence times of each root into a vector and returns all resulting vectors
  std::vector<double> getScalar(unsigned int otype=Organ::ot_organ, int stype=Plant::st_length) const; ///< Copies a scalar root parameter that is constant per root to a vector

  // std::vector<Organ*> getBaseRoots() const { return baseRoots; } ///< Base roots are tap root, basal roots, and shoot borne roots
  std::vector<int> getRootTips() const; ///< Node indices of the root tips
  std::vector<int> getRootBases() const; ///< Node indices of the root bases

  // Output Simulation results
  void write(std::string name, int otype = Organ::ot_organ) const; /// writes simulation results (type is determined from file extension in name)
  void writeRSML(std::ostream & os) const; ///< writes current simulation results as RSML
  void writeVTP(int otype, std::ostream & os) const; ///< writes current simulation results as VTP (VTK polydata file)
  void writeGeometry(std::ostream & os) const; ///< writes the current confining geometry (e.g. a plant container) as paraview python script

  std::string toString() const; ///< infos about current root system state (for debugging)

  // random stuff
  void setSeed(unsigned int seed) const ; ///< help fate (sets the seed of all random generators)
  double rand() const { return UD(gen); } ///< Uniformly distributed random number (0,1)
  double randn() const { return ND(gen); } ///< Normally distributed random number (0,1)

  int getOrganIndex() { rid++; return rid; } ///< returns next unique root id, called by the constructor of Root
  int getNodeIndex() { nid++; return nid; } ///< returns next unique node id, called by Root::addNode()
  int rsmlReduction = 5; ///< only each n-th node is written to the rsml file (to coarsely adjust axial resolution for output)

protected:

  std::vector <std::vector<OrganTypeParameter*> > organParam; ///< Parameter set for each root type

  Seed* seed;

  SignedDistanceFunction* geometry = new SignedDistanceFunction(); ///< Confining geometry (unconfined by default)

  double simtime = 0;
  int rid = -1; // unique root id counter
  int nid = -1; // unique root id counter

  int old_non=0;
  int old_nor=0;
  mutable unsigned int organs_type = -1; // type of buffered organs
  mutable std::vector<Organ*> organs = std::vector<Organ*>(); // buffer

  const int maxtypes = 20;
  const int maxorgans = 10;
  void initOTP(); // default values for organ type parameters

  void writeRSMLMeta(std::ostream & os) const;
  void writeRSMLPlant(std::ostream & os) const;

  mutable std::mt19937 gen = std::mt19937(std::chrono::system_clock::now().time_since_epoch().count());
  mutable std::uniform_real_distribution<double> UD = std::uniform_real_distribution<double>(0,1);
  mutable std::normal_distribution<double> ND = std::normal_distribution<double>(0,1);

};


































#endif /* ROOTSYSTEM_H_ */
