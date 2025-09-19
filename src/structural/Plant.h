#ifndef PLANT_H_
#define PLANT_H_

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <chrono>
#include <random>
#include <numeric>
#include <istream>

#include "Organ.h"
#include "Root.h"
#include "Seed.h"
#include "Stem.h"
#include "Leaf.h"

#include "soil.h"
#include "tropism.h"
#include "growth.h"
#include "tinyxml2.h"

namespace CPlantBox {

/**
 * Plant
 *
 * This class manages all model parameter, the simulation,
 * and stores the seed of the plant,
 * and offers utility functions for post processing
 *
 */
class Plant :public Organism
{
public:

  enum TropismTypes { tt_plagio = 0, tt_gravi = 1, tt_exo = 2, tt_hydro = 3, tt_antigravi = 4, tt_twist = 5,  tt_antigravi2gravi = 6};  ///< plant tropism types
  enum GrowthFunctionTypes { gft_negexp = 1, gft_linear = 2 , gft_CWLim = 3 }; // plant growth function

  Plant(unsigned int seednum  = 0.);
  virtual ~Plant() { };

  std::shared_ptr<Organism> copy() override; ///< deep copies the organism

  /* parameters */
  void initializeReader() override; ///< initializes XML reader
  void readParameters(std::string name, std::string  basetag = "plant", bool fromFile = true, bool verbose = true) override {this->initializeReader(); Organism::readParameters(name, basetag, fromFile, verbose); };
  void openXML(std::string name) { readParameters(name); } // old name

  /* Simulation */
  void setSoil(std::shared_ptr<SoilLookUp> soil_) { soil = soil_; } ///< optionally sets a soil for hydro tropism (call before Plant::initialize())
  void reset(); ///< resets the plant class, keeps the organ type parameters
  virtual void initializeLB(bool verbose = true); ///< creates the base roots (length based lateral emergence times), call before simulation and after setting plant and root parameters
  virtual void initializeDB(bool verbose = true); ///< creates the base roots (delay based lateral emergence times), call before simulation and after setting plant and root parameters
  void initialize(bool verbose = true) override { initializeLB(verbose); };
  void setTropism(std::shared_ptr<Tropism> tf, int organType, int subType = -1); ///< todo docme
  void simulate(); ///< simulates root system growth for the time defined in the root system parameters
  void simulate(double dt, bool verbose = false) override;
  void simulate(double dt, double maxinc, std::shared_ptr<ProportionalElongation> se, bool verbose = true); ///< simulates the plant with a maximal elongation
  void simulateLimited(double dt, double max_, std::string paramName, std::vector<double> scales, std::shared_ptr<ProportionalElongation> se, bool verbose);  ///< simulates plant with limited costs

  /* call back function creation */
  void initCallbacks(); ///< sets up callback functions for tropisms and growth functions, called by initialize()
  std::shared_ptr<Tropism> createTropismFunction(int tt, double N, double sigma, double Tage = 0.); ///< Creates the tropisms, overwrite or change this method to add more tropisms
  virtual std::shared_ptr<GrowthFunction> createGrowthFunction(int gft); ///< Creates the growth function per root type, overwrite or change this method to add more tropisms

  std::string toString() const override;

  std::vector<int> leafphytomerID = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  void abs2rel();
  void rel2abs();

protected:

  std::shared_ptr<SoilLookUp> soil; ///< callback for hydro, or chemo tropism (needs to set before initialize()) TODO should be a part of tf, or rtparam

  void initialize_(bool verbose = true); // called by initializeLB, and initializeDB
  double weightedSum(std::string paramName, std::vector<double> scales) const; // weighted sum per organ type (used by simulate_limited)

};

} // namespace CPlantBox

#endif /* PLANT_H_ */
