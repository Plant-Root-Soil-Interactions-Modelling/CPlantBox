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
#include "external/tinyxml2/tinyxml2.h"

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

  enum TropismTypes { tt_plagio = 0, tt_gravi = 1, tt_exo = 2, tt_hydro = 3, tt_antigravi = 4, tt_twist = 5};  ///< plant tropism types
  enum GrowthFunctionTypes { gft_negexp = 1, gft_linear = 2 }; // plant growth function

  Plant();
  virtual ~Plant() { };

  std::shared_ptr<Organism> copy() override; ///< deep copies the organism

  /* organs */
  std::shared_ptr<Seed> getSeed(); ///< the plant seed

  /* parameters */
  void initializeReader() override; ///< initializes XML reader
  void readParameters(std::string name, std::string  basetag = "plant") override { this->initializeReader(); Organism::readParameters(name, basetag); };
  void openXML(std::string name) { readParameters(name); } // old name

  /* Simulation */
  void setGeometry(std::shared_ptr<SignedDistanceFunction> geom) { geometry = geom; }  ///< optionally, sets a confining geometry (call before Plant::initialize())
  void reset(); ///< resets the plant class, keeps the organ type parameters
  void initialize(bool verbose = true) override; ///< creates the base roots, call before simulation and after setting the plant and root parameters
  void setTropism(std::shared_ptr<Tropism> tf, int organType, int subType = -1); ///< todo docme
  void simulate(); ///< simulates root system growth for the time defined in the root system parameters

  /* call back function creation */
  void initCallbacks(); ///< sets up callback functions for tropisms and growth functions, called by initialize()
  virtual std::shared_ptr<Tropism> createTropismFunction(int tt, double N, double sigma); ///< Creates the tropisms, overwrite or change this method to add more tropisms
  virtual std::shared_ptr<GrowthFunction> createGrowthFunction(int gft); ///< Creates the growth function per root type, overwrite or change this method to add more tropisms

  void write(std::string name) const; /// writes simulation results (type is determined from file extension in name)
  std::string toString() const override;
  void writeVTP(int otype, std::ostream & os) const;

  std::vector<int> leafphytomerID = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

protected:

  std::shared_ptr<SignedDistanceFunction> geometry = std::make_shared<SignedDistanceFunction>();  ///< Confining geometry (unconfined by default)
  std::shared_ptr<SoilLookUp> soil; ///< callback for hydro, or chemo tropism (needs to set before initialize()) TODO should be a part of tf, or rtparam

};

} // namespace CPlantBox

#endif /* ROOTSYSTEM_H_ */
