#ifndef SEED_H_
#define SEED_H_

#include "Organ.h"
#include "Organism.h"
#include "seedparameter.h"

namespace CPlantBox {

/**
 * Seed
 *
 * Creates one node, which represents the seed.
 *
 * Seed::initialize() creates the initial organs, i.e. taproot, optional: basal, shoot borne root, if plant stem, optional: tillers
 * Seed::simulate() calls the simulate method of the created organs, i.e. taproot, basal roots, shoot borne roots, stem, tillers
 *
 */
class Seed : public Organ {

  public:
    Seed(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length, bool moved = false,
         int oldNON = 0);                  ///< creates everything from scratch
    Seed(std::shared_ptr<Organism> plant); ///< used within simulation
    virtual ~Seed() {};

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> rs) override; ///< deep copies the seed

    virtual int organType() const override { return Organism::ot_seed; }

    virtual std::string toString() const override;

    void initialize(
        bool verbose = true); ///< initialization creates the initial organs, i.e. taproot, optional: basal, shoot borne root, if plant stem, optional: tillers

    std::shared_ptr<const SeedSpecificParameter> param() const { return std::static_pointer_cast<const SeedSpecificParameter>(param_); }

    std::vector<std::shared_ptr<Organ>> baseOrgans() { return children; } // created by initialize
    std::vector<std::shared_ptr<Organ>> copyBaseOrgans(std::shared_ptr<Organism> plant);

    virtual std::shared_ptr<Organ> createRoot(std::shared_ptr<Organism> plant, int type,
                                              double delay); ///< overwrite if you want to change class types used by intialize()
    virtual std::shared_ptr<Organ> createStem(std::shared_ptr<Organism> plant, int type,
                                              double delay); ///< overwrite if you want to change class types used by intialize()

    int tapType = 1; ///< defauts, may be adapted by intialize()
    int basalType = 4;
    int shootborneType = 5;
    int mainStemType = 1;
    int tillerType = 4;

    bool basalRoots() { return (param()->maxB > 0) && (param()->firstB < getMaxSimTime()); }     ///< returns true if the seed has basal roots
    bool shootBorneRoots() { return (param()->nC > 0) && (param()->firstSB < getMaxSimTime()); } ///< returns true if the seed has shoot borne roots
    bool tillers() { return (param()->maxTil > 0) && (param()->firstTil < getMaxSimTime()); }    ///< returns true if the seed has tillers

    bool basalRootDefined();      ///< is there a basal root definition in the parameter set? works only before initialization, otherwise the basal root type is
                                  ///< copied from to the tap root type
    bool shootBorneRootDefined(); ///< is there a shoot borne root definition in the parameter set? works only before initialization, otherwise the shoot borne
                                  ///< root type is copied from to the tap
    bool tillersDefined(); ///< is there a tiller definition in the parameter set? works only before initialization, otherwise the tiller type is copied from to
                           ///< the stem type

    int getNumberOfRootCrowns() const { return rootCrowns_; }            ///< returns the number of root crowns, after initialization (currently == 1)
    double getMaxSimTime() const { return param()->simtime; } ///< returns the recommended final simulation time

  protected:
    int rootCrowns_ = 0;
    int getParamSubType(int organtype, std::string str);
};

} // namespace CPlantBox

#endif /* Seed_H_ */
