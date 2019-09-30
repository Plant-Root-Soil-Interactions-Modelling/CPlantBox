#ifndef SEED_H_
#define SEED_H_

#include "Organ.h"
#include "Organism.h"
#include "seedparameter.h"

namespace CPlantBox {

/**
 * Seed
 *
 * simulate calls the simulate method of the stem, and base roots
 */
class Seed : public Organ
{
public:

    Seed(std::shared_ptr<Organism> plant) :Organ(plant, nullptr, Organism::ot_seed, 0, 0.) { };
    virtual ~Seed() { };

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> rs) override;  ///< deep copies the seed

    virtual int organType() const override { return Organism::ot_seed; }

    virtual std::string toString() const override;

    void initialize();

    std::shared_ptr<const SeedSpecificParameter> param() const { return std::static_pointer_cast<const SeedSpecificParameter>(param_); }

    int getNumberOfRootCrowns() const { return numberOfRootCrowns; } // for rootsystem initialisation
    std::vector<std::shared_ptr<Organ>> baseOrgans() { return children; } // created by initialize
    std::vector<std::shared_ptr<Organ>> copyBaseOrgans(); ///< shallow copy of the childs

    virtual std::shared_ptr<Organ> createRoot(std::shared_ptr<Organism> plant, int type, Vector3d heading, double delay); ///< overwrite if you want to change class types
    virtual std::shared_ptr<Organ> createStem(std::shared_ptr<Organism> plant, int type, Vector3d heading, double delay); ///< overwrite if you want to change class types

    // default positions (unused) (TODO) make nicer
    int tapRootType =1; // todo
    int basalType = 4;
    int shootborneType = 5;
    int mainStemType = 1; // todo
    int tillerType = 4; // todo

protected:

    int numberOfRootCrowns = 0;
    int getParamSubType(int organtype, std::string str);

};

}

#endif /* Seed_H_ */
