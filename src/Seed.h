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

    Seed(std::weak_ptr<Organism> plant) :Organ(plant, std::weak_ptr<Organ>(),Organism::ot_seed, 0, 0.) { };
    virtual ~Seed() { }; // everything is done in @Organ

    virtual int organType() const override { return Organism::ot_seed; }

    void initialize();

    std::shared_ptr<const SeedSpecificParameter> param() const { return std::static_pointer_cast<const SeedSpecificParameter>(param_); }

    virtual std::string toString() const override;

    int getNumberOfRootCrowns() const { return numberOfRootCrowns; } // for rootsystem initialisation
    std::vector<std::shared_ptr<Organ>>& baseOrgans() { return children; } // created by initialize
    std::vector<std::shared_ptr<Organ>> copyBaseOrgans(); ///< shallow copy of the childs

    // default positions (unused) (TODO) make nicer
    int tapRootType =1; // todo
    int basalType = 4;
    int shootborneType = 5;
    int mainStemType = 1; // todo
    int tillerType = 4; // todo

    virtual std::shared_ptr<Organ> createRoot(std::weak_ptr<Organism> plant, int type, Vector3d heading, double delay,
        std::weak_ptr<Organ> parent = std::weak_ptr<Organ>(), double pbl = 0., int pni = 0); // overwrite if you want to change the types
    virtual std::shared_ptr<Organ> createStem(std::weak_ptr<Organism> plant, int type, Vector3d heading, double delay,
        std::weak_ptr<Organ> parent = std::weak_ptr<Organ>(), double pbl = 0., int pni = 0); // overwrite if you want to change the types

protected:

    int numberOfRootCrowns = 0;
    int getParamSubType(int organtype, std::string str);

};

}

#endif /* Seed_H_ */
