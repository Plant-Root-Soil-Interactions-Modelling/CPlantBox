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

    Seed(Organism* plant) :Organ(plant, nullptr, Organism::ot_seed, 0, 0.) { };
    virtual ~Seed() { }; // everything is done in @Organ

    virtual int organType() const override { return Organism::ot_seed; }

    void initialize();

    SeedSpecificParameter* param() const { return (SeedSpecificParameter*)param_; }

    virtual std::string toString() const override;

    int getNumberOfRootCrowns() const { return numberOfRootCrowns; }
    std::vector<Organ*>& baseOrgans() { return children; }
    std::vector<Organ*> copyBaseOrgans(); ///< shallow copy of the childs

    // default positions (unused) (TODO) make nicer
    int tapRootType =1; // todo
    int basalType = 4;
    int shootborneType = 5;
    int mainStemType = 1; // todo
    int tillerType = 4; // todo

    virtual Organ* createRoot(Organism* rs, int type, Vector3d pheading, double delay, Organ* parent, double pbl, int pni); // overwrite if you want to change the types
    virtual Organ* createStem(Organism* rs, int type, Vector3d pheading, double delay, Organ* parent, double pbl, int pni); // overwrite if you want to change the types

protected:

    int numberOfRootCrowns = 0;
    int getParamSubType(int organtype, std::string str);

};

}

#endif /* Seed_H_ */
