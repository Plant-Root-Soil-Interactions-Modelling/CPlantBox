#ifndef SEED_H_
#define SEED_H_


#include "Organ.h"
#include "Root.h"

class Plant;
class Root;

/**
 * Seed
 *
 * the main organ of the plant,
 * simulate calls the simulate method of the stem, and base roots
 *
 */
class Seed : public Organ
{

public:

	Seed(Plant* plant);
	virtual ~Seed() { };

	virtual int organType() const override;

	virtual void initialize();
	virtual void simulate(double dt, bool silence = false) override;

	virtual std::string toString() const override;

	const int basalType = 4;
};


#endif /* Seed_H_ */
