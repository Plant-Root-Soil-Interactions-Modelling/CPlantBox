#ifndef SEED_H_
#define SEED_H_


#include "Organ.h"
#include "Root.h"

class Plant;
class Root;

/**
 * Root
 *
 * Describes a single root, by a vector of nodes representing the root.
 * The method simulate() creates new nodes of this root, and lateral roots in the root's branching zone.
 *
 */
class Seed : public Organ
{

public:

	Seed(Plant* plant);

	virtual int organType();

	virtual void initialize() { }; ///< create call backs from the parameters set TODO
	virtual void simulate(double dt, bool silence = false);

	virtual std::string toString() const;

	std::vector<Root*> roots; // tap root, and basal roots

	Organ* shoot = nullptr; // todo will be of type shoot

};


#endif /* Seed_H_ */
