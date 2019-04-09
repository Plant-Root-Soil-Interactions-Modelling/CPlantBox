#ifndef SEED_H_
#define SEED_H_


#include "Organ.h"

namespace CPlantBox {


class Plant;

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

	virtual int organType() const override { return Organ::ot_seed; };

	virtual Vector3d getRelativeOrigin() const override { return seed_pos; };
	///< the relative position within the parent organ
	Vector3d getseedPos(SeedParameter* sparam) const { return sparam->seedPos; };
	virtual void setRelativeOrigin(const Vector3d& o) override { seed_pos = o; };
	///< the relative position within the parent organ
	virtual SeedParameter* initializeparam();
	virtual void initialize(SeedParameter* sparam);
    	virtual void setRelativeHeading(const Matrix3d& m) override { this->A = m; };
	virtual std::string toString() const override;

	const int basalType = 4;
	const int tillerType = 4;
	Vector3d seed_pos;
     	Matrix3d A; // relative heading
};


} // namespace CPlantBox

#endif /* Seed_H_ */
