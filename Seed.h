#ifndef SEED_H_
#define SEED_H_


#include "Organ.h"




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

  virtual void setRelativeOrigin(const Vector3d& o) override { seed_pos = o; };
  ///< the relative position within the parent organ

  virtual void initialize();

  virtual std::string toString() const override;

  const int basalType = 4;
  Vector3d seed_pos;
};


#endif /* Seed_H_ */
