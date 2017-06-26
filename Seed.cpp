#include "Seed.h"


Seed::Seed(Plant* plant) :Organ(plant, nullptr, 0, 0)
{
}

/**
 * returns the type of the organ
 */
int Seed::organType()
{
	return Plant::ot_seed;
}

/**
 *
 */
void Seed::simulate(double dt, bool silence)
{
	// upper part
	if (shoot!=nullptr) {
		shoot->simulate(dt);
	}
	// lower part
	for (auto& r : roots)  {
		r->simulate(dt);
	}
}


/**
 * Quick info about the object for debugging
 */
std::string Seed::toString() const
{
  std::stringstream str;
  str << "Seed #"<< id <<": type "<< param->type << ", length: "<< length << ", age: " << age;
  return str.str();
}
