// *** ADDED BY HEADER FIXUP ***
#include <cassert>
#include <istream>
// *** END ***
#ifndef STEM_H_
#define STEM_H_

#include <iostream>
#include <assert.h>

#include "Organ.h"
#include "mymath.h"
#include "sdf.h"
#include "StemTropism.h"
#include "StemGrowth.h"
#include "ModelParameter.h"

class Plant;
class StemParameter;
class StemTypeParameter;

/**
* Root
*
* Describes a single root, by a vector of nodes representing the root.
* The method simulate() creates new nodes of this root, and lateral roots in the root's branching zone.
*
*/
class Stem : public Organ
{

public:





  Stem(Plant* plant, Organ* parent, int type, double delay, Vector3d isheading, int pni, double pbl); ///< typically called by constructor of RootSystem, or Root::createLaterals()
  virtual ~Stem() { }; // base class constructor is called automatically in c++

  virtual int organType() const override;

  /* simulation */
  virtual void simulate(double dt, bool silence = false) override; ///< root growth for a time span of \param dt

  /* get results */
  virtual double getScalar(int stype) const override; ///< returns an organ parameter of Plant::ScalarType

  /* exact from analytical equations */
  double getCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
  double StemgetLength(double age); ///< analytical length of the root
  double StemgetAge(double length); ///< analytical age of the root

  /* abbreviations */
  StemParameter* sParam() const { return (StemParameter*)stem_param;  } ///< type cast
  StemTypeParameter* stParam() const; // type cast
  double dx() const; ///< returns the axial resolution
  Vector3d heading() const; /// current heading of the root tip

  /* IO */
  void writeRSML(std::ostream & cout, std::string indent) const; ///< writes a RSML root tag
  std::string toString() const;

  /* nodes */
  void addNode(Vector3d n, double t); //< adds a node to the root

  /* parameters that are given per root that are constant*/
  int pni; ///< parent node index
  double pbl; ///< length [cm]

  const double smallDx = 1e-6; ///< threshold value, smaller segments will be skipped (otherwise root tip direction can become NaN)

  Vector3d initialHeading;
  Vector3d initialstemHeading;
protected:

  void createSegments(double l, bool silence); ///< creates segments of length l, called by Root::simulate()
  void createLateral(bool silence); ///< creates a new lateral, called by Root::simulate()

  int old_non = 0;

};

#endif
