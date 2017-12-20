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
* stem
*
* Describes a single stem, by a vector of nodes representing the stem.
* The method simulate() creates new nodes of this stem, and lateral stems in the stem's branching zone.
*
*/
class Stem : public Organ
{
friend class Leaf;



public:


  Stem(Plant* plant, Organ* parent, int type, double delay, Vector3d isheading, int pni, double pbl); ///< typically called by constructor of Plant::Plant, or Stem::createLaterals()
  virtual ~Stem() { }; // base class constructor is called automatically in c++

  virtual int organType() const override { return Organ::ot_stem; };

  /* simulation */
  virtual void simulate(double dt, bool silence = false) override; ///< stem growth for a time span of \param dt

  /* get results */
  virtual double getScalar(std::string name) const override; ///< returns an organ parameter of Plant::ScalarType

  /* exact from analytical equations */
  double getCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
  double StemgetLength(double age); ///< analytical length of the stem
  double StemgetAge(double length); ///< analytical age of the stem

  /* abbreviations */
  StemParameter* sParam() const { return (StemParameter*)param;  } ///< type cast
  StemTypeParameter* stParam() const; // type cast
  double dx() const; ///< returns the axial resolution
  Vector3d heading() const; /// current heading of the stem tip

  /* IO */
  void writeRSML(std::ostream & cout, std::string indent) const; ///< writes a RSML stem tag
  std::string toString() const;

  /* nodes */
  void addNode(Vector3d n, double t); //< adds a node to the stem

  /* parameters that are given per stem that are constant*/
  int pni; ///< parent node index
  double pbl; ///< length [cm]

  const double smallDx = 1e-6; ///< threshold value, smaller segments will be skipped (otherwise stem tip direction can become NaN)
  Vector3d initialHeading;/// a heading downward
  Vector3d initialStemHeading;/// heading upward

  int LeafID = 0; //declare leaf id.


protected:

  void createSegments(double l, bool silence); ///< creates segments of length l, called by stem::simulate()
  void createLateral(bool silence); ///< creates a new lateral, called by Stem::simulate()
  void LeafGrow(bool silence);
  void ShootBorneRootGrow(bool silence);

  int old_non = 0;

};

#endif
