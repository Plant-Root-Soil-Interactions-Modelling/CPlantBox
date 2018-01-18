
// *** ADDED BY HEADER FIXUP ***
#include <cassert>
#include <istream>
// *** END ***
#ifndef LEAF_H_
#define LEAF_H_

#include <iostream>
#include <assert.h>

#include "Organ.h"
#include "mymath.h"
#include "sdf.h"
#include "LeafTropism.h"
#include "LeafGrowth.h"
#include "ModelParameter.h"



class Plant;
class LeafParameter;
class LeafTypeParameter;

/**
* stem
*
* Describes a single stem, by a vector of nodes representing the stem.
* The method simulate() creates new nodes of this stem, and lateral stems in the stem's branching zone.
*
*/
class Leaf : public Organ
{

public:

  Leaf(Plant* plant, Organ* parent, int type, double delay, Vector3d ilheading, int pni, double pbl); ///< typically called by constructor of Plant::Plant, or Leaf::createLaterals()
  virtual ~Leaf() { }; // base class constructor is called automatically in c++

  virtual int organType() const override;

  /* simulation */
  virtual void simulate(double dt, bool silence = false) override; ///< stem growth for a time span of \param dt

  /* get results */
  virtual double getScalar(std::string name) const override; ///< returns an organ parameter of Plant::ScalarType

  /* exact from analytical equations */
  double getCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
  double LeafgetLength(double age); ///< analytical length of the stem
  double LeafgetAge(double length); ///< analytical age of the stem

  /* abbreviations */
  LeafParameter* lParam() const { return (LeafParameter*)param;  } ///< type cast
  LeafTypeParameter* ltParam() const; // type cast
  double dx() const; ///< returns the axial resolution
  Vector3d heading() const; /// current heading of the Leaf tip

  /* IO */
  void writeRSML(std::ostream & cout, std::string indent) const; ///< writes a RSML stem tag
  std::string toString() const;

  /* nodes */
  void addNode(Vector3d n, double t); //< adds a node to the stem

  /* parameters that are given per stem that are constant*/
  int pni; ///< parent node index
  double pbl; ///< parent base length [cm]

  const double smallDx = 1e-6; ///< threshold value, smaller segments will be skipped (otherwise stem tip direction can become NaN)
  Vector3d initialHeading;///< a heading downward
  Vector3d initialStemHeading;///< a heading upward
  Vector3d initialLeafHeading;///< leave heading direction


  void createLateral(bool silence); ///< creates a new lateral, called by Leaf::simulate()



protected:

  void createSegments(double l, bool silence); ///< creates segments of length l, called by stem::simulate()

  int old_non = 0;

};

#endif
