
#ifndef LEAF_H_
#define LEAF_H_

#include <iostream>
#include <assert.h>

#include "Organ.h"
#include "mymath.h"
#include "sdf.h"
#include "tropism.h"
#include "LeafGrowth.h"
#include "leafparameter.h"
#include "Organism.h"


namespace CPlantBox {


class Plant;

static int leafphytomerID[10]= {0};

/**
 * Stem
 *
 * Describes a single stem, by a vector of nodes representing the stem.
 * The method simulate() creates new nodes of this stem, and lateral stems in the stem's branching zone.
 *
 */
class Leaf : public Organ
{

public:

	Leaf(Organism* plant, Organ* parent, int subtype, double delay, Vector3d ilheading, int pni, double pbl); ///< typically called by constructor of Plant::Plant, or Leaf::createLaterals()
	virtual ~Leaf() { }; // base class constructor is called automatically in c++

	virtual int organType() const override;

	/* simulation */
	virtual void simulate(double dt, bool silence = false) override; ///< stem growth for a time span of \param dt

	/* get results */
	virtual double getParameter(std::string name) const override; ///< returns an organ parameter of Plant::ScalarType

	/* exact from analytical equations */
	double getCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
	double LeafGetLength(double age); ///< analytical length of the stem
	double LeafGetAge(double length); ///< analytical age of the stem

	/* abbreviations */
    LeafRandomParameter* getLeafRandomParameter() const;  ///< root type parameter of this root
    const LeafSpecificParameter* param() const; ///< root parameter
	double dx() const; ///< returns the axial resolution
	std::string name() const;
	//Vector3d relHeading() const; //< relative heading of the Leaf tip
	//Vector3d absHeading() const; //< absolute heading of the Leaf tip
	Vector3d initialLeafHeading;
	Vector3d heading() const; /// current heading of the root tip
	/* IO */
//	void writeRSML(std::ostream & cout, std::string indent) const; ///< writes a RSML stem tag
	std::string toString() const;
    Vector3d iHeading; ///< the initial heading of the root, when it was created
    double parentBaseLength; ///< length [cm]
    int parentNI; ///< parent node index
	/* nodes */
	void addNode(Vector3d n, double t); //< adds a node to the stem

	/* parameters that are given per stem that are constant*/
	int pni; ///< parent node index
	double pbl; ///< parent base length [cm]

	const double smallDx = 1e-6; ///< threshold value, smaller segments will be skipped (otherwise stem tip direction can become NaN)
	Vector3d initialHeading;///< a heading downward

//	virtual void setRelativeOrigin(const Vector3d& o) override { this->o = o; }
//	virtual Vector3d getRelativeOrigin() const override { return o;  }
//	virtual void setRelativeHeading(const Matrix3d& m) override { this->A = m; }
//	virtual Matrix3d getRelativeHeading() const override { return A; }

	void createLateral(bool silence); ///< creates a new lateral, called by Leaf::simulate()
double parent_base_length; ///< length [cm]
	int parent_ni; ///< parent node index

	Vector3d o;
	Matrix3d A; // relative heading
	    	void minusPhytomerId(int subtype) { leafphytomerID[subtype]--;  }
	int getleafphytomerID(int subtype) { return leafphytomerID[subtype]; }
	void addleafphytomerID(int subtype) { leafphytomerID[subtype]++;  }

protected:

	void createSegments(double l, bool silence); ///< creates segments of length l, called by stem::simulate()
    virtual Vector3d getIncrement(const Vector3d& p, double sdx); ///< called by createSegments, to determine growth direction
	int old_non = 0;

};

} // namespace CPlantBox

#endif
