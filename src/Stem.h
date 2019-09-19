#ifndef STEM_H_
#define STEM_H_

#include <iostream>
#include <assert.h>

#include "Organ.h"
#include "mymath.h"
#include "StemTropism.h"
#include "StemGrowth.h"
#include "stemparameter.h"
#include "Organism.h"

namespace CRootBox {

static int phytomerId[10]= {0};
/**
 * Stem
 *
 * Describes a single stem, by a vector of nodes representing the stem.
 * The method simulate() cr`eates new nodes of this stem, and lateral stems in the stem's branching zone.
 *
 */
class Stem : public Organ
{
	friend class Leaf;



public:

	Stem(Organism* plant, Organ* parent, int subtype, double delay, Vector3d isheading, int pni, double pbl); ///< typically called by constructor of Plant::Plant, or Stem::createLaterals()
	virtual ~Stem() { }; // base class constructor is called automatically in c++

	virtual int organType() const override { return Organism::ot_stem; };
    Organ* copy(Organism* rs) override;  ///< deep copies the root tree

	/* simulation */
	virtual void simulate(double dt, bool silence = false) override; ///< stem growth for a time span of \param dt

	/* get results */
	double getParameter(std::string name) const override; ///< returns an organ parameter

	/* exact from analytical equations */
	double getCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
	double StemGetLength(double age); ///< analytical length of the stem
	double StemGetAge(double length); ///< analytical age of the stem
	/* abbreviations */
StemRandomParameter* getStemRandomParameter() const;  ///< root type parameter of this root
    const StemSpecificParameter* param() const; ///< root parameter

	double dx() const; ///< returns the axial resolution
	//Vector3d relHeading() const; //< relative heading of the stem tip
	//Vector3d absHeading() const; //< absolute heading of the stem tip

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
	double pbl; ///< length [cm]

	const double smallDx = 1e-6; ///< threshold value, smaller segments will be skipped (otherwise stem tip direction can become NaN)

	int StemID = 0; //declare Stem id.
//	virtual void setRelativeOrigin(const Vector3d& o) override { this->o = o; }
//	virtual Vector3d getRelativeOrigin() const override { return o;  }
//	virtual void setRelativeHeading(const Matrix3d& m) override { this->A = m; }
//	virtual Matrix3d getRelativeHeading() const override { return A; }
	Vector3d initialStemHeading;
	Vector3d heading() const; /// current heading of the root tip
	Vector3d o;
	Matrix3d A; // relative heading
double parent_base_length; ///< length [cm]
	int parent_ni; ///< parent node index

    	void minusPhytomerId(int subtype) { phytomerId[subtype]--;  }
	int getphytomerId(int subtype) { return phytomerId[subtype]; }
	void addPhytomerId(int subtype) { phytomerId[subtype]++;  }
protected:


	void createSegments(double l, bool silence); ///< creates segments of length l, called by stem::simulate()
    virtual Vector3d getIncrement(const Vector3d& p, double sdx); ///< called by createSegments, to determine growth direction
	void createLateral(bool silence); ///< creates a new lateral, called by Stem::simulate()
	void LeafGrow(bool silence, Vector3d bud);
	void ShootBorneRootGrow(bool silence);
	int old_non = 0; // relative origin

};

} // namespace CPlantBox

#endif
