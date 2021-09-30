#ifndef LEAF_H_
#define LEAF_H_

#include "Organ.h"
#include "Organism.h"
#include "leafparameter.h"

#include <iostream>
#include <assert.h>

namespace CPlantBox {

class Plant;

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

    Leaf(int id,  std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
        Matrix3d iHeading, int pni, bool moved, int oldNON = 0);
	Leaf(std::shared_ptr<Organism> plant, int type, Matrix3d iHeading, double delay, std::shared_ptr<Organ> parent, int pni); ///< used within simulation
	virtual ~Leaf() { };

	std::shared_ptr<Organ> copy(std::shared_ptr<Organism> plant) override;   ///< deep copies the root tree

	int organType() const override { return Organism::ot_leaf; } ///< returns the organs type

	void simulate(double dt, bool silence = false) override; ///< stem growth for a time span of \param dt

    Vector3d getNode(int i) const override { return nodes.at(i); } ///< i-th node of the organ

	double getParameter(std::string name) const override; ///< returns an organ parameter of Plant::ScalarType

	/* leaf vizualisation */
    double leafLength() const { return std::max(getLength(false)-param()->lb, 0.); /* represents the leaf base*/ }; ///< leaf surface length [cm]
    double leafCenter() const { return std::max(getLength(false)-param()->la-param()->lb, 0.); }; ///< center of the radial parametrisation
    double leafArea() ; ///< returns the leaf surface area, zero if there are lateral-leafs [cm2]
	std::vector<double> getLeafVisX(int i);
	std::vector<Vector3d> getLeafVis(int i); // per node

    std::string toString() const override;

	/* exact from analytical equations */
	double calcCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
	double calcLength(double age); ///< analytical length of the stem
	double calcAge(double length); ///< analytical age of the stem

	/* abbreviations */
	std::shared_ptr<LeafRandomParameter> getLeafRandomParameter() const;  ///< root type parameter of this root
	std::shared_ptr<const LeafSpecificParameter> param() const; ///< root parameter

	/* useful */
    Vector3d heading(int n)  const override; ///< current (absolute) heading of the organs at node n
    Vector3d heading() const override {return heading( -1 ); } 
	
    void rel2abs() override;
	void abs2rel() override;
	Vector3d getiHeading() const;
	bool hasMoved() const override { return true; }; ///< have any nodes moved during the last simulate call

protected:

	Vector3d partialIHeading;

    int getleafphytomerID(int subtype);
    void minusPhytomerId(int subtype);
    void addleafphytomerID(int subtype);

    void createLateral(bool silence); ///< creates a new lateral, called by Leaf::simulate()

    Vector3d getIncrement(const Vector3d& p, double sdx, int n=-1); ///< called by createSegments, to determine growth direction
	void createSegments(double l, bool silence); ///< creates segments of length l, called by stem::simulate()

    bool nodeLeafVis(double l); ///<  leaf base (false), branched leaf (false), or leaf surface area (true)
	std::vector<double> getLeafVisX_(double l);
	bool ageDependentTropism = false;
    bool firstCall = true;
};

} // namespace CPlantBox

#endif
