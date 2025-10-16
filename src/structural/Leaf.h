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
 * Leaf
 *
 * Describes a single leaf, by a vector of nodes representing the leaf
 * The method simulate() creates new nodes of this leaf, and lateral leafs in the leaf's branching zone.
 *
 */
class Leaf : public Organ
{
public:

    Leaf(int id,  std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
        Vector3d partialIHeading_, int pni, bool moved = true, int oldNON = 0);
	Leaf(std::shared_ptr<Organism> plant, int type, double delay, std::shared_ptr<Organ> parent, int pni); ///< used within simulation
	virtual ~Leaf() { };

	std::shared_ptr<Organ> copy(std::shared_ptr<Organism> plant) override;   ///< deep copies the root tree

	int organType() const override { return Organism::ot_leaf; } ///< returns the organs type

	void simulate(double dt, bool silence = false) override; ///< leaf growth for a time span of \param dt

    Vector3d getNode(int i) const override { return nodes.at(i); } ///< i-th node of the organ

	double getParameter(std::string name) const override; ///< returns an organ parameter of Plant::ScalarType

	/* leaf vizualisation */
    double leafLength( bool realized = false) const { return std::max(getLength(realized)-param()->lb, 0.); /* represents the leaf base*/ }; ///< leaf surface length [cm]
    double leafCenter( bool realized = false) const { return std::max(getLength(realized)-param()->la-param()->lb, 0.); }; ///< center of the radial parametrisation
    double leafArea( bool realized = false, bool withPetiole = false) const; ///< returns the leaf surface area, zero if there are lateral-leafs [cm2]
	double leafAreaAtSeg(int localSegId, bool realized = false, bool withPetiole = false); //leaf area at a specific segment
	double leafVolAtSeg(int localSegId, bool realized = false, bool withPetiole = false); //leaf area at a specific segment
	double leafLengthAtSeg(int localSegId, bool withPetiole = false);
	std::vector<double> getLeafVisX(int i);
	std::vector<Vector3d> getLeafVis(int i); // per node

    std::string toString() const override;

	/* exact from analytical equations */
	double calcCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
	double calcLength(double age); ///< analytical length of the leaf
	double calcAge(double length) const; ///< analytical age of the leaf

	/* abbreviations */
	std::shared_ptr<LeafRandomParameter> getLeafRandomParameter() const;  ///< root type parameter of this root
	std::shared_ptr<const LeafSpecificParameter> param() const; ///< root parameter

	/* useful */
    bool hasMoved() const override { return true; }; ///< always need to update the coordinates of the nodes for the MappedPlant
	double orgVolume(double length_ = -1.,  bool realized = false) const override;
	double orgVolume2Length(double volume_) override;	
	bool nodeLeafVis(double l); ///<  leaf base (false), branched leaf (false), or leaf surface area (true)

protected:

	Vector3d getIncrement(const Vector3d& p, double sdx, int n = -1) override; ///< called by createSegments, to determine growth direction
    Vector3d heading(int n)  const override; ///< current (absolute) heading of the organs at node n
    int getleafphytomerID(int subtype);
    void minusPhytomerId(int subtype);
    void addleafphytomerID(int subtype);
	double beta;

	std::vector<double> getLeafVisX_(double l);
	bool ageDependentTropism = false;///< do we need to check the leaf's age to see when to update the tropism effect?, @see Leaf::rel2abs
    
};

} // namespace CPlantBox

#endif
