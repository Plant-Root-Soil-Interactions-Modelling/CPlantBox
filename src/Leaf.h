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
 * organ
 *
 * Describes a single organ, by a vector of nodes representing the organ.
 * The method simulate() creates new nodes of this organ, and lateral organs in the organ's branching zone.
 *
 */
class Leaf : public Organ
{
public:

    Leaf(int id,  std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
        Matrix3d iHeading, int pni, bool moved = true, int oldNON = 0);
	Leaf(std::shared_ptr<Organism> plant, int type, Matrix3d iHeading, double delay, std::shared_ptr<Organ> parent, int pni); ///< used within simulation
	virtual ~Leaf() { };

	std::shared_ptr<Organ> copy(std::shared_ptr<Organism> plant) override;   ///< deep copies the organ tree
	int organType() const override { return Organism::ot_leaf; } ///< returns the organs type
	void simulate(double dt, bool silence = false) override; ///< organ growth for a time span of \param dt
	double getParameter(std::string name) const override; ///< returns an organ parameter of Plant::ScalarType
	std::shared_ptr<GrowthFunction> getF_gf() override {return getLeafRandomParameter()->f_gf;}
	std::shared_ptr<Tropism> getF_tf() override {return getLeafRandomParameter()->f_tf;}

	/* leaf vizualisation */
    double leafLength() const { return std::max(getLength(false)-param()->lb, 0.); /* represents the leaf base*/ }; ///< leaf surface length [cm]
    double leafCenter() const { return std::max(getLength(false)-param()->la-param()->lb, 0.); }; ///< center of the radial parametrisation
    double leafArea() ; ///< returns the leaf surface area, zero if there are lateral-leafs [cm2]
	std::vector<double> getLeafVisX(int i);
	std::vector<Vector3d> getLeafVis(int i); // per node

    std::string toString() const override;

	
	/* abbreviations */
	std::shared_ptr<LeafRandomParameter> getLeafRandomParameter() const;  ///< organ type parameter of this organ
	std::shared_ptr<const LeafSpecificParameter> param() const; ///< organ parameter

	Vector3d getiHeading() const override;///< compute initial heading from 
	bool hasMoved() const override { return true; }; ///< always need to update the coordinates of the nodes for the MappedPlant
												   

protected:

    Vector3d partialIHeading;


    int getleafphytomerID(int subtype);
    void minusPhytomerId(int subtype);
    void addleafphytomerID(int subtype);

    void createLateral(bool silence); ///< creates a new lateral, called by Leaf::simulate()
	bool nodeLeafVis(double l); ///<  leaf base (false), branched leaf (false), or leaf surface area (true)
	std::vector<double> getLeafVisX_(double l);

};

} // namespace CPlantBox

#endif
