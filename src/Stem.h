#ifndef STEM_H_
#define STEM_H_

#include "Organ.h"
#include "Organism.h"
#include "stemparameter.h"

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
class Stem : public Organ
{
public:

    static std::vector<int> phytomerId;

    Stem(int id,  std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
    		Matrix3d iHeading, int pni, bool moved = true, int oldNON = 0);
    Stem(std::shared_ptr<Organism> plant, int type, Matrix3d iHeading, double delay, std::shared_ptr<Organ> parent, int pni); ///< used within simulation
    virtual ~Stem() { };

    std::shared_ptr<Organ> copy(std::shared_ptr<Organism> plant) override;   ///< deep copies the organ tree

    int organType() const override { return Organism::ot_stem; } ///< returns the organs type

    void simulate(double dt, bool silence = false) override; ///< stem growth for a time span of \param dt
	void internodalGrowth(double dl,double dt_, bool silence = false); ///< internodal growth of \param dl [cm]
    std::string toString() const override;

    /* abbreviations */
	std::shared_ptr<GrowthFunction> getF_gf() override {return getStemRandomParameter()->f_gf;}
	std::shared_ptr<Tropism> getF_tf() override {return getStemRandomParameter()->f_tf;}
    std::shared_ptr<StemRandomParameter> getStemRandomParameter() const;  ///< organ type parameter of this organ
    std::shared_ptr<const StemSpecificParameter> param() const; ///< organ parameter

    int shootborneType = 5;

	/* orientation */
	Vector3d getiHeading()  const override;
	
    bool hasMoved() const override { return true; }; ///< have any nodes moved during the last simulate call
																										 
protected:
	Vector3d partialIHeading;
    void minusPhytomerId(int subtype) { phytomerId[subtype]--;  }
    int getphytomerId(int subtype) { return phytomerId[subtype]; }
    void addPhytomerId(int subtype) { phytomerId[subtype]++;  }

    void createLateral(bool silence); ///< creates a new lateral, called by Stem::simulate()
    void leafGrow(bool silence);
    void shootBorneRootGrow(bool silence);


};

} // namespace CPlantBox

#endif
