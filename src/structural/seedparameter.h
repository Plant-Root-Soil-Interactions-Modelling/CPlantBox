// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef SEEDPARAMETER_H_
#define SEEDPARAMETER_H_

#include "mymath.h"
#include "organparameter.h"

/**
 * This file describes the classes SeedSpecificParameter and SeedRandomParameter.
 * SeedSpecificParameter are drawn from the SeedRandomParameter class
 */

namespace CPlantBox {

class OrganSpecificParameter;

/**
 * SeedSpecificParameter contains all plant specific parameters like planting depth and describing
 * the emergence times of basal and shoot borne roots
 *
 * The model currently rather limited, and we might replace it, if we come up with something better
 */
class SeedSpecificParameter :public OrganSpecificParameter
{

public:
    SeedSpecificParameter():SeedSpecificParameter(0,Vector3d(0.,0.,-3), 1.e9, 1.e9, 0, 0., 1.e9, 1.e9, 1.e9, 1., 0, 30.) { }; ///< Default constructor
    SeedSpecificParameter(int type, Vector3d seedPos, double fB, double dB, int mB, int nC, double fSB, double dSB, double dRC,
        double nz, int maxtil, double simtime,  double fTi = 0., double dTi = 0.):
            OrganSpecificParameter(type, 0.), seedPos(seedPos), firstB(fB), delayB(dB), maxB(mB), nC(nC), firstSB(fSB), delaySB(dSB),
            delayRC(dRC), nz(nz), maxTil(maxtil), firstTi(fTi), delayTi(dTi),  simtime(simtime) {  };
    virtual ~SeedSpecificParameter() { };

    /*
     * RootBox and PlantBox plant parameters
     */
    Vector3d seedPos;   ///< Position of the seed [cm]

    //Basal roots (nodal roots)
    double firstB;      ///< Emergence of first basal root [day]
    double delayB;      ///< Time delay between the basal roots [day]
    int maxB;           ///< Maximal number of basal roots [1]

    //Shoot borne roots (crown roots)
    int nC;             ///< Maximal number of roots per root crown [1]
    double firstSB;     ///< First emergence of a shoot borne root [day]
    double delaySB;     ///< Time delay between the shoot borne roots [day]
    double delayRC;     ///< Delay between the root crowns [day]
    double nz;          ///< Distance between the root crowns along the shoot [cm]

    //Tillers
    int maxTil;         ///< maximal number of tillers
	double firstTi;     ///< First emergence of a shoot borne root [day]
    double delayTi;     ///< Time delay between the shoot borne roots [day]

    //Simulation parameters
    double simtime;     ///< recommended final simulation time

    std::string toString() const override; ///< for debugging
};



/**
 * Contains a parameter set describing a plant
 */
class SeedRandomParameter :public OrganRandomParameter {
public:

	SeedRandomParameter(std::shared_ptr<Organism> plant); ///< default constructor
    virtual ~SeedRandomParameter() { }; ///< nothing to do

    std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant) override;

    std::shared_ptr<OrganSpecificParameter> realize() override; ///< Creates a specific plant seed from the seed parameter set

    // DEPRICATED
    void read(std::istream & cin); ///< reads a single root system parameter set
    void write(std::ostream & cout) const; ///< writes a single root system parameter set

    /* Plant parameters */
    Vector3d seedPos  = Vector3d(0.,0.,-3.); ///< Mean position of the seed [cm]
    Vector3d seedPoss = Vector3d(0.,0.,0.);  ///< Standard deviation of position  [cm]
	int delayDefinition = Organism::dd_distance; ///< how is the delay of the laterals defined

    // Basal roots (nodal roots)
    double firstB = 1.e9;  ///< Mean emergence of first basal root [day]
    double firstBs = 0.;   ///< Standard deviation of emergence of first basal root [day]
    double delayB = 1.e9;  ///< Mean time delay between the basal roots [day]
    double delayBs = 0.;   ///< Standard deviation of time delay between the basal roots [day]
    double maxB = 0.;      ///< Mean maximal number of basal roots [1]
    double maxBs = 0.;     ///< Standard deviation of maximal number of basal roots [1]

    // Shoot borne roots (crown roots)
    double nC = 0.;        ///< Mean maximal number of roots per root crown [1]
    double nCs = 0.;       ///< Standard deviation of maximal number of roots per root crown [1]
    double firstSB = 1.e9; ///< Mean first emergence of a shoot borne root [day]
    double firstSBs = 0.;  ///< Standard deviation of first emergence of a shoot borne root [day]
    double delaySB = 1.e9; ///< Mean time delay between the shoot borne roots [day]
    double delaySBs = 0.;  ///< Standard deviation of time delay between the shoot borne roots [day]
    double delayRC = 1.e9; ///< Mean delay between the root crowns [day]
    double delayRCs = 0.;  ///< Standard deviation of delay between the root crowns [day]
    double nz = 1.;        ///< Mean distance between the root crowns along the shoot [cm]
    double nzs = 0.;       ///< Standard deviation of distance between the root crowns along the shoot [cm]

    // Stem parameters
    int maxTil = 0;        ///< Maximal number of tillers
    double maxTils = 0.;   ///< Standard deviation of tillers
	double firstTi = 1.e9;  ///< Mean emergence of first basal root [day]
    double firstTis = 0.;   ///< Standard deviation of emergence of first basal root [day]
    double delayTi = 1.e9;  ///< Mean time delay between the basal roots [day]
    double delayTis = 0.;   ///< Standard deviation of time delay between the basal roots [day]

    // Simulation parameters
    double simtime = 30.;  ///< Mean recommended final simulation time
    double simtimes = 0.;  ///< Standard deviation of recommended final simulation time

protected:

    void bindParameters() override;   ///< sets up class introspection

};

} // namespace

#endif
