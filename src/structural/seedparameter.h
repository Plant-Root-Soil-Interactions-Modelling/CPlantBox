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
 */
class SeedSpecificParameter : public OrganSpecificParameter {

  public:
    Vector3d seedPos = Vector3d(0., 0., -3); ///< Position of the seed [cm]

    // Basal roots (nodal roots)
    double firstB = 1.e9; ///< Emergence of first basal root [day]
    double delayB = 1.e9; ///< Time delay between the basal roots [day]
    int maxB = 0;         ///< Maximal number of basal roots [1]

    // Shoot borne roots (crown roots)
    int nC = 0;            ///< Maximal number of roots per root crown [1]
    double firstSB = 1.e9; ///< First emergence of a shoot borne root [day]
    double delaySB = 1.e9; ///< Time delay between the shoot borne roots [day]
    double delayRC = 1.e9; ///< Delay between the root crowns [day]
    double nz = 1.;        ///< Distance between the root crowns along the shoot [cm]

    // Tillers
    int maxTil = 0;       ///< maximal number of tillers
    double firstTil = 0.; ///< First emergence of a tiller [day]
    double delayTil = 0.; ///< Time delay between the tillers [day]

    // Simulation parameters
    double simtime = 30.; ///< recommended final simulation time [day]

    std::string toString() const override; ///< for debugging
};

/**
 * Contains a parameter set describing a plant
 */
class SeedRandomParameter : public OrganRandomParameter {
  public:
    SeedRandomParameter(std::shared_ptr<Organism> plant); ///< default constructor
    virtual ~SeedRandomParameter() {};                    ///< nothing to do

    std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant) override;

    std::shared_ptr<OrganSpecificParameter> realize() override; ///< Creates a specific plant seed from the seed parameter set

    /* Plant parameters */
    Vector3d seedPos = Vector3d(0., 0., -3.); ///< Mean position of the seed [cm]
    Vector3d seedPoss = Vector3d(0., 0., 0.); ///< Standard deviation of position  [cm]
    int delayDefinitionRoot;                  ///< how is the delay of the laterals defined
    int delayDefinitionShoot;                 ///< how is the delay of the laterals defined

    // Basal roots (nodal roots)
    double firstB = 1.e9; ///< Mean emergence of first basal root [day]
    double firstBs = 0.;  ///< Standard deviation of emergence of first basal root [day]
    double delayB = 1.e9; ///< Mean time delay between the basal roots [day]
    double delayBs = 0.;  ///< Standard deviation of time delay between the basal roots [day]
    double maxB = 0.;     ///< Mean maximal number of basal roots [1]
    double maxBs = 0.;    ///< Standard deviation of maximal number of basal roots [1]

    // Shoot borne roots (crown roots)
    double nC = 7.;        ///< Mean maximal number of roots per root crown [1]
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
    int maxTil = 0;         ///< Maximal number of tillers
    double maxTils = 0.;    ///< Standard deviation of tillers
    double firstTil = 1.e9; ///< Mean emergence of first basal root [day]
    double firstTils = 0.;  ///< Standard deviation of emergence of first basal root [day]
    double delayTil = 1.e9; ///< Mean time delay between the basal roots [day]
    double delayTils = 0.;  ///< Standard deviation of time delay between the basal roots [day]

    // Simulation parameters
    double simtime = 300.; ///< Mean recommended final simulation time
    double simtimes = 0.;  ///< Standard deviation of recommended final simulation time

  protected:
    void bindParameters() override; ///< sets up class introspection
};

} // namespace CPlantBox

#endif
