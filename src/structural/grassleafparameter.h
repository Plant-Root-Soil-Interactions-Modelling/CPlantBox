// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef GRASSLEAFPARAMETER_H_
#define GRASSLEAFPARAMETER_H_

#include "mymath.h"
#include "soil.h"
#include "growth.h"
#include "organparameter.h"

namespace CPlantBox {

class Organism;

/**
 * Parameters of a specific grass leaf instance, created by GrassLeafRandomParameter::realize().
 *
 * A grass leaf consists of two zones:
 *  - the sheath, which wraps around the stem,
 *  - the blade, the flat elongated part that emerges after a delay.
 */
class GrassLeafSpecificParameter : public OrganSpecificParameter
{
public:

    GrassLeafSpecificParameter() : OrganSpecificParameter(-1, 0.) { }
    GrassLeafSpecificParameter(int subType, double a,
        double bladeAngle, double bladeWidth, double bladeLength,
        double sheathLength, double sheathDuration,
        double bladeDelay, double bladeDuration)
        : OrganSpecificParameter(subType, a),
          bladeAngle(bladeAngle), bladeWidth(bladeWidth), bladeLength(bladeLength),
          sheathLength(sheathLength), sheathDuration(sheathDuration),
          bladeDelay(bladeDelay), bladeDuration(bladeDuration) { }

    double bladeAngle = 0.;      ///< Angle between blade and stem [rad]
    double bladeWidth = 0.;      ///< Width of the blade [cm]
    double bladeLength = 0.;     ///< Length of the blade [cm]
    double sheathLength = 0.;    ///< Length of the sheath [cm]
    double sheathDuration = 0.;  ///< Duration of sheath development [day]
    double bladeDelay = 0.;      ///< Delay before blade emergence [day]
    double bladeDuration = 0.;   ///< Duration of blade growth [day]

    std::string toString() const override; ///< for debugging
};


/**
 * Contains a stochastic parameter set describing a grass leaf sub-type.
 *
 * Each parameter is described by a mean value and a standard deviation (suffix "s").
 * GrassLeafSpecificParameter instances are created by realize().
 */
class GrassLeafRandomParameter : public OrganRandomParameter
{
public:

    GrassLeafRandomParameter(std::shared_ptr<Organism> plant); ///< default constructor
    virtual ~GrassLeafRandomParameter() { }

    std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant) override;

    std::shared_ptr<OrganSpecificParameter> realize() override; ///< Creates a specific grass leaf from the parameter set

    std::string toString(bool verbose = true) const override; ///< info for debugging

    void readXML(tinyxml2::XMLElement* element, bool verbose) override; ///< reads a single sub-type organ parameter set

    /*
     * Grass leaf parameters (mean and standard deviation)
     */
    double bladeAngle = 0.3;     ///< Mean angle between blade and stem [rad]
    double bladeAngles = 0.;     ///< Standard deviation of blade angle [rad]
    double bladeWidth = 0.5;     ///< Mean width of the blade [cm]
    double bladeWidths = 0.;     ///< Standard deviation of blade width [cm]
    double bladeLength = 10.;    ///< Mean length of the blade [cm]
    double bladeLengths = 0.;    ///< Standard deviation of blade length [cm]
    double sheathLength = 5.;    ///< Mean length of the sheath [cm]
    double sheathLengths = 0.;   ///< Standard deviation of sheath length [cm]
    double sheathDuration = 10.; ///< Mean duration of sheath development [day]
    double sheathDurations = 0.; ///< Standard deviation of sheath duration [day]
    double bladeDelay = 2.;      ///< Mean delay before blade emergence [day]
    double bladeDelays = 0.;     ///< Standard deviation of blade delay [day]
    double bladeDuration = 20.;  ///< Mean duration of blade growth [day]
    double bladeDurations = 0.;  ///< Standard deviation of blade duration [day]

protected:

    void bindParameters() override; ///< sets up class introspection
};

} // end namespace CPlantBox

#endif
