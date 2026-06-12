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
        double sheathLength, double leafGrowthDuration, double bladeBend)
        : OrganSpecificParameter(subType, a),
          bladeAngle(bladeAngle), bladeWidth(bladeWidth), 
          bladeLength(bladeLength),sheathLength(sheathLength), 
          leafGrowthDuration(leafGrowthDuration), bladeBending(bladeBend) 
          { }

    double bladeAngle = 0.;        // Angle between blade and stem
    double bladeWidth = 0.;        // Width of the blade [cm]
    double bladeLength = 0.;       // Length of the blade [cm]
    double sheathLength = 0.;      // Length of the sheath [cm]
    double leafGrowthDuration = 0.;// Total duration of growth [day]
    double bladeBending = 0.05;    // Curvature of the blade

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

    GrassLeafRandomParameter(std::shared_ptr<Organism> plant); 
    virtual ~GrassLeafRandomParameter() { }

    std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant) override;
    std::shared_ptr<OrganSpecificParameter> realize() override; 
    std::string toString(bool verbose = true) const override; 
    void readXML(tinyxml2::XMLElement* element, bool verbose) override; 

    /*
     * Grass leaf parameters (mean and standard deviation)
     */
    double bladeAngle = 0.3;      // Angle between blade and stem [rad]
    double bladeAngles = 0.;      // SD of blade angle [rad]
    double bladeWidth = 0.5;      // Mean width of the blade [cm]
    double bladeWidths = 0.;      // SD of blade width [cm]
    double bladeLength = 10.;     // Mean length of the blade [cm]
    double bladeLengths = 0.;     // SD of blade length [cm]
    double sheathLength = 5.;     // Mean length of the sheath [cm]
    double sheathLengths = 0.;    // SD of sheath length [cm]
    double bladeBending = 0.05;   // Mean blade curvature  [rad/cm]
    double bladeBendings = 0.;    // SD of blade bending [rad/cm]
    double leafGrowthDuration = 20.; // Mean duration of growth [day]
    double leafGrowthDurations = 0.; // SD of duration [day]

    std::shared_ptr<SoilLookUp> f_se = std::make_shared<SoilLookUp>(); 

protected:

    void bindParameters() override; ///< sets up class introspection
};

} // end namespace CPlantBox

#endif
