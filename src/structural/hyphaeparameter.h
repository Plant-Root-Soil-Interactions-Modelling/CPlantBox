// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef HYPHAEPARAMETER_H_
#define HYPHAEPARAMETER_H_

#include "mymath.h"
#include "soil.h"
#include "growth.h"
#include "organparameter.h"

/**
 * This file describes the classes HyphaeSpecificParameter and HyphaeRandomParameter.
 * HyphaeSpecificParameter are drawn from the HyphaeRandomParameter class
 */

namespace CPlantBox {

class Organism;

/**
 * Parameters of a specific root, its created by RootRandomParameter:realize()
 */
class HyphaeSpecificParameter :public OrganSpecificParameter
{
public:

    HyphaeSpecificParameter(): HyphaeSpecificParameter(-1, 0., 0., 0., 0., 0.) { } ///< Default constructor
    HyphaeSpecificParameter(int subType, double a, double v, double b, double hlt, double theta):
            OrganSpecificParameter(subType, a),  v(v), b(b), hlt(hlt), theta(theta) { }; ///< Constructor setting all parameters

    int order = 0; // internal counter (? todo)

    double v;              ///< tip elongation rate [cm/day]
    double b;              ///< branching rate [1/day]
    double hlt;            ///< hyphal lifetime  [day]
    double theta;          ///< branching angle [rad]

    std::string toString() const override; ///< for debugging
    double getMaxLength() const; ///< returns maximal length
};


/**
 * Contains a parameter set describing a root type
 */
class HyphaeRandomParameter :public OrganRandomParameter
{

public:

    HyphaeRandomParameter(std::shared_ptr<Organism> plant); ///< default constructor
    virtual ~HyphaeRandomParameter() { };

    std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant) override;

    std::shared_ptr<OrganSpecificParameter> realize() override; ///< Creates a specific root from the root parameter set

    std::string toString(bool verbose = true) const override; ///< info for debugging

    void readXML(tinyxml2::XMLElement* element, bool verbose) override; ///< reads a single sub type organ parameter set

    void bindParameters() override; ///<sets up class introspection

    double v=0.13;              ///< tip elongation rate [cm/day] 
    double vs=0.;             ///< standard deviation of tip elongation rate [cm/day]
    double b=0.5;              ///< branching rate [1/day]
    double bs=0.01;             ///< standard deviation of branching rate [1/day]
    double hlt=10;            ///< hyphal lifetime  [day]
    double hlts=0.1;           ///< standard deviation of  hyphal lifetime  [day]
    double theta=60./180.*M_PI;          ///< branching angle [rad]
    double thetas=0.;         ///< standard deviation of branching angle  [rad]
    double distTT = 0.; ///< distance for tip tip anastomosis [cm]
    double distTH = 0.; ///< distance for tip hyphae anastomosis [cm]

    int tropismT = 2;       ///< Root tropism parameter (Type)
    double tropismN = 1.;   ///< Root tropism parameter (number of trials)
    double tropismS = 0.3;  ///< Root tropism parameter (mean value of expected changeg [1/cm]

    std::shared_ptr<SoilLookUp> v_scale = std::make_shared<SoilLookUp>(); ///< elongation rate scale
    std::shared_ptr<SoilLookUp> b_scale = std::make_shared<SoilLookUp>(); ///< scale branching rate

};

} // end namespace CPlantBox

#endif
