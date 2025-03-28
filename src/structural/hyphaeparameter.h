// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef HYPHAEPARAMETER_H_
#define HYPHAEPARAMETER_H_

#include "mymath.h"
#include "soil.h"
#include "growth.h"
#include "organparameter.h"

/**
 * This file describes the classes RootSpecificParameter and RootRandomParameter.
 * RootSpecificParameter are drawn from the RootRandomParameter class
 */

namespace CPlantBox {

class Organism;

/**
 * Parameters of a specific root, its created by RootRandomParameter:realize()
 */
class HyphaeSpecificParameter :public OrganSpecificParameter
{

public:

    HyphaeSpecificParameter(): HyphaeSpecificParameter(-1, 0., 0., std::vector<double>(0), 0., 0., 0., 0.) { } ///< Default constructor
    HyphaeSpecificParameter(int type, v, b, hlt, theta):
            OrganSpecificParameter(type, a),  lb(lb), la(la), r(r), theta(theta), rlt(rlt), ln(ln), laterals(laterals) { }; ///< Constructor setting all parameters

    double v;              ///< Basal zone [cm]
    double b;              ///< Apical zone [cm]
    double hlt;               ///< Initial growth rate [cm day-1]
    double theta;           ///< Angle between root and parent root [rad]

    std::string toString() const override; ///< for debugging

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

    double nob() const { return std::max((lmax-la-lb)/ln+1, 0.); }  ///< returns the mean maximal number of branching nodes [1]
    double nobs() const; ///< returns the standard deviation of number of branching nodes [1]

    std::string toString(bool verbose = true) const override; ///< info for debugging

    void readXML(tinyxml2::XMLElement* element, bool verbose) override; ///< reads a single sub type organ parameter set

    void bindParameters() override; ///<sets up class introspection




protected:

    double snap(double x);

};

} // end namespace CPlantBox

#endif
