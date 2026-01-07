// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "hyphaeparameter.h"

#include "Organism.h"

#include <cmath>
#include <iostream>
#include <chrono>
#include <assert.h>
#include <numeric>

namespace CPlantBox {

/**
 * @copydoc OrganParameter::toString()
 */
std::string HyphaeSpecificParameter::toString() const
{
    std::stringstream str;
    str << "HyphaeSpecificParameter: subType\t" << subType << ", radius " << a << ", " << "tip elongation rate " << v << ", " << "branching rate " << b << ", " << "hyphal lifetime " << hlt << ", branching angle " << theta << std::endl; 
    return str.str();
}


/**
 * @return Mean maximal hyphae length of this hyphae type
 */
double HyphaeSpecificParameter::getMaxLength() const {
    // double l = 0.5 / (order+1); // v*hlt is too long ie tip elongation rate * hyphal lifetime
    double l = v * hlt; // mean maximal length of a hyphae
    return l; 
}

/**
 * Default constructor sets up hashmaps for class introspection
 */
HyphaeRandomParameter::HyphaeRandomParameter(std::shared_ptr<Organism> plant) :OrganRandomParameter(plant)
{
    // base class default values
    name = "undefined";
    organType = Organism::ot_hyphae;
    subType = -1;
    f_tf = std::make_shared<Tropism>(plant);
    bindParameters();
    f_gf = std::make_shared<LinearGrowth>();
}

/**
 * @copydoc OrganTypeParameter::copy()
 */
std::shared_ptr<OrganRandomParameter> HyphaeRandomParameter::copy(std::shared_ptr<Organism> p)
{
    auto r = std::make_shared<HyphaeRandomParameter>(*this); // copy constructor breaks class introspection
    r->plant = p;
    r->bindParameters(); // fix class introspection
    r->f_tf = f_tf->copy(p); // copy call back function classes
    r->f_gf = f_gf->copy();
    return r;
}

/**
 * @copydoc OrganTypeParameter::realize()
 *
 * Creates a specific hyphae from the hyphae random parameters.
 * @return Specific hyphae parameters derived from the root type parameters
 */
std::shared_ptr<OrganSpecificParameter> HyphaeRandomParameter::realize()
{
    assert(dx > dxMin && "HyphaeRandomParameter::realize(): dxMin must be smaller than dx");
    auto p = plant.lock();

    double a_ = std::max(a + p->randn()*as, 0.); // radius
    double v_ = std::max(v + p->randn()*vs, 0.);  // tip elongation rate [cm/day]
    double b_ = std::max(b + p->randn()*bs, 0.); // branching rate [1/day]
    double hlt_ = std::max(hlt + p->randn()*hlts, 0.); // hyphal lifetime  [day]
    double theta_ = std::max(theta + p->randn()*thetas, 0.); // branching angle [rad]

     return std::make_shared<HyphaeSpecificParameter>(subType, a_, v_, b_, hlt_, theta_);
}

/**
 * @copydoc OrganTypeParameter::toString()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 */
std::string HyphaeRandomParameter::toString(bool verbose) const {

    if (verbose) {
        return OrganRandomParameter::toString(true);
    } else {
        return OrganRandomParameter::toString(false);
    }

}

/**
 * @copydoc OrganTypeParameter::readXML()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 *
 * If the parameter successor or successorP are not in the element, they are set to zero size.
 */
void HyphaeRandomParameter::readXML(tinyxml2::XMLElement* element, bool verbose)
{
    OrganRandomParameter::readXML(element, verbose);
}


/**
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void HyphaeRandomParameter::bindParameters()
{
    OrganRandomParameter::bindParameters();
    bindParameter("v", &v, "Tip elongation rate [cm/day]", &vs);
    bindParameter("b", &b, "Branching rate [1/day]", &bs);
    bindParameter("hlt", &hlt, "Hyphal lifetime  [day]", &hlts);
    bindParameter("theta", &theta, "Branching angle [rad]", &thetas);
    // bindParameter("distTT", &distTT, "Distance for tip-tip anastomosis [cm]");
    bindParameter("distTH", &distTH, "Distance for tip-hyphae anastomosis [cm]");
    bindParameter("tropismT", &tropismT, "Type of root tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)");
    bindParameter("tropismN", &tropismN, "Number of trials of root tropism");
    bindParameter("tropismS", &tropismS, "Mean value of expected change of root tropism [1/cm]");
}

} // end namespace CPlantBox
