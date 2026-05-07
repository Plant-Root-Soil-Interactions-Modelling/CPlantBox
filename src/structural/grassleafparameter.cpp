// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "grassleafparameter.h"

#include "Organism.h"
#include "tropism.h"

#include <cmath>
#include <iostream>
#include <assert.h>

namespace CPlantBox {

/**
 * @copydoc OrganSpecificParameter::toString()
 */
std::string GrassLeafSpecificParameter::toString() const
{
    std::stringstream str;
    str << "GrassLeaf" << std::endl;
    str << "subType\t"       << subType      << std::endl;
    str << "a\t"             << a            << std::endl;
    str << "bladeAngle\t"    << bladeAngle   << std::endl;
    str << "bladeWidth\t"    << bladeWidth   << std::endl;
    str << "bladeLength\t"   << bladeLength  << std::endl;
    str << "sheathLength\t"  << sheathLength << std::endl;
    str << "sheathDuration\t"<< sheathDuration << std::endl;
    str << "bladeDelay\t"    << bladeDelay   << std::endl;
    str << "bladeDuration\t" << bladeDuration << std::endl;
    return str.str();
}

/**
 * Default constructor sets up hashmaps for class introspection.
 */
GrassLeafRandomParameter::GrassLeafRandomParameter(std::shared_ptr<Organism> plant)
    : OrganRandomParameter(plant)
{
    name = "undefined";
    organType = Organism::ot_leaf;
    subType = -1;
    f_tf = std::make_shared<Tropism>(plant);
    bindParameters();
}

/**
 * @copydoc OrganRandomParameter::copy()
 */
std::shared_ptr<OrganRandomParameter> GrassLeafRandomParameter::copy(std::shared_ptr<Organism> plant)
{
    auto r = std::make_shared<GrassLeafRandomParameter>(*this); // copy constructor breaks class introspection
    r->plant = plant;
    r->bindParameters(); // fix class introspection
    r->f_tf = f_tf->copy(plant);
    r->f_gf = f_gf->copy();
    return r;
}

/**
 * @copydoc OrganRandomParameter::realize()
 *
 * Draws stochastic realisations of all grass leaf parameters and returns a
 * GrassLeafSpecificParameter instance.
 */
std::shared_ptr<OrganSpecificParameter> GrassLeafRandomParameter::realize()
{
    auto p = plant.lock();

    double a_            = std::max(a            + p->randn() * as,             0.);
    double bladeAngle_   = std::max(bladeAngle   + p->randn() * bladeAngles,    0.);
    double bladeWidth_   = std::max(bladeWidth   + p->randn() * bladeWidths,    0.);
    double bladeLength_  = std::max(bladeLength  + p->randn() * bladeLengths,   0.);
    double sheathLength_ = std::max(sheathLength + p->randn() * sheathLengths,  0.);
    double sheathDuration_ = std::max(sheathDuration + p->randn() * sheathDurations, 0.);
    double bladeDelay_   = std::max(bladeDelay   + p->randn() * bladeDelays,    0.);
    double bladeDuration_ = std::max(bladeDuration + p->randn() * bladeDurations, 0.);

    return std::make_shared<GrassLeafSpecificParameter>(
        subType, a_,
        bladeAngle_, bladeWidth_, bladeLength_,
        sheathLength_, sheathDuration_,
        bladeDelay_, bladeDuration_);
}

/**
 * @copydoc OrganRandomParameter::toString()
 */
std::string GrassLeafRandomParameter::toString(bool verbose) const
{
    return OrganRandomParameter::toString(verbose);
}

/**
 * @copydoc OrganRandomParameter::readXML()
 */
void GrassLeafRandomParameter::readXML(tinyxml2::XMLElement* element, bool verbose)
{
    OrganRandomParameter::readXML(element, verbose);
}

/**
 * Registers all scalar parameters in the introspection maps.
 * Called by the constructor; derived classes must call base first.
 */
void GrassLeafRandomParameter::bindParameters()
{
    OrganRandomParameter::bindParameters();
    bindParameter("bladeAngle",    &bladeAngle,    "Angle between blade and stem [rad]",    &bladeAngles);
    bindParameter("bladeWidth",    &bladeWidth,    "Width of the blade [cm]",               &bladeWidths);
    bindParameter("bladeLength",   &bladeLength,   "Length of the blade [cm]",              &bladeLengths);
    bindParameter("sheathLength",  &sheathLength,  "Length of the sheath [cm]",             &sheathLengths);
    bindParameter("sheathDuration",&sheathDuration,"Duration of sheath development [day]",  &sheathDurations);
    bindParameter("bladeDelay",    &bladeDelay,    "Delay before blade emergence [day]",    &bladeDelays);
    bindParameter("bladeDuration", &bladeDuration, "Duration of blade growth [day]",        &bladeDurations);
}

} // end namespace CPlantBox
