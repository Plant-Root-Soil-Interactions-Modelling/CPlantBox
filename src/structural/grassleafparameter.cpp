// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "grassleafparameter.h"

#include "Organism.h"
#include "tropism.h"

#include <cmath>
#include <iostream>

namespace CPlantBox {

/**
 * @copydoc OrganSpecificParameter::toString()
 */
std::string GrassLeafSpecificParameter::toString() const {
    std::stringstream str;
    str << "GrassLeaf" << std::endl;
    str << "subType\t" << subType << std::endl;
    str << "a\t" << a << std::endl;
    str << "bladeAngle\t" << bladeAngle << std::endl;
    str << "bladeWidth\t" << bladeWidth << std::endl;
    str << "bladeLength\t" << bladeLength << std::endl;
    str << "sheathLength\t" << sheathLength << std::endl;
    str << "bladeBending\t" << bladeBending << std::endl;
    str << "leafGrowthDuration\t" << leafGrowthDuration << std::endl;
    return str.str();
}

/**
 * Default constructor sets up hashmaps for class introspection.
 */
GrassLeafRandomParameter::GrassLeafRandomParameter(std::shared_ptr<Organism> plant) : OrganRandomParameter(plant) {
    name = "undefined";
    organType = Organism::ot_leaf;
    subType = -1;
    f_tf = std::make_shared<Tropism>(plant);
    f_gf = std::make_shared<LinearGrowth>(); // fixed to linear growth
    bindParameters();
}

/**
 * @copydoc OrganRandomParameter::copy()
 */
std::shared_ptr<OrganRandomParameter> GrassLeafRandomParameter::copy(std::shared_ptr<Organism> plant) {
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
std::shared_ptr<OrganSpecificParameter> GrassLeafRandomParameter::realize() {
    auto p = std::make_shared<GrassLeafSpecificParameter>();
    auto pl = plant.lock();
    p->subType = subType;
    p->a = std::max(a + pl->randn() * as, 0.);
    p->ldelay = std::max(ldelay + pl->randn() * ldelays, 0.);
    p->bladeAngle = std::max(bladeAngle + pl->randn() * bladeAngles, 0.);
    p->bladeWidth = std::max(bladeWidth + pl->randn() * bladeWidths, 0.);
    p->bladeLength = std::max(bladeLength + pl->randn() * bladeLengths, 0.);
    p->sheathLength = std::max(sheathLength + pl->randn() * sheathLengths, 0.);
    p->bladeBending = std::max(bladeBending + pl->randn() * bladeBendings, 0.);
    p->leafGrowthDuration = std::max(leafGrowthDuration + pl->randn() * leafGrowthDurations, 0.);
    return p;
}

/**
 * @copydoc OrganRandomParameter::toString()
 */
std::string GrassLeafRandomParameter::toString(bool verbose) const { return OrganRandomParameter::toString(verbose); }

/**
 * @copydoc OrganRandomParameter::readXML()
 */
void GrassLeafRandomParameter::readXML(tinyxml2::XMLElement *element, bool verbose) { OrganRandomParameter::readXML(element, verbose); }

/**
 * Registers all scalar parameters in the introspection maps.
 * Called by the constructor; derived classes must call base first.
 */
void GrassLeafRandomParameter::bindParameters() {
    OrganRandomParameter::bindParameters();
    bindParameter("bladeAngle", &bladeAngle, "Angle between blade and stem [rad]", &bladeAngles);
    bindParameter("bladeWidth", &bladeWidth, "Width of the blade [cm]", &bladeWidths);
    bindParameter("bladeLength", &bladeLength, "Length of the blade [cm]", &bladeLengths);
    bindParameter("sheathLength", &sheathLength, "Length of the sheath [cm]", &sheathLengths);
    bindParameter("bladeBending", &bladeBending, "Blade curvature: pitch per cm segment length [rad/cm]", &bladeBendings);
    bindParameter("leafGrowthDuration", &leafGrowthDuration, "Total duration of leaf growth [day]", &leafGrowthDurations);
}

} // end namespace CPlantBox
