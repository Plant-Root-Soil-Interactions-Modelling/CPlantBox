// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "seedparameter.h"

#include "Organism.h"
#include "tropism.h"

#include <cmath>
#include <iostream>
#include <numeric>

namespace CPlantBox {

/**
 * @copydoc OrganParameter::toString()
 */
std::string SeedSpecificParameter::toString() const {
    std::stringstream str;
    str << "subType\t" << subType << std::endl;
    str << "seedPos\t" << seedPos.toString() << std::endl;
    str << "firstB\t" << firstB << std::endl << "delayB\t" << delayB << std::endl << "maxB\t" << maxB << std::endl;
    str << "nC\t" << nC << std::endl << "firstSB\t" << firstSB << std::endl << "delaySB\t" << delaySB << std::endl;
    str << "delayRC\t" << delayRC << std::endl << "nz\t" << nz << std::endl << "maxTil\t" << maxTil << std::endl;
    str << "firstTi\t" << firstTil << std::endl << "delayTi\t" << delayTil << std::endl << "maxTil\t" << maxTil << std::endl;
    str << "simtime\t" << simtime << std::endl;
    return str.str();
}

/**
 * Default constructor sets up hashmaps for class introspection
 */
SeedRandomParameter::SeedRandomParameter(std::shared_ptr<Organism> plant) : OrganRandomParameter(plant) {
    // base class default values
    name = "undefined";
    organType = Organism::ot_seed;
    subType = 0;
    bindParameters();
    delayDefinitionRoot = Organism::dd_distance;   // for roots
    delayDefinitionShoot = Organism::dd_time_self; // for leaf or stem laterals
}

/**
 * @copydoc OrganTypeParameter::copy()
 */
std::shared_ptr<OrganRandomParameter> SeedRandomParameter::copy(std::shared_ptr<Organism> plant) {
    // std::cout << "SeedRandomParameter::copy\n"<< std::flush;
    auto s = std::make_shared<SeedRandomParameter>(*this); // copy constructor breaks class introspection
    s->plant = plant;
    s->bindParameters(); // fix class introspection
    return s;
}

/**
 * @copydoc OrganTypeParameter::realize()
 *
 * Creates a specific plant from the plant random parameters.
 * @return Specific plant parameters derived from the root random parameters
 */
std::shared_ptr<OrganSpecificParameter> SeedRandomParameter::realize() {
    auto p = std::make_shared<SeedSpecificParameter>();
    auto pl = plant.lock();
    p->subType = subType;
    p->a = std::max(a + pl->randn() * as, 0.);
    p->ldelay = std::max(ldelay + pl->randn() * ldelays, 0.);
    p->seedPos = seedPos.plus(Vector3d(pl->randn() * seedPoss.x, pl->randn() * seedPoss.y, pl->randn() * seedPoss.z));
    p->firstB = std::max(firstB + pl->randn() * firstBs, 0.);
    p->delayB = std::max(delayB + pl->randn() * delayBs, 0.);
    p->maxB = (int)(std::max(maxB + pl->randn() * maxBs, 0.) + 0.5);
    p->nC = std::max(nC + pl->randn() * nCs, 0.);
    p->firstSB = std::max(firstSB + pl->randn() * firstSBs, 0.);
    p->delaySB = std::max(delaySB + pl->randn() * delaySBs, 0.);
    p->delayRC = std::max(delayRC + pl->randn() * delayRCs, 0.);
    p->nz = std::max(nz + pl->randn() * nzs, 0.);
    p->simtime = std::max(simtime + pl->randn() * simtimes, 0.);
    p->firstTil = std::max(firstTil + pl->randn() * firstTils, 0.);
    p->delayTil = std::max(delayTil + pl->randn() * delayTils, 0.);
    p->maxTil = std::max(maxTil + pl->randn() * maxTils, 0.);
    return p;
}

/**
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void SeedRandomParameter::bindParameters() {
    OrganRandomParameter::bindParameters();
    bindParameter("seedPos.x", &seedPos.x, "X-Coordinate of seed position [cm]", &seedPoss.x);
    bindParameter("seedPos.y", &seedPos.y, "Y-Coordinate of seed position [cm]", &seedPoss.y);
    bindParameter("seedPos.z", &seedPos.z, "Z-Coordinate of seed position [cm]", &seedPoss.z);
    bindParameter("firstB", &firstB, "Emergence of first basal root [day]", &firstBs);
    bindParameter("delayB", &delayB, "Time delay between the basal roots [day]", &delayBs);
    bindParameter("maxB", &maxB, "Maximal number of basal roots [1]", &maxBs);
    bindParameter("nC", &nC, "Maximal number of roots per root crown [1]", &nCs);
    bindParameter("firstSB", &firstSB, "First emergence of a shoot borne root [day]", &firstSBs);
    bindParameter("delaySB", &delaySB, "Time delay between the shoot borne roots [day]", &delaySBs);
    bindParameter("delayRC", &delayRC, "Delay between the root crowns [day]", &delayRCs);
    bindParameter("nz", &nz, "Distance between the root crowns along the shoot [cm]", &nzs);
    bindParameter("firstTil", &firstTil, "Emergence of first tiller [day]", &firstTils);
    bindParameter("delayTil", &delayTil, "Time delay between the tillers [day]", &delayTils);
    bindParameter("maxTil", &maxTil, "Maximal number of tillers [1]", &maxTils);
    bindParameter("simulationTime", &simtime, "Recommended final simulation time  [day]", &simtimes);
    bindParameter("delayDefinitionRoot", &delayDefinitionRoot, "Method implemented to evaluate root lateral growth delay (0: distance based, 1: delay based defined by parent organ)");
    bindParameter("delayDefinitionShoot", &delayDefinitionShoot, "Method implemented to evaluate shoot lateral growth delay (0: distance based, 1: delay based defined by parent organ)");
}

} // namespace CPlantBox
