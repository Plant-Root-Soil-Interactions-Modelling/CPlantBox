// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "seedparameter.h"

#include "Organism.h"

#include <cmath>
#include <iostream>
#include <numeric>

namespace CPlantBox {

/**
 * @copydoc OrganParameter::toString()
 */
std::string SeedSpecificParameter::toString() const
{
    std::stringstream str;
    str << "subType\t" << subType << std::endl;
    str << "seedPos\t" << seedPos.toString() << std::endl;
    str << "firstB\t" << firstB << std::endl << "delayB\t" << delayB << std::endl << "maxB\t" << maxB << std::endl;
    str << "nC\t" << nC << std::endl << "firstSB\t" << firstSB << std::endl << "delaySB\t" << delaySB << std::endl;
    str << "delayRC\t" << delayRC << std::endl << "nz\t" << nz << std::endl << "maxTil\t" << maxTil << std::endl;
	str << "firstTi\t" << firstTi << std::endl << "delayTi\t" << delayTi << std::endl<< "maxTil\t" << maxTil << std::endl;
    str << "simtime\t" << simtime << std::endl;
    return str.str();
}



/**
 * Default constructor sets up hashmaps for class introspection
 */
SeedRandomParameter::SeedRandomParameter(std::shared_ptr<Organism> plant) :OrganRandomParameter(plant)
{
    // base class default values
    name = "undefined";
    organType = Organism::ot_seed;
    subType = 0;
    bindParameters();
}

/**
 * @copydoc OrganTypeParameter::copy()
 */
std::shared_ptr<OrganRandomParameter> SeedRandomParameter::copy(std::shared_ptr<Organism> plant)
{
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
std::shared_ptr<OrganSpecificParameter> SeedRandomParameter::realize()
{
    auto p = plant.lock();
    Vector3d sP = seedPos.plus(Vector3d(p->randn()*seedPoss.x, p->randn()*seedPoss.y, p->randn()*seedPoss.z));
    double fB = std::max(firstB + p->randn()*firstBs, 0.);
    double dB = std::max(delayB + p->randn()*delayBs, 0.);
	int mB = (int)(std::max(maxB + p->randn()*maxBs, 0.) +0.5);
    double nC_ = std::max(nC + p->randn()*nCs, 0.);

    double fSB = std::max(firstSB + p->randn()*firstSBs, 0.);
    double dSB = std::max(delaySB + p->randn()*delaySBs, 0.);
    double dRC = std::max(delayRC + p->randn()*delayRCs, 0.);
    double nz_ = std::max(nz , 0.);

    double st = std::max(simtime + p->randn()*simtimes, 0.);
	double fTi = std::max(firstTi + p->randn()*firstTis, 0.);
    double dTi = std::max(delayTi + p->randn()*delayTis, 0.);
    int maxtil = std::max(maxTil + p->randn()*maxTils, 0.);
    return std::make_shared<SeedSpecificParameter>(subType, sP, fB, dB, mB, nC_, fSB, dSB,dRC, nz_, maxtil, st,fTi, dTi);
}

/**
 * depricated todo no RootSystemTypeParameter yet
 */
void SeedRandomParameter::read(std::istream & cin)
{
	std::cout << "SeedRandomParameter::read is deprecated, use readXML instead \n";
    double plantingdepth;
    std::string s; // dummy
    cin  >>  s >> plantingdepth;
    cin >> s >> firstB >> s >> delayB >> s >> maxB >> s >> nC >> s >> firstSB >> s >> delaySB >> s >> delayRC >> s >> nz >> s >> simtime;
    seedPos = Vector3d(0,0,-plantingdepth);
}

/**
 * depricated todo no RootSystemTypeParameter yet
 */
void SeedRandomParameter::write(std::ostream & cout) const
{
	std::cout << "SeedRandomParameter::write is deprecated, use writeXML instead \n";
    double pd = -seedPos.z;
    cout <<  "plantingdepth\t" << pd << "\n" <<  "firstB\t" << firstB << "\n" <<  "delayB\t" << delayB << "\n"
        <<  "maxB\t" << maxB << "\n" <<  "nC\t" << nC << "\n" <<  "firstSB\t" << firstSB << "\n"
        <<  "delaySB\t" << delaySB << "\n" <<  "delayRC\t" << delayRC << "\n" <<  "nz\t" << nz << "\n" << "simulationTime\t" << simtime << "\n";
}

/**
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void SeedRandomParameter::bindParameters()
{
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
    bindParameter("nz", &nz, "Distance between the root crowns along the shoot [cm]", &nzs );
	bindParameter("firstTi", &firstTi, "Emergence of first tiller [day]", &firstTis);
    bindParameter("delayTi", &delayTi, "Time delay between the tillers [day]", &delayTis);
    bindParameter("maxTi", &maxTil, "Maximal number of tillers [1]", &maxTils);
    bindParameter("simulationTime", &simtime, "Recommended final simulation time  [day]", &simtimes );
	bindParameter("delayDefinition", &delayDefinition, "method implemented to evaluate lateral growth delay (0: distance based, 1: delay based defined by parent organ)");
}

} // namespace
