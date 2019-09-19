// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "seedparameter.h"

#include "Organism.h"

#include <cmath>
#include <iostream>

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
    str << "delayRC\t" << delayRC << std::endl << "nz\t" << nz << std::endl << "simtime\t" << simtime << std::endl;
    return str.str();
}



/**
 * Default constructor sets up hashmaps for class introspection
 */
SeedRandomParameter::SeedRandomParameter(Organism* p) :OrganRandomParameter(p)
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
OrganRandomParameter* SeedRandomParameter::copy(Organism* p)
{
	SeedRandomParameter* s = new SeedRandomParameter(*this); // copy constructor breaks class introspection
    s->plant = p;
    s->bindParameters(); // fix class introspection
    return s;
}

/**
 * @copydoc OrganTypeParameter::realize()
 *
 * Creates a specific plant from the plant random parameters.
 * @return Specific plant parameters derived from the root random parameters
 */
OrganSpecificParameter* SeedRandomParameter::realize()
{
    Vector3d sP = seedPos.plus(Vector3d(plant->randn()*seedPoss.x, plant->randn()*seedPoss.y, plant->randn()*seedPoss.z));
    double fB = std::max(firstB + plant->randn()*firstBs, 0.);
    double dB = std::max(delayB + plant->randn()*delayBs, 0.);
    int mB = (int)(std::max(maxB + plant->randn()*maxBs, 0.) +0.5);
    double nC_ = std::max(nC + plant->randn()*nCs, 0.);
    double fSB = std::max(firstSB + plant->randn()*firstSBs, 0.);
    double dSB = std::max(delaySB + plant->randn()*delaySBs, 0.);
    double dRC = std::max(delayRC + plant->randn()*delayRCs, 0.);
    double nz_ = std::max(delaySB + plant->randn()*delaySBs, 0.);
    double st = std::max(simtime + plant->randn()*simtimes, 0.);
    OrganSpecificParameter* p = new SeedSpecificParameter(subType, sP, fB, dB, mB, nC_, fSB, dSB,dRC, nz_, st);
    return p;
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
    bindParameter("organType", &organType, "Organ type (unspecified organ = 0, seed = 1, root = 2, stem = 3, leaf = 4)");
    bindParameter("subType", &subType, "Unique identifier of this sub type");
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
    bindParameter("simtime", &simtime, "Recommended final simulation time  [day]", &simtimes );
    // other parameters (descriptions only)
    description["name"]  = "Name of the sub type of the organ, e.g. small lateral";
}

} // namespace
