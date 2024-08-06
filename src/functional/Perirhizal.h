// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef PERIRHIZAL_H_
#define PERIRHIZAL_H_

#include "MappedOrganism.h"

#include <vector>

namespace CPlantBox {

/**
 * Wraps a MappedSegments (or specialisations MappedPlant, MappedRootsystem)
 * and adds functions to help modelling the perirhizal zone.
 *
 * Currently rather useless, main part is in Python (in Perirhizal.py),
 * but perfomance critical methods could be implemented in C++
 */
class Perirhizal
{
public:

    Perirhizal() { }
    Perirhizal(std::shared_ptr<MappedSegments> ms) :ms(ms) { }

    std::vector<double> segOuterRadii(int type, const std::vector<double>& vols = std::vector<double>(0)) const; ///< outer cylinder radii to match cell volume

    std::shared_ptr<MappedSegments> ms;
};

}

#endif
