// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef PERIRHIZAL_H_
#define PERIRHIZAL_H_

#include "MappedOrganism.h"

#include <vector>

namespace CPlantBox {

/**
 * Wraps a MappedSegments (or specialisations MappedPlant, MappedRootsystem)
 * and adds functions to retrieve information on the perirhizal zones of single segments.
 * See also Perirhizal.py
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
