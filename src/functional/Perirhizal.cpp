// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Perirhizal.h"

namespace CPlantBox {

/**
 * Calculates outer segment radii [cm], so that the summed segment volumes per cell equals the cell volume
 * @param type          prescribed cylinder volume proportional to 0: segment volume, 1: segment surface, 2: segment length
 * @param vols          (optional) in case of non-equidistant grids, volumes per cell must be defined
 */
std::vector<double> Perirhizal::segOuterRadii(int type, const std::vector<double>& vols) const {
    double cellVolume;
    auto& ms = this->ms;
    auto radii = ms->getEffectiveRadii(); // rename
    auto lengths =  ms->segLength();
    auto width = ms->maxBound.minus(ms->minBound);
    std::vector<double> outer_radii = std::vector<double>(ms->segments.size());
    std::fill(outer_radii.begin(), outer_radii.end(), 0.);
    for(auto iter = ms->cell2seg.begin(); iter != ms->cell2seg.end(); ++iter) {
        int cellId =  iter->first;
        if (vols.size()==0) {
            cellVolume = width.x*width.y*width.z/ms->resolution.x/ms->resolution.y/ms->resolution.z;
        } else {
            cellVolume = vols.at(cellId);
        }
        auto segs = ms->cell2seg.at(cellId);
        double v = 0.;  // calculate sum of root volumes or surfaces over cell
        for (int i : segs) {
            if (type==0) { // volume
                v += M_PI*(radii[i]*radii[i])*lengths[i];
            } else if (type==1) { // surface
                v += 2*M_PI*radii[i]*lengths[i];
            } else if (type==2) { // length
                v += lengths[i];
            }
        }
        for (int i : segs) { // calculate outer radius
            double l = lengths[i];
            double t =0.; // proportionality factor (must sum up to == 1 over cell)
            if (type==0) { // volume
                t = M_PI*(radii[i]*radii[i])*l/v;
            } else if (type==1) { // surface
                t = 2*M_PI*radii[i]*l/v;
            } else if (type==2) { // length
                t = l/v;
            }
            double targetV = t * cellVolume;  // target volume
            outer_radii[i] = std::sqrt(targetV/(M_PI*l)+radii[i]*radii[i]);
        }
    }
    return outer_radii;
}

}
