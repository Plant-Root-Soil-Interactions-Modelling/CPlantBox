// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Perirhizal.h"

namespace CPlantBox {


void Perirhizal::filter_data(
           int variable,
           int variableC,
		   const std::string& scenario,
           const std::string& pSet_str,
           double timeeval) {
			   
	const std::vector<double>* y_ptr = nullptr;
	const std::vector<double>* cdata_ptr = nullptr;

	switch (variable) {
		case 1: y_ptr = &dfcs_y; break;
		case 5: y_ptr = &dfcca_y; break;
		default: throw std::runtime_error("Invalid variable");
	}
	switch (variableC) {
		case 1: cdata_ptr = &dfcs_y; break;
		case 5: cdata_ptr = &dfcca_y; break;
		case 6: cdata_ptr = &dfccat_y; break;
		default: throw std::runtime_error("Invalid variableC");
	}

	const std::vector<double>& y = *y_ptr;
	const std::vector<double>& cdata = *cdata_ptr;
    // Block A: find closest time
    double min_diff = 1e6;
    double closest_time = times_[0];
    for (double t : times_) {
        double diff = std::abs(t - timeeval);
        if (diff < min_diff) {
            min_diff = diff;
            closest_time = t;
        }
    }

    // Block A: apply time mask
    std::vector<double> dfx, dfy, colordf;
    std::vector<std::string> pSets_filtered, scenarios_filtered;

    for (size_t i = 0; i < times_.size(); ++i) {
        if (times_[i] == closest_time) {
            dfx.push_back(x[i]);
            dfy.push_back(y[i]);
            colordf.push_back(cdata[i]);
            pSets_filtered.push_back(pSets[i]);
            scenarios_filtered.push_back(scenarios[i]);
        }
    }

    // Block B: apply metadata mask
    dfx_final.clear();
    dfy_final.clear();
    colordf_final.clear();

    for (size_t i = 0; i < pSets_filtered.size(); ++i) {
        if (pSets_filtered[i] == pSet_str && scenarios_filtered[i] == scenario) {
            dfx_final.push_back(dfx[i]);
            dfy_final.push_back(dfy[i]);
            colordf_final.push_back(colordf[i]);
        }
    }

}


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
