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
    auto& organTypes = ms->organTypes; // rename
    auto lengths =  ms->segLength();
    auto width = ms->maxBound.minus(ms->minBound);
    std::vector<double> outer_radii = std::vector<double>(ms->segments.size());
    std::fill(outer_radii.begin(), outer_radii.end(), 0.);
    for(auto iter = ms->cell2seg.begin(); iter != ms->cell2seg.end(); ++iter) {
        int cellId =  iter->first;
        auto segs = ms->cell2seg.at(cellId);
        if (cellId >= 0){
				if (vols.size()==0) {
					cellVolume = width.x*width.y*width.z/ms->resolution.x/ms->resolution.y/ms->resolution.z;
				} else {
					cellVolume = vols.at(cellId);
				}
				double v = 0.;  // calculate sum of root volumes or surfaces over cell
				for (int i : segs) {
					if(organTypes[i] == Organism::ot_root)
					{
						if (type==0) { // volume
							v += M_PI*(radii[i]*radii[i])*lengths[i];
						} else if (type==1) { // surface
							v += 2*M_PI*radii[i]*lengths[i];
						} else if (type==2) { // length
							v += lengths[i];
						}
					}
				}
				for (int i : segs) { // calculate outer radius
						if(organTypes[i] == Organism::ot_root)
						{
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
        }else{
            for (int i : segs) {
                outer_radii[i] = 1; // random value
            }
        }
    }
    return outer_radii;
}


/**
 * Sums segment fluxes over each cell
 *
 * @param segFluxes 	segment fluxes given per segment index [cm3/day]
 * @return hash map with cell indices as keys and fluxes as values [cm3/day]
 */
std::map<int,double> Perirhizal::sumSegFluxes(const std::vector<double>& segFluxes)
{
	return this->ms->sumSegFluxes(segFluxes);
}

/**
 * Splits soil fluxes per cell to the segments within the cell, so that the summed fluxes agree, @see sumSoilFluxes()
 *
 * @param soilFluxes 	cell fluxes per global index [cm3/day]
 * @param type 			split flux proportional to 0: segment volume, 1: segment surface, 2: segment length
 * @return fluxes for each segment [cm3/day]
 */
std::vector<double> Perirhizal::splitSoilFluxes(const std::vector<double>& soilFluxes, int type) const {
	return this->ms->splitSoilFluxes(soilFluxes, type);
}

void Perirhizal::redistribute_excess() {
    // Identify and handle excess values
    double extraElement = 0.0;
    for (size_t i = 0; i < val_new.size(); ++i) {
        if (val_new[i] > maxVal) {
            extraElement += (val_new[i] - maxVal) * volumes[i];
            val_new[i] = maxVal;
        }
    }

    // Redistribution loop
    int n_iter = 0;
    std::vector<double> canAdd(val_new.size(), 0.0);
    while (extraElement > 1e-20 && n_iter < 10) {
        n_iter++;
        for (size_t i = 0; i < val_new.size(); ++i) {
            canAdd[i] = std::max((maxVal - val_new[i]) * volumes[i], 0.0);
        }

        double totalCanAdd = std::accumulate(canAdd.begin(), canAdd.end(), 0.0);
        if (totalCanAdd < extraElement) {
            throw std::runtime_error("Not enough capacity to redistribute excess elements.");
        }

        if (divideEqually) {
            double canAddAbove0 = std::count_if(canAdd.begin(), canAdd.end(), [](double val) { return val > 0; });
            double addable_amount = extraElement / canAddAbove0;
            for (size_t i = 0; i < canAdd.size(); ++i) {
                double to_add = std::min(addable_amount, canAdd[i]);
                val_new[i] += to_add / volumes[i];
                extraElement -= to_add;
            }
        } else {
            for (size_t i = 0; i < canAdd.size(); ++i) {
                if (extraElement <= 0.0) break;
                double to_add = std::min(extraElement, canAdd[i]);
                val_new[i] += to_add / volumes[i];
                extraElement -= to_add;
            }
        }
    }

}

void Perirhizal::redistribute_deficit() {
    // Identify and handle deficit values
    double missingElement = 0.0;
    for (size_t i = 0; i < val_new.size(); ++i) {
        if (val_new[i] < minVal) {
            missingElement += (minVal - val_new[i]) * volumes[i];
            val_new[i] = minVal;
        }
    }

            if(verbose)
            {
                std::cout<<"Perirhizal::redistribute_deficit_init "<<std::endl;
                for (size_t i = 0; i < val_new.size(); ++i) {
                    std::cout<<val_new[i] << " ";
                } std::cout<<std::endl;
                std::cout << "missingElement "<<missingElement <<std::endl;
            }
    
    // Redistribution loop
    int n_iter = 0;
    std::vector<double> canTake(val_new.size(), 0.0);
    while (missingElement > 1e-20 && n_iter < 10) {
        n_iter++;
        
        for (size_t i = 0; i < val_new.size(); ++i) {
            canTake[i] = std::max((val_new[i] - minVal) * volumes[i], 0.0);
        }

            if(verbose)
            {
                std::cout<<"Perirhizal::redistribute_deficit_iter "<<n_iter <<std::endl;
                for (size_t i = 0; i < val_new.size(); ++i) {
                    std::cout<<val_new[i] << " ";
                } std::cout<<std::endl;
                
                std::cout << "missingElement "<<missingElement <<std::endl<<"canTake ";
                for (size_t i = 0; i < canTake.size(); ++i) {
                    std::cout<<canTake[i] << " ";
                } std::cout<<std::endl;
            }
        
        double totalCanTake = std::accumulate(canTake.begin(), canTake.end(), 0.0);
        if (totalCanTake < missingElement) {
            throw std::runtime_error("Not enough capacity to redistribute deficit elements.");
        }

        if (divideEqually) {
            double canTakeAbove0 = std::count_if(canTake.begin(), canTake.end(), [](double val) { return val > 0; });
            double removable_amount = missingElement / canTakeAbove0;
            for (size_t i = 0; i < canTake.size(); ++i) {
                double to_remove = std::min(removable_amount, canTake[i]);
                val_new[i] -= to_remove / volumes[i];
                missingElement -= to_remove;
            }
        } else {
            for (size_t i = 0; i < canTake.size(); ++i) {
                if (missingElement <= 0.0) break;
                double to_remove = std::min(missingElement, canTake[i]);
                val_new[i] -= to_remove / volumes[i];
                missingElement -= to_remove;
            }
        }
    }

            if(verbose)
            {
                std::cout<<"Perirhizal::redistribute_deficit_done "<<n_iter <<std::endl;
                for (size_t i = 0; i < val_new.size(); ++i) {
                    std::cout<<val_new[i] << " ";
                } std::cout<<std::endl;
                
                std::cout << "missingElement "<<missingElement <<std::endl<<"canTake ";
                for (size_t i = 0; i < canTake.size(); ++i) {
                    std::cout<<canTake[i] << " ";
                } std::cout<<std::endl;
            }
}

void Perirhizal::handle_excess() {
    double max_val_new = *std::max_element(val_new.begin(), val_new.end());
    if (max_val_new < maxVal) {
        if (maxVal - max_val_new < 1e-20) {
            std::transform(val_new.begin(), val_new.end(), val_new.begin(),
                  [*this](double v) { return std::min(v,this->maxVal); });
        } else {
            redistribute_excess();
        }
    }
}

void Perirhizal::handle_deficit() {
    double min_val_new = *std::min_element(val_new.begin(), val_new.end());
    if (min_val_new < minVal) {
        if (minVal - min_val_new < 1e-20) {
            std::transform(val_new.begin(), val_new.end(), val_new.begin(),
                   [*this](double v) { return std::max(v, this->minVal); });
            if(verbose)
            {
                std::cout<<"Perirhizal::handle_deficit_smalldiff "<<std::endl;
                for (size_t i = 0; i < val_new.size(); ++i) {
                    std::cout<<val_new[i] << " ";
                } std::cout<<std::endl;
            }
        } else {
            if(verbose)
            {
                std::cout<<"Perirhizal::handle_deficit_bigdiff "<<std::endl;
                for (size_t i = 0; i < val_new.size(); ++i) {
                    std::cout<<val_new[i] << " ";
                } std::cout<<std::endl;
            }
            redistribute_deficit();
        }
    }
}

    
/*
    Update the concentration to get specific total content and gradient.
    @param val_new: total concentration of the element in the cylinder cells (cm3/cm3 water or mol/cm3 solute)
    @param minVal: min acceptable content or concentration (float or list)
    @param maxVal: max acceptable content or concentration (float or list)
    @param volumes: volume of each cell of the cylinder (cm3)
    @param divideEqually : divid excess or deficit element equally (True) or incrementally from cell nearest to the root surface outwards (False)
*/
    
std::vector<double> Perirhizal::adapt_values(std::vector<double> val_new_, 
                                             double minVal_, double maxVal_, 
                                             const std::vector<double>& volumes_, 
                                             bool divideEqually_, bool verbose_ = false) {
    val_new = val_new_; volumes = volumes_;
    minVal = minVal_; maxVal = maxVal_;
    divideEqually = divideEqually_;
    verbose = verbose_;
    
    // backup to check mass balance
    double val_newBU = std::inner_product(val_new.begin(), val_new.end(), volumes.begin(), 0.0);
    int n_iter = 0;
    while ((((*std::max_element(val_new.begin(), val_new.end()) > maxVal) && (maxVal > 0.)) || 
            (*std::min_element(val_new.begin(), val_new.end()) < minVal)) && 
            n_iter <= 5) {
        n_iter++;
        if (maxVal > 0.){ handle_excess();}
        handle_deficit();
    }
    double val_newTemp = std::inner_product(val_new.begin(), val_new.end(), 
                                            volumes.begin(), 0.0);
    double diff = val_newTemp - val_newBU;
    if(std::abs(diff) > 1e-16)
    {
        std::cout<<"Perirhizal::adapt_values() "<<diff<<" "<<val_newBU<<" "<< val_newTemp<<std::endl;
        for (size_t i = 0; i < val_new.size(); ++i) {
                std::cout<<val_new[i]<<", ";// concentration to content
            }std::cout<<std::endl;
        throw std::runtime_error("Perirhizal::adapt_values() std::abs(diff) > 1e-16" );
    }
    if (maxVal > 0.){assert(*std::max_element(val_new.begin(), val_new.end()) <= maxVal);}
    assert(*std::min_element(val_new.begin(), val_new.end()) >= minVal);
    return val_new;
}
    
    
        // Function to distribute solute source across cells according to seg_values_content
    std::vector<double> Perirhizal::distributeValSolute_(
        std::vector<double> seg_values_content, 
        const std::vector<double>& volumes, 
        double source, 
        double dt) 
    {
        double total_content = std::accumulate(seg_values_content.begin(), seg_values_content.end(), 0.0);
        source = std::max(-total_content / dt, source);


        // when used for the axyssimetric segments, plant uptake rate added manually at inner cell and can lead to
        // seg_values_content[0] < 0.
        if (*std::min_element(seg_values_content.begin(), seg_values_content.end()) < 0) {
            for (size_t i = 0; i < seg_values_content.size(); ++i) {
                seg_values_content[i] /= volumes[i];// concentration to content
            }
            seg_values_content = this->adapt_values(seg_values_content, // 
                                                    0.0, // min concentration
                                                    -1., //max concentration
                                                    volumes, // cell volumes
                                                    false); // divide equally
            for (size_t i = 0; i < seg_values_content.size(); ++i) {
                seg_values_content[i] *= volumes[i];
            }
        }
        
        std::vector<double> weightVals(seg_values_content.size());
        

        if (total_content == 0.0) {// cannot use content for weighing factor, same val everywhere
            std::fill(weightVals.begin(), weightVals.end(), 1.0 / seg_values_content.size());
        } else if (source < 0) {

            double sum_content = std::accumulate(seg_values_content.begin(), seg_values_content.end(), 0.0);
            for (size_t i = 0; i < seg_values_content.size(); ++i) {
                weightVals[i] = seg_values_content[i] / sum_content;
            }
        } else {
            
            std::replace(seg_values_content.begin(), seg_values_content.end(), 0.0, 1.e-14);// avoid division by 0

            double sum_inv_content = 0.0;
            for (const auto& val : seg_values_content) {
                sum_inv_content += 1.0 / val;
            }

            for (size_t i = 0; i < seg_values_content.size(); ++i) {
                weightVals[i] = (1.0 / seg_values_content[i]) / sum_inv_content;
            }
        }

        for (size_t i = 0; i < weightVals.size(); ++i) {
            weightVals[i] *= source;
        }

        return weightVals;
    }

    // Function to distribute water
    std::vector<double> Perirhizal::distributeValWater_(
        std::vector<double> seg_values_perVol, 
        const std::vector<double>& volumes, 
        double source, 
        double dt, 
        double theta_S, 
        double theta_wilting_point) 
    {
        std::vector<double> availableSpaceOrWater(seg_values_perVol.size());
        double total_space_or_water = 0.0;

        if (source > 0) { // space available
            for (size_t i = 0; i < seg_values_perVol.size(); ++i) {
                availableSpaceOrWater[i] = (theta_S - seg_values_perVol[i]) * volumes[i];
                total_space_or_water += availableSpaceOrWater[i];
            }
            source = std::min(total_space_or_water / dt, source);
        } else { // water extraction
            for (size_t i = 0; i < seg_values_perVol.size(); ++i) {
                availableSpaceOrWater[i] = (seg_values_perVol[i] - theta_wilting_point) * volumes[i];
                total_space_or_water += availableSpaceOrWater[i];
            }
            source = std::max(-total_space_or_water / dt, source);
        }

        if (total_space_or_water > 0) {
            if (*std::min_element(availableSpaceOrWater.begin(), availableSpaceOrWater.end()) < 0) {
                availableSpaceOrWater = adapt_values(availableSpaceOrWater, //
                                                     0.0, //minVal
                                                     -1.0,//maxVal
                                                     volumes, 
                                                     false);//divideEqually
                for (size_t i = 0; i < availableSpaceOrWater.size(); ++i) {
                    availableSpaceOrWater[i] *= volumes[i];//water volume
                }
            }

            double sum_weights = std::accumulate(availableSpaceOrWater.begin(), availableSpaceOrWater.end(), 0.0);
            for (size_t i = 0; i < availableSpaceOrWater.size(); ++i) {
                availableSpaceOrWater[i] = (availableSpaceOrWater[i] / sum_weights) * source;
            }

            return availableSpaceOrWater;
        } else {
            return std::vector<double>(availableSpaceOrWater.size(), 0.0);
        }
    }
    
    std::vector<double> Perirhizal::splitSoilVals(
        const std::vector<double>& soilVals, 
        const std::vector<int>& cellIds, 
        bool isWater, 
        const std::vector<double>& seg_values, 
        const std::vector<double>& seg_volume, 
        double dt, 
        double theta_S,
        double theta_wilting_point) 
    {
        std::vector<int> organTypes = this->ms->organTypes;
        std::vector<double> splitVals(organTypes.size(), 0.0);
        std::vector<double> soilVals_(soilVals.size(), 0.0); // adapted soil vals

        for (int cellid : cellIds) 
        {
            const std::vector<int>& segIds = this->ms->cell2seg.at(cellid);
            std::vector<int> rootIds;
            for (int sid : segIds) 
            {
                if (organTypes[sid] == 2) {rootIds.push_back(sid);}
            }

            if (soilVals[cellid] != 0 && !rootIds.empty()) 
            {
                std::vector<double> seg_values_roots;
                std::vector<double> seg_values_roots_;
                std::vector<double> seg_volume_roots;
                for (int index : rootIds) 
                {
                    seg_values_roots.push_back(seg_values[index]);
                    seg_volume_roots.push_back(seg_volume[index]);
                    seg_values_roots_.push_back(seg_values[index] * seg_volume[index]);
                }
                std::vector<double> splitVals_;
                if (isWater) 
                {
                    splitVals_ = distributeValWater_(
                        seg_values_roots,
                        seg_volume_roots, 
                        soilVals[cellid], 
                        dt,
                        theta_S,
                        theta_wilting_point);
                } else 
                {
                    splitVals_ = distributeValSolute_(
                        seg_values_roots_,
                        seg_volume_roots, 
                        soilVals[cellid], 
                        dt);
                }

                soilVals_[cellid] = std::accumulate(splitVals_.begin(), splitVals_.end(), 0.0);
                for (size_t i = 0; i < rootIds.size(); ++i) {splitVals[rootIds[i]] = splitVals_[i];}
            }

            //_verify_splits(soilVals_, cellIds, splitVals, seg_values);
        }

        return splitVals;
    }
}
