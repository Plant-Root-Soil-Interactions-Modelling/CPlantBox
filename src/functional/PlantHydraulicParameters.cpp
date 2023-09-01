// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "PlantHydraulicParameters.h"

#include <algorithm>
#include <set>

namespace CPlantBox {

/**
 * Sets the radial conductivity conductivity [1 day-1] per segment (e.g. constant value per segment)
 */
void PlantHydraulicParameters::setKrValues(std::vector<double> values) {
    kr.clear();
    kr_t.clear();
    kr.push_back(values);
    kr_f = std::bind(&PlantHydraulicParameters::kr_valuePerSegment, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    std::cout << "Kr is given per segment\n";
}

/**
 * Sets the axial conductivity [cm3 day-1] per segment (e.g. constant value per segment)
 */
void PlantHydraulicParameters::setKxValues(std::vector<double> values) {
    kx.clear();
    kx_t.clear();
    kx.push_back(values);
    kx_f = std::bind(&PlantHydraulicParameters::kx_valuePerSegment, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    std::cout << "Kx is given per segment\n";
}

/**
 *  Sets the radial conductivity in [1 day-1] in case of organType or subType specific Kr, or age dependent
 *
 * @param values 		kr per organType and subType or per age [cm-1]
 * @param age 			ages if kr is per age, empty vector otherwise
 * @param kr_length_ 	exchange zone in root, where kr > 0 [cm from root tip], default = -1.0, i.e., no kr_length
 */
void PlantHydraulicParameters::setKr(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, double kr_length_) {
    kr = values;
    kr_t = age;
    if (age.size()==0) {
        if (values.size()==1) {
            if (values[0].size()==1) {
                kr_f = std::bind(&PlantHydraulicParameters::kr_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                std::cout << "Kr is constant " << values[0][0] << " 1 day-1 \n";
            } else {
                kr_f  = std::bind(&PlantHydraulicParameters::kr_perSubType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                std::cout << "Kr is constant per subtype, subtype 0 = " << values[0][0] << " 1 day-1 \n";
            }
        } else {
            if (values[0].size()==1) {
                kr_f = std::bind(&PlantHydraulicParameters::kr_perOrganType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                std::cout << "Kr is constant per organ type, organ type 2 (root) = " << values[0][0] << " 1 day-1 \n";
            } else {
				if(kr_length_ > 0.){
					std::cout << "Exchange zone in roots: kr > 0 until "<< kr_length_<<" cm from root tip"<<std::endl;
//					rs->kr_length = kr_length_; //in MappedPlant. define distance to root tipe where kr > 0 as cannot compute distance from age in case of carbon-limited growth
//					rs->calcExchangeZoneCoefs();	//computes coefficient used by XylemFlux::kr_RootExchangeZonePerType
					kr_f  = std::bind(&PlantHydraulicParameters::kr_RootExchangeZonePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
				} else {
					kr_f  = std::bind(&PlantHydraulicParameters::kr_perSubType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
				}
                std::cout << "Kr is constant per subtype of organ type, for root, subtype 0 = " << values[0][0] << " 1 day-1 \n";
            }
        }
    } else {
        kr_f  = std::bind(&PlantHydraulicParameters::kr_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        std::cout << "Kr is equal for all organs and age dependent\n";
    }
}

/**
 *  Sets the axial conductivity in [cm3 day-1] in case of organType or subType specific Kr, xor age dependent
 *
 * @param values        kx per organType and subType xor per age [cm-1]
 * @param age           ages if kx per age
 */
void PlantHydraulicParameters::setKx(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age) {
    kx = values;
    kx_t = age;
    if (age.size()==0) {
        if (values.size()==1) {
            if (values[0].size()==1) {
                kx_f = std::bind(&PlantHydraulicParameters::kx_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                std::cout << "Kx is constant " << values[0][0] << " cm3 day-1 \n";
            } else {
                kx_f  = std::bind(&PlantHydraulicParameters::kx_perSubType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                std::cout << "Kx is constant per subtype, subtype 0 = " << values[0][0] << " cm3 day-1 \n";
            }
        } else {
            if (values[0].size()==1) {
                kx_f = std::bind(&PlantHydraulicParameters::kx_perOrganType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                std::cout << "Kx is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm3 day-1 \n";
            } else {
                kx_f  = std::bind(&PlantHydraulicParameters::kx_perSubType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                std::cout << "Kx is constant per subtype of organ type, for root, subtype 0 = " << values[0][0] << " cm3 day-1 \n";
            }
        }
    } else {
        kx_f  = std::bind(&PlantHydraulicParameters::kx_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        std::cout << "Kx is equal for all organs and age dependent\n";
    }
}

/**
 *  Sets the radial conductivity in [1 day-1] age, organ type and sub-type dependent
 *
 *	@param values 			kr values for age (its linearly interpolated between these values) for each organ type and each sub-type
 *	@param age 				ages for the given values for each organ type and for each sub type
 */
void PlantHydraulicParameters::setKrTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age) {
    krs = values;
    krs_t = age;
    if (age[0].size()==1) {
        kr_f = std::bind(&PlantHydraulicParameters::kr_tablePerOrganType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        std::cout << "Kr is age dependent per organ type\n";
    }
    else{
        kr_f = std::bind(&PlantHydraulicParameters::kr_tablePerSubType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        std::cout << "Kr is age dependent per organ type and sub type\n";
    }
}

/**
 *  Sets the axial conductivity in [cm3 day-1] age, organ type and sub-type dependent
 *
 *	@param values 			kx values for age (its linearly interpolated between these values) for each organ type and each sub type
 *	@param age 				ages for the given values for each organ type and for each sub-type
 */
void PlantHydraulicParameters::setKxTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age) {
    kxs= values;
    kxs_t = age;
    if (age[0].size()==1) {
        kx_f = std::bind(&PlantHydraulicParameters::kx_tablePerOrganType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        std::cout << "Kx is age dependent per organ type\n";
    }
    else {
        kx_f = std::bind(&PlantHydraulicParameters::kx_tablePerSubType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        std::cout << "Kx is age dependent per organ type and sub type\n";
    }
}

/**
 * Returns radial conductivities per segment multiplied by segment surface for a specific simulation time (TODO numleaf is ingored)
 */
std::vector<double> PlantHydraulicParameters::getEffKr(std::shared_ptr<MappedSegments> rs, double simtime) const {
    std::vector<double> kr = std::vector<double>(rs->segments.size());
    for (int si = 0; si<rs->segments.size(); si++) {
        int i = rs->segments[si].x;
        int j = rs->segments[si].y;
        Vector3d n1 = rs->nodes[i];
        Vector3d n2 = rs->nodes[j];
        double l = (n2.minus(n1)).length();
        double a = rs->radii[si];
        int organType = rs->organTypes[si];
        double age = simtime - rs->nodeCTs[j];
        int subType = rs->subTypes[si];
        try {
            kr[si] = 2.*M_PI *a*l*kr_f(si, age, subType, organType);
        } catch(...) {
            std::cout << "\n XylemFlux::segFluxes: radial conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }
    }
    return kr;
}

/**
 * Returns radial conductivities per segment multiplied by segment surface for a specific simulation time (TODO numleaf is ingored)
 */
std::vector<double> PlantHydraulicParameters::getKr(std::shared_ptr<MappedSegments> rs, double simtime) const {
    std::vector<double> kr = std::vector<double>(rs->segments.size());
    for (int si = 0; si<rs->segments.size(); si++) {
        int j = rs->segments[si].y;
        int organType = rs->organTypes[si];
        double age = simtime - rs->nodeCTs[j];
        int subType = rs->subTypes[si];
        try {
            kr[si] = kr_f(si, age, subType, organType);
        } catch(...) {
            std::cout << "\n XylemFlux::segFluxes: radial conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }
    }
    return kr;
}


/**
 * Returns radial conductivities per segment for a specific simulation time
 */
std::vector<double> PlantHydraulicParameters::getKx(std::shared_ptr<MappedSegments> rs, double simtime) const {
    std::vector<double> kx = std::vector<double>(rs->segments.size());
    for (int si = 0; si<rs->segments.size(); si++) {
        int j = rs->segments[si].y;
        int organType = rs->organTypes[si];
        double age = simtime - rs->nodeCTs[j];
        int subType = rs->subTypes[si];
        try {
            kx[si] = kx_f(si, age, subType, organType);
        } catch(...) {
            std::cout << "\n XylemFlux::segFluxes: axial conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }
    }
    return kx;
}

/**
 * Returns soil matric potential per segment, for a given soil sx connected gy the mapper rs->seg2cell
 */
std::vector<double> PlantHydraulicParameters::getHs(std::shared_ptr<MappedSegments> rs, const std::vector<double>& sx) const {
    std::vector<double> hs = std::vector<double>(rs->segments.size());
    for (int si = 0; si<rs->segments.size(); si++) {
        int cellIndex = rs->seg2cell[si];
        if (cellIndex>=0) {
            if(sx.size()>1) {
                hs[si] = sx.at(cellIndex);
            } else {
                hs[si] = sx.at(0);
            }
        } else {
            hs[si] = psi_air;
        }
    }
    return hs;
}


} // namespace
