// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "PlantHydraulicParameters.h"

#include <algorithm>
#include <set>

namespace CPlantBox {

/**
 *  give outer water potential [cm] overloaded by @see Photosynthesis::getPsiOut
 * @param cells         sx per cell (true), or segments (false)
 * @param si            segment index
 * @param sx        [cm] soil matric potential for each cell
 */

double PlantHydraulicParameters::getPsiOut(bool cells, int si, const std::vector<double>& sx_, bool verbose) const
{
    int organType = ms->organTypes.at(si);
    double psi_s;
    if (cells) { // soil matric potential given per cell
        int cellIndex = ms->seg2cell.at(si);
        if (cellIndex>=0) {
            if((organType ==Organism::ot_leaf) && verbose){ //add a runtime error?
                std::cout<<"PlantHydraulicParameters::getPsiOut: Leaf segment n#"<<si<<" below ground. OrganType: ";
                std::cout<<organType<<" cell Index: "<<cellIndex<<std::endl;
            }
            if(sx_.size()>1) {
                psi_s = sx_.at(cellIndex);
            } else {
                psi_s = sx_.at(0);
            }
        } else {
            if((organType == Organism::ot_root) && verbose) { //add a runtime error?
                std::cout<<"PlantHydraulicParameters::getPsiOut: Root segment n#"<<si<<" aboveground. OrganType: ";
                std::cout<<organType<<" cell Index: "<<cellIndex<<std::endl;
            }
            psi_s = psi_air;
        }
    } else {
        psi_s = sx_.at(si); // j-1 = segIdx = s.y-1
    }
    return psi_s;
}



/**
 *  Sets the radial conductivity in [1 day-1]
 * TODO: make deprecated: in the examples, replace setKr[Kr] by setKr[[Kr]]
 */
//either age or type/subtype dependent
void PlantHydraulicParameters::setKr(std::vector<double> values, std::vector<double> age, bool verbose) {
    kr =  { values }; //because kr is std::vector<std::vector<double>>
    kr_t = { age };
    if (age.size()==0) {
        if (values.size()==1) {
            kr_f = std::bind(&PlantHydraulicParameters::kr_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
            if (verbose) {
                std::cout << "Kr is constant " << values[0] << " 1 day-1 \n";
            }
        } else {
            kr_f  = std::bind(&PlantHydraulicParameters::kr_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
            if (verbose) {
                std::cout << "Kr is constant per type, type 0 = " << values[0] << " 1 day-1 \n";
            }
        }
    } else {
        kr_f  = std::bind(&PlantHydraulicParameters::kr_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose) {
            std::cout << "Kr is age dependent\n";
        }
    }
}

/**
 *  Sets the radial conductivity in [1 day-1]
 * in case of organ_type specific kr
 * @param values        kr per organ pr/and organ type) or/and per age [cm-1]
 * @param age           ages if kr per age
 * @param kr_length_    exchange zone in root, where kr > 0 [cm from root tip], default = -1.0, i.e., no kr_length
 */
//either age or type/subtype dependent
void PlantHydraulicParameters::setKr(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, double kr_length_, bool verbose) {
    kr = values;
    kr_t = age;
    if (age.size()==0) {
        if (values.size()==1) {
            if (values[0].size()==1) {
                kr_f = std::bind(&PlantHydraulicParameters::kr_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose) {
                    std::cout << "Kr is constant " << values[0][0] << " 1 day-1 \n";
                }
            } else {
                kr_f  = std::bind(&PlantHydraulicParameters::kr_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose) {
                    std::cout << "Kr is constant per subtype, subtype 0 = " << values[0][0] << " 1 day-1 \n";
                }
            }
        } else {
            if (values[0].size()==1) {
                kr_f = std::bind(&PlantHydraulicParameters::kr_perOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose) {
                    std::cout << "Kr is constant per organ type, organ type 2 (root) = " << values[0][0] << " 1 day-1 \n";
                }
            } else {
                if (kr_length_ > 0.){
                    if (verbose) {
                        std::cout << "Exchange zone in roots: kr > 0 until "<< kr_length_<<"cm from root tip"<<std::endl;
                    }
                    ms->kr_length = kr_length_; //in MappedPlant. define distance to root tipe where kr > 0 as cannot compute distance from age in case of carbon-limited growth
                    ms->calcExchangeZoneCoefs();    //computes coefficient used by XylemFlux::kr_RootExchangeZonePerType
                    kr_f  = std::bind(&PlantHydraulicParameters::kr_RootExchangeZonePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                } else {
                    kr_f  = std::bind(&PlantHydraulicParameters::kr_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                }
                if (verbose) {
                    std::cout << "Kr is constant per subtype of organ type, for root, subtype 0 = " << values[0][0] << " 1 day-1 \n";
                }
            }
        }
    } else {
        kr_f  = std::bind(&PlantHydraulicParameters::kr_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose) {
            std::cout << "Kr is equal for all organs and age dependent\n";
        }
    }
}

/**
 *  Sets the axial conductivity in [cm3 day-1]
 * TODO: make deprecated: in the examples, replace setKx[Kx] by setKx[[Kx]]
 */
//either age or type/subtype dependent
void PlantHydraulicParameters::setKx(std::vector<double> values, std::vector<double> age, bool verbose) {
    kx = { values }; // because kx is std::vector<std::vector<double>>
    kx_t = { age };
    if (age.size()==0) {
        if (values.size()==1) {
            kx_f = std::bind(&PlantHydraulicParameters::kx_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
            if(verbose) {
                std::cout << "Kx is constant " << values[0] << " cm3 day-1 \n";
            }
        } else {
            kx_f  = std::bind(&PlantHydraulicParameters::kx_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
            if(verbose) {
                std::cout << "Kx is constant per subtype, subtype 0 = " << values[0] << " cm3 day-1 \n";
            }
        }
    } else {
        kx_f  = std::bind(&PlantHydraulicParameters::kx_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose) {
            std::cout << "Kx is age dependent\n";
        }
    }
}

//either age or type/subtype dependent
void PlantHydraulicParameters::setKx(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, bool verbose)
{
    kx = values;
    kx_t = age;
    if (age.size()==0) {
        if (values.size()==1) {
            if (values[0].size()==1) {
                kx_f = std::bind(&PlantHydraulicParameters::kx_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose) {
                    std::cout << "Kx is constant " << values[0][0] << " cm3 day-1 \n";
                }
            } else {
                kx_f  = std::bind(&PlantHydraulicParameters::kx_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose) {
                    std::cout << "Kx is constant per subtype, subtype 0 = " << values[0][0] << " cm3 day-1 \n";
                }
            }
        } else {
            if (values[0].size()==1) {
                kx_f = std::bind(&PlantHydraulicParameters::kx_perOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose) {
                    std::cout << "Kx is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm3 day-1 \n";
                }
            } else {
                kx_f  = std::bind(&PlantHydraulicParameters::kx_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose) {
                    std::cout << "Kx is constant per subtype of organ type, for root, subtype 0 = " << values[0][0] << " cm3 day-1 \n";
                }
            }
        }
    } else {
        kx_f  = std::bind(&PlantHydraulicParameters::kx_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose) {
            std::cout << "Kx is equal for all organs and age dependent\n";
        }
    }
}

/**
 * Sets the radial conductivity in [1 day-1]
 *  @param values           kr values for age (its linearly interpolated between these values) for each root type
 *  @param age              ages for the given values for each root type
 *
 * TODO: make deprecated (i would leave it in for now, used in pyhton_modules/root_conductivities.py)
 */
//both age and type/subtype dependent
void PlantHydraulicParameters::setKrTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, bool verbose, bool ageBased) {
    krs= { values };
    krs_t = { age };
    kr_f = std::bind(&PlantHydraulicParameters::kr_tablePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    if (verbose) {
        std::cout << "Kr is age dependent per root type\n";
    }
}

/**
 *  Sets the axial conductivity in [cm3 day-1]
 *  @param values           kx values for age (its linearly interpolated between these values) for each root type
 *  @param age              ages for the given values for each root type
 *
 * TODO: make deprecated (i would leave it in for now, used in pyhton_modules/root_conductivities.py)
 */
//both age and type/subtype dependent
void PlantHydraulicParameters::setKxTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, bool verbose) {
    kxs = {values};
    kxs_t = {age};
    kx_f = std::bind(&PlantHydraulicParameters::kx_tablePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    if (verbose) {
        std::cout << "Kx is age dependent per root type\n";
    }
}

/**
 *  Sets the radial conductivity in [1 day-1] age, organ type and sub-type dependent
 *
 *  @param values           kr values for age (its linearly interpolated between these values) for each organ type and each sub-type
 *  @param age              ages for the given values for each organ type and for each sub type
 */
void PlantHydraulicParameters::setKrTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age, bool verbose,  bool ageBased) {
    krs = values;
    krs_t = age;
    if (age[0].size()==1) {
        kr_f = std::bind(&PlantHydraulicParameters::kr_tablePerOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if (verbose) {
            std::cout << "Kr is age dependent per organ type\n";
        }
    } else {
        if (ageBased) {
            kr_f = std::bind(&PlantHydraulicParameters::kr_tablePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        } else {
            ms->kr_length = 100000.; //in MappedPlant. define distance to root tipe where kr > 0 as cannot compute distance from age in case of carbon-limited growth
            ms->calcExchangeZoneCoefs();    //computes coefficient used by XylemFlux::kr_RootExchangeZonePerType
            kr_f  = std::bind(&PlantHydraulicParameters::kr_tablePerType_distance, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        }
        if (verbose) {
            std::cout << "Kr is age dependent per organ type and sub type\n";
        }
    }
}

/**
 *  Sets the axial conductivity in [cm3 day-1] age, organ type and sub-type dependent
 *
 *  @param values           kx values for age (its linearly interpolated between these values) for each organ type and each sub type
 *  @param age              ages for the given values for each organ type and for each sub-type
 */
void PlantHydraulicParameters::setKxTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age, bool verbose) {
    kxs= values;
    kxs_t = age;
    if (age[0].size()==1) {
        kx_f = std::bind(&PlantHydraulicParameters::kx_tablePerOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose) {
            std::cout << "Kx is age dependent per organ type\n";
        }
    } else {
        kx_f = std::bind(&PlantHydraulicParameters::kx_tablePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose) {
            std::cout << "Kx is age dependent per organ type and sub type\n";
        }
    }
}

/**
 * Sets the radial conductivity conductivity [1 day-1] per segment (e.g. constant value per segment)
 */
void PlantHydraulicParameters::setKrValues(std::vector<double> values, bool verbose) {
//    assert(values.size() == rs->segments.size() && "XylemFlux::setKrValues: values size must equal number of segments");
    kr.clear();
    kr_t.clear();
    kr.push_back(values);
    kr_f = std::bind(&PlantHydraulicParameters::kr_valuePerSegment, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    if(verbose) {
        std::cout << "Kr is given per segment\n";
    }
}

/**
 * Sets the axial conductivity [cm3 day-1] per segment (e.g. constant value per segment)
 */
void PlantHydraulicParameters::setKxValues(std::vector<double> values, bool verbose) {
//    assert(values.size() == rs->segments.size() && "XylemFlux::setKxValues: values size must equal number of segments");
    kx.clear();
    kx_t.clear();
    kx.push_back(values);
    kx_f = std::bind(&PlantHydraulicParameters::kx_valuePerSegment, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    if(verbose) {
        std::cout << "Kx is given per segment\n";
    }
}

/**
 * Returns radial conductivities per segment multiplied by segment surface for a specific simulation time (TODO numleaf is ingored)
 */
std::vector<double> PlantHydraulicParameters::getEffKr(double simtime) {
    std::vector<double> kr = std::vector<double>(ms->segments.size());
    for (int si = 0; si<ms->segments.size(); si++) {
        int i = ms->segments[si].x;
        int j = ms->segments[si].y;
        Vector3d n1 = ms->nodes[i];
        Vector3d n2 = ms->nodes[j];
        double l = (n2.minus(n1)).length();
        double a = ms->radii[si];
        int organType = ms->organTypes[si];
        double age = simtime - ms->nodeCTs[j];
        int subType = ms->subTypes[si];
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
std::vector<double> PlantHydraulicParameters::getKr(double simtime) {
    std::vector<double> kr = std::vector<double>(ms->segments.size());
    for (int si = 0; si<ms->segments.size(); si++) {
        int j = ms->segments[si].y;
        int organType = ms->organTypes[si];
        double age = simtime - ms->nodeCTs[j];
        int subType = ms->subTypes[si];
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
std::vector<double> PlantHydraulicParameters::getKx(double simtime) {
    std::vector<double> kx = std::vector<double>(ms->segments.size());
    for (int si = 0; si<ms->segments.size(); si++) {
        int j = ms->segments[si].y;
        int organType = ms->organTypes[si];
        double age = simtime - ms->nodeCTs.at(j);
        int subType = ms->subTypes[si];
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
 * Returns kr of roots belowground or leaves aboveground, overwise returns 0
 */
double PlantHydraulicParameters::kr_f_wrapped(int si, double age, int subType, int organType, bool cells) const
{
    int cellIndex = ms->seg2cell.at(si);
    if (cells&&(((cellIndex>=0)&&(organType !=Organism::ot_root))||((cellIndex < 0)&&(organType ==Organism::ot_root)))) {
        return 0.;
    } else {
        return kr_f(si, age, subType, organType);
    }
}


} // namespace
