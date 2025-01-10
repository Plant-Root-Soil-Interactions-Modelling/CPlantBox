// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "PlantHydraulicParameters.h"

#include <algorithm>
#include <set>

namespace CPlantBox {

/**
 * Constructor
 */
PlantHydraulicParameters::PlantHydraulicParameters() {

    std::vector<int> organTypes = { Organism::ot_root, Organism::ot_stem, Organism::ot_leaf };
    for (int ot : organTypes) {
        std::vector<double> v0 = {0.};
        std::vector<double> v1 = {1.};
        std::vector<double> v1e_4= {1.e-4};

        // all ages are set to 0.
        std::vector<std::vector<double>> ages;
        ages.resize(maxSubTypes);
        std::fill(ages.begin(), ages.end(), v0);
        kr_ages[ot] = ages;
        kx_ages[ot] = ages;

        // kr is set to 1.e-4 for roots, 0. for stem, and 0. for leaf
        std::vector<std::vector<double>> values;
        values.resize(maxSubTypes);
        if (ot == Organism::ot_root) {
            std::fill(values.begin(), values.end(), v1e_4);
        }
        if (ot == Organism::ot_stem) {
            std::fill(values.begin(), values.end(), v0);
        }
        if (ot == Organism::ot_leaf) {
            std::fill(values.begin(), values.end(), v0);
        }
        kr_values[ot] = values;

        // kx is set to 1.e-4 for roots, 1. for stem, and 1. for leaf
        if (ot == Organism::ot_root) {
            std::fill(values.begin(), values.end(), v1e_4);
        }
        if (ot == Organism::ot_stem) {
            std::fill(values.begin(), values.end(), v1);
        }
        if (ot == Organism::ot_leaf) {
            std::fill(values.begin(), values.end(), v1);
        }
        kx_values[ot] = values;
    }
}

/**
 *  Constructor passing MappedSegments
 */
PlantHydraulicParameters::PlantHydraulicParameters(std::shared_ptr<CPlantBox::MappedSegments> ms_)
: PlantHydraulicParameters() {
    ms = ms_;
}

/**
 *  Sets the call back function according to krMode and kxMode. This is called by the Python JSON Reader,
 *  otherwise krMode and kxMode is set automatically by the setter methods.
 *
 *  krMode is const, perSegment, age, or distance
 *  kxMode is const, perSegment, age, or distance
 */
void PlantHydraulicParameters::setMode(std::string krMode, std::string kxMode) {
    this->krMode = krMode;
    this->kxMode = kxMode;
    if (krMode == "const") {
        kr_f  = std::bind(&PlantHydraulicParameters::kr_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    }
    if (kxMode == "const") {
        kx_f  = std::bind(&PlantHydraulicParameters::kx_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    }
    if (krMode == "perSegment") {
        kr_f = std::bind(&PlantHydraulicParameters::kr_valuePerSegment, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    }
    if (kxMode == "perSegment") {
        kx_f = std::bind(&PlantHydraulicParameters::kx_valuePerSegment, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    }
    if (krMode == "age") {
        kr_f = std::bind(&PlantHydraulicParameters::kr_age, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    }
    if (kxMode == "age") {
        kx_f = std::bind(&PlantHydraulicParameters::kx_age, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    }
    if (krMode == "distance") {
        kr_f = std::bind(&PlantHydraulicParameters::kr_distance, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    }
    if (kxMode == "distance") {
        kx_f = std::bind(&PlantHydraulicParameters::kx_distance, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    }
}

/**
 * Returns the radial conductivity kr of roots below ground or leaves aboveground, otherwise returns 0 (ms must be set)
 */
double PlantHydraulicParameters::kr_f_wrapped(int si, double age, int subType, int organType, bool cells) const
{
    int cellIndex = ms->seg2cell.at(si);
    if (cells&&(((cellIndex>=0)&&(organType !=Organism::ot_root))||((cellIndex < 0)&&(organType == Organism::ot_root)))) {
        return 0.;
    } else {
        return kr_f(si, age, subType, organType);
    }
}

/**
 * Sets the radial conductivity in [1 day-1] for a subType and organType
 *
 * @param v             kr per subType and organType [cm-1]
 * @param subType       organ sub type
 * @param organType     organ type (default is root)
 * @param kr_length     exchange zone in root, where kr > 0 [cm from root tip], default = -1.0, i.e., no kr_length
 */
void PlantHydraulicParameters::setKrConst(double v, int subType, int organType, double kr_length) {
    kr_values.at(organType).at(subType).at(0) = v;
    if (kr_length > 0.){
        ms->kr_length = kr_length; //in MappedPlant. define distance to root tipe where kr > 0 as cannot compute distance from age in case of carbon-limited growth
        ms->calcExchangeZoneCoefs();    //computes coefficient used by XylemFlux::kr_RootExchangeZonePerType
        kr_f  = std::bind(&PlantHydraulicParameters::kr_RootExchangeZonePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        krMode = "constExchangeZone";
    } else {
        kr_f  = std::bind(&PlantHydraulicParameters::kr_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        krMode = "const";
    }
}

/**
 * For roots the radial conductivity is multiplied by an exchange zone coefficient (that is calculated within MappedPlant per Segment)
 */
double PlantHydraulicParameters::kr_RootExchangeZonePerType(int si,double age, int subType, int organType) {
    if (organType == Organism::ot_root){
        double coef = ms->exchangeZoneCoefs.at(si);//% of segment length in the root exchange zone, see MappedPlant::simulate
        return coef * kr_values.at(organType).at(subType).at(0);
    }
    return kr_values.at(organType).at(subType).at(0);
}

/**
 * Sets the axial conductivity in [cm3 day-1] for a subType and organType
 *
 * @param v             kx per subType and organType [cm3/day]
 * @param subType       organ sub type
 * @param organType     organ type (default is root)
 */
void PlantHydraulicParameters::setKxConst(double v, int subType, int organType) {
    kx_values.at(organType).at(subType).at(0) = v;
    kx_f  = std::bind(&PlantHydraulicParameters::kx_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    kxMode = "const";
}

/**
 *  Sets the radial conductivity [day-1] in an age dependent way, organ type and sub-type dependent
 *
 * The value is linearly interpolated between the sampling points given by age, and values dependent  on age,
 * Extrapolation of the last value is used if age lies outside the domain defined by the sampling points.
 *
 * @param age           age values
 * @param v             kx values
 * @param subType       organ sub type
 * @param organType     organ type (default is root)
 */
void PlantHydraulicParameters::setKrAge(std::vector<double> age, std::vector<double> values, int subType, int organType) {
    kr_ages.at(organType).at(subType) = age;
    kr_values.at(organType).at(subType) = values;
    kr_f = std::bind(&PlantHydraulicParameters::kr_age, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    krMode = "age";
}

/**
 *  Sets the axial conductivity [cm3 day-1] in an age dependent way, organ type and sub-type dependent
 *
 * The value is linearly interpolated between the sampling points given by age, and values dependent  on age,
 * Extrapolation of the last value is used if age lies outside the domain defined by the sampling points.  *
 * @param age              age values
 * @param v                kx values
 * @param subType       organ sub type
 * @param organType     organ type (default is root)
 */
void PlantHydraulicParameters::setKxAge(std::vector<double> age, std::vector<double> values, int subType, int organType) {
    kx_ages.at(organType).at(subType) = age;
    kx_values.at(organType).at(subType) = values;
    kx_f = std::bind(&PlantHydraulicParameters::kx_age, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    kxMode = "age";
}

/**
 *  Sets the radial conductivity [day-1] in an age dependent way, organ type and sub-type dependent
 *
 * The value is linearly interpolated between the sampling points given by age, and values dependent  on age,
 * Extrapolation of the last value is used if age lies outside the domain defined by the sampling points.
 *
 * @param age           age values
 * @param v             kx values
 * @param subType       organ sub type
 * @param organType     organ type (default is root)
 */
void PlantHydraulicParameters::setKrDistance(std::vector<double> distance, std::vector<double> values, int subType, int organType) {
    kr_ages.at(organType).at(subType) = distance;
    kr_values.at(organType).at(subType) = values;
    kr_f = std::bind(&PlantHydraulicParameters::kr_distance, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    krMode = "distance";
}

/**
 * Radial conductivity per distance (distance is obtained from MappedPlant via the segment index) for roots,
 * constant values for other organs
 */
double PlantHydraulicParameters::kr_distance(int si, double age, int subType, int organType) {
    if (organType == Organism::ot_root){
        double distFromTip = ms->distanceTip.at(si);
        return Function::interp1(distFromTip, kr_ages.at(organType).at(subType), kr_values.at(organType).at(subType));
    } else {
        return kr_values.at(organType).at(subType).at(0);
    }
}

/**
 *  Sets the axial conductivity [cm3 day-1] in an age dependent way, organ type and sub-type dependent
 *
 * The value is linearly interpolated between the sampling points given by age, and values dependent  on age,
 * Extrapolation of the last value is used if age lies outside the domain defined by the sampling points.  *
 * @param age              age values
 * @param v                kx values
 * @param subType       organ sub type
 * @param organType     organ type (default is root)
 */
void PlantHydraulicParameters::setKxDistance(std::vector<double> distance, std::vector<double> values, int subType, int organType) {
    kx_ages.at(organType).at(subType) = distance;
    kx_values.at(organType).at(subType) = values;
    kx_f = std::bind(&PlantHydraulicParameters::kx_distance, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    kxMode = "distance";
}

/**
 * Axial conductivity per distance (distance is obtained from MappedPlant via the segment index) for roots,
 * constant values for other organs
 */
double PlantHydraulicParameters::kx_distance(int si, double age, int subType, int organType) {
    if (organType == Organism::ot_root){
        double distFromTip = ms->distanceTip.at(si);
        return Function::interp1(distFromTip, kx_ages.at(organType).at(subType), kx_values.at(organType).at(subType));
    } else {
        return kx_values.at(organType).at(subType).at(0);
    }
}

/**
 * Sets the radial conductivity conductivity [1 day-1] per segment (e.g. constant value per segment)
 */
void PlantHydraulicParameters::setKrValues(std::vector<double> values) {
    krValues = values;
    kr_f = std::bind(&PlantHydraulicParameters::kr_valuePerSegment, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    krMode = "perSegment";
}

/**
 * Sets the axial conductivity [cm3 day-1] per segment (e.g. constant value per segment)
 */
void PlantHydraulicParameters::setKxValues(std::vector<double> values) {
    kxValues = values;
    kx_f = std::bind(&PlantHydraulicParameters::kx_valuePerSegment, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    kxMode = "perSegment";
}

/**
 * Returns radial conductivities per segment (TODO numleaf is ingored) [1 day-1]
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
 * Returns radial conductivities per segment multiplied by segment surface for a specific simulation time (TODO numleaf is ingored) [cm2 day-1]
 */
std::vector<double> PlantHydraulicParameters::getEffKr(double simtime) {
    std::vector<double> kr = std::vector<double>(ms->segments.size());
    for (int si = 0; si<ms->segments.size(); si++) {
        int i = ms->segments[si].x;
        int j = ms->segments[si].y;
        Vector3d n1 = ms->nodes[i];
        Vector3d n2 = ms->nodes[j];
        double l = (n2.minus(n1)).length();
        double a = ms->getEffectiveRadius(si); ////// hair_zone_length, hair_eff_length
        int organType = ms->organTypes[si];
        double age = simtime - ms->nodeCTs[j];
        int subType = ms->subTypes[si];
        try {
            kr[si] = 2.*M_PI *a*l*kr_f(si, age, subType, organType);
        } catch(...) {
            std::cout << "\n PlantHydraulicParameters::getEffKr: radial conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }
    }
    return kr;
}

/**
 * Returns axial conductivities per segment for a specific simulation time [cm3 day-1]
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
 * Returns outer water potential [cm] overloaded by @see Photosynthesis::getPsiOut
 * @param cells         sx per cell (true), or segments (false)
 * @param si            segment index
 * @param sx_        [cm] soil matric potential for each cell
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


} // namespace
