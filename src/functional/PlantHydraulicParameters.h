// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef PLANT_HYDRAULIC_PARAMETERS
#define PLANT_HYDRAULIC_PARAMETERS

#include "MappedOrganism.h"

#include <map>

namespace CPlantBox {

/**
 * Plant hydraulic conductivities, i.e. axial [cm3/day] and radial [cm/day] conductivities
 *
 * defined by the functions:
 * kr_f(int segment_index, double age, int subType, int organType)
 * kx_f(int segment_index, double age, int subType, int organType)
 *
 * The parameters can be given per organType and subType:
 * - constant
 * - per segment
 * - age dependent
 * - distance dependent
 *
 * TODO needs testing and documentation
 */


class PlantHydraulicParameters
{

public:

    PlantHydraulicParameters();
    PlantHydraulicParameters(std::shared_ptr<CPlantBox::MappedSegments> ms_);

    virtual ~PlantHydraulicParameters() { }

    void setMode(std::string krMode, std::string kxMode); // sets call back functions according to krMode, and kxMode

    std::function<double(int, double, int, int)> kr_f = [](int si, double age, int type, int orgtype) { throw std::runtime_error("kr_f not implemented"); return 0.; };
    std::function<double(int, double,int,int)> kx_f = [](int si, double age, int type, int orgtype) { throw std::runtime_error("kx_f not implemented"); return 1.; };
    double kr_f_wrapped(int si, double age, int type, int orgtype, bool cells) const; ///stops transpiration if organs are not in the correct domain

    void setKrConst(double v, int subType, int organType = Organism::ot_root, double kr_length = -1);
    void setKxConst(double v, int subType, int organType = Organism::ot_root);
    void setKrAgeDependent(std::vector<double> age, std::vector<double> values, int subType, int organType = Organism::ot_root);
    void setKxAgeDependent(std::vector<double> age, std::vector<double> values, int subType, int organType = Organism::ot_root);
    void setKrDistanceDependent(std::vector<double> distance, std::vector<double> values, int subType, int organType = Organism::ot_root);
    void setKxDistanceDependent(std::vector<double> distance, std::vector<double> values, int subType, int organType = Organism::ot_root);
    void setKrValues(std::vector<double> values); ///< one value per segment
    void setKxValues(std::vector<double> values); ///< one value per segment

    std::vector<double> getKr(double simtime); // per segment [1 day-1]
    std::vector<double> getEffKr(double simtime); // per segment [cm2 day-1]
    std::vector<double> getKx(double simtime); // per segment [cm3 day-1]

    virtual double getPsiOut(bool cells, int si, const std::vector<double>& sx_, bool verbose) const; ///< get the outer water potential [cm]

    std::string krMode = "const"; // "const", "perSegment", "age", or "distance"
    std::string kxMode = "const"; // "const", "perSegment", "age", or "distance"

    std::vector<double> krValues = {0.}; // per segment [1 day-1]
    std::vector<double> kxValues = {0.}; // per segment [cm3 day-1]
    std::map<int, std::vector<std::vector<double>>> kr_ages, kr_values; // per organType, per subType [1 day-1]
    std::map<int, std::vector<std::vector<double>>> kx_ages, kx_values; // per organType, per subType [cm3 day-1]

    std::shared_ptr<CPlantBox::MappedSegments> ms; // Needs to be set for kr_RootExchangeZonePerType, and kr_tablePerType_distance

    double psi_air = -954378; // air water potential [cm] for T = 20°C and RH = 0.5

    const int maxSubTypes = 10;

protected:

    double kr_perType(int si, double age, int subType, int organType) { return kr_values.at(organType).at(subType).at(0); }
    double kr_RootExchangeZonePerType(int si,double age, int subType, int organType);
    double kx_perType(int si, double age, int subType, int organType) { return kx_values.at(organType).at(subType).at(0); }

    double kr_age(int si, double age, int subType, int organType) { return Function::interp1(age, kr_ages.at(organType).at(subType), kr_values.at(organType).at(subType)); }
    double kx_age(int si, double age, int subType, int organType) { return Function::interp1(age, kx_ages.at(organType).at(subType), kx_values.at(organType).at(subType)); }
    double kr_distance(int si, double age, int subType, int organType);
    double kx_distance(int si, double age, int subType, int organType);

    double kr_valuePerSegment(int si, double age, int subType, int organType) { return krValues.at(si); }
    double kx_valuePerSegment(int si, double age, int subType, int organType) { return kxValues.at(si); };

};

} // namespace

#endif
