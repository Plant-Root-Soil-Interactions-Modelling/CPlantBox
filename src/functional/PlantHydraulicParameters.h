// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef PLANT_HYDRAULIC_PARAMETERS
#define PLANT_HYDRAULIC_PARAMETERS

#include "MappedOrganism.h"

namespace CPlantBox {

/**
 * Plant hydraulic conductivities
 */
class PlantHydraulicParameters
{

public:

    PlantHydraulicParameters() { }

    virtual ~PlantHydraulicParameters() { }

    void setKrValues(std::vector<double> values); ///< one value per segment
    void setKxValues(std::vector<double> values); ///< one value per segment
    void setKr(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, double kr_length_ = -1.0);
    void setKx(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age);

    void setKrTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age);
    void setKxTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age);
    void setKrTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age);
    void setKxTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age);



    std::function<double(int, double, int, int)> kr_f = [](int si, double age, int subType, int organType) { throw std::runtime_error("kr_f not set"); return 0.; };
    std::function<double(int, double, int, int)> kx_f = [](int si, double age, int subType, int organType) { throw std::runtime_error("kx_f not set"); return 1.; };
    std::vector<std::vector<double>> kr, kr_t; //  [1 day-1]
    std::vector<std::vector<double>> kx, kx_t; // [cm3 day-1]
    std::vector<std::vector<std::vector<double>>> krs, krs_t; // [1 day-1]
    std::vector<std::vector<std::vector<double>>> kxs, kxs_t; // [cm3 day-1]

    std::vector<double> getEffKr(std::shared_ptr<MappedSegments> rs, double simtime) const;
    std::vector<double> getKr(std::shared_ptr<MappedSegments> rs, double simtime) const; // [1 day-1]
    std::vector<double> getKx(std::shared_ptr<MappedSegments> rs, double simtime) const; // [cm3 day-1]
    std::vector<double> getHs(std::shared_ptr<MappedSegments> rs, const std::vector<double>& sx) const;

    // additional parameters
    double psi_air = -954378; // air water potential [cm] for T = 20°C and RH = 0.5
    std::vector<double> exchangeZoneCoefs; // TODO put those values here

    // TODO what would be a good method to store and load this class

protected:

    double kr_valuePerSegment(int si, double age, int subType, int organType) { return kr.at(0).at(si); }
    double kr_const(int si,double age, int subType, int organType) { return kr.at(0).at(0); }
	double kr_perOrganType(int si,double age, int subType, int organType) { return kr.at(organType - 2).at(0); } //per organ subType (goes from 2 (root) to 4 (leaf))
    double kr_perSubType(int si,double age, int subType, int organType) { return kr.at(organType - 2).at(subType);} //per subType and organType (goes from 2 (root) to 4 (leaf))
    double kr_table(int si,double age, int subType, int organType) { return Function::interp1(age, kr_t.at(0), kr.at(0)); } //constant for all subTypes and organTypes BUT age dependent
	double kr_tablePerOrganType(int si,double age, int subType, int organType) { return Function::interp1(age, krs_t.at(organType-2).at(0), krs.at(organType-2).at(0));} //constant for all subTypes but organType and age dependent
	double kr_tablePerSubType(int si,double age, int subType, int organType) { return Function::interp1(age, krs_t.at(organType-2).at(subType), krs.at(organType-2).at(subType)); } //subtype, subType and age dependant
	double kr_RootExchangeZonePerType(int si,double age, int subType, int organType) { //when use carbon- and water-limited growth, canNOT use "kr_tablePerType" instead of this function
		if (organType == Organism::ot_root) {
			double coef = exchangeZoneCoefs.at(si);//% of segment length in the root exchange zone, see MappedPlant::simulate
			return coef * kr.at(organType - 2).at(subType);
		}
		return kr.at(organType - 2).at(subType);
	} //per subType and organType, for roots depend on distance to tip for roots

	double kx_valuePerSegment(int si, double age, int subType, int organType) { return kx.at(0).at(si); };
    double kx_const(int si,double age, int subType, int organType) { return kx.at(0).at(0); }
    double kx_perOrganType(int si,double age, int subType, int organType) { return kx.at(organType - 2)[0]; } //per organ subType (goes from 2 (root) to 4 (leaf))
    double kx_perSubType(int si,double age, int subType, int organType) { return kx.at(organType - 2).at(subType); } //per subtype and organ subType (goes from 2 (root) to 4 (leaf))
	double kx_table(int si,double age, int subType, int organType) { return Function::interp1(age, kx_t[0], kx[0]); } //constant for all subType/subtype and age dependant
    double kx_tablePerOrganType(int si,double age, int subType, int organType) { return Function::interp1(age, kxs_t.at(organType-2).at(0), kxs.at(organType-2).at(0)); } //constant for all subtype but subType and age dependant
    double kx_tablePerSubType(int si,double age, int subType, int organType) { return Function::interp1(age, kxs_t.at(organType-2).at(subType), kxs.at(organType-2).at(subType)); } //subtype, subType and age dependant

};

} // namespace

#endif
