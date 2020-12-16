// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef XYLEM_FLUX_H_
#define XYLEM_FLUX_H_

#include "MappedOrganism.h"

namespace CPlantBox {

/**
 * Hybrid solver (Meunier et al. 2017),
 * see also xylem_flux.py in CPlantBox/src/python_modules
 *
 * Units are [cm] and [day]
 *
 * Wraps a MappedSegments class (i.e. MappedRootSystem)
 */
class XylemFlux
{
public:

    XylemFlux(std::shared_ptr<CPlantBox::MappedSegments> rs);

    virtual ~XylemFlux() { }

    void linearSystem(double simTime, const std::vector<double>& sx, bool cells = true,
        const std::vector<double> soil_k = std::vector<double>()); ///< builds linear system (simTime is needed for age dependent conductivities)
    std::map<int,double> soilFluxes(double simTime, const std::vector<double>& rx, const std::vector<double>& sx,
    		bool approx = false, const std::vector<double> soil_k = std::vector<double>()); // [cm3/day]
    std::vector<double> segFluxes(double simTime, const std::vector<double>& rx, const std::vector<double>& sx,
    		bool approx = false, bool cells = false, const std::vector<double> soil_k = std::vector<double>()); // for each segment in [cm3/day]
    std::map<int,double> sumSegFluxes(const std::vector<double>& segFluxes); ///< sums segment fluxes over soil cells,  soilFluxes = sumSegFluxes(segFluxes), [cm3/day]

    std::vector<double> segSRA(double simTime, const std::vector<double>& rx, const std::vector<double>& sx, double wilting_point,
        std::function<double(double)> mfp, std::function<double(double)> imfp);
    std::vector<double> segSRAStressedFlux(const std::vector<double>& sx, double wiltingPoint, double hc,
        std::function<double(double)> mfp, std::function<double(double)> imfp, double dx = 1.e-6);
    std::vector<double> segSRAStressedAnalyticalFlux(const std::vector<double>& sx, std::function<double(double)> mfp);

    std::vector<double> segOuterRadii(int type = 0) const; ///< outer cylinder radii to match cell volume
    std::vector<double> splitSoilFluxes(const std::vector<double>& soilFluxes, int type = 0) const; ///< splits soil fluxes (per cell) into segment fluxes
    std::vector<double> segLength() const; ///< calculates segment lengths [cm]

    std::vector<int> aI; // to assemble the sparse matrix on the Python side
    std::vector<int> aJ;
    std::vector<double> aV;
    std::vector<double> aB;

    void setKr(std::vector<double> values, std::vector<double> age = std::vector<double>(0)); ///< sets a callback for kr:=kr(age,type),  [1 day-1]
    void setKx(std::vector<double> values, std::vector<double> age = std::vector<double>(0)); ///< sets a callback for kx:=kx(age,type),  [cm3 day-1]
    void setKrTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age);
    void setKxTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age);
    void setKr(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age); ///< sets a callback for kr:=kr(age,type),  [1 day-1]
    void setKx(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age); ///< sets a callback for kx:=kx(age,type),  [cm3 day-1]
    void setKrTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age);
    void setKxTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age);

    std::function<double(double,int, int, int)> kr_f = [](double age, int type, int orgtype, int numleaf) { return 0.; };
    std::function<double(double,int,int)> kx_f = [](double age, int type, int orgtype) { return 1.; };

    std::shared_ptr<CPlantBox::MappedSegments> rs;

    std::vector<std::vector<double>> kr, kr_t;
    std::vector<std::vector<double>> kx, kx_t;
    std::vector<std::vector<std::vector<double>>> krs, krs_t;
    std::vector<std::vector<std::vector<double>>> kxs, kxs_t;

    double airPressure = -1000; // static air pressure
	std::vector<double> gs;

protected:

	//type correspond to subtype or to the leaf segment number
    double kr_const(double age, int type, int orgtype, int numleaf) { if (orgtype == 4 && gs.size() > 0 ){return gs.at(numleaf);} else {return kr.at(0).at(0); } }//k constant
    double kr_perOrgType(double age, int type, int orgtype, int numleaf) { if (orgtype == 4&& gs.size() > 0 ) {return gs.at(numleaf);} else { return kr.at(orgtype - 2).at(0); }} //per organ type (goes from 2 (root) to 4 (leaf))
    double kr_perType(double age, int type, int orgtype, int numleaf) { if (orgtype == 4&& gs.size() > 0 ) {return gs.at(numleaf);} else { return kr.at(orgtype - 2).at(type); }}//per subtype and organ type (goes from 2 (root) to 4 (leaf))
    double kr_table(double age, int type, int orgtype, int numleaf) { if (orgtype == 4&& gs.size() > 0 ) {return gs.at(numleaf);} else { return interp1(age, kr_t.at(0), kr.at(0)); }} //constant for all type/subtype and age dependant
	double kr_tablePerOrgType(double age, int type, int orgtype, int numleaf){ if (orgtype == 4&& gs.size() > 0 ) {return gs.at(numleaf);} else  { return interp1(age, krs_t.at(orgtype-2).at(0), krs.at(orgtype-2).at(0)); } }//constant for all subtype but type and age dependant
	double kr_tablePerType(double age, int type, int orgtype, int numleaf) { if (orgtype == 4&& gs.size() > 0 ) {return gs.at(numleaf);} else { return interp1(age, krs_t.at(orgtype-2).at(type), krs.at(orgtype-2).at(type)); }} //subtype, type and age dependant

    double kx_const(double age, int type, int orgtype) { return kx.at(0).at(0); } //k constant
    double kx_perOrgType(double age, int type, int orgtype) { return kx.at(orgtype - 2)[0]; } //per organ type (goes from 2 (root) to 4 (leaf))
    double kx_perType(double age, int type, int orgtype) { return kx.at(orgtype - 2).at(type); } //per subtype and organ type (goes from 2 (root) to 4 (leaf))
	double kx_table(double age, int type, int orgtype) { return interp1(age, kx_t[0], kx[0]); } //constant for all type/subtype and age dependant
    double kx_tablePerOrgType(double age, int type, int orgtype) { return interp1(age, kxs_t.at(orgtype-2).at(0), kxs.at(orgtype-2).at(0)); } //constant for all subtype but type and age dependant
    double kx_tablePerType(double age, int type, int orgtype) { return interp1(age, kxs_t.at(orgtype-2).at(type), kxs.at(orgtype-2).at(type)); } //subtype, type and age dependant

    static double interp1(double ip, std::vector<double> x, std::vector<double> y);

    static double schroederStress(double r, double p, double q_out, double r_in, double r_out, std::function<double(double)> mfp, std::function<double(double)> imfp);

};

} // namespace

#endif
