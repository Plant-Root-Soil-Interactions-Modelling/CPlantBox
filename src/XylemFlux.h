// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef XYLEM_FLUX_H_
#define XYLEM_FLUX_H_

#include "MappedOrganism.h"

namespace CPlantBox {

/**
 * Hybrid solver (Meunier et al. 2017)
 *
 * Units are [cm], [g], and [day]
 *
 * Wraps a MappedSegments class (i.e. MappedRootSystem)
 */
class XylemFlux
{
public:

    XylemFlux(std::shared_ptr<CPlantBox::MappedSegments> rs) :rs(rs) { }

    virtual ~XylemFlux() { }

    void linearSystem(double simTime, const std::vector<double>& sx, bool cells = true,
    		const std::vector<double> soil_k = std::vector<double>()); ///< builds linear system (simTime is needed for age dependent conductivities)

    std::map<int,double> soilFluxes(double simTime, const std::vector<double>& rx, const std::vector<double>& sx, bool approx = false); // [cm3/day]
    std::vector<double> segFluxes(double simTime, const std::vector<double>& rx, const std::vector<double>& sx, bool approx = false); // for each segment in [cm3/day]
    std::vector<double> segFluxesSchroeder(double simTime, std::vector<double> rx, const std::vector<double>& sx, double critP, std::function<double(double)> mpf);

    std::vector<double> segOuterRadii(int type = 0) const; ///< outer cylinder radii to match cell volume
    std::map<int,double> sumSoilFluxes(const std::vector<double>& segFluxes); ///< sums segment fluxes over soil cells,  soilFluxes = sumSoilFluxes(segFluxes), [cm3/day]
    std::vector<double> splitSoilFluxes(const std::vector<double>& soilFluxes, int type = 0) const; ///< splits soil fluxes (per cell) into segment fluxes
    std::vector<double> segLength() const; ///< calculates segment lengths

    std::vector<int> aI; // to assemble the sparse matrix on the Python side
    std::vector<int> aJ;
    std::vector<double> aV;
    std::vector<double> aB;

    void setKr(std::vector<double> values, std::vector<double> age); ///< sets a callback for kr:=kr(age,type),  [1 day-1]
    void setKx(std::vector<double> values, std::vector<double> age); ///< sets a callback for kx:=kx(age,type),  [cm3 day-1]
    void setKrTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age);
    void setKxTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age);

    std::function<double(double,int)> kr_f = [](double age, int type) { return 0.; };
    std::function<double(double,int)> kx_f = [](double age, int type) { return 1.; };

    std::shared_ptr<CPlantBox::MappedSegments> rs;

    std::vector<double> kr, kr_t;
    std::vector<double> kx, kx_t;
    std::vector<std::vector<double> > krs, krs_t;
    std::vector<std::vector<double>> kxs, kxs_t;

protected:

    double kr_const(double age, int type) { return kr[0]; }
    double kr_perType(double age, int type) { return kr.at(type); }
    double kr_table(double age, int type) { return interp1(age, kr_t, kr); }
    double kr_tablePerType(double age, int type) { return interp1(age, krs_t.at(type), krs.at(type)); }

    double kx_const(double age, int type) { return kx[0]; }
    double kx_perType(double age, int type) { return kx.at(type); }
    double kx_table(double age, int type) { return interp1(age, kx_t, kx); }
    double kx_tablePerType(double age, int type) { return interp1(age, kxs_t.at(type), kxs.at(type)); }

    static double interp1(double ip, std::vector<double> x, std::vector<double> y);

};

} // namespace

#endif
