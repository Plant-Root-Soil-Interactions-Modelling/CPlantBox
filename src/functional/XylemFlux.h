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
        const std::vector<double> soil_k = std::vector<double>(), bool verbose = false); ///< builds linear system (simTime is needed for age dependent conductivities)

    std::map<int,double> soilFluxes(double simTime, const std::vector<double>& rx, const std::vector<double>& sx,
    		bool approx = false, const std::vector<double> soil_k = std::vector<double>()); // [cm3/day]
    std::vector<double> segFluxes(double simTime, const std::vector<double>& rx, const std::vector<double>& sx,
    		bool approx = false, bool cells = false, const std::vector<double> soil_k = std::vector<double>(),
			bool verbose = false) const; // for each segment in [cm3/day]
    std::map<int,double> sumSegFluxes(const std::vector<double>& segFluxes); ///< sums segment fluxes over soil cells,  soilFluxes = sumSegFluxes(segFluxes), [cm3/day]

    std::vector<double> splitSoilFluxes(const std::vector<double>& soilFluxes, int type = 0) const; ///< splits soil fluxes (per cell) into segment fluxes

    std::vector<int> aI; // to assemble the sparse matrix on the Python side
    std::vector<int> aJ;
    std::vector<double> aV;
    std::vector<double> aB;

    void setKr(std::vector<double> values, std::vector<double> age = std::vector<double>(0), bool verbose = false); ///< sets a callback for kr:=kr(age,type),  [1 day-1]
    void setKx(std::vector<double> values, std::vector<double> age = std::vector<double>(0), bool verbose = false); ///< sets a callback for kx:=kx(age,type),  [cm3 day-1]
    void setKrTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, bool verbose = false, bool ageBased = true);
    void setKxTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, bool verbose = false);
    void setKr(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, double kr_length_ = -1.0, bool verbose = false); ///< sets a callback for kr:=kr(age,type),  [1 day-1]
    void setKx(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, bool verbose = false); ///< sets a callback for kx:=kx(age,type),  [cm3 day-1]
    void setKrTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age, bool verbose, bool ageBased = true);
    void setKxTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age, bool verbose);
    void setKrValues(std::vector<double> values, bool verbose = false); ///< one value per segment
    void setKxValues(std::vector<double> values, bool verbose = false); ///< one value per segment


   std::function<double(int, double, int, int)> kr_f = [](int si, double age, int type, int orgtype){
		throw std::runtime_error("kr_f not implemented"); return 0.; };
    std::function<double(int, double,int,int)> kx_f = [](int si, double age, int type, int orgtype) {
		throw std::runtime_error("kx_f not implemented"); return 1.; };
	double kr_f_wrapped(int si, double age, int type, int orgtype, bool cells) const;///stops transpiration if organs are not in the correct domain

	virtual size_t fillVectors(size_t k, int i, int j, double bi, double cii, double cij, double psi_s) ; ///< fill the vectors aI, aJ, aV, aB
	virtual double getPsiOut(bool cells, int si, const std::vector<double>& sx_, bool verbose) const; ///< get the outer water potential [cm]
    std::vector<double> getEffKr(double simtime);
    std::vector<double> getKr(double simtime);
    std::vector<double> getKx(double simtime);
    std::vector<double> getHs(const std::vector<double>& sx);

    std::shared_ptr<CPlantBox::MappedSegments> rs;

    std::vector<std::vector<double>> kr, kr_t; //  [1 day-1]
    std::vector<std::vector<double>> kx, kx_t; // [cm3 day-1]
    std::vector<std::vector<std::vector<double>>> krs, krs_t;
    std::vector<std::vector<std::vector<double>>> kxs, kxs_t;

    double psi_air = -954378; // air water potential [cm] for T = 20°C and RH = 0.5

protected:

	//type correspond to subtype or to the leaf segment number
    double kr_const(int si,double age, int type, int organType) //k constant
	{
		return kr.at(0).at(0);
	}

	double kr_perOrgType(int si,double age, int type, int organType)
	{
		return kr.at(organType - 2).at(0);
	} //per organ type (goes from 2 (root) to 4 (leaf))
    double kr_perType(int si,double age, int type, int organType)
	{
		return kr.at(organType - 2).at(type);
	}//per subtype and organ type (goes from 2 (root) to 4 (leaf))
    double kr_table(int si,double age, int type, int organType)
	{
		double kr_ = Function::interp1(age, kr_t.at(0), kr.at(0));
		return kr_;
	} //constant for all type/subtype and age dependant


	double kr_tablePerOrgType(int si,double age, int type, int organType)
	{
		double kr_ = Function::interp1(age, krs_t.at(organType-2).at(0), krs.at(organType-2).at(0));
		return kr_;
	}//constant for all subtype but type and age dependant

		double kr_tablePerType(int si,double age, int type, int organType) {
		double kr_ = Function::interp1(age, krs_t.at(organType-2).at(type), krs.at(organType-2).at(type));
	    return kr_;
	} //subtype, type and age dependant
	double kr_valuePerSegment(int si, double age, int type, int organType)
	{
		return kr.at(0).at(si);
	}
	double kr_RootExchangeZonePerType(int si,double age, int type, int organType)//when use carbon- and water-limited growth, canNOT use "kr_tablePerType" instead of this function
	{
		if (organType == Organism::ot_root){
			double coef = rs->exchangeZoneCoefs.at(si);//% of segment length in the root exchange zone, see MappedPlant::simulate
			return coef * kr.at(organType - 2).at(type);
		}
		return kr.at(organType - 2).at(type);
	} //subtype, type and depend on distance to tip for roots

	double kr_tablePerType_distance(int si,double age, int type, int organType)//when use carbon- and water-limited growth, canNOT use "kr_tablePerType" instead of this function
	{
    if (organType == Organism::ot_root){
			//double coef = rs->exchangeZoneCoefs.at(si);//% of segment length in the root exchange zone, see MappedPlant::simulate
			double distFromTip = rs->distanceTip.at(si);//% of segment length in the root exchange zone, see MappedPlant::simulate 
			
			double kr_ = Function::interp1(distFromTip, krs_t.at(organType-2).at(type), krs.at(organType-2).at(type));
			return kr_;//coef * kr.at(organType - 2).at(type);
		}
		return krs.at(organType - 2).at(type).at(0);
	} //subtype, type and depend on distance to tip for roots

    double kx_const(int si,double age, int type, int organType) { return kx.at(0).at(0); } //k constant
    double kx_perOrgType(int si,double age, int type, int organType) { return kx.at(organType - 2)[0]; } //per organ type (goes from 2 (root) to 4 (leaf))
    double kx_perType(int si,double age, int type, int organType) { return kx.at(organType - 2).at(type); } //per subtype and organ type (goes from 2 (root) to 4 (leaf))
	double kx_table(int si,double age, int type, int organType) { return Function::interp1(age, kx_t[0], kx[0]); } //constant for all type/subtype and age dependant
    double kx_tablePerOrgType(int si,double age, int type, int organType) { return Function::interp1(age, kxs_t.at(organType-2).at(0), kxs.at(organType-2).at(0)); } //constant for all subtype but type and age dependant
    double kx_tablePerType(int si,double age, int type, int organType) { return Function::interp1(age, kxs_t.at(organType-2).at(type), kxs.at(organType-2).at(type)); } //subtype, type and age dependant
    double kx_valuePerSegment(int si, double age, int type, int organType) { return kx.at(0).at(si); };


};

} // namespace

#endif
