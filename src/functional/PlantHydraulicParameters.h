// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef PLANT_HYDRAULIC_PARAMETERS
#define PLANT_HYDRAULIC_PARAMETERS

#include "MappedOrganism.h"

namespace CPlantBox {

/**
 * Plant hydraulic conductivities, i.e. axial (kx) and radial (kr) conductivities
 *
 * TODO needs testing and documentation
 * TODO what would be a good method to store and load this class
 */


class PlantHydraulicParameters
{

public:

    PlantHydraulicParameters() { }
    PlantHydraulicParameters(std::shared_ptr<CPlantBox::MappedSegments> rs) : ms(rs) { }

    virtual ~PlantHydraulicParameters() { }

    void setKr(std::vector<double> values, std::vector<double> age = std::vector<double>(0), bool verbose = true); ///< sets a callback for kr:=kr(age,type),  [1 day-1]
    void setKx(std::vector<double> values, std::vector<double> age = std::vector<double>(0), bool verbose = true); ///< sets a callback for kx:=kx(age,type),  [cm3 day-1]
    void setKr(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, double kr_length_ = -1.0, bool verbose = false); ///< sets a callback for kr:=kr(age,type),  [1 day-1]
    void setKx(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, bool verbose = false); ///< sets a callback for kx:=kx(age,type),  [cm3 day-1]
    void setKrTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, bool verbose = false, bool ageBased = true);
    void setKxTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, bool verbose = false);
    void setKrTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age, bool verbose, bool ageBased = true);
    void setKxTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age, bool verbose);
    void setKrValues(std::vector<double> values, bool verbose = false); ///< one value per segment
    void setKxValues(std::vector<double> values, bool verbose = false); ///< one value per segment

    std::function<double(int, double, int, int)> kr_f = [](int si, double age, int type, int orgtype){
        throw std::runtime_error("kr_f not implemented"); return 0.; };
    std::function<double(int, double,int,int)> kx_f = [](int si, double age, int type, int orgtype) {
        throw std::runtime_error("kx_f not implemented"); return 1.; };
    double kr_f_wrapped(int si, double age, int type, int orgtype, bool cells) const;///stops transpiration if organs are not in the correct domain

    virtual double getPsiOut(bool cells, int si, const std::vector<double>& sx_, bool verbose) const; ///< get the outer water potential [cm]

    std::vector<std::vector<double>> kr, kr_t; //  [1 day-1]
    std::vector<std::vector<double>> kx, kx_t; // [cm3 day-1]
    std::vector<std::vector<std::vector<double>>> krs, krs_t;
    std::vector<std::vector<std::vector<double>>> kxs, kxs_t;

    std::vector<double> getEffKr(double simtime);
    std::vector<double> getKr(double simtime);
    std::vector<double> getKx(double simtime);

    std::shared_ptr<CPlantBox::MappedSegments> ms; // TODO needs to be set for kr_RootExchangeZonePerType, and kr_tablePerType_distance

    double psi_air = -954378; // air water potential [cm] for T = 20°C and RH = 0.5

protected:

    double kr_const(int si,double age, int type, int organType)  {
        return kr.at(0).at(0);
    } // type correspond to subtype or to the leaf segment number

    double kr_perOrgType(int si,double age, int type, int organType) {
        return kr.at(organType - 2).at(0);
    } // per organ type (goes from 2 (root) to 4 (leaf))

    double kr_perType(int si,double age, int type, int organType) {
        return kr.at(organType - 2).at(type);
    } // per subtype and organ type (goes from 2 (root) to 4 (leaf))

    double kr_table(int si,double age, int type, int organType) {
        double kr_ = Function::interp1(age, kr_t.at(0), kr.at(0));
        return kr_;
    } // constant for all type/subtype and age dependant

    double kr_tablePerOrgType(int si,double age, int type, int organType) {
        double kr_ = Function::interp1(age, krs_t.at(organType-2).at(0), krs.at(organType-2).at(0));
        return kr_;
    } // constant for all subtype but type and age dependant

    double kr_tablePerType(int si,double age, int type, int organType) {
        double kr_ = Function::interp1(age, krs_t.at(organType-2).at(type), krs.at(organType-2).at(type));
        return kr_;
    } // subtype, type and age dependant

    double kr_valuePerSegment(int si, double age, int type, int organType) {
        return kr.at(0).at(si);
    }

    double kr_RootExchangeZonePerType(int si,double age, int type, int organType) { //when use carbon- and water-limited growth, canNOT use "kr_tablePerType" instead of this function
        if (organType == Organism::ot_root){
            double coef = ms->exchangeZoneCoefs.at(si);//% of segment length in the root exchange zone, see MappedPlant::simulate
            return coef * kr.at(organType - 2).at(type);
        }
        return kr.at(organType - 2).at(type);
    } // subtype, type and depend on distance to tip for roots

    double kr_tablePerType_distance(int si,double age, int type, int organType) { //when use carbon- and water-limited growth, canNOT use "kr_tablePerType" instead of this function
        if (organType == Organism::ot_root){
            //double coef = rs->exchangeZoneCoefs.at(si);//% of segment length in the root exchange zone, see MappedPlant::simulate
            double distFromTip = ms->distanceTip.at(si);//% of segment length in the root exchange zone, see MappedPlant::simulate
            double kr_ = Function::interp1(distFromTip, krs_t.at(organType-2).at(type), krs.at(organType-2).at(type));
            return kr_;//coef * kr.at(organType - 2).at(type);
        }
        return krs.at(organType - 2).at(type).at(0);
    } // subtype, type and depend on distance to tip for roots

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
