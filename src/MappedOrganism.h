// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef MAPPED_ROOTSYSTEM_H_
#define MAPPED_ROOTSYSTEM_H_

#include "RootSystem.h"
#include "Plant.h"

#include <functional>
#include <vector>

namespace CPlantBox {

/**
 * Represents a connected 1d rootsystem as segments, which are mapped to a 3d soil grid.
 *
 * Holds nodes, nodeCTs, segments, radii, and types, which can be directly accessed.
 */
class MappedSegments
{
public:

    MappedSegments() { }

    MappedSegments(std::vector<Vector3d> nodes, std::vector<double> nodeCTs, std::vector<Vector2i> segs,
        std::vector<double> radii, std::vector<int> types);
    ///< for kr and kx age and type dependent

    MappedSegments(std::vector<Vector3d> nodes, std::vector<Vector2i> segs, std::vector<double> radii);
    ///< for constant kr, and kx

    void setRadius(double a); ///< sets a constant radius for all segments
    void setTypes(int t); ///< sets a constant type for all segments

    void setSoilGrid(const std::function<int(double,double,double)>& s); ///< sets the soil, resets the mappers and maps the segments
    void setSoilGrid(const std::function<int(double,double,double)>& s, Vector3d min, Vector3d max, Vector3d res); ///< sets the soil, resets the mappers and maps the segments

    void mapSegments(std::vector<Vector2i> segs);

    void removeSegments(std::vector<Vector2i> segs);
    std::vector<Vector2i> cutSegments(std::vector<Vector2i> segs) const;

    std::map<int, int> seg2cell; // root segment to soil cell mappper
    std::map<int, std::vector<int>> cell2seg; // soil cell to root segment mapper represented as two node indices n1, n2

    std::function<int(double,double,double)> soil_index = [](double x, double y, double z) { return 0; }; ///< soil cell index call back function

    std::vector<Vector3d> nodes; ///< nodes [nodes]
    std::vector<double> nodeCTs; ///< creation times [days]
    std::vector<Vector2i> segments; ///< connectivity of the nodes
    std::vector<double> radii; ///< radii [cm]
    std::vector<int> types; ///< types [1]

    Vector3d minBound;
    Vector3d maxBound;
    Vector3d resolution;
    bool rectangularGrid = false;

};


/**
 * Build MappedSegmentds sequentially from a RootSystem
 */
class MappedRootSystem : public MappedSegments, public RootSystem
{
public:

    using RootSystem::RootSystem;

    void initialize(bool verbose = true) override; ///< overridden, to map initial shoot segments,
    void initialize(int basaltype, int shootbornetype, bool verbose = true) override; ///< overridden, to map initial shoot segments,

    void simulate(double dt, bool verbose = false) override; ///< build nodes and segments sequentially

    void setRectangularGrid(Vector3d min, Vector3d max, Vector3d res); ///< sets an underlying rectangular grid, for cutting segments

};


/**
 * Hybrid solver (Meunier et al. )
 *
 * Units are [cm], [g], and [day], they are fixed by choosing g, and rho
 */
class XylemFlux
{
public:

    XylemFlux(std::shared_ptr<CPlantBox::MappedSegments> rs) :rs(rs) { }

    virtual ~XylemFlux() { }

    void linearSystem(double simTime = 0.); ///< builds linear system (simTime is needed for age dependent conductivities)
    std::vector<double> getSolution(std::vector<double> rx, std::vector<double> sx); ///< creates the solution from the homogeneous solution

    std::map<int,double> soilFluxes(double simTime, std::vector<double> rxz);
    std::map<int,double> soilFluxesApprox(double simTime, std::vector<double> rxz);

    std::vector<int> aI; // to assemble the sparse matrix on the Python side
    std::vector<int> aJ;
    std::vector<double> aV;
    std::vector<double> aB;

    void setKr(std::vector<double> values, std::vector<double> age); ///< sets a callback for kr:=kr(age,type)
    void setKx(std::vector<double> values, std::vector<double> age); ///< sets a callback for kx:=kx(age,type)
    void setKrTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age);
    void setKxTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age);

    std::function<double(double,int)> kr_f = [](double age, int type) { return 0.; };
    std::function<double(double,int)> kx_f = [](double age, int type) { return 1.; };

    double rho = 1; // [g cm-3]
    double g =  9.8065*100.*24.*3600.*24.*3600.;  // [cm day-2]

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

}

#endif

