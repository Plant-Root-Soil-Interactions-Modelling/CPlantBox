// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef MAPPED_ROOTSYSTEM_H_
#define MAPPED_ROOTSYSTEM_H_

#include "RootSystem.h"
#include "Plant.h"

#include <functional>

namespace CPlantBox {

/**
 * Represents connected 1d rootsystem segments, which are mapped to a 3d soil grid.
 *
 * Holds nodes, nodeCTs, segments, radii, and types
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
    void mapSegments(std::vector<Vector2i> segs);

    std::map<int, int> seg2cell; // root segment to soil cell mappper
    std::map<int, std::vector<int>> cell2seg; // soil cell to root segment mapper

    std::function<int(double,double,double)> soil_index = [](double x, double y, double z) { return 0; }; ///< soil cell index call back function

    std::vector<Vector3d> nodes; ///< nodes [nodes]
    std::vector<double> nodeCTs; ///< creation times [days]
    std::vector<Vector2i> segments; ///< connectivity of the nodes
    std::vector<double> radii; ///< radii [cm]
    std::vector<int> types; ///< types [1]

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

    void simulate(double dt, bool verbose = false) override; ///< build nodes and segmentssequentially

};


/**
 * Hybrid solver (Meunier et al. )
 */
class XylemFlux
{
public:

    XylemFlux(std::shared_ptr<CPlantBox::MappedSegments> rs) :rs(rs) { }

    virtual ~XylemFlux() { }

    void setKrF(const std::function<int(double,int)>&  kr) { kr_f = kr; } ///< sets a callback for kr:=kr(age,type)
    void setKxF(const std::function<int(double,int)>&  kx) { kx_f = kx; } ///< sets a callback for kx:=kx(age,type)

    void linearSystem(double simTime = 0.); ///< builds linear system (simTime is needed for age dependent conductivities)
    std::vector<double> getSolution(std::vector<double> rx, std::vector<double> sx); ///< creates the solution from the homogeneous solution

    std::vector<int> aI; // to assemble the sparse matrix on the Python side
    std::vector<int> aJ;
    std::vector<double> aV;
    std::vector<double> aB;

    std::function<int(double,int)> kr_f = [](double age, int type) { return 0.; };
    std::function<int(double,int)> kx_f = [](double age, int type) { return 1.; };

    double rho = 1; // [g cm-3]
    double g =  9.8065*100.*24.*3600.*24.*3600.;  // [cm day-2]

    std::shared_ptr<CPlantBox::MappedSegments> rs;
};

}

#endif

