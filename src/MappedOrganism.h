// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef MAPPED_ROOTSYSTEM_H_
#define MAPPED_ROOTSYSTEM_H_

#include "RootSystem.h"
#include "Plant.h"

#include <functional>

namespace CPlantBox {

/**
 * Connected 1d rootsystem segments, which are mapped to a 3d soil grid.
 */
class MappedSegments
{
public:

    MappedSegments() { }

    MappedSegments(std::vector<Vector3d> nodes, std::vector<double> nodeCTs, std::vector<Vector2i> segs,
        std::vector<double> radii, std::vector<int> types);
    ///< for kr and kx age and type dependent

    MappedSegments(std::vector<Vector3d> nodes, std::vector<Vector2i> segs, std::vector<double> radii = std::vector<double>(0));
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
 * Manages the 1d rootsystem grid, which is mapped to a 3d soil grid.
 *
 * Overwrites simulate to:
 * 1) build a 1d grid sequentially, consisting of nodes, nodeCTs, segments, radii, and types
 * 2) map segments to cell indices and vice versa
 */
// template <class BaseOrganism>
class MappedRootSystem : public MappedSegments, public RootSystem
{
public:

    using RootSystem::RootSystem;

    void initialize(bool verbose = true) override;
    void initialize(int basaltype, int shootbornetype, bool verbose = true) override;

    void simulate(double dt, bool verbose = false) override;

};


/**
 * Hybrid solver (Meunier et al. )
 */
class XylemFlux
{
public:

    XylemFlux(std::shared_ptr<CPlantBox::MappedSegments> rs)  :rs(rs) { }

    virtual ~XylemFlux() { }

    void setKrF(const std::function<int(double,int)>&  kr) { kr_f = kr; }
    void setKxF(const std::function<int(double,int)>&  kx) { kx_f = kx; }

    void linearSystem(double simTime = 0.);
    std::vector<double> getSolution(std::vector<double> rx, std::vector<double> sx);

    /**
     * Fluxes from root segments into a the soil cell with cell index cIdx [TODO]
     * (needed by the soil part)
     */
//    virtual double roots2cell(rx, sx_old)
//    {
//        if (rs->cell2seg.count(cIdx)>0) {
//            auto sIdxs = rs->cell2seg.at(cIdx);
//            double flux = 0;
//            for (int i : sIdxs) {
//                double f = 0.;
//                if (i < oldRootX.size()) {
//                    double rootP = oldRootX[i];
//                    double a = rs->radii[i];
//                    auto n1 = rs->nodes[rs->segments[i].x];
//                    auto n2 = rs->nodes[rs->segments[i].y];
//                    double l = (n2.minus(n1)).length();
//
//                    // f = 2*
//
//                }
//                flux += f;
//            }
//
//            return flux;
//        } else {
//            return 0.;
//        }
//    }

    // to assemble the sparse matrix on the pyhthon side
    std::vector<int> aI;
    std::vector<int> aJ;
    std::vector<double> aV;
    std::vector<double> aB;

    std::function<int(double,int)> kr_f = [](double age, int type) { return 0.; };
    std::function<int(double,int)> kx_f = [](double age, int type) { return 1.; };

    double rho = 1; // cm^3
    double g =  9.8065*100.*24.*3600.*24.*3600.;  // [cm/day^2]

    std::shared_ptr<CPlantBox::MappedSegments> rs;
};





}

#endif

