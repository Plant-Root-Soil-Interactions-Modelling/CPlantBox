// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef MAPPED_ORGANISM_H_
#define MAPPED_ORGANISM_H_

#include "RootSystem.h"
#include "Plant.h"

namespace CPlantBox {

/**
 * Interface to pick a soil cell
 */
class Pickable
{
public:
    virtual ~Pickable() { };
    virtual int pick(double x, double y, double z) { return 0; }; //< the cell index of a point
};

/**
 * Manages the 1d rootsystem grid, which is mapped to a 3d soil grid.
 * TODO templated solution and MappedPlant (did not work with pybind11)
 *
 * Overwrites simulate to:
 * build a 1d sequentially, consisting of nodes, nodeCTs, segments, radii, and types
 * map segments to cell indices and vice versa
 */
// template <class BaseOrganism>
class MappedRootSystem : public RootSystem
{
public:

    using RootSystem::RootSystem;

    void setSoilGrid(std::shared_ptr<Pickable> s) { soil = s; }

    void simulate(double dt, bool verbose = false) override;

    std::vector<Vector3d> nodes; ///< nodes [nodes]
    std::vector<double> nodeCTs; ///< creation times [days]
    std::vector<Vector2i> segments; ///< connectivity of the nodes
    std::vector<double> radii; ///< radii [cm]
    std::vector<double> types; ///< radii [cm]

    // index mappers
    std::map<int, int> seg2cell;
    std::map<int, std::vector<int>> cell2seg;

protected:

    std::shared_ptr<Pickable> soil;

};

//using MappedRootSystem = MappedOrganism<RootSystem>;
//using MappedPlant = MappedOrganism<Plant>;

}

#endif

