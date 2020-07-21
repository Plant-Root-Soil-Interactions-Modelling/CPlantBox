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
 *
 * Optionally, cuts segment along the boundaries of a rectangular grid
 */
class MappedSegments
{
public:

    MappedSegments() { }

    MappedSegments(std::vector<Vector3d> nodes, std::vector<double> nodeCTs, std::vector<Vector2i> segs,
        std::vector<double> radii, std::vector<int> types); ///< for kr and kx age and type dependent
    MappedSegments(std::vector<Vector3d> nodes, std::vector<Vector2i> segs, std::vector<double> radii); ///< for constant kr, and kx

    void setRadius(double a); ///< sets a constant radius for all segments
    void setTypes(int t); ///< sets a constant type for all segments

    void setSoilGrid(const std::function<int(double,double,double)>& s); ///< sets the soil, resets the mappers and maps all segments
    void setSoilGrid(const std::function<int(double,double,double)>& s, Vector3d min, Vector3d max, Vector3d res); ///< sets the soil, resets the mappers and maps all segments
    void setRectangularGrid(Vector3d min, Vector3d max, Vector3d res); ///< sets an underlying rectangular grid, and cuts all segments accordingly

    void mapSegments(std::vector<Vector2i> segs);
    void addSegments(const std::vector<Vector2i>& segs, const std::vector<double>& radii,  const std::vector<int>& types); // cut and add segments

    std::map<int, int> seg2cell; // root segment to soil cell mapper
    std::map<int, std::vector<int>> cell2seg; // soil cell to root segment mapper

    std::function<int(double,double,double)> soil_index =
        std::bind(&MappedSegments::soil_index_, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3); ///< soil cell index call back function, (care need all MPI ranks in case of dumux)

    void sort(); ///< sorts segments, each segment belongs to position s.y-1

    std::vector<Vector3d> nodes; ///< nodes [cm]
    std::vector<double> nodeCTs; ///< creation times [days]
    std::vector<Vector2i> segments; ///< connectivity of the nodes
    std::vector<double> radii; ///< radii [cm]
    std::vector<int> types; ///< types [1]

    Vector3d minBound;
    Vector3d maxBound;
    Vector3d resolution; // cells
    bool rectangularGrid = false;

    const double eps = 1.e-6;

protected:

    int soil_index_(double x, double y, double z); // default mapper to a equidistant rectangular grid
    void removeSegments(std::vector<Vector2i> segs); ///< remove segments from the mappers

};



/**
 * Build MappedSegmentds sequentially from a RootSystem
 */
class MappedRootSystem : public MappedSegments, public RootSystem
{
public:

    using RootSystem::RootSystem;

    void initialize(bool verbose = true); ///< overridden, to map initial shoot segments,
    void initializeLB(int basaltype, int shootbornetype, bool verbose = true); ///< overridden, to map initial shoot segments,

    void simulate(double dt, bool verbose = false) override; ///< build nodes and segments sequentially

    /* segments are shoot and root segments */


    std::shared_ptr<MappedSegments> mappedSegments() { return std::make_shared<MappedSegments>(*this); }  // up-cast for Python binding
    std::shared_ptr<RootSystem> rootSystem() { return std::make_shared<RootSystem>(*this); }; // up-cast for Python binding

};

}

#endif

