// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "MappedOrganism.h"

#include "SegmentAnalyser.h"

#include <algorithm>
#include <functional>
#include <cmath>

namespace CPlantBox {

/**
 * A static root system, as needed for flux computations, represented as
 *
 * @param nodes     [cm]
 * @param nodeCTs   node creation times [d]
 * @param segs      describes a segment with two node indices segs.x, and segs.y [1]
 * @param radii     [cm] segment radius
 * @param types     root type or order of the segment [1]
 */
MappedSegments::MappedSegments(std::vector<Vector3d> nodes, std::vector<double> nodeCTs, std::vector<Vector2i> segs,
    std::vector<double> radii, std::vector<int> types) : nodes(nodes), nodeCTs(nodeCTs), segments(segs), radii(radii), types(types)
{
    assert((nodes.size()==nodeCTs.size()) && "MappedSegments::MappedSegments: Unequal vector sizes nodes and nodeCTs");
    assert((segments.size()==radii.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and radii");
    assert((segments.size()==types.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and types");
}

/**
 *  A static root system, as needed for flux computations.
 *
 * @param nodes     [cm]
 * @param segs      describes a segment with two node indices segs.x, and segs.y [1]
 * @param radii     [cm] segment radius
 *
 * nodeCTs is set to 0. for all segments
 * types are set to type 0 for all segments
 */

MappedSegments::MappedSegments(std::vector<Vector3d> nodes, std::vector<Vector2i> segs, std::vector<double> radii)
:nodes(nodes), segments(segs), radii(radii) {
    nodeCTs.resize(nodes.size());
    std::fill(nodeCTs.begin(), nodeCTs.end(), 0.);
    setTypes(0);
    assert((segments.size()==radii.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and radii");
    assert((segments.size()==types.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and types");
}

/**
 * Sets the radius for all segments.
 */
void MappedSegments::setRadius(double a) {
    radii.resize(segments.size());
    std::fill(radii.begin(), radii.end(), a);
}

/**
 * Sets the type for all segments.
 */
void MappedSegments::setTypes(int t) {
    types.resize(segments.size());
    std::fill(types.begin(), types.end(), t);
}

/**
 * Sets the soil cell index call back function and a rectangular grids.
 * First cuts all segments at the grid boundaries, then resets and updates the mappers.
 *
 * The callback function takes a spatial coordinate [cm] and returns the index of the cell
 */
void MappedSegments::setSoilGrid(const std::function<int(double,double,double)>& s, Vector3d min, Vector3d max, Vector3d res) {
    soil_index = s;
    this->setRectangularGrid(min,max,res);
    this->setSoilGrid(s);
}

/**
 * Sets the soil cell index call back function. Resets and updates the mappers.
 *
 * The callback function takes a spatial coordinate [cm] and returns the index of the cell
 */
void MappedSegments::setSoilGrid(const std::function<int(double,double,double)>& s) {
    soil_index = s;
    seg2cell.clear(); // re-map all segments
    cell2seg.clear();
    mapSegments(segments);
}

/**
 * Sets a rectangular grid, and cuts all segments along the grid cells
 *
 * TODO grids are not mapped correctly after cutting
 */
void MappedSegments::setRectangularGrid(Vector3d min, Vector3d max, Vector3d res)
{
    minBound = min;
    maxBound = max;
    resolution = res;
    rectangularGrid = true;
    auto segs = segments; // copy
    auto types_ = types;
    auto radii_ = radii;
    segments.clear(); // clear
    types.clear();
    radii.clear();
    addSegments(segs, radii_, types_); // re-add
}


/**
 * Update the mappers root2cell, which maps root segment index to soil cell index, and
 * cell2seg which maps soil cell index to multiple root segments.
 *
 * @param segs      the (new) segments that need to be mapped
 */
void MappedSegments::mapSegments(std::vector<Vector2i> segs) {
    for (auto& ns : segs) {
        Vector3d mid = (nodes[ns.x].plus(nodes[ns.y])).times(0.5);
        int cellIdx = soil_index(mid.x,mid.y,mid.z);
        if (cellIdx>=0) {
            int segIdx = ns.y-1; // this is unique in a tree like structured
            seg2cell[segIdx] = cellIdx;
            if (cell2seg.count(cellIdx)>0) {
                cell2seg[cellIdx].push_back(ns.x);
                cell2seg[cellIdx].push_back(ns.y);
            } else {
                cell2seg[cellIdx] = std::vector<int>({ns.x, ns.y});
            }
        } else {
            std::cout << "MappedSegments::mapSegments: warning segment with mid " << mid.toString() << " exceeds domain, skipped segment \n";
        }
    }
}

/**
 * Adds the segments to the list.
 * Optionally, cut segments @param segs at a rectangular grid (@see MappedSegments::setSoilGrid)
 */
void MappedSegments::addSegments(const std::vector<Vector2i>& segs, const std::vector<double>& radii_,  const std::vector<int>& types_) {
    assert(segs.size()==radii_.size() && "MappedSegments::addSegments: number of segments and radii disagree!");
    assert(segs.size()==types_.size() && "MappedSegments::addSegments: number of segments and types disagree!");
    if (!rectangularGrid) {
        for (int i=0; i<segs.size(); i++ ) {
            segments.push_back(segs[i]);
            radii.push_back(radii_[i]);
            types.push_back(types_[i]);
        }
    } else {
        SDF_Cuboid sdf;
        for (int i=0; i<segs.size(); i++ ) {
            Vector2i ns = segs[i];
            Vector3d n1 = nodes[ns.x]; // since we cut at the bounds, they are often located at the bounds
            Vector3d n2 = nodes[ns.y];
            Vector3d mid = (n1.plus(n2)).times(0.5);
            Vector3d v = (n2.minus(n1)).times(eps); // never evaluate soil_index exactly at the bounds (no unique return value)
            n1 = n1.plus(v);
            n2 = n2.minus(v);
            int im = soil_index(mid.x,mid.y,mid.z); // cell indices
            int in1 = soil_index(n1.x,n1.y,n1.z);
            int in2 = soil_index(n2.x,n2.y,n2.z);

            if ((im!=in1) or (im!=in2)) { // cut
                auto width = maxBound.minus(minBound); // construct sdf
                Vector3d dx(width.x/resolution.x, width.y/resolution.y, width.z/resolution.z);
                auto mid0 = mid.minus(minBound);
                int x = std::floor(mid0.x/dx.x); // important in case of minus
                int y = std::floor(mid0.y/dx.y);
                int z = std::floor(mid0.z/dx.z);
                Vector3d minB(x*dx.x, y*dx.y, z*dx.z);
                minB = minB.plus(minBound);
                Vector3d maxB((x+1)*dx.x, (y+1)*dx.y, (z+1)*dx.z);
                maxB = maxB.plus(minBound);
                sdf.min = minB;
                sdf.max = maxB;
    //            std::cout << "Segment "<< ns.y-1 << "/" << segs.size() << ", index " << x << ", " << ", " << y << ", " << z <<
    //                "\nCuboid: "<< sdf.toString() << " at " << mid.toString() << ": dist = " << sdf.getDist(mid) << "\n";
                if (im==in1) { // is one node at mid (sort accordingly)
    //                std::cout << "im==in1, \nn2: " << n2.toString() << ", " <<  sdf.getDist(n2) << ", " << im << ", " <<in2 << "\n"  << std::flush;
    //                std::cout << "n1: " <<n1.toString() << ", " <<  sdf.getDist(n1) << ", " << im << ", " << in1 << "\n"  << std::flush;
                    auto cPoint = SegmentAnalyser::cut(mid, n2, std::make_shared<SDF_Cuboid>(sdf), eps);
                    nodes.push_back(cPoint);
                    nodeCTs.push_back(nodeCTs[ns.y]); // nn, todo: we might linearly interpolate
                } else if (im==in2) {
    //                std::cout << "im==in2, \nn1: " <<n1.toString() << ", " <<  sdf.getDist(n1) << ", " << im << ", " <<in1 << "\n"  << std::flush;
    //                std::cout << "n2: " << n2.toString() << ", " <<  sdf.getDist(n2) << ", " << im << ", " <<in2 << "\n"  << std::flush;
                    auto cPoint = SegmentAnalyser::cut(mid, n1, std::make_shared<SDF_Cuboid>(sdf), eps);
                    nodes.push_back(cPoint);
                    nodeCTs.push_back(nodeCTs[ns.x]); // nn, todo: we might linearly interpolate
                } else { // otherwise split in mid, use cutSegments on those
                    std::cout << "mid\n" << std::flush;
                    nodes.push_back(mid);
                    nodeCTs.push_back(0.5*(nodeCTs[ns.x]+nodeCTs[ns.x]));
                }
                Vector2i s1(ns.x, nodes.size()-1);
                Vector2i s2(nodes.size()-1, ns.y);
                int t = types_[i];
                double r = radii_[i];
                addSegments({s1,s2}, {r,r}, {t,t});
            } else { // dont't cut
                segments.push_back(ns);
                radii.push_back(radii_[i]);
                types.push_back(types_[i]);
            }
        }
    }
}

/**
 * Removes segments from the mappers
 */
void MappedSegments::removeSegments(std::vector<Vector2i> segs) {
    for (auto& ns : segs) {
        int cellIdx = -1;
        int segIdx = ns.y-1;
        if (seg2cell.count(segIdx)>0) {
            cellIdx = seg2cell[segIdx];
            auto it = seg2cell.find(segIdx);
            seg2cell.erase(it);
        } else {
            throw std::invalid_argument("MappedSegments::removeSegments: warning segment index "+ std::to_string(segIdx)+ " was not found in the seg2cell mapper");
        }
        if (cell2seg.count(cellIdx)>0) {
            auto& segs= cell2seg[cellIdx];
            int c = 0;
            for (int i=1; i<segs.size(); i+=2) {
                int ni = segs[i];
                if (ni == ns.y) {
                    segs.erase(segs.begin() + c -1, segs.begin() + c); // cannot be the first
                    break; // inner for
                }
                c++;
            }
        } else {
            throw std::invalid_argument("MappedSegments::removeSegments: warning cell index "+ std::to_string(cellIdx)+ " was not found in the cell2seg mapper");
        }
    }
}

/**
 * linear index
 */
int MappedSegments::soil_index_(double x, double y, double z) {
    Vector3d p(x,y,z);
    std::array<double,3>  r = { resolution.x, resolution.y, resolution.z};
    auto w = maxBound.minus(minBound);
    auto p0 = p.minus(minBound);
    std::array<double,3> i = { p0.x/w.x*r[0],p0.y/w.y*r[1],p0.z/w.z*r[2] };
    for (int k=0; k<3; k++) {
        if ((i[k] < 0) or (i[k] >= r[k])) {
            return -1;
        }
    }
    return std::floor(i[0]) * r[1] * r[2] + std::floor(i[1]) * r[1] + std::floor(i[2]); // a linear index not periodic
}



/**
 * Overridden, to map initial shoot segments (@see RootSystem::initialize).
 *
 * Shoot segments have per default radii = 0.1 cm, types = 0
 * This can be changed by directly accessing the member variables.
 */
void MappedRootSystem::initialize(int basaltype, int shootbornetype, bool verbose) {
    std::cout << "MappedRootSystem::initialize \n" << std::flush;
    RootSystem::initialize(basaltype, shootbornetype, verbose);
    segments = this->getShootSegments();
    nodes = this->getNodes();
    nodeCTs = this->getNodeCTs();
    radii.resize(segments.size());
    std::fill(radii.begin(), radii.end(), 0.1);
    types.resize(segments.size());
    std::fill(types.begin(), types.end(), 0);
    mapSegments(segments);
}

/**
 * Overridden, to map initial shoot segments (@see RootSystem::initialize).
 *
 * Shoot segments have per default radii = 0.1 cm, types = 0
 * This can be changed by directly accessing the member variables.
 */
void MappedRootSystem::initialize(bool verbose) {
    this->initialize(4,5, verbose);
}

/**
 * Simulates the development of the organism in a time span of @param dt days.
 *
 * @param dt        time step [day]
 * @param verbose   turns console output on or off
 */
void MappedRootSystem::simulate(double dt, bool verbose)
{
    if (soil_index==nullptr) {
        throw std::invalid_argument("MappedRootSystem::simulate():soil was not set, use MappedRootSystem::simulate::setSoilGrid" );
    }

    RootSystem::simulate(dt,verbose);

    auto uni = this->getUpdatedNodeIndices(); // move nodes
    auto unodes = this->getUpdatedNodes();
    assert(uni.size()==unodes.size() && "updated node indices and number of nodes must be equal");
    int c = 0;
    for (int i : uni) {
        nodes.at(i) = unodes[c];
        c++;
    }
    std::cout << "nodes moved \n" << std::flush;
    auto newnodes = this->getNewNodes(); // add nodes
    nodes.reserve(nodes.size()+newnodes.size());
    for (auto& nn : newnodes) {
        nodes.push_back(nn);
    }
    auto newnode_cts = this->getNewNodeCTs(); // add node cts
    nodeCTs.reserve(nodeCTs.size()+newnode_cts.size());
    for (auto& nct : newnode_cts) {
        nodeCTs.push_back(nct);
    }
    std::cout << "new nodes added \n" << std::flush;
    auto newsegs = this->getNewSegments(); // add segments
    segments.resize(segments.size()+newsegs.size());
    for (auto& ns : newsegs) {
        segments[ns.y-1] = ns;
    }
    std::cout << "segments added \n" << std::flush;
    auto newsegO = this->getNewSegmentOrigins(); // add radius and type
    std::cout << "resize(): << " << radii.size()+newsegO.size() << "\n"<< std::flush;
    radii.resize(radii.size()+newsegO.size());
    types.resize(types.size()+newsegO.size());
    c = 0;
    std::cout << "total length " << radii.size() << ", new " << newsegO.size() << "\n "<< std::flush;
    for (auto& so : newsegO) {
        int segIdx = newsegs[c].y-1;
        c++;
        radii[segIdx] = so->getParam()->a;
        types[segIdx] = so->getParam()->subType;
    }
    // map new segments
    this->mapSegments(newsegs);

    // update segments of moved nodes
    std::vector<Vector2i> rSegs;
    for (int i : uni) {
        int segIdx = i -1;
        int cellIdx = seg2cell[segIdx];
        auto s = segments[segIdx];
        Vector3d mid = (nodes[s.x].plus(nodes[s.y])).times(0.5);
        int newCellIdx = soil_index(mid.x,mid.y,mid.z);
        // 1. check if mid is still in same cell (otherwise, remove, and add again)
        // 2. if cut is on, check if end point is in same cell than mid point (otherwise remove and add again)
        bool remove = false;
        if (cellIdx==newCellIdx) {
            if (rectangularGrid) {
                auto endPoint = nodes[s.y];
                newCellIdx = soil_index(endPoint.x,endPoint.y,endPoint.z);
                remove = (newCellIdx!=cellIdx);
            }
        } else {
            remove = true;
        }
        if (remove) {
            rSegs.push_back(s);
        }
    }
    MappedSegments::removeSegments(rSegs);
    MappedSegments::mapSegments(rSegs);
}

} // namespace
