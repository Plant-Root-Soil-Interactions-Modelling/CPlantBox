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
 */
void MappedSegments::setRectangularGrid(Vector3d min, Vector3d max, Vector3d res, bool cut)
{
    minBound = min;
    maxBound = max;
    resolution = res;
    rectangularGrid = cut;
    cutSegments(); // re-add (for cutting)
    sort(); // todo only necessary
    seg2cell.clear(); // re-map all segments
    cell2seg.clear();
    mapSegments(segments);
}


/**
 * Update the mappers root2cell, which maps root segment index to soil cell index, and
 * cell2seg which maps soil cell index to multiple root segments.
 *
 * @param segs      the (new) segments that need to be mapped
 */
void MappedSegments::mapSegments(const std::vector<Vector2i>& segs) {
    for (auto& ns : segs) {
        Vector3d mid = (nodes[ns.x].plus(nodes[ns.y])).times(0.5);
        int cellIdx = soil_index(mid.x,mid.y,mid.z);
        if (cellIdx>=0) {
            int segIdx = ns.y-1; // this is unique in a tree like structured
            seg2cell[segIdx] = cellIdx;
            if (cell2seg.count(cellIdx)>0) {
                cell2seg[cellIdx].push_back(segIdx);
            } else {
                cell2seg[cellIdx] = std::vector<int>({segIdx});
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
void MappedSegments::cutSegments() {
    assert(segments.size()==radii.size() && "MappedSegments::addSegments: number of segments and radii disagree!");
    assert(segments.size()==types.size() && "MappedSegments::addSegments: number of segments and types disagree!");
    if (rectangularGrid) {
        int n = segments.size(); // segs.size() will change within the loop (recursive implementation)
        for (int i=0; i<n; i++ ) {
            addSegment(segments[i], radii[i], types[i], i);
        }
    }
}

/**
 * Adds a single segment at index i, appends the remaining segments, if segment is cut
 *
 * This approach may run into problems if a segment is located exactly along a face.
 */
void MappedSegments::addSegment(Vector2i ns, double r,  int t, int ii) {
    Vector3d n1 = nodes[ns.x];
    Vector3d n2 = nodes[ns.y];
    Vector3d mid = (n1.plus(n2)).times(0.5);
    int im = soil_index(mid.x,mid.y,mid.z); // cell indices
    int in1 = soil_index(n1.x,n1.y,n1.z);
    int in2 = soil_index(n2.x,n2.y,n2.z);
    if ((im!=in1) or (im!=in2)) { // cut
        // build SDF
        auto width = maxBound.minus(minBound); // construct sdf
        Vector3d dx(width.x/resolution.x, width.y/resolution.y, width.z/resolution.z);
        auto mid0 = mid.minus(minBound);
        int x = std::floor(mid0.x/dx.x);
        int y = std::floor(mid0.y/dx.y);
        int z = std::floor(mid0.z/dx.z);
        SDF_Cuboid sdf; // create a signed distance function for cutting
        Vector3d minB(x*dx.x, y*dx.y, z*dx.z);
        minB = minB.plus(minBound);
        Vector3d maxB((x+1)*dx.x, (y+1)*dx.y, (z+1)*dx.z);
        maxB = maxB.plus(minBound);
        sdf.min = minB;
        sdf.max = maxB;
        // std::cout << minB.toString() << ", " << maxB.toString() << ", width " << width.toString() << ", " << sdf.min.toString() << ", " << sdf.max.toString() << "\n";

        im = sdf.getDist(mid)>0; // redo indices, since accuracy of pickking may differ
        in1 = sdf.getDist(n1)>0;
        in2 = sdf.getDist(n2)>0;
        if ((im!=in1) or (im!=in2)) {
            Vector3d cPoint;
            if (im==in1) { // is one node at mid (sort accordingly)
                // std::cout << "n1 " << sdf.getDist(n1) << " mid " << sdf.getDist(mid) << " n2 " << sdf.getDist(n2) << ",indices "<< in1 << ", " << im << ", " << in2 << "\n";
                if (sdf.getDist(n2)<0) {
                    cPoint = SegmentAnalyser::cut(n2, mid, std::make_shared<SDF_Cuboid>(sdf), eps);
                } else {
                    cPoint = SegmentAnalyser::cut(mid, n2, std::make_shared<SDF_Cuboid>(sdf), eps);
                }
                nodeCTs.push_back(nodeCTs[ns.y]); // todo: we might linearly interpolate
            } else if (im==in2) {
                // std::cout << "n2 " << sdf.getDist(n2) << " mid " << sdf.getDist(mid) << " n1 " << sdf.getDist(n1) << ", " << in2 << ", " << im << ", " << in1 << "\n";
                if (sdf.getDist(n1)<0) {
                    cPoint = SegmentAnalyser::cut(n1, mid, std::make_shared<SDF_Cuboid>(sdf), eps);
                } else {
                    cPoint = SegmentAnalyser::cut(mid, n1, std::make_shared<SDF_Cuboid>(sdf), eps);
                }
                nodeCTs.push_back(nodeCTs[ns.x]); // todo: we might linearly interpolate
            } else { // otherwise split in mid, use cutSegments on those
                cPoint = mid;
                nodeCTs.push_back(0.5*(nodeCTs[ns.x]+nodeCTs[ns.y]));
            }
            // std::cout << "[" << n1.toString() << n2.toString() << "] -> [" << nodes[ns.x].toString() << ", " << nodes.back().toString() << "], ["<< nodes.back().toString() << ", " << n2.toString() << "], " << "\n";
            nodes.push_back(cPoint);
            Vector2i s1(ns.x, nodes.size()-1);
            Vector2i s2(nodes.size()-1, ns.y);
            if ((length(s1)<eps) ||  (length(s2)<eps)) { // if the cut segments are too small, just give up
                add(ns, r, t, ii);
                nodes.pop_back(); // remove cPoint
            } else {
                addSegment(s1, r, t , ii); // first segment replaces at index ii
                addSegment(s2, r, t , -1); // append second segment
            }
        } else { // im==in1==in2, dont't cut
            // std::cout << "ok " << ii <<": (" << ns.x <<", " << ns.y << ") [" << n1.toString() <<", "<< n2.toString() <<"\n";
            add(ns, r, t, ii);
        }
    } else { // im==in1==in2, dont't cut
        // std::cout << "ok " << ii <<": (" << ns.x <<", " << ns.y << ") [" << n1.toString() <<", "<< n2.toString() <<"\n";
        add(ns, r, t, ii);
    }
}

/**
 * Adds the segment at index i, or appends it, if i = -1
 */
void MappedSegments::add(Vector2i s, double r,  int t, int i) {
    if (i>=0) {
        segments[i] = s;
        radii[i] = r;
        types[i] = t;
    } else {
        segments.push_back(s);
        radii.push_back(r);
        types.push_back(t);
    }
}

/**
 * Length of the segment
 */
double MappedSegments::length(const Vector2i& s) const {
    return (nodes.at(s.y).minus(nodes.at(s.x))).length();
}

/**
 * Removes segments from the mappers
 */
void MappedSegments::removeSegments(const std::vector<Vector2i>& segs) {
    for (auto& ns : segs) {
        int cellIdx = -1;
        int segIdx = ns.y-1;
        if (seg2cell.count(segIdx)>0) { // remove from seg2cell
            cellIdx = seg2cell[segIdx];
            auto it = seg2cell.find(segIdx);
            seg2cell.erase(it);
        } else {
            throw std::invalid_argument("MappedSegments::removeSegments: warning segment index "+ std::to_string(segIdx)+ " was not found in the seg2cell mapper");
        }
        if (cell2seg.count(cellIdx)>0) {
            auto& csegs= cell2seg[cellIdx];
            int c = 0;
            for (int i=0; i<csegs.size(); i++) {
                if (csegs[i] == segIdx) {
                    csegs.erase(csegs.begin() + c, csegs.begin() + c);
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
 * Maps a point into a cell and return the cells linear index (for a equidistant rectangular domain)
 */
int MappedSegments::soil_index_(double x, double y, double z) { // np.array([1, 3, 5])
    Vector3d p(x,y,z);
    std::array<double,3>  r = { resolution.x, resolution.y, resolution.z};
    auto w = maxBound.minus(minBound);
    auto p0 = p.minus(minBound);
    std::array<double,3> i = { p0.x/w.x*r[0], p0.y/w.y*r[1], p0.z/w.z*r[2] };
    for (int k=0; k<3; k++) {
        if ((i[k] < 0) or (i[k] >= r[k])) {
            return -1; // point is out of domain
        }
    }
    return std::floor(i[2]) * r[0] * r[1] + std::floor(i[1]) * r[0] + std::floor(i[0]); // a linear index not periodic
}

/**
 * Sorts the segments, so that the segment index == second node index -1 (unique mapping in a tree)
 */
void MappedSegments::sort() {
    auto newSegs = segments;
    auto newRadii = radii;
    auto newTypes = types;

    for (int i=0; i<newSegs.size(); i++) {
        int ind = segments[i].y-1;
        newSegs[ind] = segments[i];
        newRadii[ind] = radii[i];
        newTypes[ind] = types[i];
    }
    segments = newSegs;
    radii = newRadii;
    types = newTypes;
}


/**
 * Overridden, to map initial shoot segments (@see RootSystem::initialize).
 *
 * Shoot segments have per default radii = 0.1 cm, types = 0
 * This can be changed by directly accessing the member variables.
 */
void MappedRootSystem::initializeLB(int basaltype, int shootbornetype, bool verbose) {
    std::cout << "MappedRootSystem::initialize \n" << std::flush;
    RootSystem::initializeLB(basaltype, shootbornetype, verbose);
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
    this->initializeLB(4,5, verbose);
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
    if (verbose) {
        std::cout << "nodes moved "<< uni.size() << "\n" << std::flush;
    }
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
    if (verbose) {
        std::cout << "new nodes added " << newnodes.size() << "\n" << std::flush;
    }
    auto newsegs = this->getNewSegments(); // add segments (TODO cutting)
    segments.resize(segments.size()+newsegs.size());
    for (auto& ns : newsegs) {
        segments[ns.y-1] = ns;
    }
    if (verbose) {
        std::cout << "segments added "<< newsegs.size() << "\n" << std::flush;
    }
    auto newsegO = this->getNewSegmentOrigins(); // to add radius and type (TODO cutting)
    radii.resize(radii.size()+newsegO.size());
    types.resize(types.size()+newsegO.size());
    c = 0;
    if (verbose) {
        std::cout << "Number of segments " << radii.size() << ", including " << newsegO.size() << " new \n"<< std::flush;
    }
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
