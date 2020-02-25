// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "MappedOrganism.h"

#include <algorithm>
#include <functional>


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
{ }

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
}

/**
 * Sets the soil cell index call back function. Resets and updates the mappers.
 * The callback function takes a spatial coordinate [cm] and returns the index of the cell
 */
void MappedSegments::setSoilGrid(const std::function<int(double,double,double)>& s, Vector3d min, Vector3d max, Vector3d res) {
    minBound = min;
    maxBound = max;
    resolution = res;
    rectangularGrid = true;
    this->setSoilGrid(s);
}

/**
 * Sets the soil cell index call back function. Resets and updates the mappers.
 * The callback function takes a spatial coordinate [cm] and returns the index of the cell
 */
void MappedSegments::setSoilGrid(const std::function<int(double,double,double)>& s) {
    soil_index = s;
    seg2cell.clear();
    cell2seg.clear();
    mapSegments(segments); // TODO use rectangularGrid
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
 * Update the mappers root2cell, which maps root segment index to soil cell index, and
 * cell2seg which maps soil cell index to multiple root segments.
 *
 * @param segs      the (new) segments that need to be mapped
 */
void MappedSegments::mapSegments(std::vector<Vector2i> segs) {
    if (rectangularGrid) {
        segs = cutSegments(segs);
    }
    for (auto& ns : segs) {
        Vector3d mid = (nodes[ns.x].plus(nodes[ns.y])).times(0.5);
        int cellIdx = soil_index(mid.x,mid.y,mid.z);
        if (cellIdx>0) {
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
 * Removes segments from the mappers todo make private
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
 * use soil_index on mid and end points, cut at rectangular grid
 * TODO
 */
std::vector<Vector2i> MappedSegments::cutSegments(std::vector<Vector2i> segs) const {
    return segs;
}
//Vector3d newnode = cut(nodes[s.x], nodes[s.y], geometry);
//nodes.push_back(newnode); // add new segment
//Vector2i newseg(s.x,nodes.size()-1);
//seg.push_back(newseg);
//sO.push_back(segO.at(i));
//std::vector<Vector2i> MappedSegments::cutSegments(std::vector<Vector2i> segs) {
//    std::vector<Vector2i> newsegs;
//    for (auto& ns : segs) {
//       auto n1 = nodes[ns.x];
//       auto n2 = nodes[ns.y];
//       int c1 = soil_index(n1.x,n1.y,n1.z);
//       int c2 = soil_index(n1.x,n1.y,n1.z);
//       if (c1==c2) { // nothing todo
//           newsegs.push_back(ns);
//       } else { // split
//
//       }
//    }
//    return newsegs;
//}



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

/**
 * TODO docme
 */
void MappedRootSystem::setRectangularGrid(Vector3d min, Vector3d max, Vector3d res)
{
   MappedSegments::minBound = min;
   MappedSegments::maxBound = max;
   MappedSegments::resolution = res;
   MappedSegments::rectangularGrid = true;
}



/**
 * Assembles the linear system as sparse matrix, given by
 * indices aI, aJ, and values aV, and load aB
 *
 * The linear system is solved to yield the homogeneous solution
 *
 * @simTime         [days] current simulation time is needed for age dependend conductivities,
 *                  to calculate the age from the creation times.
 */
void XylemFlux::linearSystem(double simTime) {

    int Ns = rs->segments.size(); // number of segments
    aI.resize(4*Ns);
    aJ.resize(4*Ns);
    aV.resize(4*Ns);
    int N = rs->nodes.size(); // number of nodes
    aB.resize(N);
    std::fill(aB.begin(), aB.end(), 0.);
    std::fill(aV.begin(), aV.end(), 0.);
    std::fill(aI.begin(), aI.end(), 0);
    std::fill(aJ.begin(), aJ.end(), 0);
    size_t k=0;
    for (int si = 0; si<Ns; si++) {

        int i = rs->segments[si].x;
        int j = rs->segments[si].y;
        auto n1 = rs->nodes[i];
        auto n2 = rs->nodes[j];

        double a = rs->radii[si]; // si is correct, with ordered and unordered segmetns
        double age = simTime - rs->nodeCTs[j];
        int type = rs->types[si];
        double kx = kx_f(age, type);
        double  kr = kr_f(age, type);

        auto v = n2.minus(n1);
        double l = v.length();
        double vz = v.z / l; // normed direction

        double c = 2.*a * M_PI * kr / kx; // Eqn (2)
        double d = std::exp(-std::sqrt(c) * l) - std::exp(std::sqrt(c) * l); // Eqn (5)
        double di = 1. / d;

        double cii = -kx * di * std::sqrt(c) * (std::exp(-std::sqrt(c) * l) + std::exp(std::sqrt(c) * l)); // Eqn 16
        double cij = 2 * kx * di * std::sqrt(c);  // Eqn 17
        double bi = kx * vz; //  # Eqn 18 (* rho * g) todo change back to kr [1/day]
        // std::cout << "cii " << cii << ", c " << c << ", d "<< d << "\n";

        aB[i] += bi;
        aI[k] = i; aJ[k]= i; aV[k] = cii;
        k += 1;
        aI[k] = i; aJ[k] = j;  aV[k] = cij;
        k += 1;

        int ii = i;
        i = j;  j = ii; // edge ji
        aB[i] -= bi; // Eqn 14 with changed sign
        aI[k] = i; aJ[k]= i; aV[k] = cii;
        k += 1;
        aI[k] = i; aJ[k] = j;  aV[k] = cij;
        k += 1;
    }
}

/**
 * Creates the inhomogeneous solution from the homogeneous one
 *
 * @param rx        root xylem solution per node
 * @param sx        soil solution per cell
 */
std::vector<double> XylemFlux::getSolution(std::vector<double> rx, std::vector<double> sx) {
    int cIdx = rs->seg2cell[0]; // the first node ends no segment attached to it
    rx[0] += sx[cIdx];
    for (int i = 1; i< rx.size(); i++) {
        int cIdx = rs->seg2cell[i-1]; // i-1 is the segment index
        rx[i] += sx[cIdx];
    }
    return rx;
}

/**
 *  Sets the radial conductivity in [1 day-1], converts to [cm2 day g-1] by dividing by rho*g
 */
void XylemFlux::setKr(std::vector<double> values, std::vector<double> age) {
    std::transform(values.begin(), values.end(), values.begin(), std::bind1st(std::multiplies<double>(),1./(rho * g)));
    kr = values;
    kr_t = age;
    if (age.size()==0) {
        if (values.size()==1) {
            kr_f = std::bind(&XylemFlux::kr_const, this, std::placeholders::_1, std::placeholders::_2);
            std::cout << "Kr is constant " << values[0] << " cm2 day g-1 \n";
        } else {
            kr_f  = std::bind(&XylemFlux::kr_perType, this, std::placeholders::_1, std::placeholders::_2);
            std::cout << "Kr is constant per type, type 0 = " << values[0] << " cm2 day g-1 \n";
        }
    } else {
        kr_f  = std::bind(&XylemFlux::kr_table, this, std::placeholders::_1, std::placeholders::_2);
        std::cout << "Kr is age dependent\n";
    }
}


/**
 *  Sets the axial conductivity in [cm3 day-1], converts to [cm5 day g-1] by dividing by rho*g
 */
void XylemFlux::setKx(std::vector<double> values, std::vector<double> age) {
    std::transform(values.begin(), values.end(), values.begin(), std::bind1st(std::multiplies<double>(),1./(rho * g)));
    kx = values;
    kx_t = age;
    if (age.size()==0) {
        if (values.size()==1) {
            kx_f = std::bind(&XylemFlux::kx_const, this, std::placeholders::_1, std::placeholders::_2);
            std::cout << "Kx is constant " << values[0] << " cm2 day g-1 \n";
        } else {
            kx_f  = std::bind(&XylemFlux::kx_perType, this, std::placeholders::_1, std::placeholders::_2);
            std::cout << "Kx is constant per type, type 0 = " << values[0] << " cm2 day g-1 \n";
        }
    } else {
        kx_f  = std::bind(&XylemFlux::kx_table, this, std::placeholders::_1, std::placeholders::_2);
        std::cout << "Kx is age dependent\n";
    }
}

/**
 *  Sets the radial conductivity in [1 day-1], converts to [cm2 day g-1] by dividing by rho*g
 */
void XylemFlux::setKrTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age) {
    krs.resize(0);
    for (auto v :values) {
        std::transform(v.begin(), v.end(), v.begin(), std::bind1st(std::multiplies<double>(),1./(rho * g)));
        krs.push_back(v);
    }
    krs_t = age;
    kr_f = std::bind(&XylemFlux::kr_tablePerType, this, std::placeholders::_1, std::placeholders::_2);
    std::cout << "Kr is age dependent per root type\n";
}

/**
 *  Sets the axial conductivity in [cm3 day-1], converts to [cm5 day g-1] by dividing by rho*g
 */
void XylemFlux::setKxTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age) {
    kxs.resize(0);
    for (auto v :values) {
        std::transform(v.begin(), v.end(), v.begin(), std::bind1st(std::multiplies<double>(),1./(rho * g)));
        kxs.push_back(v);
    }
    kxs_t = age;
    kx_f = std::bind(&XylemFlux::kx_tablePerType, this, std::placeholders::_1, std::placeholders::_2);
    std::cout << "Kx is age dependent per root type\n";
}

/**
 *
 */
double XylemFlux::interp1(double ip, std::vector<double> x, std::vector<double> y) {
    if (ip > x.back()) return y.back(); // check bounds
    if (ip < x[0]) return y[0];

    // if we are within bounds find the index of the lower bound
    const auto lookUpIndex = std::distance(x.begin(), std::lower_bound(x.begin(), x.end(), ip));
    if (lookUpIndex == 0) {
        return y[0];
    }
    double ip_ = (ip - x[lookUpIndex-1])/(x[lookUpIndex] - x[lookUpIndex-1]);
    return y[lookUpIndex-1]*(1.0 - ip_)  + y[lookUpIndex]*ip_;
}

/**
 * Fluxes from root segments into a the soil cell with cell index cIdx
 *
 * Approximation(!) TODO exact
 */
std::map<int,double> XylemFlux::soilFluxes(double simTime, std::vector<double> rx_hom)
{
    std::map<int,double> fluxes;

    for (int si = 0; si<rs->segments.size(); si++) {

        int i = rs->segments[si].x;
        int j = rs->segments[si].y;
        int segIdx = j-1;

        if (rs->seg2cell.count(segIdx)>0) {

            int cellIdx = rs->seg2cell[segIdx];

            double a = rs->radii[si]; // si is correct, with ordered and unordered segments
            double age = simTime - rs->nodeCTs[j];
            int type = rs->types[si];
            double  kr = kr_f(age, type);

            auto n1 = rs->nodes[i];
            auto n2 = rs->nodes[j];
            double l = (n2.minus(n1)).length();

            double f = - 2*a*M_PI*l*(kr*rho*g)*(-rx_hom[j]); // cm3 / day

            if (fluxes.count(cellIdx)==0) {
                fluxes[cellIdx] = f;
            } else {
                fluxes[cellIdx] += f; // sum up fluxes per cell
            }
        } else {
            std::cout << "XylemFlux::soilFluxes: Warning! unmapped segments with index " << segIdx << "\n";
        }

    }
    return fluxes;
}


}
