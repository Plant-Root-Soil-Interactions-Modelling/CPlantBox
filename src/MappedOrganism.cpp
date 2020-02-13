// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "MappedOrganism.h"

namespace CPlantBox {


MappedSegments::MappedSegments(std::vector<Vector3d> nodes, std::vector<double> nodeCTs, std::vector<Vector2i> segs,
    std::vector<double> radii, std::vector<int> types) : nodes(nodes), nodeCTs(nodeCTs), segments(segs), radii(radii), types(types)
{ }

MappedSegments::MappedSegments(std::vector<Vector3d> nodes, std::vector<Vector2i> segs, std::vector<double> radii) : nodes(nodes), segments(segs) {
    nodeCTs.resize(nodes.size());
    for (int i=0; i<segments.size(); i++) {
        nodeCTs[i] = 0;
    }
    if (radii.size()==0) {
        setRadius(0.1); // cm
    }
    setTypes(0);
}

/**
 * Sets the soil cell index call back
 */
void MappedSegments::setSoilGrid(const std::function<int(double,double,double)>& s) {
    soil_index = s;
    seg2cell.clear();
    cell2seg.clear();
    mapSegments(segments);
}

void MappedSegments::setRadius(double a) {
    radii.resize(segments.size());
    for (int i=0; i<segments.size(); i++) {
        radii[i] = a;
    }
}

void MappedSegments::setTypes(int t) {
    types.resize(segments.size());
    for (int i=0; i<segments.size(); i++) {
        types[i] = t;
    }
}

/**
 * update the mappers with the new segments
 */
void MappedSegments::mapSegments(std::vector<Vector2i> segs) {
    for (auto& ns : segs) {
        Vector3d mid = (nodes[ns.x].plus(nodes[ns.y])).times(0.5);
        int segIdx = ns.y-1;
        int cellIdx = soil_index(mid.x,mid.y,mid.z);
        seg2cell[segIdx] = cellIdx;
        if (cell2seg.count(cellIdx)>0) {
            cell2seg[cellIdx].push_back(segIdx);
        } else {
            cell2seg[cellIdx] = std::vector<int>({segIdx});
        }
    }
}



/**
 * Overridden, to map initial shoot segments,
 * shoot segments have per defautl radii = 0.1 cm, types =0
 * (can be changed by directly accessing the member variables)
 */
void MappedRootSystem::initialize(int basaltype, int shootbornetype, bool verbose) {
    std::cout << "MappedRootSystem::initialize \n" << std::flush;
    RootSystem::initialize(basaltype, shootbornetype, verbose);
    segments = this->getShootSegments();
    nodes = this->getNodes();
    nodeCTs = this->getNodeCTs();
    radii = { 0.1 };
    types = { 0 };
    mapSegments(segments);
}

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
    // update moved nodes (TODO)

}

/**
 * Assembles the linear system as sparse matrix, given by
 * indices aI, aJ, and values aV, and load aB
 *
 * The linear system is solved to yield the homogeneous solution
 */
void XylemFlux::linearSystem(double simTime) {

    int Ns = rs->segments.size(); // number of segments
    aI.resize(4*Ns);
    aJ.resize(4*Ns);
    aV.resize(4*Ns);
    int N = rs->nodes.size(); // number of nodes
    aB.resize(N);

    size_t k=0;
    for (int si = 0; si<Ns; si++) {

        int i = rs->segments[si].x;
        int j = rs->segments[si].y;
        auto n1 = rs->nodes[i];
        auto n2 = rs->nodes[j];

        double age = simTime - rs->nodeCTs[j];
        int type = rs->types[j - 1];
        double a = rs->radii[j - 1];

        double kx = kx_f(age, type);
        double kr = kr_f(age, type);

        auto v = n2.minus(n1);
        double l = v.length();
        double vz = v.z / l; // normed direction

        double c = 2.*a * M_PI * kr / kx; // Eqn (2)
        double d = std::exp(-std::sqrt(c) * l) - std::exp(std::sqrt(c) * l); // Eqn (5)
        double di = 1. / d;

        double cii = -kx * di * std::sqrt(c) * (std::exp(-std::sqrt(c) * l) + std::exp(std::sqrt(c) * l)); // Eqn 16
        double cij = 2 * kx * di * std::sqrt(c);  // Eqn 17
        double bi = kx * vz; //  # Eqn 18 (* rho * g)

        aB[k] += bi;
        aI[k] = i; aJ[k]= i; aV[k] = cii;
        k += 1;
        aI[k] = i; aJ[k] = j;  aV[k] = cij;
        k += 1;

        int ii = i;
        i = j;  j = ii; // edge ji
        aB[i] += bi;
        aI[k] = i; aJ[k]= i; aV[k] = cii;
        k += 1;
        aI[k] = i; aJ[k] = j;  aV[k] = cij;
        k += 1;

    }
}

/**
 * Creates the inhomogeneous solution from the homogeneous one
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


}
