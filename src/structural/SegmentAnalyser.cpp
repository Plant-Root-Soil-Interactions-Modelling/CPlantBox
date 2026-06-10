// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "SegmentAnalyser.h"

#include "Organ.h"
#include "Organism.h"
#include "MappedOrganism.h"
#include "PlantHydraulicParameters.h"
#include "PlantHydraulicModel.h"

#include <algorithm>
#include <iomanip>
#include <istream>
#include <iostream>
#include <string>

#include <fstream>
#include <set>
#include <math.h>

namespace CPlantBox {

/**
 * Constructs from raw vectors, initialising data fields "creationTime" and "radius".
 * Other parameters are retrieved on demand via getParameter(); missing values fall back to @p def or throw.
 *
 * @param nodes     node coordinates
 * @param segments  segment connectivity; each entry holds two node indices such that the segment runs from nodes[s.x] to nodes[s.y]
 * @param segCTs    creation time of each segment [d]
 * @param radii     radius of each segment [cm]
 */
SegmentAnalyser::SegmentAnalyser(const std::vector<Vector3d>& nodes, const std::vector<Vector2i>& segments,
    const std::vector<double>& segCTs, const std::vector<double>& radii) :nodes(nodes), segments(segments)
{
    assert((segments.size() == segCTs.size()) && "SegmentAnalyser::SegmentAnalyser(): Unequal vector sizes");
    assert((segments.size() == radii.size()) && "SegmentAnalyser::SegmentAnalyser(): Unequal vector sizes");
    data["creationTime"] = segCTs;
    data["radius"] = radii;
    segO = std::vector<std::weak_ptr<Organ>>(segments.size()); // create expired
}

/**
 * Constructs from an Organism, copying all segments and extracting per-segment
 * parameters (radius, subType, id, organType) from the organ objects.
 *
 * @param plant  organism to analyse
 */
SegmentAnalyser::SegmentAnalyser(const Organism& plant) 
{
    //std::cout << "construct from Organism\n";
    nodes = plant.getNodes();
    segments = plant.getSegments();
    auto segCTs = plant.getSegmentCTs();
    assert(segments.size()==segCTs.size() && "SegmentAnalyser::SegmentAnalyser(Organism p): Unequal vector sizes");
    data["creationTime"] = segCTs;
    auto sego = plant.getSegmentOrigins();
    segO = std::vector<std::weak_ptr<Organ>>(segments.size());
    auto radii = std::vector<double>(segments.size());
    auto subType = std::vector<double>(segments.size());
    auto id = std::vector<double>(segments.size());
    auto organType = std::vector<double>(segments.size());
    for (size_t i=0; i<segments.size(); i++) {
        segO[i] = sego[i]; // convert shared_ptr to weak_ptr
        radii[i] = segO[i].lock()->getParameter("radius");
        subType[i] = segO[i].lock()->getParameter("subType");
        id[i] = segO[i].lock()->getParameter("id");
        organType[i] = segO[i].lock()->getParameter("organType"); // = 2
    }
    data["radius"] = radii;
    data["subType"] = subType;
    data["id"] = id;
    data["organType"] = organType;
}

/**
 * Constructs from a MappedSegments object, copying nodes, segments, and per-segment
 * parameters (radius, subType, organType). Creation times are derived from nodeCTs.
 *
 * @param plant  mapped segment structure to analyse
 */
SegmentAnalyser::SegmentAnalyser(const MappedSegments& plant) :nodes(plant.nodes), segments(plant.segments)
{
    // std::cout << "construct from MappedSegments\n";
    assert((segments.size()==plant.radii.size()) && "SegmentAnalyser::SegmentAnalyser(MappedSegments p): Unequal vector sizes");
    assert((segments.size()==plant.subTypes.size()) && "SegmentAnalyser::SegmentAnalyser(MappedSegments p): Unequal vector sizes");
    assert((segments.size()==plant.organTypes.size()) && "SegmentAnalyser::SegmentAnalyser(MappedSegments p): Unequal vector sizes");
    std::vector<double> segCTs(plant.segments.size());
    std::vector<double> subTypesd(plant.subTypes.size()); // convert to double
    std::vector<double> organTypesd(plant.organTypes.size()); // convert to double
    for (size_t i=0; i<segments.size(); i++) {
        segCTs[i] = plant.nodeCTs.at(segments[i].y);
        subTypesd[i] = double(plant.subTypes[i]);
        organTypesd[i] = double(plant.organTypes[i]);
    }
    data["creationTime"] = segCTs;
    data["radius"] = plant.radii;
    data["subType"] = subTypesd;
    data["organType"] = organTypesd;
    segO = plant.segO;
}

/**
 * Appends all segments from @p plant.
 */
void SegmentAnalyser::addSegments(const Organism& plant)
{
    addSegments(SegmentAnalyser(plant));
}

/**
 * Appends all segments from analyser @p a.
 * Only data keys present in both objects are merged; keys missing in @p a are dropped.
 */
void SegmentAnalyser::addSegments(const SegmentAnalyser& a)
{
    int offset = nodes.size();
    nodes.insert(nodes.end(),a.nodes.begin(),a.nodes.end()); // copy nodes
    auto ns = a.segments;
    for (auto& s : ns) { // shift indices
        s.x += offset;
        s.y += offset;
    }
    segments.insert(segments.end(),ns.begin(),ns.end()); // copy segments
    segO.insert(segO.end(),a.segO.begin(),a.segO.end());// copy origins
    assert(segments.size() == segO.size() && "SegmentAnalyser::addSegments(): Unequal vector sizes, segments and segment origins" );
    for (auto iter = data.begin(); iter != data.end(); ++iter) {
        std::string key =  iter->first;
        if (a.data.count(key)>=0) {
            data[key].insert(data[key].end(),a.data.at(key).begin(),a.data.at(key).end());
            assert(segments.size() == data[key].size() && "SegmentAnalyser::addSegments(): Unequal vector sizes, segments and data" );
        } else {
            data.erase(key);
        }
    }
}

/**
 * Adds a single segment (e.g. an artificial shoot).
 * All other data fields are padded with -1 for the new segment.
 *
 * @param seg     segment as two node indices into the analyser's node list (note: pack() may reorder nodes)
 * @param ct      creation time [d]
 * @param radius  radius [cm]
 * @param insert  if true, prepends the segment; otherwise appends it (default: false)
 */
void SegmentAnalyser::addSegment(Vector2i seg, double ct, double radius, bool insert)
{
    if (insert) {
        segments.insert(segments.begin(),seg);
        data["creationTime"].insert(data["creationTime"].begin(),ct);
        segO.insert(segO.begin(),std::weak_ptr<Organ>()); // expired
        data["radius"].insert(data["radius"].begin(),radius);
        for (auto iter = data.begin(); iter != data.end(); ++iter) { // insert -1 for all other data
            std::string key =  iter->first;
            if (data[key].size() < segments.size()) {
                data[key].insert(data[key].begin(), -1.);
            }
        }
    } else {
        segments.push_back(seg);
        data["creationTime"].push_back(ct);
        segO.push_back(std::weak_ptr<Organ>()); // expired
        data["radius"].push_back(radius);
        for (auto iter = data.begin(); iter != data.end(); ++iter) { // push -1 for all other data
            std::string key =  iter->first;
            if (data[key].size() < segments.size()) {
                data[key].push_back(-1.);
            }
        }
    }
    for (auto iter = data.begin(); iter != data.end(); ++iter) { // check again if all data sizes are correct
        std::string key =  iter->first;
        assert(segments.size() == data[key].size() && "SegmentAnalyser::addSegment(): Unequal vector sizes, segments and data" );
    }
}

/**
 * Computes segment age as (@p simTime - creationTime), clamps to zero, and stores it as "age".
 */
void SegmentAnalyser::addAge(double simTime)
{
    std::vector<double> age(segments.size());
    for (size_t i=0; i<age.size(); i++) {
        double a = simTime - data["creationTime"].at(i);
        age.at(i) = std::max(a,0.);
    }
    this->addData("age",age);
}

/**
 * Evaluates radial (kr) and axial (kx) conductivities per segment using @p rs
 * and stores them as "kr" and "kx". Values are clamped to @p kr_max and @p kx_max.
 *
 * @param rs      hydraulic parameter object providing kr_f() and kx_f()
 * @param simTime current simulation time [d]
 * @param kr_max  upper bound for kr (default 1e6)
 * @param kx_max  upper bound for kx (default 1e6)
 */
void SegmentAnalyser::addHydraulicConductivities(const PlantHydraulicParameters& rs, double simTime, double kr_max, double kx_max)
{
    std::vector<double> kr(segments.size());
    std::vector<double> kx(segments.size());
    for (size_t i=0; i<kr.size(); i++) {
        double age = simTime - data["creationTime"].at(i);
        int subType = (int) data["subType"].at(i);
        int organType = (int) data["organType"].at(i);
        kr.at(i) = std::min(rs.kr_f(i, age, subType, organType), kr_max);
        kx.at(i) = std::min(rs.kx_f(i, age, subType, organType), kx_max);
    }
    this->addData("kr",kr);
    this->addData("kx",kx);
}


/**
 * Computes radial flux [cm³/cm²/day] and axial flux [cm³/day] per segment
 * and stores them as "radial_flux" and "axial_flux".
 * Requires addHydraulicConductivities() to be called first.
 *
 * @param rs      hydraulic model providing getRadialFluxes()
 * @param rx      xylem pressure head per node [cm]
 * @param sx      soil pressure head per segment [cm]
 * @param simTime current simulation time [d]
 */
void SegmentAnalyser::addFluxes(const PlantHydraulicModel& rs, const std::vector<double>& rx, const std::vector<double>& sx, double simTime) {

    std::vector<double> radial_flux = rs.getRadialFluxes(simTime, rx, sx, true, false); // volumetric flux, approx = true, cells = false
    std::vector<double> a = data["radius"];
    for (int i =0; i< radial_flux.size(); i++) {
        radial_flux[i] /= (2.*M_PI*a.at(i));
    }
    this->addData("radial_flux",radial_flux);
    std::cout << "added radial flux"<< "\n" << std::flush;

    //    auto& kr = data["kr"]; // use addConductivities before!
    auto& kx = data["kx"]; // use addConductivities before!
    std::vector<double> axial_flux(segments.size());
    for (size_t i=0; i<axial_flux.size(); i++) {
        //        double age = simTime - data["creationTime"].at(i);
        //        int subType = (int)data["subType"].at(i);
        //        int organType = (int)data["organType"].at(i);
        auto s = segments[i];
        auto n1 = nodes.at(s.x);
        auto n2 = nodes.at(s.y);
        auto v = n2.minus(n1);
        double l = v.length();
        v.normalize();
        //        double p_s = sx.at(i);
        double dpdz0;
        //        if (a.at(i)*kr.at(i)>1.e-16) {
        //            double tau = std::sqrt(2 * a * M_PI * kr.at(i) / kx.at(i));  // cm-2 // TODO exact value is missing, see xylem_flux.py, axial_flux
        //        } else {
        //        }
        dpdz0 = (rx[s.y] - rx[s.x]) / l;
        axial_flux.at(i) = kx.at(i) * (dpdz0 + v.z);
    }
    this->addData("axial_flux",axial_flux);
}


/**
 * Stores the soil cell index for each segment as "cell_id".
 *
 * @param plant  mapped segment structure containing the seg2cell mapping
 */
void SegmentAnalyser::addCellIds(const MappedSegments& plant)
{
    std::vector<double> cell_id(segments.size());
    for (size_t i=0; i<segments.size(); i++) {
        cell_id[i] = plant.seg2cell.at(i);
    }
    this->addData("cell_id",cell_id);
}

/**
 * Returns a named parameter per segment.
 *
 * "creationTime" and "radius" are always available. Additional data can be registered via addData().
 * "length", "surface", and "volume" are computed on the fly from geometry.
 * All other names are forwarded to the owning organ; if the organ pointer is expired,
 * @p def is returned, or an exception is thrown if @p def is NaN.
 *
 * @param name  parameter name
 * @param def   fallback value for segments without a valid organ pointer; throws if NaN
 * @return      per-segment values, one entry per segment
 */
std::vector<double> SegmentAnalyser::getParameter(std::string name, double def) const
{
    if (data.count(name)) { // first, check data
        return data.at(name);
    }
    std::vector<double> d(segments.size()); // make return vector
    if (name == "length") {
        for (size_t i=0; i<d.size(); i++) {
            d.at(i) = getSegmentLength(i);
        }
        return d;
    }
    if (name == "surface") {
        for (size_t i=0; i<d.size(); i++) {
            d.at(i) = 2*data.at("radius")[i]*M_PI*getSegmentLength(i);
        }
        return d;
    }
    if (name == "volume") {
        for (size_t i=0; i<d.size(); i++) {
            double a = data.at("radius")[i];
            d.at(i) = a*a*M_PI*getSegmentLength(i);
        }
        return d;
    }
    for (size_t i=0; i<segO.size(); i++) { // else pass to Organs
        if (!segO.at(i).expired()) {
            d.at(i) = segO.at(i).lock()->getParameter(name);
        } else { // in case the segment has no origin
            if (std::isnan(def)) {
                throw std::invalid_argument("SegmentAnalyser::getParameter: segment origin expired (segment has no owner), "
                    "for segment index " + std::to_string(i) + ", parameter name "+name);
            } else {
                d.at(i) = def;
            }
        }
    }
    return d;
}

/**
 * Returns the Euclidean length of segment @p i.
 *
 * @param i  segment index
 * @return   length [cm]
 */
double SegmentAnalyser::getSegmentLength(int i) const
{
    Vector2i s = segments.at(i);
    return (nodes.at(s.x).minus(nodes.at(s.y))).length();
}

/**
 * Crops segments to those lying inside @p geometry.
 * Segments crossing the boundary are split exactly at the intersection point.
 * Unused nodes are retained; call pack() afterwards to remove them.
 *
 * @param geometry  signed distance function defining the domain (dist ≤ 0 = inside)
 */
void SegmentAnalyser::crop(std::shared_ptr<SignedDistanceFunction> geometry)
{
    //std::cout << "cropping " << segments.size() << " segments...";
    std::vector<Vector2i> seg;
    std::vector<std::weak_ptr<Organ>> sO;
    std::map<std::string, std::vector<double>> ndata;
    for (size_t i=0; i<segments.size(); i++) {
        auto s = segments.at(i);
        bool x_ = geometry->getDist(nodes.at(s.x))<=0; // in?
        bool y_ = geometry->getDist(nodes.at(s.y))<=0; // in?
        if (x_ && y_) { //segment is inside
            seg.push_back(s);
            if (segO.size()>0) {
                sO.push_back(segO.at(i));
            }
            for(auto iter = data.begin(); iter != data.end(); ++iter) {
                std::string key =  iter->first;
                ndata[key].push_back(data[key].at(i));
            }
        } else if ((x_==false) && (y_==false)) { // segment is outside
        } else { // one node is inside, one outside
            if (x_==false) { // swap indices
                int ind = s.x;
                s.x = s.y;
                s.y = ind;
            }
            Vector3d newnode = cut(nodes[s.x], nodes[s.y], geometry);
            nodes.push_back(newnode); // add new segment
            Vector2i newseg(s.x,nodes.size()-1);
            seg.push_back(newseg);
            if (segO.size()>0) {
                sO.push_back(segO.at(i));
            }
            for(auto iter = data.begin(); iter != data.end(); ++iter) { // copy data
                std::string key =  iter->first;
                ndata[key].push_back(data[key].at(i));
            }
        }
    }
    segments = seg;
    segO  = sO;
    data = ndata;
    //std::cout << " cropped to " << segments.size() << " segments " << "\n";
}

/**
 * Crops to the box [-xx/2, xx/2] × [-yy/2, yy/2] × [-zz, 0]. @see SDF_PlantBox
 */
void SegmentAnalyser::cropDomain(double xx, double yy, double zz)
{
    crop(std::make_shared<SDF_PlantBox>(xx,yy,zz));
}

/**
 * Keeps only segments where parameter @p name is in the closed interval [@p min, @p max].
 * @see SegmentAnalyser::getParameter
 *
 * @param name  parameter name
 * @param min   lower bound (inclusive)
 * @param max   upper bound (inclusive)
 */
void SegmentAnalyser::filter(std::string name, double min, double max)
{
    std::vector<double> d_ = getParameter(name);
    std::vector<Vector2i> seg;
    std::vector<std::weak_ptr<Organ>> sO;
    std::map<std::string, std::vector<double>> ndata;
    for (size_t i=0; i<segments.size(); i++) {
        if ((d_.at(i)>=min) && (d_.at(i)<=max)) {
            seg.push_back(segments.at(i));
            if (segO.size()>0) {
                sO.push_back(segO.at(i));
            }
            for(auto iter = data.begin(); iter != data.end(); ++iter) {
                std::string key =  iter->first;
                ndata[key].push_back(data[key].at(i));
            }
        }
    }
    segments = seg;
    segO  = sO;
    data = ndata;
}

/**
 * Keeps only segments where parameter @p name equals @p value.
 * @see SegmentAnalyser::getParameter
 *
 * @param name   parameter name
 * @param value  required value
 */
void SegmentAnalyser::filter(std::string name, double value)
{
    std::vector<double> d_ = getParameter(name);
    std::vector<Vector2i> seg;
    std::vector<std::weak_ptr<Organ>> sO;
    std::map<std::string, std::vector<double>> ndata;
    for (size_t i=0; i<segments.size(); i++) {
        if (d_.at(i)==value) {
            seg.push_back(segments.at(i));
            if (segO.size()>0) {
                sO.push_back(segO.at(i));
            }
            for(auto iter = data.begin(); iter != data.end(); ++iter) {
                std::string key =  iter->first;
                ndata[key].push_back(data[key].at(i));
            }
        }
    }
    segments = seg;
    segO  = sO;
    data = ndata;
}

/**
 * Re-indexes nodes by order of first occurrence in the segment list and removes unreferenced nodes.
 * Call after crop() or filter() to free memory.
 */
void SegmentAnalyser::pack() {
    std::vector<double> ni(nodes.size());
    std::fill(ni.begin(),ni.end(), -1.);
    std::vector<Vector3d> newnodes;
    for (auto& s:segments) {
        if (ni.at(s.x) == -1.) { // the node is new
            newnodes.push_back(nodes.at(s.x));
            ni.at(s.x) = newnodes.size()-1; // set index of the new node
        }
        s.x = ni.at(s.x);
        if (ni.at(s.y) == -1.) { // the node is new
            newnodes.push_back(nodes.at(s.y));
            ni.at(s.y) = newnodes.size()-1; // set index of the new node
        }
        s.y = ni.at(s.y);
    }
    // std::cout << "pack(): nodes: " << nodes.size() << " -> " << newnodes.size() << ", " << double(newnodes.size())/double(nodes.size()) << " \n";
    nodes = newnodes; // kabum!
}


/**
 * Returns the component-wise minimum of all node coordinates (lower corner of the bounding box).
 * Result is not cached.
 */
Vector3d SegmentAnalyser::getMinBounds() {
    Vector3d min_ = Vector3d(1.e9, 1.e9, 1.e9); // much
    for (const auto& n : nodes) {
        if (n.x<min_.x) {
            min_.x = n.x;
        }
        if (n.y<min_.y) {
            min_.y = n.y;
        }
        if (n.z<min_.z) {
            min_.z = n.z;
        }
    }
    return min_;
}

/**
 * Returns the component-wise maximum of all node coordinates (upper corner of the bounding box).
 * Result is not cached.
 */
Vector3d SegmentAnalyser::getMaxBounds() {
    Vector3d max_ = Vector3d(-1.e9, -1.e9, -1.e9); // litte
    for (const auto& n : nodes) {
        if (n.x>max_.x) {
            max_.x = n.x;
        }
        if (n.y>max_.y) {
            max_.y = n.y;
        }
        if (n.z>max_.z) {
            max_.z = n.z;
        }
    }
    return max_;
}

/**
 * Numerically bisects the segment (@p in, @p out) to find its intersection with @p geometry.
 *
 * @param in        node known to be inside the domain (dist ≤ 0)
 * @param out       node known to be outside the domain (dist ≥ 0)
 * @param geometry  signed distance function of the domain boundary
 * @param eps       convergence tolerance [cm] (default 1e-6)
 * @return          intersection point
 */
Vector3d SegmentAnalyser::cut(Vector3d in, Vector3d out, const std::shared_ptr<SignedDistanceFunction>& geometry, double eps)
{
    // std::cout << out.toString() << ", " << geometry->getDist(out) << "\n";
    assert(geometry->getDist(in)<=0  && "SegmentAnalyser::cut(): in is not within domain" );
    assert(geometry->getDist(out)>=0 && "SegmentAnalyser::cut(): out is not outside domain" );
    Vector3d c =  in.plus(out).times(0.5); // mid
    if (in.minus(out).length() < eps) {
        return c;
    }
    if (geometry->getDist(c)<0) { // in
        // std::cout << geometry->getDist(c) << " c, out ";
        return cut(c, out, geometry, eps);
    } else { // out
        // std::cout  << geometry->getDist(c) << " in, c  ";
        return cut(in, c, geometry, eps);
    }
}

/**
 * Returns the sum of parameter @p name over all segments.
 */
double SegmentAnalyser::getSummed(std::string name) const {
    std::vector<double> v_ = getParameter(name);
    return std::accumulate(v_.begin(), v_.end(), 0.0);
}

/**
 * Returns the approximate sum of parameter @p name for segments whose midpoint lies inside geometry @p g.
 * For an exact result, crop to the geometry first and then call getSummed(name).
 *
 * @param name  parameter name
 * @param g     bounding geometry
 * @return      approximate sum within @p g
 */
double SegmentAnalyser::getSummed(std::string name, std::shared_ptr<SignedDistanceFunction> g) const {
    std::vector<double> data = getParameter(name);
    double v = 0;
    for (size_t i=0; i<segments.size(); i++) {
        Vector2i s = segments.at(i);
        Vector3d mid = nodes.at(s.x).plus(nodes.at(s.y)).times(0.5);
        if (g->getDist(mid)<0) {
            v += data.at(i);
        }
    }
    return v;
}

/**
 * Maps segments into the periodic domain [-xx/2, xx/2) × [-yy/2, yy/2).
 * Segments crossing the boundary are split; both resulting pieces are kept.
 *
 * @param xx  period in x direction [cm]
 * @param yy  period in y direction [cm]
 */
void SegmentAnalyser::mapPeriodic(double xx, double yy) {
    mapPeriodic_(xx, Vector3d(1.,0,0), 1.e-6);
    mapPeriodic_(yy, Vector3d(0,1.,0), 1.e-6);
}

/**
 * Maps segments into a 1-D periodic domain of width @p xx along @p axis.
 * Segments crossing the boundary are split; nodes are mapped into [-xx/2, xx/2).
 *
 * @param xx    period width [cm]
 * @param axis  unit vector of the periodic direction (e.g. Vector3d(1,0,0))
 * @param eps   boundary tolerance [cm]
 */
void SegmentAnalyser::mapPeriodic_(double xx, Vector3d axis, double eps) {
    /* 1. split segments at the boundaries */
    std::vector<Vector2i> seg;
    std::vector<std::weak_ptr<Organ>> sO;
    std::map<std::string, std::vector<double>> ndata;

    for (size_t i=0; i<segments.size(); i++) {
        auto s = segments.at(i);
        auto n1 = nodes.at(s.x);
        auto n2 = nodes.at(s.y);
        int p1 = floor((n1.times(axis)+xx/2.)/xx);
        int p2 = floor((n2.times(axis)+xx/2.)/xx);
        if (p1 == p2) { //same periodicity index, do nothing [0,xx)
            seg.push_back(s);
            if (segO.size()>0) { // if used
                sO.push_back(segO.at(i));
            }
            for(auto iter = data.begin(); iter != data.end(); ++iter) { // copy data
                std::string key =  iter->first;
                ndata[key].push_back(data[key].at(i));
            }
        } else { // otherwise split
            if (p1 > p2) { // sort
                int ind = s.x;
                s.x = s.y;
                s.y = ind;
                n1 = nodes.at(s.x);
                n2 = nodes.at(s.y);
                ind = p1;
                p1 = p2;
                p2 = ind;
            }
            // now: p1x < p2x, n1.x < n2.x
            double theta = (p2*xx - (n1.times(axis)+xx/2) )/(n2.times(axis)-n1.times(axis));
            //            std::cout << " n1.x " <<  n1.x << " n2.x " <<  n2.x << ", length: " << n2.x-n1.x << ", theta: "<< theta << " p1: "<< p1 <<", p2: " << p2 << "\n";
            auto v = n2.minus(n1);
            auto x = n1.plus(v.times(theta)); // cutting point
            auto x0 = x.minus(axis.times(eps)); // less
            auto x1 = x.plus(axis.times(eps)); // greater
            int c = 0; // segments added
            if ((n1.minus(x0)).length()>eps) { // if inside segment is large enough
                nodes.push_back(x0); // add node and new add segment
                seg.push_back(Vector2i(s.x,nodes.size()-1));
                c++;
            }
            if ((n2.minus(x1)).length()>eps) { // if outside segment is large enough
                nodes.push_back(x1); // add node and new add segment
                seg.push_back(Vector2i(nodes.size()-1, s.y));
                c++;
            }
            for (int j=0; j<c; j++) { // copy attached data
                if (segO.size()>0) { // if used
                    sO.push_back(segO.at(i));
                }
                for(auto iter = data.begin(); iter != data.end(); ++iter) {
                    std::string key =  iter->first;
                    ndata[key].push_back(data[key].at(i));
                }
            }
        }
    }
    segments = seg;
    segO  = sO;
    data = ndata;
    /* 2. map points to period [-xx/2, xx/2] */
    for (auto& n : nodes) {
        n = n.minus(axis.times(floor((n.times(axis)+xx/2.)/xx)*xx));
    }
}

/**
 * Projects nodes to the x–z plane: (x, y, z) → (√(x²+y²), 0, z).
 * The sign of the radial coordinate is taken from the first segment of each branch.
 */
void SegmentAnalyser::map2D() {

    for (size_t i=0; i<segments.size(); i++) {
        if (nodes[segments[i].x].y!=0) { // not mapped before (initial segment of base root)
            double x0 = nodes[segments[i].y].x; // first segment decides signum of branch
            double x = nodes[segments[i].x].x; // node 1
            double y = nodes[segments[i].x].y;
            nodes[segments[i].x].x = (-1*(x0<=0) + 1*(x0>0))*std::sqrt(x*x+y*y);
            nodes[segments[i].x].y = 0.;
            x = nodes[segments[i].y].x; // node 2
            y = nodes[segments[i].y].y;
            nodes[segments[i].y].x = (-1*(x0<=0) + 1*(x0>0))*std::sqrt(x*x+y*y);
            nodes[segments[i].y].y = 0.;
        } else { // map only node 2
            double x = nodes[segments[i].y].x;
            double y = nodes[segments[i].y].y;
            double r = std::sqrt(x*x+y*y);
            double x0;
            if (r>0.1) {
                x0 = nodes[segments[i].x].x;
            } else {
                x0 = nodes[segments[i].y].x;
            }
            int sgn = (-1*(x0<=0) + 1*(x0>0));
            nodes[segments[i].y].x = sgn*r;
            nodes[segments[i].y].y = 0.;
        }
    }
}

/**
 * Returns the unique set of organs that own the segments.
 * Optionally filtered by organ type @p ot (pass -1 to return all types).
 */
std::vector<std::shared_ptr<Organ>> SegmentAnalyser::getOrgans(int ot) const
{
    std::set<std::shared_ptr<Organ>> rootset;  // praise the stl
    for (auto o : segO) {
		if ((ot<0) || (ot == o.lock()->organType())) {
			rootset.insert(o.lock());
		}
    }
    return std::vector<std::shared_ptr<Organ>>(rootset.begin(), rootset.end());
}

/**
 * Returns the number of distinct organs in the segment set.
 */
int SegmentAnalyser::getNumberOfOrgans() const
{
    return getOrgans().size();
}

/**
 * Projects segments onto an image plane (perspective projection).
 * TODO: unfinished and untested.
 *
 * @param pos  camera position
 * @param ons  orthonormal system; row 0 is the optical axis, rows 1–2 span the image plane
 * @param fl   focal length [cm]
 * @return     projected segments in the x–y plane (z = 0)
 */
SegmentAnalyser SegmentAnalyser::foto(const Vector3d& pos, const Matrix3d& ons, double fl) const
{
    SegmentAnalyser f(*this); // copy
    for (auto& n : f.nodes) { // translate
        n = n.minus(pos);
    }
    Matrix3d m = ons.inverse(); // rotate
    for (auto& n : f.nodes) {
        n = m.times(n);
    }
    //	// crop to objects in front of the camera
    Vector3d o(0.,0.,0.);
    Vector3d plane(0.,0., -1);
    f.crop(std::make_shared<SDF_HalfPlane>(o,plane));
    f.pack();
    // project
    for (auto& a : f.nodes) {
        a = a.times(fl/(-plane.times(a)));
        a.z = 0;
    }
    // final image crop // TODO --> d = sqrt(2) cm
    f.crop(std::make_shared<SDF_PlantBox>(1.,1.,2.));
    f.pack();
    return f;
}

/**
 * Returns a new SegmentAnalyser containing only the segments that intersect @p plane
 * (e.g. for trench sampling). Unused nodes are removed via pack().
 *
 * @param plane  signed distance function of the cutting plane
 */
SegmentAnalyser SegmentAnalyser::cut(const SignedDistanceFunction& plane) const
{
    SegmentAnalyser f;
    f.nodes = nodes; // copy all nodes
    for (size_t i=0; i<segments.size(); i++) {
        Vector2i s = segments.at(i);
        Vector3d n1 = nodes.at(s.x);
        Vector3d n2 = nodes.at(s.y);
        double d = plane.getDist(n1)*plane.getDist(n2);
        if (d<=0) { // one is inside, one is outside
            f.segments.push_back(s);
            f.segO.push_back(segO.at(i));
            for(auto iter = data.begin(); iter != data.end(); ++iter) {
                std::string key =  iter->first;
                f.data[key].push_back(data.at(key)[i]);
            }
        }
    }
    f.pack(); // delete unused nodes
    return f;
}

/**
 * Returns a vertical distribution of parameter @p name in @p n equal layers between @p bot and @p top.
 *
 * @param name   parameter name
 * @param top    top of the domain [cm] (typically 0)
 * @param bot    bottom of the domain [cm] (e.g. -100)
 * @param n      number of layers
 * @param exact  if true, crops exactly to each layer; if false, uses segment midpoints (faster)
 * @return       vector of length @p n with the summed parameter per layer
 */
std::vector<double> SegmentAnalyser::distribution(std::string name, double top, double bot, int n, bool exact) const
{
    std::vector<double> d(n);
    double dz = (top-bot)/double(n);
    assert(dz > 0 && "SegmentAnalyser::distribution: top must be larger than bot" );
    std::shared_ptr<SignedDistanceFunction> layer = std::make_shared<SDF_PlantBox>(1e100,1e100,dz);
    for (int i=0; i<n; i++) {
        Vector3d t(0,0,top-i*dz);
        auto g = std::make_shared<SDF_RotateTranslate>(layer,t);
        if (exact) {
            SegmentAnalyser a(*this); // copy everything
            a.crop(g); // crop exactly
            d.at(i) = a.getSummed(name);
        } else {
            d.at(i) = this->getSummed(name, g);
        }
    }
    return d;
}

/**
 * Returns a vertical distribution as @p n SegmentAnalyser objects, one per layer (cropped exactly).
 *
 * @param top  top of the domain [cm]
 * @param bot  bottom of the domain [cm]
 * @param n    number of layers
 * @return     vector of length @p n
 */
std::vector<SegmentAnalyser> SegmentAnalyser::distribution(double top, double bot, int n) const
{
    std::vector<SegmentAnalyser> d(n);
    double dz = (top-bot)/double(n);
    assert(dz > 0 && "SegmentAnalyser::distribution: top must be larger than bot" );
    std::shared_ptr<SignedDistanceFunction> layer = std::make_shared<SDF_PlantBox>(1.e100, 1.e100, dz);
    for (int i=0; i<n; i++) {
        Vector3d t(0,0,top-i*dz);
        SegmentAnalyser a = SegmentAnalyser(*this); // copy everything
        a.crop(std::make_shared<SDF_RotateTranslate>(layer,t)); // crop exactly
        d.at(i) = a;
    }
    return d;
}

/**
 * Returns a 2-D (x–z) distribution of parameter @p name on an @p n × @p m grid.
 *
 * @param name   parameter name
 * @param top    top of the vertical domain [cm]
 * @param bot    bottom of the vertical domain [cm]
 * @param left   left boundary in x [cm]
 * @param right  right boundary in x [cm]
 * @param n      number of vertical layers
 * @param m      number of horizontal columns
 * @param exact  if true, crops exactly; if false, uses segment midpoints
 * @return       n × m matrix of summed parameter values
 */
std::vector<std::vector<double>> SegmentAnalyser::distribution2(std::string name, double top, double bot, double left, double right, int n, int m, bool exact) const
{
    std::vector<std::vector<double>> d(n);
    double dz = (top-bot)/double(n);
    assert(dz > 0 && "SegmentAnalyser::distribution2: top must be larger than bot" );
    double dx = (right-left)/double(m);
    assert(dx > 0 && "SegmentAnalyser::distribution2: right must be larger than left" );
    std::shared_ptr<SignedDistanceFunction> layer = std::make_shared<SDF_PlantBox>(dx, 1.e9, dz);
    for (int i=0; i<n; i++) {
        std::vector<double> row(m); // m columns
        for (int j=0; j<m; j++) {
            Vector3d t(left+(j+0.5)*dx,0.,top-i*dz); // box is [-x/2,-y/2,0] - [x/2,y/2,-z]
            auto g = std::make_shared<SDF_RotateTranslate>(layer,t);
            if (exact) {
                SegmentAnalyser a(*this); // copy everything
                a.crop(g); // crop exactly
                row.at(j) = a.getSummed(name);
            } else {
                row.at(j) = this->getSummed(name, g);
            }
        }
        d.at(i)=row; // store the row (n rows)
    }
    return d;
}

/**
 * Returns a 2-D (x–z) grid of SegmentAnalyser objects, one per cell (cropped exactly).
 *
 * @param top    top of the vertical domain [cm]
 * @param bot    bottom of the vertical domain [cm]
 * @param left   left boundary in x [cm]
 * @param right  right boundary in x [cm]
 * @param n      number of vertical layers
 * @param m      number of horizontal columns
 * @return       n × m matrix of SegmentAnalyser objects
 */
std::vector<std::vector<SegmentAnalyser>> SegmentAnalyser::distribution2(double top, double bot, double left, double right, int n, int m) const
{
    std::vector<std::vector<SegmentAnalyser>> d(n);
    double dz = (top-bot)/double(n);
    assert(dz > 0 && "SegmentAnalyser::distribution2: top must be larger than bot" );
    double dx = (right-left)/double(m);
    assert(dx > 0 && "SegmentAnalyser::distribution2: right must be larger than left" );
    std::shared_ptr<SignedDistanceFunction> layer = std::make_shared<SDF_PlantBox>(dx, 1.e4, dz);
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            Vector3d t(left+(j+0.5)*dx,0.,top-i*dz); // box is [-x/2,-y/2,0] - [x/2,y/2,-z]
            SegmentAnalyser a(*this); // copy everything
            a.crop(std::make_shared<SDF_RotateTranslate>(layer,t)); // crop exactly
            d.at(i).push_back(a);
        }
    }
    return d;
}

/**
 * Attaches named per-segment data accessible via getParameter() and written by writeVTP().
 * Node-sized vectors are automatically converted to segment data using the tip node index.
 *
 * @param name    parameter name
 * @param values  per-segment or per-node values; size must match segments or nodes
 */
void SegmentAnalyser::addData(std::string name, std::vector<double> values)
{
    if (values.size()== segments.size()) {
        data[name] = values;
    } else if (values.size()==nodes.size()) { // convert node to segment data
        std::vector<double> d;
        d.reserve(segments.size());
        for (int i = 0; i<segments.size(); i++) {
            d.push_back(values.at(segments[i].y));
        }
        data[name] = d;
    } else {
        throw std::invalid_argument("SegmentAnalyser::addData: parameter values has wrong size.");
    }
}

/**
 * Writes simulation results to @p name; format is determined by the file extension.
 * Supported: ".vtp" (VTK PolyData), ".txt" (RootBox format), ".dgf" (DuMux format).
 * pack() is called before writing.
 *
 * @param name   output file path (extension must be lower-case)
 * @param types  cell-data fields to include in VTP output
 */
void SegmentAnalyser::write(std::string name, std::vector<std::string> types)
{
    this->pack(); // a good idea before writing any file
    std::ofstream fos;
    fos.open(name.c_str());
    std::string ext = name.substr(name.size()-3,name.size()); // pick the right writer
    if (ext.compare("vtp")==0) {
        std::cout << "writing VTP: " << name << "\n" << std::flush;
        this->writeVTP(fos, types);
    } else if (ext.compare("txt")==0)  {
        std::cout << "writing text file for Matlab import: "<< name << "\n"<< std::flush;
        writeRBSegments(fos);
    } else if (ext.compare("dgf")==0)  {
        std::cout << "writing dgf file: "<< name << "\n"<< std::flush;
        writeDGF(fos);
    } else {
        throw std::invalid_argument("SegmentAnalyser::write: Unknown file type");
    }
    fos.close();
}

/**
 * Writes a VTK PolyData (.vtp) file with the specified cell-data fields.
 *
 * @param os     output stream
 * @param types  names of per-segment data fields to include as CellData
 */
void SegmentAnalyser::writeVTP(std::ostream & os, std::vector<std::string> types) const
{
    os << "<?xml version=\"1.0\"?>";
    os << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    os << "<PolyData>\n";
    os << "<Piece NumberOfLines=\""<< segments.size() << "\" NumberOfPoints=\""<< nodes.size()<< "\">\n";
    // data (CellData)
    os << "<CellData Scalars=\" CellData\">\n";
    for (auto name : types) {
        std::vector<double> data = getParameter(name, -1.);
        os << "<DataArray type=\"Float32\" Name=\"" << name << "\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
        for (const auto& t : data) {
            os << t << " ";
        }
        os << "\n</DataArray>\n";
    }
    os << "\n</CellData>\n";
    // nodes (Points)
    os << "<Points>\n"<<"<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\" >\n";
    for (const auto& n : nodes) {
        os << n.x << " "<< n.y <<" "<< n.z<< " ";
    }
    os << "\n</DataArray>\n"<< "</Points>\n";
    // segments (Lines)
    os << "<Lines>\n"<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
    for (const auto& s : segments) {
        os << s.x << " " << s.y << " ";
    }
    os << "\n</DataArray>\n"<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
    for (size_t i=0; i<segments.size(); i++) {
        os << 2*i+2 << " ";
    }
    os << "\n</DataArray>\n";
    os << "\n</Lines>\n";
    //
    os << "</Piece>\n";
    os << "</PolyData>\n" << "</VTKFile>\n";
}

/**
 * Writes segments in RootBox text format, sorted by creation time.
 * Mimics the Matlab getSegments() script from RootBox.
 *
 * @param os  output stream
 */
void SegmentAnalyser::writeRBSegments(std::ostream & os) const
{
    auto ctime = getParameter("creationTime");
    std::vector< std::tuple<double,  // creation time 0
    int,  // node1ID 1
    int,  // node2ID 2
    int,  // branchID 3
    double, // x1 4
    double, // y1 5
    double, // z1 6
    double, // x2 7
    double, // y2 8
    double, // z2 9
    double, // radius 10
    double, // age 11
    int,    // subtype 12
    int,    //organ 13
    int    //organid 14
    > > ctime_seg;
    os << "node1ID node2ID branchID x1 y1 z1 x2 y2 z2 radius time age type organ parent_node_id \n";
    for (size_t i=0; i<segments.size(); i++) {
        std::shared_ptr<Organ> o = segO.at(i).lock();
        int organ=o->organType();
        int branchnumber = o->getId();
        Vector2i s = segments.at(i);
        Vector3d n1 = nodes.at(s.x);
        Vector3d n2 = nodes.at(s.y);
        int x = s.x;
        int y = s.y;
        double radius = o->getParameter("radius");
        double time = ctime[i];
        double age = o->getParameter("age");
        int subType = o->getParameter("subType");
        int id = o->getParameter("id");
        ctime_seg.push_back( std::make_tuple(time, x, y, branchnumber, n1.x, n1.y, n1.z, n2.x, n2.y, n2.z, radius, age, subType, organ, id  ));
        //        auto parent_id = std::find_if(segments.begin(), segments.end(), [y](const Vector2i s_){return s_.x == y; });
        //        if (parent_id != segments.end())
        //          int index = std::distance(begin(segments), parent_id);
        //        else
        //           int index = 0;
        //         std::string parent_node_id = std::to_string(index);

    }
    std::sort(ctime_seg.begin(), ctime_seg.end()); // sort based on the time
    for (size_t i=0; i<segments.size(); i++) { // write output
        //int parent_node_id = std::find(segments.y, segments.y, s.at(i).x);
        os << std::fixed << std::setprecision(4)<< std::get<1>(ctime_seg[i])
                                                     << " " << std::get<2>(ctime_seg[i])
                                                     << " " << std::get<3>(ctime_seg[i])
                                                     << " " << std::get<4>(ctime_seg[i])
                                                     << " " << std::get<5>(ctime_seg[i])
                                                     << " " << std::get<6>(ctime_seg[i])
                                                     << " " << std::get<7>(ctime_seg[i])
                                                     << " " << std::get<8>(ctime_seg[i])
                                                     << " " << std::get<9>(ctime_seg[i])
                                                     << " " << std::get<10>(ctime_seg[i])
                                                     << " " << std::get<0>(ctime_seg[i])
                                                     << " " << std::get<11>(ctime_seg[i])
                                                     << " " << std::get<12>(ctime_seg[i])
                                                     << " " << std::get<13>(ctime_seg[i])
                                                     << " " << std::get<14>(ctime_seg[i])<< " " <<" \n";
        //for debug//        os << std::fixed << std::setprecision(4)<< std::get<1>(ctime_seg[i]) << " " << std::get<2>(ctime_seg[i]) << " " << std::get<0>(ctime_seg[i]) <<" "  << "\n"; << branchnumber << " " << n1.x << " " << n1.y << " " << n1.z << " " << n2.x << " " << n2.y << " " << n2.z << " " << radius << " " << std::get<1>(ctime_seg[i])<< " " << age<<" " <<subType<< " " <<organ << " " <<" \n";

    }
}

/**
 * Writes segments in DGF format for use with DuMux.
 *
 * Column layout (IBG-3 convention):
 *  - 0: order (branching order, starting from 0; -1 for artificial shoots)
 *  - 1: brnID (unique organ id)
 *  - 2: surface [cm²]
 *  - 3: length [cm]
 *  - 4: radius [cm]
 *  - 5: kz [cm⁴ hPa⁻¹ d⁻¹] (set to 0; unknown to CPlantBox)
 *  - 6: kr [cm hPa⁻¹ d⁻¹] (set to 0; unknown to CPlantBox)
 *  - 7: emergence time [d]
 *  - 8: subType
 *  - 9: organType (seed=1, root=2, stem=3, leaf=4; -1 for artificial shoots)
 *
 * @param os  output stream
 */
void SegmentAnalyser::writeDGF(std::ostream & os) const
{
    os << "DGF" << std::endl;
    os << "Vertex" << std::endl;
    for (auto& n : nodes) {
        os << n.x/100 << " " << n.y/100 << " " << n.z/100  << std::endl;
    }
    os << "#" << std::endl;
    os << "SIMPLEX" << std::endl;
    // node1ID, node2ID, order, brnID, surf [cm2], length [cm], radius [cm],
    // kz [cm4 hPa-1 d-1], kr [cm hPa-1 d-1], emergence time [d], subType, organType
    os << "parameters 10 # id0, id1, order, branchId, surf[cm2], length[cm], radius[cm], "
        "kz[cm4 hPa-1 d-1], kr[cm hPa-1 d-1], emergence time [d], subType, organType " << std::endl;
    std::vector<double> radius = getParameter("radius");
    std::vector<double> length  = getParameter("length");
    std::vector<double> surface = getParameter("surface");
    std::vector<double> ctime = getParameter("creationTime");
    std::vector<double> subType = getParameter("subType", -1.); // -1 if no origin organ
    std::vector<double> organType = getParameter("organType", Organism::ot_stem); // artificial stem
    std::vector<double> order = getParameter("order", -1.); // -1 if no origin organ
    std::vector<double> organId =  getParameter("id", -1.); // -1 if no origin organ
    for (size_t i=0; i<segments.size(); i++) {
        Vector2i s = segments.at(i);
        os << s.x << " " << s.y << " " << order.at(i) <<  " " << organId.at(i) << " " << surface.at(i) << " "
            << length.at(i) << " " << radius.at(i) << " " << 0. << " " << 0. << " " << ctime.at(i) << " "
            << subType.at(i) << " " << organType.at(i) << std::endl;
    }
    os << "#" << std::endl;
    os << "BOUNDARYDOMAIN " << std::endl;
    os << "default 1" << std::endl;
    os << "#" << std::endl;
}

} // end namespace CPlantBox
