// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "SegmentAnalyser.h"

#include "Organ.h"
#include "Organism.h"
#include "MappedOrganism.h"
#include "XylemFlux.h"
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
 * Creates segments, with data "creationTime" and "radius"
 *
 * The rest of parameters, that might be needed are either a default value, or
 * if not defined an exception is thrown, see SegmentAnalyser::getParameter.
 *
 * @param nodes     a list of nodes
 * @param segment   each segment is represented by two node indices the segment, the line segment is given by nodes[segment.x] to nodes[segment.y].
 * @param segCT     creation time of each segment
 * @param radii     the segment radii
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
 * Copies the line segments representing the plant to the analysis class
 *
 * @param plant     the the organism that is analysed
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
 * Copies the line segments representing the plant to the analysis class
 *
 * @param plant     the the organism that is analysed
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
 * Adds all line segments from @param plant to the analysis
 */
void SegmentAnalyser::addSegments(const Organism& plant)
{
    addSegments(SegmentAnalyser(plant));
}

/**
 * Adds all line segments from the analyser @param a to this analysis.
 *
 * Ignores any user data from @param a.
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
 * Adds a single segment (e.g. an artificial shoot),
 *
 * @param seg       a single segment represented of two node indices of the analyser's node list.
 *                  care, some methods modify the node list order (e.g. pack)
 * @param ct        segment creation time
 * @param radius    segment radius
 * @param insert    inserts the segment at the first index (true), or appends it to the segment list (false)
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

void SegmentAnalyser::addAge(double simTime)
/**
 * adds "age" for vizualisation (cuts off negative values)
 */
{
    std::vector<double> age(segments.size());
    for (size_t i=0; i<age.size(); i++) {
        double a = simTime - data["creationTime"].at(i);
        age.at(i) = std::max(a,0.);
    }
    this->addData("age",age);
}

/**
 * Adds kr and kx to the user data for vizualisation ("kr", "kx"),
 *
 * @param rs    XylemFlux for determination of radial and axial conductivities (kr, and kx)
 */
void SegmentAnalyser::addConductivities(const XylemFlux& rs, double simTime, double kr_max, double kx_max)
{
    // std::cout << "creationTime " << data["creationTime"].size() << ", " << data["subType"].size() << ", " << data["organType"].size() << "\n";

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
 * Evaluates kr and kx using PlantHydraulicParameters
 * and adds kr and kx per segment to the user data for vizualisation ("kr", "kx"),
 *
 * @param rs    XylemFlux for determination of radial and axial conductivities (kr, and kx)
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
 * Adds radial and axial fluxes to the user data for vizualisation,
 * ("radial_flux" [cm3/cm2 / day], and "axial_flux" [cm3/day])
 *
 * use addConductivities before!
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


void SegmentAnalyser::addCellIds(const MappedSegments& plant)
{
    std::vector<double> cell_id(segments.size());
    for (size_t i=0; i<segments.size(); i++) {
        cell_id[i] = plant.seg2cell.at(i);
    }
    this->addData("cell_id",cell_id);
}

/**
 * Returns a specific parameter per root segment.
 *
 * The parameters "creationTime" and "radius" are stored in. Additional parameters can be added by using addData.
 *
 * Other parameters are retrieved from the segment's organ.
 * If the pointer of the organ is expired the default value is returned, or if default is nan an exception is thrown.
 *
 * @param name  parameter name
 * @param def	default parameter, if organ's origin is expired. If nan an exception is thrown.
 * @return      vector containing parameter value per segment
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
 * Returns the length of a segment
 *
 * @param i 	index of the segment
 * @return 		the length of segment i
 */
double SegmentAnalyser::getSegmentLength(int i) const
{
    Vector2i s = segments.at(i);
    return (nodes.at(s.x).minus(nodes.at(s.y))).length();
}

/**
 * Crops the segments with some geometry.
 * This is done exact, i.e. segments are cut in two at the geometry border.
 *
 * All nodes are kept, use pack() to remove unused nodes.
 *
 * @param geometry      signed distance function of the geometry
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
 *  Crops the segments to the domain [-xx/2, -yy/2, -zz] - [xx/2, yy/2, 0.]
 */
void SegmentAnalyser::cropDomain(double xx, double yy, double zz)
{
    crop(std::make_shared<SDF_PlantBox>(xx,yy,zz));
}

/**
 * Filters the segments to the ones, where the parameter is within the interval [min,max], @see SegmentAnalyser::getParameter,
 * i.e. all other segments are deleted.
 *
 * @param name  parameter name
 * @param min   minimal value
 * @param max   maximal value
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
 * Filters the segments to the ones, where parameter equals a specific value, @see SegmentAnalyser::getParameter,
 * i.e. all other segments are deleted.
 *
 * @param name      parameter name
 * @param value     parameter value of the segments that are kept
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
 * Sorts nodes by occurrence in the segment list, and deletes unused nodes.
 *
 * This can save a lot of memory, since SegmentAnalyser::crop and SegmentAnalyser::filter
 * only delete segments, not unused nodes
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
 * Calculates the minimum of node coordinates
 * (e.g. minimum corner of bounding box)
 * value not cached
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
 * Calculates the maximum of node coordinates
 * (e.g. maximum corner of bounding box)
 * value not cached
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
 * Numerically computes the intersection point
 *
 * @param in       the node within the domain
 * @param out      the node outside of the domain
 * @param geometry signed distance function of the geometry
 * @param eps      accuracy of intersection (default = 1.e-6 cm)
 * @return         the intersection point
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
 * @return The sum of parameter @param name
 */
double SegmentAnalyser::getSummed(std::string name) const {
    std::vector<double> v_ = getParameter(name);
    return std::accumulate(v_.begin(), v_.end(), 0.0);
}

/**
 * Returns the sum parameter of parameter @param name, within geometry @param g based
 * on the segment mid point (i.e. not exact).
 *
 * To sum up exactly, crop to the geometry, and then run SegmentAnalyser::getSummed(name).
 *
 * @return Approximated sum of parameter @param name within the geometry @param g
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
 * A unconfined rootsystem is mapped into a periodic domain with period xx.
 *
 * Segments crossing the periodic boundary are split, nodes are mapped into [ [-xx/2, xx/2), [-yy/2, yy/2) ]
 *
 * @param xx    period in x direction
 * @param xx    period in y direction
 */
void SegmentAnalyser::mapPeriodic(double xx, double yy) {
    mapPeriodic_(xx, Vector3d(1.,0,0), 1.e-6);
    mapPeriodic_(yy, Vector3d(0,1.,0), 1.e-6);
}

/**
 * A unconfined rootsystem is mapped into a periodic domain with period xx.
 *
 * Segments crossing the periodic boundary are split, nodes are mapped into [-xx/2, xx/2)
 *
 * @param xx    width
 * @param axis  unit vector in the periodic direction (e.g. Vector3d(1,0,0), or Vector(0,1,0)).
 * @param eps   accuracy at boundary
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
 * Maps the 3d coordinates to the x-z plan (sqrt(x2+y2), 0., z)
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
 * @return The origin's of the segments, i.e. the organ's where the segments are part of (unique, no special ordering)
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
 * @return The number of different organs
 */
int SegmentAnalyser::getNumberOfOrgans() const
{
    return getOrgans().size();
}

/**
 * Projects the segments to an image plane (todo verify this code, its unfinished)
 *
 * @param pos       position of camera
 * @param ons       orthonormal system, row 1 is orthogonal to the image plane given by [row 2,row 3]
 * @param fl        focal length, alpha = 2*arctan(d/(2*fl)), were alpha is the angle of field, and d the image diagonal
 *
 * @return The image segments in the x-y plane (z=0)
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
 * Keeps the segments that intersect with a plane  (e.g. simulating trenches).
 *
 * @param plane 	half plane
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
 * Creates a vertical distribution of the parameter @param name.
 *
 * @param name      parameter name
 * @param top       vertical top position (cm) (normally = 0)
 * @param bot       vertical bot position (cm) (e.g. = -100 cm)
 * @param n         number of layers, each with a height of (top-bot)/n
 * @param exact     calculates the intersection with the layer boundaries (true), only based on segment midpoints (false)
 *                  (TODO) intersections do not work with user data
 * @return          vector of size @param n containing the summed parameter in this layer
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
 *  Creates a vertical distribution, of the geometry.
 *
 * @param top       vertical top position (cm) (normally = 0)
 * @param bot       vertical bot position (cm) (e.g. = -100 cm)
 * @param n         number of layers, each with a height of (top-bot)/n
 * @return          vector of size @param n containing an Analysis object of the layers (cropped exact)
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
 *  Creates a two-dimensional distribution of the parameter @param name.
 *
 * @param name      parameter name
 * @param top       vertical top position (cm) (normally = 0)
 * @param bot       vertical bot position (cm) (e.g. = -100 cm)
 * @param left      left along x-axis (cm)
 * @param right     right along x-axis (cm)
 * @param n         number of layers, each with a height of (top-bot)/n
 * @param m 		number of horizontal grid elements (each with length of (right-left)/m)
 * @param exact     calculates the intersection with the layer boundaries (true), only based on segment midpoints (false)
 * @return          vector of size @param n containing the summed parameter in this layer
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
 *  Creates a vertical distribution
 *
 * @param top       vertical top position (cm) (normally = 0)
 * @param bot       vertical bot position (cm) (e.g. = -100 cm)
 * @param left      left along x-axis (cm)
 * @param right     right along x-axis (cm)
 * @param n         number of layers, each with a height of (top-bot)/n
 * @param m 		number of horizontal grid elements (each with length of (right-left)/m)
 * @return          vector of size @param n containing the summed parameter in this layer
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
 * Adds user data that can be accessed by SegmentAnalyser::getParameter, and that can be written to the VTP file
 * (e.g. used to add simulation results like xylem pressure to the output).
 *
 * @param name      parameter name
 * @param values    segment or node data, node data are converted to segment data
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
 * Exports the simulation results with the type from the file extension in name (must be lower case)
 * Currently its possible to write "vtp", "txt", or "dgf" files.
 *
 * SegmentAnalyser::pack() is called before writing the file, i.e. nodes are sorted.
 *
 * @param name      file name e.g. output.vtp
 * @param types 	Optionally, for vtp we can determine the cell data by a vector of parameter names
 *                  (default = { "radius", "subType", "creationTime", "organType" })
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
 * Writes a VTP file with @param types data per segment.
 *
 * @param os        a file out stream
 * @param types     parameter names of the cell data  (default = { "radius", "subType", "creationTime", "organType" })
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
 * Writes the (line)segments of the root system, and
 * mimics the Matlab script getSegments() of RootBox
 *
 * @param os      a file out stream
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
 * Writes the (line)segments of the root system in dgf format used by DuMux
 *
 * For parameters, the IBG-3 default is:
 * 0 order (starting from 0)
 * 1 brnID (the unique organ id, for RootSystem starting at 1 because seed has id 0) ,
 * 2 surf [cm2],
 * 3 length [cm],
 * 4 radius [cm],
 * 5 kz [cm4 hPa-1 d-1],  axial root conductivity (unknown to CPlantBox) is set to 0.
 * 6 kr [cm hPa-1 d-1],   radial root conductivity (unknown to CPlantBox) is set to 0.
 * 7 emergence time [d],
 * 8 subType, (normally, starting from 1, type is set in the parameter xml files)
 * 9 organType (seed=1, root=2 ,stem=3, leaf=4)
 *
 * For artificial shoot segments (that can be added with SegmentAnalyser::addSegment)
 * order = -1
 * subType = -1
 * organType = -1.
 *
 * @param os      typically a file out stream
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
