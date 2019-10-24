// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef SEGMENTANALYSER_H_
#define SEGMENTANALYSER_H_

#include "sdf.h"
#include "Organ.h"
#include "organparameter.h"

#include <memory>

namespace CPlantBox {

class Organism;

/**
 * Meshfree analysis of the root system based on signed distance functions.
 */
class SegmentAnalyser
{

public:

    SegmentAnalyser() { }; ///< creates an empty object (use AnalysisSDF::addSegments)
    SegmentAnalyser(std::vector<Vector3d> nodes, std::vector<Vector2i> segments, std::vector<double> segCTs, std::vector<double> radii); ///< everything from scratch
    SegmentAnalyser(const Organism& plant); ///< creates an analyser object containing the segments from the root system
    SegmentAnalyser(const SegmentAnalyser& a) : nodes(a.nodes), segments(a.segments), segCTs(a.segCTs), radii(a.radii), segO(a.segO) { } ///< copy constructor, does not copy user data
    virtual ~SegmentAnalyser() { }; ///< nothing to do here

    // merge segments
    void addSegments(const Organism& plant); ///< adds the segments
    void addSegments(const SegmentAnalyser& a); ///< adds the segments

    // reduce number of segments
    void crop(SignedDistanceFunction* geometry); ///< crops the data to a geometry
    void filter(std::string name, double min, double max); ///< filters the segments to the data @see AnalysisSDF::getParameter
    void filter(std::string name, double value); ///< filters the segments to the data @see AnalysisSDF::getParameter
    void pack(); ///< sorts the nodes and deletes unused nodes

    // some things we might want to know
    std::vector<double> getParameter(std::string name) const; ///< Returns a specific parameter per segment @see RootSystem::ScalarType
    double getSegmentLength(int i) const; ///< returns the length of a segment
    double getSummed(std::string name) const; ///< Sums up the parameter
    double getSummed(std::string name, SignedDistanceFunction* geometry) const; ///< Sums up the parameter within the geometry
    std::vector<double> distribution(std::string name, double top, double bot, int n, bool exact=false) const; ///< vertical distribution of a parameter
    std::vector<SegmentAnalyser> distribution(double top, double bot, int n) const; ///< vertical distribution of a parameter
    std::vector<std::vector<double>> distribution2(std::string name, double top, double bot, double left, double right, int n, int m, bool exact=false) const; ///< 2d distribution (x,z) of a parameter
    std::vector<std::vector<SegmentAnalyser>> distribution2(double top, double bot, double left, double right, int n, int m) const; ///< 2d distribution (x,z) of a parameter
    // todo distribution3

    // rather specialized things we want to know
    std::vector<std::shared_ptr<Organ>> getOrgans() const; ///< segment origins
    int getNumberOfOrgans() const; ///< number of different organs
    SegmentAnalyser foto(const Vector3d& pos, const Matrix3d& ons, double height) const; ///< takes a picture TODO unfinished, untested
    SegmentAnalyser cut(const SDF_HalfPlane& plane) const; ///< returns the segments intersecting with a plane (e.g. for trenches)

    // User data for export or distributions
    void addUserData(std::vector<double> data, std::string name); ///< adds user data that are written into the VTP file, @see SegmentAnalyser::writeVTP
    void clearUserData() { userData.clear(); userDataNames.clear(); } ///< resets the user data

    // some exports
    void write(std::string name, std::vector<std::string>  types = { "radius", "subType", "creationTime", "organType" }); ///< writes simulation results (type is determined from file extension in name)
    void writeVTP(std::ostream & os, std::vector<std::string>  types = { }) const; ///< writes a VTP file
    void writeRBSegments(std::ostream & os) const; ///< Writes the segments of the root system, mimics the Matlab script getSegments()
    void writeDGF(std::ostream & os) const; ///< Writes the segments of the root system in DGF format used by DuMux

    // auxiliary
    static Vector3d cut(Vector3d in, Vector3d out, SignedDistanceFunction* geometry); ///< intersects a line with  the geometry

    std::vector<Vector3d> nodes; ///< nodes
    std::vector<Vector2i> segments; ///< connectivity of the nodes
    std::vector<double> segCTs; ///< creation times of the segments
    std::vector<double> radii; ///< segment radii
    std::vector<std::weak_ptr<Organ>> segO; ///< to look up things

protected:

    std::vector<std::vector<double>> userData; ///< user data attached to the segments (for vtp file), e.g. flux, pressure, etc.
    std::vector<std::string> userDataNames; ///< names of the data added, e.g. "Flux", "Pressure", etc.

};

} // end namespace CPlantBox

#endif
