// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef SEGMENTANALYSER_H_
#define SEGMENTANALYSER_H_

#include "sdf.h"
#include "Organ.h"
#include "organparameter.h"

#include <memory>
#include <limits>
#include <tuple>

namespace CPlantBox {

class Organism;
class MappedSegments;
class XylemFlux;
class PlantHydraulicParameters;
class PlantHydraulicModel;

/**
 * Meshfree analysis of the root system based on signed distance functions.
 */
class SegmentAnalyser
{

public:

    SegmentAnalyser() { }; ///< creates an empty object (use AnalysisSDF::addSegments)
    SegmentAnalyser(const std::vector<Vector3d>& nodes, const std::vector<Vector2i>& segments,
    		const std::vector<double>& segCTs, const std::vector<double>& radii); ///< everything from scratch
    SegmentAnalyser(const Organism& plant); ///< creates an analyser object containing the segments from the organism @param plant
    SegmentAnalyser(const MappedSegments& plant); ///< creates an analyser object containing the segments from @param plant
    virtual ~SegmentAnalyser() { }; ///< nothing to do here

    // merge segments
    void addSegments(const Organism& plant); ///< adds the segments
    void addSegments(const SegmentAnalyser& a); ///< adds the segments
    void addSegment(Vector2i seg, double ct, double radius, bool insert = false); ///< adds a single segment

    // to user data to later visualize results
    void addAge(double simtime);  // "age"
    void addConductivities(const XylemFlux& xylem, double simtime, double kr_max = 1.e6, double kx_max = 1.e6); // "kr", "kx"
    void addHydraulicConductivities(const PlantHydraulicParameters& xylem, double simtime, double kr_max = 1.e6, double kx_max = 1.e6); // "kr", "kx"
    void addFluxes(const PlantHydraulicModel& rs, const std::vector<double>& rx, const std::vector<double>& sx, double simTime); // "axial_flux", "radial_flux"
    void addCellIds(const MappedSegments& plant); // "cell_id"

    // reduce number of segments
    void crop(std::shared_ptr<SignedDistanceFunction> geometry); ///< crops the data to a geometry
    void cropDomain(double xx, double yy, double zz); // crops to the domain @see SDF_PlantBox
    void filter(std::string name, double min, double max); ///< filters the segments to the data @see AnalysisSDF::getParameter
    void filter(std::string name, double value); ///< filters the segments to the data @see AnalysisSDF::getParameter
    void pack(); ///< sorts the nodes and deletes unused nodes
    Vector3d getMinBounds(); ///< get minimum of node coordinates (e.g. minimum corner of bounding box)
    Vector3d getMaxBounds(); ///< get maximum of node coordinates (e.g. maximum corner of bounding box)

    // some things we might want to know
    std::vector<double> getParameter(std::string name, double def = std::numeric_limits<double>::quiet_NaN()) const; ///< Returns a specific parameter per segment @see RootSystem::ScalarType
    double getSegmentLength(int i) const; ///< returns the length of a segment
    double getSummed(std::string name) const; ///< Sums up the parameter
    double getSummed(std::string name, std::shared_ptr<SignedDistanceFunction> geometry) const; ///< Sums up the parameter within the geometry (e.g. for length or surface)
    std::vector<double> distribution(std::string name, double top, double bot, int n, bool exact=false) const; ///< vertical distribution of a parameter
    std::vector<SegmentAnalyser> distribution(double top, double bot, int n) const; ///< vertical distribution
    std::vector<std::vector<double>> distribution2(std::string name, double top, double bot, double left, double right, int n, int m, bool exact=false) const; ///< 2d distribution (x,z) of a parameter
    std::vector<std::vector<SegmentAnalyser>> distribution2(double top, double bot, double left, double right, int n, int m) const; ///< 2d distribution (x,z)
    // todo distribution3

    // rather specialized things we want to know
    void mapPeriodic(double xx, double yy); /// maps into a periodic domain, splits up intersecting segments
    void map2D(); ///< maps the 3d coordinates to the x-z plan (sqrt(x2+y2), 0., z)
    std::vector<std::shared_ptr<Organ>> getOrgans(int ot = -1) const; ///< segment origins
    int getNumberOfOrgans() const; ///< number of different organs
    SegmentAnalyser foto(const Vector3d& pos, const Matrix3d& ons, double height) const; ///< takes a picture TODO unfinished, untested
    SegmentAnalyser cut(const SignedDistanceFunction& plane) const; ///< returns the segments intersecting with a plane (e.g. for trenches)


    // User data for export or distributions
    void addData(std::string name, std::vector<double> data); ///< adds user data that are written into the VTP file, @see SegmentAnalyser::writeVTP

    // some exports
    void write(std::string name, std::vector<std::string>  types = { "radius", "subType", "creationTime", "organType" }); ///< writes simulation results (type is determined from file extension in name)
    void writeVTP(std::ostream & os, std::vector<std::string>  types = { "radius", "subType", "creationTime", "organType"  }) const; ///< writes a VTP file
    void writeRBSegments(std::ostream & os) const; ///< Writes the segments of the root system, mimics the Matlab script getSegments()
    void writeDGF(std::ostream & os) const; ///< Writes the segments of the root system in DGF format used by DuMux

    // auxiliary
    static Vector3d cut(Vector3d in, Vector3d out, const std::shared_ptr<SignedDistanceFunction>& geometry, double eps = 1.e-6); ///< intersects a line with  the geometry

    std::vector<Vector3d> nodes; ///< nodes
    std::vector<Vector2i> segments; ///< connectivity of the nodes
    std::vector<std::weak_ptr<Organ>> segO; ///< to look up things
    std::map<std::string, std::vector<double>> data; ///< user data attached to the segments (for vtp file), e.g. flux, pressure, etc.	

protected:

    void mapPeriodic_(double xx, Vector3d axis, double eps);

};

} // end namespace CPlantBox

#endif
