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
class PlantHydraulicParameters;
class PlantHydraulicModel;

/**
 * @brief Mesh-free analysis of plant segments based on signed distance functions.
 *
 * Stores nodes, segments, and named per-segment data (e.g. radius, creationTime).
 * Provides methods to add, filter, crop, and export segments, as well as to
 * compute spatial distributions and hydraulic quantities.
 */
class SegmentAnalyser
{

public:

    SegmentAnalyser() { }; ///< Creates an empty object; populate with addSegments().
    SegmentAnalyser(const std::vector<Vector3d>& nodes, const std::vector<Vector2i>& segments,
    		const std::vector<double>& segCTs, const std::vector<double>& radii); ///< Constructs from raw node, segment, creation-time, and radius vectors.
    SegmentAnalyser(const Organism& plant); ///< Constructs from an Organism, copying its segments and per-segment parameters.
    SegmentAnalyser(const MappedSegments& plant); ///< Constructs from a MappedSegments object.
    virtual ~SegmentAnalyser() { };

    /** @name Add segments */
    ///@{
    void addSegments(const Organism& plant); ///< Appends all segments from @p plant.
    void addSegments(const SegmentAnalyser& a); ///< Appends all segments from analyser @p a.
    void addSegment(Vector2i seg, double ct, double radius, bool insert = false); ///< Appends (or prepends if @p insert) a single segment.
    ///@}

    /** @name Add per-segment data */
    ///@{
    void addAge(double simtime);  ///< Computes segment age (= @p simtime - creationTime) and stores it as "age".
    void addHydraulicConductivities(const PlantHydraulicParameters& xylem, double simtime, double kr_max = 1.e6, double kx_max = 1.e6); ///< Evaluates and stores radial ("kr") and axial ("kx") conductivities.
    void addFluxes(const PlantHydraulicModel& rs, const std::vector<double>& rx, const std::vector<double>& sx, double simTime); ///< Computes and stores "axial_flux" [cm³/day] and "radial_flux" [cm³/cm²/day].
    void addCellIds(const MappedSegments& plant); ///< Stores the soil cell index for each segment as "cell_id".
    ///@}

    /** @name Filter and reduce */
    ///@{
    void crop(std::shared_ptr<SignedDistanceFunction> geometry); ///< Crops segments to those inside @p geometry, splitting intersecting segments exactly.
    void cropDomain(double xx, double yy, double zz); ///< Crops to the box [-xx/2, xx/2] x [-yy/2, yy/2] x [-zz, 0]. @see SDF_PlantBox
    void filter(std::string name, double min, double max); ///< Keeps only segments where parameter @p name is in [@p min, @p max].
    void filter(std::string name, double value); ///< Keeps only segments where parameter @p name equals @p value.
    void pack(); ///< Re-indexes nodes and removes unused ones (call after crop/filter to save memory).
    Vector3d getMinBounds(); ///< Returns the component-wise minimum of all node coordinates.
    Vector3d getMaxBounds(); ///< Returns the component-wise maximum of all node coordinates.
    ///@}

    /** @name Query */
    ///@{
    std::vector<double> getParameter(std::string name, double def = std::numeric_limits<double>::quiet_NaN()) const; ///< Returns a named parameter per segment; falls back to @p def (or throws) for expired organ pointers.
    double getSegmentLength(int i) const; ///< Returns the Euclidean length of segment @p i.
    double getSummed(std::string name) const; ///< Returns the sum of parameter @p name over all segments.
    double getSummed(std::string name, std::shared_ptr<SignedDistanceFunction> geometry) const; ///< Returns the approximate sum of @p name for segments whose midpoint lies inside @p geometry.
    std::vector<double> distribution(std::string name, double top, double bot, int n, bool exact=false) const; ///< Returns a vertical distribution of @p name in @p n layers between @p bot and @p top.
    std::vector<SegmentAnalyser> distribution(double top, double bot, int n) const; ///< Returns @p n SegmentAnalyser objects for vertical layers between @p bot and @p top.
    std::vector<std::vector<double>> distribution2(std::string name, double top, double bot, double left, double right, int n, int m, bool exact=false) const; ///< Returns a 2-D (x–z) distribution of @p name on an @p n × @p m grid.
    std::vector<std::vector<SegmentAnalyser>> distribution2(double top, double bot, double left, double right, int n, int m) const; ///< Returns a 2-D (x–z) grid of SegmentAnalyser objects.
    // todo distribution3
    ///@}

    /** @name Geometry transforms */
    ///@{
    void mapPeriodic(double xx, double yy); ///< Maps segments into the periodic domain [-xx/2, xx/2) x [-yy/2, yy/2), splitting boundary-crossing segments.
    void map2D(); ///< Projects nodes to the x–z plane: (x,y,z) → (√(x²+y²), 0, z).
    std::vector<std::shared_ptr<Organ>> getOrgans(int ot = -1) const; ///< Returns the unique organs that own the segments (optionally filtered by organ type @p ot).
    int getNumberOfOrgans() const; ///< Returns the number of distinct organs in the segment set.
    SegmentAnalyser foto(const Vector3d& pos, const Matrix3d& ons, double height) const; ///< Projects segments onto an image plane (unfinished/untested).
    SegmentAnalyser cut(const SignedDistanceFunction& plane) const; ///< Returns the segments that intersect @p plane (e.g. for trench sampling).
    ///@}

    /** @name User data */
    ///@{
    void addData(std::string name, std::vector<double> data); ///< Attaches named per-segment (or per-node) data for use in getParameter() and writeVTP(). @see writeVTP
    ///@}

    /** @name Export */
    ///@{
    void write(std::string name, std::vector<std::string>  types = { "radius", "subType", "creationTime", "organType" }); ///< Writes to file; format determined by extension (.vtp, .txt, .dgf).
    void writeVTP(std::ostream & os, std::vector<std::string>  types = { "radius", "subType", "creationTime", "organType"  }) const; ///< Writes a VTK PolyData (.vtp) file with the specified cell-data fields.
    void writeRBSegments(std::ostream & os) const; ///< Writes segments in RootBox format (mimics the Matlab getSegments() script).
    void writeDGF(std::ostream & os) const; ///< Writes segments in DGF format for use with DuMux.
    ///@}

    /** @name Auxiliary */
    ///@{
    static Vector3d cut(Vector3d in, Vector3d out, const std::shared_ptr<SignedDistanceFunction>& geometry, double eps = 1.e-6); ///< Numerically bisects a segment to find its intersection with @p geometry.
    ///@}

    std::vector<Vector3d> nodes; ///< Node coordinates.
    std::vector<Vector2i> segments; ///< Segment connectivity (pairs of node indices).
    std::vector<std::weak_ptr<Organ>> segO; ///< Weak pointers to the organ each segment belongs to.
    std::map<std::string, std::vector<double>> data; ///< Named per-segment scalar data (e.g. "radius", "creationTime", "flux").

protected:

    void mapPeriodic_(double xx, Vector3d axis, double eps); ///< Helper: maps segments into a 1-D periodic domain along @p axis.

};

} // end namespace CPlantBox

#endif
