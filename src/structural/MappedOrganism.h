// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef MAPPED_ROOTSYSTEM_H_
#define MAPPED_ROOTSYSTEM_H_

#include "RootSystem.h"
#include "Plant.h"

#include <functional>
#include <vector>
#include <tuple>

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
        std::vector<double> radii, std::vector<int> subTypes, std::vector<int> organTypes); ///< for kr and kx age and type dependent
    MappedSegments(std::vector<Vector3d> nodes, std::vector<double> nodeCTs, std::vector<Vector2i> segs,
        std::vector<double> radii, std::vector<int> subTypes); ///< for kr and kx age and type dependent
    MappedSegments(std::vector<Vector3d> nodes, std::vector<Vector2i> segs, std::vector<double> radii); ///< for constant kr, and kx

    virtual ~MappedSegments() { }

	std::map<int,double> sumSegFluxes(const std::vector<double>& segFluxes); ///< sums segment fluxes over soil cells,  soilFluxes = sumSegFluxes(segFluxes), [cm3/day]
    std::vector<double> splitSoilFluxes(const std::vector<double>& soilFluxes, int type = 0) const; ///< splits soil fluxes (per cell) into segment fluxes


    void setSoilGrid(const std::function<int(double,double,double)>& s, bool noChanges = false); ///< sets the soil, resets the mappers, and maps all segments
    void setRectangularGrid(Vector3d min, Vector3d max, Vector3d res, bool cut = true, bool noChanges = false); ///< sets an underlying rectangular grid, and cuts all segments accordingly
    void mapSegments(const std::vector<Vector2i>& segs);
    void cutSegments(); // cut and add segments
    int getNumberOfMappedSegments() const { return segments.size(); };  // for the python binding, != getNumberOfSegments (because of shoot roots or cutting)
    int getNumberOfMappedNodes() const { return nodes.size(); }; // might contain extra nodes
    std::vector<int> getSegmentMapper() const;  // seg2cell mapper as vector

    void sort(); ///< sorts segments, each segment belongs to position s.y-1

    std::vector<double> segLength() const; ///< calculates segment lengths [cm]
    std::vector<double> getHs(const std::vector<double> sx) const; // return the potential per segment that is given per soil cell
    std::vector<double> getSegmentZ() const; // z-coordinate of segment mid
    std::vector<double> matric2total(const std::vector<double> sx) const;
    std::vector<double> total2matric(const std::vector<double> sx) const;
    virtual std::vector<double> getEffectiveRadii();


    virtual double getEffectiveRadius(int si) { return this->radii.at(si); };
    virtual double getPerimeter(int si_, double l_) { return 2 * M_PI * radii[si_];} ///< Perimeter of the segment [cm] overloaded by @see MappedPlant::getPerimeter
    virtual int getSegment2leafId(int si_);
    virtual void calcExchangeZoneCoefs() { throw std::runtime_error("calcExchangeZoneCoefs used on MappedSegment instead of MappedPlant object"); }; // calcExchangeZoneCoefs() only usefull for carbon-limited growth i.e., with a MappedPlant

    Vector3d getMinBounds();
    void setRadius(double a); ///< sets a constant radius for all segments
    void setSubTypes(int t); ///< sets a constant sub type for all segments

    // nodes
    std::vector<Vector3d> nodes; ///< nodes [cm]
    std::vector<double> nodeCTs; ///< creation times [days]

    // segments
    std::vector<Vector2i> segments; ///< connectivity of the nodes
    std::vector<double> radii; ///< radii [cm]
    std::vector<int> subTypes; ///< types [1]
    std::vector<int> organTypes; ///< types of the organ[1]
    std::vector<double> leafBladeSurface; //leaf blade area per segment (to define water radial flux assume no radial flux in petiole)
    std::vector<double> segVol; //segment volume <= needed for MappedPlant as leaf does not have cylinder shape necessarally only do segLeaf to have shorter vector?
    std::vector<double> bladeLength;//blade length <= needed for MappedPlant as leaf does not have cylinder shape necessarally only do segLeaf to have shorter vector?
    std::vector<std::weak_ptr<Organ>> segO; ///< for SegmentAnalyser

    // calculated by calcExchangeZoneCoefs()
    std::vector<double> distanceTip;// save the distance between root segment and root tip (for location-dependent kr)
    std::vector<double> exchangeZoneCoefs;

    Vector3d minBound; // grid bounds
    Vector3d maxBound;
    Vector3d resolution; // cells
    bool cutAtGrid = false;
    bool constantLoc = false;// the roots remain in the soil voxel they appear in

    std::map<int, int> seg2cell; // root segment to soil cell mapper
    std::map<int, std::vector<int>> cell2seg; // soil cell to root segment mapper

    std::function<int(double,double,double)> soil_index =
        std::bind(&MappedSegments::soil_index_, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3); ///< soil cell index call back function, (care need all MPI ranks in case of dumux)

    const double eps = 1.e-5; // for cutting

	double kr_length = -1.0; //define distance to root tipe where kr > 0 as cannot compute distance from age in case of carbon-limited growth
    //% of segment length in the root exchange zone, see MappedPlant::simulate.
    //only needed if carbon- and water-limited growth (i.e., for plants with phloem module)


protected:

    void addCutSegment(Vector2i ns, double a, int st, int ot, double lbsurf, double vols, double blen, std::weak_ptr<Organ> segO, int i); // adds a single segment at index i, appends the rest if cutted
    void addSegment(Vector2i ns, double a,  int st, int ot, double lbsurf, double vols, double blen, std::weak_ptr<Organ> segO, int i); // adds without cutting, at index i, or appends if i = -1
    double length(const Vector2i& s) const;

    int soil_index_(double x, double y, double z); // default mapper to a equidistant rectangular grid
    void unmapSegments(const std::vector<Vector2i>& segs); ///< remove segments from the mappers

};


/**
 * Build MappedSegments sequentially from a Plant
 */
class MappedPlant : public MappedSegments, public Plant
{
public:

	MappedPlant(unsigned int randomSeed = 0): MappedSegments(), Plant(randomSeed) { }; ///< constructor
    virtual ~MappedPlant() { }; ///< destructor

    std::shared_ptr<MappedSegments> mappedSegments() { return std::make_shared<MappedSegments>(*this); }  // up-cast for Python binding
    std::shared_ptr<Plant> plant() { return std::make_shared<Plant>(*this); }; // up-cast for Python binding

    void disableExtraNode() { extraNode = 0; }
    void enableExtraNode() { extraNode = 1; } // extra collar node for collar BC easier

    // I made all initializer functions virtual and having only the verbose as argument to avoid confusion, stochasity can be set by Organism::setStochastic
    void initializeLB(bool verbose = true) override { initialize_(verbose,  true); }; ///< overridden, length based initialization
	void initializeDB(bool verbose = true) override { initialize_(verbose,  false); }; ///< overridden, delay based based initialization
	void initialize(bool verbose = true) override { initializeLB(verbose); }; ///< overridden, to map initial nodes, segs
	void simulate(double dt, bool verbose) override ; ///< build nodes and segments sequentially

    void printNodes(); ///< print information
	void mapSubTypes();

	//for photosynthesis and phloem module:
	void calcExchangeZoneCoefs() override;
	std::vector<int> getSegmentIds(int ot = -1) const;//needed in phloem module
	std::vector<int> getNodeIds(int ot = -1) const;	//needed in phloem module
	double getEffectiveRadius(int si) override;
	double getPerimeter(int si_, double l_) override; ///< Perimeter of the segment [cm] overloaded by @see MappedPlant::getPerimeter
	int getSegment2leafId(int si_) override; ///< fill segment2Leaf vector

//	bool stochastic = true;//< whether or not to implement stochasticity, usefull for test files @see test_relative_coordinates.py
//	virtual double rand() override {if(stochastic){return UD(gen);} else {return 0.5; } }  ///< uniformly distributed random number (0,1)
//	virtual double randn() override {if(stochastic){return std::min(std::max(ND(gen),-1.),1.);} else {return 0.5; } }  ///< normally distributed random number (0,1)
//	// TODO randn() is cut, in Plant it is not, this is confusing

	std::map<std::tuple<int, int>, int > st2newst; // replace subtypes with other int nummer, so that the N subtypes of one organ type go from 0 to N-1
    std::vector<int> segment2leafIds;///< to go from vector of size segment to vector of size leaf_segment


 protected:

    int extraNode = -1; // -1 .. choose automatic, 0.. False, 1.. True

	bool rootHairs = true; // todo: determine from parameters, and set within constructor

	void initialize_(bool verbose = true, bool lengthBased = true);
	void getSegment2leafIds(); ///< fill segment2Leaf vector

};

}

#endif

