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

    void setRadius(double a); ///< sets a constant radius for all segments
    void setSubTypes(int t); ///< sets a constant sub type for all segments

    void setSoilGrid(const std::function<int(double,double,double)>& s); ///< sets the soil, resets the mappers, and maps all segments
    void setSoilGrid(const std::function<int(double,double,double)>& s, Vector3d min, Vector3d max, Vector3d res, bool cut = true); ///< sets the soil, resets the mappers, cuts and maps all segments
    void setRectangularGrid(Vector3d min, Vector3d max, Vector3d res, bool cut = true, bool noChanges = false); ///< sets an underlying rectangular grid, and cuts all segments accordingly

    void mapSegments(const std::vector<Vector2i>& segs);
    void cutSegments(); // cut and add segments

    void sort(); ///< sorts segments, each segment belongs to position s.y-1

    std::vector<double> segOuterRadii(int type = 0, const std::vector<double>& vols = std::vector<double>(0)) const; ///< outer cylinder radii to match cell volume
    std::vector<double> segLength() const; ///< calculates segment lengths [cm]
    std::vector<double> getHs(const std::vector<double> sx) const; // return the potential per segment that is given per soil cell
    std::vector<double> getSegmentZ() const; // z-coordinate of segment mid
    std::vector<double> matric2total(const std::vector<double> sx) const;
    std::vector<double> total2matric(const std::vector<double> sx) const;

    int getNumberOfMappedSegments() const { return segments.size(); };  // for the python binding, != getNumberOfSegments (because of shoot roots or cutting)
    std::vector<int> getSegmentMapper() const;  // seg2cell mapper as vector


    std::map<int, int> seg2cell; // root segment to soil cell mapper
    std::map<int, std::vector<int>> cell2seg; // soil cell to root segment mapper

    std::function<int(double,double,double)> soil_index =
        std::bind(&MappedSegments::soil_index_, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3); ///< soil cell index call back function, (care need all MPI ranks in case of dumux)

    std::vector<Vector3d> nodes; ///< nodes [cm]
    std::vector<double> nodeCTs; ///< creation times [days]
    std::vector<Vector2i> segments; ///< connectivity of the nodes
    std::vector<double> radii; ///< radii [cm]
    std::vector<int> subTypes; ///< types [1]
    std::vector<int> organTypes; ///< types of the organ[1]

    Vector3d minBound;
    Vector3d maxBound;
    Vector3d resolution; // cells
    bool cutAtGrid = false;
  bool constantLoc = false;// the roots remain in the soil voxel they appear in


	virtual double getPerimeter(int si_, double l_){return 2 * M_PI * radii[si_];} ///< Perimeter of the segment [cm] overloaded by @see MappedPlant::getPerimeter
	virtual int getSegment2leafId(int si_);

    const double eps = 1.e-5;
    std::array<std::map<int, std::shared_ptr<OrganRandomParameter>>, 5> plantParam;
	double kr_length = -1.0; //define distance to root tipe where kr > 0 as cannot compute distance from age in case of carbon-limited growth
	//% of segment length in the root exchange zone, see MappedPlant::simulate.
	//only needed if carbon- and water-limited growth (i.e., for plants with phloem module)
    std::vector<bool> isRootTip;
	std::vector<double> exchangeZoneCoefs;
  std::vector<double> distanceTip;// save the distance between root segment and root tip (for location-dependent kr)
	std::vector<double> leafBladeSurface; //leaf blade area per segment to define water radial flux. assume no radial flux in petiole
	std::vector<double> segVol; //segment volume <= needed for MappedPlant as leaf does not have cylinder shape necessarally only do segLeaf to have shorter vector?
	std::vector<double> bladeLength;//blade length <= needed for MappedPlant as leaf does not have cylinder shape necessarally only do segLeaf to have shorter vector?
	Vector3d getMinBounds();
		// calcExchangeZoneCoefs() only usefull for carbon-limited growth i.e., with a MappedPlant
	virtual void calcExchangeZoneCoefs(){throw std::runtime_error("calcExchangeZoneCoefs used on MappedSegment instead of MappedPlant object");};
	virtual void calcIsRootTip(){throw std::runtime_error("calcIsRootTip used on MappedSegment instead of MappedPlant object");};

protected:

    void addSegment(Vector2i ns, double radius,  int st,int ot, int i); // adds a single segment at index i, appends the rest if cutted
    void add(Vector2i ns, double radius,  int st, int ot, int i); // adds without cutting, at index i, or appends if i = -1
    double length(const Vector2i& s) const;

    int soil_index_(double x, double y, double z); // default mapper to a equidistant rectangular grid
    void unmapSegments(const std::vector<Vector2i>& segs); ///< remove segments from the mappers

};



/**
 * Build MappedSegmentds sequentially from a RootSystem
 */
class MappedRootSystem : public MappedSegments, public RootSystem
{
public:

    using RootSystem::RootSystem;
	MappedRootSystem(unsigned int seednum = 0): RootSystem(seednum){}; ///< constructor
    virtual ~MappedRootSystem() { }; ///< destructor

    void initialize(bool verbose = true) override { initializeLB(4, 5, verbose); }; ///< overridden, to map initial shoot segments,
    void initializeLB(int basaltype, int shootbornetype, bool verbose = true) override { bool LB = true; initialize_(basaltype, shootbornetype, verbose, LB); }; ///< overridden, to map initial shoot segments,
	void initializeDB(int basaltype, int shootbornetype, bool verbose = true) override { bool LB = false; initialize_(basaltype, shootbornetype, verbose, LB); }; ///< overridden, to map initial shoot segments,

    void simulate(double dt, bool verbose = false) override; ///< build nodes and segments sequentially

    /* segments are shoot and root segments */


    std::shared_ptr<MappedSegments> mappedSegments() { return std::make_shared<MappedSegments>(*this); }  // up-cast for Python binding
    std::shared_ptr<RootSystem> rootSystem() { return std::make_shared<RootSystem>(*this); }; // up-cast for Python binding

protected:
	void initialize_(int basaltype, int shootbornetype, bool verbose = true, bool LB = true);


};

/**
 * Build MappedSegmentds sequentially from a Plant
 */
class MappedPlant : public MappedSegments, public Plant
{
public:

    using Plant::Plant;
	MappedPlant(unsigned int seednum = 0): Plant(seednum){}; ///< constructor
    virtual ~MappedPlant() { }; ///< destructor

    void initializeLB(bool verbose = true, bool stochastic = true) { bool LB = true; initialize_(verbose, stochastic, LB); }; ///< overridden, to map initial shoot segments,
	void initializeDB(bool verbose = true, bool stochastic = true) { bool LB = false; initialize_(verbose, stochastic, LB); }; ///< overridden, to map initial shoot segments,
	void initialize(bool verbose = true, bool stochastic = true) { initializeLB(verbose, stochastic); }; ///< overridden, to map initial shoot segments,
    void simulate(double dt, bool verbose) override ; ///< build nodes and segments sequentially
    void printNodes(); ///< print information
	void mapSubTypes();

    std::shared_ptr<MappedSegments> mappedSegments() { return std::make_shared<MappedSegments>(*this); }  // up-cast for Python binding
    std::shared_ptr<Plant> plant() { return std::make_shared<Plant>(*this); }; // up-cast for Python binding

	std::map<std::tuple<int, int>, int > st2newst; // replace subtypes with other int nummer, so that the N subtypes of one organ type go from 0 to N-1

    virtual double rand() override {if(stochastic){return UD(gen);} else {return 0.5; } }  ///< uniformly distributed random number (0,1)
	virtual double randn() override {if(stochastic){return std::min(std::max(ND(gen),-1.),1.);} else {return 0.5; } }  ///< normally distributed random number (0,1)
	bool stochastic = true;//< whether or not to implement stochasticity, usefull for test files @see test_relative_coordinates.py
	//for photosynthesis and phloem module:
	void calcExchangeZoneCoefs() override;
    void calcIsRootTip() override;
	std::vector<int> getSegmentIds(int ot = -1) const;//needed in phloem module
	std::vector<int> getNodeIds(int ot = -1) const;	//needed in phloem module
	double getPerimeter(int si_, double l_) override; ///< Perimeter of the segment [cm] overloaded by @see MappedPlant::getPerimeter

	int getSegment2leafId(int si_) override; ///< fill segment2Leaf vector
	std::vector<int> segment2leafIds;///< to go from vector of size segment to vectoer of size leaf_segment

 protected:
	void initialize_(bool verbose = true, bool stochastic = true, bool LB = true);
	void getSegment2leafIds(); ///< fill segment2Leaf vector
};

}

#endif

