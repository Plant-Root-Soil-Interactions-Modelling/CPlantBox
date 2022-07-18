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


    void setRadius(double a); ///< sets a constant radius for all segments
    void setSubTypes(int t); ///< sets a constant sub type for all segments

    void setSoilGrid(const std::function<int(double,double,double)>& s); ///< sets the soil, resets the mappers, and maps all segments
    void setSoilGrid(const std::function<int(double,double,double)>& s, Vector3d min, Vector3d max, Vector3d res, bool cut = true); ///< sets the soil, resets the mappers, cuts and maps all segments
    void setRectangularGrid(Vector3d min, Vector3d max, Vector3d res, bool cut = true); ///< sets an underlying rectangular grid, and cuts all segments accordingly

    void mapSegments(const std::vector<Vector2i>& segs);
    void cutSegments(); // cut and add segments

    void sort(); ///< sorts segments, each segment belongs to position s.y-1

    std::vector<double> segOuterRadii(int type = 0, const std::vector<double>& vols = std::vector<double>(0)) const; ///< outer cylinder radii to match cell volume
    std::vector<double> segLength() const; ///< calculates segment lengths [cm]

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

    const double eps = 1.e-5;
    std::array<std::map<int, std::shared_ptr<OrganRandomParameter>>, 5> plantParam;
	double kr_length = -1.0; //define distance to root tipe where kr > 0 as cannot compute distance from age in case of carbon-limited growth
	//% of segment length in the root exchange zone, see MappedPlant::simulate. 
	//only needed if carbon- and water-limited growth (i.e., for plants with phloem module)
	std::vector<double> exchangeZoneCoefs; 
	std::vector<double> leafBladeSurface; //leaf blade area per segment to define water radial flux. assume no radial flux in petiole
	std::vector<double> segVol; //segment volume <= needed for MappedPlant as leaf does not have cylinder shape necessarally only do segLeaf to have shorter vector?
	std::vector<double> bladeLength;//blade length <= needed for MappedPlant as leaf does not have cylinder shape necessarally only do segLeaf to have shorter vector?
	Vector3d getMinBounds();		
		// calcExchangeZoneCoefs() only usefull for carbon-limited growth i.e., with a MappedPlant
	virtual void calcExchangeZoneCoefs(){throw std::runtime_error("calcExchangeZoneCoefs used on MappedSegment instead of MappedPlant object");};		

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

    void initialize(bool verbose = true); ///< overridden, to map initial shoot segments,
    void initializeLB(int basaltype, int shootbornetype, bool verbose = true); ///< overridden, to map initial shoot segments,

    void simulate(double dt, bool verbose = false) override; ///< build nodes and segments sequentially

    /* segments are shoot and root segments */


    std::shared_ptr<MappedSegments> mappedSegments() { return std::make_shared<MappedSegments>(*this); }  // up-cast for Python binding
    std::shared_ptr<RootSystem> rootSystem() { return std::make_shared<RootSystem>(*this); }; // up-cast for Python binding

};

/**
 * Build MappedSegmentds sequentially from a Plant
 */
class MappedPlant : public MappedSegments, public Plant
{
public:

    using Plant::Plant;
	MappedPlant(double seednum = 0): Plant(seednum){}; ///< constructor
    virtual ~MappedPlant() { }; ///< destructor
	
    void initialize(bool verbose = true, bool stochastic = true); ///< overridden, to map initial shoot segments,
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
	std::vector<int> getSegmentIds(int ot = -1) const;//needed in phloem module
	std::vector<int> getNodeIds(int ot = -1) const;	//needed in phloem module			  
};

}

#endif

