// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef LEAFPARAMETER_H_
#define LEAFPARAMETER_H_

#include "mymath.h"
#include "soil.h"
#include "organparameter.h"
#include "growth.h"

#include <vector>

namespace CPlantBox {

class Organism;

/**
 * Parameters of a single leaf (created by LeafSpecificParameter)
 */
class LeafSpecificParameter : public OrganSpecificParameter
{
public:

	LeafSpecificParameter() :OrganSpecificParameter(-1, 0.) { };
	LeafSpecificParameter(int subType, double lb, double la, 
	const std::vector<double>& ln, double r, double a, double theta, 
	double rlt, double leafArea, bool laterals, double Width_blade, double Width_petiole):
		OrganSpecificParameter(subType, a) , lb(lb), la(la), r(r), 
		theta(theta), rlt(rlt), areaMax(leafArea), laterals(laterals), 
		ln(ln), Width_blade(Width_blade), Width_petiole(Width_petiole)  { }; ///< Constructor setting all parameters

	/*
	 * Parameters per leaf
	 */
	double lb = 0.; 		///< Basal zone of leaf (leaf-stem) [cm]
	double la = 0.;			///< Apical zone of leaf vein [cm];
	double r = 0.;			///< Initial growth rate [cm day-1]
	double theta = 0.; 		///< Branching angle between veins [rad]
	double rlt = 0.;		///< Leaf life time [day]
	double areaMax = 0.; 	///< Leaf area [cm2]
	bool laterals = false;  ///< Indicates if lateral leaves exist
	std::vector<double> ln = std::vector<double>(); ///< Inter-lateral distances (if laterals) or mid for radial parametrisation (if there are no laterals) [cm]
	double Width_blade = 0.;		///< width of leafe blade (cm) = length - lb zone. define later a width growth rate?
	double Width_petiole = 0.;		///< width of leafe petiole (cm) = lb zone. define later a width growth rate?
	int nob() const { return ln.size() + laterals; } //number of laterals = number of phytomers + 1
	double getK() const; ///< Returns the exact maximal leaf length (including leaf stem) of this realization [cm]
	double leafLength() const { return getK()-lb; }; ///< Returns the exact maximal leaf length (excluding leaf stem) of this realization [cm]

	std::string toString() const override;

};



/**
 * A parameter set describing a leaf type
 */
class LeafRandomParameter : public OrganRandomParameter
{
public:
	enum shapeTypes { shape_cylinder = 0, shape_cuboid = 1, shape_2D = 2}; ///< how is the shape of the leaf defined?, see @orgVolume and @orgVolume2Length
    LeafRandomParameter(std::shared_ptr<Organism> plant); ///< default constructor
	virtual ~LeafRandomParameter() { };

    void createLeafGeometry(std::vector<double> y, std::vector<double> l, int N); // create normalized leaf geometry
    void createLeafRadialGeometry(std::vector<double> phi, std::vector<double> l, int N); // create normalized leaf geometry from a radial parameterization

	std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant) override;

	std::shared_ptr<OrganSpecificParameter> realize() override; ///< Creates a specific leaf from the leaf parameter set

	double nob() const { return std::max((lmax-la-lb)/ln+1, 1.); }  ///< returns the mean number of branching points [1]
    double nobs() const; ///< returns the standard deviation of number of branches [1]
    double leafLength() { return lmax-lb; }; // lb represents the leaf base
    double leafMid() { return lmax-la-lb; }; //

	std::string toString(bool verbose = true) const override; ///< writes parameter to a string

    void readXML(tinyxml2::XMLElement* element, bool verbose) override; ///< reads a single sub type organ parameter set
    tinyxml2::XMLElement* writeXML(tinyxml2::XMLDocument& doc, bool comments = true) const override; ///< writes a organ leaf parameter set

	/*
	 * Parameters per leaf type
	 */
	double lb = 0.; 	///< Basal zone [cm]
	double lbs = 0.;  	///< Standard deviation basal zone [cm]
	double la = 10.;	///< Apical zone [cm];
	double las = 0.;	///< Standard deviation apical zone [cm];
	double ln = 1.; 	///< Inter-lateral distance [cm]
	double lns = 0.;  	///< Standard deviation inter-lateral distance [cm]
	int lnf = 0; 		///< type of inter-branching distance (0 homogeneous, 1 linear inc, 2 linear dec, 3 exp inc, 4 exp dec)
    double lmax = 0.;       ///< Maximal stem length [cm]
    double lmaxs = 0.;      ///< Standard deviation of maximal stem length [cm]
    double areaMax = 10.; 	///< maximal leaf area (reached when stem length reaches lmax) [cm2]
    double areaMaxs = 0.; 	///< Standard deviation of maximal leaf area [cm2]
    double r = 1.;			///< Initial growth rate [cm day-1]
	double rs = 0.;			///< Standard deviation initial growth rate [cm day-1]
	double rotBeta = 0.6;	///< Radial rotation (roll) (rad)
	double betaDev = 0.2;	///< Deviation of radial rotation (rad)
	double initBeta = 0.2;	///< Initial radial rotation (rad)
	int tropismT = 1;		///< Leaf tropism parameter (Type)
	double tropismN = 1.;	///< Leaf tropism parameter (number of trials)
	double tropismS = 0.2;	///< Leaf tropism parameter (mean value of expected changeg) [1/cm]
	double tropismAge = 0.;	///< Leaf tropism parameter (age when switch tropism)
	double tropismAges = 0.;///< Leaf tropism parameter (age when switch tropism, standard deviation)
	double theta = 1.22;	///< Angle between leafvein and parent leafvein (rad)
	double thetas = 0.; 	///< Standard deviation angle between leafvein and parent leafvein (rad)
	double rlt = 1.e9;		///< Leaf life time (days)
	double rlts = 0.;		///< Standard deviation of leaf life time (days)
	double Width_blade = 0.;		///< width of leafe blade (cm) = length - lb zone. define later a width growth rate?
	double Width_blades = 0.;		///< Standard deviation of leaf blade width (cm)
	double Width_petiole = 0.;		///< width of leafe petiole (cm) = lb zone. define later a width growth rate?
	double Width_petioles = 0.;		///< Standard deviation of leaf petiole width (cm)
	int gf = 1;				///< Growth function (1=negative exponential, 2=linear)
	int isPseudostem = 0;				///< do the leaf sheaths make a pseudostem? (0 false, 1 true)
	std::vector<int> successor = {};			///< Lateral types [1]
	std::vector<double> successorP = {}; 	///< Probabiltities of lateral type to emerge (sum of values == 1) [1]

	/* describes the plant geometry */
	std::vector<double> leafGeometryPhi= {}; //2D shape
	std::vector<double> leafGeometryX= {};//2D shape
	int parametrisationType = 0; // 2D shape type : 0 .. radial, 1..along main axis
	//how is the shape of the leaf deined? cylinder (a = radius), cuboid (a = thickness, Width_blade, Width_petiole), 2D (leafGeometryPhi, leafGeometryX, areaMax)
	int shapeType = 2; //default: 2D shape
	/* call back functions */
    std::shared_ptr<SoilLookUp> f_se = std::make_shared<SoilLookUp>(); ///< scale elongation function
    std::shared_ptr<SoilLookUp> f_sa = std::make_shared<SoilLookUp>(); ///< scale angle function
    std::shared_ptr<SoilLookUp> f_sbp = std::make_shared<SoilLookUp>(); ///< scale branching probability function

    std::vector<std::vector<double>> leafGeometry; // normalized x - coordinates per along the normalized mid vein
	int geometryN = 100; // leaf geometry resolution (not in XML)

protected:

	
    void bindParameters() override; ///<sets up class introspectionbindParameters
    std::vector<double> intersections(double y, std::vector<double> phi, std::vector<double> l); ///< returns the intersection of a horizontal line at y-coordinate with the leaf geometry
    void normalizeLeafNodes(); ///< scales leaf area to 1

};

} // end namespace CPlantBox

#endif
