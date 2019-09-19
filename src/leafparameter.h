// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef LEAFPARAMETER_H_
#define LEAFPARAMETER_H_

#include "mymath.h"
#include "soil.h"
#include "organparameter.h"
#include "tropism.h"
#include "growth.h"

#include <vector>

namespace CPlantBox {

/**
 * Parameters of a single leaf (created by LeafSpecificParameter)
 */
class LeafSpecificParameter : public OrganSpecificParameter
{
public:

	LeafSpecificParameter() :OrganSpecificParameter() { };
	LeafSpecificParameter(int type, double lb, double la, const std::vector<double>& ln, double r, double a, double theta, double rlt):
		OrganSpecificParameter() , lb(lb), la(la), r(r), a(a), theta(theta), rlt(rlt), ln(ln){  };

	/*
	 * Parameters per leaf
	 */
	double lb = 0.; 		///< Basal zone of leaf vein [cm]
	double la = 0.;			///< Apical zone of leaf vein [cm];
	double r = 0.;			///< Initial growth rate [cm day-1]
	double a = 0.; 			///< Leaf width [cm]
	double theta = 0.; 		///< Branching angle between veins [rad]
	double rlt = 0.;		///< Leaf life time [day]
	std::vector<double> ln = std::vector<double>();    ///< Inter-lateral distances [cm]

	double getK() const; ///< Returns the exact maximal leaf length of this realization [cm]

	std::string toString() const override;

};



/**
 * A parameter set describing a leaf type
 */
class LeafRandomParameter : public OrganRandomParameter
{
public:

	LeafRandomParameter(Organism* plant); ///< default constructor
	virtual ~LeafRandomParameter();

    OrganRandomParameter* copy(Organism* plant_) override;


	OrganSpecificParameter* realize() override; ///< Creates a specific leaf from the leaf parameter set
	int getLateralType(const Vector3d& pos); ///< Choose (dice) lateral type based on leaf parameter set
	double getK() const { return std::max(nob-1,double(0))*ln+la+lb; }  ///< returns the mean maximal leaf length [cm]

	std::string toString(bool verbose = true) const override; ///< writes parameter to a string

    void readXML(tinyxml2::XMLElement* element) override; ///< reads a single sub type organ parameter set
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
	double k = 10.;			///< Maximal leaf length [cm]
	double ks =0.;			///< Standard deviation maxial leaf length [cm]
	double nob = 0.; 		///< Number of branches [1]
	double nobs = 0.; 		///<  TODO get rid of nobs
	double r = 1.;			///< Initial growth rate [cm day-1]
	double rs = 0.;			///< Standard deviation initial growth rate [cm day-1]
	double a = 0.1; 		///< Leaf width [cm]
	double as = 0.; 		///< Standard deviation leaf width [cm]
	double RotBeta = 0.6;	///< Radial rotation (roll) (rad)
	double BetaDev = 0.2;	///< Deviation of radial rotation (rad)
	double InitBeta = 0.2;	///< Initial radial rotation (rad)
	int tropismT = 1;		///< Leaf tropism parameter (Type)
	double tropismN = 1.;	///< Leaf tropism parameter (number of trials)
	double tropismS = 0.2;	///< Leaf tropism parameter (mean value of expected changeg) [1/cm]
	double dx = 0.25; 		///< Maximal segment size [cm]
	double theta = 1.22;	///< Angle between leafvein and parent leafvein (rad)
	double thetas = 0.; 	///< Standard deviation angle between leafvein and parent leafvein (rad)
	double rlt = 1.e9;		///< Leaf life time (days)
	double rlts = 0.;		///< Standard deviation of leaf life time (days)
	int gf = 1;				///< Growth function (1=negative exponential, 2=linear)
	std::vector<int> successor;			///< Lateral types [1]
	std::vector<double> successorP; 	///< Probabiltities of lateral type to emerge (sum of values == 1) [1]

	/* call back functions */
    Tropism* f_tf;  ///< tropism function ( = new Tropism(plant) )
    GrowthFunction* f_gf = new ExponentialGrowth(); ///< growth function
    SoilLookUp* f_se = new SoilLookUp(); ///< scale elongation function
    SoilLookUp* f_sa = new SoilLookUp(); ///< scale angle function
    SoilLookUp* f_sbp = new SoilLookUp(); ///< scale branching probability functiongrowth

protected:

    void bindParameters(); ///<sets up class introspectionbindParameters

};

} // end namespace CPlantBox

#endif
