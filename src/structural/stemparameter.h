// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef STEMPARAMETER_H_
#define STEMPARAMETER_H_

#include "mymath.h"
#include "soil.h"
#include "growth.h"
#include "organparameter.h"

/**
 * This file describes the classes StemSpecificParameter and StemRandomParameter.
 * StemSpecificParameter are drawn from the StemRandomParameter class
 */

namespace CPlantBox {

class Organism;

/**
 * Parameters of a specific stem, its created by StemRandomParameter:realize()
 */
class StemSpecificParameter :public OrganSpecificParameter
{

public:

    StemSpecificParameter(): StemSpecificParameter(-1,0.,0.,std::vector<double>(0),0,0.,0.,0.,0.) { } ///< Default constructor
    StemSpecificParameter(int type, double lb, double la, const std::vector<double>& ln, double r, double a, double theta, double rlt,
	bool laterals= false, int nodalGrowth = 0, double delayNGStart = 0.,double delayNGEnd = 0., double delayLat = 0.):
        OrganSpecificParameter(type, a),  lb(lb), la(la), r(r), theta(theta), rlt(rlt), ln(ln),
		laterals(laterals), nodalGrowth(nodalGrowth), delayNGStart(delayNGStart), delayNGEnd(delayNGEnd), delayLat(delayLat) { } ///< Constructor setting all parameters

    /*
     * Stem parameters per single stem
     */
    double lb;              ///< Basal zone [cm]
    double la;              ///< Apical zone [cm]
    double r;               ///< Initial growth rate [cm day-1]
    double theta;           ///< Angle between stem and parent stem [rad]
    double rlt;             ///< Stem life time [day]
    std::vector<double> ln; ///< Inter-lateral distances [cm]
	bool laterals = false;
	int nodalGrowth;			///< whether to implement the internodal growth [1] (see @stem::simulate)
	double delayNGStart;
	double delayNGEnd;
	double delayLat;
    int nob() const { return ln.size() + laterals; } ///< return the maximal number of branching points [1]
    double getK() const; ///< Returns the exact maximal stem length of this realization [cm]

    std::string toString() const override; ///< for debugging

};



/**
 * Contains a parameter set describing a stem type
 */
class StemRandomParameter :public OrganRandomParameter
{

public:

    StemRandomParameter(std::shared_ptr<Organism> plant); ///< default constructor
    virtual ~StemRandomParameter() { };

    std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant) override;

    std::shared_ptr<OrganSpecificParameter> realize() override; ///< Creates a specific stem from the stem parameter set

    double nob() const { return std::max((lmax-la-lb)/ln +1, 1.); }  ///< returns the mean number of branches [1]
    double nobs() const; ///< returns the standard deviation of number of branches [1]


    /*
     * Parameters per stem type
     */
    double lb = 0.; 	    ///< Basal zone [cm]
    double lbs = 0.;        ///< Standard deviation basal zone [cm]
    double la = 10.;	    ///< Apical zone [cm];
    double las = 0.;    	///< Standard deviation apical zone [cm];
    double ln = 1; 		    ///< Inter-lateral distance [cm]
    double lns = 0.;    	///< Standard deviation inter-lateral distance [cm]
    int lnf = 0;            ///< type of inter-branching distance (0 homogeneous, 1 linear inc, 2 linear dec, 3 exp inc, 4 exp dec)
    double lmax = 0.;       ///< Maximal stem length [cm]
    double lmaxs = 0.;      ///< Standard deviation of maximal stem length [cm]
    double r = 1;		    ///< Initial growth rate [cm day-1]
    double rs = 0.;	    	///< Standard deviation initial growth rate [cm day-1]
    double a = 0.1; 		///< Stem radius [cm]
    double as = 0.; 		///< Standard deviation stem radius [cm]
    double k = 0.;          ///< Maximal stem length [cm]
    double ks = 0.;         ///< Maximal stem length deviation [cm]
    double rotBeta = 0.6;	///< Revrotation
    double betaDev = 0.2;	///< Deviation of RevRotation
    double initBeta = 0.2;	///< Initial RevRotation
    int tropismT = 1;	    ///< Stem tropism parameter (Type)
    double tropismN = 1.;   ///< Stem tropism parameter (number of trials)
    double tropismS = 0.2;  ///< Stem tropism parameter (mean value of expected changeg) [1/cm]
	double tropismAge = 0.;	///< Leaf tropism parameter (age when switch tropism)
	double tropismAges = 0.;///< Leaf tropism parameter (age when switch tropism, standard deviation)
    double theta = 1.22; 	///< Angle between stem and parent stem (rad)
    double thetas= 0.; 	    ///< Standard deviation angle between stem and parent stem (rad)
    double rlt = 1e9;		///< Stem life time (days)
    double rlts = 0.;	    ///< Standard deviation stem life time (days)
    int gf = 1;			    ///< Growth function (1=negative exponential, 2=linear)
	int nodalGrowth = 0;		///< whether to implement the internodal growth (see @stem::simulate)
	double delayNGStart = 0.;		///< delay between stem creation and start of nodal growth [day]
	double delayNGStarts = 0.;		///< delay between stem creation and start of nodal growth, deviation [day]
	double delayNGEnd = 0.;		///< delay between stem creation and start of nodal growth [day]
	double delayNGEnds = 0.;		///< delay between stem creation and start of nodal growth, deviation [day]

    /*
     * Callback functions for the Stem (set up by the class StemSystem)
     */
    std::shared_ptr<SoilLookUp> f_se = std::make_shared<SoilLookUp>(); ///< scale elongation function
    std::shared_ptr<SoilLookUp> f_sa = std::make_shared<SoilLookUp>(); ///< scale angle function
    std::shared_ptr<SoilLookUp> f_sbp = std::make_shared<SoilLookUp>(); ///< scale branching probability functiongrowth

protected:

    void bindParameters() override; ///<sets up class introspection

};

} // end namespace CStemBox

#endif
