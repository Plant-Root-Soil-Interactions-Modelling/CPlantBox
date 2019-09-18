// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef STEMPARAMETER_H_
#define STEMPARAMETER_H_

#include "mymath.h"
#include "soil.h"
#include "growth.h"
#include "organparameter.h"
#include "tropism.h"

/**
 * This file describes the classes StemSpecificParameter and StemRandomParameter.
 * StemSpecificParameter are drawn from the StemRandomParameter class
 */

namespace CRootBox {

class Organism;

/**
 * Parameters of a specific stem, its created by StemRandomParameter:realize()
 */
class StemSpecificParameter :public OrganSpecificParameter
{

public:

    StemSpecificParameter(): StemSpecificParameter(-1,0.,0.,std::vector<double>(0),0,0.,0.,0.,0.) { } ///< Default constructor
    StemSpecificParameter(int type, double lb, double la, const std::vector<double>& ln, int nob, double r, double a, double theta, double rlt):
        OrganSpecificParameter(),  lb(lb), la(la), nob(nob), r(r), a(a), theta(theta), rlt(rlt), ln(ln) { subType = type; } ///< Constructor setting all parameters

    /*
     * StemBox parameters per single stem
     */
    double lb;              ///< Basal zone [cm]
    double la;              ///< Apical zone [cm]
    int nob;                ///< Number of branches [1], ln.size()== nob-1 for nob>0
    double r;               ///< Initial growth rate [cm day-1]
    double a;               ///< Stem radius [cm]
    double theta;           ///< Angle between stem and parent stem [rad]
    double rlt;             ///< Stem life time [day]
    std::vector<double> ln; ///< Inter-lateral distances [cm]

    double getK() const; ///< Returns the exact maximal stem length of this realization [cm]

    std::string toString() const override; ///< for debugging

};



/**
 * Contains a parameter set describing a stem type
 */
class StemRandomParameter :public OrganRandomParameter
{

public:

    StemRandomParameter(Organism* plant); ///< default constructor
    virtual ~StemRandomParameter();

    OrganRandomParameter* copy(Organism* plant_) override;

    OrganSpecificParameter* realize() override; ///< Creates a specific stem from the stem parameter set
    int getLateralType(const Vector3d& pos); ///< Choose (dice) lateral type based on stem parameter set
    double getK() const { return std::max(nob-1,double(0))*ln+la+lb; }  ///< returns the mean maximal stem length [cm]

    std::string toString(bool verbose = true) const override; ///< info for debugging

    void readXML(tinyxml2::XMLElement* element) override; ///< reads a single sub type organ parameter set
    tinyxml2::XMLElement* writeXML(tinyxml2::XMLDocument& doc, bool comments = true) const override; ///< writes a organ stem parameter set

    /*
     * Parameters per stem type
     */
    double lb = 0.; 	    ///< Basal zone [cm]
    double lbs = 0.;        ///< Standard deviation basal zone [cm]
    double la = 10.;	    ///< Apical zone [cm];
    double las = 0.;    	///< Standard deviation apical zone [cm];
    double ln = 1; 		    ///< Inter-lateral distance [cm]
    double lns = 0.;    	///< Standard deviation inter-lateral distance [cm]
    int lnf = 0;            ///< Inter-lateral type [type]
    double nob = 0.;    	///< Number of branches [1]
    double nobs = 0.;   	///< Standard deviation number of branches [1]
    double r = 1;		    ///< Initial growth rate [cm day-1]
    double rs = 0.;	    	///< Standard deviation initial growth rate [cm day-1]
    double a = 0.1; 		///< Stem radius [cm]
    double as = 0.; 		///< Standard deviation stem radius [cm]
    double k = 0.;          ///< Maximal stem length [cm]
    double RotBeta = 0.6;	///< Revrotation
    double BetaDev = 0.2;	///< Deviation of RevRotation
    double InitBeta = 0.2;	///< Initial RevRotation
    int tropismT = 1;	    ///< Stem tropism parameter (Type)
    double tropismN = 1.;   ///< Stem tropism parameter (number of trials)
    double tropismS = 0.2;  ///< Stem tropism parameter (mean value of expected changeg) [1/cm]
    double dx = 0.25; 		///< Maximal segment size [cm]
    double theta = 1.22; 	///< Angle between stem and parent stem (rad)
    double thetas= 0.; 	    ///< Standard deviation angle between stem and parent stem (rad)
    double rlt = 1e9;		///< Stem life time (days)
    double rlts = 0.;	    ///< Standard deviation stem life time (days)
    int gf = 1;			    ///< Growth function (1=negative exponential, 2=linear)
    std::vector<int> successor = std::vector<int>(0);			///< Lateral types [1]
    std::vector<double> successorP = std::vector<double>(0);  	///< Probabilities of lateral type to emerge (sum of values == 1) [1]

    /*
     * Callback functions for the Stem (set up by the class StemSystem)
     */
    Tropism* f_tf;  ///< tropism function ( = new Tropism(plant) )
    GrowthFunction* f_gf = new ExponentialGrowth(); ///< growth function
    SoilLookUp* f_se = new SoilLookUp(); ///< scale elongation function
    SoilLookUp* f_sa = new SoilLookUp(); ///< scale angle function
    SoilLookUp* f_sbp = new SoilLookUp(); ///< scale branching probability function

protected:

    void bindParmeters(); ///<sets up class introspection

};

} // end namespace CStemBox

#endif
