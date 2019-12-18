// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ROOTPARAMETER_H_
#define ROOTPARAMETER_H_

#include "mymath.h"
#include "soil.h"
#include "growth.h"
#include "organparameter.h"
#include "tropism.h"

/**
 * This file describes the classes RootSpecificParameter and RootRandomParameter.
 * RootSpecificParameter are drawn from the RootRandomParameter class
 */

namespace CPlantBox {

class Organism;

/**
 * Parameters of a specific root, its created by RootRandomParameter:realize()
 */
class RootSpecificParameter :public OrganSpecificParameter
{

public:

    RootSpecificParameter(): RootSpecificParameter(-1,0.,0.,std::vector<double>(0),0,0.,0.,0.,0.) { } ///< Default constructor
    RootSpecificParameter(int type, double lb, double la, const std::vector<double>& ln, int nob, double r, double a, double theta, double rlt):
        OrganSpecificParameter(type, a),  lb(lb), la(la), r(r), theta(theta), rlt(rlt), ln(ln) { }; ///< Constructor setting all parameters

    /*
     * RootBox parameters per single root
     */
    double lb;              ///< Basal zone [cm]
    double la;              ///< Apical zone [cm]
    double r;               ///< Initial growth rate [cm day-1]
    double theta;           ///< Angle between root and parent root [rad]
    double rlt;             ///< Root life time [day]
    std::vector<double> ln; ///< Inter-lateral distances [cm]

    int nob() const { return ln.size(); } ///< return the maximal number of lateral branches [1]
    double getK() const; ///< Returns the exact maximal root length of this realization [cm]

    std::string toString() const override; ///< for debugging

};



/**
 * Contains a parameter set describing a root type
 */
class RootRandomParameter :public OrganRandomParameter
{

public:

    RootRandomParameter(std::shared_ptr<Organism> plant); ///< default constructor
    virtual ~RootRandomParameter() { };

    std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant) override;

    std::shared_ptr<OrganSpecificParameter> realize() override; ///< Creates a specific root from the root parameter set

    int getLateralType(const Vector3d& pos); ///< Choose (dice) lateral type based on root parameter set
    double nob() const { return std::max((lmax-la-lb)/ln+1, 1.); }  ///< returns the mean number of branches [1]
    double nobs() const; ///< returns the standard deviation of number of branches [1]

    std::string toString(bool verbose = true) const override; ///< info for debugging

    void readXML(tinyxml2::XMLElement* element) override; ///< reads a single sub type organ parameter set
    tinyxml2::XMLElement* writeXML(tinyxml2::XMLDocument& doc, bool comments = true) const override; ///< writes a organ root parameter set

    // DEPRICATED
    void read(std::istream & cin); ///< reads a single root parameter set
    void write(std::ostream & cout) const; ///< writes a single root parameter set

    void bindParameters() override; ///<sets up class introspection


    /*
     * RootBox parameters per root type
     */
    double lb = 0.; 	    ///< Basal zone [cm]
    double lbs = 0.;        ///< Standard deviation basal zone [cm]
    double la = 10.;	    ///< Apical zone [cm];
    double las = 0.;    	///< Standard deviation apical zone [cm];
    double ln = 1; 		    ///< Inter-lateral distance [cm]
    double lns = 0.;    	///< Standard deviation inter-lateral distance [cm]
    double lmax = 0.;    	///< Maximal root length [cm]
    double lmaxs = 0.;   	///< Standard deviation of maximal root length [cm]
    double r = 1;		    ///< Initial growth rate [cm day-1]
    double rs = 0.;	    	///< Standard deviation initial growth rate [cm day-1]
    double colorR = 0.6;	///< Root color (red)
    double colorG = 0.2;	///< Root color (green)
    double colorB = 0.2;	///< Root color (blue)
    int tropismT = 1;	    ///< Root tropism parameter (Type)
    double tropismN = 1.;   ///< Root tropism parameter (number of trials)
    double tropismS = 0.2;  ///< Root tropism parameter (mean value of expected changeg) [1/cm]
    double dx = 0.25; 		///< Maximal segment size [cm]
    double theta = 1.22; 	///< Angle between root and parent root (rad)
    double thetas= 0.; 	    ///< Standard deviation angle between root and parent root (rad)
    double rlt = 1e9;		///< Root life time (days)
    double rlts = 0.;	    ///< Standard deviation root life time (days)
    int gf = 1;			    ///< Growth function (1=negative exponential, 2=linear)
    std::vector<int> successor = std::vector<int>(0);			///< Lateral types [1]
    std::vector<double> successorP = std::vector<double>(0);  	///< Probabilities of lateral type to emerge (sum of values == 1) [1]

    /*
     * Callback functions for the Root (set up by the class RootSystem)
     */
    std::shared_ptr<Tropism> f_tf;  ///< tropism function (defined in constructor as new Tropism(plant))
    std::shared_ptr<GrowthFunction> f_gf = std::make_shared<ExponentialGrowth>(); ///< growth function
    std::shared_ptr<SoilLookUp> f_se = std::make_shared<SoilLookUp>(); ///< scale elongation function
    std::shared_ptr<SoilLookUp> f_sa = std::make_shared<SoilLookUp>(); ///< scale angle function
    std::shared_ptr<SoilLookUp> f_sbp = std::make_shared<SoilLookUp>(); ///< scale branching probability functiongrowth

};

} // end namespace CPlantBox

#endif
