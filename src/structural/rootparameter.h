// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ROOTPARAMETER_H_
#define ROOTPARAMETER_H_

#include "mymath.h"
#include "soil.h"
#include "growth.h"
#include "organparameter.h"

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

    RootSpecificParameter(): RootSpecificParameter(-1, 0., 0., std::vector<double>(0), 0., 0., 0., 0.) { } ///< Default constructor
    RootSpecificParameter(int type, double lb, double la, const std::vector<double>& ln, double r, double a,
    		double theta, double rlt, bool laterals = false, double a_growth = 0., double rlt_winter = 1e10):
        OrganSpecificParameter(type, a),  lb(lb), la(la), ln(ln), r(r), theta(theta), rlt(rlt), laterals(laterals), a_gr(a_growth), rlt_winter(rlt_winter) { }; ///< Constructor setting all parameters

    /*
     * RootBox parameters per single root
     */
    double lb;              ///< Basal zone [cm]
    double la;              ///< Apical zone [cm]
    std::vector<double> ln; ///< Inter-lateral distances [cm]
    double r;               ///< Initial growth rate [cm day-1]
    double theta;           ///< Angle between root and parent root [rad]
    double rlt;             ///< Root life time [day]

    bool laterals = false;
    double a_gr;            ///< Secondary growth rate [cm day-1]
    int nob() const { return ln.size()+ laterals; } ///< return the maximal number of lateral branching nodes [1]
    double getK() const; ///< Returns the exact maximal root length of this realization [cm]
	double rlt_winter = 1.e10;

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

    double nob() const { if(ln>0){ return std::max((lmax-la-lb)/ln+1, 0.);}else{return 0.;} }  ///< returns the mean maximal number of branching nodes [1]
    double nobs() const; ///< returns the standard deviation of number of branching nodes [1]

    std::string toString(bool verbose = true) const override; ///< info for debugging

    void readXML(tinyxml2::XMLElement* element, bool verbose) override; ///< reads a single sub type organ parameter set

    // DEPRICATED
    void read(std::istream & cin); ///< reads a single root parameter set
    void write(std::ostream & cout) const; ///< writes a single root parameter set

    void bindParameters() override; ///<sets up class introspection

    /*
     * RootBox parameters per root type
     */
    double lb = 0.;         ///< Basal zone [cm]
    double lbs = 0.;        ///< Standard deviation basal zone [cm]
    double la = 10.;        ///< Apical zone [cm];
    double las = 0.;        ///< Standard deviation apical zone [cm];
    double ln = 1;          ///< Inter-lateral distance [cm]
    double lns = 0.;        ///< Standard deviation inter-lateral distance [cm]
    double lmax = 0.;       ///< Maximal root length [cm]
    double lmaxs = 0.;      ///< Standard deviation of maximal root length [cm]
    double r = 1;           ///< Initial growth rate [cm day-1]
    double rs = 0.;         ///< Standard deviation initial growth rate [cm day-1]
    double a_gr = 0.;        ///< Secondary growth rate [cm day-1]
    double a_grs = 0.;      ///< Standard deviation Secondary growth rate [cm day-1]
    int tropismT = 1;       ///< Root tropism parameter (Type)
    //int tropismT1 = 1;       ///< Root tropism parameter (Type)
    //int tropismT2 = 1;       ///< Root tropism parameter (Type)
    double tropismW1 = 1;       ///< Root tropism parameter (Type)
    double tropismW2 = 0;       ///< Root tropism parameter (Type)
    double tropismN = 1.;   ///< Root tropism parameter (number of trials)
    double tropismS = 0.2;  ///< Root tropism parameter (mean value of expected changeg) [1/cm]
    double theta = 1.22;    ///< Angle between root and parent root (rad)
    double thetas= 0.;      ///< Standard deviation angle between root and parent root (rad)
    double rlt = 1e9;       ///< Root life time (days)
    double rlts = 0.;       ///< Standard deviation root life time (days)
    int gf = 1;             ///< Growth function (1=negative exponential, 2=linear)

    // new
    double lnk = 0.;        ///< Slope of inter-lateral distances [1]
	int is_fine_root = 0; // fine [1] or long-lived [0] root? maybe better way to define that.
	double k_survive = 1.88;
	double lambda_survive = 7.53;

    // root hair
    double hairsElongation = 0.; // zone behind the tip without root hairs [cm]
    double hairsZone = 0.; // root hair zone [cm]
	double hairsLength = 0.; // root hair length [cm]

    /*
     * Callback functions for the Root (set up by the class RootSystem)
     */
    std::shared_ptr<SoilLookUp> f_se = std::make_shared<SoilLookUp>(); ///< scale elongation function
    std::shared_ptr<SoilLookUp> f_sa = std::make_shared<SoilLookUp>(); ///< scale angle function
    std::shared_ptr<SoilLookUp> f_sbp = std::make_shared<SoilLookUp>(); ///< scale branching probability function growth

protected:

    double snap(double x);

};

} // end namespace CPlantBox

#endif
