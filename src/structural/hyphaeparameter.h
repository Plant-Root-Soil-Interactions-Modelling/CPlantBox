// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef HYPHAEPARAMETER_H_
#define HYPHAEPARAMETER_H_

#include "mymath.h"
#include "soil.h"
#include "growth.h"
#include "organparameter.h"

/**
 * This file describes the classes HyphaeSpecificParameter and HyphaeRandomParameter.
 * HyphaeSpecificParameter are drawn from the HyphaeRandomParameter class
 */

namespace CPlantBox {

class Organism;

/**
 * Parameters of a specific root, its created by RootRandomParameter:realize()
 */
class HyphaeSpecificParameter :public OrganSpecificParameter
{
public:

    HyphaeSpecificParameter(): HyphaeSpecificParameter(-1,0., 0., 0., std::vector<double>(0), 0., 0., 0., 0., false) { } ///< Default constructor
    HyphaeSpecificParameter(int subType, double a, double lb, double la, const std::vector<double>& ln,  double v, double b, double hlt, double theta, bool laterals = true):
            OrganSpecificParameter(subType, a), lb(lb), la(la), ln(ln), v(v), b(b), hlt(hlt), theta(theta), laterals(laterals) { }; ///< Constructor setting all parameters

    // int order = 0; // internal counter (? todo)
    double lb;              ///< Basal zone [cm]
    double la;              ///< Apical zone [cm]
    std::vector<double> ln;  ///< Inter-lateral distance [cm]
    double v;              ///< tip elongation rate [cm/day]
    double b;              ///< branching rate [1/day]
    double hlt;            ///< hyphal lifetime  [day]
    double theta;          ///< branching angle [rad]
    bool laterals;         ///< indicates if hyphae has laterals

    int nob() const { return ln.size()+ laterals; } ///< return the maximal number of lateral branching nodes [1]
    std::string toString() const override; ///< for debugging
    double getMaxLength() const; ///< returns maximal length
};


/**
 * Contains a parameter set describing a root type
 */
class HyphaeRandomParameter :public OrganRandomParameter
{

public:

    HyphaeRandomParameter(std::shared_ptr<Organism> plant); ///< default constructor
    virtual ~HyphaeRandomParameter() { };

    std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant) override;

    std::shared_ptr<OrganSpecificParameter> realize() override; ///< Creates a specific root from the root parameter set

    std::string toString(bool verbose = true) const override; ///< info for debugging

    double nob() const { if(ln>0){ return std::max((lmax-la-lb)/ln+1, 0.);}else{return 0.;} }  ///< returns the mean maximal number of branching nodes [1]
    double nobs() const; ///< returns the standard deviation of number of branching nodes [1]

    void readXML(tinyxml2::XMLElement* element, bool verbose) override; ///< reads a single sub type organ parameter set

    void bindParameters() override; ///<sets up class introspection

    double lb = 0.0001;         ///< Basal zone [cm]
    double lbs = 0.;        ///< Standard deviation basal zone [cm]
    double la = 0.003;        ///< Apical zone [cm];
    double las = 0.;        ///< Standard deviation apical zone [cm];
    double lmax = 10.;       ///< Maximal length of the hyphae [cm]
    double lmaxs = 0.;       ///< Standard deviation of maximal length of the hyphae [cm]
    double ln = 0.005;          ///< Inter-lateral distance [cm]
    double lns = 0.;        ///< Standard deviation inter-lateral distance [cm]
    double v=0.13;              ///< tip elongation rate [cm/day] 
    double vs=0.01;             ///< standard deviation of tip elongation rate [cm/day]
    double b=0.5;              ///< branching rate [1/day]
    double bs=0.05;             ///< standard deviation of branching rate [1/day]
    double hlt=10;            ///< hyphal lifetime  [day]
    double hlts=0.1;           ///< standard deviation of  hyphal lifetime  [day]
    double theta=30./180.*M_PI;          ///< branching angle [rad]
    double thetas=0./180.*M_PI;         ///< standard deviation of branching angle  [rad]
    // double distTT = 0.; ///< distance for tip tip anastomosis [cm] 
    double distTH = 0.; ///< distance for tip hyphae anastomosis [cm]
    double ana= 1.0; ///< Probability of anastomosis occuring if tip is close enough

    double lnk = 0.; //< Slope of inter-lateral disntances [1]

    int tropismT = 2;       ///< Hypha tropism parameter (Type)
    double tropismN = 1.;   ///< Hypha tropism parameter (number of trials)
    double tropismS = 0.3;  ///< Hypha tropism parameter (mean value of expected changeg [1/cm]

    std::shared_ptr<SoilLookUp> v_scale = std::make_shared<SoilLookUp>(); ///< elongation rate scale
    std::shared_ptr<SoilLookUp> b_scale = std::make_shared<SoilLookUp>(); ///< scale branching rate

    protected:

    double snap(double x) const; ///< snap value to grid defined by dxMin
};

} // end namespace CPlantBox

#endif
