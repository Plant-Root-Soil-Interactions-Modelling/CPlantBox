// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "hyphaeparameter.h"

#include "Organism.h"
#include "tropism.h"

#include <cmath>
#include <iostream>
#include <chrono>
#include <assert.h>
#include <numeric>

namespace CPlantBox {

/**
 * @copydoc OrganParameter::toString()
 */
std::string HyphaeSpecificParameter::toString() const
{
    std::stringstream str;
    str << "HyphaeSpecificParameter: subType\t" << subType << ", radius " << a << ", " << "tip elongation rate " << v << ", " << "branching rate " << b << ", " << "hyphal lifetime " << hlt << ", branching angle " << theta << std::endl; 
    return str.str();
}


/**
 * @return Mean maximal hyphae length of this hyphae type
 */
double HyphaeSpecificParameter::getMaxLength() const {
    // double l = 0.5 / (order+1); // v*hlt is too long ie tip elongation rate * hyphal lifetime
    double l = v * hlt; // mean maximal length of a hyphae
    return l; 
}

/**
 * Default constructor sets up hashmaps for class introspection
 */
HyphaeRandomParameter::HyphaeRandomParameter(std::shared_ptr<Organism> plant) :OrganRandomParameter(plant)
{
    // base class default values
    name = "undefined";
    organType = Organism::ot_hyphae;
    subType = -1;
    f_tf = std::make_shared<Tropism>(plant);
    bindParameters();
    f_gf = std::make_shared<LinearGrowth>();
}

/**
 * @copydoc OrganTypeParameter::copy()
 */
std::shared_ptr<OrganRandomParameter> HyphaeRandomParameter::copy(std::shared_ptr<Organism> p)
{
    auto r = std::make_shared<HyphaeRandomParameter>(*this); // copy constructor breaks class introspection
    r->plant = std::weak_ptr<Organism>(p);
    r->plant = p;
    r->bindParameters(); // fix class introspection
    r->f_tf = f_tf->copy(p); // copy call back function classes
    r->f_gf = f_gf->copy();
    return r;
}

/**
 * @copydoc OrganTypeParameter::realize()
 *
 * Creates a specific hyphae from the hyphae random parameters.
 * @return Specific hyphae parameters derived from the root type parameters
 */
std::shared_ptr<OrganSpecificParameter> HyphaeRandomParameter::realize()
{
    assert(dx > dxMin && "HyphaeRandomParameter::realize(): dxMin must be smaller than dx");
    auto p = plant.lock();

    double a_ = std::max(a + p->randn()*as, 0.); // radius
    double lb_ = std::max(lb + p->randn()*lbs, 0.); // basal zone [cm]
    double la_ = std::max(la + p->randn()*las, 0.); // apical zone [cm]
    double v_ = std::max(v + p->randn()*vs, 0.);  // tip elongation rate [cm/day]
    double b_ = std::max(b + p->randn()*bs, 0.); // branching rate [1/day]
    double hlt_ = std::max(hlt + p->randn()*hlts, 0.); // hyphal lifetime  [day]
    double theta_ = std::max(theta + p->randn()*thetas, 0.); // branching angle [rad]
    std::vector<double> ln_; // stores the inter-distances
    double nob_sd = p->randn()*nobs();
    int nob_real = round(std::max(nob() + nob_sd, 0.)); // real maximal number of branching points
    bool hasLaterals = (successorST.size()>0) && (nob_real>0);

    if (!hasLaterals) { // no laterals
        lb_ = 0;
        la_ = std::max(lmax + p->randn()*lmaxs, 0.); // la, and lb is ignored
        la_ = snap(la_);
    } else { // laterals
        lb_ = snap(std::max(lb + p->randn()*lbs, 0.)); // length of basal zone
        la_ = snap(std::max(la + p->randn()*las, 0.)); // length of apical zone
        double ln_mean = ln;
        if(ln < dxMin && ln != 0) { // limit to minimum resolution
            ln_mean = dxMin;
        }
        double nob1 = std::max((lmax-la_-lb_)/ln_mean+1, 0.); // use new la_, lb_ and ln_mean
        int nob_ = std::min(std::max(round(nob1 + nob_sd), 0.), double(nob_real)); // maximal number of branches +1
        int latMissing = nob_real - nob_;
        assert((latMissing >= 0) && "RootRandomParameter::realize(): latMissing < 0");
        int latExtraMean = floor(latMissing/nob_real); // mean number of extra laterals per branching point to keep correct number
        int latExtra = latMissing - latExtraMean*(nob_);
        for (int j = 0; j<latExtraMean; j++) { //at end of basal zone
            ln_.push_back(0);
        }
        if (latExtra> 0) { //at end of basal zone
            ln_.push_back(0);
            latExtra--;
        }
        double sum_ln = nob_*ln_mean; // mean length of lateral zone

        for (int i = 0; i<nob_-1; i++) { // create inter-root distances
            double z = ((double)i+0.5)*ln_mean; // regular position along root lateral zone
            double f = lnk*(z-sum_ln/2.); // evaluate slope lnk f(mid) = 0
            double pf = (ln_mean + f) / ln_mean; // we scale lns by the change in percentage
            double d = std::max(ln_mean + f + pf*p->randn()*lns, 1.e-5); // miminum is 1.e-5
            d = snap(d);
            ln_.push_back(d);
            for (int j = 0; j<latExtraMean; j++) {
                ln_.push_back(0);
            }
            if (latExtra> 0) {
                ln_.push_back(0);
                latExtra--;
            }
        }
    }

     return std::make_shared<HyphaeSpecificParameter>(subType, a_, lb_, la_, ln_, v_, b_, hlt_, theta_, hasLaterals);
}

/**
 * @copydoc OrganTypeParameter::toString()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 */
std::string HyphaeRandomParameter::toString(bool verbose) const {

    if (verbose) {
        return OrganRandomParameter::toString(true);
    } else {
        return OrganRandomParameter::toString(false);
    }

}

/**
 * @copydoc OrganTypeParameter::readXML()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 *
 * If the parameter successor or successorP are not in the element, they are set to zero size.
 */
void HyphaeRandomParameter::readXML(tinyxml2::XMLElement* element, bool verbose)
{
    OrganRandomParameter::readXML(element, verbose);
}


/**
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void HyphaeRandomParameter::bindParameters()
{
    OrganRandomParameter::bindParameters();
    bindParameter("v", &v, "Tip elongation rate [cm/day]", &vs);
    bindParameter("b", &b, "Branching rate [1/day]", &bs);
    bindParameter("hlt", &hlt, "Hyphal lifetime  [day]", &hlts);
    bindParameter("theta", &theta, "Branching angle [rad]", &thetas);
    // bindParameter("distTT", &distTT, "Distance for tip-tip anastomosis [cm]");
    bindParameter("distTH", &distTH, "Distance for tip-hyphae anastomosis [cm]");
    bindParameter("ana", &ana, "Probability of anastomosis occuring if distance is long enough");
    bindParameter("tropismT", &tropismT, "Type of root tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)");
    bindParameter("tropismN", &tropismN, "Number of trials of root tropism");
    bindParameter("tropismS", &tropismS, "Mean value of expected change of root tropism [1/cm]");
}

/**
 * snaps to the grid, to make it compatible with dx() and dxMin() copied from RootRandomPArameter
 */
double HyphaeRandomParameter::snap(double x) const
{
    double res = x - floor(x / dx)*dx; // res < dx
    if ((res < dxMin) && (res != 0)) { //make ln compatible with dx() and dxMin().
        if(res <= dxMin/2){
            x -= res;
        } else {
            x = floor(x / dx)*dx + dxMin;
        }
    }
    return x;
}
double HyphaeRandomParameter::nobs() const
{
    double nobs = 0;
    if(ln >0)
    {
        nobs = (lmaxs/lmax - lns/ln)*lmax/ln; // error propagation
        if (la>0) {
            nobs -= (las/la - lns/ln)*la/ln;
        }
        if (lb>0) {
            nobs -= (lbs/lb - lns/ln)*lb/ln;
        }
    }
    return std::max(nobs,0.);
}

} // end namespace CPlantBox
