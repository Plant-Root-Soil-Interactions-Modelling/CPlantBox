// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "rootparameter.h"

#include "Organism.h"
#include "tropism.h"

#include <assert.h>
#include <cmath>
#include <iostream>
#include <numeric>

namespace CPlantBox {

/**
 * @return Mean maximal root length of this root type
 */
double RootSpecificParameter::getK() const {
    double l = std::accumulate(ln.begin(), ln.end(), 0.);
    return l + la + lb;
}

/**
 * @copydoc OrganParameter::toString()
 */
std::string RootSpecificParameter::toString() const {
    std::stringstream str;
    str << "subType\t" << subType << std::endl;
    str << "lb\t" << lb << std::endl << "la\t" << la << std::endl;
    str << "nob\t" << nob() << std::endl << "r\t" << r << std::endl << "a\t" << a << std::endl;
    str << "theta\t" << theta << std::endl << "rlt\t" << rlt << std::endl;
    str << "ln\t";
    for (int i = 0; i < ln.size(); i++) {
        str << ln[i] << " ";
    }
    str << std::endl;
    return str.str();
}

/**
 * Default constructor sets up hashmaps for class introspection
 */
RootRandomParameter::RootRandomParameter(std::shared_ptr<Organism> plant) : OrganRandomParameter(plant) {
    // base class default values
    name = "undefined";
    organType = Organism::ot_root;
    subType = -1;
    f_tf = std::make_shared<Tropism>(plant);
    bindParameters();
}

/**
 * @copydoc OrganTypeParameter::copy()
 */
std::shared_ptr<OrganRandomParameter> RootRandomParameter::copy(std::shared_ptr<Organism> pl) {
    // std::cout << "RootRandomParameter::copy() into plant "<< (p->plantId) <<  "\n"<< std::flush;
    auto r = std::make_shared<RootRandomParameter>(*this); // copy constructor breaks class introspection
    r->plant = std::weak_ptr<Organism>();                  // necessary?
    r->plant = pl;
    r->bindParameters();      // fix class introspection
    r->f_tf = f_tf->copy(pl); // copy call back function classes
    r->f_gf = f_gf->copy();
    r->f_se = f_se; // for carbon limited grow we want the same reference
    r->f_sa = f_sa->copy();
    r->f_sbp = f_sbp->copy();
    return r;
}

/**
 * @copydoc OrganTypeParameter::realize()
 *
 * Creates a specific root from the root type parameters.
 * @return Specific root parameters derived from the root type parameters
 */
std::shared_ptr<OrganSpecificParameter> RootRandomParameter::realize() {
    assert(dx > dxMin && "RootRandomParameter::realize(): dxMin must be smaller than dx");
    auto p = std::make_shared<RootSpecificParameter>();
    auto pl = plant.lock();
    p->subType = subType;
    p->a = std::max(a + pl->randn() * as, 0.);
    p->ldelay = std::max(ldelay + pl->randn() * ldelays, 0.);
    std::vector<double> ln_; // stores the inter-distances
    double nob_sd = pl->randn() * nobs();
    int nob_real = round(std::max(nob() + nob_sd, 0.)); // real maximal number of branching points
    bool hasLaterals = (successorST.size() > 0) && (nob_real > 0);
    if (!hasLaterals) { // no laterals
        p->lb = 0;
        double la_ = std::max(lmax + pl->randn() * lmaxs, 0.); // la, and lb is ignored
        p->la = snap(la_);                                     // lmax is moved to la
    } else {                                                   // laterals
        p->lb = snap(std::max(lb + pl->randn() * lbs, 0.));    // length of basal zone
        p->la = snap(std::max(la + pl->randn() * las, 0.));    // length of apical zone
        double ln_mean = ln;
        if (ln < dxMin && ln != 0) { // limit to minimum resolution
            ln_mean = dxMin;
        }
        double nob1 = std::max((lmax - p->la - p->lb) / ln_mean + 1, 0.);          // use new la_, lb_ and ln_mean
        int nob_ = std::min(std::max(round(nob1 + nob_sd), 0.), double(nob_real)); // maximal number of branches +1
        int latMissing = nob_real - nob_;
        assert((latMissing >= 0) && "RootRandomParameter::realize(): latMissing < 0");
        int latExtraMean = floor(latMissing / nob_real); // mean number of extra laterals per branching point to keep correct number
        int latExtra = latMissing - latExtraMean * (nob_);
        for (int j = 0; j < latExtraMean; j++) { // at end of basal zone
            ln_.push_back(0);
        }
        if (latExtra > 0) { // at end of basal zone
            ln_.push_back(0);
            latExtra--;
        }
        double sum_ln = nob_ * ln_mean;                                       // mean length of lateral zone
        for (int i = 0; i < nob_ - 1; i++) {                                  // create inter-root distances
            double z = ((double)i + 0.5) * ln_mean;                           // regular position along root lateral zone
            double f = lnk * (z - sum_ln / 2.);                               // evaluate slope lnk f(mid) = 0
            double pf = (ln_mean + f) / ln_mean;                              // we scale lns by the change in percentage
            double d = std::max(ln_mean + f + pf * pl->randn() * lns, 1.e-5); // miminum is 1.e-5
            d = snap(d);
            ln_.push_back(d);
            for (int j = 0; j < latExtraMean; j++) {
                ln_.push_back(0);
            }
            if (latExtra > 0) {
                ln_.push_back(0);
                latExtra--;
            }
        }
    }
    p->laterals = hasLaterals;
    p->ln = ln_;
    p->r = std::max(r + pl->randn() * rs, 0.);             // initial elongation
    p->theta = std::max(theta + pl->randn() * thetas, 0.); // initial elongation
    p->rlt = std::max(rlt + pl->randn() * rlts, 0.);       // root life time
    return p;
}

/**
 * snaps to the grid, to make it compatible with dx() and dxMin()
 */
double RootRandomParameter::snap(double x) {
    double res = x - floor(x / dx) * dx; // res < dx
    if ((res < dxMin) && (res != 0)) {   // make ln compatible with dx() and dxMin().
        if (res <= dxMin / 2) {
            x -= res;
        } else {
            x = floor(x / dx) * dx + dxMin;
        }
    }
    return x;
}

/**
 * @return Mean maximal number of branching nodes
 */
double RootRandomParameter::nob() const {
    if (ln > 0) {
        return std::max((lmax - la - lb) / ln + 1, 0.);
    }
    return 0.;
}

/**
 * The standard deviation of number of branches is calculated form root maximal length lmax, and
 * its standard deviation lmaxs
 *
 * todo I have no idea why this holds...
 */
double RootRandomParameter::nobs() const {
    double nobs = 0;
    if (ln > 0) {
        nobs = (lmaxs / lmax - lns / ln) * lmax / ln; // error propagation
        if (la > 0) {
            nobs -= (las / la - lns / ln) * la / ln;
        }
        if (lb > 0) {
            nobs -= (lbs / lb - lns / ln) * lb / ln;
        }
    }
    return std::max(nobs, 0.);
}

/**
 * @copydoc OrganTypeParameter::toString()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 */
std::string RootRandomParameter::toString(bool verbose) const {

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
void RootRandomParameter::readXML(tinyxml2::XMLElement *element, bool verbose) { OrganRandomParameter::readXML(element, verbose); }

/**
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void RootRandomParameter::bindParameters() {
    OrganRandomParameter::bindParameters();
    bindParameter("lb", &lb, "Basal zone [cm]", &lbs);
    bindParameter("la", &la, "Apical zone [cm]", &las);
    bindParameter("ln", &ln, "Inter-lateral distance [cm]", &lns);
    bindParameter("lmax", &lmax, "Maximal root length [cm]", &lmaxs);
    bindParameter("r", &r, "Initial growth rate [cm day-1]", &rs);
    bindParameter("tropismT", &tropismT, "Type of root tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)");
    bindParameter("tropismN", &tropismN, "Number of trials of root tropism");
    bindParameter("tropismS", &tropismS, "Mean value of expected change of root tropism [1/cm]");
    bindParameter("theta", &theta, "Angle between root and parent root [rad]", &thetas);
    bindParameter("rlt", &rlt, "Root life time [day]", &rlts);
    bindParameter("gf", &gf, "Growth function number [1]", &rlts);
    bindParameter("lnk", &lnk, "Slope of inter-lateral distances [1]");
    bindParameter("hairsElongation", &hairsElongation, "Zone behind the tip without root hairs  [cm]");
    bindParameter("hairsZone", &hairsZone, "Length of the root hair zone [cm]");
    bindParameter("hairsLength", &hairsLength, "Root hair length [cm]");
}

} // end namespace CPlantBox
