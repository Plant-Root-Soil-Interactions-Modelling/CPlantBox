// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "stemparameter.h"

#include "Organism.h"
#include "Seed.h"
#include "tropism.h"

#include <assert.h>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>

namespace CPlantBox {

/**
 * @return Mean maximal stem length of this stem type
 */
double StemSpecificParameter::getK() const {
    double l = std::accumulate(ln.begin(), ln.end(), 0.);
    return l + la + lb;
}

/**
 * @copydoc OrganParameter::toString()
 */
std::string StemSpecificParameter::toString() const {
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
StemRandomParameter::StemRandomParameter(std::shared_ptr<Organism> plant) : OrganRandomParameter(plant) {
    // base class default values
    name = "undefined";
    organType = Organism::ot_stem;
    subType = -1;
    f_tf = std::make_shared<Tropism>(plant);
    bindParameters();
}

/**
 * @copydoc OrganTypeParameter::copy()
 */
std::shared_ptr<OrganRandomParameter> StemRandomParameter::copy(std::shared_ptr<Organism> plant) {
    auto r = std::make_shared<StemRandomParameter>(*this); // copy constructor breaks class introspection
    r->plant = plant;
    r->bindParameters();         // fix class introspection
    r->f_tf = f_tf->copy(plant); // copy call back classes
    r->f_gf = f_gf->copy();
    r->f_se = f_se->copy();
    r->f_sa = f_sa->copy();
    r->f_sbp = f_sbp->copy();
    return r;
}

/**
 * @copydoc OrganTypeParameter::realize()
 *
 * Creates a specific stem from the stem type parameters.
 * @return Specific stem parameters derived from the stem type parameters
 */
std::shared_ptr<OrganSpecificParameter> StemRandomParameter::realize() {
    auto p = std::make_shared<StemSpecificParameter>();
    auto pl = plant.lock();
    p->subType = subType;
    p->a = std::max(a + pl->randn() * as, 0.);
    p->ldelay = std::max(ldelay + pl->randn() * ldelays, 0.);
    std::vector<double> ln_; // stores the inter-distances
    double res;
    double nob_sd = pl->randn() * nobs();
    int nob_real = round(std::max(nob() + nob_sd, 0.)); // real maximal number of branching points
    bool hasLaterals = (successorST.size() > 0);
    if (dx <= dxMin) {
        std::cout << "dx <= dxMin, dxMin set to dx/2" << std::endl;
        this->dxMin = dx / 2;
    }
    if (!hasLaterals) { // no laterals

        p->lb = 0;
        p->la = std::max(lmax + pl->randn() * lmaxs, 0.); // la, and lb is ignored
        res = p->la - floor(p->la / dx) * dx;
        if (res < dxMin && res != 0) {
            if (res <= dxMin / 2) {
                p->la -= res;
            } else {
                p->la = floor(p->la / dx) * dx + dxMin;
            }
        } // make la_ compatible with dx() and dxMin()

    } else {
        p->lb = std::max(lb + pl->randn() * lbs, 0.); // length of basal zone
        p->la = std::max(la + pl->randn() * las, 0.); // length of apical zone
        res = p->lb - floor(p->lb / dx) * dx;
        if ((res < dxMin) && (res != 0)) {
            if (res <= dxMin / 2) {
                p->lb -= res;
            } else {
                p->lb = floor(p->lb / dx) * dx + dxMin;
            }
        }

        bool hasSeed = (pl->baseOrgans.size() > 0) && (pl->baseOrgans.at(0)->organType() == Organism::ot_seed);
        if (hasSeed && (p->lb < dxMin * 2) && (pl->getSeed()->param()->nC > 0)) {
            p->lb = dxMin * 2;
        }

        res = p->la - floor(p->la / dx) * dx;
        if (res < dxMin && res != 0) {
            if (res <= dxMin / 2) {
                p->la -= res;
            } else {
                p->la = floor(p->la / dx) * dx + dxMin;
            }
        }
        double ln_mean = ln;
        if (ln < dxMin * 0.99 && ln != 0) {
            std::cout << "\nStemRandomParameter::realize inter-lateral distance (ln) " << ln << " below minimum resolution (dxMin) " << dxMin << ". ln set to dxMin" << std::endl;
            ln_mean = dxMin;
        }

        // adapt number of laterals by branching point to keep same total number of lats
        // in spite of dxMin
        int nob1 = std::max((lmax - p->la - p->lb) / ln_mean + 1, 1.);                           // use new la_, lb_ and ln_mean
        int nob_ = std::min(std::max(round(nob1 + pl->randn() * nobs()), 1.), double(nob_real)); // maximal number of branches
        int latMissing = nob_real - nob_;
        int latExtra1 = floor(latMissing / nob_); // mean number of extra laterals per branching point to keep correct number
        int latExtra2 = latMissing - latExtra1 * (nob_);
        int latExtra2_ = latExtra2;

        // at end of basal zone
        for (int j = 0; j < latExtra1; j++) {
            ln_.push_back(0);
        }
        if (latExtra2_ > 0) {
            ln_.push_back(0);
            latExtra2_--;
        }

        switch (lnf) {
        case 0:                                                          // homogeneously distributed stem nodes
            for (int i = 0; i < nob_ - 1; i++) {                         // create inter-stem distances
                double d = std::max(ln_mean + pl->randn() * lns, 1.e-5); // Normal function of equal internode distance
                res = d - floor(d / dx) * dx;
                if (res < dxMin && res != 0) {
                    if (res <= dxMin / 2) {
                        d -= res;
                    } else {
                        d = floor(d / dx) * dx + dxMin;
                    }

                } // make ln compatible with dx() and dxMin().

                ln_.push_back(d);
                for (int j = 0; j < latExtra1; j++) {
                    ln_.push_back(0);
                }
                if (latExtra2_ > 0) {
                    ln_.push_back(0);
                    latExtra2_--;
                }
            };
            break;
        case 1:                                                                    //  nodes distance increase linearly
            for (int i = 0; i < nob_ * 2 - 1; i++) {                               // create inter-stem distances
                double d = std::max(ln_mean * (1 + i) + pl->randn() * lns, 1.e-5); // std::max(  );//ln +p->randn()*lns,1e-9);
                res = d - floor(d / dx) * dx;
                if (res < dxMin && res != 0) {
                    if (res <= dxMin / 2) {
                        d -= res;
                    } else {
                        d = floor(d / dx) * dx + dxMin;
                    }

                } // make ln compatible with dx() and dxMin().

                ln_.push_back(d);
                for (int j = 0; j < latExtra1; j++) {
                    ln_.push_back(0);
                }
                if (latExtra2_ > 0) {
                    ln_.push_back(0);
                    latExtra2_ -= 0.5;
                }
                ln_.push_back(0);
                for (int j = 0; j < latExtra1; j++) {
                    ln_.push_back(0);
                }
                if (latExtra2_ > 0) {
                    ln_.push_back(0);
                    latExtra2_ -= 0.5;
                }
            };
            break;
        case 2:                                                                    // nodes distance decrease linearly
            for (int i = 0; i < nob_ - 1; i++) {                                   // create inter-stem distances
                double d = std::max(ln_mean * (1 + i) + pl->randn() * lns, 1.e-5); // std::max(  );//ln +p->randn()*lns,1e-9);
                res = d - floor(d / dx) * dx;
                if (res < dxMin && res != 0) {
                    if (res <= dxMin / 2) {
                        d -= res;
                    } else {
                        d = floor(d / dx) * dx + dxMin;
                    }

                } // make ln compatible with dx() and dxMin().

                ln_.push_back(d);
                for (int j = 0; j < latExtra1; j++) {
                    ln_.push_back(0);
                }
                if (latExtra2_ > 0) {
                    ln_.push_back(0);
                    latExtra2_--;
                }
            };
            break;
        case 3:                                                          // nodes distance decrease exponential
            for (int i = 0; i < nob_ - 1; i++) {                         // create inter-stem distances
                double d = std::max(ln_mean + pl->randn() * lns, 1.e-5); // std::max(  );//ln +p->randn()*lns,1e-9);
                res = d - floor(d / dx) * dx;
                if (res < dxMin && res != 0) {
                    if (res <= dxMin / 2) {
                        d -= res;
                    } else {
                        d = floor(d / dx) * dx + dxMin;
                    }

                } // make ln compatible with dx() and dxMin().

                ln_.push_back(d);
                for (int j = 0; j < latExtra1; j++) {
                    ln_.push_back(0);
                }
                if (latExtra2_ > 0) {
                    ln_.push_back(0);
                    latExtra2_--;
                }
            };
            break;

        case 4:                                                                    // nodes distance decrease exponential
            for (int i = 0; i < nob_ * 2 - 1; i++) {                               // create inter-stem distances
                double d = std::max(ln_mean / (1 + i) + pl->randn() * lns, 1.e-5); // std::max(  );//ln +p->randn()*lns,1e-9);
                res = d - floor(d / dx) * dx;
                if (res < dxMin && res != 0) {
                    if (res <= dxMin / 2) {
                        d -= res;
                    } else {
                        d = floor(d / dx) * dx + dxMin;
                    }

                } // make ln compatible with dx() and dxMin().

                ln_.push_back(d);
                for (int j = 0; j < latExtra1; j++) {
                    ln_.push_back(0);
                }
                if (latExtra2_ > 0) {
                    ln_.push_back(0);
                    latExtra2_ -= 0.5;
                }
                ln_.push_back(0);
                for (int j = 0; j < latExtra1; j++) {
                    ln_.push_back(0);
                }
                if (latExtra2_ > 0) {
                    ln_.push_back(0);
                    latExtra2_ -= 0.5;
                }
            };
            break;
        case 5:                                                                    // nodes distance decrease exponential
            for (int i = 0; i < nob_ * 2 - 1; i++) {                               // create inter-stem distances
                double d = std::max(ln_mean / (1 + i) + pl->randn() * lns, 1.e-5); // std::max(  );//ln +p->randn()*lns,1e-9);
                res = d - floor(d / dx) * dx;
                if (res < dxMin && res != 0) {
                    if (res <= dxMin / 2) {
                        d -= res;
                    } else {
                        d = floor(d / dx) * dx + dxMin;
                    }

                } // make ln compatible with dx() and dxMin().

                ln_.push_back(d);
                for (int j = 0; j < latExtra1; j++) {
                    ln_.push_back(0);
                }
                if (latExtra2_ > 0) {
                    ln_.push_back(0);
                    latExtra2_--;
                }
            };
            break;
        default:
            throw std::runtime_error("StemRandomParameter::realize type of inter-branching distance not recognized");
        }
    }
    p->laterals = hasLaterals;
    p->ln = ln_;
    p->r = std::max(r + pl->randn() * rs, 0.); // initial elongation
    p->ldelay = std::max(ldelay + pl->randn() * ldelays, 0.);
    p->theta = std::max(theta + pl->randn() * thetas, 0.); // initial elongation
    p->rlt = std::max(rlt + pl->randn() * rlts, 0.);       // stem life time
    p->delayNGStart = std::max(delayNGStart + pl->randn() * delayNGStarts, 0.);
    p->delayNGEnd = std::max(delayNGEnd + pl->randn() * delayNGEnds, 0.);
    if (p->delayNGEnd < p->delayNGStart) {
        std::cout << "StemRandomParameter::realize() : delayNGEnd_ < delayNGStart_ \n";
        std::cout << "set delayNGEnd_ = delayNGStart_ = " << p->delayNGStart << std::endl;
        p->delayNGEnd = p->delayNGStart;
    }
    return p;
}

/**
 * @return Mean maximal number of branching nodes
 */
double StemRandomParameter::nob() const {
    if (ln > 0) {
        return std::max((lmax - la - lb) / ln + 1, 0.);
    }
    return 0.;
}

/**
 * todo docme
 *
 * todo I have no idea why this holds...
 */
double StemRandomParameter::nobs() const {
    double nobs = 0;
    if (ln > 0) {
        double nobs = (lmaxs / lmax - lns / ln) * lmax / ln; // error propagation
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
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void StemRandomParameter::bindParameters() {
    OrganRandomParameter::bindParameters();
    bindParameter("lb", &lb, "Basal zone, hypocotyl for main stem [cm]", &lbs);
    bindParameter("la", &la, "Apical zone, length before shoot apical meristem [cm]", &las);
    bindParameter("ln", &ln, "Inter-lateral length, e.g. internode length [cm]", &lns);
    bindParameter("lnf", &lnf, "Type of inter-lateral distance (0 homogeneous, 1 linear inc, 2 linear dec, 3 exp inc, 4 exp dec)");
    bindParameter("lmax", &lmax, "Maximal stem length [cm]", &lmaxs);
    bindParameter("r", &r, "Initial growth rate [cm day-1]", &rs);
    bindParameter("a", &a, "Stem radius [cm]", &as);
    bindParameter("RotBeta", &rotBeta, "RevRotation of the stem");
    bindParameter("BetaDev", &betaDev, "RevRotation deviation");
    bindParameter("InitBeta", &initBeta, "Initial RevRotation");
    bindParameter("tropismT", &tropismT, "Type of stem tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)");
    bindParameter("tropismN", &tropismN, "Number of trials of stem tropism");
    bindParameter("tropismS", &tropismS, "Mean value of expected change of stem tropism [1/cm]");
    bindParameter("tropismAge", &tropismAge, "Age at which organ switch tropism", &tropismAges);
    bindParameter("theta", &theta, "Angle between stem and parent stem [rad]", &thetas);
    bindParameter("rlt", &rlt, "Stem life time [day]", &rlts);
    bindParameter("gf", &gf, "Growth function number [1]", &rlts);
    bindParameter("nodalGrowth", &nodalGrowth, "nodal growth function (sequential = 0, equal = 0)");
    bindParameter("delayNGStart", &delayNGStart, "delay between stem creation and start of nodal growth", &delayNGStarts);
    bindParameter("delayNGEnd", &delayNGEnd, "delay between stem creation and start of nodal growth", &delayNGEnds);
}

} // end namespace CPlantBox
