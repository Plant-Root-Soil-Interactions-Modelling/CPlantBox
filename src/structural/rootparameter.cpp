// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "rootparameter.h"

#include "Organism.h"

#include <cmath>
#include <iostream>
#include <chrono>
#include <assert.h>
#include <numeric>

namespace CPlantBox {

/**
 * @return Mean maximal root length of this root type
 */
double RootSpecificParameter::getK() const {
    double l = std::accumulate(ln.begin(), ln.end(), 0.);
    return l+la+lb;
}

/**
 * @copydoc OrganParameter::toString()
 */
std::string RootSpecificParameter::toString() const
{
    std::stringstream str;
    str << "subType\t" << subType << std::endl;
    str << "lb\t" << lb << std::endl << "la\t" << la << std::endl;
    str << "nob\t" << nob() << std::endl << "r\t" << r << std::endl << "a\t" << a << std::endl;
    str << "theta\t" << theta << std::endl << "rlt\t" << rlt << std::endl;
    str << "ln\t";
    for (int i=0; i<ln.size(); i++) {
        str << ln[i] << " ";
    }
    str << std::endl;
    return str.str();
}

/**
 * Default constructor sets up hashmaps for class introspection
 */
RootRandomParameter::RootRandomParameter(std::shared_ptr<Organism> plant) :OrganRandomParameter(plant)
{
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
std::shared_ptr<OrganRandomParameter> RootRandomParameter::copy(std::shared_ptr<Organism> p)
{
    // std::cout << "RootRandomParameter::copy\n"<< std::flush;
    auto r = std::make_shared<RootRandomParameter>(*this); // copy constructor breaks class introspection
    r->plant = p;
    r->bindParameters(); // fix class introspection
    r->f_tf = f_tf->copy(p); // copy call back function classes
    r->f_gf = f_gf->copy();
    r->f_se = f_se->copy();
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
std::shared_ptr<OrganSpecificParameter> RootRandomParameter::realize()
{
    assert(dx > dxMin && "RootRandomParameter::realize(): dxMin must be smaller than dx");
    auto p = plant.lock();
    double lb_; //define the parameters outside of the if functions:
    double la_;
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

    double r_ = std::max(r + p->randn()*rs, 0.); // initial elongation
    double a_ = std::max(a + p->randn()*as, 0.001); // radius
    double theta_ = std::max(theta + p->randn()*thetas, 0.); // initial elongation
    double rlt_ = std::max(rlt + p->randn()*rlts, 0.); // root life time

    return std::make_shared<RootSpecificParameter>(subType,lb_,la_,ln_,r_,a_,theta_,rlt_, hasLaterals);
}

/**
 * snaps to the grid, to make it compatible with dx() and dxMin()
 */
double RootRandomParameter::snap(double x)
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



/**
 * The standard deviation of number of branches is calculated form root maximal length lmax, and
 * its standard deviation lmaxs
 *
 * todo I have no idea why this holds...
 */
double RootRandomParameter::nobs() const
{
    double nobs = (lmaxs/lmax - lns/ln)*lmax/ln; // error propagation
    if (la>0) {
        nobs -= (las/la - lns/ln)*la/ln;
    }
    if (lb>0) {
        nobs -= (lbs/lb - lns/ln)*lb/ln;
    }
    return std::max(nobs,0.);
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
void RootRandomParameter::readXML(tinyxml2::XMLElement* element, bool verbose)
{
    OrganRandomParameter::readXML(element, verbose);
}


/**
 * CPlantBox parameter reader (DEPRICATED)
 */
void RootRandomParameter::read(std::istream & cin)
{
    std::cout << "RootRandomParameter::read is deprecated, use readXML instead \n";
    char ch[256]; // dummy
    cin.getline(ch,256);
    std::string s; // dummy
    cin >> s >> subType >> s >> name >> s >> lb >> lbs >> s >> la >> las >> s >> ln >> lns >> s >> lmax >> lmaxs;
    cin >> s >> r >> rs >> s >> a >> as >> s >> colorR >> colorG >> colorB >> s >> tropismT >> tropismN >> tropismS >> s >> dx;
    int n;
    cin  >> s >> n;
    successorST.clear();
    successorOT.clear();
    successorWhere.clear();
    successorNo.clear();
    int is;
    for (int i=0; i<n; i++) {
        cin >> is;
        successorST.push_back(std::vector<int>(is));
        successorOT.push_back(std::vector<int>(2));
        successorNo.push_back(1);
        successorWhere.push_back(std::vector<double>(0,0));
    }
    cin >> s >> n;
    successorP.clear();
    double ds;
    for (int i=0; i<n; i++) {
        cin >> ds;
        successorP.push_back(std::vector<double>(ds));
    }
    cin >> s >> theta >> thetas >> s >> rlt >> rlts >> s >> gf >> s;
}

/**
 * CPlantBox parameter write (DEPRICATED)
 */
void RootRandomParameter::write(std::ostream & cout) const {
    std::cout << "RootRandomParameter::write is deprecated, use writeXML instead \n";
    cout << "# Root type parameter for " << name << "\n";
    cout << "type\t" << subType << "\n" << "name\t" << name << "\n" << "lb\t"<< lb <<"\t"<< lbs << "\n" << "la\t"<< la <<"\t"<< las << "\n"
        << "ln\t" << ln << "\t" << lns << "\n" << "lmax\t"<< lmax <<"\t"<< lmaxs << "\n" << "r\t"<< r <<"\t"<< rs << "\n" <<
        "a\t" << a << "\t" << as << "\n" << "color\t"<< colorR <<"\t"<< colorG << "\t" << colorB << "\n"
        << "tropism\t"<< tropismT <<"\t"<< tropismN << "\t" << tropismS << "\n" << "dx\t" << dx << "\n" << "successor\t" << successorST.size() << "\t";
    for (size_t i=0; i<successorST.size(); i++) {
        for (int j=0; i<successorST.at(i).size(); j++) {
				cout << successorST.at(i).at(j) << "\t";
		}
    }
    cout << "\n" << "successorP\t" << successorP.size() <<"\t";
    for (size_t i=0; i<successorP.size(); i++) {
        for (int j=0; i<successorP.at(i).size(); j++) {
				cout << successorP.at(i).at(j) << "\t";
		}
    }
    cout << "\n" << "theta\t" << theta << "\t" << thetas << "\n" << "rlt\t" << rlt << "\t" << rlts << "\n" << "gf\t" << gf << "\n";
}

/**
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void RootRandomParameter::bindParameters()
{
    OrganRandomParameter::bindParameters();
    bindParameter("lb", &lb, "Basal zone [cm]", &lbs);
    bindParameter("la", &la, "Apical zone [cm]", &las);
    bindParameter("ln", &ln, "Inter-lateral distance [cm]", &lns);
    bindParameter("lmax", &lmax, "Maximal root length [cm]", &lmaxs);
    bindParameter("r", &r, "Initial growth rate [cm day-1]", &rs);
    bindParameter("colorR", &colorR, "Root color, red component [0.-1.]");
    bindParameter("colorG", &colorG, "Root color, green component [0.-1.]");
    bindParameter("colorB", &colorB, "Root color, blue component [0.-1.]");
    bindParameter("tropismT", &tropismT, "Type of root tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)");
    bindParameter("tropismN", &tropismN, "Number of trials of root tropism");
    bindParameter("tropismS", &tropismS, "Mean value of expected change of root tropism [1/cm]");
    bindParameter("theta", &theta, "Angle between root and parent root [rad]", &thetas);
    bindParameter("rlt", &rlt, "Root life time [day]", &rlts);
    bindParameter("gf", &gf, "Growth function number [1]", &rlts);
    // NEW
    bindParameter("lnk", &lnk, "Slope of inter-lateral distances [1]");
    bindParameter("ldelay", &ldelay, "Lateral root emergence delay [day]", &ldelays);
}

} // end namespace CPlantBox
