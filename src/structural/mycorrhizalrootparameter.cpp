#include "mycorrhizalrootparameter.h"
#include "rootparameter.h"
#include "Organism.h"

namespace CPlantBox {

void MycorrhizalRootRandomParameter::bindParameters(){
    // std::cout << "MycorrhizalRootRandomParameter::bindParameters called" << std::endl;
    RootRandomParameter::bindParameters();
    bindParameter("lambda", &lambda, "Rate of primary infection for dispersed inoculum [1/(cm day)]");
    bindParameter("minAge", &minAge, "Minimal infectious age of a root segment [day]");
    bindParameter("maxAge", &maxAge, "Maximal infection age of a root segment [day]");
    bindParameter("vi", &vi, "Rate of internal infection [cm / day]");
    bindParameter("maxInfection", &maxInfection, "Percentage of maximal infection");

    // TODO
}

// std::shared_ptr<OrganSpecificParameter> MycorrhizalRootRandomParameter::realize() {
//     // std::c// std::shared_ptr<OrganSpecificParameter> MycorrhizalRootRandomParameter::realize() {
//     // std::cout << "MycorrhizalRootRandomParameter::realize called" << std::endl;
//     assert(dx > dxMin && "MycorrhizalRootRandomParameter::realize(): dxMin must be smaller than dx");
//     auto p = plant.lock();

//     double lambda_ = std::max(lambda + p->randn()*lambdas, 0.); // rate of primary infection for dispersed inoculum
//     double vi_ = std::max(vi + p->randn()*vis, 0.);  // speed of node to node infection
//     // std::cout<< "MycorrhizalRootRandomParameter::realize called" << std::endl;

//     // // Use base class realization to get all randomized root parameters
//     // auto baseroot = this -> RootRandomParameter::realize();
//     // std::cout << baseroot->toString() << std::endl;
//     // auto castedroot = std::dynamic_pointer_cast<RootSpecificParameter>(baseroot);

//     // // Use the realized parameters from the base class
//     // auto lb_ = castedroot->lb;
//     // auto la_ = castedroot->la;
//     // auto ln_ = castedroot->ln;
//     // auto r_ = castedroot->r;
//     // auto a_ = castedroot->a;
//     // auto theta_ = castedroot->theta;
//     // auto rlt_ = castedroot->rlt;
//     // auto hasLaterals = castedroot->laterals;

//     double lb_; //define the parameters outside of the if functions:
//     double la_;
//     std::vector<double> ln_; // stores the inter-distances
//     double nob_sd = p->randn()*nobs();
//     int nob_real = round(std::max(nob() + nob_sd, 0.)); // real maximal number of branching points
//     bool hasLaterals = (successorST.size()>0) && (nob_real>0);

//     if (!hasLaterals) { // no laterals
//         lb_ = 0;
//         la_ = std::max(lmax + p->randn()*lmaxs, 0.); // la, and lb is ignored
//         la_ = snap(la_);
//     } else { // laterals
//         lb_ = snap(std::max(lb + p->randn()*lbs, 0.)); // length of basal zone
//         la_ = snap(std::max(la + p->randn()*las, 0.)); // length of apical zone
//         double ln_mean = ln;
//         if(ln < dxMin && ln != 0) { // limit to minimum resolution
//             ln_mean = dxMin;
//         }
//         double nob1 = std::max((lmax-la_-lb_)/ln_mean+1, 0.); // use new la_, lb_ and ln_mean
//         int nob_ = std::min(std::max(round(nob1 + nob_sd), 0.), double(nob_real)); // maximal number of branches +1
//         int latMissing = nob_real - nob_;
//         assert((latMissing >= 0) && "RootRandomParameter::realize(): latMissing < 0");
//         int latExtraMean = floor(latMissing/nob_real); // mean number of extra laterals per branching point to keep correct number
//         int latExtra = latMissing - latExtraMean*(nob_);
//         for (int j = 0; j<latExtraMean; j++) { //at end of basal zone
//             ln_.push_back(0);
//         }
//         if (latExtra> 0) { //at end of basal zone
//             ln_.push_back(0);
//             latExtra--;
//         }
//         double sum_ln = nob_*ln_mean; // mean length of lateral zone

//         for (int i = 0; i<nob_-1; i++) { // create inter-root distances
//             double z = ((double)i+0.5)*ln_mean; // regular position along root lateral zone
//             double f = lnk*(z-sum_ln/2.); // evaluate slope lnk f(mid) = 0
//             double pf = (ln_mean + f) / ln_mean; // we scale lns by the change in percentage
//             double d = std::max(ln_mean + f + pf*p->randn()*lns, 1.e-5); // miminum is 1.e-5
//             d = snap(d);
//             ln_.push_back(d);
//             for (int j = 0; j<latExtraMean; j++) {
//                 ln_.push_back(0);
//             }
//             if (latExtra> 0) {
//                 ln_.push_back(0);
//                 latExtra--;
//             }
//         }
//     }

//     double r_ = std::max(r + p->randn()*rs, 0.); // initial elongation
//     double a_ = std::max(a + p->randn()*as, 0.01); // radius
//     double theta_ = std::max(theta + p->randn()*thetas, 0.); // initial elongation
//     double rlt_ = std::max(rlt + p->randn()*rlts, 0.); // root life time
//     // std::cout<< this->toString() << std::endl;

//     return std::make_shared<MycorrhizalRootSpecificParameter>(
//         subType, lb_, la_, ln_, r_, a_, theta_, rlt_, hasLaterals, p_, vi_);
// }

MycorrhizalRootRandomParameter::MycorrhizalRootRandomParameter(std::shared_ptr<Organism> plant) :RootRandomParameter(plant) {
    // std::cout << "MycorrhizalRootRandomParameter::MycorrhizalRootRandomParameter called" << std::endl;
    bindParameters();
}

std::shared_ptr<OrganRandomParameter> MycorrhizalRootRandomParameter::copy(std::shared_ptr<Organism> p) {
    // std::cout << "MycorrhizalRootRandomParameter::copy called" << std::endl;
    auto r = std::make_shared<MycorrhizalRootRandomParameter>(*this);
    r->plant = p;
    r->bindParameters();
    r->f_inf = f_inf;
    r->f_tf = f_tf->copy(p); // copy call back function classes
    r->f_gf = f_gf->copy();
    r->f_se = f_se; // for carbon limited grow we want the same reference
    r->f_sa = f_sa->copy();
    r->f_sbp = f_sbp->copy();
    return r;
}

std::string MycorrhizalRootRandomParameter::toString(bool verbose) const {
    // std::cout << "MycorrhizalRootRandomParameter::toString called" << std::endl;
    if (verbose) {
        return OrganRandomParameter::toString(true);
    } else {
        return OrganRandomParameter::toString(false);
    }
}

}
