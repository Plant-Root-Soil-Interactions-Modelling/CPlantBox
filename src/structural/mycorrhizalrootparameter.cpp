#include "mycorrhizalrootparameter.h"
#include "rootparameter.h"
#include "Organism.h"

namespace CPlantBox {

void MycorrhizalRootRandomParameter::bindParameters(){
    // std::cout << "MycorrhizalRootRandomParameter::bindParameters called" << std::endl;
    RootRandomParameter::bindParameters();
    bindParameter("p", &p, "Probability of primary infection for dispersed inoculum [1/(cm day)]");
    bindParameter("minAge", &minAge, "Minimal infectious age of a root segment [day]");
    bindParameter("maxAge", &maxAge, "Maximal infection age of a root segment [day]");
    bindParameter("vi", &vi, "Rate of internal infection [cm / day]");
    bindParameter("maxInfection", &maxInfection, "Percentage of maximal infection");
    bindParameter("infradius", &infradius, "Radius of the localized infection front");

    // TODO
}

// std::shared_ptr<OrganSpecificParameter> MycorrhizalRootRandomParameter::realize() {
//     // std::cout << "MycorrhizalRootRandomParameter::realize called" << std::endl;
//     assert(dx > dxMin && "MycorrhizalRootRandomParameter::realize(): dxMin must be smaller than dx");
//     auto pl = plant.lock();

//     double p_ = std::max(p + pl->randn()*ps, 0.); // rate of primary infection for dispersed inoculum
//     double vi_ = std::max(vi + pl->randn()*vis, 0.);  // speed of node to node infection
//     std::cout<< "MycorrhizalRootRandomParameter::realize called" << std::endl;

//     // Use base class realization to get all randomized root parameters
//     auto baseroot = this -> RootRandomParameter::realize();
//     std::cout << baseroot->toString() << std::endl;
//     auto castedroot = std::dynamic_pointer_cast<RootSpecificParameter>(baseroot);

//     // Use the realized parameters from the base class
//     auto lb_ = castedroot->lb;
//     auto la_ = castedroot->la;
//     auto ln_ = castedroot->ln;
//     auto r_ = castedroot->r;
//     auto a_ = castedroot->a;
//     auto theta_ = castedroot->theta;
//     auto rlt_ = castedroot->rlt;
//     auto hasLaterals = castedroot->laterals;

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
