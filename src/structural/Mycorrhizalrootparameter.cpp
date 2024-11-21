#include "Mycorrhizalrootparameter.h"

#include "Organism.h"

namespace CPlantBox {
    void MycorrhizalRootRandomParameter::bindParameters(){
        RootRandomParameter::bindParameters();
        bindParameter("p", &p, "Probability of primary infection for dispersed inoculum [1/(cm day)]");
        bindParameter("minAge", &minAge, "Minimal infectious age of a root segment [day]");
        bindParameter("maxAge", &maxAge, "Maximal infection age of a root segment [day]");
        bindParameter("vi", &vi, "Rate of internal infection front [cm / day]");
        bindParameter("maxInfection", &maxInfection, "Percentage of maximal infection");
    }

    // std::string MycorrhizalRootRandomParameter::toString(bool verbose) const {

    //     if (verbose) {
    //         return RootRandomParameter::toString(true);
    //     } else {
    //         return RootRandomParameter::toString(false);
    //     }
    // }

    MycorrhizalRootRandomParameter::MycorrhizalRootRandomParameter(std::shared_ptr<Organism> plant) :RootRandomParameter(plant) {
        bindParameters();
    }

    std::shared_ptr<OrganRandomParameter> MycorrhizalRootRandomParameter::copy(std::shared_ptr<Organism> p) {
        // std::cout << "MycorrhizalRootRandomParameter::copy\n" << std:.flush;
        auto r = std::make_shared<MycorrhizalRootRandomParameter>(*this);
        r->plant = p;
        r->bindParameters();
        // r->f_tf = f_tf->copy(p); // copy call back function classes
        // r->f_gf = f_gf->copy();
        // r->f_se = f_se; // for carbon limited grow we want the same reference
        // r->f_sa = f_sa->copy();
        // r->f_sbp = f_sbp->copy();
        return r;
    }
}
