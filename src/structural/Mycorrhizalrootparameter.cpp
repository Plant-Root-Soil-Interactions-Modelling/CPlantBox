#include "Mycorrhizalrootparameter.h"

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
    bindParameter("nEntryP", &nEntryP, "Number of entry points for external hyphae");
}


MycorrhizalRootRandomParameter::MycorrhizalRootRandomParameter(std::shared_ptr<Organism> plant) :RootRandomParameter(plant) {
    // std::cout << "MycorrhizalRootRandomParameter::MycorrhizalRootRandomParameter called" << std::endl;
    bindParameters();
}

std::shared_ptr<OrganRandomParameter> MycorrhizalRootRandomParameter::copy(std::shared_ptr<Organism> p) {
    // std::cout << "MycorrhizalRootRandomParameter::copy called" << std::endl;
    auto r = std::make_shared<MycorrhizalRootRandomParameter>(*this);
    r->plant = p;
    r->bindParameters();
    r->f_inf = f_inf; // why no ->copy() here?
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
