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


    MycorrhizalRootRandomParameter::MycorrhizalRootRandomParameter(std::shared_ptr<Organism> plant) :RootRandomParameter(plant) {
        bindParameters();
    }

    std::shared_ptr<OrganRandomParameter> MycorrhizalRootRandomParameter::copy(std::shared_ptr<Organism> p) {
        // std::cout << "MycorrhizalRootRandomParameter::copy\n" << std:.flush;
        auto r = std::make_shared<MycorrhizalRootRandomParameter>(*this);
        r->plant = p;
        r->bindParameters();
        return r;
    }
}