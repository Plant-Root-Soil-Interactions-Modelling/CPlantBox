#include "Mycorrhizalrootparameter.h"

#include "Organism.h"

namespace CPlantBox {
    void MycorrhizalRootRandomParameter::bindParameters(){
        RootRandomParameter::bindParameters();
        bindParameter("p", &p, "Probability of primary infection for dispersed inoculum [1/(cm day)]");
        bindParameter("minAge", &minAge, "Minimal Infectious age of a root segment [day]");
        bindParameter("maxAge", &maxAge, "Maximal Infection age of a root segment [day]");
        bindParameter("vi", &vi, "rate of internal infection front [cm / day]");
        bindParameter("maxInfection", &maxInfection, "Percentage of maximal infection");
    }

    std::string MycorrhizalRootRandomParameter::toString(bool verbose) const {

        if (verbose) {
            return OrganRandomParameter::toString(true);
        } else {
            return OrganRandomParameter::toString(false);
        }
    }
}
