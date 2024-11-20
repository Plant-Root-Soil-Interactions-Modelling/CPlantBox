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
            return RootRandomParameter::toString(true);
        } else {
            return RootRandomParameter::toString(false);
        }
    }

    MycorrhizalRootRandomParameter::MycorrhizalRootRandomParameter(std::shared_ptr<Organism> plant) :RootRandomParameter(plant) {
        name = "undefined";
        organType = Organism::ot_root;// SOLLTE NICHT WAS ANDERES SEIN?
        subType = -1;
        f_tf = std::make_shared<Tropism>(plant);
        bindParameters();
    }
}
