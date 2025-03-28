#ifndef MYCORRHIZIALROOTPARAMETER_H_
#define MYCORRHIZIALROOTPARAMETER_H_

#include "rootparameter.h"
#include "organparameter.h"
#include "soil.h"


namespace CPlantBox {


    class MycorrhizalRootSpecificParameter :public RootSpecificParameter{
        public:
        MycorrhizalRootSpecificParameter(): MycorrhizalRootSpecificParameter(-1, 0., 0., std::vector<double>(0), 0., 0., 0., 0.) { }

        MycorrhizalRootSpecificParameter(int type, double lb, double la, const std::vector<double>& ln, double r, double a, double theta, double rlt, int infected = false, bool laterals = false):
            RootSpecificParameter(type,lb,la,ln,r,a,theta,rlt,laterals) {};

    };

    class MycorrhizalRootRandomParameter :public RootRandomParameter {
        public:

        MycorrhizalRootRandomParameter(std::shared_ptr<Organism> plant);

        virtual ~MycorrhizalRootRandomParameter() {};

        std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant) override;
        std::string toString(bool verbose = true) const override;

        //void readXML(tinyxml2::XMLElement* element, bool verbose) override;


        void bindParameters() override;
        /*
            Internal AMF Infection Parameters
        */
        double p = 0.15;        ///< Probability of primary infection for dispersed inoculum [1/(cm day)]
        double minAge = 0;      ///< Minimal Infectious age of a root segment [day]
        double maxAge = 32;     ///< Maximal Infection age of a root segment [day]
        double vi = 0.13;       ///< speed of internal infection [cm / day]
        double maxInfection = 1;    ///< Percentage of maximal infection
        double infradius = 1;     ///< Radius of the localized inoculum
        // double nEntryP = 0; //< verbindung zu externen hyphen

        std::shared_ptr<SoilLookUpSDF> f_inf = std::make_shared<SoilLookUpSDF>();

    };
}

#endif
