#ifndef MYCORRHIZIALROOTPARAMETER_H_
#define MYCORRHIZIALROOTPARAMETER_H_

#include "rootparameter.h"
#include "organparameter.h"


namespace CPlantBox {


    class MycorrhizalRootSpecificParameter :public RootSpecificParameter{
        public:
        MycorrhizalRootSpecificParameter(): MycorrhizalRootSpecificParameter(-1, 0., 0., std::vector<double>(0), 0., 0., 0., 0.) { }

        MycorrhizalRootSpecificParameter(int type, double lb, 
        double la, const std::vector<double>& ln, double r, double a,
        double theta, double rlt, int infected = false, bool laterals = false): 
        RootSpecificParameter(type,lb,la,ln,r,a,theta,rlt,laterals), infected(infected) {};

        int infected;  ///< Status of infection with AMF


    };

    class MycorrhizalRootRandomParameter :public RootRandomParameter {
        public:

        MycorrhizalRootRandomParameter(std::shared_ptr<Organism> plant);
        
        virtual ~MycorrhizalRootRandomParameter() {};

        std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant) override;
        // realize
        std::string toString(bool verbose = true) const override;
        
        //void readXML(tinyxml2::XMLElement* element, bool verbose) override;


        void bindParameters() override;
        /*
            Internal AMF Infection Parameters
        */
        double p = 0.15;        ///< Probability of primary infection for dispersed inoculum [1/(cm day)]
        double minAge = 0;      ///< Minimal Infectious age of a root segment [day]
        double maxAge = 32;     ///< Maximal Infection age of a root segment [day]
        double vi = 0.13;       ///< rate of internal infection front [cm / day]
        double maxInfection = 1;    ///< Percentage of maximal infection
        double posX = 0;    ///< x Position of the localized inoculum
        double posY = 0;      ///< y Position of the localized inoculum
        double posZ = -3;       ///< z Position of the localized inoculum
        double infradius = 1;     ///< Radius of the localized inoculum
        // double nEntryP = 0; //< verbindung zu externen hyphen 
        int infected = 0;  ///< status of AMF infection
        // add parameter for localized infection here? like f_tf that can be modified from python binding

    };
}

#endif
