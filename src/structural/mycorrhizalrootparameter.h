#ifndef MYCORRHIZIALROOTPARAMETER_H_
#define MYCORRHIZIALROOTPARAMETER_H_

#include "rootparameter.h"
#include "organparameter.h"
#include "soil.h"


namespace CPlantBox {


    class MycorrhizalRootSpecificParameter :public RootSpecificParameter{
        public:
        MycorrhizalRootSpecificParameter(): MycorrhizalRootSpecificParameter(-1, 0., 0., std::vector<double>(0), 0., 0.,0.15, 0.13, false, 0., 0.) { }

        MycorrhizalRootSpecificParameter(int type, double lb, double la, const std::vector<double>& ln, double r, double a, double theta, double rlt, 
				bool laterals, double lmbd, double vi, int infected = false):
            RootSpecificParameter(type,lb,la,ln,r,a,theta,rlt,laterals), lmbd(lmbd), vi(vi) {};

        double lmbd = 0.15;        ///< Rate of primary infection for dispersed inoculum [1/(cm day)]
        double vi = 0.13;          ///< speed of internal infection [cm / day]

    };

    class MycorrhizalRootRandomParameter :public RootRandomParameter {
        public:

        MycorrhizalRootRandomParameter(std::shared_ptr<Organism> plant);

        virtual ~MycorrhizalRootRandomParameter() {};

        std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant) override;
        std::shared_ptr<OrganSpecificParameter> realize() override; ///< Creates a specific root from the root parameter set
        std::string toString(bool verbose = true) const override;

        void bindParameters() override;
        /*
            Internal AMF Infection Parameters
        */
        double lmbd = 0.15;        ///< Rate of primary infection for dispersed inoculum [1/(cm day)]
        double lmbds = 0.;         ///< Standard deviation of primary infection [1/(cm day)]
        double minAge = 0;         ///< Minimal Infectious age of a root segment [day]
        double maxAge = 32;        ///< Maximal Infection age of a root segment [day]
        double vi = 0.13;          ///< speed of internal infection [cm / day]
        double vis = 0.;           ///< Standard deviation of speed of internal infection [cm / day]
        double maxInfection = 1;   ///< Percentage of maximal infection
        
        double hyphalDelay = 0.;      ///< Dummy delay time before hyphal emergence after infection [day] there is none.
        double highresolution = 1; ///< If true, a hypha is created at every infected node, otherwise hyphae are created based on hyphalEmergenceDensity
        double dx_inf = 0.1;       ///< Segment length for infected root segments only in high resolution case [cm]

        double hyphalEmergenceDensity = 0; //< [1 / cm]

        std::shared_ptr<SoilLookUp> f_inf = std::make_shared<SoilLookUp>();    
		void readXML(tinyxml2::XMLElement* element, bool verbose) override; ///< reads a single sub type organ parameter set

    };
}

#endif
