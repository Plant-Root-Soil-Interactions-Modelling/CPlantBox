#ifndef MYCORRHIZIALROOTPARAMETER_H_
#define MYCORRHIZIALROOTPARAMETER_H_

#include "rootparameter.h"
#include "organparameter.h"


namespace CPlantBox {

    class Organsim;

    class MycorrhizalRootSpecificParameter : public RootSpecificParameter{
        public:
        MycorrhizalRootSpecificParameter(): MycorrhizalRootSpecificParameter(-1-1, 0., 0., std::vector<double>(0), 0., 0., 0., 0.,0.,0.,0.,0.) {};

        MycorrhizalRootSpecificParameter(int type, double lb, 
        double la, const std::vector<double>& ln, double r, double a,
        double theta, double rlt, bool laterals = false, double p,
        double minAge = 0, double maxAge = 32, double vi = 0.13): 
        RootSpecificParameter(type,lb,la,ln,r,a,theta,rlt,laterals), p(p),minAge(minAge),
        maxAge(maxAge), vi(vi) {};

        /*
         * AMF Infection Parameters
         */
        double p;        ///< Probability of primary infection for dispersed inoculum [1/(cm day)]
        double minAge;      ///< Minimal Infectious age of a root segment [day]
        double maxAge;     ///< Maximal Infection age of a root segment [day]
        double vi;       ///< rate of internal infection front [cm / day]
    };

    class MycorrhizalRootRandomParameter :public RootRandomParameter {
        public:

        MycorrhizalRootRandomParameter(std::shared_ptr<Organism> plant); // plant oder mycplant???
        virtual ~MycorrhizalRootRandomParameter() {};

        void bindParameters() override;
        /*
            Internal AMF Infection Parameters
        */
        double p = 0.15;        ///< Probability of primary infection for dispersed inoculum [1/(cm day)]
        double minAge = 0;      ///< Minimal Infectious age of a root segment [day]
        double maxAge = 32;     ///< Maximal Infection age of a root segment [day]
        double vi = 0.13;       ///< rate of internal infection front [cm / day]
        double maxInfection = 1;    ///< Percentage of maximal infection 

    };
}


#endif