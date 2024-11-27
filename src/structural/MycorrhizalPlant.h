#ifndef MYCORRHIZAL_PLANT_H_
#define MYCORRHIZAL_PLANT_H_

#include "Plant.h"
#include "MycorrhizalRoot.h"

namespace CPlantBox {
    class MycorrhizalPlant :public Plant {
        public:
        MycorrhizalPlant(unsigned int seednum = 0.);
        virtual ~MycorrhizalPlant() {};

        //copy
        void initializeReader() override; ///< initializes XML reader
        //readParameters neuer basetag?
        //void initializeLB(bool verbose = true);
        //void initializeDB(bool verbose = true);
        //simulate ?
        //toString    
    };
}

#endif