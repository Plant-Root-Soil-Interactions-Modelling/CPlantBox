#ifndef MYCORRHIZAL_PLANT_H_
#define MYCORRHIZAL_PLANT_H_

#include <iostream>

#include "Plant.h"
#include "MycorrhizalRoot.h"

namespace CPlantBox {
    class MycorrhizalPlant :public Plant {
        MycorrhizalPlant(unsigned int seednum = 0.);
        virtual ~MycorrhizalPlant() {};

        //copy
        void initializeReader() override; ///< initializes XML reader
        //readParameters
        //initialize LB / DB
        //simulate ?
        //toString    
    };
}

#endif