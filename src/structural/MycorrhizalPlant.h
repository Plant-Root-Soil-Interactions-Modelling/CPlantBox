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
        void readParameters(std::string name, std::string basetag = "plant", bool fromFile = true, bool verbose = true) override {this -> initializeReader(); Organism::readParameters(name, basetag, fromFile, verbose);};
        //readParameters neuer basetag?
        //void initializeLB(bool verbose = true);
        //void initializeDB(bool verbose = true);
        // simulate
        // für visualisierung was rausholen
        // plant.getSegmentCTs()  als vorlage für infektion auch zum fehler suchen !!! in den segmentanalyzer dazu tun
        // getSegmentsInfection()
        // getNodeInfection()
        //toString    
    };
}

#endif