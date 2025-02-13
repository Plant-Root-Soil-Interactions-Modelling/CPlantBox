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
        void readParameters(std::string name, std::string basetag = "plant", bool fromFile = true, bool verbose = true) override {std::cout<< "MycorrhizalPlant::readParameter called"<<std::endl; this -> initializeReader(); Organism::readParameters(name, basetag, fromFile, verbose);};
        
        std::shared_ptr<Organ> createRoot(std::shared_ptr<Organism> plant, int type, double delay);
        void initialize(bool verbose = true) override {std::cout<< "MycorrhizalPlant::initialize called" <<std::endl; initializeLB(verbose);};
        void initializeLB(bool verbose = true);

        // virtual std::vector<int> getNodeInfections(int ot) const;
        // virtual std::vector<int> getSegmentInfections(int ot) const;

        // für visualisierung was rausholen
        // plant.getSegmentCTs()  als vorlage für infektion auch zum fehler suchen !!! in den segmentanalyzer dazu tun
        //toString    
    };
}

#endif