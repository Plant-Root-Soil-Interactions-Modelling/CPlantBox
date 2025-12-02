#ifndef MYCORRHIZAL_PLANT_H_
#define MYCORRHIZAL_PLANT_H_

#include "Plant.h"
#include "MycorrhizalRoot.h"
// #include "soil.h"
#include "sdf_rs.h" 
#include "aabbcc/AABB.h"
// #include "sdf.h"
//#include "Hyphae.h"

namespace CPlantBox {
    class MycorrhizalPlant :public Plant {
        public:
        MycorrhizalPlant(unsigned int seednum = 0.);
        virtual ~MycorrhizalPlant() {};

        //copy
        void initializeReader() override; ///< initializes XML reader
        void readParameters(std::string name, std::string basetag = "plant", bool fromFile = true, bool verbose = true) override {this -> initializeReader(); Organism::readParameters(name, basetag, fromFile, verbose);};

        // std::shared_ptr<Organ> createRoot(std::shared_ptr<Organism> plant, int type, double delay);
        void initialize(bool verbose = true) override {initializeLB(verbose);};
        void initializeLB(bool verbose = true);

        void simulate(double dt, bool verbose) override;
        void simulatePrimaryInfection(double dt, bool verbose);
        void simulateSecondaryInfection(double dt, bool verbose);
        void simulateHyphalGrowth(double dt);
        void simulateAnastomosis();
        void simulateAnastomosisTree();

        virtual std::vector<int> getNodeInfections(int ot) const; // returns Infections
        virtual std::vector<double> getNodeInfectionTime(int ot) const; // returns Infection Time

        // void setInfectionSoil(std::shared_ptr<Soil> soil); //?? set a soil here
        void initCallbacks() override;
        void addTree(); // AABB tree

        std::vector<double> stopTime; // time when root stopped growing, 0 if it has not
        std::vector<Vector3d> tips;
        std::vector<SDF_RootSystem> sdfs; // direction from tip towards root base

        aabb::Tree tree = aabb::Tree(); // aabb tree for anastomosis
        void buildAnastomosisTree();
        void updateAnastomosisTree(double dt);
        unsigned int getDistTree(unsigned int p,Vector3d  tip, double dist);

    };
}

#endif
