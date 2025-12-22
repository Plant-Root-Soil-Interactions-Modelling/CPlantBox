#ifndef MYCORRHIZAL_PLANT_H_
#define MYCORRHIZAL_PLANT_H_

#include "Plant.h"
#include "sdf_rs.h"


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

    virtual std::vector<int> getNodeInfections(int ot) const; // returns Infections
    virtual std::vector<double> getNodeInfectionTime(int ot) const; // returns Infection Time
    virtual std::vector<Vector3d> getAnastomosisPoints(int ot) const; // returns Anastomosis Points

    // void setInfectionSoil(std::shared_ptr<Soil> soil); //?? set a soil here
    void initCallbacks() override;

    int getNextHyphalTreeIndex() { hyphalTreeIndex++; return hyphalTreeIndex; }

    std::shared_ptr<SDF_RootSystem> sdf; // direction from tip towards root base

    std::vector<int> localNodes;
    std::vector<std::shared_ptr<Hyphae>> localHyphae;
    int hyphalTreeIndex = -1;
    // void updateAnastomosisTree(double dt);

    };

}

#endif // namespace


