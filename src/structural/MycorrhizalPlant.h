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

    void initialize(bool verbose = true) override {initializeLB(verbose);};
    void initializeLB(bool verbose = true);

    void simulate(double dt, bool verbose) override;
    void simulateColonization(double dt, bool verbose);
    void simulatePrimaryColonization(double dt, bool verbose);
    void simulateSecondaryColonization(double dt, bool verbose);
    void simulateHyphalGrowth(double dt, bool verbose);
    void simulateAnastomosis(double dt, bool verbose);
    void simulateHyphae(double dt, bool verbose);

    virtual std::vector<int> getNodeColonizations(int ot) const; // returns Colonizations
    virtual std::vector<double> getNodeColonizationTime(int ot) const; // returns Colonization Time
    virtual std::vector<int> getAnastomosisPoints(int ot) const; // returns Anastomosis Points
    std::vector<int> getNodeTips(int ot) const;

    void initCallbacks() override;

    void changeGeometry(int ot,std::shared_ptr<SignedDistanceFunction> geom); // changes the geometry for a specific organ type, call before simulation

    int getNextHyphalTreeIndex() { hyphalTreeIndex++; return hyphalTreeIndex; }
	void turnOffSidePetriDish(double minnodex, double maxnodez, double maxnodey, double minnodey) ;

    std::shared_ptr<SDF_RootSystem> sdf; // direction from tip towards root base

    std::vector<int> localNodes;
    std::vector<std::shared_ptr<Hyphae>> localHyphae;
    int hyphalTreeIndex = -1;

    };

}

#endif // namespace


