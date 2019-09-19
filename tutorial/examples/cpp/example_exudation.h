#include "../external/gauss_legendre/gauss_legendre.h"
#include "exudation.h"

#include <functional>

namespace CRootBox {


/**
 *
 */
void example_exudation()
{
    RootSystem rootsystem;
    std::string name = "Zea_mays_1_Leitner_2010";

    /*
     * Plant and root parameter from a file
     */
    rootsystem.openFile(name);
    rootsystem.writeParameters(std::cout);

    // THETA by hand
    rootsystem.getRootRandomParameter(2)->theta = M_PI/4.;
    rootsystem.getRootRandomParameter(2)->thetas = 0; // no std

    /*
     * Set geometry
     */
    //creates a square 27*27 cm containter with height 1.5 cm (used in parametrisation experiment)
    SDF_PlantBox box(900,900,900);
    rootsystem.setGeometry(&box);

    /*
     * Initialize
     */
    rootsystem.initialize();

    /*
     * Simulate
     */
    double simtime = 10.; // 20, 40, 60 days
    double dt = 1.; // try other values here
    int N = round(simtime/dt);
    for (int i=0; i<N; i++) {
        rootsystem.simulate(dt);
    }

    std::cout << "Number of roots " << rootsystem.getRoots().size() << "\n";

    /*
     * Export final result (as vtp)
     */
    rootsystem.write(name+".vtp");

//    /*
//     * Exudation
//     */
//    ExudationParameters params;
//
//    std::vector<double> allC = getExudateConcentration(rootsystem, params, 10, 10, 100, 30, 50);
//
//    /*
//     * write as txt file
//     */
//    std::ofstream fos;
//    fos.open(name+".txt");
//    for (size_t i=0; i<allc.size(); i++) {
//        fos << allc[i] << "\n";
//    }
//    fos.close();
//
//    std::cout << "fin \n";

}

} // end namespace CRootBox




