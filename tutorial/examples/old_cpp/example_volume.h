#include "RootSystem.h"
#include "SegmentAnalyser.h"

/**
 * Example 1
 *
 * 1) Opens plant and root parameters from a file
 * 2) Simulates root growth
 * 3) Outputs a VTP (for vizualisation in ParaView)
 *    In Paraview: use tubePLot.py script for fancy visualisation (Macro/Add new macro...), apply after opening file
 *
 *  Additionally, exports the line segments as .txt file to import into Matlab for postprocessing
 */
namespace CPlantBox {

void example_volume()
{
    auto rs = std::make_shared<RootSystem>();

    std::string path = "../../../modelparameter/rootsystem/";
    std::string name = "maize_p1_zero_std.xml"; // "maize_p1_zero_std", "maize_p2_zero_std", "maize_p3_zero_std"
    rs->readParameters(path+name);

    rs->initialize();

    /*
     * Simulate
     */
    double target_volume = 90; // cm^3
    double simtime = 365;  // maximal simulation time
    double dt = 0.5; // days

    double t = 0;
    double vol = 0;
    while ((vol<target_volume) && (t<simtime)) {
        // simulate
        rs->simulate(dt);
        t += dt;
        // calculate root system volume
        std::vector<double> vols = rs->getParameter("volume");
        vol = std::accumulate(vols.begin(), vols.end(), 0.);
    }
    std::cout << "\nfinished with " << vol << " cm^3 after "<< t << " days for " << name  << "\n\n";

    /*
     * Export final result (as vtp)
     */
    rs->write(name +".vtp");

    auto analysis = SegmentAnalyser(*rs);
    double l = analysis.getSummed("volume");
    std::cout << "Root system volume " << l << " cm^3 \n";
    std::cout << "Finished with a total of " << rs->getNumberOfNodes()<< " nodes\n";
}

} // end namespace CRootBox


