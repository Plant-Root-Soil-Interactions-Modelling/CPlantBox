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
using namespace std;

void example1()
{
    Plant plant;

    string name = "CPlantBox_test_Xiaoran";

    /*
     * Open plant and root parameter from a file
     */
    plant.openFile(name);
    plant.writeAlltoXML(name);
    plant.openXML(name);
//     plant.writeParameters(std::cout);

    /*
     * Initialize
     */
    plant.initialize();
    std::cout << "finished initialize\n";
    /*
     * Simulate
     */
    double simtime = 60; // 20, 40, 60 days
    double dt = 60; // try other values here
    int N = round(simtime/dt);

    for (int i=0; i<N; i++) {
        plant.simulate(dt);
    }
    cout << "fin (with " << plant.getNumberOfNodes() << " nodes) \n";
    auto o_ = plant.getOrgans(Organ::ot_organ);
    cout << o_.size() << " organs \n";

    /*
     * Export final result (as vtp)
     */
    auto t1 = std::chrono::high_resolution_clock::now();
    plant.write("plant_write_rootsystem.vtp");
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "test function took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";



    /*
     * Export segments for Matlab analysis
     */
    SegmentAnalyser analysis(plant);
    analysis.write("rootsystem.txt");
    analysis.write("analysis_write_rootsystem.vtp");

    /*
     * Export dgf file
     */
    // analysis.write("rootsystem.dgf");

    /*
      Total length and surface
     */
//    double l = analysis.getSummed(Plant::st_length);
//    std::cout << "Root system length " << l << " cm \n";

    cout << "Finished with a total of " << plant.getNumberOfNodes()<< " nodes\n";

}
