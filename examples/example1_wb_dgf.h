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

namespace CRootBox {

void example1_wb_dgf()
{
    using namespace std;
    RootSystem rootsystem;

    string name = "Anagallis_femina_Leitner_2010";

    /*
     * Plant and root parameter from a file
     */
    rootsystem.openFile(name);
    rootsystem.writeParameters(std::cout);

    /*
     * Set geometry
     */
    //creates a box
    SDF_PlantBox box(900,900,900);
    rootsystem.setGeometry(&box);
    /*
     * Initialize
     */
    rootsystem.initialize();

    /*
     * Simulate
     */
    // double simtime = 10; // 20, 40, 60 days
    // double dt = 1; // try other values here
    // int N = round(simtime/dt);

    // for (int i=0; i<N; i++) {
        rootsystem.simulate(30);
    // }

    /*
     * Export final result (as vtp)
     */
    rootsystem.write(name+".vtp");

    /*
     * Export segments in RSML format
     */
    rootsystem.write(name+".rsml");

    /*
     * Export dgf format
     */
    SegmentAnalyser analysis(rootsystem);
    analysis.write(name+".dgf");

    /*
      Total length and surface
     */
    double l = analysis.getSummed("length");
    std::cout << "Visible Length " << l << " cm \n";


    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes\n";
    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes\n";

}

} // end namespace CRootBox
