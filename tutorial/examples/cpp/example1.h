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

void example1()
{
    using namespace std;

    RootSystem rootsystem;

    string name = "Anagallis_femina_Leitner_2010";

    /*
     * Open plant and root parameter from a file
     */
    rootsystem.openFile(name);
    // rootsystem.writeParameters(std::cout);

    /*
     * Initialize
     */
    rootsystem.initialize();

    /*
     * Simulate
     */
    double simtime = 60; // 20, 40, 60 days
    double dt = 60; // try other values here
    int N = round(simtime/dt);

    for (int i=0; i<N; i++) {
        rootsystem.simulate(dt);
    }
    cout << "fin\n";

    /*
     * Export final result (as vtp)
     */
    rootsystem.write("rootsystem.vtp");

    /*
     * Export segments in RSML format
     */
    rootsystem.write("rootsystem.rsml");

    /*
     * Export segments for Matlab analysis
     */
    SegmentAnalyser analysis(rootsystem);
    analysis.write("rootsystem.txt");

    /*
     * Export dgf file
     */
    SegmentAnalyser analysis_dgf(rootsystem);
    analysis_dgf.write(name+".dgf");

    /*
      Total length and surface
     */
    double l = analysis.getSummed("length");
    std::cout << "Root system length " << l << " cm \n";

    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes\n";

//	// Root System copy example
//	RootSystem rs = RootSystem();
//	string name = "Anagallis_femina_Leitner_2010";
//	rs.openFile(name);
//	rs.initialize();
//
//	cout << "Simulate... \n";
//
//	rs.simulate(20); // for a bit
//
//	auto rs2 = RootSystem(rs);
//
//	cout << "Nodes... \n";
//
//	auto nodes1 = rs.getNodes();
//	auto nodes2 = rs2.getNodes();
//
//	cout << "nodes 1: " << nodes1.size() << "\n";
//	cout << "nodes 2: " << nodes2.size() << "\n";
//
//	cout << "Simulate both \n";
//
//	rs.simulate(10);
//	rs2.simulate(10);
//
//	nodes1 = rs.getNodes();
//	nodes2 = rs2.getNodes();
//
//	cout << "New nodes... \n";
//	cout << "nodes 1: " << nodes1.size() << "\n";
//	cout << "nodes 2: " << nodes2.size() << "\n";
//

}

} // end namespace CRootBox
