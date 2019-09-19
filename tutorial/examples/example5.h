/**
 * Example 5
 *
 * Hydrotropism ,
 * proof of concept, with a static soil water content
 *
 */

namespace CRootBox {

void example5()
{
    using namespace std;

    RootSystem rootsystem;

    string name = "Zea_mays_5_Leitner_2014";

    /*
     * Plant and root parameter from a file
     */
    rootsystem.openFile(name);

    // "manually" set tropism to hydrotropism
    for (int i=1; i<7; i++) {
        rootsystem.getRootRandomParameter(i)->tropismT = RootSystem::tt_hydro;
        rootsystem.getRootRandomParameter(i)->tropismN = 1; //N
        rootsystem.getRootRandomParameter(i)->tropismS = 0.4; //sigma
    }

    /**
     * Static soil property
     */
    SDF_PlantBox sideBox(10,20,50);
    SDF_RotateTranslate left(&sideBox, Vector3d(-4.99,0,0));
    SDF_RotateTranslate right(&sideBox, Vector3d(4.99,0,0));
    SDF_Union leftright(&left,&right);
    rootsystem.setGeometry(&leftright);  //for vizualisation

    double maxS = 0.7; // maximal saturation
    double minS = 0.1; // minimal saturation
    double slope = 20; // [cm] linear gradient between min and max
    SoilLookUpSDF soilprop(&left, maxS, minS, slope);

    // set the soil properties before calling initialize
    rootsystem.setSoil(&soilprop);

    /**
     * Initialize
     */
    rootsystem.initialize(); //it is important to call initialize() after setGeometry()

    /*
     * Simulate
     */
    double simtime = 30; // e.g. 30 or 60 days
    double dt = 1;
    int N = round(simtime/dt);
    for (int i=0; i<N; i++) {
        rootsystem.simulate(dt);
    }

    /**
     * Export results (as vtp)
     */
    rootsystem.write(name+".vtp");

    /**
     * Export container geometry as Paraview Python script (run file in Paraview by Tools->Python Shell, Run Script)
     */
    rootsystem.write(name+".py");

    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes\n";
}

} // end namespace CRootBox
