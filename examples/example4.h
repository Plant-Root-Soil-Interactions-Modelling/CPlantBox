/**
 * Example 4
 *
 * More complex geometries
 *
 * 1. a split pot experiment
 * 2. rhizotubes as obstacles (from Example_Rhizotubes.m)
 *
 * Additionally, exports the confining geometry as paraview pyhton script
 * (run file in Paraview by Tools->Python Shell, Run Script)
 */
using namespace std;

void example4()
{
    RootSystem rootsystem;

    string name = "Zea_mays_5_Leitner_2014";

    /*
     * Plant and root parameter from a file
     */
    rootsystem.openFile(name);

    /*
     * 1. A split pot experiment
     */
    SDF_PlantBox topBox(22,20,5);
    SDF_PlantBox sideBox(10,20,35);

    SDF_RotateTranslate left(&sideBox, Vector3d(-6,0,-5));
    SDF_RotateTranslate right(&sideBox, Vector3d(6,0,-5));

    vector<SignedDistanceFunction*> box_;
    box_.push_back(&topBox);
    box_.push_back(&left);
    box_.push_back(&right);
    SDF_Union splitBox(box_);

    /*
     * 2. Rhizotubes as obstacles
     */
    // Box
    double boxX = 96.;
    double boxY = 126.;
    double boxZ = 130.;
    SDF_PlantBox box(boxX,boxY,boxZ);

    // A single Rhizotube
    double r = 2.*3.2; // cm
    SDF_PlantContainer rhizotube(r,r,96.,false);
    SDF_RotateTranslate rhizoX(&rhizotube, 90., SDF_RotateTranslate::yaxis, Vector3d(boxX/2.,0.,0));

    // The experimental setting
    vector<SignedDistanceFunction*> rhizotubes_;
    const int tubeN=6;
    double y_[tubeN] = { -30, -18, -6, 6, 18, 30 };
    double z_[tubeN]= { -10, -20, -40, -60, -80, -120 };
    for (int i=0; i<tubeN; i++) {
        rhizotubes_.push_back(new SDF_RotateTranslate(&rhizoX, 0, SDF_RotateTranslate::xaxis, Vector3d(0,y_[i],z_[i])));
    }

    // Final geometry
    SDF_Union rhizotubes(rhizotubes_);
    SDF_Difference rhizoTubes(&box, &rhizotubes);

    rootsystem.setGeometry(&rhizoTubes); // &splitBox, or &rhizoTubes

    /**
     * Initialize
     */
    rootsystem.initialize(); //it is important to call initialize() after setGeometry()

    /*
     * Simulate
     */
    double simtime = 100; // e.g. 30 or 60 days
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
