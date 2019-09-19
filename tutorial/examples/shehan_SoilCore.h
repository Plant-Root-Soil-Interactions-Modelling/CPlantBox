/**
 * Shehans experimental set-up
 *
 * (A) Virtual soil core analysis
 *
 * Simulates the growing root systems (6*37) in the field
 * uses the soil coring method to analyse the root systems.
 *
 */

#include <iostream>
#include <iomanip>

namespace CRootBox {

/**
 * Creates a field of multiply root system (6*37)
 *
 * @param name     filename of the parameter file
 * @param geom     confining geometry
 */
std::vector<RootSystem*> initializeRootSystems(const std::string& name, SignedDistanceFunction* geom = new SignedDistanceFunction())
{
    using namespace std;
    auto gen = mt19937(chrono::system_clock::now().time_since_epoch().count());
    auto UID = std::uniform_int_distribution<unsigned int>(); // random stuff, does it work now?
    int M=6;
    int N=37;
    double dist1 = 18; // [cm] (M-1)*18 = 90
    double dist2 = 3; // [cm] (N-1)*3 = 108
    vector<RootSystem*> allRS;
    for (int i=0; i<M; i++) {
        for (int j=0;j<N; j++) {
            RootSystem* rs = new RootSystem();
            allRS.push_back(rs);
            rs->openFile(name);
            rs->setGeometry(geom);
            rs->getRootRandomParameter(4)->theta = 80./180.*M_PI; // fix insertion angle of the basal roots
            rs->getRootSystemParameter()->seedPos = Vector3d(dist1*i,dist2*j,-3); // set position of seed [cm]
            unsigned int s = UID(gen);
            rs->setSeed(s); // randomly select a seed
            rs->initialize();
        }
    }
    return allRS;
}

/**
 * Simulates all root systems and stores the results in one analyser class per time step
 *
 * @param times     simulation times we want to simulate (starting with 0)
 * @param allRS     root systems to simulate
 */
void simulateRS(std::vector<double> times, const std::vector<RootSystem*>& allRS)
{
    using namespace std;
    int N_ = allRS.size();
    for (size_t i=0; i<times.size()-1; i++) {
        std::cout << "\nTIME " << times.at(i+1) <<"\n\n";
        double dt = times.at(i+1)-times.at(i);
        for (int r=0; r<N_; r++) { // simulate all
            allRS.at(r)->simulate(dt,true);
        }
    }
}

/**
 * Retrieves a the state of the root systems at a certain time in a SegmentAnalyser object
 */
SegmentAnalyser getResult(const std::vector<RootSystem*>& allRS, double time) {
  using namespace std;
  SegmentAnalyser a = SegmentAnalyser();
  for (const auto& rs : allRS) { // merge all into one analyser object
      //cout << "Analyser " << a.segments.size()<< ", " << a.nodes.size() << "\n";
      //cout << "Root system: " << rs->getNumberOfNodes() <<"\n";
      auto news = SegmentAnalyser(*rs);
      news.filter("creation_time",0,time);
      news.pack(); // delete unused nodes
      a.addSegments(news);
  }
  return a;
}

/**
 * The geometry of the 15 soil cores with fixed location in the field
 *
 * @param r     radius of the soil core
 * @param h     height of the soil core
 */
std::vector<SignedDistanceFunction*> soilCores(double r,double h)
{
    using namespace std;
    SDF_PlantContainer* core =  new SDF_PlantContainer(r,r,h,false);
    vector<double> x = {27, 27, 27, 27, 27, 45, 45, 45, 45, 45, 63, 63, 63, 63, 63};
    vector<double> y = {30, 42, 54, 66, 78, 30, 42, 54, 66, 78, 30, 42, 54, 66, 78};
    vector<SignedDistanceFunction*> cores_;
    for (size_t i=0; i<x.size(); i++) {
        cores_.push_back(new SDF_RotateTranslate(core, 0, SDF_RotateTranslate::xaxis, Vector3d(x[i],y[i],0)));
    }
    // SDF_Union* cores = new SDF_Union(cores_);
    return cores_;
}



/**
 * Soil core analysis in the virtual field experiment
 */
void shehan_SoilCore(const std::string& name = "wheat", bool exportVTP = false)
{
    using namespace std;
    /*
     * Initialize
     */
    auto allRS = initializeRootSystems(name);

    /*
     * Simulate
     */
    vector<double> times = {0, 30, 60, 90};
    simulateRS(times, allRS);

    /*
     * Analysis: vertical distribution within the soil cores
     */
    double r = 2.1; // radius of the soil core (cm)
    double h = 160; // depth of the soil core (cm)
    double dz = 5; // width of layer (cm)
    vector<SignedDistanceFunction*> coregeometry = soilCores(r,h);

    vector<vector<double>> finalmatrix(15*(times.size()-1));
    for (size_t i=1; i<times.size(); i++) {
        for (size_t j=0; j<15; j++) {

            std::cout << "\nANALYSE TIME " << times.at(i) <<"\n\n";
            SegmentAnalyser coreanalyser = getResult(allRS, times.at(i));
            std::cout << "crop\n";
            coreanalyser.crop(coregeometry[j]); // throw segments away
            std::cout << "pack\n";
            coreanalyser.pack(); // throw unused nodes away
            std::cout << "get distribution\n";
            vector<double> tl = coreanalyser.distribution("length",0,h,round(h/dz),true);  // vertical distribution
            for (double& d :tl) {
                d = d/(15*r*r*M_PI*dz); // 15 soil cores
            }
            finalmatrix.at((i-1)*15+j) = tl; // store length density in the finalmatrix
            if (exportVTP) {
                string vtpname = name + "_core_cropped"+ std::to_string(i)+".vtp"; // export cropped segments for vizualisaten
                coreanalyser.write(vtpname);
            }

        }
    }

    /*
     * Export the final matrix
     */
    int n = finalmatrix[0].size(); // h/dz = 32
    int m = finalmatrix.size();
    string matname = name+"_core_matrix.txt";
    std::ofstream fos;
    fos.open(matname.c_str());
    for (int j=0; j<n; j++) {
        for (int i=0; i<m; i++) {
            std::cout << std::fixed << std::setprecision(4)<< finalmatrix[i][j] << "\t";
            fos << std::fixed << std::setprecision(4)<< finalmatrix[i][j] << "\t";
        }
        std::cout << "\n";
        fos << "\n";
    }
    fos.close();

//	/*
//	 * Export rootsystems (around 4.5GB)
//	 */
//	SegmentAnalyser analyser;
//	for (const auto& rs : allRS) {
//	    analyser.addSegments(*rs);
//	}
//	string vtpname = name +".vtp";
//	analyser.write(vtpname);


    /**
     * Export core geometry
     */
    if (exportVTP) {
        string gname = name + "_core.py";
        SDF_Union* cores = new SDF_Union(coregeometry);
        allRS[0]->setGeometry(cores); // just for writing
        allRS[0]->write(gname);
    }

}

} // end namespace CRootBox
