#include <sstream>
#include <fstream>

/**
 * Example for Shehan
 *
 * Scale elongation rate according to some 1d tabular data
 * including carbon limitation
 *
 * see shehan_ScaleElongation.h for definition of class ScaleElongation and function readCSV
 *
 * proof of concept, please adjust
 */

namespace CRootBox {

/**
 * The maximal possible total root system length increment per day (cm/day) at time t
 * (based on carbon, or leaf area)
 */
double length_increment(double t, Grid1D& grid) {
    double lai = grid.getValue(Vector3d(0., 0., t), nullptr); // leaf area index
    double li = 2 * lai; // TODO calculate increment based on lai with some magic formula
    return li; //cm per day
}

/**
 * Example how to implement varying elongation rates
 */
void shehan_ScaleElongation_CL()
{
    using namespace std;

    RootSystem rootsystem;

    string name = "Zea_mays_5_Leitner_2014";

    /*
     * Plant and root parameter from a file
     */
    rootsystem.openFile(name);

    /*
     * Create the tables for water content and temperature and the scale elongation function
     */
    // grid coordinates
    int n = 7;
    std::vector<double> z_{0., 15., 25., 50., 70., 100., 140.};

    // LAI data
    std::cout << "reading lai data ...\n";
    auto field_lai = readCSV("/home/daniel/workspace/CRootBox/examples/lai.csv", ',', 0, 0);
    std::vector<double> time = field_lai.at(0); // n time data points (eg.; 0, 1, 2, 3 days)
    std::vector<double> lai_data = field_lai.at(1); // n-1 lai data points (value 1 is between day 0-1, value 2, between day 1-2, and so on)
    lai_data.pop_back(); // the last value is thrown away (i.e. we give n lai_data points in the file, but only use n-1)
    Grid1D lai_grid = Grid1D(time.size(), time, lai_data);

    // temperature data
    std::cout << "reading temperature data ...\n";
    auto field_temp = readCSV("/home/daniel/workspace/CRootBox/examples/TEMP.csv", ';', 2, 1);
    std::vector<double> temp = field_temp.at(0);

    // water content data
    std::cout << "reading water content data ...";
    auto field_wc = readCSV("/home/daniel/workspace/CRootBox/examples/WC.csv", ',', 2, 1);
    std::vector<double> wc = field_wc.at(0);

    // create scale elongation function
    Grid1D water_content = Grid1D(n, z_, wc);
    Grid1D temperature = Grid1D(n, z_, temp );
    ScaleElongation se = ScaleElongation(&water_content, &temperature);

    // create carbon scaling (on top)
    ProportionalElongation pe = ProportionalElongation();
    pe.setBaseLookUp(&se);

    // "manually" set the scale elongation function
    for (int i = 1; i < 7; i++) {
        rootsystem.getRootRandomParameter(i)->f_se = &pe;
    }

    /**
     * Initialize
     */
    rootsystem.initialize(); //it is important to call initialize() after setGeometry()

    /*
     * Simulate
     */
    double simtime = 1; // 7*30; // days
    double dt = 0.5 * 1./24.;  // dt is half an hour
    size_t N = round(simtime/dt);

    assert(N<=field_temp.size()); // check if enough data are available

    string matname = "rld.csv";
    std::ofstream fos;
    fos.open(matname.c_str());
    double h = 160; // depth of the soil core (cm)
    double dz = 5; // width of layer (cm)


    for (size_t i=0; i<N; i++) {

        rootsystem.simulate(dt, length_increment(i * dt, lai_grid), &pe, false);

        // update field data:
        temperature.data =  field_temp.at(i);
        water_content.data = field_wc.at(i);

        // compute and write RLD
        SegmentAnalyser ana = SegmentAnalyser(rootsystem);
        vector<double> tl = ana.distribution("length", 0, h, round(h / dz), true);  // vertical distribution
        for (size_t i = 0; i < tl.size(); i++) {
            fos << tl[i] << ", ";
        }
        fos << "\n";

    }

    fos.close();

    /**
     * Export results (as vtp)
     */
    rootsystem.write(name+".vtp");

    double tl = rootsystem.getSummed("length");
    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes, " << tl << " cm total length \n";
}

} // end namespace CRootBox
