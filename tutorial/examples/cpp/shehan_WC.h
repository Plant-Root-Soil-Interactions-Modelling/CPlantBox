#include <sstream>
#include <fstream>
#include <algorithm>

/**
 * Example for Shehan
 *
 * Scale elongation rate according to some 1d tabular data
 *
 * proof of concept, please adjust
 */

namespace CRootBox {

/*
 *  reads a csv into a vector^2
 */
std::vector<std::vector<double>> readCSV(std::string name, char delimeter, int ignore_rows, int ignore_cols)
        {
    std::ifstream file(name);
    std::vector<std::vector<double>> data(0);
    int rc = 0;
    std::string line = "";
    while (getline(file, line)) {
        rc++;
        if (rc>ignore_rows) {
            std::stringstream rss(line);
            std::vector<double> row(0);
            int cc = 0;
            std::string val;
            while (std::getline(rss, val, delimeter)) {
                cc++;
                if (cc>ignore_cols) {
                    std::replace(val.begin(), val.end(), ',', '.'); // replace all ',' to '.'
                    std::stringstream convert(val);
                    double d;
                    convert >> d;
                    row.push_back(d);
                    // std::cout << d << ", ";
                }
            }
            // std::cout << "\n";
            data.push_back(row);
        }
    }
    return data;
}

/**
 * The elongation rate is calculated in dependence of water content and temperature,
 * The empirical relationship is located in the private function f.
 *
 */
class ScaleElongation :public SoilLookUp
{

public:

    ScaleElongation(Grid1D* wc, Grid1D* temp) {
        water_content = wc;
        temperature = temp;
    }

    virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const { ///< Returns a scalar property of the soil, 1. per default
        double wc = water_content->getValue(pos,root);
        double temp = temperature->getValue(pos,root);
        double y = f(wc, temp);
        // std::cout << "Root type : " << root->param.type << " pos " << pos.z << ", wc " << wc << ", temp " << temp << ", scaling " << y << "\n";
        return y;
    }

    virtual std::string toString() const { return "ScaleElongation based on water content and temperature"; } ///< Quick info about the object for debugging

    Grid1D* water_content;
    Grid1D* temperature;

private:

    // this function describes the dependency the elongation scale to
    double f(double wc, double temp) const {
        if (wc<tp) {
            return 0;
        }
        if (wc>ts) {
            return 1;
        }

        double Sp;
        double impt;
        double P;

        P = ((Pb - Pbmin) / (Pbmax - Pbmin));
        Sp = ((wc - tp) / (ts - tp));
        double PR = exp(1.5 + (2 * P) - (4 * Sp)); // scaled temperature
        impt = (1 - (PR / PR1) / 100);
        //std::cout << impt << " \t";
        return impt;
    }

    double Pb = 1.55;// bulk density
    double Pbmin = 1.50;// minimum soil bulk density
    double Pbmax = 1.60;//maximum soil bulk density
    double tp = 0.10; //permanent wilting point
    double ts = 0.45;//saturation
    double PR1 = 0.223; //penetrometer resistant at saturation(Mpa)


};

/**
 * Example how to implement varying elongation rates
 */
void shehan_WC()
{
    using namespace std;

    RootSystem rootsystem;

    string name = "wheat";

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

    // temperature data
    std::cout << "reading temperature data ...";
    auto field_temp = readCSV("/home/daniel/workspace/CRootBox/examples/TEMP.csv", ';', 2, 1);
    std::vector<double> temp = field_temp.at(0);

    // water content data
    std::cout << "reading water content data ...";
    auto field_wc = readCSV("/home/daniel/workspace/CRootBox/examples/WC.csv", ';', 1, 1);
    std::vector<double> wc = field_wc.at(0);

    // create scale elongation function
    Grid1D water_content = Grid1D(n, z_, wc);
    Grid1D temperature = Grid1D(n, z_, temp );
    ScaleElongation se = ScaleElongation(&water_content, &temperature);
    for (int i = 1; i < 7; i++) {  // "manually" set the scale elongation function
        rootsystem.getRootRandomParameter(i)->f_se = &se;
    }

    /**
     * Initialize
     */
    rootsystem.initialize(); //it is important to call initialize() after setGeometry()

    /*
     * Simulate
     */
    double simtime = 8 * 30; // days
    double dt = 0.5 * 1./24.;  // dt is half an hour
    size_t N = round(simtime/dt);

    assert(N<=field_temp.size()); // check if enough data are available

    for (size_t i = 0; i < N; i++) {

        rootsystem.simulate(dt);

        // update field data:
        temperature.data =  field_temp.at(i);
        water_content.data = field_wc.at(i);

        //        auto rl = rootsystem.getScalar(RootSystem::st_length);
        //        double tl = std::accumulate(rl.begin(), rl.end(), 0);
        //        cout << "Time step i " << tl << " cm total length \n";

    }

    /**
     * Export results (as vtp)
     */
    rootsystem.write(name+".vtp");

    double tl = rootsystem.getSummed("length");
    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes, " << tl << " cm total length \n";
}

} // end namespace CRootBox
