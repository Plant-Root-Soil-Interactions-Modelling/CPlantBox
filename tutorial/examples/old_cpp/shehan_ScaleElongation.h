#include <sstream>
#include <fstream>

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
                    std::stringstream convert(val);
                    double d;
                    convert >> d;
                    row.push_back(d);
                    //std::cout << d << "\n";
                }
            }
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

    ScaleElongation(Grid1D* wc, Grid1D* temp)
{
        water_content = wc;
        temperature = temp;
}

    virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const ///< Returns a scalar property of the soil, 1. per default
    {
        // std::cout << "type : " << root->param.type << "\n";
        double wc = water_content->getValue(pos,root);
        double temp = temperature->getValue(pos,root);
        return f(wc,temp);
    }

    virtual std::string toString() const { return "ScaleElongation based on water content and temperature"; } ///< Quick info about the object for debugging

    Grid1D* water_content;
    Grid1D* temperature;

private:

    // this function describes the dependency the elongation scale to
    double f(double wc, double temp) const
    {
        if (temp<minT) {
            return 0;
        }
        if (temp>maxT) {
            return 0;
        }
        double sigma;
        double impt;
        if (optT < 0.5*(minT+maxT)) {
            sigma = log(0.5) /log((optT-minT)/(maxT-minT));
            double t = (temp-minT)/(maxT-minT); // scaled temperature
            impt = pow(sin(M_PI*t),sigma);
        } else {
            sigma = log(0.5) /log((optT-maxT)/(minT-maxT));
            double t = (temp-maxT)/(minT-maxT); // scaled temperature
            impt = pow(sin(M_PI*t),sigma);
        }
        // std::cout << impt << " \n";
        return impt;
    }

    double minT = 0;
    double maxT = 30;
    double optT = 20;

};

/**
 * Example how to implement varying elongation rates
 */
void shehan_ScaleElongation()
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

    // temperature data
    std::cout << "reading temperature data ...";
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
    // "manually" set the scale elongation function
    for (int i=1; i<7; i++) {
        rootsystem.getRootRandomParameter(i)->se = &se;
    }

    /**
     * Initialize
     */
    rootsystem.initialize(); //it is important to call initialize() after setGeometry()

    /*
     * Simulate
     */
    double simtime = 7*30; // days
    double dt = 0.5 * 1./24.;  // dt is half an hour
    size_t N = round(simtime/dt);

    assert(N<=field_temp.size()); // check if enough data are available

    for (size_t i=0; i<N; i++) {

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

    auto rl = rootsystem.getParameter("length");
    double tl = std::accumulate(rl.begin(), rl.end(), 0);
    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes, " << tl << " cm total length \n";
}

} // end namespace CRootBox
