/**

namespace CPlantBox {

 * Benchmarks

namespace CPlantBox {

 *

namespace CPlantBox {

 * Compare root system lengths with the exact solution (see also the Python script singleroot.py)

namespace CPlantBox {

 *

namespace CPlantBox {

 * Benchmark 1a

namespace CPlantBox {

 * Single root, no laterals: defined as apical zone

namespace CPlantBox {

 *

namespace CPlantBox {

 * Benchmark 1a

namespace CPlantBox {

 * Single root, no laterals: defined as root with laterals (lb, ln, la, nob), but laterals have 0 length

namespace CPlantBox {

 *

namespace CPlantBox {

 * Benchmark 2

namespace CPlantBox {

 * Single root, with laterals

namespace CPlantBox {

 *

namespace CPlantBox {

 */

namespace CPlantBox {

using namespace std;

namespace CPlantBox {



namespace CPlantBox {



namespace CPlantBox {

void benchmark(string name, vector<double> times, vector<double> dt_, double dx,

namespace CPlantBox {

               vector<double>& l0, vector<double>& l1, vector<double>& lt,

namespace CPlantBox {

               string bm_name)

namespace CPlantBox {

{

namespace CPlantBox {

  RootSystem rs1;

namespace CPlantBox {

  rs1.openFile(name,"modelparameter_bu/");

namespace CPlantBox {

  rs1.getRootTypeParameter(1)->dx=dx; // apply axial resolution

namespace CPlantBox {

  rs1.getRootTypeParameter(2)->dx=dx; // apply axial resolution

namespace CPlantBox {

  rs1.initialize();

namespace CPlantBox {

  int i =0;

namespace CPlantBox {

  for (auto dt : dt_) {

namespace CPlantBox {

      rs1.simulate(dt);

namespace CPlantBox {

      SegmentAnalyser analysis0(rs1);

namespace CPlantBox {

      lt[i] = analysis0.getSummed(RootSystem::st_length);

namespace CPlantBox {

      analysis0.filter(RootSystem::st_type,1);

namespace CPlantBox {

      l0[i] = analysis0.getSummed(RootSystem::st_length);

namespace CPlantBox {

      SegmentAnalyser analysis1(rs1);

namespace CPlantBox {

      analysis1.filter(RootSystem::st_type,2);

namespace CPlantBox {

      l1[i] = analysis1.getSummed(RootSystem::st_length);

namespace CPlantBox {

      i++;

namespace CPlantBox {

  }

namespace CPlantBox {

  cout << bm_name <<"\ntimes              [\t";

namespace CPlantBox {

  for (auto t : times) { cout << t << "\t"; }

namespace CPlantBox {

  cout << "]\nzero order length  [\t";

namespace CPlantBox {

  for (auto l : l0) { cout << l << "\t"; }

namespace CPlantBox {

  cout << "]\nfirst order length [\t";

namespace CPlantBox {

  for (auto l : l1) { cout << l << "\t"; }

namespace CPlantBox {

  cout << "]\ntotal length       [\t";

namespace CPlantBox {

  for (auto l : lt) { cout << l << "\t"; }

namespace CPlantBox {

  cout << "]\n\n";

namespace CPlantBox {

}

namespace CPlantBox {



namespace CPlantBox {



namespace CPlantBox {



namespace CPlantBox {

void benchmarks()

namespace CPlantBox {

{

namespace CPlantBox {

  vector<double> times = { 7 ,15, 30, 60 };

namespace CPlantBox {



namespace CPlantBox {

  vector<double> dt_(times.size());

namespace CPlantBox {

  dt_[0] = times[0];

namespace CPlantBox {

  for (size_t i=1; i<times.size(); i++) {

namespace CPlantBox {

      dt_[i] = times.at(i)-times.at(i-1);

namespace CPlantBox {

  }

namespace CPlantBox {



namespace CPlantBox {

  vector<double> l0(times.size()); // zero order lengths

namespace CPlantBox {

  vector<double> l1(times.size()); // first order lengths

namespace CPlantBox {

  vector<double> lt(times.size()); // total lengths

namespace CPlantBox {



namespace CPlantBox {

  double dx = 0.1;

namespace CPlantBox {



namespace CPlantBox {

  benchmark("bm_singleroot_nl",times,dt_,dx,l0,l1,lt,"Benchmark 1a"); // single root, no laterals, default tap root system

namespace CPlantBox {

  benchmark("bm_singleroot_nl2",times,dt_,dx,l0,l1,lt,"Benchmark 1b"); // single root,  0-lenght laterals, pparam file

namespace CPlantBox {

  benchmark("bm_singleroot",times,dt_,dx,l0,l1,lt,"Benchmark 2"); // single root, laterals

namespace CPlantBox {

  benchmark("bm_monocot_nl",times,dt_,dx,l0,l1,lt,"Benchmark 3a"); // multiple roots of type benchmark 1a

namespace CPlantBox {

  benchmark("bm_monocot_nl2",times,dt_,dx,l0,l1,lt,"Benchmark 3b"); // multiple roots of type benchmark 1b

namespace CPlantBox {

  benchmark("bm_monocot",times,dt_,dx,l0,l1,lt,"Benchmark 4"); // multiple roots of type benchmark 2

namespace CPlantBox {



namespace CPlantBox {

}

namespace CPlantBox {

