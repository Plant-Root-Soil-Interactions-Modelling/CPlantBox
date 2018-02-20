/**
 * Benchmarks
 *
 * Compare root system lengths with the exact solution (see also the Python script singleroot.py)
 *
 * Benchmark 1a
 * Single root, no laterals: defined as apical zone
 *
 * Benchmark 1a
 * Single root, no laterals: defined as root with laterals (lb, ln, la, nob), but laterals have 0 length
 *
 * Benchmark 2
 * Single root, with laterals
 *
 */
using namespace std;


void benchmark(string name, vector<double> times, vector<double> dt_, double dx,
               vector<double>& l0, vector<double>& l1, vector<double>& lt,
               string bm_name)
{
  RootSystem rs1;
  rs1.openFile(name,"modelparameter_bu/");
  rs1.getRootTypeParameter(1)->dx=dx; // apply axial resolution
  rs1.getRootTypeParameter(2)->dx=dx; // apply axial resolution
  rs1.initialize();
  int i =0;
  for (auto dt : dt_) {
      rs1.simulate(dt);
      SegmentAnalyser analysis0(rs1);
      lt[i] = analysis0.getSummed(RootSystem::st_length);
      analysis0.filter(RootSystem::st_type,1);
      l0[i] = analysis0.getSummed(RootSystem::st_length);
      SegmentAnalyser analysis1(rs1);
      analysis1.filter(RootSystem::st_type,2);
      l1[i] = analysis1.getSummed(RootSystem::st_length);
      i++;
  }
  cout << bm_name <<"\ntimes              [\t";
  for (auto t : times) { cout << t << "\t"; }
  cout << "]\nzero order length  [\t";
  for (auto l : l0) { cout << l << "\t"; }
  cout << "]\nfirst order length [\t";
  for (auto l : l1) { cout << l << "\t"; }
  cout << "]\ntotal length       [\t";
  for (auto l : lt) { cout << l << "\t"; }
  cout << "]\n\n";
}



void benchmarks()
{
  vector<double> times = { 7 ,15, 30, 60 };

  vector<double> dt_(times.size());
  dt_[0] = times[0];
  for (size_t i=1; i<times.size(); i++) {
      dt_[i] = times.at(i)-times.at(i-1);
  }

  vector<double> l0(times.size()); // zero order lengths
  vector<double> l1(times.size()); // first order lengths
  vector<double> lt(times.size()); // total lengths

  double dx = 0.1;

  benchmark("bm_singleroot_nl",times,dt_,dx,l0,l1,lt,"Benchmark 1a"); // single root, no laterals, default tap root system
  benchmark("bm_singleroot_nl2",times,dt_,dx,l0,l1,lt,"Benchmark 1b"); // single root,  0-lenght laterals, pparam file
  benchmark("bm_singleroot",times,dt_,dx,l0,l1,lt,"Benchmark 2"); // single root, laterals
  benchmark("bm_monocot_nl",times,dt_,dx,l0,l1,lt,"Benchmark 3a"); // multiple roots of type benchmark 1a
  benchmark("bm_monocot_nl2",times,dt_,dx,l0,l1,lt,"Benchmark 3b"); // multiple roots of type benchmark 1b
  benchmark("bm_monocot",times,dt_,dx,l0,l1,lt,"Benchmark 4"); // multiple roots of type benchmark 2

}
