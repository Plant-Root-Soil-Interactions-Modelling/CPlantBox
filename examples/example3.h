/**
 * Example 3
 *
 * Creates N*N root systems at different locations, no confining geometry
 *
 */
using namespace std;

void example3()
{
  auto gen = default_random_engine(chrono::system_clock::now().time_since_epoch().count());
  auto UD = uniform_real_distribution<double>(0,1); // random stuff

  string name = "Zea_mays_5_Leitner_2014";

  /*
   * Creates N*N root systems
   */
  int N=3;
  double dist = 40; // [cm]
  vector<Plant*> allRS;
  for (int i=0; i<N; i++) {
      for (int j=0;j<N; j++) {
          Plant* rs = new Plant();
          allRS.push_back(rs);
          rs->openFile(name);
          rs->getRootSystemParameter()->seedPos = Vector3d(dist*i,dist*j,-3); // set position of seed [cm]
      }
  }

  /*
   * Simulate
   */
  double simtime = 120; // days
  for (auto rs : allRS) { // simulate all
      double s = UD(gen);
      rs->setSeed(s); // randomly select a seed
      rs->initialize(4,5);
      rs->simulate(simtime);
      cout << "Finished with a total of " << rs->getNumberOfNodes()<< " nodes\n";
  }

  /*
   * Export results as single vtp files
   */
  int c=0;
  for (auto rs : allRS) {
      c++; // root system number
      string vtpname = "results/"+name+std::to_string(c)+".vtp";
      rs->write(vtpname);
  }

}
