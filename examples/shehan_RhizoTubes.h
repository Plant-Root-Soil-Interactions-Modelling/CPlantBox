/**

namespace CPlantBox {

 * Shehans experimental set-up

namespace CPlantBox {

 *

namespace CPlantBox {

 * (B) Virtual analysis using rhizo tubes

namespace CPlantBox {

 *

namespace CPlantBox {

 * simulates the growing root systems (6*37) in the field

namespace CPlantBox {

 * uses trenches to analyse the number of roots in a coarse grid

namespace CPlantBox {

 * *

namespace CPlantBox {

 */

namespace CPlantBox {

using namespace std;

namespace CPlantBox {



namespace CPlantBox {

/*

namespace CPlantBox {

 * Creates the geometry of the rhizo tubes

namespace CPlantBox {

 *

namespace CPlantBox {

 * @param r     radius of the rhizo tubes

namespace CPlantBox {

 * \return union of the rhizo tubes

namespace CPlantBox {

 */

namespace CPlantBox {

SignedDistanceFunction* fieldRhizoTubes(double r)

namespace CPlantBox {

{

namespace CPlantBox {

  double l = 90; // cm tube length

namespace CPlantBox {

  SDF_PlantContainer* rhizotube = new SDF_PlantContainer(r,r,l,false);

namespace CPlantBox {

  SDF_RotateTranslate* rhizoX = new SDF_RotateTranslate(rhizotube, l, SDF_RotateTranslate::yaxis, Vector3d(l,0.,0.));

namespace CPlantBox {

  vector<SignedDistanceFunction*> rhizotubes_;

namespace CPlantBox {

  const int tubeN=6;

namespace CPlantBox {

  double y_[tubeN] = { 35, 45, 55, 65, 75, 85 };

namespace CPlantBox {

  double z_[tubeN] = { -10, -20, -40, -60, -80, -120 };

namespace CPlantBox {

  for (int i=0; i<tubeN; i++) {

namespace CPlantBox {

      rhizotubes_.push_back(new SDF_RotateTranslate(rhizoX, 0, SDF_RotateTranslate::xaxis, Vector3d(0,y_[i],z_[i])));

namespace CPlantBox {

  }

namespace CPlantBox {

  SDF_Union* rhizotubes = new SDF_Union(rhizotubes_);

namespace CPlantBox {

  return rhizotubes;

namespace CPlantBox {

}

namespace CPlantBox {



namespace CPlantBox {

/*

namespace CPlantBox {

 * The positions in the rhizotubes were fotos are taken

namespace CPlantBox {

 *

namespace CPlantBox {

 * \return Positions of the fotos

namespace CPlantBox {

 */

namespace CPlantBox {

vector<Vector3d> fotoPos() {

namespace CPlantBox {

  double l = 90; // cm

namespace CPlantBox {

  vector<Vector3d> fP;

namespace CPlantBox {

  const int tubeN=6;

namespace CPlantBox {

  double y_[tubeN] = { 35, 45, 55, 65, 75, 85 };

namespace CPlantBox {

  double z_[tubeN]= { -10, -20, -40, -60, -80, -120 };

namespace CPlantBox {

  int nf = 7; // number of foto positions per tube

namespace CPlantBox {

  double dist = l/double(nf+2);  // = 10 cm

namespace CPlantBox {

  for (int i=0; i<tubeN; i++) {

namespace CPlantBox {

      for (int j=0; j<nf; j++) {

namespace CPlantBox {

          fP.push_back(Vector3d((j+1)*dist,y_[i],z_[i]));

namespace CPlantBox {

      }

namespace CPlantBox {

  }

namespace CPlantBox {

  return fP;

namespace CPlantBox {

}

namespace CPlantBox {



namespace CPlantBox {

/**

namespace CPlantBox {

 * Rhizo tube analysis in the virtual field experiment

namespace CPlantBox {

 */

namespace CPlantBox {

void shehan_RhizoTubes(string name = "wheat", bool exportVTP = false)

namespace CPlantBox {

{

namespace CPlantBox {

  /*

namespace CPlantBox {

   * Initialize

namespace CPlantBox {

   */

namespace CPlantBox {

  double r = 3.2;

namespace CPlantBox {

  SignedDistanceFunction* geometry = new SDF_Complement(fieldRhizoTubes(r));

namespace CPlantBox {

  auto allRS = initializeRootSystems(name,geometry);

namespace CPlantBox {



namespace CPlantBox {

  /*

namespace CPlantBox {

   * Simulate

namespace CPlantBox {

   */

namespace CPlantBox {

  vector<double> times = {0, 30, 60, 90, 120, 150, 180, 210, 240};

namespace CPlantBox {

  //vector<double> times = {0, 7, 14, 21, 30, 37, 44, 51, 60, 67, 74, 81, 90, 97, 104, 111, 120, 127, 134, 141, 150, 157, 164, 171, 180, 187, 194, 201, 210,

namespace CPlantBox {

  //		217, 223, 231, 240}; // test...

namespace CPlantBox {

  simulateRS(times, allRS);

namespace CPlantBox {



namespace CPlantBox {

  /*

namespace CPlantBox {

   * Analysis: crop to the segments visible from the rhizotubes

namespace CPlantBox {

   */

namespace CPlantBox {

  double dx = 0.3; //cm

namespace CPlantBox {

  SignedDistanceFunction* visible = fieldRhizoTubes(3.2+dx);

namespace CPlantBox {



namespace CPlantBox {

  vector<SegmentAnalyser> croppedTubes;

namespace CPlantBox {

  for (size_t i=1; i<times.size(); i++) {

namespace CPlantBox {

      std::cout << "\nANALYSE TIME " << times.at(i) <<"\n\n";

namespace CPlantBox {

      SegmentAnalyser analyser = getResult(allRS,times.at(i));

namespace CPlantBox {

      analyser.crop(visible); // throw segments away

namespace CPlantBox {

      analyser.pack(); // throw nodes away

namespace CPlantBox {

      if (exportVTP) {

namespace CPlantBox {

          string vtpname = name + "_tube_cropped"+ std::to_string(i)+".vtp"; // export cropped segments for visualization

namespace CPlantBox {

          analyser.write(vtpname);

namespace CPlantBox {

      }

namespace CPlantBox {

      croppedTubes.push_back(analyser); // copies into the vector

namespace CPlantBox {

  }

namespace CPlantBox {



namespace CPlantBox {

  /*

namespace CPlantBox {

   * Foto

namespace CPlantBox {

   */

namespace CPlantBox {

  vector<vector<double>> finalmatrixC(times.size()-1);

namespace CPlantBox {

  vector<vector<double>> finalmatrixL(times.size()-1);

namespace CPlantBox {



namespace CPlantBox {

  vector<Vector3d> fpos = fotoPos();

namespace CPlantBox {



namespace CPlantBox {

  vector<SegmentAnalyser> fotos; // not real fotos, but slices of the croopedTube

namespace CPlantBox {



namespace CPlantBox {

  double xx = 1.4; // cm slice thickness

namespace CPlantBox {

  double a1 = 8./180.*M_PI;

namespace CPlantBox {

  double a2 = 22./180.*M_PI;

namespace CPlantBox {

  double r_ = r+dx; // 3.5

namespace CPlantBox {

  double zz = sin(a1)*r_ + sin(a2)*r_;

namespace CPlantBox {

  cout << " foto box " << xx << ", 100, "<<zz<< "\n";

namespace CPlantBox {

  SDF_PlantBox fotoBox(xx,100.,zz);

namespace CPlantBox {



namespace CPlantBox {

  vector<SignedDistanceFunction*> allcameras;

namespace CPlantBox {



namespace CPlantBox {

  for (size_t t=0; t<times.size()-1; t++) {

namespace CPlantBox {



namespace CPlantBox {

      int c = 0;

namespace CPlantBox {

      cout << "Time: " << times.at(t) << "\n";

namespace CPlantBox {

      SegmentAnalyser empty = SegmentAnalyser();

namespace CPlantBox {

      fotos.push_back(empty); // fill in this time step

namespace CPlantBox {

      SegmentAnalyser f = croppedTubes.at(t);

namespace CPlantBox {



namespace CPlantBox {

      for (auto& p :fpos) {

namespace CPlantBox {



namespace CPlantBox {

          int j = floor(double(c)/7.);

namespace CPlantBox {



namespace CPlantBox {

          if (int(finalmatrixC[t].size())<=j) {

namespace CPlantBox {

              finalmatrixC[t].push_back(0.);

namespace CPlantBox {

              finalmatrixL[t].push_back(0.);

namespace CPlantBox {

          }

namespace CPlantBox {



namespace CPlantBox {

          // foto 1 + 2

namespace CPlantBox {

          SDF_RotateTranslate* foto = new SDF_RotateTranslate(&fotoBox,

namespace CPlantBox {

                                                              Vector3d(p.x,p.y,p.z+r_*sin(a2)));

namespace CPlantBox {

          SegmentAnalyser foto12 = f;

namespace CPlantBox {

          foto12.crop(foto);

namespace CPlantBox {

          foto12.pack();

namespace CPlantBox {



namespace CPlantBox {

          finalmatrixC[t][j] += double(foto12.getNumberOfRoots())/7.;

namespace CPlantBox {

          // number of segments is not feasible, since a root might grow around the tube

namespace CPlantBox {

          finalmatrixL[t][j] += double(foto12.getSummed(RootSystem::st_length))/7.;

namespace CPlantBox {

//          finalmatrixC[t][j] += 1./7.; // for debugging

namespace CPlantBox {

//          finalmatrixL[t][j] += 1./7.;

namespace CPlantBox {



namespace CPlantBox {

          std::cout << "#" << c << " foto at location: "<< p.toString() << " tube " << j  << "\n";

namespace CPlantBox {



namespace CPlantBox {

          allcameras.push_back(foto); // save geometry for later vizualisation

namespace CPlantBox {

          fotos.back().addSegments(foto12); // save segments for vizualisation

namespace CPlantBox {

          c++;

namespace CPlantBox {

      }

namespace CPlantBox {

  }

namespace CPlantBox {



namespace CPlantBox {

  /*

namespace CPlantBox {

   * export fotos

namespace CPlantBox {

   */

namespace CPlantBox {

  int c=0;

namespace CPlantBox {

  for (auto& ft : fotos) {

namespace CPlantBox {

    string vtpname = name +  "_tube_fotos"+std::to_string(++c) +".vtp";

namespace CPlantBox {

    ft.write(vtpname);

namespace CPlantBox {

  }

namespace CPlantBox {



namespace CPlantBox {

  int n = finalmatrixC[0].size(); // 32

namespace CPlantBox {

  int m = times.size()-1; // 8

namespace CPlantBox {

  string matname = name+"_tube_matrix1.txt";

namespace CPlantBox {

  std::ofstream fos;

namespace CPlantBox {

  fos.open(matname.c_str());

namespace CPlantBox {

  for (int j=0; j<n; j++) {

namespace CPlantBox {

      for (int i=0; i<m; i++) {

namespace CPlantBox {

          std::cout << std::fixed << std::setprecision(4)<< finalmatrixC[i][j] << "\t";

namespace CPlantBox {

          fos << std::fixed << std::setprecision(4)<< finalmatrixC[i][j] << "\t";

namespace CPlantBox {

      }

namespace CPlantBox {

      std::cout << "\n";

namespace CPlantBox {

      fos << "\n";

namespace CPlantBox {

  }

namespace CPlantBox {

  fos.close();

namespace CPlantBox {



namespace CPlantBox {

  n = finalmatrixL[0].size(); // 32

namespace CPlantBox {

  m = times.size()-1; // 8

namespace CPlantBox {

  matname = name+"_tube_matrix2.txt";

namespace CPlantBox {

  fos.open(matname.c_str());

namespace CPlantBox {

  for (int j=0; j<n; j++) {

namespace CPlantBox {

      for (int i=0; i<m; i++) {

namespace CPlantBox {

          std::cout << std::fixed << std::setprecision(4)<< finalmatrixL[i][j] << "\t";

namespace CPlantBox {

          fos << std::fixed << std::setprecision(4)<< finalmatrixL[i][j] << "\t";

namespace CPlantBox {

      }

namespace CPlantBox {

      std::cout << "\n";

namespace CPlantBox {

      fos << "\n";

namespace CPlantBox {

  }

namespace CPlantBox {

  fos.close();

namespace CPlantBox {



namespace CPlantBox {

  /**

namespace CPlantBox {

   * Export rhizotube geometry

namespace CPlantBox {

   */

namespace CPlantBox {

  if (exportVTP) {

namespace CPlantBox {

      string gname = name + "_tube.py";

namespace CPlantBox {

      allRS[0]->write(gname);

namespace CPlantBox {

      string gname2 = name + "_tubeCamera.py";

namespace CPlantBox {

      SDF_Union cam(allcameras);

namespace CPlantBox {

      allRS[0]->setGeometry(&cam);

namespace CPlantBox {

      allRS[0]->write(gname2);

namespace CPlantBox {

  }

namespace CPlantBox {

}

namespace CPlantBox {

