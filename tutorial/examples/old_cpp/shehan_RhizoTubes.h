/**
 * Shehans experimental set-up
 *
 * (B) Virtual analysis using rhizo tubes
 *
 * simulates the growing root systems (6*37) in the field
 * uses trenches to analyse the number of roots in a coarse grid
 * *
 */
namespace CRootBox {

/*
 * Creates the geometry of the rhizo tubes
 *
 * @param r     radius of the rhizo tubes
 * \return union of the rhizo tubes
 */
SignedDistanceFunction* fieldRhizoTubes(double r)
{
  using namespace std;
  double l = 90; // cm tube length
  SDF_PlantContainer* rhizotube = new SDF_PlantContainer(r,r,l,false);
  SDF_RotateTranslate* rhizoX = new SDF_RotateTranslate(rhizotube, l, SDF_RotateTranslate::yaxis, Vector3d(l,0.,0.));
  vector<SignedDistanceFunction*> rhizotubes_;
  const int tubeN = 6;
  double y_[tubeN] = { 35, 45, 55, 65, 75, 85 };
  double z_[tubeN] = { -10, -20, -40, -60, -80, -120 };
  for (int i=0; i<tubeN; i++) {
      rhizotubes_.push_back(new SDF_RotateTranslate(rhizoX, 0, SDF_RotateTranslate::xaxis, Vector3d(0,y_[i],z_[i])));
  }
  SDF_Union* rhizotubes = new SDF_Union(rhizotubes_);
  return rhizotubes;
}

/*
 * The positions in the rhizotubes were fotos are taken
 *
 * \return Positions of the fotos
 */
std::vector<Vector3d> fotoPos() {
  using namespace std;
  double l = 90; // cm
  vector<Vector3d> fP;
  const int tubeN=6;
  double y_[tubeN] = { 35, 45, 55, 65, 75, 85 };
  double z_[tubeN]= { -10, -20, -40, -60, -80, -120 };
  int nf = 7; // number of foto positions per tube
  double dist = l/double(nf+2);  // = 10 cm
  for (int i=0; i<tubeN; i++) {
      for (int j=0; j<nf; j++) {
          fP.push_back(Vector3d((j+1)*dist,y_[i],z_[i]));
      }
  }
  return fP;
}

/**
 * Rhizo tube analysis in the virtual field experiment
 */
void shehan_RhizoTubes(const std::string& name = "wheat", bool exportVTP = false)
{
  using namespace std;
  /*
   * Initialize
   */
  double r = 3.2;
  SignedDistanceFunction* geometry = new SDF_Complement(fieldRhizoTubes(r));
  auto allRS = initializeRootSystems(name,geometry);

  /*
   * Simulate
   */
  vector<double> times = {0, 30, 60, 90, 120, 150, 180, 210, 240};
  //vector<double> times = {0, 7, 14, 21, 30, 37, 44, 51, 60, 67, 74, 81, 90, 97, 104, 111, 120, 127, 134, 141, 150, 157, 164, 171, 180, 187, 194, 201, 210,
  //		217, 223, 231, 240}; // test...
  simulateRS(times, allRS);

  /*
   * Analysis: crop to the segments visible from the rhizotubes
   */
  double dx = 0.3; //cm
  SignedDistanceFunction* visible = fieldRhizoTubes(3.2+dx);

  vector<SegmentAnalyser> croppedTubes;
  for (size_t i=1; i<times.size(); i++) {
      std::cout << "\nANALYSE TIME " << times.at(i) <<"\n\n";
      SegmentAnalyser analyser = getResult(allRS,times.at(i));
      analyser.crop(visible); // throw segments away
      analyser.pack(); // throw nodes away
      if (exportVTP) {
          string vtpname = name + "_tube_cropped"+ std::to_string(i)+".vtp"; // export cropped segments for visualization
          analyser.write(vtpname);
      }
      croppedTubes.push_back(analyser); // copies into the vector
  }

  /*
   * Foto
   */
  vector<vector<double>> finalmatrixC(times.size()-1);
  vector<vector<double>> finalmatrixL(times.size()-1);

  vector<Vector3d> fpos = fotoPos();

  vector<SegmentAnalyser> fotos; // not real fotos, but slices of the croopedTube

  double xx = 1.4; // cm slice thickness
  double a1 = 8./180.*M_PI;
  double a2 = 22./180.*M_PI;
  double r_ = r+dx; // 3.5
  double zz = sin(a1)*r_ + sin(a2)*r_;
  cout << " foto box " << xx << ", 100, "<<zz<< "\n";
  SDF_PlantBox fotoBox(xx,100.,zz);

  vector<SignedDistanceFunction*> allcameras;

  for (size_t t=0; t<times.size()-1; t++) {

      int c = 0;
      cout << "Time: " << times.at(t) << "\n";
      SegmentAnalyser empty = SegmentAnalyser();
      fotos.push_back(empty); // fill in this time step
      SegmentAnalyser f = croppedTubes.at(t);

      for (auto& p :fpos) {

          int j = floor(double(c)/7.);

          if (int(finalmatrixC[t].size())<=j) {
              finalmatrixC[t].push_back(0.);
              finalmatrixL[t].push_back(0.);
          }

          // foto 1 + 2
          SDF_RotateTranslate* foto = new SDF_RotateTranslate(&fotoBox,
                                                              Vector3d(p.x,p.y,p.z+r_*sin(a2)));
          SegmentAnalyser foto12 = f;
          foto12.crop(foto);
          foto12.pack();

          finalmatrixC[t][j] += double(foto12.getNumberOfOrgans())/7.;
          // number of segments is not feasible, since a root might grow around the tube
          finalmatrixL[t][j] += double(foto12.getSummed("length"))/7.;
//          finalmatrixC[t][j] += 1./7.; // for debugging
//          finalmatrixL[t][j] += 1./7.;

          std::cout << "#" << c << " foto at location: "<< p.toString() << " tube " << j  << "\n";

          allcameras.push_back(foto); // save geometry for later vizualisation
          fotos.back().addSegments(foto12); // save segments for vizualisation
          c++;
      }
  }

  /*
   * export fotos
   */
  int c=0;
  for (auto& ft : fotos) {
    string vtpname = name +  "_tube_fotos"+std::to_string(++c) +".vtp";
    ft.write(vtpname);
  }

  int n = finalmatrixC[0].size(); // 32
  int m = times.size()-1; // 8
  string matname = name+"_tube_matrix1.txt";
  std::ofstream fos;
  fos.open(matname.c_str());
  for (int j=0; j<n; j++) {
      for (int i=0; i<m; i++) {
          std::cout << std::fixed << std::setprecision(4)<< finalmatrixC[i][j] << "\t";
          fos << std::fixed << std::setprecision(4)<< finalmatrixC[i][j] << "\t";
      }
      std::cout << "\n";
      fos << "\n";
  }
  fos.close();

  n = finalmatrixL[0].size(); // 32
  m = times.size()-1; // 8
  matname = name+"_tube_matrix2.txt";
  fos.open(matname.c_str());
  for (int j=0; j<n; j++) {
      for (int i=0; i<m; i++) {
          std::cout << std::fixed << std::setprecision(4)<< finalmatrixL[i][j] << "\t";
          fos << std::fixed << std::setprecision(4)<< finalmatrixL[i][j] << "\t";
      }
      std::cout << "\n";
      fos << "\n";
  }
  fos.close();

  /**
   * Export rhizotube geometry
   */
  if (exportVTP) {
      string gname = name + "_tube.py";
      allRS[0]->write(gname);
      string gname2 = name + "_tubeCamera.py";
      SDF_Union cam(allcameras);
      allRS[0]->setGeometry(&cam);
      allRS[0]->write(gname2);
  }
}

} // end namespace CRootBox
