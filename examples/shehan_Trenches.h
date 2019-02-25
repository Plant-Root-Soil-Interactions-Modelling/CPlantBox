/**

namespace CPlantBox {

 * Shehans experimental set-up

namespace CPlantBox {

 *

namespace CPlantBox {

 * (C) Virtual analysis of trenches

namespace CPlantBox {

 *

namespace CPlantBox {

 * simulates the growing root systems (6*37) in the field

namespace CPlantBox {

 * uses trenches to analyse the number of roots in a coarse grid

namespace CPlantBox {

 *

namespace CPlantBox {

 */

namespace CPlantBox {

using namespace std;

namespace CPlantBox {



namespace CPlantBox {



namespace CPlantBox {

/**

namespace CPlantBox {

 * Creates the geometry of the trenches

namespace CPlantBox {

 *

namespace CPlantBox {

 * \return a vector of 9 infinitely large trenches

namespace CPlantBox {

 */

namespace CPlantBox {

vector<SDF_HalfPlane*> fieldTrenches() {

namespace CPlantBox {

	int N = 9; // number of trenches

namespace CPlantBox {

	vector<SDF_HalfPlane*> trenches;

namespace CPlantBox {

	double dist1 = 18; // cm horizontal distance between plants

namespace CPlantBox {

	double dist2 = 6; 	// cm vertical distance between the trenches

namespace CPlantBox {

	for (int i=0; i<N; i++) {

namespace CPlantBox {

		Vector3d o(1.5*dist1, 5.*dist2+i*dist2,  -160);

namespace CPlantBox {

		std::cout << o.toString() << "\n";

namespace CPlantBox {

		Vector3d p1(1.5*dist1, 5.*dist2+i*dist2,  0);

namespace CPlantBox {

		Vector3d p2(3.5*dist1, 5.*dist2+i*dist2,  -160);

namespace CPlantBox {

		trenches.push_back(new SDF_HalfPlane(o,p1,p2));

namespace CPlantBox {

	}

namespace CPlantBox {

	return trenches;

namespace CPlantBox {

}

namespace CPlantBox {



namespace CPlantBox {

/**

namespace CPlantBox {

 * Creates a field of multiply root system relevant for the trenches (4*21)

namespace CPlantBox {

 *

namespace CPlantBox {

 * @param name     filename of the parameter file

namespace CPlantBox {

 * @param geom     confining geometry

namespace CPlantBox {

 */

namespace CPlantBox {

vector<RootSystem*> initializeRootSystems2(string name, SignedDistanceFunction* geom = new SignedDistanceFunction())

namespace CPlantBox {

						{

namespace CPlantBox {

	// auto gen = mt19937(chrono::system_clock::now().time_since_epoch().count());

namespace CPlantBox {

	// auto UD = uniform_real_distribution<double>(0,1); // random stuff, does it work now?

namespace CPlantBox {

	int M=6;

namespace CPlantBox {

	int N=37;

namespace CPlantBox {

	double dist1 = 18; // [cm] (M-1)*18 = 90

namespace CPlantBox {

	double dist2 = 3; // [cm] (N-1)*3 = 108

namespace CPlantBox {

	vector<RootSystem*> allRS;

namespace CPlantBox {

	for (int i=1; i<(M-1); i++) {

namespace CPlantBox {

		for (int j=8;j<(N-8); j++) {

namespace CPlantBox {

			RootSystem* rs = new RootSystem();

namespace CPlantBox {

			// double s = UD(gen);

namespace CPlantBox {

			// rs->setSeed(s); // randomly select a seed

namespace CPlantBox {

			allRS.push_back(rs);

namespace CPlantBox {

			rs->openFile(name);

namespace CPlantBox {

			rs->setGeometry(geom);

namespace CPlantBox {

			rs->getRootTypeParameter(4)->theta = 80./180.*M_PI; // fix insertion angle of the basal roots

namespace CPlantBox {

			auto pos = Vector3d(dist1*i,dist2*j,-3);

namespace CPlantBox {

			std::cout << pos.toString() << "\n";

namespace CPlantBox {

			rs->getRootSystemParameter()->seedPos = pos; // set position of seed [cm]

namespace CPlantBox {

			rs->initialize();

namespace CPlantBox {

		}

namespace CPlantBox {

	}

namespace CPlantBox {

	return allRS;

namespace CPlantBox {

}

namespace CPlantBox {



namespace CPlantBox {



namespace CPlantBox {



namespace CPlantBox {

/**

namespace CPlantBox {

 * Trench analysis in the virtual field experiment

namespace CPlantBox {

 */

namespace CPlantBox {

void shehan_Trenches(string name = "wheat", bool exportVTP = false)

namespace CPlantBox {

{

namespace CPlantBox {



namespace CPlantBox {

	/*

namespace CPlantBox {

	 * Initialize

namespace CPlantBox {

	 */

namespace CPlantBox {

	auto allRS = initializeRootSystems2(name); // no geometry, trenches are destructive and created after roots are grown

namespace CPlantBox {



namespace CPlantBox {

	/*

namespace CPlantBox {

	 * Simulate

namespace CPlantBox {

	 */

namespace CPlantBox {

	vector<double> times = {0, 30}; //, 60, 90, 120, 150, 180, 210, 240};

namespace CPlantBox {

	simulateRS(times, allRS);

namespace CPlantBox {



namespace CPlantBox {

	/*

namespace CPlantBox {

	 * Analysis: cut the roots along the trenches

namespace CPlantBox {

	 */

namespace CPlantBox {

	double dist1 = 18; // cm

namespace CPlantBox {

	int m = 4;

namespace CPlantBox {

	int n = 32;

namespace CPlantBox {

	SDF_PlantBox box_(36.,1.e9,160.);

namespace CPlantBox {

	//SDF_RotateTranslate box(&box_, Vector3d(18.+ 27., 0., 0.) );

namespace CPlantBox {

	vector<SDF_HalfPlane*> trenches = fieldTrenches();

namespace CPlantBox {

	vector<vector<vector<double>>> finalMT;

namespace CPlantBox {



namespace CPlantBox {

	for (size_t t=0; t<times.size()-1; t++) {

namespace CPlantBox {



namespace CPlantBox {

		std::cout<<"create empty matrix";

namespace CPlantBox {

		// create empty n*m matrix

namespace CPlantBox {

		vector<vector<double>> finalmatrix(n);

namespace CPlantBox {

		for (int i=0; i<n; i++) {

namespace CPlantBox {

			for (int j=0; j<m; j++) {

namespace CPlantBox {

				finalmatrix[i].push_back(0.);

namespace CPlantBox {

			}

namespace CPlantBox {

		}

namespace CPlantBox {



namespace CPlantBox {

		// loop over trenches

namespace CPlantBox {

		int ii=0;

namespace CPlantBox {

		for (auto& tr : trenches) {

namespace CPlantBox {

		//auto& tr = trenches[0];

namespace CPlantBox {

			std::string s;

namespace CPlantBox {

			std::ostringstream os;

namespace CPlantBox {

			os << "trench_" << ++ii;

namespace CPlantBox {

			std::cout << os.str() <<"\n";

namespace CPlantBox {



namespace CPlantBox {

			// cut along the trench

namespace CPlantBox {

			std::cout << "cut \n";

namespace CPlantBox {

			SegmentAnalyser analyser = getResult(allRS,times.at(t+1));

namespace CPlantBox {

			SegmentAnalyser cut = analyser.cut(*tr);

namespace CPlantBox {

			// cut.crop(&box); // cut with bounding box

namespace CPlantBox {

			// cut.write(os.str()+".vtp");

namespace CPlantBox {



namespace CPlantBox {

			// split into grid

namespace CPlantBox {

			std::cout <<"grid \n";

namespace CPlantBox {

			//			vector<vector<SegmentAnalyser>> anamatrix1 = cut.distribution2(0,160, 1.5*dist1, 2.5*dist1,n,m);

namespace CPlantBox {

			//			vector<vector<SegmentAnalyser>> anamatrix2 = cut.distribution2(0,160, 2.5*dist1, 3.5*dist1,n,m);

namespace CPlantBox {

			vector<vector<SegmentAnalyser>> anamatrix = cut.distribution2(0,160, 1.5*dist1, 3.5*dist1,n,2*m); // grid of 8 * 32

namespace CPlantBox {



namespace CPlantBox {

			// save root count into matrix

namespace CPlantBox {

			std::cout<<"count \n";

namespace CPlantBox {

			for (int i=0; i<n; i++) {

namespace CPlantBox {

				for (int j=0; j<m; j++) {

namespace CPlantBox {

					//					finalmatrix[i][j] += anamatrix[i][j].segments.size();

namespace CPlantBox {

					//					finalmatrix[i][j] += anamatrix[i][j+m].segments.size();

namespace CPlantBox {

					finalmatrix[i][j] += anamatrix[i][j].getSummed(RootSystem::st_length);

namespace CPlantBox {

					finalmatrix[i][j] += anamatrix[i][j+m].getSummed(RootSystem::st_length);;

namespace CPlantBox {

				}

namespace CPlantBox {

			}

namespace CPlantBox {

		}

namespace CPlantBox {



namespace CPlantBox {

		finalMT.push_back(finalmatrix); // a lot of finals we have here (one matrix per time step)

namespace CPlantBox {

	}

namespace CPlantBox {



namespace CPlantBox {

	/*

namespace CPlantBox {

	 * Export the final matrix (one file per matrix)

namespace CPlantBox {

	 */

namespace CPlantBox {

	double N = double(trenches.size())*2.*(4.5*5.);

namespace CPlantBox {

	// to calculate the mean number of roots per cm^2 over the trenches

namespace CPlantBox {

	for (int j=0; j<m; j++) {

namespace CPlantBox {

		string matname = name+"_trench_matrix"+std::to_string(j+1)+".txt"; // grid number

namespace CPlantBox {

		std::ofstream fos;

namespace CPlantBox {

		fos.open(matname.c_str());

namespace CPlantBox {

		for (int i=0; i<n ; i++) {

namespace CPlantBox {

			for (size_t t=0; t<times.size()-1; t++) {

namespace CPlantBox {

				std::cout << std::fixed << std::setprecision(4)<< finalMT.at(t)[i][j]/N << "\t";

namespace CPlantBox {

				fos << std::fixed << std::setprecision(4)<< finalMT.at(t)[i][j]/N << "\t";

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

		std::cout << "\n\n";

namespace CPlantBox {

	}

namespace CPlantBox {



namespace CPlantBox {

	/*

namespace CPlantBox {

	 * Export rootsystems (around 4.5GB)

namespace CPlantBox {

	 */

namespace CPlantBox {

	if (exportVTP) {

namespace CPlantBox {

		for (size_t i=0; i<times.size()-1; i++) {

namespace CPlantBox {

			string vtpname = name + std::to_string(i+1)+".vtp";

namespace CPlantBox {

			getResult(allRS,times.at(i+1)).write(vtpname);

namespace CPlantBox {

		}

namespace CPlantBox {



namespace CPlantBox {

		// export trench geometry

namespace CPlantBox {

		vector<SignedDistanceFunction*> tr_;

namespace CPlantBox {

		for (const auto& tr : trenches) { // copy vector

namespace CPlantBox {

			tr_.push_back(tr);

namespace CPlantBox {

		}

namespace CPlantBox {

		string gname = name + "_trench.py";

namespace CPlantBox {

		SDF_Union trenchgeometry = SDF_Union(tr_);

namespace CPlantBox {

		allRS[0]->setGeometry(&trenchgeometry); // just for writing

namespace CPlantBox {

		allRS[0]->write(gname);

namespace CPlantBox {



namespace CPlantBox {

//		allRS[0]->setGeometry(&box); // just for writing

namespace CPlantBox {

//		allRS[0]->write("boundingbox.py");

namespace CPlantBox {



namespace CPlantBox {

	}

namespace CPlantBox {



namespace CPlantBox {

}

namespace CPlantBox {



namespace CPlantBox {



namespace CPlantBox {

