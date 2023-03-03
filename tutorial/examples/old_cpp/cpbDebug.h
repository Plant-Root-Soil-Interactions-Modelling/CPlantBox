#include <iostream>
#include "Plant.h"
#include "mymath.h"
#include <cmath>

// sepcialized
#include "MappedOrganism.h"
#include "XylemFlux.h"

#include <iostream>
#include <vector>
#include "PiafMunch/runPM.h"
#include <limits>

std::shared_ptr<PhloemFlux> phloem_; 

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

void print_vector(std::vector<double> vec)
{
  std::cout << "size: " << vec.size() << std::endl;
  for (double d : vec)
    std::cout << d << " ";
  std::cout << std::endl;
}

void print_map(std::map<std::string, double>  m)
{
    for (auto const &pair: m) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}
 

double theta2H(std::vector<double> vg, double theta)//(-) to cm
{
    double thetar =vg[0];// 0.059
    double thetas = vg[1];//0.445
    double alpha = vg[2];//0.00644
    double n = vg[3];//1.503
    double nrev = 1/(1-1/n);
    double H =-(pow((pow( (thetas - thetar)/(theta - thetar), nrev)) - 1,(1/n)))/alpha;
    return H;
}

double sinusoidal(double t){return (sin(M_PI*t*2)+1)/2;}

double qair2rh(double qair, double es_, double press)
{
    double e =qair * press / (0.378 * qair + 0.622);
    double rh = e / es_;
    rh=max(min(rh, 1.),0.);
    return rh;
}

std::map<std::string, double> weather(double simDuration)
{
    std::vector<double> vgSoil{0.059, 0.45, 0.00644, 1.503, 1.};
    
    double Qmin = 0; double Qmax = 960e-6 ;
    double Tmin = 15.8; double Tmax = 22;
    double specificHumidity = 0.0097;
    double Pair = 1010.00 ;//hPa
    double thetaInit = 30./100.;

    double coefhours = sinusoidal(simDuration);
    double TairC_ = Tmin + (Tmax - Tmin) * coefhours;
    double Q_ = Qmin + (Qmax - Qmin) * coefhours;
    double cs = 850e-6; //co2 paartial pressure at leaf surface (mol mol-1)
    //RH = 0.5 # relative humidity
    double es =  6.112 * exp((17.67 * TairC_)/(TairC_ + 243.5));
    double RH = qair2rh(specificHumidity, es, Pair);
    
    double pmean = theta2H(vgSoil, thetaInit);
    std::map<std::string, double> weatherVar;
	weatherVar["TairC"] 	= TairC_;
	weatherVar["Qlight"] 	= Q_;
	weatherVar["cs"] 		= cs;
	weatherVar["RH"] 		= RH;
	weatherVar["p_mean"] 	= pmean;
	//weatherVar["vg"] 		= vgSoil;
    std::cout<<"Env variables at "<<floor(simDuration)<<"d ";
	std::cout<<floor((simDuration- floor(simDuration))*24)<<"hrs :\n";
	print_map(weatherVar);
    assert(std::isfinite(pmean)&&"pmean not a number");
    return weatherVar;
}


void setKrKx_xylem(double TairC, double RH)
{
    double hPa2cm = 1/0.9806806;
    //number of vascular bundles
    double VascBundle_leaf = 32.;
    double VascBundle_stem = 52.;
    double VascBundle_root = 1. ;//valid for all root type
	
	double dEauPure = (999.83952 + TairC * (16.952577 + TairC * 
		(- 0.0079905127 + TairC * (- 0.000046241757 + TairC * 
		(0.00000010584601 + TairC * (- 0.00000000028103006)))))) /  (1 + 0.016887236 * TairC);
    double siPhi = (30. - TairC) / (91. + TairC);
    double siEnne=0;
    double mu =  pow(10., (- 0.114 + (siPhi * (1.1 + 43.1 * pow(siEnne, 1.25) )))) ;
    mu = mu /(24.*60.*60.)/100./1000.; //mPa s to hPa d, 1.11837e-10 hPa d for pure water at 293.15K
    mu = mu * hPa2cm ;//hPa d to cmh2o d 

            
    //radius of xylem type^4 * number per bundle
    double rad_x_l_1   = (pow(0.0015 ,4.)) * 2.; double rad_x_l_2   = pow(0.0005 ,4.) * 2.  ; 
    double rad_x_s_1   = pow(0.0017 ,4.) * 3.; double rad_x_s_2   = pow(0.0008 ,4.) * 1.   ;  
    double rad_x_r0_1  = pow(0.0015 ,4.) * 4.    ;
    double rad_x_r12_1 = pow(0.00041,4.) * 4.; double rad_x_r12_2 = pow(0.00087,4.) * 1.;
    //double rad_x_r3_1  = pow(0.00068,4.) * 1.      ;
	
	
    //axial conductivity [cm^3/day]        
    double kz_l  = VascBundle_leaf *(rad_x_l_1 + rad_x_l_2)    *M_PI /(mu * 8.)  ;
    double kz_s  = VascBundle_stem *(rad_x_s_1 + rad_x_s_2)    *M_PI /(mu * 8.) ;
    double kz_r0 = VascBundle_root * rad_x_r0_1                *M_PI /(mu * 8.)  ;
    double kz_r1 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*M_PI /(mu * 8.) ;
    double kz_r2 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*M_PI /(mu * 8.)  ;
    //double kz_r3 = VascBundle_root * rad_x_r3_1                *M_PI /(mu * 8.) ;// 4.32e-1

    //radial conductivity [1/day],
    double kr_l  = 3.83e-4 * hPa2cm;// init: 3.83e-4 cm/d/hPa
    double kr_s  = 0. ;//1.e-20  * hPa2cm ;// set to almost 0
    double kr_r0 = 6.37e-5 * hPa2cm ;
    double kr_r1 = 7.9e-5  * hPa2cm ;
    double kr_r2 = 7.9e-5  * hPa2cm ; 
    //double kr_r3 = 6.8e-5  * hPa2cm ;
	double l_kr = 0.8 ;
	
	std::vector<std::vector<double>> kr_x{
		{kr_r0,kr_r1,kr_r2,kr_r0},{kr_s,kr_s},{kr_l}};
	std::vector<std::vector<double>> kz_x{
		{kz_r0,kz_r1,kz_r2,kz_r0},{kz_s,kz_s },{kz_l}};
	phloem_->setKr(kr_x, std::vector<std::vector<double>>(0), l_kr);
	phloem_->setKx(kz_x, std::vector<std::vector<double>>(0));
	
	double Rgaz=8.314 ;//#J K-1 mol-1 = cm^3*MPa/K/mol
    double rho_h2o = dEauPure/1000;//#g/cm3
    double Mh2o = 18.05;// #g/mol
    double MPa2hPa = 10000;
    //log(-) * (cm^3*MPa/K/mol) * (K) *(g/cm3)/ (g/mol) * (hPa/MPa) * (cm/hPa) =  cm                      
    double p_a = log(RH) * Rgaz * rho_h2o * (TairC + 273.15)/Mh2o * MPa2hPa * hPa2cm;

    phloem_->psi_air = p_a;//*MPa2hPa #used only with xylem
	
}
//sudo make install
//export CPLUS_INCLUDE_PATH=/home/rbtlm640/CPB_PMKleaf/report3/PMinCPB/src
// export LD_LIBRARY_PATH=/home/rbtlm640/CPB_PMKleaf/report3/PMinCPB/src
//g++ -O2 -pg -Wall -I. -o pmtry launchPM.cpp -lPMinCPB1
//gprof pmtry gmon.out > analysis.txt
//gprof pmtry gmon.out | less > analysis.txt

namespace CPlantBox {
void tryLaunch(double dt = 1./24.)//, double maxT = -1.) 
{
	
    //auto plant = std::make_shared<CPlantBox::PhloemFlux>();
	auto plant = std::make_shared<CPlantBox::MappedPlant>(2.);
    std::string path = "/home/rbtlm640/CPB2506/CPlantBox/modelparameter/plant/";
    std::string name = "Triticum_aestivum_adapted_2021.xml";//"testGrmax.xml";//
	plant->readParameters(path+name);
	
	//""" soil """
	double depth=60.;// double layers = 5.; // Soil core analysis
	double Inf = std::numeric_limits<double>::max();
	auto sdf = std::make_shared<CPlantBox::SDF_PlantBox>(Inf, Inf, depth );
	plant->setGeometry(sdf); // creates soil space to stop roots from growing out of the soil
	
    plant->initialize(true, false);
	double simInit = 7.;
	double simDuration = simInit;
	double simMax = simInit + 7.;
	
	plant->simulate(simInit, false);
	
	
	
	while(simDuration < simMax){
		plant->simulate(dt, false);
		simDuration += dt;
	}
	
	
}
}

//in source folder:
//cmake. or to have debug : cmake -DCMAKE_BUILD_TYPE=Debug .
//sudo make install
//installed chared PiafMunch library using cmake install

//g++ -o pmtry launchPM.cpp -lPIAFMunch
//g++ -o pmtry launchPM.cpp -lPIAFMunch
//g++ -o pmtry launchPM.cpp -lCPlantBox
//need to change the name of the two "Main"files

//to debug:
//g++ -g launchPM.cpp -lPIAFMunch
//g++ -g launchPM.cpp -lPIAFMunch

//g++ -D decoupled -Wl,-V -o pmtry launchPM.cpp -lPIAFMunch0712
//g++ -D decoupled -o pmtry launchPM.cpp -lPIAFMunch0712

//valgrind --leak-check=full --log-file="valgrind_output.txt" ./test
//valgrind --leak-check=yes --log-file="valgrind_output.txt" ./a.out                                    
//or
//valgrind --leak-check=yes ./PIAFMunch
//*************************************************
//bin/not used

//gcc -fpic -shared launchPM.cpp -lPIAFMunch
//lib not foud: sudo /sbin/ldconfig -v
//For C++ includes use 
//g++ -o pmtry launchPM.cpp -lPIAFMunch

//export CPLUS_INCLUDE_PATH=/home/rbtlm640/CPB_PMKleaf/report3/PMinCPB/src
//export CPLUS_INCLUDE_PATH
//echo $CPLUS_INCLUDE_PATH

/*I've got profiling working with a custom kernel. The missing piece is the CONFIG_HIGH_RES_TIMERS=y kernel option. I got this hint through the notice on the pprof documentation which links to #13841 where the option is referenced. For anyone who is curios how to get this working on WSL2:

clone https://github.com/microsoft/WSL2-Linux-Kernel
add option CONFIG_HIGH_RES_TIMERS=y to Microsoft/config-wsl
follow README-Microsoft.WSL2 (notice a missing apt within the dependency installation)
copy vmlinux within the root folder of the repo over to windows
add %USERPROFILE%/.wslconfig which sets kernel to vmlinux copied earlier
see Release Notes for file format
reboot computer or restart wsl with wsl --shutdown
I will contact the WSL team to see if it's possible to get the flag into the official kernel.
*/