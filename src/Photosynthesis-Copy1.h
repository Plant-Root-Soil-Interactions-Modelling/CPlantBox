// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef PHOTOSYNTHESIS_H_
#define PHOTOSYNTHESIS_H_

#include "MappedOrganism.h"
#include "XylemFlux.h"
#include <map>
#include <iostream>
#include <fstream>
namespace CPlantBox {

/**
 * Photosynthesis computation according to Tuzet (2003), Dewar 2002 and Lu 2020
 * see also photosynthesis.py in CPlantBox/src/python_modules
 *
 * Units are [cm] and [day]
 *
 * Wraps a MappedSegments class (i.e. MappedPlant)
 */
class Photosynthesis: public XylemFlux
{
public:
    Photosynthesis(std::shared_ptr<CPlantBox::MappedPlant> plant_, double psiXylInit = -500., double ciInit = 350e-6);
	virtual ~Photosynthesis() {}
	
    std::shared_ptr<CPlantBox::MappedPlant> plant;
	std::vector<std::map<int,double>> waterLimitedGrowth(double t);
	void solve_photosynthesis(double ea_,double es_, double sim_time_ =1., std::vector<double> sxx_= std::vector<double>(1,-200.0) , 
				bool cells_ = true, std::vector<double> soil_k_ = std::vector<double>(),
				bool doLog_ = false, int verbose_ = 0, 
				double TairC_=25, std::string outputDir_=""); ///< main function, makes looping until convergence
	void linearSystemSolve(double simTime_, const std::vector<double>& sxx_, bool cells_, 
				const std::vector<double> soil_k_);///< main function, solves the flux equations
	
	void loopCalcs(double simTime); ///<solves photosynthesis/stomatal opening equations 
	void getAg4Phloem(); ///< Converts An [mol CO2 m-2 s-1] to Ag4Phloem [mmol Suc d-1]
	void getError(double simTime);///< Computes error % for each segment for each of the variables of interestes.
	
	void doAddGravity(); ///< add gravitational wat. pot to total wat. pot. (used in phloem module)
	//void r_forPhloem(double lightTimeRatio, int ot);
    void toFile(std::string fileName,  std::vector<double>& descriptorsValues);
	
	//		intermediary results and outputs
    std::vector<double>  tauv, fv,  lengths, deltagco2, dv; 
	double Ko, Kc, delta, J, Rd;
	double es,ea;// 
	std::vector<double> psiXyl; //saves the wat. pot. values of xylem for photosynthesis and phloem modules
	std::vector<double> Ag4Phloem;//gross assimilation rate per segment for phloem module [mmol Suc d-1]
	std::vector<double> An; //gross assimilation rate assimilation per unit of surface [mol CO2 m-2 s-1]
	std::vector<double> Vc; //gross assimilation rate per unit of surface [mol CO2 m-2 s-1]
	std::vector<double> Vcmax; //gross assimilation rate per unit of surface [mol CO2 m-2 s-1]
	std::vector<double> Vj; //gross assimilation rate per unit of surface [mol CO2 m-2 s-1]
	std::vector<double> ci;
	std::vector<double> Jw;
	std::vector<double> Ev;
	std::vector<double> PVD;
	std::vector<double> EAL;
	std::vector<double> hrelL;
	std::vector<double> pg;//leaf guard cell water potential [cm]
	std::vector<double> outputFlux;
	std::vector<double> fw;
    std::vector<double> fwmesophyll;
	std::vector<double> psiXyl4Phloem; //sum of psiXyl + gravitational wat. pot.
	std::vector<double> gco2;
	std::vector<double> gtotOx;
    std::vector<double> kr_meso;
    //std::vector<double> k_stomatas;
    std::vector<double> Vcrefmax;std::vector<double> Jrefmax;
    float gm = 0.05; //mesophyll resistance 
    
    
	//		to evaluate convergence, @see Photosynthesis::getError
	bool doLog = false; int verbose_photosynthesis = 0;
    bool oldciEq = true;bool followTrace = false;
    bool usePg4Fw = true; //keep false or it does not converge I think//
    bool DoSun2012 = false; bool doWiggle = false;
    bool alternativeAn = false;
    bool useVc = true; bool useVj = true;
	int maxLoop = 1000; int minLoop = 1;std::string outputDir="";
    int loop;double limMaxErr = 1e-4;
	double maxMaxErr;
	std::vector<double> maxErr= std::vector<double>(9, 0.);	
	std::vector<double> outputFlux_old;
	std::vector<double> psiXyl_old;
	std::vector<double> An_old ;
	std::vector<double> gco2_old ;
	std::vector<double> ci_old ;
	std::vector<double> pg_old ;
	std::vector<double> k_stomatas_old ;
	std::vector<double> fw_old ;
	std::vector<double> fw_very_old ;
	
	
	//			initial guesses
	double psiXylInit; //initial guess for xylem total wat. pot. [cm]
	double ciInit; //initial guess for leaf internal [CO2] [mol mol-1]
	
	//			default value of env variable, to re-set at runtime
	double Patm = 1013.15;//default [hPa]
	double cs = 350e-6; //example from Dewar2002 [mol mol-1]
	double TleafK; //[K]
	double TairC = 20; //[°C]
	double Qlight = 900e-6;//mean absorbed photon irradiance per leaf segment [mol photons m-2 s-1]  
	std::vector<double>  Chl = std::vector<double>(1,55.); 
	double oi = 210e-3;//leaf internal [O2] [mol mol-1]
	double g_bl = 2.8;//leaf boundary molar conductance [mol CO2 m-2 s-1]
	double g_canopy = 6.1;//aerodynamic molar conductance [mol CO2 m-2 s-1]
	double g_air = 11.4;//aerodynamic molar conductance [mol CO2 m-2 s-1]
	
	//			parameter to re-parametrise , put in phloem files
	//water stress factor, parametrised from data of Corso2020
    double fwr = 9.308e-2; //residual opening when water stress parametrised with data from corso2020 [-]
	double sh = 4e-4;//3.765e-4;//sensibility to water stress
	double p_lcrit = -15000/2;//-0.869;//min psiXil for stomatal opening [Mpa]
    double fwrmesophyll = 9.308e-2; //residual opening when water stress parametrised with data from corso2020 [-]
	double shmesophyll = 4e-4;//3.765e-4;//sensibility to water stress
	double p_lcritmesophyll = -15000/2;//-0.869;//min psiXil for stomatal opening [Mpa]
	//influence of N contant, to reparametrise!, 
	double VcmaxrefChl1 = 1.28/2;//otherwise value too high, original: 1.31
	double VcmaxrefChl2 = 8.33/2; //otherwise value too high, original: 8,52
	//for Vc, Vj, Ag, to reparametrise!
    double Theta_ag = 0.98;//empirical parameter that governs the transition between Vc, Vj
	double a1=4.; //g0+ fw[i] * a1 *( An[i] + Rd)/(ci[i] - deltagco2[i]);//tuzet2003
	double a3 = 1.7;//Jrefmax = Vcrefmax * a3 ;//Eq 25
	double alpha = 0.2; //or 0.44 , coefb = -(alpha * Qlight + Jmax);, alpha * Qlight * Jmax;
	double theta = 0.9;//or 0.67 coefa = theta;
	double gamma0 = 28e-6; double gamma1 = 0.0509; double gamma2 = 0.001;
	double g0 = 0.3e-3;//residual stomatal opening to CO2, Tuzet 2003 [mol CO2 m-2 s-1]
	
	//			physicall constant (no need to parametrise)
	double a2_stomata = 1.6; //gco2 * a2 = gh2o
	double a2_bl = 1.37;//gco2 * a2 = gh2o
	double a2_canopy = 1;//gco2 * a2 = gh2o
	double a2_air = 1;//gco2 * a2 = gh2o
	//(de)activate parameters
    double Eac = 59430; double Eaj = 37000; double Eao = 36000;//mJ mmol-1
    double Eard = 53000;double Eav = 58520;
    double Edj =  220000;double Edv = 220000;
	//ref value at T = T_ref
	double Kc_ref = 302e-6; double Ko_ref = 256e-3; //mmol mmol-1
	double Rd_ref = 0.32e-6;double Tref = 293.2;//K
	double S = 700;//enthropy mJ mmol-1 K-1
	double Mh2o = 18;//molar weight, g mol-1, mg mmol-1
	double R_ph = 83.143;//perfect gas constant, hPa cm3K−1mmol−1
	double rho_h2o = 1000;//water density mg cm-3
	double M_Chla = 893.51;//chlorophyle a molar mass (g / mol)
    
	//void setKr_meso(std::vector<double> values); ///< sets a callback for kx_suc:=kx_suc(ot,type), 
    //std::function<double(int)> kr_meso_f = [](int type) {
		//throw std::runtime_error("kr_meso_f not implemented"); 
		//return 0.; };
    
protected:
	//Compute variables which do not vary during one "solve_photosynthesis " computation
	void initCalcs(double sim_time_);
	void initStruct(double sim_time_);
	void initVcVjRd();
	
	
	//			other parameters and runtime variables
	std::vector<std::shared_ptr<Organ>> orgsVec;
	std::vector<int> seg_leaves_idx;
    bool stop = false;
    
    
	//double kr_meso_const(  int type) { return kr_meso.at(0); }  //constant
    //double kr_meso_perOrgType( int type) { return kr_meso.at(type ); } //per organ type (goes from 2 (root) to 4 (leaf))
    

};

} // namespace

#endif