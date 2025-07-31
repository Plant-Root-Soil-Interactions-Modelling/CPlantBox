// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef PHOTOSYNTHESIS_H_
#define PHOTOSYNTHESIS_H_

#include "MappedOrganism.h"
#include "PlantHydraulicModel.h"
#include <map>
#include <iostream>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace CPlantBox {

/**
 * Photosynthesis computation according to Tuzet (2003), Dewar 2002 and Lu 2020
 * see also photosynthesis.py in CPlantBox/src/python_modules
 *
 * Units are [cm] and [day]
 *
 * Wraps a MappedSegments class (i.e. MappedPlant)
 */
class Photosynthesis: public PlantHydraulicModel
{
public:
	enum PhotoTypes { C3 = 0, C4 = 1}; ///< add other later (CAM?)
    Photosynthesis(std::shared_ptr<CPlantBox::MappedPlant> plant_, std::shared_ptr<CPlantBox::PlantHydraulicParameters> params, double psiXylInit = -500., double ciInit = 350e-6);
	virtual ~Photosynthesis() {}
	
    std::shared_ptr<CPlantBox::MappedPlant> plant;
	std::vector<std::map<int,double>> waterLimitedGrowth(double t);
	void solve_photosynthesis(double sim_time_ , std::vector<double> sxx_, 
				double ea_,double es_, std::vector<double> TleafK_, 
				bool cells_ = true, std::vector<double> soil_k_ = std::vector<double>(),
				bool doLog_ = false, int verbose_ = 0, 
				std::string outputDir_=""); ///< main function, makes looping until convergence
	void linearSystemSolve(double simTime_, const std::vector<double>& sxx_, bool cells_, 
				const std::vector<double> soil_k_);///< main function, solves the flux equations
	
	void loopCalcs(double simTime, std::vector<double> sxx_, bool cells_); ///<solves photosynthesis/stomatal opening equations 
	void photoC4_loop(int i);
	void photoC3_loop(int i);
	//Compute variables which do not vary during one "solve_photosynthesis " computation
	void initCalcs(double sim_time_);
	void initStruct(double sim_time_);
	void initVcVjRd();
	void photoC4_init(int i);
	void photoC3_init(int i);
	
	void getAg4Phloem(); ///< Converts An [mol CO2 m-2 s-1] to Ag4Phloem [mmol Suc d-1]
	void getError(double simTime);///< Computes error % for each segment for each of the variables of interestes.
	
	void doAddGravity(); ///< add gravitational wat. pot to total wat. pot. (used in phloem module)
	
	double Q10f(int index){		return std::pow(Q10_photo,(getMeanOrSegData(TleafK, index) - Tref)/10);		}
	double thermalBreakdown(int index, double Ed);
	double Arrhenius(int index, double Ea);
	
	double getMeanOrSegData(std::vector<double> data, int index)
	{
		if(data.size() != seg_leaves_idx.size()){return data.at(0);
		}else{return data.at(index);}
	}
	
	
	double getPsiOut(bool cells, int si, const std::vector<double>& sx_) const override;
	size_t fillVectors(size_t k, int i, int j, double bi, double cii, double cij, double psi_s) override ; ///< fill the vectors aI, aJ, aV, aB
	double kr_f(int si, double age, int type, int orgtype);
		
		
	//void r_forPhloem(double lightTimeRatio, int ot);
	
	//					intermediary results and outputs
	// 				C3 and C4
    std::vector<double>  tauv, fv,  lengths, deltagco2, dv; 
	double es,ea;// 
	std::vector<double> delta;
	std::vector<double> Rd_ref;
	std::vector<double> Rd;
	std::vector<double> psiXyl; //saves the wat. pot. values of xylem for photosynthesis and phloem modules
	std::vector<double> Ag4Phloem;//gross assimilation rate per segment for phloem module [mmol Suc d-1]
	std::vector<double> An; //gross assimilation rate assimilation per unit of surface [mol CO2 m-2 s-1]
	std::vector<double> Vc; //gross assimilation rate per unit of surface [mol CO2 m-2 s-1]
	std::vector<double> Vcmax; //gross assimilation rate per unit of surface [mol CO2 m-2 s-1]
    std::vector<double> Vcrefmax;
	std::vector<double> Vj; //gross assimilation rate per unit of surface [mol CO2 m-2 s-1]
	std::vector<double> ci;
	std::vector<double> Jw;//transpiration [cm3 cm-2 d-1]
	std::vector<double> Ev; //transpiration [cm3/day]
	std::vector<double> PVD;
	std::vector<double> EAL;
	std::vector<double> hrelL;
	std::vector<double> pg;//leaf guard cell water potential [cm]
	std::vector<double> outputFlux;
	std::vector<double> fw;
	std::vector<double> psiXyl4Phloem; //sum of psiXyl + gravitational wat. pot.
	std::vector<double> gco2;
	std::vector<double> gtotOx;
	// 				C3 only
	std::vector<double> Jrefmax;
	std::vector<double> J;
	std::vector<double> Jmax;
	// 				C3 MOSTLY (very small effect on C4)
	std::vector<double> Ko;
	std::vector<double> Kc;
	// 				C4 only
	std::vector<double> Vp; //PEP carboxylase-limited (or CO2-limited) rate [mol CO2 m-2 s-1]
	std::vector<double> kp; //initial slope of the CO2 response curve for Vp [mol CO2 m-2 s-1]
	std::vector<double> kp25; //initial slope of the CO2 response curve for Vp [mol CO2 m-2 s-1] at the reference temperature
	
	//		to evaluate convergence, @see Photosynthesis::getError
	int maxLoop = 50000; int minLoop = 1;std::string outputDir="";
	bool doLog = false; int verbose_photosynthesis = 0;
    int loop;									   
    std::vector<double> maxErrAbs= std::vector<double>(9, 0.);	
	std::vector<double> maxErr= std::vector<double>(9, 0.);	
    //psi_x, An, gco2, ci, pg, outputFlux, 0, sum(F), sum(An)
    std::vector<double> maxErrAbsLim = { 1.,1e-6 , 0.002, 1e-6, 1,1,1e-6,1,1e-6 };	
	std::vector<double> maxErrLim  =  std::vector<double>(9, 1e-4);  	
    bool canStop();
	std::vector<double> outputFlux_old;
	std::vector<double> psiXyl_old;
	std::vector<double> An_old ;
	std::vector<double> gco2_old ;
	std::vector<double> ci_old ;
	std::vector<double> pg_old ;
	
	//			initial guesses
	double psiXylInit; //initial guess for xylem total wat. pot. [cm]
	double ciInit; //initial guess for leaf internal [CO2] [mol mol-1]
	
	//___________
	//			env variable, potentially re-set at each time step
	double Patm = 1013.15;//default [hPa]
	std::vector<double>  cs = std::vector<double>(1, 350e-6); //example from Dewar2002 [mol mol-1]
	std::vector<double>  TleafK; //[K]
	std::vector<double>  Qlight;//mean absorbed photon irradiance per leaf segment [mol photons m-2 s-1]  
	std::vector<double>  g_bl = std::vector<double>(1,2.8);//leaf boundary molar conductance [mol CO2 m-2 s-1]
	std::vector<double>  g_canopy = std::vector<double>(1,6.1);//aerodynamic molar conductance [mol CO2 m-2 s-1]
	std::vector<double>  g_air = std::vector<double>(1,11.4);//aerodynamic molar conductance [mol CO2 m-2 s-1]
	//___________
	
	
	//_______
	//			plant parameters to re-parametrise 
	// 				C3 and C4
	int PhotoType = C3;	
    float gm = 0.03; //mesophyll resistance 
	//water stress factor, parametrised from data of Corso2020
    double fwr = 9.308e-2; //residual opening when water stress parametrised with data from corso2020 [-]
	double fw_cutoff = 0;// to make it easier to get fw
	double sh = 3.765e-4;//sensibility to water stress
	double p_lcrit = -15000/2;//min psiXil for stomatal opening [Mpa]
	//influence of N contant, to reparametrise!, 
	std::vector<double>  Chl = std::vector<double>(1,55.); // mug/cm^-2
	double VcmaxrefChl1 = 1.28/2;//otherwise value too high, original: 1.31
	double VcmaxrefChl2 = 8.33/2; //otherwise value too high, original: 8,52
	double alpha = 0.2; //or 0.44 , coefb = -(alpha * Qlight + Jmax);, alpha * Qlight * Jmax;
	double a1=4.; //g0+ fw[i] * a1 *( An[i] + Rd)/(ci[i] - deltagco2[i]);//tuzet2003
	double g0 = 0.3e-3;//residual stomatal opening to CO2, Tuzet 2003 [mol CO2 m-2 s-1]	
	// 				C3 only
	double a3 = 1.7;//Jrefmax = Vcrefmax * a3 ;//Eq 25
	double theta = 0.9;//or 0.67 coefa = theta;
	double oi = 210e-3;//leaf internal [O2] [mol mol-1]	
	// 				C4 only
	double Q10_photo = 2;
	//___________
	
	
	//___________
	//			physical constant (no need to reparametrise?)
	// 				C3 and C4
	double a2_stomata = 1.6; //gco2 * a2 = gh2o
	double a2_bl = 1.37;//gco2 * a2 = gh2o
	double a2_canopy = 1;//gco2 * a2 = gh2o
	double a2_air = 1;//gco2 * a2 = gh2o
	double Mh2o = 18;//molar weight, g mol-1, mg mmol-1
	double R_ph = 83.143;//perfect gas constant, hPa cm3K−1mmol−1
	double rho_h2o = 1000;//water density mg cm-3
	double M_Chla = 893.51;//chlorophyle a molar mass (g / mol)
	double Tref = 293.2;//K
	// 				MOSTLY C3 (very small effect on C4)
	//(de)activate parameters
    double Eac = 59430; double Eaj = 37000; double Eao = 36000;//mJ mmol-1
    double Eard = 53000;double Eav = 58520;double Ead = 37830;
    double Edj =  220000;double Edv = 220000;double Edrd;// = 53000;
	double S = 700;//enthropy mJ mmol-1 K-1
	//ref value at T = T_ref
	double Kc_ref = 302e-6; double Ko_ref = 256e-3; //mmol mmol-1
    double delta_ref= 42.75e-6;	
	// 				C4 only
	double s1 = 0.3; //K-1	Bonan2019Chap11: temperature
	double s2 = 313.15;//K	Bonan2019Chap11
	double s3 = 0.2; //K-1	Bonan2019Chap11
	double s4 = 288.15;//K	Bonan2019Chap11
	double s5 = 1.3; //K-1	Bonan2019Chap11
	double s6 = 328.15;//K	Bonan2019Chap11
	//___________
	
	
	//			other parameters and runtime variables
	std::vector<std::shared_ptr<Organ>> orgsVec;
	std::vector<int> seg_leaves_idx;
    bool stop = false;
protected:
	//for Photosynthesis::linearSystemSolve
	std::vector<Eigen::Triplet<double>> tripletList;
	Eigen::VectorXd b;
	

};

} // namespace

#endif