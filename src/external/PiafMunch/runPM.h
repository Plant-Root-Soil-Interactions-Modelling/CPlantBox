/*
* PiafMunch (v.2) -- Implementing and Solving the Munch model of phloem sap flow in higher plants
*
* Copyright (C) 2004-2019 INRA
*
* Author: A. Lacointe, UMR PIAF, Clermont-Ferrand, France
*
* File: Main.cpp
*
* This file is part of PiafMunch. PiafMunch is free software: you can redistribute it and/or
* modify it under the terms of the GNU General Public License version 3.0 as published by
* the Free Software Foundation and appearing in the file LICENSE.GPL included in the
* packaging of this file. Please  review the following information to ensure the GNU
* General Public License version 3.0  requirements will be met:
* http://www.gnu.org/copyleft/gpl.html.
*
* PiafMunch is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
* without even the implied warranty of FITNESS FOR A PARTICULAR PURPOSE.
* See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with PiafMunch.
* If not, see <http://www.gnu.org/licenses/>.
*
-----------------------------------------------------------------------------------------------------------------------------------*/
#ifndef runPM_H_
#define runPM_H_
// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

#include "Plant.h"
#include "mymath.h"

//#include "PiafMunch2.h"
// sepcialized
#include "MappedOrganism.h"
#include "Photosynthesis.h"
#include "CPB_to_PM.h"


#include <math.h>
#include <vector>
#include <list>
#include <locale>
#include <algorithm>
#include "odepack.h"
//#include "PiafMunch2.h"
#include <fstream>

#include <stdio.h>
#include "PM_arrays.h"

#include <sundials/sundials_types.h>    /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>
#include <sundials/sundials_pcg.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_spbcgs.h>
#include <sundials/sundials_spfgmr.h>
#include <sundials/sundials_spgmr.h>
#include <sundials/sundials_sptfqmr.h>
#include <nvector/nvector_serial.h>     /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
#include <sunlinsol/sunlinsol_klu.h>    /* access to KLU sparse direct solver   */
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include <cvode/cvode.h>                /* prototypes for CVODE fcts., consts.  */
#include <cvode/cvode_bandpre.h>
#include <cvode/cvode_diag.h>
#include <cvode/cvode_direct.h>
#include <cvode/cvode_spils.h>
#include <cvode/cvode_impl.h>

#include <arkode/arkode.h>                /* prototypes for arkODE fcts., consts.  */
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_butcher.h>

//take those functions out
extern int Jac_(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
void TimeSegmentConfig() ; // User-editable (implemented in 'PiafMunch2.cpp')
void OutputSettings(bool doTroubleshooting_ = false) ;
void BreakpointSharpParameterChanges(int s, double t) ; // User-editable (implemented in 'PiafMunch2.cpp') ; s = # of integration segment (first = 1) ; t = time
void aux(double t, double * y) ;	


/**
 * Phloem flow based on PiafMunch code of LAcointe et al. 2019
 * see also phloem_flow.py in CPlantBox/src/python_modules
 *
 * Units are [hPa] and [day]
 * CplantBox object making link with PiafMunh
 * Wraps a Photosynthesis class
 */
class PhloemFlux: public CPlantBox::Photosynthesis, public std::enable_shared_from_this<PhloemFlux>
{
	public:
	PhloemFlux(std::shared_ptr<CPlantBox::MappedPlant> plant_, std::shared_ptr<CPlantBox::PlantHydraulicParameters> params,
	double psiXylInit = -500., double ciInit = 350e-6): 
		CPlantBox::Photosynthesis(plant_, params,psiXylInit, ciInit){};
	void initialize(){cpb_2_pm = std::make_shared<CPB_to_PM>(plant, shared_from_this());};
    std::shared_ptr<PhloemFlux> Phloem() { return shared_from_this(); }; // up-cast for Python binding
		
	virtual ~PhloemFlux() { }
	int startPM(double StartTime ,double EndTime, int OutputStep,double TairK, bool verbose = true , 
		std::string filename= "outpm.txt");///< main function called from python
	void computeOrgGrowth(double t);///< returns max sucrose need for growth per segment
	std::shared_ptr<CPB_to_PM> cpb_2_pm;								 
	


	//		for python post-processing and checks
    std::vector<int> orgTypes;
	std::vector<double> Q_Rmmax1v;
	std::vector<double> a_STv;
	vector<double> Q_initOutv;//Y0 vector at end of simulation
	vector<double> Q_init;//Y0 vector at beginning of simulation
	vector<double> Q_initOut_dotv;
	vector<double> C_STv;//Suc_ST for python byding
	vector<double> vol_STv;//void(*aux)(double,double*),
	vector<double> r_ST_refv; 
	vector<double> r_STv; 
	//length and volume incrase per segment and organs
	std::vector<double> delta_suc_org;
	std::vector<double> delta_ls;
	std::vector<double> delta_ls_org_i;
	std::vector<double> delta_ls_org;
	std::vector<double> delta_ls_org_imax;
	std::vector<double> delta_ls_org_max;
	std::vector<double> delta_vol_org_i, delta_vol_org, delta_vol_org_imax, delta_vol_org_max;
	std::vector<double> delta_vol_node_i, delta_vol_node, delta_vol_node_imax, delta_vol_node_max;
	
	//all in (mmol Suc d-1)
	std::vector<double> Agv;//assimilation (mmol Suc d-1)
	std::vector<double> Q_Grmaxv;//maximal sucrose sink for growth (mmol Suc d-1)
	std::vector<double> Q_Exudmaxv ;//maximal exudatoin rate, (mmol Suc d-1)
	std::vector<double> Q_Rmmaxv ;//maximal sucrose usage for maintenance, (mmol Suc d-1)
	std::vector<double> Flv ;//sucrose flow from mesophyll to sieve tube, (mmol Suc d-1)
	std::vector<double> vol_Mesov;//volume of mesophyll (same as leaf blade volume), (cm3)
	std::vector<double> JW_STv;//sieve tube water flow, (cm3 d-1)
	std::vector<double> JS_STv;//sieve tube sucrose flow, (mmol d-1)
	std::vector<double> Fpsi;//water scarcity factor for growth, (-)
	
	//		To calibrate
	double Q10 = 2.; double TrefQ10 = 20;//to compute effect of T on growth (see CN-wheat, residual respiration @Barillot 2016, appendix)
	double psiMax = 0; double psiMin = -2000*(1/0.9806806);//limit wat. pot. in xylem for water-limited growth, [cm]
	//double KMgr = 0.16; //@see C_fluxes,Michaelis menten coef for growth, not implemented
	double KMfu = 0.2; //@see C_fluxes,Michaelis menten coef for active sucrose usage
	//double k_meso = 1e-4;//conductivity if implement ohm analogy for Fl, not implemented
	double CsoilDefault =1e-4;//dummy value for soil concentration so that we always have ((Exud==0)||(Gr*Rm>0))
	std::vector<double> Csoil_seg;
	std::vector<double> Csoil_node;
	//used if sameVolume_meso_st == false, sameVolume_meso_seg == false
	double surfMeso =0.01 ;//cross sectinnal area of mesophyll (cm2). 
	//double Cobj_ST = 1.;// ==> when doing loading or unloading with objective value. not implemented
	double Vmaxloading = 0.019872;//mmol cm-1 d-1 for leaf blade 1cm wide
	double CSTimin = 0.4;//minimum CST value below which there is no sink of sucrose
	double beta_loading = 1;//@see C_fluxes, feedback effect of C_ST on Q_FL
	double Mloading = 0.2;//@see C_fluxes,Michaelis menten coef for Fl
	double Gr_Y = 0.75;//growth efficiency
	double k_S_ST = 0.;
	double C_targ = 0.4;
	double atol_double = 1e-017;//max absolute error
	double rtol_double = 1e-023;//max realtive error
	double initValST = 0.8;//initial concentration in sieve tube
	double initValMeso = 0.9;//initial concentration in mesophyll
	
	//		boolean choices
    bool doTroubleshooting = false;
	bool update_viscosity_ = true;
	bool usePsiXyl = true;//use PsiXyl value of xylem tissue
	bool sameVolume_meso_seg = true; //use same volume for mesophyll and leaf blade compartment?
	bool sameVolume_meso_st = true; //use same volume for mesophyll and leaf st compartment?
	bool withInitVal = false;//use initValST and initValMeso
	int solver = 1;//which solver to use
	int expression = 1;//if implement several possible expression in C_fluxes
	bool useStemTip = true;
	int growthType = 0; //0: standard; 1: threshold value for start; 2: 
	bool krFromLen = true;
    bool canStartActivating = true;
	double CSTthreshold = 0.3;
    bool leafAsIAASource = false;
    
    
	//		Auxin
    double StopAt_ = 0;
    double stopAt =-1;
	double auxin_threshold = 0.3;
	double auxin_P = 0.3;//production rate
	double auxin_D = 0.3;//decay
	double auxin_alpha = 0.3;//there is only one way down so we should be fine
    //Fortran_vector JAuxin_ST1;
    //Fortran_vector JAuxin_ST2;
	double initValAuxin = 0.9;//initial concentration in active tip
    std::vector<double> C_Auxinv;//Suc_ST for python byding
    std::vector<double> Delta_JA_STv;
    std::vector<double> C_AuxinOutv;
    std::vector<double> JAuxin_ST2v;
    bool deleteAtRootTip = false;
    std::vector<bool> isRootTip;
    bool burnInTime = false;
    std::vector<double> SucSTLost;
    std::vector<double> SucMesoLost;
    std::vector<double> AuxinLost;
    std::vector<double> manualAddST;
    std::vector<double> manualAddMeso;
    std::vector<double> manualAddAux;
    bool StopLoss = false;
    double L_dead_threshold = 2.;
    double auxin_init_mean;
    double BerthLim = -1;
    int useLength = 0;
    
	//internal PiafMunch functions but cannot protect
	void initialize_carbon(vector<double> vecIn) ;							// initializes carbon system parameters & constants (implemented in 'initialize.cpp')
	void initialize_hydric() ;							// initializes hydric system parameters & constants (implemented in 'initialize.cpp')
	void initializePM_(double dt,  double TairK); //copmutes PiafMunch input data from CPlantBox data
	void f(double t, double *y, double *y_dot) ;	//function launched by CVODE-SUNDIALS
	void aux(double t, double * y);
	void update_viscosity() ;
	void C_fluxes(double t, int Nt) ; // in  PiafMunch2.cpp
	
	protected:
	//internal parameters
	double TairK_phloem;//temperature in K for phloem tissue
	int Nt_old = 0; //BU old seg size
    Fortran_vector Q_GrowthtotBU ;
	Fortran_vector Q_GrmaxBU ;
	//bool hayErrores = false;
	int errorID = -1;
	int neq_coef = 10;//number of variables solved by PiafMunch. n# eq = num nodes * neq_coef
	
	
};



#endif