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

// ,"P_ST (MPa)",
	// "Resp_Maint_rate (mmol / d)","Resp_Maint (mmol)",
	// "Resp_Gr_rate (mmol / d)","Resp_Gr (mmol)",
	// "Exud_rate (mmol / d)","Exud (mmol)",
	// "Input (mmol / d)",
	// "Q_meso (mmol)", "Q_meso_dot (mmol/d)",
	 // "JS_ST (mmol / d)", "JW_ST (ml / d)"
// } ;

extern int Jac_(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
extern string Upper(string &str) ;
extern void LoadSettings(string ini_name) ;
//void Launch() ;
void TimeSegmentConfig() ; // User-editable (implemented in 'PiafMunch2.cpp')
void OutputSettings(string pos = "" ) ;
void BreakpointSharpParameterChanges(int s, double t) ; // User-editable (implemented in 'PiafMunch2.cpp') ; s = # of integration segment (first = 1) ; t = time
//void UpdateResistances(double t) ;
void aux(double t, double * y) ;	

//#ifndef decoupled
//
class PhloemFlux: public CPlantBox::Photosynthesis
{
	public:
	PhloemFlux(std::shared_ptr<CPlantBox::MappedPlant> plant_, double psiXylInit = -500., double ciInit = 350e-6): 
		CPlantBox::Photosynthesis(plant_, psiXylInit, ciInit){std::cout<<"new PM "<<std::endl;};//, plant(plant_)
    std::shared_ptr<PhloemFlux> Phloem() { return std::make_shared<PhloemFlux>(*this); }; // up-cast for Python binding
	void initialize_carbon(vector<double> vecIn) ;							// initializes carbon system parameters & constants (implemented in 'initialize.cpp')
	void initialize_hydric() ;							// initializes hydric system parameters & constants (implemented in 'initialize.cpp')
	void computeOrgGrowth(double t);
	double TairK_phloem;
	int Nt_old = 0; //BU old seg size
    virtual ~PhloemFlux() { }
	int startPM(double StartTime ,double EndTime, int OutputStep,double TairK,string dir_output = "outPM.txt", bool verbose = true );
	void Launch() ;
    //std::shared_ptr<CPlantBox::MappedPlant> plant;
	void initializePM_(double dt,  double TairK); 
	void prepareOutput(std::string dir_output) ;
	void f(double t, double *y, double *y_dot) ;
	//void *ff(double t, double *y, double *y_dot);//{f(t, y, y_dot);} ;
	void aux(double t, double * y);
	void update_viscosity() ;
	void C_fluxes(double t, int Nt) ; // in  PiafMunch2.cpp
	
	
	//		from plant shape
	std::vector<double> rhoSucrose{0,0,1.,0.8,0.6};//0.07*0.4*0.4/12/12 *1000; //mmol Sucrose /cm3 fresh matter 
    std::vector<std::vector<double>> kr_suc; //  [mmol Suc /hPa day-1]
    std::vector<std::vector<double>> kx_suc; // [mmol Suc day-1]
    std::vector<std::vector<double>> Across_STs; // [cm2]
	double kx_suc_f(int type, int organType) {return kx_suc.at(organType - 2).at(type); }//tot kx_st in segment per subtype and organ type (goes from 2 (root) to 4 (leaf))
    double kr_suc_f(int type, int organType) {return kr_suc.at(organType - 2).at(type); }//tot kz_st in segment per subtype and organ type (goes from 2 (root) to 4 (leaf))
    double get_Across_ST( int type, int organType){return Across_STs.at(organType - 2).at(type); }//cross-sectional area of all the sieve tubes in segment
	//std::vector<std::vector<double>> a_ST;
    //double get_a_ST( int type, int organType){return a_ST.at(organType - 2).at(type-1); }
	
	
	//		for python post-processing and checks
	vector<double> Q_outv;
	vector<double> Q_init;
	vector<double> Q_out_dotv;
	vector<double> C_STv;
	vector<double> vol_STv;//void(*aux)(double,double*),
	vector<double> r_ST_refv; 
	vector<double> r_STv; 
	std::vector<double> delta_suc_org;
	std::vector<double> delta_ls;
	std::vector<double> delta_ls_org_i;
	std::vector<double> delta_ls_org;
	std::vector<double> delta_ls_org_imax;
	std::vector<double> delta_ls_org_max;
	std::vector<double> delta_vol_org_i, delta_vol_org, delta_vol_org_imax, delta_vol_org_max;
	std::vector<double> delta_vol_node_i, delta_vol_node, delta_vol_node_imax, delta_vol_node_max;
	std::vector<double> Agv;
	std::vector<double> Q_Grmaxv;
	std::vector<double> Q_Exudmaxv ;
	std::vector<double> Q_Rmmaxv ;
	std::vector<double> Flv ;
	std::vector<double> vol_Mesov;
	std::vector<double> JW_STv;
	std::vector<double> a_STv;
	std::vector<std::map<int,double>> deltaVolOrgNode_;
	
	
	//		To calibrate
	bool update_viscosity_ = true;
	vector<double> krm1v = vector<double>(1, 2e-4*1000);
	vector<double> krm2v = vector<double> (1, 0.);
	vector<double> exud_kv = vector<double> (1, 0.1);
	double KMgr = 0.16;
	double KMrm = 0.2;//0.01;
	double k_gr = 1.;
	double k_meso = 1e-4;
	double Csoil =1e-4;//dummy value for soil concentration so that we always have ((Exud==0)||(Gr*Rm>0))
	double kout = 1e-4;
	double surfMeso =0.01 ;//cm2
	double Cobj_ST = 1.;
	double Vmaxloading = 0.019872;//mmol cm-1 d-1 for leaf blade 1cm wide
	double CSTi_objectiv = 1.;
	double CSTimin = 0.4;
	double beta_loading = 1;
	double Mloading = 0.2;
	int neq_coef = 9;
	int errorID = -1;
	int expression = 1;
	double k_gr2 =1.;
	double Gr_Y = 0.75;
	double atol_double = 1e-017;
	double rtol_double = 1e-023;
	double initValST = 0.8;
	double initValMeso = 0.9;
	//		boolean choices
	bool usePsiXyl = true;
	bool sameVolume_meso_st = true;
	bool withInitVal = false;
	bool sameVolume_meso_st = true;
	bool withInitVal = false;
	int solver = 1;
	
	protected
	//internal parameters
	Fortran_vector Q_ParApoBU ;
	Fortran_vector Q_GrmaxBU ;
	bool hayErrores = false;
};

//#else
//: public XylemFlux
/*
namespace CPlantBox {
class PhloemFlux
{
	public:
	//PhloemFlux(std::shared_ptr<CPlantBox::MappedPlant> plant_):XylemFlux(std::shared_ptr<CPlantBox::MappedSegments>(plant_)), plant(plant_){cout<<"alive creation "<<endl;};
	//PhloemFlux():XylemFlux(std::shared_ptr<CPlantBox::MappedSegments>(NULL)){cout<<"alive creation empty"<<endl;};
	PhloemFlux(){cout<<"alive creation "<<endl;};

    virtual ~PhloemFlux() { }
	int startPM(string post = "");
	void Launch() ;
	void initializePM_();
    double kst_soil; 
    double kst; 
    double kmeso; 
};
}
#endif
*/



#endif