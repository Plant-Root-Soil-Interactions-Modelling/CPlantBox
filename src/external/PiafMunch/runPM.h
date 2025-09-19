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

//take those functions out
extern int Jac_(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
void TimeSegmentConfig() ; // User-editable (implemented in 'PiafMunch2.cpp')
void OutputSettings() ;
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
	PhloemFlux(std::shared_ptr<CPlantBox::MappedPlant> plant_,std::shared_ptr<CPlantBox::PlantHydraulicParameters> params,
	double psiXylInit = -500., double ciInit = 350e-6): 
		CPlantBox::Photosynthesis(plant_, params, psiXylInit, ciInit){};
    std::weak_ptr<PhloemFlux> Phloem() { return shared_from_this(); }; // up-cast for Python binding
	virtual ~PhloemFlux() { }
	int startPM(double StartTime ,double EndTime, int OutputStep,double TairK, bool verbose = true , 
		std::string filename= "outpm.txt");///< main function called from python
	void computeOrgGrowth(double t);///< returns max sucrose need for growth per segment
	
	//		from plant shape	
	void setKr_st(std::vector<std::vector<double>> values, double kr_length_, bool verbose = false); ///< sets a callback for kr_suc:=kr_suc(ot,type), 
	void setKx_st(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type), 
	void setRmax_st(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type), 
	void setAcross_st(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type), 
	void setPerimeter_st(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type),  
	void setRhoSucrose(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type),  
	void setKrm1(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type), 
	void setKrm2(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type), 
	
    std::function<double(int, int, int)> kr_st_f = []( int type, int orgtype, int si) {
		throw std::runtime_error("kr_st_f not implemented"); 
		return 0.; };
    std::function<double(int,int)> kx_st_f = [](int type, int orgtype) {
		throw std::runtime_error("kx_st_f not implemented"); 
		return 0.; };
    std::function<double(int, int)> Across_st_f = [](int type, int orgtype) {//cross-sectional area of all the sieve tubes in segment
		throw std::runtime_error("get_Across_st not implemented"); 
		return 0.; };
    std::function<double(int, int)> Perimeter_st_f = [](int type, int orgtype) {//cross-sectional area of all the sieve tubes in segment
		throw std::runtime_error("get_Across_st not implemented"); 
		return 0.; };
    std::function<double(int,int)> Rmax_st_f = []( int type, int orgtype){//maximum initial growth rate use by phloem module
		throw std::runtime_error("get_Rmax not implemented"); 
		return 0.; };
    std::function<double(int,int)> rhoSucrose_f = []( int type, int orgtype){//maximum initial growth rate use by phloem module
		throw std::runtime_error("get_rhoSucrose not implemented"); 
		return 0.; };
    std::function<double(int,int)> krm1_f = []( int type, int orgtype){//maximum initial growth rate use by phloem module
		throw std::runtime_error("get_rhoSucrose not implemented"); 
		return 0.; };
    std::function<double(int,int)> krm2_f = []( int type, int orgtype){//maximum initial growth rate use by phloem module
		throw std::runtime_error("get_rhoSucrose not implemented"); 
		return 0.; };


	//		for python post-processing and checks
	std::vector<double> a_STv;
	vector<double> Q_outv;//Y0 vector at end of simulation
	vector<double> Q_init;//Y0 vector at beginning of simulation
	vector<double> Q_out_dotv;
	vector<double> C_STv;//Suc_ST for python byding
	vector<double> vol_STv;//void(*aux)(double,double*),
	vector<double> r_ST_refv; 
	vector<double> r_STv; 
	
	//all in (mmol Suc d-1)
	std::vector<double> Agv;//assimilation (mmol Suc d-1)
	std::vector<double> Q_Grmaxv;//maximal sucrose sink for growth (mmol Suc d-1)
	std::vector<double> Q_Exudmaxv ;//maximal exudatoin rate, (mmol Suc d-1)
	std::vector<double> Q_Rmmaxv ;//maximal sucrose usage for maintenance, (mmol Suc d-1)
	std::vector<double> Flv ;//sucrose flow from mesophyll to sieve tube, (mmol Suc d-1)
	std::vector<double> vol_Mesov;//volume of mesophyll (same as leaf blade volume), (cm3)
	std::vector<double> JW_STv;//sieve tube water flow, (cm3 d-1)
    //for growth
	std::vector<double> Fpsi;//water scarcity factor for growth, (-)
    std::vector<double> Flen;
    std::vector<int> GrowthZone;
    std::vector<int> GrowthZoneLat;
	std::vector<std::map<int,double>> deltaSucOrgNode_;//maximal sucrose need for growth per node, (mmol Suc d-1)
    std::vector<double> psi_p_symplasm;
	
	// intermediate outputs
	bool do_gf_warning = true; // get a warning in case we have an unexpected growth function (gives warning only once)
		
	//		To calibrate	
	// soil (mmol Suc cm-3)
	std::vector<double> Csoil_seg;
	std::vector<double> Csoil_node;
	double CsoilDefault = 1e-4;//dummy value for soil concentration so that we always have ((Exud==0)||(Gr*Rm>0))
	
	// initial values
	double initValST = 0.8;//initial concentration in sieve tube
	double initValMeso = 0.9;//initial concentration in mesophyll
	bool withInitVal = false;//use initValST and initValMeso
	
	// growth
	double psi_osmo_proto = -4000*1.0197;
    double psiMin = 2000*1.0197;//limit wat. pot. in xylem for water-limited growth, [cm]
    double leafGrowthZone = 2;//cm
	double Gr_Y = 0.75;//growth efficiency
    bool StemGrowthPerPhytomer = true;
	bool useCWGr; //use water- and carbon- limited growth?
		
	// sieve tube
	double Vmaxloading = 0.019872;//mmol cm-1 d-1 for leaf blade 1cm wide
	double CSTimin = 0.4;//minimum CST value below which there is no sink of sucrose
	double beta_loading = 1;//@see C_fluxes, feedback effect of C_ST on Q_FL
	double Mloading = 0.2;//@see C_fluxes,Michaelis menten coef for Fl	
	double C_targ = 0;//(mmol Suc cm-3 )
	double Q10 = 2.; double TrefQ10 = 20;// effect of T on respiration (see CN-wheat, residual respiration @Barillot 2016, appendix)
	double KMfu = 0.2; //@see C_fluxes,Michaelis menten coef for active sucrose usage
	double k_mucil = 0;//(d-1)
	std::vector<double> k_mucil_;
	double Vmax_S_ST = 0;//(mmol Suc d-1 cm-3 )
	double kM_S_ST = 0;//(mmol Suc cm-3)
	double kHyd_S_ST = 0;//(d-1)
	double k_S_ST = 0;//(d-1)
	bool update_viscosity_ = true;
	bool usePsiXyl = true;//use PsiXyl value of xylem tissue
	
	// mesophyll 
	double C_targMesophyll = 0;//(mmol Suc cm-3 )
	double Vmax_S_Mesophyll = 0;//(mmol Suc d-1 cm-3 )
	double kM_S_Mesophyll = 0;//(mmol Suc cm-3)
	double kHyd_S_Mesophyll = 0;//(d-1)
	double k_S_Mesophyll = 0;//(d-1)
	double surfMeso =0.01 ;//cross sectinnal area of mesophyll (cm2). 
	bool sameVolume_meso_seg = true; //use same volume for mesophyll and leaf blade compartment?
	bool sameVolume_meso_st = true; //use same volume for mesophyll and leaf st compartment?
	
	// per type
    std::vector<std::vector<double>> kr_st;//  [mmol hPa-1 day-1]
	std::vector<std::vector<double>> kx_st; //  [cm3 hPa-1 day-1]
    std::vector<std::vector<double>> Across_st; // [cm2]
    std::vector<std::vector<double>> Perimeter_st; // [cm]
	std::vector<std::vector<double>> Rmax_st; // [cm day-1]
	std::vector<std::vector<double>> rhoSucrose;
	std::vector<std::vector<double>> krm1v; 
	std::vector<std::vector<double>> krm2v;
	
	// solver and other
	double atol_double = 1e-017;//max absolute error
	double rtol_double = 1e-023;//max realtive error
	int solver = 1;//which solver to use
	bool doTroubleshooting =false; //do extra printing
	
	// not implemented																								  
	//used if sameVolume_meso_st == false, sameVolume_meso_seg == false
	//double Cobj_ST = 1.;// ==> when doing loading or unloading with objective value. not implemented
	//double k_meso = 1e-4;//conductivity if implement ohm analogy for Fl, not implemented
	//double KMgr = 0.16; //@see C_fluxes,Michaelis menten coef for growth, not implemented
	
    	//internal PiafMunch functions but cannot protect
	void initialize_carbon(vector<double> vecIn) ;							// initializes carbon system parameters and constants (implemented in 'initialize.cpp')
	void initialize_hydric() ;							// initializes hydric system parameters and constants (implemented in 'initialize.cpp')
	void initializePM_(double dt,  double TairK); //copmutes PiafMunch input data from CPlantBox data
	void f(double t, double *y, double *y_dot) ;	//function launched by CVODE-SUNDIALS
	void aux(double t, double * y);
	void update_viscosity() ;
	void C_fluxes(double t, int Nt) ; // in  PiafMunch2.cpp
	
	std::vector<std::map<int,double>> waterLimitedGrowth(double t);
	double adaptDt(std::shared_ptr<CPlantBox::Organ> org, double dt);
	void assertUsedCReserves(std::shared_ptr<CPlantBox::Organ> org);
	double getMaxVolumicGrowth(std::shared_ptr<CPlantBox::Organ> org, double t);
	std::vector<int> getGrowingNodes(std::shared_ptr<CPlantBox::Organ> org);
	void computeFpsi();
	void computeFlen(std::shared_ptr<CPlantBox::Organ> org, std::vector<int> growingNodesId);
	
	protected:
	//internal parameters
	double TairK_phloem;//temperature in K for phloem tissue
	std::size_t Nt_old = 0; //BU old seg size
    Fortran_vector Q_GrowthtotBU ;
	Fortran_vector Q_GrmaxBU ;
	//bool hayErrores = false;
	int errorID = -1;
	std::size_t neq_coef = 10;//number of variables solved by PiafMunch. n# eq = num nodes * neq_coef
	std::map<int, double> BackUpMaxGrowth;//to check at runtime if growth is correct. is a map because max(orgID) > len(orgs)
	
	//retrieve tissue-specific parameters
	double kr_st_const(  int type, int organType, int si) { return kr_st.at(0).at(0); }  //constant
    double kr_st_perOrgType( int type, int organType, int si) { return kr_st.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double kr_st_perType( int type, int organType, int si) {return kr_st.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	double kr_st_RootExchangeZonePerType(int type, int organType, int si)
	{ 
		if (organType == CPlantBox::Organism::ot_root){
			double coef = plant->exchangeZoneCoefs.at(si);//% of segment length in the root exchange zone, see MappedPlant::simulate
			return coef * kr_st.at(organType - 2).at(type); 
		}
		return kr_st.at(organType - 2).at(type);  
	} //subtype, type and depend on distance to tip for roots
	
	double kx_st_const( int type, int organType) { return kx_st.at(0).at(0); }  //constant
    double kx_st_perOrgType( int type, int organType) { return kx_st.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double kx_st_perType( int type, int organType) {return kx_st.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	
	double Across_st_const( int type, int organType) { return Across_st.at(0).at(0); }  //constant
    double Across_st_perOrgType( int type, int organType) { return Across_st.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double Across_st_perType( int type, int organType) {return Across_st.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
    
	double Perimeter_st_const( int type, int organType) { return Perimeter_st.at(0).at(0); }  //constant
    double Perimeter_st_perOrgType( int type, int organType) { return Perimeter_st.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double Perimeter_st_perType( int type, int organType) {return Perimeter_st.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	
	double Rmax_st_const( int type, int organType) { return Rmax_st.at(0).at(0); }  //constant
    double Rmax_st_perOrgType( int type, int organType) { return Rmax_st.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double Rmax_st_perType( int type, int organType) {return Rmax_st.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	
	double rhoSucrose_const( int type, int organType) { return rhoSucrose.at(0).at(0); }  //constant
    double rhoSucrose_perOrgType( int type, int organType) { return rhoSucrose.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double rhoSucrose_perType( int type, int organType) {return rhoSucrose.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	
	double krm1_const( int type, int organType) { return krm1v.at(0).at(0); }  //constant
    double krm1_perOrgType( int type, int organType) { return krm1v.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double krm1_perType( int type, int organType) {return krm1v.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	
	double krm2_const( int type, int organType) { return krm2v.at(0).at(0); }  //constant
    double krm2_perOrgType( int type, int organType) { return krm2v.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double krm2_perType( int type, int organType) {return krm2v.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	
	
};



#endif