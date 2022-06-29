/*
* PiafMunch (v.2) -- Implementing and Solving the Munch model of phloem sap flow in higher plants
*
* Copyright (C) 2004-2019 INRA
*
* Author: A. Lacointe, UMR PIAF, Clermont-Ferrand, France
*
* File: PiafMunch2.cpp -- main Model Code file for PiafMunch v.2
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
#include "runPM.h"
//#include "PiafMunch2.h"
string GnuplotPath = "C:/gnuplot/bin/gnuplot.exe" ; // Path to GnuPlot executable, including executable file name
//Fortran_vector vol_ST_seg;
//Fortran_vector C_ST_seg, Q_ST_seg, Q_ST_seg_new, Q_ST_newFuFl, grad_Q_ST_seg;
Fortran_vector Length;
Fortran_vector  Q_Exud, Q_Gr, Q_Rm, Q_Fl; 
Fortran_vector  Ag, Q_Grmax, Q_Rmmax, Q_Exudmax, exud_k, krm2, len_leaf;


// Architecture :
extern int N1L ; // number of nodes of type # 1L  'leaf end' : conn.order = 1 ; has imposed transpiration stream as a limit condition
extern int N1R ; // number of nodes of type # 1R 'root end' : conn.order = 1 ; has imposed soil water potential as a limit condition
extern int Nt, Nc, N[] ; // N[o = 1..8] : number of nodes of conn.order o (N[1] = N1L + N1R) ; Nt : Tot. nb of nodes = N[1] + N[2] + N[3] + ... + N[8] ; Nc : tot. nb of connectors
extern Index_vector i_[] ; // i_[o = 0..8][k = 1..No] = id# of kth node in (1-based list)  list of all nodes of conn.ord. co(o). For o=0: No=N1L ; o=1: No=N1R ; o=2..8: No=N[o] :
extern Index_vector &RootEnds, &LeafEnds ; // (= i_[1] and i_[0], resp.) : label which nodes are network ends = nodes of co=1 , which are either 'leaf' or 'root' ends :
extern vector<int> I_Upflow, I_Downflow ; // I_Upflow(resp.I_Downflow)[jf=1..Nc] = id# du noeud amont (resp. aval) : jf = JF(i,i2) > 0 si I_Upflow[(abs(jf)]==i, i.e. si I_Downflow[(abs(jf)]==i2
extern vector<Index_vector> Connect ; // Connect[i][j] = Id# of (j)th node (as ranked in input .ini file) connected to node #i (i = 1..Nt ; j = 1..co(i))
#ifndef co
#define co(i) Connect[i].size() // connectivity order of node# i (Id# i = 1..Nt)
#endif
extern int i, nbv ; extern bool LogScale ; // for output settings  in function ' parameter_and_boundary_conditions(t)'  at the end of this file
// Tracer-related parameters :
double TracerHalfLife ; // (h) set in GUI (default = 0.33967 h  for 13C)
double TracerDecay_k ; // (h-1) exponential tracer decay coeff = ln(2) / TracerHalfLife
// param. for sympl. volume changes, elastic or irreversible (Lockhart) :
Fortran_vector k_Lockhart, P_thr, vol_Sympl_max ;

int solver = 1 ; // (int. in the range [1,35] -- see Main.cpp) solver config.# : normally the smaller the better -- but may change in specific configs. (very complex architecture and/or equations...). Try diff. value if calc. fails or unstable or too slow !

/******************************************  Environmental Variables / boundary conditions : *********************************************/
// The following 4 imposed (but possibly changing) external conditions are set/updated by the user in function  'parameter_and_boundary_conditions()' below :
double T ; // (K) absolute temperature : set in GUI, but may be updated anytime
Fortran_vector Transpirat		; // (mmol / h)  Leaf transpiration rate
Fortran_vector PsiSoil			; // (MPa)  Soil water potential at root end
double * vol_Sympl_dot = NULL ; // (ml / h) Variation rate of lateral parenchyma symplasm volume -- to implement pressure-dependent volume variation (reversible elastic reservoir, irreversible Lockhart growth...)

/******************* VARIABLES INVOLVED in CARBON METABOLISM AND FLUXES *******************************  */
// as components of function f argument double* y in the solving procedure, the first 3 are stored as double* ; all others are Fortran_vectors
double* Q_out = NULL;
double* Q_Sympl = NULL						; // Amount of sugar in parenchyma symplasm										(mmol)
double* Q_ST = NULL						; // Amount of sugar in sieve tubes	= C_TC * Vol_ST						(mmol)
double* Starch = NULL						; // Amount of starch in parenchyma										(mmol)
double* Q_PhlApo = NULL                 ; // amount of sugar in phloem apoplasm   (mmol)
double* Q_ParApo = NULL                 ; // amount of sugar in lateral parenchyma apoplasm   (mmol)
Fortran_vector JS_ST						; // (mmol / h)  : Axial phloem sugar flux
Fortran_vector JS_PhlMb				; // Phloem cross-membrane sugar fluxes into sieve tubes from apoplasm				(mmol / h)
Fortran_vector RespMaint					; // Maintenance respiration rate										(mmol / h)
Fortran_vector StarchSyn						; // Michaelis-Menten rate of starch synthesis from sugar substrate						(mmol sug.eq./ h)
Fortran_vector C_amont					; //  (mmol / ml) : ST Sugar concentration at upflow node
Fortran_vector Input					; // e.g. Leaf photosynthetic assimilation rate, but may occur at any node, whether 'leaf' or not (boundary condition)			(mmol / h)
Fortran_vector C_Sympl						; // Concentration of sugar in parenchyma = Q_Par / vol_Sympl			(mmol / ml)
Fortran_vector C_ST							; // Concentration of sugar in sieve tubes								(mmol / ml solution))
double* vol_Sympl = NULL							; // total volume of symplasm	(ml), which is a reservoir of variable size
Fortran_vector JS_Sympl						; // Symplasmic flux of sugar from Lateral parenchyma to phloem ST (mmol / h)
Fortran_vector JS_ParMb				; // Lateral parenchyma cross-membrane sugar flux into symplasm from apoplasm				(mmol / h)
Fortran_vector JS_Apo						; // apoplasmic sugar flux from phloem to Lateral parenchyma  (mmol / h)
extern Fortran_vector JW_ParMb, JW_Apo, JW_Sympl ; // water fluxes corresponding to above 3 sugar fluxes   (ml / h)
Fortran_vector C_SymplUpflow					; // upflow concentration (mmol / ml) for JS_Sympl
Fortran_vector C_PhlApo, C_ParApo, C_ApoUpflow ; // (mmol / ml) apoplasmic sugar conc., resp. in phloem and lat.parenchyma, and upflow conc. for JS_Apo
Fortran_vector Delta_JS_ST ; // sera la composante purement phloémienne de Q_TC_dot[ ]							(mmol / h)
extern Fortran_vector P_Sympl		; // Lateral parenchyma symplasmic turgor pressure 								(MPa)
double Starch_dot_alt ; // for an alternate, target-oriented,  expression of starch variation rate

/******** TRACER-RELATED VARIABLES ***********************************************************/
double* TracerQ_Sympl = NULL						; // Amount of tracer in parenchyma										(MBq)
double* TracerQ_ST = NULL						; // Amount of soluble tracer in sieve tubes	= TracerC_ST * Vol_ST						(MBq)
double* TracerStarch = NULL						; // Amount of insoluble tracer in parenchyma										(MBq)
double *TracerQ_PhlApo=NULL, *TracerQ_ParApo=NULL ; // Amount of soluble tracer in apoplasm	= TracerC_xxxApo * Vol_xxxApo			(MBq)
Fortran_vector TracerJS_ST						; // TracerJS_ST[i] (MBq / h)  = TracerJS_ST_[i-1], i = 2..N  : Axial phloem soluble tracer flux ; TracerJS_ST[1]= NA, whereas TracerJS_ST_[1] = TracerJS_ST[2]
Fortran_vector TracerJS_PhlMb						; // Lateral (Apoplasmic) soluble tracer fluxes into sieve tubes from parenchyma				(MBq / h)
Fortran_vector TracerRespMaint					; // Tracer Maintenance respiration rate										(MBq / h)
Fortran_vector TracerStarchSyn						; // Michaelis-Menten rate of tracer starch synthesis from tracer sugar substrate						(MBq / h)
Fortran_vector TracerInput					; // Leaf photosynthetic Tracer Assimilation rate (boundary condition)			(MBq / h)
Fortran_vector TracerC_Sympl						; // Concentration of soluble tracer in parenchyma symplasm= TracerQ_Sympl / vol_Sympl			(MBq / ml)
Fortran_vector TracerC_ST							; // Concentration of soluble tracer in sieve tubes								(MBq / ml solution))
Fortran_vector TracerC_PhlApo, TracerC_ParApo ; // Concentration of soluble tracer in apoplasms					(MBq / ml solution))
Fortran_vector TracerRatioSympl ; //   TracerQ_Sympl / Q_Sympl = TracerC_Sympl / C_Sympl  (MBq / mmol)
Fortran_vector TracerRatioStarch ; //   = TracerStarch / Starch (MBq / mmol)
Fortran_vector Delta_TracerJS_ST ; // sera la composante purement phloémienne de TracerQ_ST_dot[ ]
Fortran_vector TracerJS_Sympl, TracerJS_Apo, TracerJS_ParMb ; // Lateral (Sympl., Apopl.; cross-membr...) soluble tracer fluxes FROM sieve tubes INTO parenchyma							(MBq / h)     !!! ATTENTION : sens positif opposé à JS_Trsv !!!
Fortran_vector TracerC_SymplUpflow					;  // upflow tracer concentration (mmol / ml) for TracerJS_Sympl
Fortran_vector TracerC_ApoUpflow ;				; // upflow tracer concentration (mmol / ml) for TracerJS_Apo

/******* variation rate of any variable X is noted: X_dot = dX/dt : *******/
double * Q_Sympl_dot = NULL, * Q_ST_dot = NULL, * Starch_dot = NULL, *Q_PhlApo_dot = NULL, *Q_ParApo_dot = NULL,  *Q_out_dot = NULL  ;
double * TracerQ_Sympl_dot=NULL, * TracerQ_ST_dot=NULL, * TracerStarch_dot=NULL, *TracerQ_PhlApo_dot=NULL, *TracerQ_ParApo_dot=NULL ;
// Next 2 variables are not considered as such, but as possible inputs to compute vol_Sympl_dot :
Fortran_vector P_ST_dot, P_Sympl_dot			; //  dP_ST/dt , dP_Sympl/dt					(MPa h / h)    -- for elasticity...

/********************* C-FLUXES-RELATED BIOPHYSICAL & PHYSIOLOGICAL PARAMETERS **********************************/
Fortran_vector kML						; // kinetic parameter / Michaelis - phloem loading					(mmol / ml)
Fortran_vector vML					; // kinetic parameter / phloem loading								(mmol /h)
Fortran_vector kMU						; // kinetic parameter / Michaelis - phloem unloading					(mmol / ml)
Fortran_vector vMU					; // kinetic parameter / phloem unloading								(mmol /h)
Fortran_vector kMParMb			; // kinetic parameter / Michaelis - parenchyma crossmembrane C flux					(mmol / ml)
Fortran_vector vMParMb			; // kinetic parameter / Michaelis - parenchyma crossmembrane C flux					(mmol /h)
Fortran_vector kM						; // kinetic parameter / Michaelis - starch Synthesis					(mmol / ml)
Fortran_vector Vmax					; // kinetic parameter / starch Synthesis								(mmol ml-1 h-1)
Fortran_vector C_targ			  	; // kinetic parameter / starch/sugar equilibrium. (regul. par. sugar conc.) 			(mmol / ml)
Fortran_vector kHyd					; // kinetic parameter / starch hydrolysis								(h-1)
Fortran_vector k1					; // kinetic parameter / maintenance respiration 						(h-1)
Fortran_vector k2					; // kinetic parameter / maintenance respiration						(ml mmol-1 h-1)
Fortran_vector k3						; // kinetic parameter / starch/sugar equilibrium. (regul. par. sugar conc.) 			(h-1)
Fortran_vector StructC				; // structural C subject to maintenance respiration					(mmol sug. eq.)
Fortran_vector vol_ST, vol_ParApo, vol_PhlApo ; // (ml) total volume of, resp. : sieve tubes, and  parenchyma and phloem apoplasm
//Fortran_vector	isTip ,Ag,  Lmax_org, Rmax_org, krm1, krm2 , exud_k, Q_Grmax_; 
//double krm2;
/*******************************   LIQUID-FLUXES-RELATED BIOPHYSICAL PARAMETERS  **********************************/
extern Fortran_vector r_Xyl, r_Trsv, r_ST, r_ST_ref, r_Sympl, r_Apo, r_PhlMb, r_Par_Mb ; // Hydraulic resistances (might change in response to embolism, cold block, aquaporin function, etc...)
// constants to express hydro resistance changes :
#define NONE 0
#define XYL 1
#define TRSV 2
#define APO 4
#define SYMPL 8
#define PHL 16
#define PHLMB 32
#define PARMB 64
#define RABS 128
extern int resistance_changed ;
Fortran_vector  radius_ST ; // vol_Sympl is considered a variable, driven by its variation rate -- see function  Smooth_Parameter_and_BoundaryConditions_Changes()  below


void PhloemFlux::C_fluxes(double t, int Nt)  {
// May be edited, including redifinition of GUI-defined parameters (kML, vML, k1, k2, kM, vM, vMU, k3, etc...).
// Furthermore, additional parameters can be defined and used, after sizing and assigning as in following example :
// Fortran_vector AddParam(Nt, 2.) ; AddParam[2] = 3. ;  // in this example, AddParam is set to 2. for all nodes except for node #3
    //int vecI;
	//at 2020-10-13T12:02:17+02:00,urn:fzj:tereno:samples:SE_Y_046,SE_Y_046-L1-000200
	//Leach10cmSmpConcentrationTOC => 14,82 mg L-1 => 1,235 mmol C L => 1,235e-3 mmol C ml => 0,0001 mmol Suc ml
	
	
	for (int i = 1 ; i <= Nt ; i++) { // edit (make different loops) to enter specific equations for specific nodes or conn.orders
		//Input[i] = Ag[i];
		
		double CSTi = max(0.,C_ST[i]);// From A.Lacointe: solver may try C<0 even if actual C never does
		double CSTi_delta = max(0.,C_ST[i]-Csoil);
		//vecI = i;1.51357e-006
		double Cmeso = max(0.,Q_Sympl[i]/vol_ParApo[i]);
		Q_Fl[i] = k_meso*max(Cmeso - CSTi, 0.);//Ag[i];////Ag[i] * (1-1/(std::exp(10*(-CSTi+1))+1));
		 
		double Q_Rmmax_ ;double Q_Exudmax_;double Fu_lim;
		Q_out_dot[i] = 0;
		switch(expression) {
		case 0:{
			Q_Exud[i] = exud_k[i]*CSTi_delta*Q_Exudmax[i];//exud_k[i]*CSTi_delta*(2*3.14*Length[i]*(radius_ST[i]*1e2));
			Q_Gr[i] = (CSTi/(CSTi + KMgr)) * Q_Grmax[i];//
			Q_Rm[i] = (CSTi/(CSTi + KMrm)) * (Q_Rmmax[i] + krm2[i] * CSTi);//krm1[i]*StructC[i];//(krm2[i]*CSTi + krm1[i])*StructC[i];
			//Q_ST_dot:
			Q_ST_dot[i] = Q_Fl[i] - Q_Exud[i] - Q_Gr[i] - Q_Rm[i]+ Delta_JS_ST[i];			
			//Q_meso_dot:
			Q_Sympl_dot[i] = Ag[i] -Q_Fl[i];//0.;////Q_Gr[i];
			Input[i] = Q_Fl[i];
			//Resp_maint:
			Starch_dot[i] = Q_Rm[i];
			//Growth:
			Q_ParApo_dot[i] = Q_Gr[i];
			//Exudation:
			Q_PhlApo_dot[i] = Q_Exud[i] ;
			//Resp_maintmax:
			TracerQ_ST_dot[i] = Q_Rmmax[i] + krm2[i] * CSTi;//Q_Rm[i];
			//Growthmax:
			TracerQ_ParApo_dot[i] = Q_Grmax[i];//Q_Gr[i];
			//Exudationmax:
			TracerQ_PhlApo_dot[i] = Q_Exudmax_;}//Q_Exud[i] ;
			break;
		case 1:{
			//Q_Rmmax[i] = krm1 * StructSucrose
			Q_Rmmax_ = Q_Rmmax[i] + krm2[i] * CSTi;
			//std::cout<<i<<" "<<Q_Rmmax[i]<<" "<<krm2[i]<<std::endl;
			Q_Exudmax_ = exud_k[i]*CSTi_delta*Q_Exudmax[i];
			//std::cout<<exud_k[i]<<" "<<Q_Exudmax[i]<<" "<<Csoil<<std::endl;
			Fu_lim = (Q_Rmmax_  + Q_Grmax[i])* (CSTi/(CSTi + KMrm));
			Q_ST_dot[i] = Q_Fl[i] - Fu_lim -Q_Exudmax_ + Delta_JS_ST[i];
			//Q_meso_dot:
			Q_Sympl_dot[i] = Ag[i] -Q_Fl[i];//0.;////Q_Gr[i];
			Input[i] = Q_Fl[i];
			//Resp_maint:
			Starch_dot[i] = min(Fu_lim, Q_Rmmax_);//Q_Rm[i];
			//Growth:
			Q_ParApo_dot[i] = min(Fu_lim - Starch_dot[i], Q_Grmax[i]);//Q_Gr[i];
			//Exudation:
			Q_PhlApo_dot[i] =  Q_Exudmax_;//min(Fu_lim - Starch_dot[i] - Q_ParApo_dot[i], Q_Exudmax_);//Q_Exud[i] ;
			
			
			//Resp_maintmax:
			TracerQ_ST_dot[i] = Q_Rmmax_ +krm2[i] * CSTi;//Q_Rm[i];
			//Growthmax:
			TracerQ_ParApo_dot[i] = Q_Grmax[i];//Q_Gr[i];
			//Exudationmax:
			TracerQ_PhlApo_dot[i] = Q_Exudmax_;//Q_Exud[i] ;
			// if(C_ST[i] > 3){
			//std::cout<<"test "<<i<<" "<<C_ST[i]<<" "<<CSTi<<std::endl;//<<" "<<(CSTi/(CSTi + KMrm))<<std::endl;
			 // std::cout<<Q_Sympl[i] <<" "<<vol_ST[i]<<" "<<Q_Fl[i]<<std::endl;
			 // std::cout<<Q_Sympl_dot[i] <<" "<<Ag[i]<<" "<<t<<std::endl;
			// assert(Q_Fl[i]==0&&"eend test");
			// std::cout<<(Q_Rmmax_ + Q_Exudmax_ + Q_Grmax[i])<<" "<<Fu_lim<<std::endl;
			// std::cout<<Q_Rmmax_<<" "<<Starch_dot[i]<<std::endl;
			// std::cout<<Q_Grmax[i]<<" "<<Q_ParApo_dot[i]<<std::endl;
			// std::cout<<Q_Exudmax_<<" "<<Q_PhlApo_dot[i]<<std::endl;
			// std::cout<<Delta_JS_ST[i]<<std::endl;
			// }
			if((Starch_dot[i]<0)||(Q_ParApo_dot[i]<0)||(Q_PhlApo_dot[i]<0)){
				std::cout<<"neg Sink "<<Starch_dot[i]<<" "<<Q_ParApo_dot[i]<<" "<<Q_PhlApo_dot[i]<<std::endl;
				std::cout<<Fu_lim<<std::endl;
				assert(false&&"neg Sink ");
			}
			if((Q_ST[i]<= 0) &&(Q_ST_dot[i] < 0)){
				std::ofstream outfile;
				outfile.open("errors.txt", std::ios_base::app); // append instead of overwrite
				outfile<< std::endl<<"C_fluxes "<<t<<" "<<i<<" "<<CSTi<<" "<<Q_ST[i]<<" qdot: "<<Q_ST_dot[i]<<", Fu: "<<Fu_lim <<std::flush;
				outfile<<" rm:"<<Starch_dot[i]<<" gr:"<<Q_ParApo_dot[i]<<" ex:"<<Q_PhlApo_dot[i]<<" js:"<<Delta_JS_ST[i]<<std::flush;
				outfile<<" vol_ST:"<<vol_ST[i]<<std::flush;
				outfile<<std::endl<<std::flush;
				errorID = i;
			}
			}break;
		case 2:{
			//Q_Rmmax[i] = krm1 * StructSucrose
			Q_Rmmax_ = Q_Rmmax[i] + krm2[i] * CSTi;
			Q_Exudmax_ = exud_k[i]*CSTi_delta*Q_Exudmax[i];
			Fu_lim = (Q_Rmmax_ + Q_Exudmax_ + Q_Grmax[i])* (CSTi/(CSTi + KMrm));
			Q_ST_dot[i] = Q_Fl[i] - Fu_lim + Delta_JS_ST[i];
			//Q_meso_dot:
			Q_Sympl_dot[i] = Ag[i] -Q_Fl[i];//0.;////Q_Gr[i];
			Input[i] = Q_Fl[i];
			//Resp_maint:
			Starch_dot[i] = min(Fu_lim, Q_Rmmax_);//Q_Rm[i];
			//Growth:
			Q_ParApo_dot[i] = min(Fu_lim - Starch_dot[i], Q_Grmax[i]);//Q_Gr[i];
			//Exudation:
			Q_PhlApo_dot[i] = min(Fu_lim - Starch_dot[i] - Q_ParApo_dot[i], Q_Exudmax_);//Q_Exud[i] ;
			
			
			//Resp_maintmax:
			TracerQ_ST_dot[i] = Q_Rmmax_;//Q_Rm[i];
			//Growthmax:
			TracerQ_ParApo_dot[i] = Q_Grmax[i];//Q_Gr[i];
			//Exudationmax:
			TracerQ_PhlApo_dot[i] = Q_Exudmax_;//Q_Exud[i] ;
			}break;
		case 3:{
		
			//Q_Rmmax[i] = krm1 * StructSucrose
			double Couti = Q_out[i]/vol_ST[i];
			double Couti_delta = max(0.,Couti-Csoil);
			Q_Rmmax_ = Q_Rmmax[i] + krm2[i] * Couti;
			Q_Exudmax_ = exud_k[i]*Couti_delta*Q_Exudmax[i];
			Fu_lim = kout*(CSTi - Couti) ;
			double Fu_lim2 =(Q_Rmmax_  + Q_Grmax[i])* (Couti/(Couti + KMrm));
			Q_ST_dot[i] = Q_Fl[i] - Fu_lim + Delta_JS_ST[i];
			//Q_meso_dot:
			Q_Sympl_dot[i] = Ag[i] -Q_Fl[i];//0.;////Q_Gr[i];
			Input[i] = Q_Fl[i];
			//Resp_maint:
			Starch_dot[i] = min(Fu_lim2, Q_Rmmax_);//Q_Rm[i];
			//Growth:
			Q_ParApo_dot[i] = min(Fu_lim2 - Starch_dot[i], Q_Grmax[i]);//Q_Gr[i];
			//Exudation:
			Q_PhlApo_dot[i] =Q_Exudmax_;// min(Fu_lim - Starch_dot[i] - Q_ParApo_dot[i], Q_Exudmax_);//Q_Exud[i] ;
			//Qout
			Q_out_dot[i] = Fu_lim - Fu_lim2 - Q_Exudmax_;
			
			//Resp_maintmax:
			TracerQ_ST_dot[i] = Q_Rmmax_;//Q_Rm[i];
			//Growthmax:
			TracerQ_ParApo_dot[i] = Q_Grmax[i];//Q_Gr[i];
			//Exudationmax:
			TracerQ_PhlApo_dot[i] = Q_Exudmax_;//Q_Exud[i] ;
			}break;
			
		case 4:{
			// from https://www.frontiersin.org/articles/10.3389/fpls.2022.787837/full
			
			Q_Fl[i] = (Vmaxloading *len_leaf[i])* Cmeso/(Mloading + Cmeso) * exp(-CSTi* beta_loading);//k_meso*min(max(Cobj_ST - CSTi, 0.), Cmeso);//
			Q_Rmmax_ = Q_Rmmax[i] + krm2[i] * CSTi;
			//std::cout<<i<<" "<<Q_Rmmax[i]<<" "<<krm2[i]<<std::endl;
			Q_Exudmax_ = exud_k[i]*CSTi_delta*Q_Exudmax[i];
			//std::cout<<exud_k[i]<<" "<<Q_Exudmax[i]<<" "<<Csoil<<std::endl;
			Fu_lim = (Q_Rmmax_  + Q_Grmax[i])* (CSTi/(CSTi + KMrm));
			Q_ST_dot[i] = Q_Fl[i] - Fu_lim -Q_Exudmax_ + Delta_JS_ST[i];
			//Q_meso_dot:
			Q_Sympl_dot[i] = Ag[i] -Q_Fl[i];//0.;////Q_Gr[i];
			Input[i] = Q_Fl[i];
			//Resp_maint:
			Starch_dot[i] = min(Fu_lim, Q_Rmmax_);//Q_Rm[i];
			//Growth:
			//add max(X,0.) in case of issues with rounding
			Q_ParApo_dot[i] = max(min(Fu_lim - Starch_dot[i], Q_Grmax[i]),0.);//Q_Gr[i];
			//Exudation:
			Q_PhlApo_dot[i] =  Q_Exudmax_;//min(Fu_lim - Starch_dot[i] - Q_ParApo_dot[i], Q_Exudmax_);//Q_Exud[i] ;
			
			
			//Resp_maintmax:
			TracerQ_ST_dot[i] = Q_Rmmax_ ;//+krm2[i] * CSTi;//Q_Rm[i];
			//Growthmax:
			TracerQ_ParApo_dot[i] = Q_Grmax[i];//Q_Gr[i];
			//Exudationmax:
			TracerQ_PhlApo_dot[i] = Q_Exudmax_;//Q_Exud[i] ;
			// if((Starch_dot[i]<0)||(Q_ParApo_dot[i]<0)||(Q_PhlApo_dot[i]<0)){
				// std::cout<<"neg Sink "<<Starch_dot[i]<<" "<<Q_ParApo_dot[i]<<" "<<Q_PhlApo_dot[i]<<std::endl;
				// std::cout<<Fu_lim<<std::endl;
				// assert(false&&"neg Sink ");
			// }
			if(((Q_ST[i]<= 0) &&(Q_ST_dot[i] < 0))||((Starch_dot[i]<0)||(Q_ParApo_dot[i]<0)||(Q_PhlApo_dot[i]<0))){
				std::ofstream outfile;
				outfile.open("errors.txt", std::ios_base::app); // append instead of overwrite
				outfile<< std::endl<<"C_fluxes "<<t<<" "<<i<<" "<<CSTi<<" "<<Q_ST[i]<<" qdot: "<<Q_ST_dot[i]<<", Fu: "<<Fu_lim <<std::flush;
				outfile<<" rm:"<<Starch_dot[i]<<" maxrm: "<<TracerQ_ST_dot[i]<<std::flush;
				outfile<<" gr:"<<Q_ParApo_dot[i]<<" grmax: "<<TracerQ_ParApo_dot[i]<<" exud:"<<Q_PhlApo_dot[i]<<" js:"<<Delta_JS_ST[i]<<std::flush;
				outfile<<" vol_ST:"<<vol_ST[i]<<std::flush;
				outfile<<std::endl<<std::flush;
				errorID = i;
				assert(false&&"neg in piafmunch, expression 4");
			}
			}break;
		case 5:{
			double Vmaxloading_ = max(0.,CSTi_objectiv-CSTi) *vol_ST[i];//in mmol
			Q_Fl[i] = Vmaxloading_ *Cmeso/(Mloading + Cmeso);// * exp(-CSTi* beta_loading);//k_meso*min(max(Cobj_ST - CSTi, 0.), Cmeso);//
			Q_Rmmax_ = Q_Rmmax[i] + krm2[i] * CSTi;
			//std::cout<<i<<" "<<Q_Rmmax[i]<<" "<<krm2[i]<<std::endl;
			Q_Exudmax_ = exud_k[i]*CSTi_delta*Q_Exudmax[i];
			//std::cout<<exud_k[i]<<" "<<Q_Exudmax[i]<<" "<<Csoil<<std::endl;
			Fu_lim = (Q_Rmmax_  + Q_Grmax[i])* (CSTi/(CSTi + KMrm));
			Q_ST_dot[i] = Q_Fl[i] - Fu_lim -Q_Exudmax_ + Delta_JS_ST[i];
			//Q_meso_dot:
			Q_Sympl_dot[i] = Ag[i] -Q_Fl[i];//0.;////Q_Gr[i];
			Input[i] = Q_Fl[i];
			//Resp_maint:
			Starch_dot[i] = min(Fu_lim, Q_Rmmax_);//Q_Rm[i];
			//Growth:
			//add max(X,0.) in case of issues with rounding
			Q_ParApo_dot[i] = max(min(Fu_lim - Starch_dot[i], Q_Grmax[i]),0.);//Q_Gr[i];
			//Exudation:
			Q_PhlApo_dot[i] =  Q_Exudmax_;//min(Fu_lim - Starch_dot[i] - Q_ParApo_dot[i], Q_Exudmax_);//Q_Exud[i] ;
			
			
			//Resp_maintmax:
			TracerQ_ST_dot[i] = Q_Rmmax_ ;//+krm2[i] * CSTi;//Q_Rm[i];
			//Growthmax:
			TracerQ_ParApo_dot[i] = Q_Grmax[i];//Q_Gr[i];
			//Exudationmax:
			TracerQ_PhlApo_dot[i] = Q_Exudmax_;//Q_Exud[i] ;
			// if((Starch_dot[i]<0)||(Q_ParApo_dot[i]<0)||(Q_PhlApo_dot[i]<0)){
				// std::cout<<"neg Sink "<<Starch_dot[i]<<" "<<Q_ParApo_dot[i]<<" "<<Q_PhlApo_dot[i]<<std::endl;
				// std::cout<<Fu_lim<<std::endl;
				// assert(false&&"neg Sink ");
			// }
			if((Q_Fl[i]<0)||((Q_ST[i]<= 0) &&(Q_ST_dot[i] < 0))||((Starch_dot[i]<0)||(Q_ParApo_dot[i]<0)||(Q_PhlApo_dot[i]<0))){
				std::ofstream outfile;
				outfile.open("errors.txt", std::ios_base::app); // append instead of overwrite
				outfile<< std::endl<<"C_fluxes "<<t<<" "<<i<<" "<<CSTi<<" "<<Q_ST[i]<<" qdot: "<<Q_ST_dot[i]<<", Fu: "<<Fu_lim <<std::flush;
				outfile<<" rm:"<<Starch_dot[i]<<" maxrm: "<<TracerQ_ST_dot[i]<<std::flush;
				outfile<<" gr:"<<Q_ParApo_dot[i]<<" grmax: "<<TracerQ_ParApo_dot[i]<<" exud:"<<Q_PhlApo_dot[i]<<" js:"<<Delta_JS_ST[i]<<std::flush;
				outfile<<" vol_ST:"<<vol_ST[i]<<" fl "<<Q_Fl[i]<<" vmaxloading "<<Vmaxloading_<<std::flush;
				outfile<<std::endl<<std::flush;
				errorID = i;
				assert(false&&"neg in piafmunch, expression 5");
			}
			}break;
		case 6:{
			Q_Fl[i] = (Vmaxloading *len_leaf[i])* Cmeso/(Mloading + Cmeso) * exp(-CSTi* beta_loading);//k_meso*min(max(Cobj_ST - CSTi, 0.), Cmeso);//
			CSTi = max(0., CSTi-CSTimin);
			Q_Rmmax_ = Q_Rmmax[i] + krm2[i] * CSTi;
			//std::cout<<i<<" "<<Q_Rmmax[i]<<" "<<krm2[i]<<std::endl;
			Q_Exudmax_ = exud_k[i]*CSTi_delta*Q_Exudmax[i];
			//std::cout<<exud_k[i]<<" "<<Q_Exudmax[i]<<" "<<Csoil<<std::endl;
			Fu_lim = (Q_Rmmax_  + Q_Grmax[i])* (CSTi/(CSTi + KMrm));
			Q_ST_dot[i] = Q_Fl[i] - Fu_lim -Q_Exudmax_ + Delta_JS_ST[i];
			//Q_meso_dot:
			Q_Sympl_dot[i] = Ag[i] -Q_Fl[i];//0.;////Q_Gr[i];
			Input[i] = Q_Fl[i];
			//Resp_maint:
			Starch_dot[i] = min(Fu_lim, Q_Rmmax_);//Q_Rm[i];
			//Growth:
			//add max(X,0.) in case of issues with rounding
			Q_ParApo_dot[i] = max(min(Fu_lim - Starch_dot[i], Q_Grmax[i]),0.);//Q_Gr[i];
			//Exudation:
			Q_PhlApo_dot[i] =  Q_Exudmax_;//min(Fu_lim - Starch_dot[i] - Q_ParApo_dot[i], Q_Exudmax_);//Q_Exud[i] ;
			
			
			//Resp_maintmax:
			TracerQ_ST_dot[i] = Q_Rmmax_ ;//+krm2[i] * CSTi;//Q_Rm[i];
			//Growthmax:
			TracerQ_ParApo_dot[i] = Q_Grmax[i];//Q_Gr[i];
			//Exudationmax:
			TracerQ_PhlApo_dot[i] = Q_Exudmax_;//Q_Exud[i] ;
			// if((Starch_dot[i]<0)||(Q_ParApo_dot[i]<0)||(Q_PhlApo_dot[i]<0)){
				// std::cout<<"neg Sink "<<Starch_dot[i]<<" "<<Q_ParApo_dot[i]<<" "<<Q_PhlApo_dot[i]<<std::endl;
				// std::cout<<Fu_lim<<std::endl;
				// assert(false&&"neg Sink ");
			// }
			if(((Q_ST[i]<= 0) &&(Q_ST_dot[i] < 0))||((Starch_dot[i]<0)||(Q_ParApo_dot[i]<0)||(Q_PhlApo_dot[i]<0))){
				std::ofstream outfile;
				outfile.open("errors.txt", std::ios_base::app); // append instead of overwrite
				outfile<< std::endl<<"C_fluxes "<<t<<" "<<i<<" "<<CSTi<<" "<<CSTimin<<" "<<Q_ST[i]<<" qdot: "<<Q_ST_dot[i]<<", Fu: "<<Fu_lim <<std::flush;
				outfile<<" rm:"<<Starch_dot[i]<<" maxrm: "<<TracerQ_ST_dot[i]<<std::flush;
				outfile<<" gr:"<<Q_ParApo_dot[i]<<" grmax: "<<TracerQ_ParApo_dot[i]<<" exud:"<<Q_PhlApo_dot[i]<<" js:"<<Delta_JS_ST[i]<<std::flush;
				outfile<<" vol_ST:"<<vol_ST[i]<<std::flush;
				outfile<<std::endl<<std::flush;
				errorID = i;
				assert(false&&"neg in piafmunch, expression 4");
			}
			}break;
             
		}
		// kgr = km2 * (((1-eps)/eps)**2)
		
		//checks:
			//assert((Ag[i]>=Q_Fl[i])&&"Rmmax == 0 outside of seed node");
			assert((Q_Rmmax[i] >= 0.)&&"Rmmax < 0 ");
			
			if(Q_Grmax[i ]<0.){
				std::cout<<"gr, node "<<i<<" "<<Q_Grmax[i]<<std::endl;
			}
			if(Q_Exudmax[i]<0.){
				std::cout<<"exud, node "<<i<<" "<<Q_Exudmax[i]<<std::endl;
			}
			assert((Q_Grmax[i]>=0. )&& (Q_Exudmax[i]>=0.)&& "carbon sink max < 0 ");
			assert(((Q_Rm[i]>=0. )&&(Q_Gr[i]>=0. )&& (Q_Exud[i]>=0.))&& "carbon sink < 0 ");
			if(Q_Rm[i]==0.){
				assert(((Q_Gr[i]==0.)&&(Q_Exud[i]==0.)) && "Starch_dot[i]==0. && ((Q_Sympl_dot[i]>0.)||(Q_PhlApo_dot[i]>0.))");
			}
			if(Q_Grmax[i]>0){
				assert((!(((Q_Rm[i]/Q_Rmmax_)<0.8)&&((Q_Gr[i]/Q_Grmax[i])>0.2)))&&"Rm/Rmmax < 0.8 && Gr/Grmax >0.2");
				if(Q_Exud[i]>0){assert((!(((Q_Gr[i]/Q_Grmax[i])<0.8)&&((Q_Exud[i]/Q_Exudmax[i])>0.2)))&&"Gr/Grmax < 0.8 && exud/exudmax>0.2");}
			}else{
				if(Q_Exud[i]>0){
					assert((!(((Q_Rm[i]/Q_Rmmax_)<0.8)&&((Q_Exud[i]/Q_Exudmax[i])>0.2)))&&"Rm/Rmmax < 0.8 && Gr/Grmax >0.2");
				}
			}
		
  }
}
// 


//C_ST[i]cout<<"var Q_node "<<i<<" "<<Q_ST[i]<<" "<<Q_Fl[vecI]<<" "<<Q_Exud[vecI] <<" "<< Delta_JS_ST[i]<<" "<<Q_ST_dot[i]<<endl; 
        //if(i==1)//seed node: no sink or source
		//{assert((Q_Rmmax[i]==0. &&Q_Grmax[i]==0. &&Q_Exudmax[i]==0.) && "carbon sink != 0 at seed node");
		//assert((Q_Rm[i]==0. &&Q_Gr[i]==0. &&Q_Exud[i]==0. )&& "carbon sink != 0 at seed node");
		//}else{
		//	if(Ag[i] >0){std::cout<<"ag "<<Ag[i]<<" "<<Q_Fl[i]<<" "<<(Q_Fl[i]/Ag[i])<<std::endl;}
		/*if(!(C_ST[i]>0 || Q_ST_dot[i]>=0)){
			std::cout<<"C_fluxes "<<t<<" "<<i<<" "<<CSTi<<" qdot: "<<Q_ST_dot[i]<<", ag: "<<Ag[i]<<std::flush;
			std::cout<<" rm:"<<Q_Rm[i]<<" gr:"<<Q_Gr[i]<<" ex:"<<Q_Exud[i]<<" js:"<<Delta_JS_ST[i]<<std::flush;
			std::cout<<", rmmax:"<<Q_Rmmax[i]<<" grmax:"<<Q_Grmax[i]<<" exudmax:"<<Q_Exudmax[i]<<std::flush;
			std::cout<<", ratio: "<<(Q_Rm[i]/Q_Rmmax[i])<<" "<<(Q_Gr[i]/Q_Grmax[i])<<" "<<(Q_Exud[i]/Q_Exudmax[i])<<std::flush;
			std::cout<<std::endl<<std::flush;
		//}
		assert((C_ST[i]>0 || Q_ST_dot[i]>=0) &&"C_ST< 0 && Q_ST_dot <0");
		assert((C_ST[i]>0 || Delta_JS_ST[i]>=0) &&"C_ST< 0 && Delta_JS_ST <0");*/
		
		/*
		
		//C_fluxes 0.00033001 7 qdot: 0.339211, ag: 0.339292 0
		//3.97348e-06 7.18033e-08 Delta_JS_ST[i]<<-7.73505e-05, 
		//Max: 1.22173e-07 1.91492e-05 0 
		
		
		
		Q_ST_dot[i] = Q_Fl[i] - Q_Exud[i] - Q_Gr[i] - Q_Rm[i]+ Delta_JS_ST[i];//Q_Fl_node[i] + Delta_JS_ST[i] - Q_Exud_node[i] - Q_Gr_node[i] - Q_Rm_node[i]; 
		
		
		if((Q_ST[i]<=0)&&(Q_ST_dot[i]<0)) { // control negative concentrations
            
			cout << "at t=" << t << ", seg#" << i << ": Q_ST<= 0 and Q_ST_dot < 0  =>  Q_Fu reduced by " ;
            if(Q_Exud[i] > 0.) {
                cout << "reducing exudation, " << endl ;
				// double QSTdi = Q_ST_dot[i];
				Q_Exud[i] = max(0., Q_Exud[i] + Q_ST_dot[i]);
                Q_ST_dot[i] = Q_Fl[i] - Q_Exud[i] - Q_Gr[i] - Q_Rm[i]+ Delta_JS_ST[i];
            } 
			if(((Q_ST[i]<=0.)&&(Q_ST_dot[i]<0.)) && (Q_Gr[i] > 0.)) {
                cout << "reducing growth and growth respiration, "  << endl ;
				// double QSTdi = Q_ST_dot[i];
				Q_Gr[i] = max(0., Q_Gr[i] +  Q_ST_dot[i]);
                Q_ST_dot[i] = Q_Fl[i] - Q_Exud[i] - Q_Gr[i] - Q_Rm[i]+ Delta_JS_ST[i];
            } 
			if(((Q_ST[i]<=0.)&&(Q_ST_dot[i]<0.)) && (Q_Rm[i] > 0.)) {
                cout << "reducing growth and growth respiration, ->" <<Q_Rm[i] <<" "<<  Q_ST_dot[i]<<" ";
				//cout<<(Q_Rm[i] +  Q_ST_dot[i]) <<"<-"<< endl ;
				// double QSTdi = Q_ST_dot[i];
				Q_Rm[i] =  Q_Fl[i] - Q_Exud[i] - Q_Gr[i]+ Delta_JS_ST[i];//Q_Rm[i] +  Q_ST_dot[i];
                Q_ST_dot[i] = Q_Fl[i] - Q_Exud[i] - Q_Gr[i] - Q_Rm[i]+ Delta_JS_ST[i];
				assert(Q_Rm[i]>=0);
            } 
			
			if((Q_ST[i]<=0)&&(Q_ST_dot[i]<0.))
				{
					cout<<endl<<Q_Rm[i]<<" "<<Q_Rmmax[i]<<" "<<" "<<Q_Exudmax[i]<<std::flush;
					cout<<" "<<exud_k[i]<<" "<<Length[i]<<" "<<radius_ST[i]<<" "<<Ag[i]<<std::endl<<std::flush;
					cout<<"var Q_node "<<i<<" "<<C_ST[i]<<" "<<Q_ST[i]<<std::flush;
					cout<<" "<<Q_Fl[i]<<" "<<Q_Exud[i] <<" "<< Delta_JS_ST[i]<<" "<<Q_ST_dot[i]<<std::endl<<std::flush; 
					throw std::runtime_error("(Q_ST_dot[i] < 0.) even after reducing unloading");}
        }
		
		
		
			if((Q_ST[i] <= 0.) && (Q_ST_dot[i] < 0.)) { // control negative concentrations 
			    cout << "at t=" << t << ", node#" << i << ": Q_ST <= 0 and Q_T_dot < 0  =>  Q_ST_dot set to zero by " ;
			    if(Q_PhlApo_dot[i] > 0.) {
                cout << "reducing exudation, " << endl ;
				// double QSTdi = Q_ST_dot[i];
				Q_PhlApo_dot[i] = max(0., Q_PhlApo_dot[i] + Q_ST_dot[i]);
				Fu_lim = Q_PhlApo_dot[i] + Starch_dot[i] + Q_ParApo_dot[i];
                Q_ST_dot[i] = Q_Fl[i] - Q_Exud[i] - Q_Gr[i] - Q_Rm[i]+ Delta_JS_ST[i];
            } 
			if(((Q_ST[i]<=0.)&&(Q_ST_dot[i]<0.)) && (Q_Gr[i] > 0.)) {
                cout << "reducing growth and growth respiration, "  << endl ;
				// double QSTdi = Q_ST_dot[i];
				Q_Gr[i] = max(0., Q_Gr[i] +  Q_ST_dot[i]);
                Q_ST_dot[i] = Q_Fl[i] - Q_Exud[i] - Q_Gr[i] - Q_Rm[i]+ Delta_JS_ST[i];
            } 
			if(((Q_ST[i]<=0.)&&(Q_ST_dot[i]<0.)) && (Q_Rm[i] > 0.)) {
                cout << "reducing growth and growth respiration, ->" <<Q_Rm[i] <<" "<<  Q_ST_dot[i]<<" ";
				//cout<<(Q_Rm[i] +  Q_ST_dot[i]) <<"<-"<< endl ;
				// double QSTdi = Q_ST_dot[i];
				Q_Rm[i] =  Q_Fl[i] - Q_Exud[i] - Q_Gr[i]+ Delta_JS_ST[i];//Q_Rm[i] +  Q_ST_dot[i];
                Q_ST_dot[i] = Q_Fl[i] - Q_Exud[i] - Q_Gr[i] - Q_Rm[i]+ Delta_JS_ST[i];
				assert(Q_Rm[i]>=0);
            } 
		    }
		*/
		//todo