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
double* Q_Mucil = NULL;
double* Q_S_ST = NULL;
double* Q_Mesophyll = NULL						; // Amount of sugar in parenchyma symplasm										(mmol)
double* Q_ST = NULL						; // Amount of sugar in sieve tubes	= C_TC * Vol_ST						(mmol)
double* Q_RespMaint = NULL						; // Amount of starch in parenchyma										(mmol)
double* Q_Exudation = NULL                 ; // amount of sugar in phloem apoplasm   (mmol)
double* Q_Growthtot = NULL                 ; // amount of sugar in lateral parenchyma apoplasm   (mmol)
Fortran_vector JS_ST						; // (mmol / h)  : Axial phloem sugar flux
Fortran_vector JS_PhlMb				; // Phloem cross-membrane sugar fluxes into sieve tubes from apoplasm				(mmol / h)
Fortran_vector RespMaint					; // Maintenance respiration rate										(mmol / h)
Fortran_vector Q_RespMaintSyn						; // Michaelis-Menten rate of starch synthesis from sugar substrate						(mmol sug.eq./ h)
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
Fortran_vector Delta_JS_ST ; // sera la composante purement phloemienne de Q_TC_dot[ ]							(mmol / h)
extern Fortran_vector P_Sympl		; // Lateral parenchyma symplasmic turgor pressure 								(MPa)
double Q_Rm_dot_alt ; // for an alternate, target-oriented,  expression of starch variation rate

/******** TRACER-RELATED VARIABLES ***********************************************************/
double* TracerQ_Mesophyll = NULL						; // Amount of tracer in parenchyma										(MBq)
double* Q_RespMaintmax = NULL						; // Amount of soluble tracer in sieve tubes	= TracerC_ST * Vol_ST						(MBq)
double* TracerQ_RespMaint = NULL						; // Amount of insoluble tracer in parenchyma										(MBq)
double *Q_S_Mesophyll=NULL, *Q_Growthtotmax=NULL ; // Amount of soluble tracer in apoplasm	= TracerC_xxxApo * Vol_xxxApo			(MBq)
Fortran_vector TracerJS_ST						; // TracerJS_ST[i] (MBq / h)  = TracerJS_ST_[i-1], i = 2..N  : Axial phloem soluble tracer flux ; TracerJS_ST[1]= NA, whereas TracerJS_ST_[1] = TracerJS_ST[2]
Fortran_vector TracerJS_PhlMb						; // Lateral (Apoplasmic) soluble tracer fluxes into sieve tubes from parenchyma				(MBq / h)
Fortran_vector TracerRespMaint					; // Tracer Maintenance respiration rate										(MBq / h)
Fortran_vector TracerQ_RespMaintSyn						; // Michaelis-Menten rate of tracer starch synthesis from tracer sugar substrate						(MBq / h)
Fortran_vector TracerInput					; // Leaf photosynthetic Tracer Assimilation rate (boundary condition)			(MBq / h)
Fortran_vector TracerC_Sympl						; // Concentration of soluble tracer in parenchyma symplasm= TracerQ_Mesophyll / vol_Sympl			(MBq / ml)
Fortran_vector TracerC_ST							; // Concentration of soluble tracer in sieve tubes								(MBq / ml solution))
Fortran_vector TracerC_PhlApo, TracerC_ParApo ; // Concentration of soluble tracer in apoplasms					(MBq / ml solution))
Fortran_vector TracerRatioSympl ; //   TracerQ_Mesophyll / Q_Mesophyll = TracerC_Sympl / C_Sympl  (MBq / mmol)
Fortran_vector TracerRatioQ_RespMaint ; //   = TracerQ_RespMaint / Q_RespMaint (MBq / mmol)
Fortran_vector Delta_TracerJS_ST ; // sera la composante purement phloemienne de Q_Rmmax_dot[ ]
Fortran_vector TracerJS_Sympl, TracerJS_Apo, TracerJS_ParMb ; // Lateral (Sympl., Apopl.; cross-membr...) soluble tracer fluxes FROM sieve tubes INTO parenchyma							(MBq / h)     !!! ATTENTION : sens positif oppose a JS_Trsv !!!
Fortran_vector TracerC_SymplUpflow					;  // upflow tracer concentration (mmol / ml) for TracerJS_Sympl
Fortran_vector TracerC_ApoUpflow ;				; // upflow tracer concentration (mmol / ml) for TracerJS_Apo

/******* variation rate of any variable X is noted: X_dot = dX/dt : *******/
double * Q_Mesophyll_dot = NULL, * Q_ST_dot = NULL, * Q_Rm_dot = NULL, *Q_Exud_dot = NULL, *Q_Gtot_dot = NULL,  *Q_S_ST_dot = NULL , *Q_Mucil_dot = NULL ;
double * TracerQ_Mesophyll_dot=NULL, * Q_Rmmax_dot=NULL, * TracerQ_Rm_dot=NULL, *Q_S_Mesophyll_dot=NULL, *Q_Gtotmax_dot=NULL ;
// Next 2 variables are not considered as such, but as possible inputs to compute vol_Sympl_dot :
Fortran_vector P_ST_dot, P_Sympl_dot			; //  dP_ST/dt , dP_Sympl/dt					(MPa h / h)    -- for elasticity...

/********************* C-FLUXES-RELATED BIOPHYSICAL and PHYSIOLOGICAL PARAMETERS **********************************/
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
Fortran_vector vol_ST, vol_ParApo, vol_PhlApo,vol_Seg ; // (ml) total volume of, resp. : sieve tubes, and  parenchyma and phloem apoplasm
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


void PhloemFlux::C_fluxes(double t, int Nt)  
{
	TairC = TairK_phloem - 273.15;
	for (int i = 1 ; i <= Nt ; i++) 
	{ // edit (make different loops) to enter specific equations for specific nodes or conn.orders
		int cpp_id = i -1;// o go from Fortran_vector numeration to cpp vector numeration
		double CSTi = max(0.,C_ST[i]);// From A.Lacointe: solver may try C<0 even if actual C never does
		//double CSTi_exud = max(0.,C_ST[i]);// From A.Lacointe: solver may try C<0 even if actual C never does
		double Cmeso = max(0.,Q_Mesophyll[i]/vol_ParApo[i]);//concentration in meosphyll compartment
		//Q_Fl[i] = k_meso*max(Cmeso - CSTi, 0.);//flux from mesophyll to sieve tube
		 
		double Q_Rmmax_ ;double Q_Exudmax_;double Fu_lim;
        
		//starch		
			//ST	
		//Q_S_ST_dot[i] = 0;
		double denominator = kM_S_ST + CSTi;
		double StarchSyn = 0;
		if(denominator != 0)
		{
			StarchSyn = Vmax_S_ST * max(0.,Q_ST[i]) / (denominator) ; //  (Vmax and kHyd below) or k3 (below),  or all three, should be zero
		}
		Q_Mucil_dot[i] = std::max(0.,k_mucil_[cpp_id] *   Q_S_ST[i]);//mucilage exudation
		//an alternate, target-oriented,  expression of starch variation rate, mutually exclusive of (AmSyn - kHyd * Amid), so that...
        double Starch_dot_alt = k_S_ST * (CSTi - C_targ) * vol_ST[i] ;	
		Q_S_ST_dot[i] = StarchSyn + Starch_dot_alt - kHyd_S_ST *std::max(0., Q_S_ST[i]) ; // (Vmax and kHyd) or k3 (below),  or all three, should be zero
		
        if((Q_S_ST[i] <= 0.) && (Q_S_ST_dot[i] < 0.)) { // control negative starch concentrations (remove 4-lines-block if not relevant)
            //cout << "at t=" << t << ", node#" << i << ": Starch <= 0 and Starch_dot < 0  =>  Starch_dot set to zero" << endl ;
            Q_S_ST_dot[i] = 0. ;
			Q_Mucil_dot[i] = 0. ;
        }
        double st_2_starch = Q_S_ST_dot[i] ;//+ Q_Mucil_dot[i];
		Q_S_ST_dot[i] -=  Q_Mucil_dot[i];
		
        	//MEsophyll
		Q_S_Mesophyll_dot[i] = 0.;
        if(vol_ParApo[i]>0) {
            denominator = kM_S_Mesophyll + Cmeso;
            StarchSyn = 0;
            if(denominator != 0)
            {
                StarchSyn = Vmax_S_Mesophyll * max(0.,Q_Mesophyll[i]) /  (denominator) ; //  (Vmax and kHyd below) or k3 (below),  or all three, should be zero
            }
            //an alternate, target-oriented,  expression of starch variation rate, mutually exclusive of (AmSyn - kHyd * Amid), so that...
            Starch_dot_alt = k_S_Mesophyll * (Cmeso - C_targMesophyll) * vol_ParApo[i] ;	
            Q_S_Mesophyll_dot[i] = StarchSyn - kHyd_S_Mesophyll * Q_S_Mesophyll[i] + Starch_dot_alt ; // (Vmax and kHyd) or k3 (below),  or all three, should be zero
            if((Q_S_Mesophyll[i] <= 0.) && (Q_S_Mesophyll_dot[i] < 0.)) { // control negative starch concentrations (remove 4-lines-block if not relevant)
                //cout << "at t=" << t << ", node#" << i << ": Starch <= 0 and Starch_dot < 0  =>  Starch_dot set to zero" << endl ;
                Q_S_Mesophyll_dot[i] = 0. ;
            }
        }
		
		Q_Fl[i] = (Vmaxloading *len_leaf[i])* Cmeso/(Mloading + Cmeso) * exp(-CSTi* beta_loading);//phloem loading. from Stanfield&Bartlett_2022
		CSTi = max(0., CSTi-CSTimin); //if CSTi < CSTimin, no sucrose usage
		CSTi_exud.at(cpp_id) = max(0., max(0.,CSTi)-CSTimin_exud); //if CSTi < CSTimin, no sucrose usage
		Crsi_exud.at(cpp_id) = max(0.,Csoil_node[cpp_id]-CSTimin_exud); //if CSTi < CSTimin, no sucrose usage
		
		CSTi_delta.at(cpp_id) = max(0.,CSTi_exud.at(cpp_id)-Crsi_exud.at(cpp_id)); //concentration gradient for passive exudation. TODO: take Csoil from dumux 
		Q_Rmmax_ = (Q_Rmmax[i] + krm2[i] * CSTi) * pow(Q10,(TairC - TrefQ10)/10);//max maintenance respiration rate
		
		Q_Exudmax_ = CSTi_delta.at(cpp_id)*Q_Exudmax[i];//max exudation rate
		Fu_lim = (Q_Rmmax_  + Q_Grmax[i])* (CSTi/(CSTi + KMfu));//active transport of sucrose out of sieve tube			
		Q_ST_dot[i] = Q_Fl[i] - Fu_lim -Q_Exudmax_ + Delta_JS_ST[i] - st_2_starch;//variation of sucrose content in node
		
		//Q_meso_dot:
		Q_Mesophyll_dot[i] = Ag[i] -Q_Fl[i] - Q_S_Mesophyll_dot[i] ;//variaiton of sucrose content in mesophyll compartment 
		
		Input[i] = Q_Fl[i];//phloem loading
		
		//Q_Rm_dot:
		Q_Rm_dot[i] = min(Fu_lim, Q_Rmmax_);//realized rate of maintenance respiration 
		
		//Growth:
		//add max(X,0.) in case of issues with rounding
		Q_Gtot_dot[i] = max(min(Fu_lim - Q_Rm_dot[i], Q_Grmax[i]),0.);//realized rate of sucrose usage for growth + growth respiration
		//Exudation:
		Q_Exud_dot[i] =  Q_Exudmax_;//realized rate of exudation
		
		////save maximum sucrose usage rate for post-processing
		//Resp_maintmax:
		Q_Rmmax_dot[i] = Q_Rmmax_ ;
		//Growthmax:
		Q_Gtotmax_dot[i] = Q_Grmax[i];
		
		if(doTroubleshooting){
			std::cout<<"C_fluxes "<<i<<" "<<vol_ST[i]<<" "<<vol_ParApo[i]<<" "<<vol_Seg[i]<<" CSTimin "<<CSTimin<<std::endl;
			std::cout<<"max(0.,C_ST[i]) "<<max(0.,C_ST[i])<<std::endl;
			std::cout<<" C_ST[i] "<<C_ST[i]<<" Q_ST[i] "<<Q_ST[i]<<" "<<Q_Fl[i]<<" "<<CSTi<<" "<<Cmeso<<" "<<len_leaf[i]<<" max(0., CSTi-CSTimin) "<< max(0., CSTi-CSTimin)<<std::endl;

			std::cout<<Q_Rmmax_<<" "<<Q_Rmmax[i]<<" "<< krm2[i]<<" "<<CSTi_delta.at(cpp_id)<<std::endl;

			std::cout<<Q_Exudmax_<<" Fu_lim "<<Fu_lim<<" Q_ST_dot "<<Q_ST_dot[i]<<" "<<Q_Mesophyll_dot[i]<<" "<<Input[i]<<" "<<Q_Rm_dot[i]<<std::endl;
			std::cout<<"Qgri "<<Q_Gtot_dot[i] <<" Q_Exudmax_ "<<Q_Exud_dot[i]<<" Q_Rmmax_ "<<Q_Rmmax_dot[i]<<" Qgrmaxi "<<Q_Gtotmax_dot[i]<<std::endl;
			std::cout<<"Qmeso "<<Q_Mesophyll[i]<<" "<<Ag[i]<<std::endl;
		}
		
		//check if error
		if(((Q_ST[i]<= 0) &&(Q_ST_dot[i] < 0))||((Q_Rm_dot[i]<0)||(Q_Gtot_dot[i]<0)||(Q_Exud_dot[i]<0))){
			std::cout<<"error, see file errors.txt"<<std::endl;
			std::ofstream outfile;
			outfile.open("errors.txt", std::ios_base::app); // append instead of overwrite
			outfile<< std::endl<<"C_fluxes "<<t<<" "<<i<<" "<<CSTi<<" "<<CSTimin<<" "<<Q_ST[i]<<" qdot: "<<Q_ST_dot[i]<<", Fu: "<<Fu_lim <<std::flush;
			outfile<<" rm:"<<Q_Rm_dot[i]<<" maxrm: "<<Q_Rmmax_dot[i]<<std::flush;
			outfile<<" gr:"<<Q_Gtot_dot[i]<<" grmax: "<<Q_Gtotmax_dot[i]<<" exud:"<<Q_Exud_dot[i]<<" js:"<<Delta_JS_ST[i]<<std::flush;
			outfile<<" vol_ST:"<<vol_ST[i]<<std::flush;
			outfile<<std::endl<<std::flush;
			errorID = i; //will stop computation in PhloemFlux::f
		}
             
	}
		
}


