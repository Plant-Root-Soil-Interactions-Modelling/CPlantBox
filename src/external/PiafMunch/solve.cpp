/*
* PiafMunch (v.2) -- Implementing and Solving the Munch model of phloem sap flow in higher plants
*
* Copyright (C) 2004-2019 INRA
*
* Author: A. Lacointe, UMR PIAF, Clermont-Ferrand, France
*
* File: solve.cpp
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

#include <math.h>
#include <vector>
#include "PM_arrays.h"
#include "runPM.h"

/**********************************************  Network Architecture: *********************************************/

//extern Fortran_vector Q_FuFl_seg;
//extern Fortran_vector Q_FuFl_node;
Fortran_vector Q_ST_seg_init;

extern Fortran_vector  Q_Exud_seg, Q_Gr_seg, Q_Rm_seg, Q_Fl_seg;
extern Fortran_vector C_ST_seg, Q_ST_seg, Q_ST_seg_new, Q_ST_newFuFl, grad_Q_ST_seg;
extern Fortran_vector vol_ST_seg, Length;
vector<int> Deg ;
int Nt, Nc ; // Nt = Total number of nodes, Nc = total number of connections in network.
int N1R ; // Number of  'root tip's, i.e. nodes of conn.order 1 that have imposed soil water potential as a limit condition.
extern int N1L ; // Organ type # 1L : 'leaf end' ; has imposed transpiration stream as a limit condition
extern Index_vector i_[] ; // i_[o = 0..8][k = 1..No] = id# of kth node in (1-based list)  list of all nodes of conn.ord. co(o). For o=0: No=N1L ; o=1: No=N1R ; o=2..8: No=N[o] :
extern Index_vector &RootEnds, &LeafEnds ; // (= i_[1] and i_[0], resp.) : label which nodes are network ends = nodes of co=1 , which are either 'leaf' or 'root' ends :
vector<int> jf_RootEnds ; // indices to label end nodes = those of conn.order 1
extern vector<int> I_Upflow, I_Downflow ; // I_Upflow(resp.I_Downflow)[jf=1..Nc] = id# du noeud amont (resp. aval) : jf = JF(i,i2) > 0 si I_Upflow[(abs(jf)]==i, i.e. si I_Downflow[(abs(jf)]==i2
SpUnit_matrix Delta2 ; // describe hydraulic architecture (topology)
SpUnit_matrix Delta2abs ; // describe hydraulic architecture (topology)
Sparse_matrix Deltaabs ; // describe hydraulic architecture (topology)
Sparse_matrix Delta2W;
Sparse_matrix X ;	// intermediaires de calcul ;
int** ipiv_ptr ; void**TM_ptr ; // id. -- initialise dans inialize_hydric()
#ifdef Full_Matrix // solve full linear system for all hydraulic variables, chained as YY below :
Fortran_vector YY, SM ; // resp., inconnue et second membre de l'eq. matricielle
#else // reduce linear system to only one unknown hydraulic variable (P_ST) => significant speed up and reduce required memory -- needs pre-solving recalculation if local model is changed
SpUnit_matrix Delta ; // =  - Transpose(Delta2) : describe hydraulic architecture (topology)
Sparse_matrix M1, Delta_rxyl, Delta_rphl ;	// intermediaires de calcul ;
Fortran_vector Km, rG, O ; //  intermediaires de calcul (Km : rien a voir avec Michaelis !)
Fortran_vector inv_Km, rs_rG, inv_rPhlM, inv_rG ; //  intermediaires de calcul, derives ou inverses des precedents
# endif
double C, _0089_099803 = 0.089/0.99803 ; // intermediaire de celcul de la molalite
extern int jf, i ; int is ; // is = # of current integration time segment (first = 1)

/******************************************  Environmental Variables: *********************************************/
extern double T ; // (K) absolute temperature : set in GUI, but may be updated anytime in function 'parameter_and_boundary_conditions()' at the end of file 'PiafMunch2.cpp'.
// The following 2 environment data are set by the user in function  'parameter_and_boundary_conditions()'  at the end of file 'PiafMunch2.cpp' :
extern Fortran_vector Transpirat		; // Leaf transpiration rate (used in water-fluxes calc.)				(mmol / h)
extern Fortran_vector PsiSoil			; // Soil water potential at root end (used in water-fluxes calc.)		(MPa)

/******************************************  Constants and Parameters: *********************************************/
#define R 83.14//0.0083143	// constante des gaz parfaits -											(MPa ml K-1 mmol-1)
double TdC, dEauPure,  siPhi, newPhi ; // pour visc. calc. par  www.seas.upenn.edu et/ou NonLinPsi
bool Adv_BioPhysics ; // true if  non-zero sugar specific volume, osmotic pot.=non-linear function of molality (Thompson and Holbrook), and viscosity changes with C_TC (Thompson and Holbrook ; Seas, Flanagan)  ; set in IntroDialogBox

// Hydro Resistance Parameters (set in GUI) :
Fortran_vector r_Xyl			; //(MPa h / ml) : xylem water resistance
Fortran_vector r_ST			; // (MPa h / ml) : axial phloem water resistance
Fortran_vector r_ST_ref 	; // r_ST_ ref. values for C_TC =  0.5 mmol / ml sap solution
Fortran_vector r_abs			; // soil - root  water resistances (may be changed by user)			(MPa h / ml)
Fortran_vector r_Trsv			; // Transverse xylem to phloem apoplasm water resistance										(MPa h / ml)
Fortran_vector r_PhlMb		; // Transmembrane Phloem (ST,CC) Apoplasm to Symplasm water resistance		(MPa h / ml)
Fortran_vector r_ParMb		; // Transmembrane Parenchyma Apoplasm to Symplasm water resistance		(MPa h / ml)
Fortran_vector r_Apo			;//  Lateral parenchyma to phloem Apoplastic (MPa h / ml)
Fortran_vector r_Sympl			; // Lateral parenchyma to phloem ST Symplasmic water resistance (MPa h / ml)

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

int resistance_changed = NONE ; // = NONE a priori
void UpdateResistances(double t) ; // in  initialize.cpp
void Smooth_Parameter_and_BoundaryConditions_Changes(int s, double t) ; // User-editable (implemented in 'PiafMunch2.cpp') ; s = # of integration segment
//void C_fluxes(double t, int Nt) ; // in  PiafMunch2.cpp
void vector_init(double t, double *y, double *y_dot);

// ******************  mere C-fluxes related variables or parameters ********** :
extern double *Q_ST, *Q_Mesophyll, *Q_RespMaint, *Q_Exudation, *Q_Growthtot, *Q_S_ST, *Q_Mucil ;		  // components of vector y as used in diff. system f()...
extern double *Q_ST_dot, *Q_Mesophyll_dot, *Q_Rm_dot, *Q_Exud_dot, *Q_Gtot_dot, *Q_S_ST_dot, *Q_Mucil_dot ; //... and its derivatives.  ;
extern double *vol_Sympl ;
extern Fortran_vector JS_ST, C_amont, JS_Sympl, JS_Apo, RespMaint ;
extern Fortran_vector vol_ST, vol_PhlApo, vol_ParApo ;
extern Fortran_vector r_abs, PsiSoil, Transpirat, Input, Q_RespMaintSyn ;
extern Fortran_vector P_ST_dot, P_Sympl_dot ; /******* variation rate of any variable X is noted: X_dot = dX/dt *******/
//Index_vector i_amont					; // i_amont[i] =Id# of  *true*  upflow node
Fortran_vector i_amont					;
extern Fortran_vector C_SymplUpflow, C_ApoUpflow ; // to derive JS_Sympl, and maybe JS_Apo, from resp. water fluxes JW_xxx

// ******************  auxiliary sugar- or anatomy- related, non-hydraulic variables or parameters that may be involved in water-system equations ********** :
extern Fortran_vector C_ST, C_PhlApo		; // Concentration of sugar in phloem  sieve tubes  and  apoplasm, resp.							(mmol / ml)   (ml of solution)
extern Fortran_vector  JS_PhlMb, JS_ParMb ; //  (if Adv_BioPhysics)  to compute NZS and NZSP (see next line)
// NZS : optional non-zero volume sugar flow (not a distinct variable)   (ml / h) : NZS = JS_PhlMb * PartMolalVol (=0.2155 in Thompson and Holbrook -- 0.214 might be more accurate)
extern Fortran_vector C_Sympl, C_ParApo		; // Concentration of sugar in lateral tissue symplasm  and  apoplasm, resp.								(mmol / ml)   (ml of solution)
extern Fortran_vector Delta_JS_ST ; // sera la composante purement phloemienne de Q_TC_dot[ ]							(mmol / h)

// tracer-specific add-ins :
extern double *Q_RespMaintmax, *TracerQ_Mesophyll, *TracerQ_RespMaint, *Q_S_Mesophyll, *Q_Growthtotmax ;		  // components of vector y as used in diff. system f()...
extern double *Q_Rmmax_dot, *TracerQ_Mesophyll_dot, *TracerQ_Rm_dot, *Q_S_Mesophyll_dot, *Q_Gtotmax_dot ; //... and its derivatives.  ;
extern Fortran_vector TracerJS_ST, TracerC_Sympl, TracerC_ST, TracerJS_Sympl, TracerJS_Apo, TracerJS_ParMb, TracerJS_PhlMb, TracerC_PhlApo, TracerC_ParApo ;
extern Fortran_vector TracerQ_RespMaintSyn, TracerInput, TracerRespMaint, TracerC_SymplUpflow, TracerC_ApoUpflow ;
extern Fortran_vector TracerRatioSympl, TracerRatioQ_RespMaint ;
extern Fortran_vector Delta_TracerJS_ST ; // sera la composante purement phloemienne de Q_TC_dot[ ]							(mmol / h)

/*************************** VARIABLES INVOLVED IN HYDRIC SYSTEM (Water fluxes): ******************************* */
extern double* vol_Sympl_dot ; // Rate of symplasmic vol. change (eq. 8)
Fortran_vector P_Xyl		    	; // Xylem Pressure (= water potential since there is no solute in this version)					(MPa)
Fortran_vector Psi_Xyl		    	; // Xylem Pressure (= water potential since there is no solute in this version)					(MPa)
Fortran_vector Psi_ST	 		; // Phloem sieve-tube (+CC) water potential												(MPa)
Fortran_vector Psi_PhlApo		; // Phloem Apoplam water potential												(MPa)
Fortran_vector Psi_ParApo		; // Lateral parenchyma apoplasmic water potential												(MPa)
Fortran_vector P_PhlApo		; // Phloem Apoplam pressure												(MPa)
Fortran_vector P_ParApo		; // Lateral parenchyma apoplasmic pressure												(MPa)
Fortran_vector Psi_Sympl		; // Lateral parenchyma symplasmic water potential												(MPa)
Fortran_vector P_ST			; // Phloem sieve-tube (+CC) turgor pressure 								(MPa)
Fortran_vector P_Sympl		; // Lateral parenchyma symplasmic turgor pressure 								(MPa)
Fortran_vector Absorb			; // Root end water absorbtion flux										(ml / h)
Fortran_vector JW_Xyl			; // Axial xylem water flux
Fortran_vector JW_ST			; // Axial phloem sieve-tube liquid flux
Fortran_vector JW_Trsv			; // Transverse xylem to phloem apoplasm water flux									(ml / h)
Fortran_vector JW_PhlMb		; // Transmembrane Phloem (ST,CC) Apoplasm to Symplasm water flux
Fortran_vector JW_ParMb		; // Transmembrane Parenchyma Apoplasm to Symplasm water flux
Fortran_vector JW_Apo			; // Lateral parenchyma to phloem Apoplastic water flux									(ml / h)
Fortran_vector JW_Sympl		; // Lateral parenchyma to phloem ST Symplasmic liquid flux									(ml / h)

/********************* WATER SYSTEM EQUATIONS  ***************************/


double  PartMolalVol =0;
//comes from Genotelle, J. Expression de la viscosite des solutions sucrees. Ind. Aliment. Agric. 1978, 95, 747-755
//prooved to hold by https://doi.org/10.1021/ie000266e
//taken here as presented in "Sucrose Properties and Applications" for pure sucrose solution
//see 10.1007/978-1-4615-2676-6_6
//Mathlouthi, M.; Reiser, P., 1995
//eq. 6.29
void PhloemFlux::update_viscosity() { // called if (Adv_BioPhysics)
	static double T_old(-9999.);//V_ref(-9999.) ; // memory of previous T (do not compute dEauPure, siPhi, etc. if T unchanged), and V_ref
	double d, siEnne,  mu ; 
	if (T_old != TairK_phloem) {
		TdC = TairK_phloem - 273.15;
		//in g/L or mg/cm3
		dEauPure = (999.83952 + TdC * (16.952577 + TdC * (- 0.0079905127 + TdC * (- 0.000046241757 + TdC * (0.00000010584601 + TdC * (- 0.00000000028103006)))))) / (1 + 0.016887236 * TdC); 
		siPhi = (30 - TdC) / (91 + TdC) ;  T_old = TairK_phloem ;//  R.Gilli 1997, after Mathlouthi and Genotelle 1995 - valid for any T :
	}
   for (int i=1 ; i <= Nc ; i++) { // this loop computes the viscosity profile, temporarily stored as vector r_ST_ before scaling to actual r_ST values:
		C = C_amont[i] ; // (mmol / ml solution)
		//  R.Gilli 1997, after Mathlouthi and Genotelle 1995 - valid for any T :
		if (C < 0.) C = 0. ; // fix any artefact from solver (may try C<0 even if actual C never does)
		//342.3 g/mol or mg/mmol
		double PartMolalVol_ = 0;//0.2155;
		d = C * 342.3 + (1 - C * PartMolalVol_) * dEauPure ;//in mg/cm3
		siEnne = (100 * 342.30 * C) / d ; // actually this is sc = sucrose content (g.suc. % g.solution) ; 342.30 = molar mass of sacch.
		siEnne /= 1900 - (18 * siEnne) ;
		//mPa s
		mu =  pow(10, ((22.46 * siEnne) - 0.114 + (siPhi * (1.1 + 43.1 * pow(siEnne, 1.25) )))) ; // peut atteindre des valeurs > 1.e200 !! (sans signification evidemment -- le sucre doit precipiter bien avant !!)
		mu = mu /(24*60*60)/100/1000; //mPa s to hPa d, 1.11837e-10 hPa d for pure water at 293.15K
		r_ST[i] = mu*r_ST_ref[i];
   }
}




double R_ = 83.14;//cm^3 hPa K-1 mmol-1
void PhloemFlux::f(double t, double *y, double *y_dot) { // the function to be processed by the solver
	
	double RT = R_*TairK_phloem ; // int k ; T may be changed anytime, in function 'parameter_and_boundary_conditions(t)  in PiafMunch2.cpp
	static Fortran_vector dummy(Nt);
	Q_ST = y ;							// note: zero_indices  Q_ST[0],..., Q_RespMaint[0]... are ignored
	Q_Mesophyll = Q_ST + Nt ; 
	Q_RespMaint = Q_Mesophyll + Nt ; 
	Q_Exudation = Q_RespMaint + Nt ; 
	Q_Growthtot = Q_Exudation + Nt ; 
	
	Q_RespMaintmax = Q_Growthtot + Nt ; 
	Q_Growthtotmax = Q_RespMaintmax + Nt ; 
	Q_S_Mesophyll = Q_Growthtotmax + Nt ; 
	
	//if delete, lower neq
	Q_S_ST = Q_S_Mesophyll + Nt;
	Q_Mucil = Q_S_ST + Nt ;
	
	for (int i=1; i<=Nt; i++)  {
		double volSTi = vol_ST[i];
		double QSTi = Q_ST[i];
		
		if (QSTi < 0.){ QSTi = 0. ; Q_ST[i] =0;}// fix any artefact from solver (may try C<0 even if actual C never does)
		if (Q_Mesophyll[i] < 0.){ Q_Mesophyll[i] =0;}// fix any artefact from solver (may try C<0 even if actual C never does)
		if (Q_S_ST[i] < 0.){ Q_S_ST[i] =0;}// fix any artefact from solver (may try C<0 even if actual C never does)
		if (Q_S_Mesophyll[i] < 0.){ Q_S_Mesophyll[i] =0;}// fix any artefact from solver (may try C<0 even if actual C never does)
		C_ST[i] = QSTi / volSTi ; // Concentration of sugar in sieve tubes		(mmol / ml)
	
	}
	Q_ST_dot = y_dot ; 
	Q_Mesophyll_dot = Q_ST_dot + Nt ; 
	Q_Rm_dot = Q_Mesophyll_dot + Nt ; 
	Q_Exud_dot = Q_Rm_dot + Nt ;
	Q_Gtot_dot = Q_Exud_dot + Nt ;
	
	Q_Rmmax_dot = Q_Gtot_dot + Nt ; 
	Q_Gtotmax_dot = Q_Rmmax_dot + Nt ; 
	Q_S_Mesophyll_dot = Q_Gtotmax_dot + Nt ; 
	
	
	Q_S_ST_dot = Q_S_Mesophyll_dot + Nt ;
	Q_Mucil_dot = Q_S_ST_dot + Nt ;
	
	//Add later
	/*if (Adv_BioPhysics) {
		dummy.set(1.) ; dummy.sub_mult(PartMolalVol, C_ST) ; // dummy = 1 - PartMolalVol * C_ST
		Psi_ST.set_elediv(C_ST, dummy) ; // temporary buffer for  molarity  m_  =  C_ST / (1 - PartMolalVol * C_ST) ; true Psi_ST value will be set later
		Psi_ST *= _0089_099803 * Psi_ST + 0.998 ; // temp. buffer for  m_ * (0.998 + 0.089/0.99803 * m_) =  M_ST ; true Psi_ST value will be set later
		Psi_ST *= - RT ;
	}	*/
	
	P_ST.set(C_ST);
	P_ST *= RT;
	if(usePsiXyl){P_ST += Psi_Xyl ;}
	JW_ST.set_matmult(Delta, P_ST) ; 
	for(int j = 1 ; j <= Nc ; j ++) {
		//int i_aval; double C_aval;
		if(JW_ST[j] > 0)
		{  
			i_amont[j] = I_Upflow[j] ; 
		}else{  
			i_amont[j] = I_Downflow[j] ; 
		}
		C_amont[j] = C_ST[i_amont[j]] ; 
		
		if((errorID == I_Upflow[j])||(errorID == I_Downflow[j]))//found an error in last run of C_fluxes function
		{
			std::ofstream outfile;
			outfile.open("errors.txt", std::ios_base::app); // append instead of overwrite
			outfile<< "JW_ST "<<JW_ST[j]<<" r_ST "<<r_ST[j]<<" ide "<<errorID<<" ids "<<I_Upflow[j]<<" "<< I_Downflow[j]<<std::flush;
			outfile<<" cst: "<<C_ST[I_Upflow[j]]<<" "<<C_ST[I_Downflow[j]] <<std::flush;
			outfile<<" pst: "<<P_ST[I_Upflow[j]]<<" "<<P_ST[I_Downflow[j]] <<std::flush;
			outfile<<" vol: "<<vol_ST[I_Upflow[j]]<<" "<<vol_ST[I_Downflow[j]] <<std::flush;
			outfile<<std::endl<<std::flush;
			assert(false&&"PhloemFlux::C_fluxes : negative maximal flux of sucrose concentration");
			
		}
		if (C_amont[j] < 0.) C_amont[j] = 0. ; // fix any artefact from solver (may try C<0 even if actual C never does)
	}
	if(update_viscosity_){update_viscosity() ;}
	JW_ST /= r_ST;
	
	JS_ST.set_elemult(JW_ST, C_amont) ;	// 	i.e.   JS_ST = JW_ST * C_amont			   (eq. 11)
	Delta_JS_ST.set_matmult(Delta2, JS_ST) ;
	
	C_fluxes(t, Nt) ; //see PiafMunch2.cpp
	
}
