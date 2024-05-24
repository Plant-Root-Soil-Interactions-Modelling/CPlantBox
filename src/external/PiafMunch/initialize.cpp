/*
* PiafMunch (v.2) -- Implementing and Solving the Munch model of phloem sap flow in higher plants
*
* Copyright (C) 2004-2019 INRA
*
* Author: A. Lacointe, UMR PIAF, Clermont-Ferrand, France
*
* File: initialize.cpp
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
extern Fortran_vector vol_ST_seg;//, Length;
extern Fortran_vector vol_ST, vol_ParApo;
extern vector<int> Deg ;
extern Sparse_matrix Delta2W;
extern SpUnit_matrix Delta2		; // describe hydraulic architecture (topology)
extern Sparse_matrix Deltaabs		; // describe hydraulic architecture (topology)
#ifdef Full_Matrix // solve full linear system for all hydraulic variables, chained as YY below :
extern Fortran_vector YY, SM ; // resp., inconnu et second membre de l'eq. matricielle
#else // reduce linear system to only one unknown hydraulic variable (P_ST) => significant speed up and reduce required memory -- needs pre-solving recalculation if local model is changed
extern SpUnit_matrix Delta		; // describe hydraulic architecture (topology)
extern Sparse_matrix M1, Delta_rxyl, Delta_rphl ; // intermediary complex matrices to speed up calculations
Sparse_matrix A22, A22_I, Delta2rpDelta, Delta2rxDelta, D6, I, X1 ; // intermediary complex matrices to speed up calculations
Fortran_vector dummy, dummy_ ; // temporary buffers to speed up calculations(not have to create new objects at each call)
extern Fortran_vector Km, rG, O ; // (def. in hydric_solve.cpp) intermediaires de calcul des matrices, cf. fichier Excel (Km : rien a voir avec Michaelis !)
extern Fortran_vector inv_Km, rs_rG, inv_rPhlM, inv_rG ; //  intermediaires de calcul, derives ou inverses des precedents
#endif
extern Sparse_matrix X ;
extern int** ipiv_ptr; extern void**TM_ptr; // buffers to speed up linear solving
extern int jf, i, j ;

/**********************************************  Network Architecture: *********************************************/
extern int Nt, Nc ; // Nt = Total number of nodes, Nc = total number of de connections in network.
extern int N1L ; // Number of 'leaf end's, i.e. nodes of conn.order 1 that have imposed transpiration stream as a limit condition.
extern int N1R ; // Number of  'root tip's, i.e. nodes of conn.order 1 that have imposed soil water potential as a limit condition (N[1] = N1L + N1R).
extern Index_vector i_[] ; // i_[o = 0..8][k = 1..No] = id# of kth node in (1-based list)  list of all nodes of conn.ord. co(o). For o=0: No=N1L ; o=1: No=N1R ; o=2..8: No=N[o] :
extern Index_vector &RootEnds, &LeafEnds ; // (= i_[1] and i_[0], resp.) : label which nodes are network ends = nodes of co=1 , which are either 'leaf' or 'root' ends :
extern vector<int> jf_RootEnds ; // to index those nodes that are root- or leaf-ends, resp.
extern vector<vector<int> > JF_ ; // JF_[i][j] = Id# du (j)eme flux (dans l'ordre declare dans le fichier) issu du noeud #i :> 0 ou < 0 par convention dans le sens (I_Upflow[ ] --> I_Downflow[ ]) :
extern vector<int> I_Upflow, I_Downflow ; // I_Upflow(resp.I_Downflow)[jf=1..Nc] = id# du noeud amont (resp. aval) par convention : jf > 0 dans le sens (I_Upflow[(abs(jf)] --> I_Downflow[(abs(jf)], < 0 dans l'autre sens

/******************* VARIABLES INVOLVED in CARBON METABOLISM AND FLUXES *******************************  */
extern Fortran_vector JS_ST						; // (mmol / h)  : Axial phloem sugar flux
extern Fortran_vector JS_PhlMb				; // Phloem cross-membrane sugar fluxes into sieve tubes from apoplasm				(mmol / h)
extern Fortran_vector RespMaint					; // Maintenance respiration rate										(mmol / h)
extern Fortran_vector Q_RespMaintSyn						; // Michaelis-Menten rate of starch synthesis from sugar substrate						(mmol sug.eq./ h)
extern Fortran_vector C_amont					; //  (mmol / ml) : ST Sugar concentration at upflow node
extern Fortran_vector Input					; // Local C input (photosynthetic assimilation rate in leaves, but may be defined anywhere) (boundary condition)			(mmol / h)
extern Fortran_vector C_Sympl						; // Concentration of sugar in parenchyma = Q_Par / vol_Sympl			(mmol / ml)
extern Fortran_vector C_ST							; // Concentration of sugar in sieve tubes								(mmol / ml solution))
//extern Index_vector i_amont					; // i_amont[i] = Id# of upflow node  ; def. ONLY FOR i > 1  (i_amont[1] = NA)
extern Fortran_vector i_amont					;
extern Fortran_vector JS_Sympl						; // Symplasmic flux of sugar from Lateral parenchyma to phloem ST (mmol / h)
extern Fortran_vector JS_ParMb				; // Lateral parenchyma cross-membrane sugar flux into symplasm from apoplasm				(mmol / h)
extern Fortran_vector JS_Apo						; // apoplasmic sugar flux from phloem to Lateral parenchyma  (mmol / h)
extern Fortran_vector JW_ParMb, JW_Apo, JW_Sympl ; // water fluxes corresponding to above 3 sugar fluxes   (ml / h)
extern Fortran_vector C_SymplUpflow					; // upflow concentration (mmol / ml) for JS_Sympl
extern Fortran_vector C_PhlApo, C_ParApo, C_ApoUpflow ; // (mmol / ml) apoplasmic sugar conc., resp. in phloem and lat.parenchyma, and upflow conc. for JS_Apo
extern Fortran_vector Delta_JS_ST ; // sera la composante purement phloemienne de Q_TC_dot[ ]							(mmol / h)

/******** TRACER-RELATED VARIABLES ***********************************************************/
extern double TracerDecay_k, TracerHalfLife ;
extern Fortran_vector TracerJS_ST						; // TracerJS_ST[i] (MBq / h)  = TracerJS_ST_[i-1], i = 2..N  : Axial phloem soluble tracer flux ; TracerJS_ST[1]= NA, whereas TracerJS_ST_[1] = TracerJS_ST[2]
extern Fortran_vector TracerJS_PhlMb						; // Lateral (Apoplasmic) soluble tracer fluxes into sieve tubes from parenchyma				(MBq / h)
extern Fortran_vector TracerRespMaint					; // Tracer Maintenance respiration rate										(MBq / h)
extern Fortran_vector TracerQ_RespMaintSyn						; // Michaelis-Menten rate of tracer starch synthesis from tracer sugar substrate						(MBq / h)
extern Fortran_vector TracerInput					; // Tracer input (photosynthetic Tracer Assimilation rate in leaves, but more general) (boundary condition)			(MBq / h)
extern Fortran_vector TracerC_Sympl						; // Concentration of soluble tracer in parenchyma symplasm= TracerQ_Mesophyll / vol_Sympl			(MBq / ml)
extern Fortran_vector TracerC_ST							; // Concentration of soluble tracer in sieve tubes								(MBq / ml solution))
extern Fortran_vector TracerC_PhlApo, TracerC_ParApo ; // Concentration of soluble tracer in apoplasms					(MBq / ml solution))
extern Fortran_vector TracerRatioSympl ; //   TracerQ_Mesophyll / Q_Mesophyll = TracerC_Sympl / C_Sympl  (MBq / mmol)
extern Fortran_vector TracerRatioQ_RespMaint ; //   = TracerQ_RespMaint / Q_RespMaint (MBq / mmol)
extern Fortran_vector Delta_TracerJS_ST ; // sera la composante purement phloemienne de Q_Rmmax_dot[ ]
extern Fortran_vector TracerJS_Sympl, TracerJS_Apo, TracerJS_ParMb ; // Lateral (Sympl., Apopl.; cross-membr...) soluble tracer fluxes FROM sieve tubes INTO parenchyma							(MBq / h)     !!! ATTENTION : sens positif oppose a JS_Trsv !!!
extern Fortran_vector TracerC_SymplUpflow					;  // upflow tracer concentration (mmol / ml) for TracerJS_Sympl
extern Fortran_vector TracerC_ApoUpflow ;				; // upflow tracer concentration (mmol / ml) for TracerJS_Apo
extern Fortran_vector Delta_TracerJS_ST ; // sera la composante purement phloemienne de Q_Rmmax_dot[ ]

extern Fortran_vector P_Xyl			; // Xylem water potential = pressure (no solute)												(MPa)
extern Fortran_vector Psi_Xyl			; // Xylem water potential = pressure (no solute)												(MPa)
extern Fortran_vector Psi_ST			; // Phloem sieve-tube (+CC) water potential												(MPa)
extern Fortran_vector Psi_PhlApo		; // Phloem Apoplam water potential												(MPa)
extern Fortran_vector Psi_ParApo		; // Lateral parenchyma apoplasmic water potential												(MPa)
extern Fortran_vector P_PhlApo		; // Phloem Apoplam pressure												(MPa)
extern Fortran_vector P_ParApo		; // Lateral parenchyma apoplasmic pressure												(MPa)
extern Fortran_vector Psi_Sympl		; // Lateral parenchyma symplasmic water potential												(MPa)
extern Fortran_vector P_ST			; // Phloem sieve-tube (+CC) turgor pressure 								(MPa)
extern Fortran_vector P_Sympl		; // Lateral parenchyma symplasmic turgor pressure 								(MPa)
extern Fortran_vector Absorb			; // Root end water absorbtion flux										(ml / h)
extern Fortran_vector JW_Xyl			; // Axial xylem water flux
extern Fortran_vector JW_ST			; // Axial phloem sieve-tube liquid flux
extern Fortran_vector JW_Trsv			; // Transverse xylem to phloem apoplasm water flux									(ml / h)
extern Fortran_vector JW_PhlMb		; // Transmembrane Phloem (ST,CC) Apoplasm to Symplasm water flux
extern Fortran_vector JW_ParMb		; // Transmembrane Parenchyma Apoplasm to Symplasm water flux
extern Fortran_vector JW_Apo			; // Phloem to Lateral parenchyma Apoplastic water flux									(ml / h)
extern Fortran_vector JW_Sympl		; // Lateral parenchyma to phloem ST Symplasmic liquid flux									(ml / h)
// Next 2 variables are not considered as such, but as possible inputs to compute vol_Sympl_dot :
extern Fortran_vector P_ST_dot, P_Sympl_dot			; //  dP_ST/dt , dP_Sympl/dt					(MPa h / h)    -- for elasticity...
// NZS : optional non-zero volume sugar flow (not a distinct variable)   (ml / h) : NZS = JS_Trsv * PartMolalVol (=0.2155 in Thompson and Holbrook -- 0.214 might be more accurate)

/******************************************  Constants and Parameters: *********************************************/
extern double TdC, dEauPure, PartMolalVol, siPhi, newPhi ; // pour visc. calc. par  www.seas.upenn.edu
extern bool Adv_BioPhysics ; // true if  non-zero sugar specific volume, osmotic pot.=non-linear function of molality (Thompson and Holbrook), and viscosity changes with C_TC (Thompson and Holbrook ; Seas, Flanagan)  ; set in IntroDialogBox
extern Fortran_vector Q_ST_seg_init;
// Hydro Resistance Parameters (set in GUI) :
extern Fortran_vector r_Xyl		; // (MPa h / ml)  : axial xylem water resistance
extern Fortran_vector r_ST			; // (MPa h / ml) : axial phloem water resistance
extern Fortran_vector r_ST_ref 	; // r_ST_ ref. values for C_TC =  0.5 mmol / ml sap solution
extern Fortran_vector r_abs			; // soil - root  water resistances (may be changed by user)			(MPa h / ml)
extern Fortran_vector r_Trsv			; // Transverse xylem to phloem apoplasm water resistance										(MPa h / ml)
extern Fortran_vector r_PhlMb		; // Transmembrane Phloem (ST,CC) Apoplasm to Symplasm water resistance		(MPa h / ml)
extern Fortran_vector r_ParMb		; // Transmembrane Parenchyma Apoplasm to Symplasm water resistance		(MPa h / ml)
extern Fortran_vector r_Apo			;// Phloem to Lateral parenchyma Apoplastic water resistence (MPa h / ml)
extern Fortran_vector r_Sympl			; // Lateral parenchyma to phloem ST Symplasmic water resistance (MPa h / ml)

extern int resistance_changed ;
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


extern Fortran_vector  Q_Exud, Q_Gr, Q_Rm, Q_Fl; 
extern Fortran_vector  Q_Exudmax, Q_Grmax, Q_Rmmax;//, Q_Fl_seg;
//extern Fortran_vector C_ST_seg, Q_ST_seg, Q_ST_seg_new,Q_ST_newFuFl, grad_Q_ST_seg;

/******************************************  Environmental Variables: *********************************************/
extern double T ;	// absolute temperature													(K)
// The following 2 environment data are set by the user in function  'hydric_boundary_conditions()'  at the end of file 'water-fluxes.cpp' :
extern Fortran_vector Transpirat		; // Leaf transpiration rate (used in water-fluxes calc.)				(mmol / h)
extern Fortran_vector PsiSoil			; // Soil water potential at root end (used in water-fluxes calc.)		(MPa)
void UpdateResistances(double t) ;
extern Fortran_vector Y0	; // (set in GUI) initial condition vector is made of Q_ST_0, Q_Mesophyll_0, Q_RespMaint_0, Q_Exudation_0, Q_Growthtot_0, and homologous tracer init values, and vol_Sympl_0 :




void PhloemFlux::initialize_carbon(vector<double> vecIn) {
	
	JS_ST = Fortran_vector(Nc, 0.) ; // (mmol / h)   Axial phloem sugar flux
    //JS_PhlMb = Fortran_vector(Nt, 0.)			; // CrossMembrane phloem sugar fluxes from apoplasm into sieve tubes (mmol / h)
    RespMaint = Fortran_vector(Nt, 0.)			; // Maintenance respiration rate										(mmol / h)
    //Q_RespMaintSyn = Fortran_vector(Nt, 0.)				; // Rate of starch synthesis from sugar substrate						(mmol sug.eq./ h)
    Input = Fortran_vector(Nt, 0.)		; // External sugar input (may be photosynthetic Assimilation rate, but not restricted to leaves)	(mmol / h)
    i_amont = Fortran_vector(Nc, 0.)	; //  Index_vector(Nc)	; // true upflow node : sera I_Upflow[j]  ou  I_Downflow[j] suivant le sens reel du flux
    C_amont = Fortran_vector(Nc, 0.)	; //  (mmol / ml) : ST Sugar concentration at true upflow node
	C_ST = Fortran_vector(Nt, 0.);
	if(doTroubleshooting){
		std::cout<<"initial size of vector: "<<vecIn.size()<<" nodes: "<<Nt<<" connections "<<Nc<<" "<<std::endl;
	}
    if(vecIn.size() == (Nt*neq_coef)){ //gave input vector with starting values 
		if(doTroubleshooting){std::cout<<"setup full y0 "<<std::endl;}
		Y0= Fortran_vector(vecIn); 
		
		//Y0.display();
		Fortran_vector Q_GrowthtotBU_temp = Fortran_vector(Nt - Nt_old, 0.) ;
		Q_GrowthtotBU.append(Q_GrowthtotBU_temp);
		Fortran_vector Q_GrmaxBU_temp = Fortran_vector(Nt - Nt_old, 0.) ;
		Q_GrmaxBU.append(Q_GrmaxBU_temp);
	}else{
		if(vecIn.size() > 0){// plant grew since last phloem flow computatoin
			if(doTroubleshooting){std::cout<<"complete y0 "<<std::endl;}
			Y0 =  Fortran_vector(Nt*neq_coef, 0.) ;
			Y0.sequentialFill(vecIn, Nt_old, Nt);
			
			Fortran_vector Q_GrowthtotBU_temp = Fortran_vector(Nt - Nt_old, 0.) ;
			Q_GrowthtotBU.append(Q_GrowthtotBU_temp);
			Fortran_vector Q_GrmaxBU_temp = Fortran_vector(Nt - Nt_old, 0.) ;
			Q_GrmaxBU.append(Q_GrmaxBU_temp);
		}else{//first phloem flow computation ==> vecIn is empty
			if(doTroubleshooting){std::cout<<"setup empty y0 "<<std::endl;}
			Y0 =  Fortran_vector(Nt*neq_coef, 0.) ;
			Q_GrowthtotBU = Fortran_vector(Nt, 0.) ;
			Q_GrmaxBU = Fortran_vector(Nt, 0.) ;
			if(withInitVal){//initial value for mesophyll and sieve tube
				for(int z = 0; z < Nt;z++){
					Y0[z + 1] = initValST * vol_ST[z + 1]; //conz to content
					Y0[z + 1  + Nt * 8] = initValST * vol_ST[z + 1]; //starch
				}
				for(int z = 0; z < Nt;z++){
					if(vol_ParApo[z + 1] > 0){
						Y0[z + 1 + Nt ] = initValMeso * vol_ParApo[z + 1];
					Y0[z + 1  + Nt * 7] = initValMeso * vol_ParApo[z + 1]; //starch in mesophil
					}
				}
				this->Q_init = Y0.toCppVector(); //for post processing
			}
		}
	}
	
	if(doTroubleshooting){cout<<"Y0_STinit "<<Y0[1]<<" "<<Nc<<" "<<Nt<<" "<<Nt_old<<endl;}
	Nt_old = Nt; //BU Nt
	
	Q_Exud = Fortran_vector(Nt, 0.)			; 
    Q_Gr = Fortran_vector(Nt, 0.)			; 
    Q_Rm = Fortran_vector(Nt, 0.)			; 
	
    Q_Fl = Fortran_vector(Nt, 0.)			; 
	
    }

void PhloemFlux::initialize_hydric() {
	int j;//, k ;
	//jf_RootEnds =vector<int>(N1R + 1) ; for(k = 1 ; k <= N1R ; k ++)  jf_RootEnds[k] = JF_[RootEnds[k]][1] ; // donne le n# (positif si part du noeud considere, negatif s'il y arrive) du flux reliant le noeud a son voisin
	P_Xyl = Fortran_vector(Nt, 0.)	; // Xylem water potential (initial value in relation to resistance update criterion...) (MPa)
	Psi_ST = Fortran_vector(Nt, 0.)		; // (id.) Phloem water potential	
	//									(MPa)
	P_ST = Fortran_vector(Nt, 0.)			; // Sieve tube turgor pressure											(MPa)
	JW_Trsv = Fortran_vector(Nt, 0.)		; // Transverse xylem-phloem water flux									(ml / h)
	JW_Xyl = Fortran_vector(Nc, 0.)		; // (ml / h)  :  Axial xylem water flow
	JW_ST = Fortran_vector(Nc, 0.)		; // (ml / h) : Axial phloem liquid flow
	
	JW_PhlMb = Fortran_vector(Nt, 0.)		; // (ml / h) : Transmembrane Apo->(ST,CC) phloem water flow
	JW_ParMb = Fortran_vector(Nt, 0.)		; // (ml / h) : Transmembrane parenchyma (Apo->Sympl) water flow
	Psi_Sympl = Fortran_vector(Nt, 0.)			; // Parenchyma Symplasmic water potential												(MPa)
	P_Sympl = Fortran_vector(Nt, 0.)					; // Parenchyma Symplasmic turgor pressure 								(MPa)
	JW_Apo = Fortran_vector(Nt, 0.)				; // (Transverse) Apoplasmic (phloem -> parenchyma) water flux									(ml / h)
	JW_Sympl = Fortran_vector(Nt, 0.)				; // (Transverse) Symplasmic (parenchyma -> phloem) water flux									(ml / h)
	Psi_PhlApo = Fortran_vector(Nt, 0.)		; // (id.) Phloem Apoplstic water potential										(MPa)
	Psi_ParApo = Fortran_vector(Nt, 0.)		; // (id.) Parenchyma Apoplstic water potential										(MPa)
	P_PhlApo = Fortran_vector(Nt, 0.)		; // (id.) Phloem Apoplstic pressure										(MPa)
	P_ParApo = Fortran_vector(Nt, 0.)		; // (id.) Parenchyma Apoplstic pressure										(MPa)
	
	P_ST_dot = Fortran_vector(Nt, 0.) ;   // ajout pour elasticite (dP/dt)
	P_Sympl_dot = Fortran_vector(Nt, 0.) ;   // ajout pour elasticite (dP_Par/dt)
	r_ST = r_ST_ref		; //(MPa h / ml) :  phloem water resistance
	//add later
	//if (Adv_BioPhysics) PartMolalVol = 0.2155 ; else 
	PartMolalVol = 0. ; // (L / mol) valeur Thompson et Holbrook ; en fait serait plutôt 0.214, independant (a 0.1% pres) a la fois de T et de C (cf. SucroseViscosity.xls)
	// pour visc. calc. par  www.seas.upenn.edu :
	TdC = T - 273.15;
	dEauPure = (999.83952 + 16.952577 * TdC - 7.9905127 * (0.001) * (TdC*TdC) - 46.241757 * (0.000001) * (TdC*TdC*TdC) + 105.84601 * (0.000000001) * (TdC*TdC*TdC*TdC) - 281.03006 * (0.000001*0.000001) * (TdC*TdC*TdC*TdC*TdC)) / (1 + 16.887236 * (0.001) * TdC); // g/L
	siPhi = (30 - TdC) / (91 + TdC);
	newPhi=( - 0.114 + (siPhi *1.1));
	Delta_JS_ST = Fortran_vector(Nt, 0.) ;
	Delta_TracerJS_ST = Fortran_vector(Nt, 0.) ;
	// les signes ci-dessous sont donnes en coherence avec eq. (1) a (4), i.e. considerant que j=JW_ST (et sens oppose a JW_Xyl) :
	Sparse_matrix* Delta2_ = new Sparse_matrix(Nt, Nc, 2*Nc) ;//IM shape : ligne shape = Nt; column.size = Nc, 
	Sparse_matrix* Delta2_abs = new Sparse_matrix(Nt, Nc, 2*Nc) ;//IM shape : ligne shape = Nt; column.size = Nc, 
	
	for(j = 1 ; j <= Nc ; j ++) { 
		Delta2_->set_(I_Upflow[j], j, -1.) ; 
		Delta2_->set_(I_Downflow[j], j, 1.) ; 
		
		Delta2_abs->set_(I_Upflow[j], j, 0.5) ; 
		Delta2_abs->set_(I_Downflow[j], j, 0.5) ; 	
		}
	Delta2 = SpUnit_matrix(*Delta2_) ; delete Delta2_ ;
	Sparse_matrix Delta2abs = *Delta2_abs ; delete Delta2_abs ;

	Delta = -transpose(Delta2) ;
	Deltaabs = transpose(Delta2abs) ;
	Delta_rphl = Delta / r_ST_ref ;
			
	ipiv_ptr = NULL ; TM_ptr = NULL;
	JW_Xyl.set(0.) ; //JW_Xyl.sub_matmult(Delta_rxyl, P_Xyl) ; // JW_Xyl = - (1/r_Xyl * Delta) · P_Xyl (eq. 4)
	Absorb.set(0.) ;  P_PhlApo.set(0.) ; //.set(P_Xyl) ; P_PhlApo.sub_elemult(r_Trsv, JW_Trsv) ; // P_PhA = PX - rT * JT (eq. 6)
    Psi_PhlApo.set(0.) ; //.set(P_PhlApo) ; Psi_PhlApo.sub_mult(RT, C_PhlApo) ; // PsiPhA = P_PhA - RT C_PhA  (eq. 6')
    JW_PhlMb.set(0.) ; //.set_sub(Psi_PhlApo, Psi_ST) ; JW_PhlMb *= inv_rPhlM ; // (ST,CC) transmembrane flow (eq. 7)
    JW_Apo.set(0.) ; //.set_sub(JW_Trsv, JW_PhlMb) ;  JW_Apo.zero(1.e-10 * JW_Trsv) ; // (loi des noeuds (eq. 14) -- et elimination des erreurs d'arrondi a la soustraction
    P_ParApo.set(0.) ; //.set(P_PhlApo) ; P_ParApo.sub_elemult(r_Apo, JW_Apo) ; // eq. (9)
    Psi_ParApo.set(0.) ; //.set(P_ParApo) ; Psi_ParApo.sub_mult(RT, C_ParApo) ; // eq. (9')
    JW_ParMb.set(0.) ; //.set(JW_Apo) ;  // (eq. 13)
    Psi_Sympl.set(0.) ; //.set(Psi_ParApo) ; Psi_Sympl.sub_elemult(r_ParMb, JW_ParMb) ; // eq. (10)
    P_Sympl.set(0.) ; // += Psi_Sympl ; // (eq. 12)  --   'P_Sympl' valait provisoirement  (RT * M_Sympl)
	JW_Sympl.set(0.) ; //.set_sub(JW_ParMb, dummy) ; // (eq. 8) dummy vaut encore  (volSymplDot - NZSP).
	
}

void UpdateResistances(double t) { // called if (resistance_changed)
	//eventually fill it to represent changes of resistances when adding compartment
	resistance_changed = NONE ;
}
