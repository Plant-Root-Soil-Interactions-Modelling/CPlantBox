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

#include <PiafMunch/runPM.h>

vector<int> I_Upflow, I_Downflow ; //  I_Upflow(resp.I_Downflow)[jf=1..Nc] = id# of conv. upflow (resp. conv. downflow) node of internode flux #jf
int i,j;
float rs = 1.1f ; // (if relevant) rate of geometric increment of output step -- to be adjusted to yield 100 steps between t0 and tf
extern bool Adv_BioPhysics ; 
extern double dEauPure;
string term = "qt" ; // name of terminal for gnuplot (default "wxt" bugs with Windows)
extern string GnuplotPath ; // (defined in PiafMunch2.cpp) Path to GnuPlot executable, including executable file name
extern Fortran_vector r_ST_ref 	; // r_ST_ ref. values for C_TC =  0.5 mmol / ml sap solution

// Network Architecture :
extern int N1L ; // Organ type # 1L : 'leaf end' ; has imposed transpiration stream as a limit condition
extern int N1R ; // Organ type # 1R : 'root tip' ; has imposed soil water potential as a limit condition
// N[o = 1..8] = number of topological nodes of connectivity order o : N[1] = N1L + N1R ; N[0] = N[1] + N[2] + N[3] + ... + N[8] = Nt :
extern int Nt, Nc, N[] ; // Nt = total number of nodes, Nc = total number of internode connections [= (1*N[1] + 2*N[2] + 3*N[3] + ... + 8*N[8]) / 2]
extern Index_vector i_[] ; // i_[o = 0..8][k = 1..No] = id# of kth node in (1-based list)  list of all nodes of conn.ord. co(o). For o=0: No=N1L ; o=1: No=N1R ; o=2..8: No=N[o] :
extern Index_vector &RootEnds, &LeafEnds ; // (= i_[1] and i_[0], resp.) : label which nodes are network ends = nodes of co=1 , which are either 'leaf' or 'root' ends :
extern std::vector<int> I_Upflow, I_Downflow ; // I_Upflow(resp.I_Downflow)[jf=1..Nc] = id# du noeud amont (resp. aval) : jf = JF(i,i2) > 0 si I_Upflow[(abs(jf)]==i, i.e. si I_Downflow[(abs(jf)]==i2
extern Fortran_vector kML						; // kinetic parameter / Michaelis - phloem loading					(mmol / ml)
extern Fortran_vector vMU					; // kinetic parameter / phloem unloading								(mmol /h)
//extern Fortran_vector	isTip ,Ag,  Lmax_org, Rmax_org, krm1, krm2 , StructC, exud_k; 
extern Fortran_vector radius_ST ,  Length, vol_Seg;// Q_Grmax, Q_Rmmax, Q_Exudmax; 

extern Fortran_vector Ag, Q_Grmax, Q_Rmmax, Q_Exudmax, exud_k, krm2, len_leaf; 
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

std::weak_ptr<PhloemFlux> phloem_; 
vector<double> extra_output_times, breakpoint_times ; // may contain already existing output time points
extern int is, solver ; // (see below) solver config.# -- default = 1 = cvode (SPILS: SPFGMR, MODIFIED_GS, PREC_NONE... Try diff. value if calc. fails or too slow !
void auxout(double t, double * y){phloem_.lock()->aux(t, y);} ;	// launch auxiliary calculations to store dynamics of temporary variables
void fout(double t, double *y, double *y_dot){phloem_.lock()->f(t, y, y_dot);} ; // the function to be processed by the solver  (implemented in 'solve.cpp')

// Full list of printable/savable variables -- to be updated in case of hard-updating the model (number and nature of local compartments within each archit. element) :
int NumAllNodeVariablesNames = 62 ; // 62 node variables :
string AllNodeVariablesNamesQList[] = {
	"Absorb (ml / h)", "Q_RespMaint (mmol eq. Glu)", "Q_RespMaintSyn (mmol eq. Glu / h)", "C_ST (mmol / ml)", "C_Sympl (mmol / ml)", "JS_PhlMb (mmol / h)", "JW_Trsv (ml / h)", "JW_PhlMb (ml / h)", "JW_ParMb (ml / h)", "Psi_ST (MPa)",
	"P_Xyl (MPa)", "Q_Mesophyll (mmol)", "Q_ST (mmol)", "Resp_Maint (mmol / h)", "C_SymplUpflow (mmol / ml)", "JS_Sympl (mmol / h)", "JW_Apo (ml / h)", "JW_Sympl (ml / h)", "P_ST (MPa)", "P_Sympl (MPa)",
	"Psi_PhlApo (MPa)", "Psi_ParApo (MPa)", "Psi_Sympl (MPa)", "vol_Sympl (ml)", "vol_Sympl_dot (ml / h)", "P_PhlApo (MPa)", "P_ParApo (MPa)", "Q_Exudation (mmol)", "Q_Growthtot (mmol)",
	"C_PhlApo (mmol / ml)", "C_ParApo (mmol / ml)", "Q_Exud_dot (mmol / h)", "Q_Gtot_dot (mmol / h)", "JS_Apo (mmol / h)", "JS_ParMb (mmol / h)", "C_ApoUpflow (mmol / ml)",
	"P_ST_dot (MPa / h)", "P_Sympl_dot (MPa / h)", "Q_Rm_dot (mmol eq. Glu / h)", "Q_ST_dot (mmol / h)", "Q_Mesophyll_dot (mmol / h)", "Input (mmol / h)", "Transpirat (ml / h)", "PsiSoil (MPa)",
// Tracer-specific add-ins :
	"TracerQ_RespMaint (MBq)", "TracerQ_RespMaintSyn (MBq / h)", "TracerC_ST (MBq / ml)", "TracerC_Sympl (MBq / ml)", "TracerJS_PhlMb (MBq / h)", "TracerQ_Mesophyll (MBq)", "Q_RespMaintmax (MBq)", "TracerRespMaint (MBq / h)",
	"TracerInput (MBq / h)", "TracerJS_Sympl (MBq / h)", "TracerJS_ParMb (MBq / h)", "Q_S_Mesophyll (MBq)", "Q_Growthtotmax (MBq)", "TracerC_PhlApo (MBq / ml)",
	"TracerC_ParApo (MBq / ml)", "TracerJS_Apo (MBq / h)", "TracerC_SymplUpflow (MBq / ml)", "TracerC_ApoUpflow (MBq / ml)"
} ;
int NumAllConnVariablesNames = 5 ; // 5 internode connector flux variables :
string AllConnVariablesNamesQList[] = { "C_Upflow (mmol / ml)", "JS_ST (mmol / h)", "JW_ST (ml / h)", "JW_Xyl (ml / h)", "TracerJS_ST (MBq / h)" } ;

extern double *Q_ST, *Q_Mesophyll, *Q_RespMaint, *Q_Exudation, *Q_Growthtot, *Q_S_ST, *Q_Mucil ;		  // components of vector y as used in diff. system f()...
extern double *Q_ST_dot, *Q_Mesophyll_dot, *Q_Rm_dot, *Q_Exud_dot, *Q_Gtot_dot, *Q_S_ST_dot, *Q_Mucil_dot ; //... and its derivatives.  ;
extern double *Q_RespMaintmax, *TracerQ_Mesophyll, *TracerQ_RespMaint, *Q_S_Mesophyll, *Q_Growthtotmax ;		  // components of vector y as used in diff. system f()...
extern double *Q_Rmmax_dot, *TracerQ_Mesophyll_dot, *TracerQ_Rm_dot, *Q_S_Mesophyll_dot, *Q_Gtotmax_dot ; //... and its derivatives.  ;
extern double *vol_Sympl, *vol_Sympl_dot ;
extern Fortran_vector Psi_Xyl, JW_ST, JS_ST, C_Sympl, C_ST, C_amont, JS_Sympl, JS_Apo, JS_ParMb, JS_PhlMb, Psi_Sympl, C_PhlApo, C_ParApo, P_ST, Absorb, RespMaint ;
extern Fortran_vector JW_Trsv, JW_Xyl, JW_PhlMb, JW_ParMb, JW_Apo, JW_Sympl, P_Sympl, Psi_ST, P_Xyl, Psi_PhlApo, Psi_ParApo, P_PhlApo, P_ParApo ;
extern Fortran_vector r_Xyl, r_ST, r_ST_ref, r_Trsv, r_PhlMb, r_ParMb, r_Apo, r_Sympl, vol_ST, vol_PhlApo, vol_ParApo ;
extern Fortran_vector r_abs, PsiSoil, Transpirat, Input, Q_RespMaintSyn ;
extern Fortran_vector P_ST_dot, P_Sympl_dot ; // variation rate of any variable X is noted: X_dot = dX/dt
extern Fortran_vector C_SymplUpflow, C_ApoUpflow ; // to derive JS_Sympl, and maybe JS_Apo, from resp. water fluxes JW_xxx
// tracer-specific add-ins :
extern Fortran_vector TracerJS_ST, TracerC_Sympl, TracerC_ST, TracerJS_Sympl, TracerJS_Apo, TracerJS_ParMb, TracerJS_PhlMb, TracerC_PhlApo, TracerC_ParApo ;
extern Fortran_vector TracerQ_RespMaintSyn, TracerInput, TracerRespMaint, TracerC_SymplUpflow, TracerC_ApoUpflow ;
extern Fortran_vector TracerRatioSympl, TracerRatioQ_RespMaint ;

// variables pointers are of type 'Fortran_vector*' (or NULL for those that are components of double* y in definition of function f).
// Check for consistency between these lists AllxxxxVariablesRefsQList hereafter and AllxxxxVariablesNamesQList[] above:
Fortran_vector* AllNodeVariablesRefsQList[] = {  // 62 node variables :
	&Absorb, NULL, &Q_RespMaintSyn, &C_ST, &C_Sympl, &JS_PhlMb, &JW_Trsv, &JW_PhlMb, &JW_ParMb, &Psi_ST,
	&P_Xyl, NULL, NULL, &RespMaint, &C_SymplUpflow, &JS_Sympl, &JW_Apo, &JW_Sympl , &P_ST, &P_Sympl,
	&Psi_PhlApo, &Psi_ParApo, &Psi_Sympl, NULL, NULL, &P_PhlApo, &P_ParApo, NULL, NULL,
	&C_PhlApo, &C_ParApo, NULL, NULL, &JS_Apo, &JS_ParMb, &C_ApoUpflow,
	&P_ST_dot, &P_Sympl_dot,  NULL, NULL, NULL, &Input, &Transpirat, &PsiSoil,
	// tracer-specific add-ons :
	NULL, &TracerQ_RespMaintSyn, &TracerC_ST, &TracerC_Sympl, &TracerJS_PhlMb, NULL, NULL, &TracerRespMaint,
	&TracerInput, &TracerJS_Sympl, &TracerJS_ParMb, NULL, NULL, &TracerC_PhlApo,
	&TracerC_ParApo, &TracerJS_Apo , &TracerC_SymplUpflow, &TracerC_ApoUpflow
} ;
Fortran_vector* AllConnVariablesRefsQList[] = {&C_amont, &JS_ST,  &JW_ST, &JW_Xyl, &TracerJS_ST} ; // 5 internode-connector variables

// temporary storage of results for output :
Fortran_matrix Q_ST_t, Q_Mesophyll_t, Q_RespMaint_t, Q_Exudation_t, Q_Growthtot_t, Q_ST_dot_t, Q_Mesophyll_dot_t, Q_Rm_dot_t, Q_Exud_dot_t, Q_Gtot_dot_t, vol_Sympl_t, vol_Sympl_dot_t ;
Fortran_matrix JW_ST_t, JS_ST_t, C_Sympl_t, C_ST_t, JS_Sympl_t, JS_Apo_t, JS_ParMb_t, JS_PhlMb_t, Psi_Sympl_t, C_PhlApo_t, C_ParApo_t, P_ST_t, Absorb_t, RespMaint_t ;
Fortran_matrix JW_Trsv_t, JW_Xyl_t, JW_PhlMb_t, JW_ParMb_t, JW_Apo_t, JW_Sympl_t, P_Sympl_t, Psi_ST_t, P_Xyl_t, Psi_PhlApo_t, Psi_ParApo_t, P_PhlApo_t, P_ParApo_t ;
Fortran_matrix PsiSoil_t, Transpirat_t, Input_t, Q_RespMaintSyn_t, P_ST_dot_t, P_Sympl_dot_t, C_SymplUpflow_t, C_ApoUpflow_t, C_amont_t ;
// tracer-specific add-ins :
Fortran_matrix Q_RespMaintmax_t, TracerQ_Mesophyll_t, TracerQ_RespMaint_t, Q_S_Mesophyll_t, Q_Growthtotmax_t ;		  // components of vector y as used in diff. system f()...
// Fortran_matrix Q_Rmmax_dot_t, TracerQ_Mesophyll_dot_t, TracerQ_Rm_dot_t, Q_S_Mesophyll_dot_t, Q_Gtotmax_dot_t ;		  // components of vector y as used in diff. system f()...
Fortran_matrix TracerJS_ST_t, TracerC_Sympl_t, TracerC_ST_t, TracerJS_Sympl_t, TracerJS_Apo_t, TracerJS_ParMb_t, TracerJS_PhlMb_t, TracerC_PhlApo_t, TracerC_ParApo_t, TracerC_SymplUpflow_t, TracerC_ApoUpflow_t ;
Fortran_matrix TracerPhotSynth_t, TracerQ_RespMaintSyn_t, TracerInput_t, TracerRespMaint_t ;
// Fortran_matrix TracerRatioSympl_t, TracerRatioQ_RespMaint_t ;

Fortran_matrix* AllNodeMatricesRefsQList[] = {
	&Absorb_t, &Q_RespMaint_t, &Q_RespMaintSyn_t, &C_ST_t, &C_Sympl_t, &JS_PhlMb_t, &JW_Trsv_t, &JW_PhlMb_t, &JW_ParMb_t, &Psi_ST_t,
	&P_Xyl_t, &Q_Mesophyll_t, &Q_ST_t, &RespMaint_t, &C_SymplUpflow_t, &JS_Sympl_t, &JW_Apo_t, &JW_Sympl_t, &P_ST_t, &P_Sympl_t,
	&Psi_PhlApo_t, &Psi_ParApo_t, &Psi_Sympl_t, &vol_Sympl_t, &vol_Sympl_dot_t, &P_PhlApo_t, &P_ParApo_t, &Q_Exudation_t, &Q_Growthtot_t,
	&C_PhlApo_t, &C_ParApo_t, &Q_Exud_dot_t, &Q_Gtot_dot_t, &JS_Apo_t, &JS_ParMb_t, &C_ApoUpflow_t,
	&P_ST_dot_t, &P_Sympl_dot_t, &Q_Rm_dot_t, &Q_ST_dot_t, &Q_Mesophyll_dot_t, &Input_t, &Transpirat_t, &PsiSoil_t,
	// tracer-specific add-ons :
	&TracerQ_RespMaint_t, &TracerQ_RespMaintSyn_t, &TracerC_ST_t, &TracerC_Sympl_t, &TracerJS_PhlMb_t, &TracerQ_Mesophyll_t, &Q_RespMaintmax_t, &TracerRespMaint_t,
	&TracerInput_t, &TracerJS_Sympl_t, &TracerJS_ParMb_t, &Q_S_Mesophyll_t, &Q_Growthtotmax_t, &TracerC_PhlApo_t,
	&TracerC_ParApo_t, &TracerJS_Apo_t, &TracerC_SymplUpflow_t, &TracerC_ApoUpflow_t
};
Fortran_matrix* AllConnMatricesRefsQList[] = {&C_amont_t, &JS_ST_t,  &JW_ST_t, &JW_Xyl_t, &TracerJS_ST_t} ;
Index_vector Breakpoint_index(1, 1) ; // aux. for output time settings
double * y_dot = NULL ;	// used in f() as called by aux() which is called by odesolve() -- contains all derivatives of vector Y (whose starting values are Y0 below) below.
Fortran_vector Y0	; // (set in GUI) initial condition vector is made of Q_ST_0, Q_Mesophyll_0, Q_RespMaint_0, Q_Exudation_0 & Q_Growthtot_0, and homologous tracer init values, and vol_Sympl_0 :
Fortran_vector atol_, rtol ; // integration accuracy parameters (set in GUI)
double t0, tf, pas;		// solve for t = [t0, tf], save result at each time step (hours)
extern double T;
float t1 ; int nbv ; // pas = initial output step = between first 2 output steps t0 and t1
bool LogScale = false ; // set to 'true' if log scale wanted (in 'PiafMunch2.cpp')
Fortran_vector OutputTimes, SegmentTimes, truc ; // will be vector of, resp.: all output times including user-added extra- or breakpoint times ; output times for current integration time segment ; and an aux.


extern int resistance_changed ;
extern int jf, i, j ; double x ; bool to_store(true) ; FILE* file = NULL ; extern string S ;
int nl = 0 ; // # of current line in matrices currently under construction in aux() -- one for each element SegmentTimes[nl] of time vector SegmentTimes
extern std::vector<int> k_ ; // k_[i = 1 : Nt] = index of node #i in list i_[o(i)]   where o(i) = 0 for LeafEnds (k_[i] = - k then) ; o(i) = co(i) else  -- for the moment, filled only for (o = 0 or 1)   --   k_[i] = 0 else

// Selected nodes/connectors and variables for plotting or saving output :
int nsp, nvp, nss, nvs, fsp, fvp, fss, fvs, rsp(0), lsp(0), xsp ;
list<int> SelectedSaveNodesQList ; list<string> SelectedSaveNodeVariablesNamesQList ;
list<int> SelectedSaveConnsQList ; list<string> SelectedSaveConnVariablesNamesQList ;
list<int> SelectedPlotNodesQList ; list<string> SelectedPlotNodeVariablesNamesQList ;
list<int> SelectedPlotConnsQList ; list<string> SelectedPlotConnVariablesNamesQList ;
list<Fortran_vector*> SelectedPlotNodeVariablesRefsQList ; list<Fortran_vector*> SelectedPlotConnVariablesRefsQList ;
list<int> SelectedPlotNodeVariablesNamesIndicesInQList ; list<int> SelectedPlotConnVariablesNamesIndicesInQList ;
list<Fortran_vector*> SelectedSaveNodeVariablesRefsQList ; list<Fortran_vector*> SelectedSaveConnVariablesRefsQList ;
list<int> SelectedSaveNodeVariablesNamesIndicesInQList ; list<int> SelectedSaveConnVariablesNamesIndicesInQList ;
// following matrices will store dynamics of variables to plot   -- filled in aux() :
list<Fortran_matrix*> SelectedPlotNodeMatricesRefsQList, SelectedPlotConnMatricesRefsQList ;
list<int> SelectedPlotNetworkEnds_k, SelectedSaveNetworkEnds_k ;

extern string IniTitle ; // Path/name of .ini configuration file
list<int>::iterator it, itk ; list<Fortran_vector*>::iterator itv ; list<string>::iterator its ; list<Fortran_matrix*>::iterator itm ;
Fortran_vector* FV ; Fortran_matrix* FM ; string QS ;
bool AutoQuit(false); // will be TRUEd if launched with '-q' (which will also disable Plotting)


int PhloemFlux::startPM(double StartTime, double EndTime, int OutputStep,double TairK, bool verbose, std::string filename) {
	//std::string nametroubleshootingfile = "outputPM"+"_"+ std::to_string(simTime)+".txt"
	std::ofstream out(filename);
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	if(doTroubleshooting){
		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
		cout<<"extainr TairK "<<TairK<<" "<<StartTime<<" "<<EndTime<<" "<<OutputStep<<" verbose "<<verbose<<std::endl;
	}
	
    t0  = StartTime; tf  = EndTime;	nbv = OutputStep; T = TairK; TairK_phloem = TairK;
	
	nl = 0;
	if(doTroubleshooting){
		std::cout <<"set out put vec "<<t0<<" "<<t1<<" "<<tf<<" "<<nbv<<std::endl;
	}
	initializePM_(tf-t0, TairK);
	initialize_carbon(this->Q_outv) ;		// sizes C-fluxes-related variable vectors
    initialize_hydric() ;		// sizes water-fluxes-related variable vectors and sets up hydric system
	
	// building up output times vector OutputTimes according to GUI seetings :
    pas = (tf-t0)/double(nbv);
	t1 = t0 + pas ; //nbv = (int)((tf-t0)/pas) ; // pas = initial output step = between first 2 output steps t0 and t1 as set in GUI
	if(doTroubleshooting){
		std::cout <<"set out put vec "<<t0<<" "<<t1<<" "<<tf<<" "<<pas<<" "<<nbv<<std::endl;}
	
	OutputSettings() ; // sets up file and graph outputs, including possible breakpoint- and/or extra-output- times as stated above

    // *************** SOLVING THE DIFFERENTIAL EQUATION SYSTEM ************************************* :
    std::size_t neq = Y0.size() ;							// number of differential eq. = problem size (= 8*Nt apres ajouts FAD)
    if(doTroubleshooting){cout<<"neq "<<neq<<" "<<Nc<<" "<<Nt<<endl;}
	
	assert((Nt == (neq/neq_coef))&&"Wrong seg and node number");
	assert((Nc == (Nt -1))&&"Wrong seg and node number");
	y_dot = new double[1 + neq] ;			// for use in aux()
    
    // -------------  To compute P_dot = dP/dt  and  P_symp_dot = dP_Sympl/dt  and make them available to all modules in real time : ---------------------
    // 1#) Oversize both 'Var_integrale' and 'Var_derivee' pointers and initialize all of them to NULL :
    Fortran_vector** Var_integrale = new Fortran_vector*[100] ; Fortran_vector** Var_derivee = new Fortran_vector*[100] ;
    for (i = 0 ; i < 100 ; i ++) {Var_integrale[i] = Var_derivee[i] = NULL ;}
    // 2#) initialize pointers to derivatives (and corresponding integrals) that are actually used (or may be so), in this case for...
    Var_integrale[0] = &P_Sympl ; Var_derivee[0] = &P_Sympl_dot ; //...elastic changes of symplastic volume in function (PiafMunch2.cpp)Smooth_Parameter_and_BoundaryConditions_Changes()
    Var_integrale[1] = &P_ST ; Var_derivee[1] = &P_ST_dot ;
		//throw std::runtime_error("Breakpoint_index.size() ");

	if(doTroubleshooting){
		std::cout <<"add pointer "<<std::endl;
	}
	phloem_ = this->Phloem();
    for(is = 1 ; is < Breakpoint_index.size() ; is ++) { // allows several integration segments in relation to breakpoints (if any -- OK if none)
        SegmentTimes = subvector(OutputTimes, Breakpoint_index(is), Breakpoint_index(is+1))  ; // time segment between 2 breakpoints (Breakpoint_index(1) = 1 ; Breakpoint_index(Breakpoint_index.size()) = 1 + nbv)
        if(doTroubleshooting){
			cout << endl << "starting integration on time segment #" << is << " = [" << SegmentTimes[1] << ", " << SegmentTimes[SegmentTimes.size()] << "] "<< Breakpoint_index.size()<<" "<<SegmentTimes.size()<< endl ;
        }
		// ***** the following solver configs are ranked from the most efficient (in most tested situations) to the least (in most tested situations). Can change in different situations ! *****
        //   see  SUNDIALS  documentation  for  cvode  solver options (SPxxxx, xxxx_GS, PREC_xxxx, BAND, etc.)
		if(doTroubleshooting){cout <<"solver, Y0: "<<solver <<endl;}
		switch (solver) {
			case 1: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, MODIFIED_GS, PREC_NONE, 2, Var_integrale, Var_derivee); break; // STALD = true, verbose = true, rootfind = NULL
			case 2: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, MODIFIED_GS, PREC_NONE, 2, Var_integrale, Var_derivee); break; // STALD = true, verbose = true, rootfind = NULL
			case 3: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, CLASSICAL_GS, PREC_NONE, 2, Var_integrale, Var_derivee); break; // le + rapide : best for large N (>= 1000) ; break ; sinon, un peu instable avec VarVisc
			case 4 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, CLASSICAL_GS, PREC_NONE, 2, Var_integrale, Var_derivee) ; break ; // le + rapide : best for large N (>= 1000) ; break ; sinon, un peu instable avec VarVisc
				case 5 : j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, DIAG, 2, Var_integrale, Var_derivee) ; break ; // TB pour Thompson, meme avec vol_Sympl_dot ...
			case 6 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPBCGS, 1, PREC_NONE, 2, Var_integrale, Var_derivee) ; break ;
				case 7: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, BAND, 2, Var_integrale, Var_derivee, true, true, NULL, 0, neq / 30, neq / 30); break; // mu = ml = neq/30 : TB
			case 8: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, PCG, 1, PREC_NONE, 2, Var_integrale, Var_derivee); break; //
			case 9 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPTFQMR, 1, PREC_NONE, 2, Var_integrale, Var_derivee) ; break ; //
			case 10 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, MODIFIED_GS, PREC_LEFT, 2, Var_integrale, Var_derivee) ; break ; // STALD = true, ... (id. ci-dessus) :  bon choix en general, mais pas avec vol_Sympl_dot !
			case 11 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, CLASSICAL_GS, PREC_LEFT, 2, Var_integrale, Var_derivee) ; break ; //
			case 12: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, CLASSICAL_GS, PREC_RIGHT, 2, Var_integrale, Var_derivee); break; // le + rapide : best for large N (>= 1000) ; break ; sinon, un peu instable avec VarVisc
			case 13: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, MODIFIED_GS, PREC_BOTH, 2, Var_integrale, Var_derivee); break; // STALD = true, verbose = true, rootfind = NULL
			case 14: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, CLASSICAL_GS, PREC_LEFT, 2, Var_integrale, Var_derivee); break; // le + rapide : best for large N (>= 1000) ; break ; sinon, un peu instable avec VarVisc
			case 15 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, MODIFIED_GS, PREC_RIGHT, 2, Var_integrale, Var_derivee) ; break ; // STALD = true, ... (id. ci-dessus)  :  bon choix en general, mais pas avec vol_Sympl_dot
			case 16 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, CLASSICAL_GS, PREC_RIGHT, 2, Var_integrale, Var_derivee) ; break ; //
			case 17: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, MODIFIED_GS, PREC_RIGHT, 2, Var_integrale, Var_derivee); break; // STALD = true, verbose = true, rootfind = NULL
			case 18: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, CLASSICAL_GS, PREC_BOTH, 2, Var_integrale, Var_derivee); break; // le + rapide : best for large N (>= 1000) ; break ; sinon, un peu instable avec VarVisc
			case 19: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, MODIFIED_GS, PREC_BOTH, 2, Var_integrale, Var_derivee); break; // STALD = true, verbose = true, rootfind = NULL
			case 20: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, CLASSICAL_GS, PREC_BOTH, 2, Var_integrale, Var_derivee); break; // le + rapide : best for large N (>= 1000) ; break ; sinon, un peu instable avec VarVisc
			case 21: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, MODIFIED_GS, PREC_LEFT, 2, Var_integrale, Var_derivee); break; // STALD = true, verbose = true, rootfind = NULL
			case 22: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, PCG, 1, PREC_BOTH, 2, Var_integrale, Var_derivee); break; //
			case 23: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPBCGS, 1, PREC_LEFT, 2, Var_integrale, Var_derivee); break; //
				case 24: j = arkode(fout, Y0, SegmentTimes, auxout, atol_, rtol, 2, Var_integrale, Var_derivee); break; //
			case 25: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPBCGS, 1, PREC_BOTH, 2, Var_integrale, Var_derivee); break;
			case 26 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPBCGS, 1, PREC_RIGHT, 2, Var_integrale, Var_derivee) ; break ; //
			case 27 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPTFQMR, 1, PREC_RIGHT, 2, Var_integrale, Var_derivee) ; break ; //
				case 28: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, BAND, 2, Var_integrale, Var_derivee, true, true, NULL, 0, neq / 100, neq / 100); break; // mu = ml = neq/100 : marche tres bien meme si neq < 100 !
				case 29: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, BAND, 2, Var_integrale, Var_derivee, true, true, NULL, 0, neq / 10, neq / 10); break; // mu = ml = neq/10 : OK
				case 30: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, BAND, 2, Var_integrale, Var_derivee, true, true, NULL, 0, neq / 3, neq / 3); break; // pas de rootfind, ; break ; mu = ml = neq/3 : OK
				case 31: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, DENSE, 2, Var_integrale, Var_derivee); break; // solver = cvode DENSE, STALD = true
				case 32: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, KLU, 2, Var_integrale, Var_derivee); break; //
			case 33 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPTFQMR, 1, PREC_LEFT, 2, Var_integrale, Var_derivee) ; break ; //
				case 34: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, BAND, 2, Var_integrale, Var_derivee, true, true, NULL, 0, neq / 300, neq / 300); break; //
			case 35: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPTFQMR, 1, PREC_BOTH, 2, Var_integrale, Var_derivee); break; //
			default : cout << endl << "!! Error !! solver # must be within the range [1 , 35] !!" << endl ; j = -1 ; exit(-1) ;
        }		
		if (j < 0) return (-1); // solver init error
		if (j > 0) break ; // simulation run-time error
        
    }
    // ********************* OUTPUT *********************************
	if(doTroubleshooting){
		std::cout<<"MEMORY LIBERATIONS"<<std::endl;
	}
    // MEMORY LIBERATIONS:
    delete [] y_dot;
	//for python:
	this->Q_outv = Y0.toCppVector();//att: now has several outputs
	
	this->Q_out_dotv = std::vector<double>(Y0.size(),0);
	for(int j = 1 ; j <= Y0.size() ; j ++) 
	{
		this->Q_out_dotv[j-1] = y_dot[j] ; 
	}
	this->C_STv = C_ST.toCppVector();
	this->vol_STv = vol_ST.toCppVector();
	this->Q_Grmaxv = Q_Grmax.toCppVector();
	this->Q_Exudmaxv = Q_Exudmax.toCppVector();
	this->Q_Rmmaxv = Q_Rmmax.toCppVector() ;
	this->Flv = Input.toCppVector() ;
	this->r_STv = r_ST.toCppVector() ;
	this->JW_STv = JW_ST.toCppVector() ;
	if(doTroubleshooting){
		std::cout<<"computeOrgGrowth"<<std::endl;
	}
	computeOrgGrowth(tf-t0);
	if(doTroubleshooting){
		std::cout<<"computeOrgGrowth	DONE"<<std::endl;
	}
	if(doTroubleshooting){std::cout.rdbuf(coutbuf);} //reset to standard output again
	return(1) ;
	
}

void OutputSettings()  {
	OutputTimes = Fortran_vector(nbv + 1) ;
	OutputTimes[1]=t0 ; if(nbv > 0) OutputTimes[2] = t1 ;
	for(i=2 ; i <= nbv ; i++) {
		if(LogScale) pas *= rs ;
		OutputTimes[i+1] = OutputTimes[i] + pas ;
	}
	OutputTimes[nbv + 1] = tf;
	if(Breakpoint_index[Breakpoint_index.size()]  != OutputTimes.size())  Breakpoint_index.append(Index_vector(1, OutputTimes.size())) ;
	//cout << "Output times :" << endl ; for(i=0 ; i < OutputTimes.size() ; i ++) cout << " " <<  OutputTimes[i+1] ; // if newly inserted extra output time (just after breakpoint) is too close, it may print as the same but actually is not
	nbv = OutputTimes.size() - 1 ; // at this point,  output times vector truc = OutputTimes ; and Breakpoint_index = {1, [...,] nbv+1}  with  [...,] = extra breakpoints indices (empty if none)
	
}



void PhloemFlux::aux(double t, double * y) {	// launch auxiliary calculations to store dynamics of temporary variables, , shared_ptr<PhloemFlux> PhloemObject
	f(t, y, y_dot) ;			// Update all variables from Y as returned by solver
	if(nbv > 1)
	{
		for (std::size_t it = 1 ; it <= Nt  ; it++){
			if(Q_ST[j]<0.)
			{
				std::cout<< "at t = " << t << " : Y0.size() = " << Y0.size()<<std::endl;
				std::cout<<"negative Q_ST value at j = "<<j<<", Q_ST[j] = "<< Q_ST[j] <<std::endl;
				assert(false);
			}
		}
	}
	
	for(int j = 1 ; j <= Y0.size() ; j ++) 
	{
		Y0[j] = y[j] ; 
		if((j <= Nt)&&(Y0[j] < 0.)){assert(false && "negative C_ST value");}
	}// update Y0
	if(doTroubleshooting){
		std::cout << "at t = " << t << " : Y0.size() = " << Y0.size();
		std::cout<<std::endl;
	}
}


void PhloemFlux::initializePM_(double dt, double TairK){	
	//all is in mmol/ml (=> M) and in d-1
	if(doTroubleshooting){cout<<"initializePM_1 "<<endl;}
	Adv_BioPhysics = false;
	vector<CPlantBox::Vector2i> segmentsPlant = plant->segments;
	Nc = segmentsPlant.size(); Nt = plant->nodes.size();
	assert((Nc == (Nt -1))&&"Wrong seg and node number");
	atol_= Fortran_vector(Nt*neq_coef, atol_double);//1e-017
	r_ST_ref = Fortran_vector(Nc)	; 
	I_Upflow = I_Downflow = std::vector<int>(Nc +1) ;
	vector<int> orgTypes = plant->organTypes;//per seg
	vector<int> subTypes = plant->subTypes;
	vector<double> Radii = plant->radii;
	vector<double> Lengthvec = plant->segLength();
	rtol =  Fortran_vector(1, rtol_double);
	double a_seg, l;double mu =0.;
	int ot, st;
	Ag =Fortran_vector(Nt, 0.);
	a_STv.resize(Nc , 0.)  ;//for postprocessing
	len_leaf =Fortran_vector(Nt, 0.);
	vol_ParApo =Fortran_vector(Nt, 0.);
	Q_Grmax =Fortran_vector(Nt, 0.);
	Q_Exudmax=Fortran_vector(Nt, 0.);
	Q_Rmmax =Fortran_vector(Nt, 0.);
	vol_ST=Fortran_vector(Nt, 0.);
	vol_Seg=Fortran_vector(Nt, 0.);//for postprocessing	
	exud_k=Fortran_vector(Nt, 0.);
	krm2=Fortran_vector(Nt, 0.);
	Csoil_node.resize(Nt,0.);
	deltaSucOrgNode_ = waterLimitedGrowth(dt);//water limited growth
	if(doTroubleshooting){cout<<"initializePM_new "<<Nc<<" "<<Nt<<endl;}
	int nodeID;
	double StructSucrose;//double deltaStructSucrose;
	double cmH2O_to_hPa = 0.980638	;//cm to hPa
	double krm1;
	k_mucil_.resize(Nt , 0.);
	if(psiXyl4Phloem.size() == Nt){Psi_Xyl = Fortran_vector(psiXyl4Phloem,cmH2O_to_hPa); //computed by xylem object
	}else{Psi_Xyl = Fortran_vector(Nt, 0.);}
	
	if(! update_viscosity_)
	{	

		double TdC = TairK - 273.15;
		//in g/L or mg/cm3
		dEauPure = (999.83952 + TdC * (16.952577 + TdC * (- 0.0079905127 + TdC * (- 0.000046241757 + TdC * (0.00000010584601 + TdC * (- 0.00000000028103006)))))) / (1 + 0.016887236 * TdC); 
		double siPhi = (30 - TdC) / (91 + TdC) ; // T_old = TairK_phloem ;//  R.Gilli 1997, after Mathlouthi and Genotelle 1995 - valid for any T :
				
		double C = 0.5 ; // (mmol / ml solution)
		//  R.Gilli 1997, after Mathlouthi and Genotelle 1995 - valid for any T :
		//342.3 g/mol or mg/mmol
		double PartMolalVol_ =0;// 0.2155;
		double d = C * 342.3 + (1 - C * PartMolalVol_) * dEauPure ;//in mg/cm3
		double siEnne = (100 * 342.30 * C) / d ; // actually this is sc = sucrose content (g.suc. % g.solution) ; 342.30 = molar mass of sacch.
		siEnne /= 1900 - (18 * siEnne) ;
		//mPa s
		mu =  pow(10, ((22.46 * siEnne) - 0.114 + (siPhi * (1.1 + 43.1 * pow(siEnne, 1.25) )))) ; // peut atteindre des valeurs > 1.e200 !! (sans signification evidemment -- le sucre doit precipiter bien avant !!)
		mu = mu /(24*60*60)/100/1000; //mPa s to hPa d, 1.11837e-10 hPa d for pure water at 293.15K
		
	}
	
	for ( std::size_t k=1 ;k <= Nc;k++ ) // seg id. == node.y id
	{
			//the PiafMunch Fortran_vector are numbered from 1 (and not from 0)
			//so we need to account for that...
			// element 0 contains the vector size i think
			ot = orgTypes[k-1]; 
			st = plant->st2newst[std::make_tuple(ot,subTypes[k-1])];
			a_seg = Radii[k-1];
			l = Lengthvec[k-1];
			I_Upflow[k] = segmentsPlant[k-1].x +1;// +1 to go from C++ numeration to pseudo-Fortran numeration
			I_Downflow[k] = segmentsPlant[k-1].y +1;
			nodeID = I_Downflow[k];
			int nodeID_cpb = segmentsPlant[k-1].y;
			len_leaf[nodeID] = l;assert((len_leaf[nodeID]>0)&&"len_seg[nodeID] <=0");
			
			double Across_ST =  Across_st_f(st,ot);assert((Across_ST>0)&&"Across_ST <=0");
			a_STv[k-1] = Across_ST;
			
			vol_ST[nodeID] = Across_ST * l;
			
		//general
			r_ST_ref[k] = 1/kx_st_f(  st,ot)*l;
			if(! update_viscosity_)
			{
				r_ST_ref[k] = mu*r_ST_ref[k];
			}assert((r_ST_ref[k]>0)&&"r_ST_ref[k] <=0");
			
		//Fl
			Ag[nodeID] = Ag4Phloem[segmentsPlant[k-1].y] ;//can be negative at night
		//Fu
			//Exudation
			exud_k[nodeID] =kr_st_f(st,ot, k-1);
			Q_Exudmax[nodeID ]=0;
			if((ot==2)&&(this->plant->seg2cell[k-1] >=0)){//root segment belowground
				Q_Exudmax[nodeID ] =   2 * M_PI * a_seg * l*exud_k[nodeID] ;
				if (exud_k[nodeID] > 0.)
				{
					// k_mucil is a std::vecto not a fotran vector
					k_mucil_[segmentsPlant[k-1].y ] = k_mucil;
				}
				if(Csoil_seg.size() > 0.)
				{
					// Csoil_node is a std::vecto not a fotran vector
					Csoil_node.at(segmentsPlant[k-1].y ) = Csoil_seg.at(k -1);					
				}else{
					Csoil_node.at(segmentsPlant[k-1].y ) = CsoilDefault;
				}
			}
			if(doTroubleshooting){std::cout<<"QexudMax "<<Q_Exudmax[nodeID ]<<std::endl;}
			
			
			//Rm, maintenance respiration 
			krm1 =  krm1_f(st,ot);
			krm2[nodeID] =  krm2_f(st,ot);
			
			vol_Seg[nodeID] = plant->segVol[k-1];
			double l_blade = plant->bladeLength[k-1];//area in leaf blade
			vol_ParApo[nodeID] = -1;
			
			if( l_blade>0){
				if(sameVolume_meso_st){
						vol_ParApo[nodeID] = vol_ST[nodeID];
					}else
					{
						if(sameVolume_meso_seg){
							vol_ParApo[nodeID] = vol_Seg[nodeID];
						}else{
							vol_ParApo[nodeID] = l_blade * surfMeso;
						}
					}
			}		
			StructSucrose = rhoSucrose_f(st,ot) * vol_Seg[nodeID]; 
			if(doTroubleshooting){std::cout<<"LeafShape "<<(k-1)<<" "<<l_blade<<" "<<vol_ParApo[nodeID]<<" "<<surfMeso<<std::endl;}
			
			Q_Rmmax[nodeID] = krm1 * StructSucrose;
			if(doTroubleshooting){
				std::cout<<"forRmmax "<<nodeID<<" "<< vol_Seg[nodeID]<<" "<<krm1<<" "<<krm2[nodeID]<<" "<<rhoSucrose_f(st,ot)<<" "<<Q_Rmmax[nodeID]<<std::endl;
			}
						
			Q_Grmax[nodeID] = deltaSucOrgNode_.at(nodeID_cpb).at(-1)/Gr_Y/dt;
			if(doTroubleshooting){
				std::cout<<"forGrmax"<<nodeID<<" "<< deltaSucOrgNode_.at(k).at(-1)<<" "<<Q_Grmax[nodeID]<<std::endl;
			}
		
		//Test
			if(exud_k[nodeID]<0.){
				std::cout<<"exud_k[nodeID]: loop n#"<<k<<", node "<<nodeID<<" "<<exud_k[nodeID]<<" "<<(exud_k[nodeID]<0.);
				std::cout<<" "<<" "<<(exud_k[nodeID]==0.)<<" "<<" "<<(exud_k[nodeID]>0.)<<std::endl;
				assert(false);
			}
			if(krm2[nodeID]<0.){
				std::cout<<"krm2: loop n#"<<k<<", node "<<nodeID<<" "<<krm2[nodeID]<<std::endl;
				assert(false);
			}
			if(Q_Grmax[nodeID ]<0.){
				std::cout<<"gr: loop n#"<<k<<", node "<<nodeID<<" "<<Q_Grmax[nodeID]<<" "<< deltaSucOrgNode_.at(k).at(-1)<<std::endl;
				assert(false);
			}
			if(Q_Exudmax[nodeID ]<0.){
				std::cout<<"exud: loop n#"<<k<<", node "<<nodeID<<" "<<Q_Exudmax[nodeID]<<" "<< l<<" "<<Radii[k-1]<<std::endl;
				assert(false);
			}
			if(Q_Rmmax[nodeID ]< 0.){
				std::cout<<"rm: loop n#"<<k<<", node "<<nodeID<<" "<<Q_Rmmax[nodeID]<<" "<< krm1<<" "<<StructSucrose<<std::endl;
				assert(false);
			}
			//
		}
		
	//for seed node:
	a_STv.at(0) = a_STv.at(1);
	len_leaf[1] = len_leaf[2];
	vol_ParApo[1] = vol_ParApo[2];
	vol_ST[1] = vol_ST[2];
	Q_Rmmax[1] =0;
	Q_Grmax[1] =  deltaSucOrgNode_.at(0).at(-1)/Gr_Y/dt;
	vol_Seg[1] = vol_Seg[2];
	krm2[1]=0;
	exud_k[1]=0;Q_Exudmax[1]=0;
	
	//for post processing in python
	this->Q_Grmaxv = Q_Grmax.toCppVector();
	this->r_ST_refv = r_ST_ref.toCppVector();
	this->Agv = Ag.toCppVector();
	this->vol_STv = vol_ST.toCppVector();
	this->vol_Mesov = vol_ParApo.toCppVector();
}



void rootfind(double t, double* y, double* g) {return ; } // required but not used by cvode_xxxxx()

void PhloemFlux::computeOrgGrowth(double t){
	
	int nNode = plant->getNumberOfNodes();
	auto orgs = plant->getOrgans(-1, true); //get all organs, even small ones
	int nOrg = orgs.size();
	std::map<int, double> cWGrRoot; //set cWGr in this function instead of in Mappedorganism.h : cWGr is then reset to empy every tim efunction is called
	std::map<int, double> cWGrStem; // + no need for mapped organism to keep cWGr in memory. just has to be in growth function
	std::map<int, double> cWGrLeaf;
	
	double delta_volOrg;
	double delta_volOrgmax;
	int orgID2 = 0;//needed as max(orgID) can be > plant->getOrgans(-1, true).size()
	for(auto org: orgs){
		int orgID = org->getId();
		std::vector<int> nodeIds_;
		int ot = org->organType();
		int stold = org->getParameter("subType");
		int st = plant->st2newst[std::make_tuple(ot,stold)];
		int nNodes = org->getNumberOfNodes();
		
		if(ot == 2){
			nodeIds_.push_back(org->getNodeId(nNodes-1));			
			}else{nodeIds_ = org->getNodeIds();}
		
		delta_volOrg = 0;
		delta_volOrgmax = 0;
		std::vector<int> allnodeids = org->getNodeIds();
		// TODO: Y not jsut ocmpute again the growing nodes ID?
		for(std::size_t k = 0; k < nodeIds_.size(); k++)//take  first node as might carry small laterals
		{
			double deltaSucmax_1 = 0;
			if(deltaSucOrgNode_.at(nodeIds_.at(k)).count(orgID) !=0){
				deltaSucmax_1 = deltaSucOrgNode_.at(nodeIds_.at(k)).at(orgID);//maxsuc needed for Gr
				if(doTroubleshooting){std::cout<<"deltaSucmax_1 "<<orgID<<" "<<k<<" "<<deltaSucmax_1<<std::endl;}
			}
			double deltaSucmax_ = deltaSucmax_1 /Gr_Y;//max suc needed for Gtot
			double deltaSuc = std::min(Q_Growthtot[nodeIds_.at(k) +1] - Q_GrowthtotBU[nodeIds_.at(k) +1], deltaSucmax_);
			if((deltaSucmax_ > Q_Grmax[nodeIds_.at(k) +1]))
			{
				std::cout<<"error Q_Grmax "<<k<<" "<<(nodeIds_.at(k)+1)<<" "<<deltaSucmax_<<" "<<Q_Grmax[nodeIds_.at(k) +1]<<std::endl;
				assert(false);
			}
			if((deltaSuc < 0.)&&(deltaSuc > -1e-5)){deltaSuc=0.;}
			if((deltaSuc < 0.)){
				assert((deltaSuc >= 0.) &&"negative deltaSuc");
			}
			delta_volOrg += deltaSuc*Gr_Y/rhoSucrose_f(st,ot);
			delta_volOrgmax += deltaSucmax_*Gr_Y/rhoSucrose_f(st,ot);
			
			
			Q_GrowthtotBU[nodeIds_.at(k) +1] += deltaSuc;//sucrose is spent
			Q_GrmaxBU[nodeIds_.at(k) +1] += deltaSucmax_;//sucrose is spent
		}
		
		double vol = org->orgVolume(-1, false);
		double l_ =  org->getLength(false);
		
		double vol_check = org->orgVolume(l_, false);
		double l_check = org->orgVolume2Length(vol);
		
		if((abs(l_ - l_check)>1e-5) || (abs(vol - vol_check)>1e-5)){
			std::cout<<"(abs(l_ - l_check)<1e-5) || (abs(vol - vol_check)<1e-5)"<<std::endl;
			std::cout<<ot<<" "<<(abs(l_ - l_check)<1e-5)<<" "<< (abs(vol - vol_check)<1e-5)<<std::endl;
			std::cout<<l_<<" "<<l_check<<" "<<vol<<" "<<vol_check<<std::endl;
			assert(false);
		}
		double newl = l_;
		double newl_max = l_;
		if(delta_volOrg > 0)//if vol + delta_volOrg == 0 can lead to nan
		{
			newl =org->orgVolume2Length(vol + delta_volOrg);
		}
		if(delta_volOrgmax > 0)
		{
			newl_max = org->orgVolume2Length(vol + delta_volOrgmax);
		}
		double orgGr = newl - l_;
		double delta_lmax = newl_max - l_;
		
		if((orgGr < 0.)&&(orgGr > -1e-10)){orgGr=0.;}
			assert((orgGr >= 0.) &&"negative orgGr");
		if((delta_lmax < 0.)&&(delta_lmax > -1e-10)){delta_lmax=0.;}
			assert((delta_lmax >= 0.) &&"negative delta_lmax");	
		
		if(ot == 2){//each organ type has it s own growth function (and thus CW_Gr map)
			cWGrRoot.insert(std::pair<int, double>(orgID, orgGr));
		}
		if(ot == 3){
			cWGrStem.insert(std::pair<int, double>(orgID, orgGr));
		}
		if(ot == 4){
			cWGrLeaf.insert(std::pair<int, double>(orgID,orgGr));
		}
		if((BackUpMaxGrowth[org->getId()] - newl)< -1e-10)
		{
			std::cout <<ot<<" "<<st<<" "<<orgID<<" "<<orgID2<<" "<<newl<<" "<<BackUpMaxGrowth[org->getId()] <<" ";
			std::cout<<(BackUpMaxGrowth[org->getId()]- newl)<<std::endl;
			throw std::runtime_error("PhloemFlux::computeOrgGrowth: new target length of organ too high");
		}
		orgID2 ++;
		
	}
	for (auto orp : plant->getOrganRandomParameter(2)) { // each maps is copied for each sub type; or (todo), we could pass a pointer, and keep maps in this class
		if(orp!= NULL) {orp->f_gf->CW_Gr = cWGrRoot;}
	}
	
	for (auto orp : plant->getOrganRandomParameter(3)) {
		if(orp!= NULL) {orp->f_gf->CW_Gr = cWGrStem;}
	}
	for (auto orp : plant->getOrganRandomParameter(4)) {
		if(orp!= NULL) {orp->f_gf->CW_Gr = cWGrLeaf;}
	}
}

double PhloemFlux::adaptDt(std::shared_ptr<CPlantBox::Organ> org, double dt)
{
	double age = org->getAge();
	int ot = org->organType(); 
	double orgLT = org->getParameter("rlt");
	
	if (age+dt>orgLT) { // root life time
		dt=orgLT-age; // remaining life span
	}
	
	//no probabilistic branching models and no other scaling via getRootRandomParameter()->f_se->getValue(nodes.back(), shared_from_this());
	if(ot == CPlantBox::Organism::ot_stem){
		
		double delayNGStart = org->getParameter("delayNGStart");
		double delayNGEnd = org->getParameter("delayNGEnd");
		
		if((age+dt) > delayNGStart){//simulation ends after start of growth pause
			if((age+dt)  < delayNGEnd){dt = 0;//during growth pause
			}else{
				double beforePause = 0.; double afterPause =0.;
				if(age < delayNGStart){
					beforePause = std::min(dt,delayNGStart - age); //reaches start of growth pause and then stopes
				}
				if((age+dt)  > delayNGEnd){
					afterPause = std::min(dt, std::max(age +dt  - delayNGEnd, 0.)); //part of the simulation after end of pause
				}
				dt = beforePause + afterPause;//part of growth pause during simulation				
			}
		}
	}
	if(dt < 0){
		throw std::runtime_error("PhloemFlux::waterLimitedGrowth: dt <0");
	}
	return dt;
	
}

void PhloemFlux::assertUsedCReserves(std::shared_ptr<CPlantBox::Organ> org)
{
	if((!((org->getOrganRandomParameter()->f_gf->CW_Gr.empty()) || 
					(org->getOrganRandomParameter()->f_gf->CW_Gr.count(org->getId()) ==0) ||
					(org->getOrganRandomParameter()->f_gf->CW_Gr.find(org->getId())->second<0.)))&&
					org->isActive()&&useCWGr)
	{
		std::cout<<org->getId()<<" "<<org->getOrganRandomParameter()->f_gf->CW_Gr.find(org->getId())->second<<std::endl;
		std::cout<<org->calcLength(1)<<" "<< org->organType() <<" "<<org->getAge()<<std::endl;
		throw std::runtime_error("PhloemFlux::waterLimitedGrowth: sucrose for growth has not been used at last time step");
	}
}

double PhloemFlux::getMaxVolumicGrowth(std::shared_ptr<CPlantBox::Organ> org, double t)
{
	int ot = org->organType(); 
	int st = plant->st2newst[std::make_tuple(ot,org->getParameter("subType"))];
	double rmax = Rmax_st_f(st,ot);
	double Linit = org->getLength(false);//theoretical length 
			
	if(org->getParameter("gf") != 3)
	{
		std::cout<<"PhloemFlux::waterLimitedGrowth: organ(s) do(es) not use carbon-limited growth"<<std::endl;
		do_gf_warning = false;
	}
	
	auto f_gf =  plant->createGrowthFunction(1);
	double age_ = f_gf->getAge(Linit, rmax, org->getParameter("k"), org->shared_from_this());			
	assertUsedCReserves(org);			
	double dt = adaptDt(org, t);				
	//	params to compute growth
	double LinitTemp = f_gf->getLength(age_ , rmax, org->getParameter("k"), org->shared_from_this());
	double targetlength = f_gf->getLength(age_ + dt , rmax, org->getParameter("k"), org->shared_from_this());
	double e = std::max(0.,targetlength-LinitTemp);// unimpeded elongation in time step dt
	
	if((e + org->getLength(false))> org->getParameter("k") + 1e-10){
		std::cout<<"PhloemFlux::getMaxLengthGrowth: target length too high "<<e<<" "<<dt<<" "<<Linit;
		std::cout<<" "<<org->getParameter("k")<<" "<<org->getId()<<std::endl;
		assert(false);
	}
	BackUpMaxGrowth[org->getId()] = Linit + e; // to compare with final growth length in @see PhloemFlux::computeOrgGrowth
	double deltavol = std::max(0.,org->orgVolume(Linit + e, false) - org->orgVolume(Linit, false));//volume from theoretical length
	return deltavol;
}

/*
	//roots grow at their tip
	//leaves grow in their groth zones
	//stems grow at each phytomere nodes
*/

std::vector<int> PhloemFlux::getGrowingNodes(std::shared_ptr<CPlantBox::Organ> org)
{
	int ot = org->organType(); 
	int nNodes = org->getNumberOfNodes();
	std::vector<int> nodeLocalIds;

	if ((nNodes==1)||(ot == 2)) { //organ not represented because below dx limit or is root
		nodeLocalIds.push_back(nNodes - 1);
	}else{
		if(ot == 4)//only for nodes in the leaf growth zone
		{
			for(int k = 0; k < nNodes; k++)
			{
				nodeLocalIds.push_back(k);
				if(org->getLength(k) > this->leafGrowthZone)
				{
					break;
				}
			}
		}
		if(ot == 3)
		{
			int nn1 = 0;
			if((StemGrowthPerPhytomer)&&(org->getNumberOfChildren() > 0))
			{
				// only for nodes in the growing phytomers. Phytomer ends each time we have a lateral.
				auto stemParams = std::static_pointer_cast<CPlantBox::Stem>(org)->param();
				bool foundPhytoIdx = false;
				int PhytoIdx = 0; // to remember PhytoIdx outside of loop
				for(; ((PhytoIdx < stemParams->ln.size())&(!foundPhytoIdx));PhytoIdx++)
				{
					double maxPhytoLen = stemParams->ln.at(PhytoIdx);
					double currentPhytoLen = org->getLength(org->getChild(PhytoIdx+1)->parentNI) - org->getLength(org->getChild(std::max(0,PhytoIdx))->parentNI);
					foundPhytoIdx = (maxPhytoLen - currentPhytoLen > 1e-10);//first still growing phytomere
					if(maxPhytoLen - currentPhytoLen < -1e-10){throw std::runtime_error("maxPhytoLen - currentPhytoLen < -1e10;");}
					
				}
				assert((nn1 < nNodes)&&"nn1 > nNodes");
			}
			for (int i = nn1; i < nNodes; ++i) {
				nodeLocalIds.push_back(i);
			}
		}
	}
	return nodeLocalIds;
}

void PhloemFlux::computeFpsi()
{
	Fpsi = std::vector<double>(Nt,1.); 
	psi_p_symplasm = std::vector<double>(Nt,1.); 
	for(std::size_t nodeId = 0; nodeId < psiXyl.size(); nodeId++)
	{
		psi_p_symplasm.at(nodeId) =  psiXyl.at(nodeId) - psi_osmo_proto;
		if(-psi_osmo_proto - psiMin>1e-5){
			Fpsi.at(nodeId) = std::max((std::max(psi_p_symplasm.at(nodeId), psiMin) - psiMin)/(-psi_osmo_proto - psiMin),0.);
		}else{Fpsi.at(nodeId) = 0.;}
		
		if((Fpsi.at(nodeId) > (1+1e-6))||(Fpsi.at(nodeId) < -1e-6))
		{
			throw std::runtime_error("Fpsi < 0 or Fpsi > 1");
		}
	}
}

void PhloemFlux::computeFlen(std::shared_ptr<CPlantBox::Organ> org, std::vector<int> growingNodesId)
{
	double Flen_tot = 0.;
	if (growingNodesId.size() == 1)
	{
		Flen.at(org->getNodeId(growingNodesId.at(0))) = 1.;
	}else{
		double L_n0 = org->getLength(growingNodesId.at(0)); 
		double L_growth_length = org->getLength(growingNodesId.at(growingNodesId.size() - 1)) - L_n0;
		for(std::size_t k = 1; k < growingNodesId.size(); k++)//int nodeId : growingNodesId)
		{
			int nodeGlobalId = org->getNodeId(growingNodesId.at(k));
			auto nodei = plant->nodes.at(org->getNodeId(growingNodesId.at(k - 1)));
			auto nodej = plant->nodes.at(nodeGlobalId);
			double Lseg = nodej.minus(nodei).length();
			Flen.at(nodeGlobalId) = (Lseg/L_growth_length);
			Flen_tot += Flen.at(nodeGlobalId);
		}
		if((std::abs(Flen_tot - 1.)>1e-10)||(std::abs(Flen_tot - 1.)<-1e-10))
		{
			throw std::runtime_error("PhloemFlux::computeFlen: (Flen_tot != 1.)&&wrong tot Flen");
		}
	}
}

	/*water limited deltaSuc per node*/
std::vector<std::map<int,double>> PhloemFlux::waterLimitedGrowth(double t)
{
	std::map<int,double> initMap = { { -1,0. } };//-1 : total need per node
	std::vector<std::map<int,double>> deltaSucOrgNode(Nt, initMap);//, 0.);
	bool allOrgs = true;
	auto orgs = plant->getOrgans(-1,allOrgs);//also org with length - Epsilon == 0	
	BackUpMaxGrowth = {};//std::vector<double>(orgs.size(), 0.); //for checking in @see PhloemFlux::computeOrgGrowth
	Flen = std::vector<double>(Nt,0.); 
	computeFpsi();
	for(auto org: orgs)
	{
		assert(org != nullptr);
		double dt = t; //copy
		if(((org->getAge()+dt)>0)&&(org->isAlive())&&(org->isActive()))
		{ //organ alive and active at this time step
			int ot = org->organType(); 
			int stold = org->getParameter("subType");
			int st = plant->st2newst[std::make_tuple(ot,stold)];
			std::vector<int> localGrowingNodesId = getGrowingNodes(org); // localIds			
			double deltaVol = getMaxVolumicGrowth(org, dt); // cm3				
			computeFlen(org, localGrowingNodesId);
			double deltaSucGrowth_tot = 0.;
			double rhoSucrose_double = rhoSucrose_f(st,ot);
			for(int localNodeId : localGrowingNodesId)
			{
				int nodeId = org->getNodeId(localNodeId);
				double deltaSucGrowth_per_node = deltaVol * Flen.at(nodeId) * Fpsi.at(nodeId) * rhoSucrose_double;
				
				if((deltaSucGrowth_per_node<0.)||(deltaSucGrowth_per_node != deltaSucGrowth_per_node)){
					// could be error of pressision (if l = Lmax)
					// or that, because of nodal growth and dxMin, org->getEpsilon() <0
					if(((deltaSucGrowth_per_node>-1e-5)||(deltaVol>0.))
							&&(deltaSucGrowth_per_node == deltaSucGrowth_per_node)){
						deltaSucGrowth_per_node=0.;	//within margin of error
					}else{
						throw std::runtime_error("(deltaSucGrowth_per_node<0.)||(deltaSucGrowth_per_node != deltaSucGrowth_per_node)");
					}					
				}				
				assert(nodeId >= 0 && nodeId < (int)deltaSucOrgNode.size());
				deltaSucOrgNode.at(nodeId)[org->getId()] = deltaSucGrowth_per_node;
				//deltaSucOrgNode.at(nodeId).insert(std::pair<int, double>(org->getId(), deltaSucGrowth_per_node));
				deltaSucOrgNode.at(nodeId).at(-1) += deltaSucGrowth_per_node;								
				deltaSucGrowth_tot += deltaSucGrowth_per_node;
			}
			assert(((deltaSucGrowth_tot * rhoSucrose_double - deltaVol)<1e-10)&&"deltaSucGrowth_tot too high");//deltaSucGrowth_tot * rhoSucrose_double <= deltaVol
		}else{
			BackUpMaxGrowth[org->getId()] = org->getLength(false);
		}
	}	
	return deltaSucOrgNode;
}


/**
 *  Sets the radial conductivity in [1 day-1]
 * in case of organ_type specific kr
 * @param values 		kr per organ pr/and organ type) or/and per age [cm-1]
 * @param age 			ages if kr per age
 * @param kr_length_ 	exchange zone in root, where kr > 0 [cm from root tip], default = -1.0, i.e., no kr_length
 */
//type/subtype dependent
void PhloemFlux::setKr_st(std::vector<std::vector<double>> values, double kr_length_, bool verbose) {
    kr_st = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			kr_st_f = std::bind(&PhloemFlux::kr_st_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
            if (verbose)
            {
                std::cout << "Kr_st is constant " << values[0][0] << " 1 day-1 \n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			kr_st_f = std::bind(&PhloemFlux::kr_st_perOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
			if (verbose)
            {
                std::cout << "Kr_st is constant per organ type, organ type 2 (root) = " << values[0][0] << " 1 day-1 \n";
            }
		} else {
			if (verbose)
            {
                std::cout << "Exchange zone in roots: kr_st > 0 until "<< kr_length_<<"cm from root tip "<<(kr_length_ > 0)<<" "<<(kr_length_ > 0.)<<std::endl;
            }
			if(kr_length_ > 0.){
				if (verbose)
            {
                std::cout << "Exchange zone in roots: kr > 0 until "<< kr_length_<<"cm from root tip"<<std::endl;
                }
				plant->kr_length = kr_length_; //in MappedPlant. define distance to root tipe where kr > 0 as cannot compute distance from age in case of carbon-limited growth
				plant->calcExchangeZoneCoefs();	
				kr_st_f  = std::bind(&PhloemFlux::kr_st_RootExchangeZonePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
			}else{
				kr_st_f  = std::bind(&PhloemFlux::kr_st_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
			}
			if (verbose)
            {
                std::cout << "Kr_st is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << " 1 day-1 \n";
            }
		}
	}
    
}	
// type/subtype dependent
void PhloemFlux::setKx_st(std::vector<std::vector<double>> values, bool verbose) {
    kx_st = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			kx_st_f = std::bind(&PhloemFlux::kx_st_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Kx_st is constant " << values[0][0] << " cm3 day-1 \n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			kx_st_f = std::bind(&PhloemFlux::kx_st_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Kx_st is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm3 day-1 \n";
            }
		} else {
			kx_st_f  = std::bind(&PhloemFlux::kx_st_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Kx_st is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << " cm3 day-1 \n";
            }
		}
	}
}	


//type/subtype dependent
void PhloemFlux::setAcross_st(std::vector<std::vector<double>> values, bool verbose) {
    Across_st = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			Across_st_f = std::bind(&PhloemFlux::Across_st_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Across_st is constant " << values[0][0] << " cm2\n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			Across_st_f = std::bind(&PhloemFlux::Across_st_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Across_st is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm2 \n";
            }
		} else {
			Across_st_f  = std::bind(&PhloemFlux::Across_st_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Across_st is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << " cm2 \n";
            }
		}
	}
}	


//type/subtype dependent ==> to compute Rhat Fhat
void PhloemFlux::setPerimeter_st(std::vector<std::vector<double>> values, bool verbose) {
    Perimeter_st = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			Perimeter_st_f = std::bind(&PhloemFlux::Perimeter_st_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Perimeter_st is constant " << values[0][0] << " cm\n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			Perimeter_st_f = std::bind(&PhloemFlux::Perimeter_st_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Perimeter_st is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm\n";
            }
		} else {
			Perimeter_st_f  = std::bind(&PhloemFlux::Perimeter_st_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Perimeter_st is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << " cm\n";
            }
		}
	}
}	

// type/subtype dependent
void PhloemFlux::setRmax_st(std::vector<std::vector<double>> values, bool verbose) {
    Rmax_st = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			Rmax_st_f = std::bind(&PhloemFlux::Rmax_st_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Rmax_st is constant " << values[0][0] << " cm day-1 \n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			Rmax_st_f = std::bind(&PhloemFlux::Rmax_st_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Rmax_st is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm day-1 \n";
            }
		} else {
			Rmax_st_f  = std::bind(&PhloemFlux::Rmax_st_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Rmax_st is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << " cm day-1 \n";
            }
		}
	}
}	


//either age or type/subtype dependent
void PhloemFlux::setRhoSucrose(std::vector<std::vector<double>> values, bool verbose) {
    rhoSucrose = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			rhoSucrose_f = std::bind(&PhloemFlux::rhoSucrose_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "rhoSucrose is constant " << values[0][0] << " mmol cm-3\n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			rhoSucrose_f = std::bind(&PhloemFlux::rhoSucrose_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "rhoSucrose is constant per organ type, organ type 2 (root) = " << values[0][0] << " mmol cm-3\n";}
		} else {
			rhoSucrose_f  = std::bind(&PhloemFlux::rhoSucrose_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "rhoSucrose is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << " mmol cm-3\n";}
		}
	}
}	


//either age or type/subtype dependent
void PhloemFlux::setKrm1(std::vector<std::vector<double>> values, bool verbose) {
    krm1v = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			krm1_f = std::bind(&PhloemFlux::krm1_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "krm1 is constant " << values[0][0] << " -\n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			krm1_f = std::bind(&PhloemFlux::krm1_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "krm1 is constant per organ type, organ type 2 (root) = " << values[0][0] << " -\n";
            }
		} else {
			krm1_f  = std::bind(&PhloemFlux::krm1_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "krm1 is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << "-\n";
            }
		}
	}
}	

//either age or type/subtype dependent
void PhloemFlux::setKrm2(std::vector<std::vector<double>> values, bool verbose) {
    krm2v = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			krm2_f = std::bind(&PhloemFlux::krm2_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "krm2 is constant " << values[0][0] << " -\n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			krm2_f = std::bind(&PhloemFlux::krm2_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			std::cout << "krm2 is constant per organ type, organ type 2 (root) = " << values[0][0] << " -\n";
		} else {
			krm2_f  = std::bind(&PhloemFlux::krm2_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "krm2 is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << "-\n";
            }
		}
	}
}	
