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
#ifndef decoupled
#include <PiafMunch/runPM.h>
#else
#include <runPM.h>
#endif

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
extern vector<int> I_Upflow, I_Downflow ; // I_Upflow(resp.I_Downflow)[jf=1..Nc] = id# du noeud amont (resp. aval) : jf = JF(i,i2) > 0 si I_Upflow[(abs(jf)]==i, i.e. si I_Downflow[(abs(jf)]==i2
extern Fortran_vector kML						; // kinetic parameter / Michaelis - phloem loading					(mmol / ml)
extern Fortran_vector vMU					; // kinetic parameter / phloem unloading								(mmol /h)
//extern Fortran_vector	isTip ,Ag,  Lmax_org, Rmax_org, krm1, krm2 , StructC, exud_k; 
extern Fortran_vector radius_ST ,  Length, vol_ST_seg;// Q_Grmax, Q_Rmmax, Q_Exudmax; 

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

std::shared_ptr<PhloemFlux> phloem_; 
vector<double> extra_output_times, breakpoint_times ; // may contain already existing output time points
extern int is, solver ; // (see below) solver config.# -- default = 1 = cvode (SPILS: SPFGMR, MODIFIED_GS, PREC_NONE... Try diff. value if calc. fails or too slow !
void auxout(double t, double * y){phloem_->aux(t, y);} ;	// launch auxiliary calculations to store dynamics of temporary variables
void fout(double t, double *y, double *y_dot){phloem_->f(t, y, y_dot);} ; // the function to be processed by the solver  (implemented in 'solve.cpp')

// Full list of printable/savable variables -- to be updated in case of hard-updating the model (number and nature of local compartments within each archit. element) :
int NumAllNodeVariablesNames = 62 ; // 62 node variables :
string AllNodeVariablesNamesQList[] = {
	"Absorb (ml / h)", "Starch (mmol eq. Glu)", "StarchSyn (mmol eq. Glu / h)", "C_ST (mmol / ml)", "C_Sympl (mmol / ml)", "JS_PhlMb (mmol / h)", "JW_Trsv (ml / h)", "JW_PhlMb (ml / h)", "JW_ParMb (ml / h)", "Psi_ST (MPa)",
	"P_Xyl (MPa)", "Q_Sympl (mmol)", "Q_ST (mmol)", "Resp_Maint (mmol / h)", "C_SymplUpflow (mmol / ml)", "JS_Sympl (mmol / h)", "JW_Apo (ml / h)", "JW_Sympl (ml / h)", "P_ST (MPa)", "P_Sympl (MPa)",
	"Psi_PhlApo (MPa)", "Psi_ParApo (MPa)", "Psi_Sympl (MPa)", "vol_Sympl (ml)", "vol_Sympl_dot (ml / h)", "P_PhlApo (MPa)", "P_ParApo (MPa)", "Q_PhlApo (mmol)", "Q_ParApo (mmol)",
	"C_PhlApo (mmol / ml)", "C_ParApo (mmol / ml)", "Q_PhlApo_dot (mmol / h)", "Q_ParApo_dot (mmol / h)", "JS_Apo (mmol / h)", "JS_ParMb (mmol / h)", "C_ApoUpflow (mmol / ml)",
	"P_ST_dot (MPa / h)", "P_Sympl_dot (MPa / h)", "Starch_dot (mmol eq. Glu / h)", "Q_ST_dot (mmol / h)", "Q_Sympl_dot (mmol / h)", "Input (mmol / h)", "Transpirat (ml / h)", "PsiSoil (MPa)",
// Tracer-specific add-ins :
	"TracerStarch (MBq)", "TracerStarchSyn (MBq / h)", "TracerC_ST (MBq / ml)", "TracerC_Sympl (MBq / ml)", "TracerJS_PhlMb (MBq / h)", "TracerQ_Sympl (MBq)", "TracerQ_ST (MBq)", "TracerRespMaint (MBq / h)",
	"TracerInput (MBq / h)", "TracerJS_Sympl (MBq / h)", "TracerJS_ParMb (MBq / h)", "TracerQ_PhlApo (MBq)", "TracerQ_ParApo (MBq)", "TracerC_PhlApo (MBq / ml)",
	"TracerC_ParApo (MBq / ml)", "TracerJS_Apo (MBq / h)", "TracerC_SymplUpflow (MBq / ml)", "TracerC_ApoUpflow (MBq / ml)"
} ;
int NumAllConnVariablesNames = 5 ; // 5 internode connector flux variables :
string AllConnVariablesNamesQList[] = { "C_Upflow (mmol / ml)", "JS_ST (mmol / h)", "JW_ST (ml / h)", "JW_Xyl (ml / h)", "TracerJS_ST (MBq / h)" } ;

extern double *Q_ST, *Q_Sympl, *Starch, *Q_PhlApo, *Q_ParApo, *Q_out ;		  // components of vector y as used in diff. system f()...
extern double *Q_ST_dot, *Q_Sympl_dot, *Starch_dot, *Q_PhlApo_dot, *Q_ParApo_dot, *Q_out_dot ; //... and its derivatives.  ;
extern double *TracerQ_ST, *TracerQ_Sympl, *TracerStarch, *TracerQ_PhlApo, *TracerQ_ParApo ;		  // components of vector y as used in diff. system f()...
extern double *TracerQ_ST_dot, *TracerQ_Sympl_dot, *TracerStarch_dot, *TracerQ_PhlApo_dot, *TracerQ_ParApo_dot ; //... and its derivatives.  ;
extern double *vol_Sympl, *vol_Sympl_dot ;
extern Fortran_vector Psi_Xyl, JW_ST, JS_ST, C_Sympl, C_ST, C_amont, JS_Sympl, JS_Apo, JS_ParMb, JS_PhlMb, Psi_Sympl, C_PhlApo, C_ParApo, P_ST, Absorb, RespMaint ;
extern Fortran_vector JW_Trsv, JW_Xyl, JW_PhlMb, JW_ParMb, JW_Apo, JW_Sympl, P_Sympl, Psi_ST, P_Xyl, Psi_PhlApo, Psi_ParApo, P_PhlApo, P_ParApo ;
extern Fortran_vector r_Xyl, r_ST, r_ST_ref, r_Trsv, r_PhlMb, r_ParMb, r_Apo, r_Sympl, vol_ST, vol_PhlApo, vol_ParApo ;
extern Fortran_vector r_abs, PsiSoil, Transpirat, Input, StarchSyn ;
extern Fortran_vector P_ST_dot, P_Sympl_dot ; // variation rate of any variable X is noted: X_dot = dX/dt
extern Fortran_vector C_SymplUpflow, C_ApoUpflow ; // to derive JS_Sympl, and maybe JS_Apo, from resp. water fluxes JW_xxx
// tracer-specific add-ins :
extern Fortran_vector TracerJS_ST, TracerC_Sympl, TracerC_ST, TracerJS_Sympl, TracerJS_Apo, TracerJS_ParMb, TracerJS_PhlMb, TracerC_PhlApo, TracerC_ParApo ;
extern Fortran_vector TracerStarchSyn, TracerInput, TracerRespMaint, TracerC_SymplUpflow, TracerC_ApoUpflow ;
extern Fortran_vector TracerRatioSympl, TracerRatioStarch ;

// variables pointers are of type 'Fortran_vector*' (or NULL for those that are components of double* y in definition of function f).
// Check for consistency between these lists AllxxxxVariablesRefsQList hereafter and AllxxxxVariablesNamesQList[] above:
Fortran_vector* AllNodeVariablesRefsQList[] = {  // 62 node variables :
	&Absorb, NULL, &StarchSyn, &C_ST, &C_Sympl, &JS_PhlMb, &JW_Trsv, &JW_PhlMb, &JW_ParMb, &Psi_ST,
	&P_Xyl, NULL, NULL, &RespMaint, &C_SymplUpflow, &JS_Sympl, &JW_Apo, &JW_Sympl , &P_ST, &P_Sympl,
	&Psi_PhlApo, &Psi_ParApo, &Psi_Sympl, NULL, NULL, &P_PhlApo, &P_ParApo, NULL, NULL,
	&C_PhlApo, &C_ParApo, NULL, NULL, &JS_Apo, &JS_ParMb, &C_ApoUpflow,
	&P_ST_dot, &P_Sympl_dot,  NULL, NULL, NULL, &Input, &Transpirat, &PsiSoil,
	// tracer-specific add-ons :
	NULL, &TracerStarchSyn, &TracerC_ST, &TracerC_Sympl, &TracerJS_PhlMb, NULL, NULL, &TracerRespMaint,
	&TracerInput, &TracerJS_Sympl, &TracerJS_ParMb, NULL, NULL, &TracerC_PhlApo,
	&TracerC_ParApo, &TracerJS_Apo , &TracerC_SymplUpflow, &TracerC_ApoUpflow
} ;
Fortran_vector* AllConnVariablesRefsQList[] = {&C_amont, &JS_ST,  &JW_ST, &JW_Xyl, &TracerJS_ST} ; // 5 internode-connector variables

// temporary storage of results for output :
Fortran_matrix Q_ST_t, Q_Sympl_t, Starch_t, Q_PhlApo_t, Q_ParApo_t, Q_ST_dot_t, Q_Sympl_dot_t, Starch_dot_t, Q_PhlApo_dot_t, Q_ParApo_dot_t, vol_Sympl_t, vol_Sympl_dot_t ;
Fortran_matrix JW_ST_t, JS_ST_t, C_Sympl_t, C_ST_t, JS_Sympl_t, JS_Apo_t, JS_ParMb_t, JS_PhlMb_t, Psi_Sympl_t, C_PhlApo_t, C_ParApo_t, P_ST_t, Absorb_t, RespMaint_t ;
Fortran_matrix JW_Trsv_t, JW_Xyl_t, JW_PhlMb_t, JW_ParMb_t, JW_Apo_t, JW_Sympl_t, P_Sympl_t, Psi_ST_t, P_Xyl_t, Psi_PhlApo_t, Psi_ParApo_t, P_PhlApo_t, P_ParApo_t ;
Fortran_matrix PsiSoil_t, Transpirat_t, Input_t, StarchSyn_t, P_ST_dot_t, P_Sympl_dot_t, C_SymplUpflow_t, C_ApoUpflow_t, C_amont_t ;
// tracer-specific add-ins :
Fortran_matrix TracerQ_ST_t, TracerQ_Sympl_t, TracerStarch_t, TracerQ_PhlApo_t, TracerQ_ParApo_t ;		  // components of vector y as used in diff. system f()...
// Fortran_matrix TracerQ_ST_dot_t, TracerQ_Sympl_dot_t, TracerStarch_dot_t, TracerQ_PhlApo_dot_t, TracerQ_ParApo_dot_t ;		  // components of vector y as used in diff. system f()...
Fortran_matrix TracerJS_ST_t, TracerC_Sympl_t, TracerC_ST_t, TracerJS_Sympl_t, TracerJS_Apo_t, TracerJS_ParMb_t, TracerJS_PhlMb_t, TracerC_PhlApo_t, TracerC_ParApo_t, TracerC_SymplUpflow_t, TracerC_ApoUpflow_t ;
Fortran_matrix TracerPhotSynth_t, TracerStarchSyn_t, TracerInput_t, TracerRespMaint_t ;
// Fortran_matrix TracerRatioSympl_t, TracerRatioStarch_t ;

Fortran_matrix* AllNodeMatricesRefsQList[] = {
	&Absorb_t, &Starch_t, &StarchSyn_t, &C_ST_t, &C_Sympl_t, &JS_PhlMb_t, &JW_Trsv_t, &JW_PhlMb_t, &JW_ParMb_t, &Psi_ST_t,
	&P_Xyl_t, &Q_Sympl_t, &Q_ST_t, &RespMaint_t, &C_SymplUpflow_t, &JS_Sympl_t, &JW_Apo_t, &JW_Sympl_t, &P_ST_t, &P_Sympl_t,
	&Psi_PhlApo_t, &Psi_ParApo_t, &Psi_Sympl_t, &vol_Sympl_t, &vol_Sympl_dot_t, &P_PhlApo_t, &P_ParApo_t, &Q_PhlApo_t, &Q_ParApo_t,
	&C_PhlApo_t, &C_ParApo_t, &Q_PhlApo_dot_t, &Q_ParApo_dot_t, &JS_Apo_t, &JS_ParMb_t, &C_ApoUpflow_t,
	&P_ST_dot_t, &P_Sympl_dot_t, &Starch_dot_t, &Q_ST_dot_t, &Q_Sympl_dot_t, &Input_t, &Transpirat_t, &PsiSoil_t,
	// tracer-specific add-ons :
	&TracerStarch_t, &TracerStarchSyn_t, &TracerC_ST_t, &TracerC_Sympl_t, &TracerJS_PhlMb_t, &TracerQ_Sympl_t, &TracerQ_ST_t, &TracerRespMaint_t,
	&TracerInput_t, &TracerJS_Sympl_t, &TracerJS_ParMb_t, &TracerQ_PhlApo_t, &TracerQ_ParApo_t, &TracerC_PhlApo_t,
	&TracerC_ParApo_t, &TracerJS_Apo_t, &TracerC_SymplUpflow_t, &TracerC_ApoUpflow_t
};
Fortran_matrix* AllConnMatricesRefsQList[] = {&C_amont_t, &JS_ST_t,  &JW_ST_t, &JW_Xyl_t, &TracerJS_ST_t} ;
Index_vector Breakpoint_index(1, 1) ; // aux. for output time settings
double * y_dot = NULL ;	// used in f() as called by aux() which is called by odesolve() -- contains all derivatives of vector Y (whose starting values are Y0 below) below.
Fortran_vector Y0	; // (set in GUI) initial condition vector is made of Q_ST_0, Q_Sympl_0, Starch_0, Q_PhlApo_0 & Q_ParApo_0, and homologous tracer init values, and vol_Sympl_0 :
Fortran_vector atol_, rtol ; // integration accuracy parameters (set in GUI)
double t0, tf, pas;		// solve for t = [t0, tf], save result at each time step (hours)
extern double T;
float t1 ; int nbv ; // pas = initial output step = between first 2 output steps t0 and t1
bool LogScale = false ; // set to 'true' if log scale wanted (in 'PiafMunch2.cpp')
Fortran_vector OutputTimes, SegmentTimes, truc ; // will be vector of, resp.: all output times including user-added extra- or breakpoint times ; output times for current integration time segment ; and an aux.


extern int resistance_changed ;
extern int jf, i, j ; double x ; bool to_store(true) ; FILE* file = NULL ; extern string S ;
int nl = 0 ; // # of current line in matrices currently under construction in aux() -- one for each element SegmentTimes[nl] of time vector SegmentTimes
extern vector<int> k_ ; // k_[i = 1 : Nt] = index of node #i in list i_[o(i)]   where o(i) = 0 for LeafEnds (k_[i] = - k then) ; o(i) = co(i) else  -- for the moment, filled only for (o = 0 or 1)   --   k_[i] = 0 else

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

//#ifndef decoupled
//std::shared_ptr<CPlantBox::MappedPlant> plant__; 
//#endif
//extern Fortran_vector Q_ST_seg_init;
std::vector<std::string> outputVar{
	"Q_ST (mmol)", "Q_ST_dot (mmol/d)", 
	"C_ST (mmol / ml)"};


int PhloemFlux::startPM(double StartTime, double EndTime, int OutputStep,double TairK,
		std::string dir_output, bool verbose) {//, bool PMcompare, std::string file_dir 
	//std::cout<<"in startpm "<<std::endl;
	std::ofstream out("outputPM.txt");
	//cout<<"alive0 "<<endl;
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	//cout<<"alive1 "<<endl;
	if(!verbose){
		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
	}
	cout<<"extainr TairK "<<TairK<<" "<<StartTime<<" "<<EndTime<<" "<<OutputStep<<" "<<verbose<<std::endl;
    t0  = StartTime; tf  = EndTime;	nbv = OutputStep; T = TairK; TairK_phloem = TairK;
	
	cout<<"alive3 "<<endl;
	phloem_ = this->Phloem();
	
	
	cout<<"alive4 "<<endl;
	//#ifndef decoupled
	//plant__ = this->plant;
	//#endif
	
	nl = 0;
	to_store = true;
    //Launch() ; // will load settings file
	//if(PMcompare){
		
	//}
	initializePM_(tf-t0, TairK);
	initialize_carbon(this->Q_outv) ;		// sizes C-fluxes-related variable vectors
    initialize_hydric() ;		// sizes water-fluxes-related variable vectors and sets up hydric system
	prepareOutput(dir_output);
	
	// building up output times vector OutputTimes according to GUI seetings :
    pas = (tf-t0)/double(nbv);
	t1 = t0 + pas ; //nbv = (int)((tf-t0)/pas) ; // pas = initial output step = between first 2 output steps t0 and t1 as set in GUI
	std::cout <<"set out put vec "<<t0<<" "<<t1<<" "<<tf<<" "<<pas<<" "<<nbv<<std::endl;
	//t1 = tf;
	//nbv = 1;
	//TimeSegmentConfig() ; // User-editable (implemented in 'PiafMunch2.cpp') : set of Log Scale and/or addition of extra output- or breakpoint- times
    OutputSettings("test") ; // sets up file and graph outputs, including possible breakpoint- and/or extra-output- times as stated above

    // *************** SOLVING THE DIFFERENTIAL EQUATION SYSTEM ************************************* :
    int neq = Y0.size() ;							// number of differential eq. = problem size (= 8*Nt après ajouts FAD)
    cout<<"neq "<<neq<<" "<<Nc<<" "<<Nt<<endl;
	
	assert((Nt == (neq/neq_coef))&&"Wrong seg and node number");
	//assert((Nc == (neq /3-1))&&"Wrong seg and node number");
	assert((Nc == (Nt -1))&&"Wrong seg and node number");
	y_dot = new double[1 + neq] ;			// for use in aux()
    
    // -------------  To compute P_dot = dP/dt  and  P_symp_dot = dP_Sympl/dt  and make them available to all modules in real time : ---------------------
    // 1°) Oversize both 'Var_integrale' and 'Var_derivee' pointers and initialize all of them to NULL :
    Fortran_vector** Var_integrale = new Fortran_vector*[100] ; Fortran_vector** Var_derivee = new Fortran_vector*[100] ;
    for (i = 0 ; i < 100 ; i ++) {Var_integrale[i] = Var_derivee[i] = NULL ;}
    // 2°) initialize pointers to derivatives (and corresponding integrals) that are actually used (or may be so), in this case for...
    Var_integrale[0] = &P_Sympl ; Var_derivee[0] = &P_Sympl_dot ; //...elastic changes of symplastic volume in function (PiafMunch2.cpp)Smooth_Parameter_and_BoundaryConditions_Changes()
    Var_integrale[1] = &P_ST ; Var_derivee[1] = &P_ST_dot ;
		//throw std::runtime_error("Breakpoint_index.size() ");

    for(is = 1 ; is < Breakpoint_index.size() ; is ++) { // allows several integration segments in relation to breakpoints (if any -- OK if none)
        SegmentTimes = subvector(OutputTimes, Breakpoint_index(is), Breakpoint_index(is+1))  ; // time segment between 2 breakpoints (Breakpoint_index(1) = 1 ; Breakpoint_index(Breakpoint_index.size()) = 1 + nbv)
        cout << endl << "starting integration on time segment #" << is << " = [" << SegmentTimes[1] << ", " << SegmentTimes[SegmentTimes.size()] << "] "<< Breakpoint_index.size()<<" "<<SegmentTimes.size()<< endl ;
        double t = SegmentTimes[1] ; //BreakpointSharpParameterChanges(is, t) ;
		resistance_changed = XYL | TRSV | APO | SYMPL | PHL |  PHLMB | PARMB | RABS ;
		
		//UpdateResistances(t) ;
        
		resistance_changed = NONE ;
        // NB : Y0 is updated at each output time point (hence, at each breakpoint) in function aux()

        // ***** the following solver configs are ranked from the most efficient (in most tested situations) to the least (in most tested situations). Can change in different situations ! *****
        //   see  SUNDIALS  documentation  for  cvode  solver options (SPxxxx, xxxx_GS, PREC_xxxx, BAND, etc.)
		cout <<"solver, Y0: "<<solver <<endl;
		//for (i = 1 ; i <= Y0.size() ; i++) {cout<<" "<<Y0[i];}//Y0: all the data in 1 vector, so far so good
		//cout<<endl<<"end Y0"<<endl;//f, aux,
		//atol_.display(); rtol.display();
		//throw std::runtime_error("avant cvode ");
        switch (solver) {
		//j = cvode_spils(fout, Y0, SegmentTimes,auxout,   atol_, rtol, SPFGMR, MODIFIED_GS, PREC_NONE, 2, Var_integrale, Var_derivee); 
			case 1: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, MODIFIED_GS, PREC_NONE, 2, Var_integrale, Var_derivee); break; // STALD = true, verbose = true, rootfind = NULL
			case 2: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, MODIFIED_GS, PREC_NONE, 2, Var_integrale, Var_derivee); break; // STALD = true, verbose = true, rootfind = NULL
			case 3: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, CLASSICAL_GS, PREC_NONE, 2, Var_integrale, Var_derivee); break; // le + rapide : best for large N (>= 1000) ; break ; sinon, un peu instable avec VarVisc
			case 4 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, CLASSICAL_GS, PREC_NONE, 2, Var_integrale, Var_derivee) ; break ; // le + rapide : best for large N (>= 1000) ; break ; sinon, un peu instable avec VarVisc
				case 5 : j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, DIAG, 2, Var_integrale, Var_derivee) ; break ; // TB pour Thompson, même avec vol_Sympl_dot ...
			case 6 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPBCGS, 1, PREC_NONE, 2, Var_integrale, Var_derivee) ; break ;
				case 7: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, BAND, 2, Var_integrale, Var_derivee, true, true, NULL, 0, neq / 30, neq / 30); break; // mu = ml = neq/30 : TB
			case 8: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, PCG, 1, PREC_NONE, 2, Var_integrale, Var_derivee); break; //
			case 9 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPTFQMR, 1, PREC_NONE, 2, Var_integrale, Var_derivee) ; break ; //
			case 10 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, MODIFIED_GS, PREC_LEFT, 2, Var_integrale, Var_derivee) ; break ; // STALD = true, ... (id. ci-dessus) :  bon choix en général, mais pas avec vol_Sympl_dot !
			case 11 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, CLASSICAL_GS, PREC_LEFT, 2, Var_integrale, Var_derivee) ; break ; //
			case 12: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, CLASSICAL_GS, PREC_RIGHT, 2, Var_integrale, Var_derivee); break; // le + rapide : best for large N (>= 1000) ; break ; sinon, un peu instable avec VarVisc
			case 13: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, MODIFIED_GS, PREC_BOTH, 2, Var_integrale, Var_derivee); break; // STALD = true, verbose = true, rootfind = NULL
			case 14: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPFGMR, CLASSICAL_GS, PREC_LEFT, 2, Var_integrale, Var_derivee); break; // le + rapide : best for large N (>= 1000) ; break ; sinon, un peu instable avec VarVisc
			case 15 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPGMR, MODIFIED_GS, PREC_RIGHT, 2, Var_integrale, Var_derivee) ; break ; // STALD = true, ... (id. ci-dessus)  :  bon choix en général, mais pas avec vol_Sympl_dot
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
				case 28: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, BAND, 2, Var_integrale, Var_derivee, true, true, NULL, 0, neq / 100, neq / 100); break; // mu = ml = neq/100 : marche très bien même si neq < 100 !
				case 29: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, BAND, 2, Var_integrale, Var_derivee, true, true, NULL, 0, neq / 10, neq / 10); break; // mu = ml = neq/10 : OK
				case 30: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, BAND, 2, Var_integrale, Var_derivee, true, true, NULL, 0, neq / 3, neq / 3); break; // pas de rootfind, ; break ; mu = ml = neq/3 : OK
				case 31: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, DENSE, 2, Var_integrale, Var_derivee); break; // solver = cvode DENSE, STALD = true
				case 32: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, KLU, 2, Var_integrale, Var_derivee); break; //
			case 33 : j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPTFQMR, 1, PREC_LEFT, 2, Var_integrale, Var_derivee) ; break ; //
				case 34: j = cvode_direct(fout, Y0, SegmentTimes, auxout, atol_, rtol, BAND, 2, Var_integrale, Var_derivee, true, true, NULL, 0, neq / 300, neq / 300); break; //
			case 35: j = cvode_spils(fout, Y0, SegmentTimes, auxout, atol_, rtol, SPTFQMR, 1, PREC_BOTH, 2, Var_integrale, Var_derivee); break; //
			default : cout << endl << "!! Error !! solver # must be within the range [1 , 35] !!" << endl ; j = -1 ; exit(-1) ;
        }		
		//j = cvode_spils(fout, Y0, SegmentTimes,auxout,   atol_, rtol, SPFGMR, MODIFIED_GS, PREC_NONE, 2, Var_integrale, Var_derivee); 
		//j = cvode_spils(f, Y0, SegmentTimes,aux,   atol_, rtol, SPFGMR, MODIFIED_GS, PREC_NONE, 2, Var_integrale, Var_derivee); 
		if (j < 0) return (-1); // solver init error
		if (j > 0) break ; // simulation run-time error
        //to_store = false ; // do not store next datapoint, i.e. the first datapoint of next integration segment (if any), as this will be the same as the last one in this integration segment
    std::cout << std::endl << "ending integration on time segment [" <<std::endl;//<< SegmentTimes[1] << ", " << SegmentTimes[SegmentTimes.size()] << "] : ii = " << ii << " ; Y0 = " ; Y0.display() ;
    }
    // ********************* OUTPUT *********************************
    //if (file != NULL) {fclose(file) ; file = NULL;}
	
	//std::cout<<"output0 "<<Y0.size()<<std::endl;
    // MEMORY LIBERATIONS:
    delete [] y_dot;
	//if(!AutoQuit) system("pause") ;
	//cout.rdbuf(coutbuf); //reset to standard output again
	//std::cout<<"output1 "<<std::endl;
	//Y0.display();
	//C_ST.display();
	// Q_Sympl.display();
	// Starch.display();
	// Q_PhlApo.display();
	// Q_ParApo.display();
	//std::cout<<"output2 "<<std::endl;
			
	this->Q_outv = Y0.toCppVector();//att: now has several outputs
	
	this->Q_out_dotv = std::vector<double>(Y0.size(),0);
	for(int j = 1 ; j <= Y0.size() ; j ++) 
	{
		this->Q_out_dotv[j-1] = y_dot[j] ; 
	}
	//std::cout<<"output3 "<<std::endl;
	this->C_STv = C_ST.toCppVector();//att: now has several outputs
	//std::cout<<"output4 "<<std::endl;
	this->vol_STv = vol_ST.toCppVector();
	//std::cout<<"output5 "<<std::endl;
	//this->Q_Grmaxv = Q_Grmax.toCppVector();
	this->Q_Exudmaxv = Q_Exudmax.toCppVector();
	this->Q_Rmmaxv = Q_Rmmax.toCppVector() ;
	this->Flv = Input.toCppVector() ;
	this->r_STv = r_ST.toCppVector() ;
	this->JW_STv = JW_ST.toCppVector() ;
	//std::cout<<"computeOrgGrowth"<<std::endl;
	//C_ST.display();
	computeOrgGrowth(tf-t0);//set Gr sink to 0
		//std::cout<<"aliveendend "<<std::endl;
	if(!verbose){std::cout.rdbuf(coutbuf);} //reset to standard output again
		//std::cout<<"aliveendend2 "<<std::endl;
	return(1) ;
	
}

void OutputSettings(string pos)  {
	//int k ; 
	OutputTimes = Fortran_vector(nbv + 1) ;
	OutputTimes[1]=t0 ; if(nbv > 0) OutputTimes[2] = t1 ;
	for(i=2 ; i <= nbv ; i++) {
		if(LogScale) pas *= rs ;
		OutputTimes[i+1] = OutputTimes[i] + pas ;
	}
	OutputTimes[nbv + 1] = tf;
	//if (OutputTimes[nbv + 1] < tf) { // may happen even if step was set as 'pas = (tf-t0)/nbv',  because of rounding errors
		
		//nbv ++ ;
		//OutputTimes.append(Fortran_vector(1, tf)) ;
	//}
	//else  OutputTimes[nbv + 1] = tf ; // fixes rounding errors
		
	if(Breakpoint_index[Breakpoint_index.size()]  != OutputTimes.size())  Breakpoint_index.append(Index_vector(1, OutputTimes.size())) ;
	cout << "Output times :" << endl ; for(i=0 ; i < OutputTimes.size() ; i ++) cout << " " <<  OutputTimes[i+1] ; // if newly inserted extra output time (just after breakpoint) is too close, it may print as the same but actually is not
	cout << endl << endl << "breakpoint times :" << endl ; for(i=1 ; i <= breakpoint_times.size() ; i ++) cout << " " <<  OutputTimes[Breakpoint_index[i+1]] ;
	cout << endl << "AutoQuit "<<AutoQuit<<" "<<nsp<<" "<<nvp<<" "<<fsp<<" "<<fvp<<endl;
	nbv = OutputTimes.size() - 1 ; // at this point,  output times vector truc = OutputTimes ; and Breakpoint_index = {1, [...,] nbv+1}  with  [...,] = extra breakpoints indices (empty if none)
	
	cout<<endl<<"prepared outputnames "<<endl;
	//const time_t current = time(NULL) ; strftime(name_, 50, "_%Y-%m-%d_%Hh-%Mmn-%Ss_output.txt", localtime(&current)) ;
}

void PhloemFlux::prepareOutput(std::string dir_output)  {
	
	if(nbv > 1)
	{
		file = fopen(dir_output.c_str(), "w") ; assert(file) ;
		//std::cout<<"dir_output "<<dir_output<<std::endl;
		fprintf(file, "time") ;
		
		//std::cout<<"expression "<<expression<<std::endl;
		switch(expression) {
			case 1:
				std::vector<std::string> outputVar_ { 
					"Q_ST (mmol)", "Q_ST_dot (mmol/d)", 
					"C_ST (mmol / ml)","Q_Sympl (mmol)",
					"Input (mmol / h)", "Starch (mmol eq. Glu)",
					 "Q_ParApo (mmol)","Q_PhlApo (mmol)"};
				outputVar = outputVar_;
			break;
		}
				
		for(const auto& outVar : outputVar ){
			//std::cout<<outVar<<" "<<Nt<<std::endl;
			for (int its = 0 ; its < Nt ; its++) {
				QS = "\t" + outVar  + "["+std::to_string(its)+"]" ;
				//std::cout<<QS<<std::endl;
				fprintf(file, QS.c_str());
			}
		}
		fprintf(file," \n");
	}
	//throw std::runtime_error( "Root::simulate: p.ln.at(i) - length < dxMin");
	
}

void PhloemFlux::aux(double t, double * y) {	// launch auxiliary calculations to store dynamics of temporary variables, , shared_ptr<PhloemFlux> PhloemObject
/*list of data to print:
	C_ST, JS_ST, JW_ST

*/
	//cout<< t <<" t == t0: " << bool(t == t0) <<endl;
	//int j;//i,  k, ik ;
	//if (t == t0) { // true simulation start point --> reset all primary variables that are not direct components of Y0 to make sure the system is balanced :
		//cout << "aux() : initializing at t = t0 : nl = " << nl << endl ;
		//JS_PhlMb.set(0.) ; C_amont.set(C_ST[1]) ; // arbitrary, for Adv_BioPhysics...
		//f(t, y, y_dot) ;
	//}
	f(t, y, y_dot) ;			// Update all variables from Y as returned by solver
	//std::cout<<"in aux "<<Nt<<std::endl;
	if(nbv > 1)
	{
		fprintf(file, "%g", t) ;
		for (int it = 1 ; it <= Nt  ; it++){
			fprintf(file, "\t%g", Q_ST[it]) ;
			if(Q_ST[j]<0.)
			{
				std::cout<< "at t = " << t << " : Y0.size() = " << Y0.size()<<std::endl;
				Y0.display();
				std::cout<<"negative Q_ST value at j = "<<j<<", Q_ST[j] = "<< Q_ST[j] <<std::endl;
				assert(false);
			}
		}
		//std::cout<<"expression "<<expression<<std::endl;
		switch(expression) {
			case 1:
				for (int it = 1 ; it <= Nt  ; it++){
					fprintf(file, "\t%g", Q_ST_dot[it]) ;
				}
				for (int it = 1 ; it <= Nt  ; it++){
					fprintf(file, "\t%g", C_ST[it]) ;
				}
				for (int it = 1 ; it <= Nt  ; it++){
					fprintf(file, "\t%g", Q_Sympl[it]) ;
				}
				for (int it = 1 ; it <= Nt  ; it++){
					fprintf(file, "\t%g", Input[it]) ;
				}
				for (int it = 1 ; it <= Nt  ; it++){
					fprintf(file, "\t%g", Starch[it]) ;
				}
				for (int it = 1 ; it <= Nt  ; it++){
					fprintf(file, "\t%g", Q_ParApo[it]) ;
				}
				for (int it = 1 ; it <= Nt  ; it++){
					fprintf(file, "\t%g", Q_PhlApo[it]) ;
				}
			break;
		}
		
		fprintf(file," \n");
	}
	
	for(int j = 1 ; j <= Y0.size() ; j ++) 
	{
		Y0[j] = y[j] ; 
		if((j <= Nt)&&(Y0[j] < 0.)){assert(false && "negative C_ST value");}
	}// update Y0
	std::cout << "at t = " << t << " : Y0.size() = " << Y0.size();// << " ; r_ST =" ; //for(i=1 ; i <= Nc ; i ++) cout << " " <<  r_ST[i] ; cout << endl ; // to be edited according to parameter change at breakpoint
	std::cout<<std::endl;
	//assert(false);
}


//still to do:
//1) make that input come from plant
//2) add evapotranspiration in f() function
//3) to exteroir loop?
//4) add Gr and Rm
//5) 
void PhloemFlux::initializePM_(double dt, double TairK){	
/*	____________________________________________________________________________________________
	
						cell-specific				cell-shape, type			set constant
	____________________________________________________________________________________________
	param exud: |	exud_k					|	rad, len, ot		|
	------------|---------------------------|-----------------------|---------------------------
	param Rm: 	|	krm1, krm2, rhoC		|		volST			|
	------------|---------------------------|-----------------------|---------------------------
	param Gr: 	|	Rmax_org, Lmax_org		|	len, ot				|	Y
				|	 isTip					|						|
	------------|---------------------------|-----------------------|---------------------------
	source: 	|	Ag, k_meso				|	rad, len, ot		|	 
				|							|						|
	
*/
	//Att! all is in mmol/ml (=> M) and in d-1
	cout<<"initializePM_1 "<<endl;
	Adv_BioPhysics = false;
	vector<CPlantBox::Vector2i> segmentsPlant = plant->segments;
	Nc = segmentsPlant.size(); Nt = plant->nodes.size();
	assert((Nc == (Nt -1))&&"Wrong seg and node number");
	atol_= Fortran_vector(Nt*neq_coef, atol_double);//1e-017
	r_ST_ref = Fortran_vector(Nc)	; 
	I_Upflow = I_Downflow = vector<int>(Nc +1) ;
	a_STv.resize(Nc , 0.)  ;
	//StructSucrose.resize(Nc + 1, 0.)  ;
	//vector<int> tips = plant->get_organ_seg_tips();//seg id
	vector<int> orgTypes = plant->organTypes;//per seg
	vector<int> subTypes = plant->subTypes;
	vector<double> Radii = plant->radii;
	vector<double> dist2tip = plant->dist2tips;
	//std::cout<<"plant->segLength() "<<std::endl;
	vector<double> Lengthvec = plant->segLength();
	//vector<double> sideSurface = Flen();
	rtol =  Fortran_vector(1, rtol_double);//1e-023
	double a_seg, l;double mu ;
	int ot, st;
	Ag =Fortran_vector(Nt, 0.);
	len_leaf =Fortran_vector(Nt, 0.);
	vol_ParApo =Fortran_vector(Nt, 0.);
	Q_Grmax =Fortran_vector(Nt, 0.);
	Q_Exudmax=Fortran_vector(Nt, 0.);
	Q_Rmmax =Fortran_vector(Nt, 0.);
	vol_ST=Fortran_vector(Nt, 0.);
	exud_k=Fortran_vector(Nt, 0.);
	krm2=Fortran_vector(Nt, 0.);//0.04
	//rmaxPotv = rmaxSeg(dt, k_gr);
	deltaVolOrgNode_ = calcDeltaVolOrgNode(dt, k_gr);
	//rmaxPot = Fortran_vector(rmaxPotv);
	//vol_ST.display(); Q_Rmmax.display();
	cout<<"initializePM_new "<<Nc<<" "<<Nt<<endl;
	int nodeID;
	
	/*straw bulk density: mean(24, 111) kg m3 => 0.07 g cm-3 #http://biomasslogistics.org/Publications/22lam.pdf
    // mol C/cm3 fresh plant: [0.07 g/cm3freshmatter] [0.4 gC/gdryplant (cf elemental analysis ackelysimeter)]*[0.4 gDryplant/gwetplant (cf erntetabelle ackelysimeter)]\
    // *[1/12 molC/gC]content (mol) per unit of volume
	// *[1/12 molSuc/molC]*1000 mmol/mol*/
	double StructSucrose, deltaStructSucrose;//, sideSurface_seg;//deltaVol,
			
	double cmH2O_to_hPa = 0.980638	;//cm to hPa
	double krm1;
	int numleaf = 0;
	//std::cout<<"before psiXyl.size()"<<std::endl;
	if(psiXyl4Phloem.size() == Nt){Psi_Xyl = Fortran_vector(psiXyl4Phloem,cmH2O_to_hPa);
	}else{Psi_Xyl = Fortran_vector(Nt, 0.);}// (id.) Phloem water potential	
	//std::cout<<"after psiXyl.size()"<<std::endl;
	
	if(not update_viscosity_)
	{	

		double TdC = TairK - 273.15;
		//in g/L or mg/cm3
		dEauPure = (999.83952 + TdC * (16.952577 + TdC * (- 0.0079905127 + TdC * (- 0.000046241757 + TdC * (0.00000010584601 + TdC * (- 0.00000000028103006)))))) / (1 + 0.016887236 * TdC); 
		double siPhi = (30 - TdC) / (91 + TdC) ; // T_old = TairK_phloem ;//  R.Gilli 1997, after Mathlouthi & Génotelle 1995 - valid for any T :
				
		double C = 0.5 ; // (mmol / ml solution)
		//  R.Gilli 1997, after Mathlouthi & Génotelle 1995 - valid for any T :
		//342.3 g/mol or mg/mmol
		double PartMolalVol_ =0;// 0.2155;
		double d = C * 342.3 + (1 - C * PartMolalVol_) * dEauPure ;//in mg/cm3
		double siEnne = (100 * 342.30 * C) / d ; // actually this is sc = sucrose content (g.suc. % g.solution) ; 342.30 = molar mass of sacch.
		siEnne /= 1900 - (18 * siEnne) ;
		//mPa s
		mu =  pow(10, ((22.46 * siEnne) - 0.114 + (siPhi * (1.1 + 43.1 * pow(siEnne, 1.25) )))) ; // peut atteindre des valeurs > 1.e200 !! (sans signification évidemment -- le sucre doit précipiter bien avant !!)
		//std::cout<<"visco water "<<TdC<<" "<<mu<<std::endl;
		mu = mu /(24*60*60)/100/1000; //mPa s to hPa d, 1.11837e-10 hPa d for pure water at 293.15K
		
	}
	for ( int k=1 ;k <= Nc;k++ ) 
	{
			ot = orgTypes[k-1]; st = subTypes[k-1];
			a_seg = Radii[k-1];
			l = Lengthvec[k-1];
			I_Upflow[k] = segmentsPlant[k-1].x +1;
			I_Downflow[k] = segmentsPlant[k-1].y +1;
			//if(ot==4){st += 2;}	
			nodeID = I_Downflow[k];
			len_leaf[nodeID] = l;
			//std::cout<<"ot st "<<ot<<" "<<st<<" "<< k <<" "<<orgTypes[k-2]<<" "<<orgTypes[k]<<" st "<<subTypes[k-2]<<" "<<subTypes[k]<<std::endl;
			double Across_ST = get_Across_ST(st,ot);
			//std::cout<<"a_ST_eqs "<<a_ST_eqs.size()<<//std::endl;
			vol_ST[nodeID] = Across_ST * l;// * M_PI;//a_ST[k-1] * a_ST[k-1] * l * M_PI;
			Across_STv[k-1] = get_Across_ST(st,ot);
			//std::cout<<"a_st_seq "<<a_ST_eq<<std::endl;
		//general
			//double waterViscosity = 2.414e-5 * std::pow(10,247.8/(TairK- 140))/(3600*24)/100; //in [hPa d] https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
			//1.7mPa s => 1.7 /(3600*24)/1000 /100 = 1.967593e-10 hPa d
			//mu =  1.967593e-10;//mean  viscosity in work of XZ//hPa d-1
			
			//double surface_ST = M_PI * a_ST[k-1] *a_ST[k-1];//cm2
			r_ST_ref[k] = 1/kx_suc_f(  st,ot)*l;//(1/kst)*(l/surface_ST);//320;//1405;//mu/kst*l/a;
			if(not update_viscosity_)
			{
				r_ST_ref[k] = mu*r_ST_ref[k];
			}
			//std::cout<<"r_ST "<<kx_suc_f(  st,ot)<<" "<<mu<<" "<<r_ST_ref[k]<<std::endl;//" "<<M_PI<<" "<<surface_ST<<" "<<l<<" "<<kst<<" "<<mu<<" "<<r_ST_ref[k]<<std::endl;
			//assert(false);
			//if(add_gravity){Psi_Xyl[k] += vz*g;} //add gravity potential
		//Fl
			/*mean val: [10, 20] micromol C m-2 s-1 => 15 /12 * 24*60*60 *0,001 /10000 mmol Sucrose cm-2 d-1*/
			//in mmol Suc d-1
			Ag[nodeID] = Ag4Phloem[nodeID-1] ;//* sideSurface[nodeID-1];//for leaf: ratio of leaf area
			//6.40314e-006 * (ot == 4);
			/*mol CO2 m-2 s-1 => mol C m-2 s-1
			//6.174359149662036e-05
			//X /12 * 24*60*60 *1000/10000 = X2 mmol Suc cm-2 d-1
			//0.04445539*/
			
		//Fu
			//Exudation
			// if(exud_kv.size()==Nt)
			// {exud_k[nodeID] = exud_kv[nodeID-1]*(ot == 2)+0.;
			// }else{exud_k[nodeID] = exud_kv[0]*(ot == 2)+0.;}
			exud_k[nodeID] =kr_suc_f(st,ot);
			//double sideSurface_st = 2*M_PI*l*a_ST[k-1]; //pour exud predre surface de root
			//std::cout<<"exud "<<exud_k[nodeID]<<std::endl;
			Q_Exudmax[nodeID ]=0;
			if((ot==2)&&(dist2tip[k-1] < l_kr)){ //in root and in length for exud?
				double l_exud = std::min(l,l_kr - dist2tip[k-1]);
				Q_Exudmax[nodeID ]= 2*M_PI*l_exud*a_seg;
			}//sideSurface_st*(ot == 2);//for root: cylinder side surface
			//std::cout<<"QexudMax "<<Q_Exudmax[nodeID ]<<std::endl;
			
			//Rm
			if(krm1v.size()==Nt)
			{krm1 = krm1v[nodeID-1];
			}else{krm1 = krm1v[0];}
			
			if(krm2v.size()==Nt)
			{krm2[nodeID] = krm2v[nodeID-1];
			}else{krm2[nodeID] = krm2v[0];}
			double segVol;
			vol_ParApo[nodeID] = -1;//dummy value to avoid division by 0 when cimputing cnncentration in mesophyle
			if(ot == 4){
				segVol = plant->segVolLeaf[numleaf];//a_seg;
				
				vol_ParApo[nodeID] = l * surfMeso;
				// if(segVol == 0)//basal zone
				// {
					// segVol = l*(a_seg*a_seg)*M_PI;
				// }
			}else{segVol = l*(a_seg*a_seg)*M_PI;}
			if(sameVolume_meso_st){
					vol_ParApo[nodeID] = vol_ST[nodeID];
				}
			StructSucrose = rhoSucrose[ot] * segVol;//* 
			//std::cout<<"segvol "<<segVol<<" "<<l<<" "<<a_seg<<" "<<StructSucrose<<std::endl;
			
			Q_Rmmax[nodeID] = krm1 * StructSucrose;
			
			//Gr
			//deltaVol = rmaxPotv[nodeID -1]* a_seg * a_seg * M_PI ;
			//deltaStructSucrose = deltaVol * rhoSucrose;
			//in delta_cm3 to mmol Suc d-1
			deltaStructSucrose = deltaVolOrgNode_[k][-1]* rhoSucrose[ot];//total sucrose rate needed for growth on that node
			//std::cout<<deltaVolOrgNode_[k].find(-1)->second<<" "<<deltaVolOrgNode_[k][-1]<<std::endl;
			// "/dt" to go from max change in cm3 to change rate in cm3 d-1
			Q_Grmax[nodeID] = deltaStructSucrose/Gr_Y/dt;//*k_gr2;
		
		//Test
			if(Q_Grmax[nodeID ]<0.){
				std::cout<<"gr: loop n°"<<k<<", node "<<nodeID<<" "<<Q_Grmax[nodeID]<<" "<< deltaVolOrgNode_[k][-1]<<std::endl;
				assert(false);
			}
			if(Q_Exudmax[nodeID ]<0.){
				std::cout<<"exud: loop n°"<<k<<", node "<<nodeID<<" "<<Q_Exudmax[nodeID]<<" "<< l<<" "<<Radii[k-1]<<std::endl;
				assert(false);
			}
			if(Q_Rmmax[nodeID ]<=0.){
				std::cout<<"rm: loop n°"<<k<<", node "<<nodeID<<" "<<Q_Rmmax[nodeID]<<" "<< krm1<<" "<<StructSucrose<<std::endl;
				assert(false);
			}
			//
		}
	//for first node:
	//a_ST[0] = a_ST[1];
	len_leaf[1] = len_leaf[2];
	vol_ParApo[1] = vol_ParApo[2];
	vol_ST[1] = vol_ST[2];//Length[1] = Length[2]; radius_ST[1] =radius_ST[2];
	Q_Rmmax[1] = Q_Rmmax[2];
	Q_Grmax[1] =  deltaVolOrgNode_[0][-1]* rhoSucrose[orgTypes[0]]/Gr_Y/dt;
	//krm2[1]=0;
	//vol_ST.display();
	//r_ST_ref.display();
	//vol_ParApo.display();
	// assert(false);
	this->Q_Grmaxv = Q_Grmax.toCppVector();
	this->r_ST_refv = r_ST_ref.toCppVector();
	this->Agv = Ag.toCppVector();
	this->vol_STv = vol_ST.toCppVector();
	this->vol_Mesov = vol_ParApo.toCppVector();
}



void rootfind(double t, double* y, double* g) {return ; } // required but not used by cvode_xxxxx()

void PhloemFlux::computeOrgGrowth(double t){
	//int nOrg = plant->getNumberOfOrgans();
	int nNode = plant->getNumberOfNodes();
	//std::vector<double> orgGr(nOrg, 0.);
	auto orgs = plant->getOrgans(-1, true); //get all organs, even small ones
	int nOrg = orgs.size();
	//double orgGr;
	std::map<int, double> cWGrRoot; //set cWGr in this function instead of in Mappedorganism.h : cWGr is then reset to empy every tim efunction is called
	std::map<int, double> cWGrStem; // + no need for mapped organism to keep cWGr in memory. just has to be in growth function
	std::map<int, double> cWGrLeaf;
	
	delta_suc_org.resize(plant->getNumberOfOrgans(),0.);
	//delta_suc_org_i = std::vector<double>(nOrg, 0.);
	//delta_suc_org_max.resize(nOrg,0.);
	//delta_suc_org_imax = std::vector<double>(nOrg, 0.);
	
	delta_ls_org.resize(nOrg,0.);
	delta_ls_org_i = std::vector<double>(nOrg, 0.);
	delta_ls_org_max.resize(nOrg,0.);
	delta_ls_org_imax = std::vector<double>(nOrg, 0.);
		
	delta_vol_org.resize(nOrg,0.);
	delta_vol_org_max.resize(nOrg,0.);
	delta_vol_node.resize(nNode,0.);
	delta_vol_node_i = std::vector<double>(nNode, 0.);
	delta_vol_node_max.resize(nNode,0.);
	delta_vol_node_imax = std::vector<double>(nNode, 0.);
	double delta_volOrg;
	double delta_volOrgmax;
	int orgID2 = 0;
	for(auto org: orgs){
		//orgGr = 0;
		int orgID = org->getId();
		//double a = org->getParameter("a");//radius
		//double crossSurface = M_PI * a * a;
		std::vector<int> nodeIds_;
		double ot = org->organType();
		int nNodes = org->getNumberOfNodes();
		//std::cout<<"org "<<orgID<<" "<<a<<" "<<ot<<" "<<nNodes<<std::endl;
		if(ot == 2){
			//nodeIds_.push_back(-1); not needed anymore as take into account seed node
			//std::cout<<"root id schemw "<<std::endl;
			nodeIds_.push_back(org->getNodeId(nNodes-1));			
			}else{nodeIds_ = org->getNodeIds();}
		
		delta_volOrg = 0;
		delta_volOrgmax = 0;
		std::vector<int> allnodeids = org->getNodeIds();
		//std::cout<<"nodeids1: "<<allnodeids.size()<<" "<< nNodes <<" "<<nodeIds_.size()<<" "<< (org->getNodeId(nNodes-1)) <<std::endl;
		
		//for(int k=0;k< allnodeids.size();k++)//take  first node as might carry small laterals
		//{
			//std::cout<< allnodeids[k] <<" ";
		//}
		//std::cout<<std::endl<<"id 2 "<<std::endl;
		
		//for(int k=0;k< nodeIds_.size();k++)//take  first node as might carry small laterals
		//{
			//std::cout<< nodeIds_[k] <<" ";
		//}
		//std::cout<<std::endl;
		for(int k=0;k< nodeIds_.size();k++)//take  first node as might carry small laterals
		{
			////if(ot==3){std::cout<<"output growth "<<k<<" "<<nodeIds_[k]<<" "<< orgID <<std::endl;}
			double deltaVolmax_ = 0;
			if(deltaVolOrgNode_[nodeIds_[k]].count(org->getId()) !=0){
				deltaVolmax_ = deltaVolOrgNode_[nodeIds_[k]].find(org->getId())->second;
			}
			//if(ot==3){std::cout<<deltaVolOrgNode_[nodeIds_[k]].find(org->getId())->second<<" "<<deltaVolOrgNode_[nodeIds_[k]][org->getId()]<<std::endl;}
			//if(ot==3){std::cout<<deltaVolOrgNode_[nodeIds_[k]].find(-1)->second<<" "<<deltaVolOrgNode_[nodeIds_[k]][-1]<<std::endl;}
			double deltaSucmax_ = deltaVolmax_ * rhoSucrose[ot]/Gr_Y;//to differentiate between parent organs and its laterals to small to be represented
			//if(ot==3){std::cout<<"		"<<deltaVolmax_<<" "<<deltaSucmax_<<std::endl;}
			//if(ot==3){std::cout<<"		"<<Q_ParApo[nodeIds_[k] +1]<<" "<< Q_ParApoBU[nodeIds_[k] +1]<<std::endl;}
			double deltaSuc = min(Q_ParApo[nodeIds_[k] +1] - Q_ParApoBU[nodeIds_[k] +1], deltaSucmax_);
			//double delta_Sucmax = min(Q_Grmax[nodeIds_[k] +1]*t, deltaSucmax_);
			//if(ot==3){std::cout<<"deltaSuc "<<deltaSuc<<" "<<delta_volOrg<<std::endl;}
			assert((deltaSucmax_ <= Q_Grmax[nodeIds_[k] +1])&&"error Q_Grmax");
			if((deltaSuc < 0.)&&(deltaSuc > -1e-5)){deltaSuc=0.;}
			if((deltaSuc < 0.)){
				//std::cout<<"deltasuc "<<deltaSuc<<std::flush<<std::endl;
				assert((deltaSuc >= 0.) &&"negative deltaSuc");
			}
			delta_volOrg += deltaSuc*Gr_Y/rhoSucrose[ot];
			//if(ot==3){std::cout<<"deltaSuc "<<deltaSuc<<" "<<delta_volOrg<<" "<<Gr_Y<<" "<<rhoSucrose[ot]<<std::endl;}
			delta_volOrgmax += deltaSucmax_*Gr_Y/rhoSucrose[ot];
			
			
			Q_ParApoBU[nodeIds_[k] +1] += deltaSuc;//Q_ParApo[nodeIds_[k] +1];//sucrose is spent
			
			//orgGr += delta_l; //[orgID]
			Q_GrmaxBU[nodeIds_[k] +1] += deltaSucmax_;//sucrose is spent
			
			delta_vol_node_i[nodeIds_[k]] += deltaSuc*Gr_Y/rhoSucrose[ot];
			delta_vol_node[nodeIds_[k]] += deltaSuc*Gr_Y/rhoSucrose[ot];
			
			delta_vol_node_imax[nodeIds_[k]] += deltaSucmax_*Gr_Y/rhoSucrose[ot];
			delta_vol_node_max[nodeIds_[k]] += deltaSucmax_*Gr_Y/rhoSucrose[ot];
			
			delta_suc_org[orgID] += deltaSuc*Gr_Y;
			
		}
		
		double vol = org->orgVolume(-1, false);
		double l_ =  org->getLength(false);
		//if(ot==3){std::cout<<"vol new vol "<<ot<<" "<<vol<<" delta_volOrg: "<<delta_volOrg<<" "<<l_<<std::endl;}
		double vol_check = org->orgVolume(l_, false);
		double l_check = org->orgVolume2Length(vol);
		//if(ot==3){std::cout<<"vol, l, checks: "<<ot<<" "<<l_<<" "<<l_check<<" "<<(l_ - l_check) <<" "<<vol<<" "<<vol_check<<" "<<(vol - vol_check)<<std::endl;}
		if((abs(l_ - l_check)>1e-5) || (abs(vol - vol_check)>1e-5)){
			std::cout<<"(abs(l_ - l_check)<1e-5) || (abs(vol - vol_check)<1e-5)"<<std::endl;
			std::cout<<ot<<" "<<(abs(l_ - l_check)<1e-5)<<" "<< (abs(vol - vol_check)<1e-5)<<std::endl;
			std::cout<<l_<<" "<<l_check<<" "<<vol<<" "<<vol_check<<std::endl;
			assert(false);
		}
		double newl = org->orgVolume2Length(vol + delta_volOrg);
		double newl_max = org->orgVolume2Length(vol + delta_volOrgmax);
		double orgGr = newl - l_;
		double delta_lmax = newl_max - l_;
		//if(ot==3){std::cout<<"deltavolorg "<<delta_volOrg<<" "<<delta_volOrgmax<<std::endl;}
		//if(ot==3){std::cout<<"vol increase "<<vol<<" "<<l_<<" "<<newl<<" "<<newl_max<<" "<<orgGr<<" "<<delta_lmax<<std::endl;}
		if((orgGr < 0.)&&(orgGr > -1e-10)){orgGr=0.;}
			assert((orgGr >= 0.) &&"negative orgGr");
		if((delta_lmax < 0.)&&(delta_lmax > -1e-10)){delta_lmax=0.;}
			assert((delta_lmax >= 0.) &&"negative delta_lmax");	
		//orgGr += delta_l; //[orgID]
		if(ot == 2){//each organ type has it s own growth function (and thus CW_Gr map)
			cWGrRoot.insert(std::pair<int, double>(orgID, orgGr));
		}
		if(ot == 3){
			cWGrStem.insert(std::pair<int, double>(orgID, orgGr));
		}
		if(ot == 4){
			cWGrLeaf.insert(std::pair<int, double>(orgID,orgGr));
		}
		// std::cout<<"aliveend1 "<<orgGr<<std::endl;
		delta_ls_org[orgID2] += orgGr;
		// std::cout<<"aliveend2 "<<std::endl;
		delta_ls_org_i[orgID2] = orgGr;
		//std::cout<<"orgID "<<orgID<<" "<<orgID2<<" "<< orgGr<<std::endl;
		//std::cout<<"orgID "<<delta_ls_org[orgID2]<<" "<< delta_ls_org_i[orgID2]<<std::endl;
		// std::cout<<"aliveend3 "<<delta_lmax<<std::endl;
		delta_ls_org_max[orgID2] += delta_lmax;
		// std::cout<<"aliveend4 "<<std::endl;
		delta_ls_org_imax[orgID2] = delta_lmax;
		// std::cout<<"aliveend5 "<<std::endl;
		orgID2 ++;
		
	}
	
	for (auto orp : plant->getOrganRandomParameter(2)) { // each maps is copied for each sub type; or (todo), we could pass a pointer, and keep maps in this class
		orp->f_gf->CW_Gr = cWGrRoot;
	}
	for (auto orp : plant->getOrganRandomParameter(3)) {
		orp->f_gf->CW_Gr = cWGrStem;
	}
	for (auto orp : plant->getOrganRandomParameter(4)) {
		orp->f_gf->CW_Gr = cWGrLeaf;
	}
		//std::cout<<"aliveend6 "<<std::endl;
	//Q_Grmax.display();
}

		