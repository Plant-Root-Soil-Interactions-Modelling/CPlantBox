/*
* PiafMunch (v.2) -- Implementing and Solving the Munch model of phloem sap flow in higher plants
*
* Copyright (C) 2004-2019 INRA
*
* Author: A. Lacointe, UMR PIAF, Clermont-Ferrand, France
*
* File: odepack.cpp
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

#include <PiafMunch/odepack.h>
#include <time.h>

// La fonction suivante "habille" un Fortran_vector en N_Vector, i.e. "cree" un N_Vector
// qui PARTAGE LES MEMES DONNEES (atributs v_ et NV_DATA_S, resp.)  que le Fortran_vector d'origine.
// Ceci permet d'appliquer "in-place" a un Fortran_vector les fonctions appliquables a un N_Vector :
//
N_Vector InPlace_NVector(Fortran_vector &v) {
         double* v_ = InPlace_Array(v);
		return(N_VMake_Serial((sunindextype)v_[0], (realtype*)(v_ + 1)));
}

int check_flag(void *flagvalue, string funcname_, int opt);   // utilise en 'verbose' dans cvode
const time_t current = time(NULL) ;
Fortran_vector y ;
void* cvode_mem ;        // espace de travail du solveur; eventuellement utilise par la fonction aux()
void* arkode_mem;        // espace de travail du solveur; eventuellement utilise par la fonction aux()

void (*ff)(double t ,double* y, double* ydot) ;

inline int ffff(realtype t, N_Vector yy, N_Vector yydot, void *f_data) {
// Forme de f compatible avec cvode: noter decalage des indices entre double* y et NV_DATA_S(N_Vector y)
    ff(t, NV_DATA_S(yy) - 1, NV_DATA_S(yydot) - 1) ;
	return 0;
}

// les 2 fonctions suivantes sont a reformuler avec fonctions 'synonymes' a traiter comme f
// pour devenir independantes du nom d'appel ('rootfind'...)
extern void rootfind(double t, double* y, double* g);       // (optionnel) definit d'eventuelles equations g(t,y)=0 a resoudre

inline int gg(realtype t, N_Vector yy, realtype *g, void *g_data) {  // Idem pour rootfind / gg :
    rootfind(t, NV_DATA_S(yy) - 1, g-1) ; // retourne le nb d'equations dans la var. n
    return 0;
}

SUNLinearSolver LS(NULL);
/* Other Constants pour calcul KLU_DQ_Jac : */
#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

static int Jac_(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
	CVodeMem cv_mem;
	sunindextype *colptrs = SUNSparseMatrix_IndexPointers(J);
	sunindextype *rowvals = SUNSparseMatrix_IndexValues(J);
	realtype *data = SUNSparseMatrix_Data(J);
	N_Vector  ftemp; // stockera temporairement  f(y + inc) = f(y) + (df/dy[j]) * dy[j] = valeur de f quand y[j] est decale de inc=dy[j]
	realtype fnorm, minInc, inc, yjsaved, srur, conj;
	realtype *y_data, *fy_data, *ftemp_data, *ewt_data, *cns_data;
	sunindextype i, j, N, NNZ, NNZ0, retval, npnz, ntnz = 0; // ntnz sera le nombre de NZ effectivement trouve, y compris les elements diagonaux meme si nuls
	realtype J_ij ; bool nnz_str_changed(false) ;
	cv_mem = (CVodeMem)cvode_mem;
	SUNMatZero(J);
	y_data = N_VGetArrayPointer(y);
	retval = 0;
	npnz = NNZ = NNZ0 = SM_NNZ_S(J);// npnz de Sparse_matrix
	// access matrix dimension  et  verif. coherence J, y, fy :
	N = SM_COLUMNS_S(J);
	assert(SM_ROWS_S(J) == N);
	assert(NV_LENGTH_S(y) == N);
	assert(NV_LENGTH_S(fy) == N);
	sunindextype** KLU_Ai = NULL; realtype ** KLU_Ax = NULL; // serviront eventuellement de tampon de stockage si nnz_eff depasse NNZ
	sunindextype n_KLU_ptrs_1(-1), n_KLU_ptrs_2(N - 1); //   n_KLU_ptrs_1 =  nombre de tampons de taille N deja utilises  - 1
	/* Rename work vector for readibility */
	ftemp = tmp1;
	/* Obtain pointers to the data for ewt, y */
	ftemp_data = N_VGetArrayPointer(ftemp);
	ewt_data = N_VGetArrayPointer(cv_mem->cv_ewt);
	y_data = N_VGetArrayPointer(y);
	fy_data = N_VGetArrayPointer(fy);
	if (cv_mem->cv_constraints != NULL)
		cns_data = N_VGetArrayPointer(cv_mem->cv_constraints);
	/* Set minimum increment based on uround and norm of f */
	srur = SUNRsqrt(cv_mem->cv_uround);
	fnorm = N_VWrmsNorm(fy, cv_mem->cv_ewt);
	minInc = (fnorm != ZERO) ?
		(MIN_INC_MULT * SUNRabs(cv_mem->cv_h) * cv_mem->cv_uround * N * fnorm) : ONE;
	for (j = 0; j < N; j++) { // Generate the jth col of J(tn,y)
		if (colptrs[j] != ntnz) {
			colptrs[j] = ntnz;
			nnz_str_changed = true;
		}
		yjsaved = y_data[j];
		inc = SUNMAX(srur*SUNRabs(yjsaved), minInc / ewt_data[j]); // inc = dy[j]
		// Adjust sign(inc) if y_j has an inequality constraint :
		if (cv_mem->cv_constraints != NULL) {
			conj = cns_data[j];
			if (SUNRabs(conj) == ONE) { if ((yjsaved + inc)*conj < ZERO)  inc = -inc; }
			else if (SUNRabs(conj) == TWO) { if ((yjsaved + inc)*conj <= ZERO) inc = -inc; }
		}
		y_data[j] += inc; // = y[j] + dy[j]
		retval = cv_mem->cv_f(t, y, ftemp, cv_mem->cv_user_data); // ftemp = y_dot(t, y + inc)     (y+inc ne differe de y qu par sa jeme composante)
		if (retval != 0) break;
		y_data[j] = yjsaved; // restaure y
		for (i = 0; i < N; i++) { // ligne i de la colonne j
			if ((ftemp_data[i] != fy_data[i]) || (i == j)) {
				if (ftemp_data[i] != fy_data[i])   J_ij = (ftemp_data[i] - fy_data[i]) / inc;   else  J_ij = ZERO;
				if (ntnz >= NNZ0) { // avec celui-ci en plus, on va depasser l'espace reserve pour J
					if (ntnz == NNZ0) cout << "ntnz = "  << ntnz << " va depasser NNZ = " << NNZ0 << " => creation des tampons :" << endl;
					if (!KLU_Ai) {
						assert(!KLU_Ax);
						KLU_Ai = new sunindextype*[N]; KLU_Ax = new realtype*[N];
//						cout << "KLU_Ai , _Ax = new ...type*[N]" << endl ;
					}
					n_KLU_ptrs_2++;
					if (n_KLU_ptrs_2 == N) {
						n_KLU_ptrs_2 = 0; // de la ligne KLU_Ai[1 + n_KLU_ptrs_1], qui n'est pas encore utilisee, et donc n'a pas encore ete initialisee :
						n_KLU_ptrs_1++; assert(n_KLU_ptrs_1 < N);
						KLU_Ai[n_KLU_ptrs_1] = new sunindextype[N]; KLU_Ax[n_KLU_ptrs_1] = new realtype[N];
	//					cout << "KLU_Ai[...ptrs1] , _Ax[...ptrs1] = new ...type[N]" << endl;
					}
					KLU_Ai[n_KLU_ptrs_1][n_KLU_ptrs_2] = i; KLU_Ax[n_KLU_ptrs_1][n_KLU_ptrs_2] = J_ij;
				}
				else {
					if (rowvals[ntnz] != i) {
						rowvals[ntnz] = i;
						nnz_str_changed = true;
					}
					data[ntnz] = J_ij;
				}
				ntnz++;
			}
		}
	}
	if (ntnz > NNZ0) { // on a du stocker des valeurs dans les tampons KLU_A.[0 .. n_KLU_ptrs_1 ][ ]
//		cout << "ntnz > NNZ : reallocation matrice J ; utilisation et destruction des tampons :" << endl;
		npnz = ntnz / N ; // (provisoirement) nombre de lignes  pleines de longueur N, pouvant contenir tous les ntnz 'zeros'
		assert(npnz * N <= ntnz) ;
		if (npnz * N < ntnz) npnz ++ ; // on dimensionne le nouveau 'NNZ' = npnz  en multiples entiers de N :
		npnz = N * npnz ; assert(npnz >= ntnz) ; // le vrai nouveau npnz a reserver
		SUNSparseMatrix_Reallocate(J, npnz) ;	// nouveau NNZ officiel = SM_NNZ_S(J) = npnz de Sparse_matrix ; conserve les valeurs d'index < NNZ
		rowvals = SUNSparseMatrix_IndexValues(J);
		data = SUNSparseMatrix_Data(J);
		colptrs = SUNSparseMatrix_IndexPointers(J); // n'a pas du changer...
		for (i = 0; i < n_KLU_ptrs_1; i++) { // recopie les lignes KLU_A.[ ] completes a coup sur (s'il y en a, i.e. si n_KLU_ptrs_1 >= 1)
			for (j = 0; j < N; j++) {
				rowvals[NNZ] = KLU_Ai[i][j]; data[NNZ] = KLU_Ax[i][j];
				NNZ++;
			}
			delete[] KLU_Ai[i]; delete[] KLU_Ax[i];
//			cout << "delete[] KLU_Ai[i] , _Ax[i] ; " ;
		}
		for (j = 0; j <= n_KLU_ptrs_2; j++) { // la derniere ligne, eventuellement incomplete
			rowvals[NNZ] = KLU_Ai[n_KLU_ptrs_1][j]; data[NNZ] = KLU_Ax[n_KLU_ptrs_1][j];
			NNZ++;
		}
		delete[] KLU_Ai[n_KLU_ptrs_1]; delete[] KLU_Ax[n_KLU_ptrs_1];
//		cout << "delete[] KLU_Ai[...ptrs1] , _Ax[...ptrs1] ; " ;
		delete[] KLU_Ai; delete[] KLU_Ax;
//		cout << "delete[] KLU_Ai , _Ax" << endl;
		assert(NNZ == ntnz); // attention, ce n'est plus le NNZ officiel ! -- qui reste stocke en NNZ0
	}
	if (colptrs[N] != ntnz) {
		colptrs[N] = ntnz;
		nnz_str_changed = true;
	}
	if (nnz_str_changed) {
		if (LS) { // le solveur a deja ete initialise par cvode_direct( ), donc une reinit. PARTIAL suffit :
			SUNLinSol_KLUReInit(LS, J, ntnz, SUNKLU_REINIT_PARTIAL);
		}
	}
	return(0);
}


int cvode_direct(void(*f)(double,double*,double*), Fortran_vector& y, Fortran_vector &T, void(*aux)(double,double*), Fortran_vector& atol, Fortran_vector& rtol, int solver,
	     int nbVar_dot, Fortran_vector** Var_primitive, Fortran_vector** Var_dot, bool verbose, bool STALD, void(*rootfind)(double, double*, double*), int nrootfns, int mu, int ml) {
	ff = f ;
  int itol = 1 ;
/* 	itol = 1 (= CV_SS) : atol and rtol tous 2 scalaires ;
		itol = 2 (= CV_SV) : atol vecteur, rtol scalaire ;
		itol = 3 (= CV_WF) : ewt[i] defini par une fonction utilisateur (non encore implemente ; cf. doc CVODE p.29) */
  if (atol.size() > 1) itol ++ ;
  if (rtol.size() > 1) itol += 2 ;
  if (itol == 3) {_LogMessage("Erreur itol = CV_WF : cas non encore traite par cvode") ; return -1;}
  if (itol == 4) {_LogMessage("Erreur itol = 4 (atol et rtol tous deux vecteurs): non traite par cvode") ; return -1 ;}
  int neq = y.size();							 // taille du probleme = nb d'equations = nb de variables y[i]
  int nbt = T.size() ;							 // nombre de temps ou l'on veut la solution, y compris t0 = T[1]
  int i, j, flag, flagr ;
  double t, tout, t_sav, dt ;
  N_Vector abstol(InPlace_NVector(atol)) ;
  N_Vector yy(InPlace_NVector(y));              // "habille" le Fortran_vector y en N_Vector sans occuper de memoire suppl.
  SUNMatrix A;
  double* y_ = InPlace_Array(y);	// y_ = y.v_
	t = t_sav = T[1] ;
   Fortran_vector** Var_primitive_sav ;
   if (aux != NULL) aux(T[1], y_) ; // instant t0 = T[1]
  Update_Output(true) ;
  if (nbVar_dot) {
	  assert(Var_primitive) ;
	  assert(Var_dot) ;
	  for (j = 0 ; j < nbVar_dot ; j ++) {
		  assert(Var_primitive[j]) ;
		  assert(Var_dot[j]) ;
	  }
	  assert(Var_primitive[nbVar_dot] == NULL) ; assert(Var_dot[nbVar_dot] == NULL) ; // previendra l'utilisateur si on a declare moins de derivees que l'on n'en a dimensionnees
	  Var_primitive_sav = new Fortran_vector*[nbVar_dot] ;
	  for (j = 0 ; j < nbVar_dot ; j ++) 	  Var_primitive_sav[j] = new Fortran_vector(Var_primitive[j]->size()) ;
  }
  cvode_mem = CVodeCreate(CV_BDF);    // init. le solveur avec la methode d'integr. BDF (et la resol. de type Newton)
    if (check_flag(cvode_mem, "CVodeCreate", 0)) {_LogMessage("erreur CVodeCreate") ; return -1 ; }
  flag = CVodeInit(cvode_mem, ffff, T[1], yy); // alloue l'espace de travail du solveur
    if (check_flag(&flag, "CVodeInit", 1)) {_LogMessage("erreur CVodeInit") ; return -1 ; }

/* Call CVodeSVtolerances to specify the scalar relative tolerance
* and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, rtol[1], abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

  if (STALD) {                // algorithme de detection/correction de stabilite aux ordres > 2 (integration par methode BDF)
      flag = CVodeSetStabLimDet(cvode_mem, SUNTRUE);    // FALSE par defaut
        if (check_flag(&flag, "CVodeSetStabLimDet", 1)) {_LogMessage("erreur CVodeSetStabLimDet") ; return -1 ; }
  }
  int* rootsfound = NULL ;
  if (rootfind != NULL) {
      if(nrootfns < 1) {_LogMessage("erreur RootFind : nrootfns =< 0 eq. a resoudre !") ; return -1 ; }
      rootsfound = new int[nrootfns] ;
	  flag = CVodeRootInit(cvode_mem, nrootfns, gg);
	  if (check_flag(&flag, "CVodeRootInit", 1)) {_LogMessage("erreur CVodeRootInit") ; return -1 ; }
  }
  if (solver == DIAG) {
	  flag = CVDiag(cvode_mem);
	  if (check_flag(&flag, "CVDiag", 1)) { _LogMessage("erreur CVDiag"); return -1; }
  } else {
	  if (solver == DENSE) {
		  /* Create dense SUNMatrix for use in linear solves */
		  A = SUNDenseMatrix(neq, neq);
		  if (!(SM_ROWS_S(A) == neq)) { _LogMessage("erreur SunDenseMatrix"); return -1; }
		  /* Create dense SUNLinearSolver object for use by CVode */
		  LS = SUNLinSol_Dense(yy, A);
		  if (!LS) { _LogMessage("erreur SUNLinSol_Dense"); return -1; }
	  } else {
		  if (solver == BAND) {
			  A = SUNBandMatrix(neq, mu, ml);
			  if (!(SM_ROWS_S(A) == neq)) { _LogMessage("erreur SunBandMatrix"); return -1; }
			  /* Create banded SUNLinearSolver object for use by CVode */
			  LS = SUNLinSol_Band(yy, A);
			  if (!LS) { _LogMessage("erreur SUNLinSol_Band"); return -1; }
		  } else {
				if (solver == KLU) {// cree les 3 espaces memoire (KLU::Ax, Ai, Ap) dimensionnes ; on pourra les redimensionner ulterieurement...
					A = SUNSparseMatrix(neq, neq, neq, CSC_MAT);// 3eme arg : nnz = neq = le + petit a priori (la diag.)
					if (check_flag((void *)A, "SUNSparseMatrix", 0)) return(1);
					/* Create KLU solver object for use by CVode */
					LS = SUNLinSol_KLU(yy, A);
					if (check_flag((void *)LS, "SUNLinSol_KLU", 0)) return(1);
				}
				else { (void)sprintf(message, "%d : nom de solveur inconnu", solver); _LogMessage(message); return -1; }
		  }
	  }
	  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
	  flag = CVodeSetLinearSolver(cvode_mem, LS, A);
	  if (check_flag(&flag, "CVodeSetLinearSolver", 1)) return(1);
	  if (solver == KLU) {
		  flag = CVodeSetJacFn(cvode_mem, Jac_);
		  if (check_flag(&flag, "CVodeSetJacFn", 1)) return(1); // Non-zero = erreur
	  }
  }
  flag = CVodeSetMaxConvFails(cvode_mem, 100) ; // pour prevenir l'erreur de non-convergence (error code = -4), sauf si la convergence n'est effectivement jamais atteinte !
  for (i = 2 ; i <= nbt ; i++) {						// pour chaque instant ou l'on souhaite la solution
	tout=T[i];
	if (verbose) {
			strftime(message, 50, "%H:%M:%S", localtime(&current)) ; cout <<  "at " << message << " :  starting step n#" << i-1 << " (tf = " << tout << ")" << endl ;
			Update_Output(i == 2) ;
	}
	while (t < tout) {
		flag = -9999 ; // bidon
		while (flag != CV_SUCCESS) {
			  flag = CVode(cvode_mem, tout, yy, &t, CV_ONE_STEP); // appel au solveur en mode 'CV_ONE_STEP' // attention : cette instruction MAJ  t, qui peut donc maintenant etre > tout !!
			  if (flag == CV_ROOT_RETURN) {
				 flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
				   if (check_flag(&flagr, "CVodeGetRootInfo", 1)) { _LogMessage(" Erreur CVodeGetRootInfo") ;
						//strftime(message, 50, "%H:%M:%S", localtime(&current)) ; cout <<  "at " << message << " :  exiting solver" << endl ; 
						Update_Output() ; return i ;
				   }
				 for (j = 0 ; j < nrootfns ; j++) {
					 if (rootsfound[j] == 1) {
						(void)sprintf(message, "Eq.[%d] : Root Found at t = %g :", j+1, t) ; _LogMessage(message) ;
						Update_Output() ;
					 }
				 }
				 aux(t, y_) ;
			  }
			  else {
				   if (flag != CV_SUCCESS) {
				  //    if (check_flag(&flag, "CVode", 1)) break ; // routine a  reprendre...
					   (void)sprintf(message, "error-flag CVode = %d", flag) ; _LogMessage(message) ;
						//strftime(message, 50, "%H:-%M:%S", localtime(&current)) ; cout <<  "at " << message << " :  exiting solver" << endl ; 
						Update_Output(true) ;
					   return i ;
				  }
			  }
		}
		if (nbVar_dot) {
			dt = t - t_sav;
			for (j = 0; j < nbVar_dot; j++) {
				Var_dot[j]->set((*(Var_primitive[j]) - *(Var_primitive_sav[j])) / dt);
				Var_primitive_sav[j]->set(*(Var_primitive[j]));
			}
			if (t < tout)  t_sav = t;
		}
	} // fin boucle  while(t < tout) : les 8 lignes suivantes sont donc executees si (t > tout) :
		if (aux != NULL) {		// calculs auxiliaires, par ex. sauvegarder var. intermediaires
								// a ce stade t_curr est juste superieur a tout : on interpole y a y_out = y(tout) :
			flag = CVodeGetDky(cvode_mem, tout, 0, yy);
			check_flag(&flag, "CVodeGetDky", 1);
			aux(T[i], y_);
		}
		t_sav = t;
	} // fin boucle i (T[i])
  if (verbose) {  // Print some final statistics :
      long int nst, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
      flag = CVodeGetNumSteps(cvode_mem, &nst);
            check_flag(&flag, "CVodeGetNumSteps", 1);
      flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
            check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
      flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
            check_flag(&flag, "CVodeGetNumErrTestFails", 1);
      flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
           check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
      flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
           check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
      flag = CVodeGetNumJacEvals(cvode_mem, &nje);
           check_flag(&flag, "CVodeGetNumJacEvals", 1);
      flag = CVodeGetNumRhsEvals(cvode_mem, &nfeLS);
           check_flag(&flag, "CVodeGetNumRhsEvals", 1);
      flag = CVodeGetNumGEvals(cvode_mem, &nge);
           check_flag(&flag, "CVodeGetNumGEvals", 1);
      _LogMessage("\nFinal Statistics:");
      (void)sprintf(message, "nst = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n", nst, nsetups, nfeLS, nje);
           _LogMessage(message) ;
      (void)sprintf(message, "nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n", nni, ncfn, netf, nge); _LogMessage(message) ;
  }
  /* Free integrator memory : */
  CVodeFree(&cvode_mem);
  N_VDestroy_Serial(yy) ; N_VDestroy_Serial(abstol) ; // heureusement, ne desallouent pas leurs 'double* NV_DATA_S()' car issus de 'N_VMake_Serial'
  if (rootfind != NULL) delete [] rootsfound ;
  if (nbVar_dot) {
	  for (j = 0 ; j < nbVar_dot ; j ++) 	 delete Var_primitive_sav[j] ;
	delete[] Var_primitive_sav ;
  }
//  _strtime_s(message, 100); cout << "at " << message << " :  exiting solver" << endl ; Update_Output() ;
   //strftime(message, 50, "%H:%M:%S", localtime(&current)) ; cout <<  "at " << message << " :  exiting solver" << endl ; 
   Update_Output() ;
	return 0 ;
}


int cvode_spils(void(*f)(double, double*, double*), Fortran_vector &y, Fortran_vector &T, void(*aux)(double, double*), Fortran_vector& atol, Fortran_vector& rtol,
	int solver, int GSType, int prectype, int nbVar_dot, Fortran_vector** Var_primitive, Fortran_vector** Var_dot,
	bool verbose, bool STALD, void(*rootfind)(double, double*, double*), int nrootfns, int mu, int ml, int maxl) {
	ff = f;
	int itol = 1;
	// 		itol = 1 (= CV_SS) : atol and rtol tous 2 scalaires ;
	//		itol = 2 (= CV_SV) : atol vecteur, rtol scalaire ;
	//		itol = 3 (= CV_WF) : ewt[i] defini par une fonction utilisateur (non encore implemente ; cf. doc CVODE p.29)
	if (atol.size() > 1) itol++;
	if (rtol.size() > 1) itol += 2;
	if (itol == 3) { _LogMessage("Erreur itol = CV_WF : cas non encore traite par cvode"); return -1; }
	if (itol == 4) { _LogMessage("Erreur itol = 4 (atol et rtol tous deux vecteurs): non traite par cvode"); return -1; }
	int neq = y.size();							 // taille du probleme = nb d'equations = nb de variables y[i]
	int nbt = T.size();							 // nombre de temps ou l'on veut la solution, y compris t0 = T[1]
	int i, j, flag, flagr;
	double t, tout, t_sav, dt;
	N_Vector yy(InPlace_NVector(y));              // "habille" le Fortran_vector y en N_Vector sans occuper de memoire suppl.
	N_Vector abstol(InPlace_NVector(atol));
	double* y_ = NV_DATA_S(yy) - 1;
	t = t_sav = T[1];
	Fortran_vector** Var_primitive_sav;
	if (aux != NULL) aux(T[1], y_); // instant t0 = T[1]
	Update_Output(true);
	if (nbVar_dot) {
		assert(Var_primitive);
		assert(Var_dot);
		for (j = 0; j < nbVar_dot; j++) {
			assert(Var_primitive[j]);
			assert(Var_dot[j]);
		}
		assert(Var_primitive[nbVar_dot] == NULL); assert(Var_dot[nbVar_dot] == NULL); // previendra l'utilisateur si on a declare moins de derivees que l'on n'en a dimensionnees
		Var_primitive_sav = new Fortran_vector*[nbVar_dot];
		for (j = 0; j < nbVar_dot; j++) 	  Var_primitive_sav[j] = new Fortran_vector(Var_primitive[j]->size());
	}
	cvode_mem = CVodeCreate(CV_BDF);    // init. le solveur avec la methode d'integr. BDF et la resol. de type Newton
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) { _LogMessage("erreur CVodeCreate"); return -1; }

	flag = CVodeInit(cvode_mem, ffff, T[1], yy); // alloue l'espace de travail du solveur
	if (check_flag(&flag, "CVodeInit", 1)) { _LogMessage("erreur CVodeInit"); return -1; }

	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	* and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, rtol[1], abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1))
	{
		std::cout<<"check_flag(&flag, CVodeSVtolerances, 1) "<< flag <<std::endl;
		return(1);
	}

	if (STALD) {                // algorithme de detection/correction de stabilite aux ordres > 2 (integration par methode BDF)
		flag = CVodeSetStabLimDet(cvode_mem, SUNTRUE);    // FALSE par defaut
		if (check_flag(&flag, "CVodeSetStabLimDet", 1)) { _LogMessage("erreur CVodeSetStabLimDet"); return -1; }
	}
	int* rootsfound = NULL;
	if (rootfind != NULL) {
		if (nrootfns < 1) { _LogMessage("erreur RootFind : nrootfns =< 0 eq. a resoudre !"); return -1; }
		rootsfound = new int[nrootfns];
		flag = CVodeRootInit(cvode_mem, nrootfns, gg);
		if (check_flag(&flag, "CVodeRootInit", 1)) { _LogMessage("erreur CVodeRootInit"); return -1; }
	}
	if (solver == SPGMR) {
		LS = SUNLinSol_SPGMR(yy, prectype, maxl);
		if (!LS) { _LogMessage("erreur SUNLinSol_spgmr"); return -1; }
		flag = SUNLinSol_SPGMRSetGSType(LS, GSType); // Gram-Schmidt orthogonalisation method
		if (check_flag(&flag, "SUNLinSol_SPGMRSetGSType", 1)) { _LogMessage("erreur SUNLinSol_SPGMRSetGSType"); return -1; }
	}
	else {
		if (solver == SPFGMR) {
			LS = SUNLinSol_SPFGMR(yy, prectype, maxl);
			if (!LS) { _LogMessage("erreur SUNLinSol_spfgmr"); return -1; }
			flag = SUNLinSol_SPFGMRSetGSType(LS, GSType); // Gram-Schmidt orthogonalisation method
			if (check_flag(&flag, "SUNLinSol_SPFGMRSetGSType", 1)) { _LogMessage("erreur SUNLinSol_SPFGMRSetGSType"); return -1; }
		}
		else {
			if (solver == SPBCGS) {
				LS = SUNLinSol_SPBCGS(yy, prectype, maxl);
				if (!LS) { _LogMessage("erreur SUNLinSol_spbcgs"); return -1; }
			}
			else {
				if (solver == SPTFQMR) {
					LS = SUNLinSol_SPTFQMR(yy, prectype, maxl);
					if (!LS) { _LogMessage("erreur SUNLinSol_sptfqmr"); return -1; }
				}
				else {
					if (solver == PCG) {
						LS = SUNLinSol_PCG(yy, prectype, maxl);
						if (!LS) { _LogMessage("erreur SUNLinSol_PCG"); return -1; }
					}
					else { (void)sprintf(message, "%d : nom de solveur inconnu", solver); _LogMessage(message); return -1; }
				}
			}
		}
	}
	flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
	if (check_flag(&flag, "CVodeSetLinearSolver", 1)) return -1;
	if (prectype != PREC_NONE) {
		flag = CVBandPrecInit(cvode_mem, neq, mu, ml);
		cout << "CVBandPrecInit flag = " << flag << endl;
		if (check_flag(&flag, "CVBandPrecInit", 1)) { _LogMessage("erreur CVBandPrecInit"); return -1; }
	}
  flag = CVodeSetMaxConvFails(cvode_mem, 100) ; // pour prevenir l'erreur de non-convergence (error code = -4), sauf si la convergence n'est effectivement jamais atteinte !
  for (i = 2 ; i <= nbt ; i++) {						// pour chaque instant ou l'on souhaite la solution
	tout=T[i];
	if (verbose) {
			strftime(message, 50, "%H:%M:%S", localtime(&current)) ; cout <<  "at " << message << " :  starting step n#" << i-1 << " (tf = " << tout << ")" << endl ;
			Update_Output(i == 2) ;
	}
	while (t < tout) {
		flag = -9999 ; // bidon
		while (flag != CV_SUCCESS) {
			  flag = CVode(cvode_mem, tout, yy, &t, CV_ONE_STEP); // appel au solveur
			  if (flag == CV_ROOT_RETURN) {
				 flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
				   if (check_flag(&flagr, "CVodeGetRootInfo", 1)) { _LogMessage(" Erreur CVodeGetRootInfo") ;
		       			//strftime(message, 50, "%H:%M:%S", localtime(&current)) ; cout <<  "at " << message << " :  exiting solver" << endl ; 
						Update_Output() ; return i ;
				   }
				 for (j = 0 ; j < nrootfns ; j++) {
					 if (rootsfound[j] == 1) {
						(void)sprintf(message, "Eq.[%d] : Root Found at t = %g :", j+1, t) ; _LogMessage(message) ;
						Update_Output() ;
					 }
				 }
				 
				std::cout<<"aux for t < tout"<<std::endl;
				aux(t, y_) ;
			  }
			  else {
				  if (flag != CV_SUCCESS) {
				  //    if (check_flag(&flag, "CVode", 1)) break ; // routine a  reprendre...
					   (void)sprintf(message, "error-flag CVode = %d", flag) ; _LogMessage(message) ;
					   	//strftime(message, 50, "%H:%M:%S", localtime(&current)) ; cout <<  "at " << message << " :  exiting solver" << endl ; 
						Update_Output() ;
					   return i ;
				  }
			  }
		}
		if (nbVar_dot) {
			dt = t - t_sav;
			for (j = 0; j < nbVar_dot; j++) {
				Var_dot[j]->set((*(Var_primitive[j]) - *(Var_primitive_sav[j])) / dt);
				Var_primitive_sav[j]->set(*(Var_primitive[j]));
			}
			if (t < tout)  t_sav = t;
		}
	  } // fin boucle  while(t < tout) : les 8 lignes suivantes sont donc executees si (t > tout) :
	  if (aux != NULL) {		// calculs auxiliaires, par ex. sauvegarder var. intermediaires
		//if (true){//aux != NULL) {		// calculs auxiliaires, par ex. sauvegarder var. intermediaires
								// a ce stade t_curr est juste superieur a tout : on interpole y a y_out = y(tout) :
		  flag = CVodeGetDky(cvode_mem, tout, 0, yy);
		  check_flag(&flag, "CVodeGetDky", 1);
		  aux(T[i], y_);
	  }
	  t_sav = t;
	} // fin boucle i (T[i])
	if (verbose) {  // Print some final statistics :
      long int nst, nfe, nsetups, nni, ncfn, netf, nge;
      flag = CVodeGetNumSteps(cvode_mem, &nst);
            check_flag(&flag, "CVodeGetNumSteps", 1);
      flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
            check_flag(&flag, "CVodeGetNumRhsEvals", 1);
      flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
            check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
      flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
            check_flag(&flag, "CVodeGetNumErrTestFails", 1);
      flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
           check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
      flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
           check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
      flag = CVodeGetNumGEvals(cvode_mem, &nge);
           check_flag(&flag, "CVodeGetNumGEvals", 1);
      _LogMessage("\nFinal Statistics:");
      (void)sprintf(message, "nst (num steps) = %-6ld nfe  (num call to f)= %-6ld nsetups (call to lin solver setup func)= %-6ld\n", nst, nfe, nsetups);
           _LogMessage(message) ;
      (void)sprintf(message, "nni (iter of nonlinear solver) = %-6ld ncfn (non linsolver conv fail)= %-6ld netf (num err test fail) = %-6ld nge (call to root function) = %ld\n \n", nni, ncfn, netf, nge); _LogMessage(message) ;
  }
  // Free integrator memory :
  CVodeFree(&cvode_mem);
  N_VDestroy_Serial(yy) ; N_VDestroy_Serial(abstol) ; // heureusement, ne desallouent pas leurs 'double* NV_DATA_S()' car issus de 'N_VMake_Serial'
  if (rootfind != NULL) delete [] rootsfound ;
  if (nbVar_dot) {
	  for (j = 0 ; j < nbVar_dot ; j ++) 	 delete Var_primitive_sav[j] ;
	delete[] Var_primitive_sav ;
  }
//  _strtime_s(message, 100); cout << "at " << message << " :  exiting solver" << endl ; Update_Output();
	// strftime(message, 50, "%H:%M:%S", localtime(&current)) ; cout <<  "at " << message << " :  exiting solver" << endl ; 
	Update_Output();

	return 0 ;
}
/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */
	int check_flag(void *flagvalue, string funcname_, int opt) {
  int *errflag;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  const char * funcname = funcname_.c_str();
  if (opt == 0 && flagvalue == NULL) {
    (void)sprintf(message, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname); _LogMessage(message) ;
//    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return(1);
  }
  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      (void)sprintf(message, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag); _LogMessage(message) ;
  //    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
      return(1);
    }
  }
  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    (void)sprintf(message, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname); _LogMessage(message) ;
//    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return(1);
  }
  return(0);
}


int arkode(void(*f)(double, double*, double*), Fortran_vector& y, Fortran_vector &T, void(*aux)(double, double*), Fortran_vector& atol, Fortran_vector& rtol,
	int nbVar_dot, Fortran_vector** Var_primitive, Fortran_vector** Var_dot, bool verbose,	void(*rootfind)(double, double*, double*), int nrootfns) {
	ff = f;
	int itol = 1;
	/* 	itol = 1 (= CV_SS) : atol and rtol tous 2 scalaires ;
	itol = 2 (= CV_SV) : atol vecteur, rtol scalaire ;
	itol = 3 (= CV_WF) : ewt[i] defini par une fonction utilisateur (non encore implemente ; cf. doc CVODE p.29) */
	if (atol.size() > 1) itol++;
	if (rtol.size() > 1) itol += 2;
	if (itol == 3) { _LogMessage("Erreur itol = CV_WF : cas non encore traite par cvode"); return -1; }
	if (itol == 4) { _LogMessage("Erreur itol = 4 (atol et rtol tous deux vecteurs): non traite par cvode"); return -1; }
	//int neq = y.size();							 // taille du probleme = nb d'equations = nb de variables y[i]
	int nbt = T.size();							 // nombre de temps ou l'on veut la solution, y compris t0 = T[1]
	int i, j, flag, flagr;
	double t, tout, t_sav, dt;
	N_Vector abstol(InPlace_NVector(atol));
	N_Vector yy(InPlace_NVector(y));              // "habille" le Fortran_vector y en N_Vector sans occuper de memoire suppl.
	double* y_ = InPlace_Array(y);	// y_ = y.v_
	t = t_sav = T[1];
	Fortran_vector** Var_primitive_sav;
	if (aux != NULL) aux(T[1], y_); // instant t0 = T[1]
	Update_Output(true);
	if (nbVar_dot) {
		assert(Var_primitive);
		assert(Var_dot);
		for (j = 0; j < nbVar_dot; j++) {
			assert(Var_primitive[j]);
			assert(Var_dot[j]);
		}
		assert(Var_primitive[nbVar_dot] == NULL); assert(Var_dot[nbVar_dot] == NULL); // previendra l'utilisateur si on a declare moins de derivees que l'on n'en a dimensionnees
		Var_primitive_sav = new Fortran_vector*[nbVar_dot];
		for (j = 0; j < nbVar_dot; j++) 	  Var_primitive_sav[j] = new Fortran_vector(Var_primitive[j]->size());
	}
	arkode_mem = ERKStepCreate(ffff, T[1], yy);    // init. le solveur
	if (check_flag(arkode_mem, "ERKStepCreate", 0)) { _LogMessage("erreur ERKStepCreate"); return -1; }

	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	* and vector absolute tolerances */
	flag = ERKStepSVtolerances(arkode_mem, rtol[1], abstol);
	if (check_flag(&flag, "ERKStepSVtolerances", 1)) return(1);
	int* rootsfound = NULL;
	if (rootfind != NULL) {
		if (nrootfns < 1) { _LogMessage("erreur RootFind : nrootfns =< 0 eq. a resoudre !"); return -1; }
		rootsfound = new int[nrootfns];
		flag = ERKStepRootInit(arkode_mem, nrootfns, gg);
		if (check_flag(&flag, "ERKStepRootInit", 1)) { _LogMessage("erreur ERKStepRootInit"); return -1; }
	}
	for (i = 2; i <= nbt; i++) {						// pour chaque instant ou l'on souhaite la solution
		tout = T[i];
		if (verbose) {
			strftime(message, 50, "%H:%M:%S", localtime(&current)) ; cout <<  "at " << message << " :  starting step n#" << i - 1 << " (tf = " << tout << ")" << endl;
			Update_Output(i == 2);
		}
		while (t < tout) {
			flag = -9999; // bidon
			while (flag != ARK_SUCCESS) {
				flag = ERKStepEvolve(arkode_mem, tout, yy, &t, ARK_ONE_STEP); // appel au solveur en mode 'ARK_ONE_STEP' // attention : cette instruction MAJ  t, qui peut donc maintenant etre > tout !!
				if (flag == ARK_ROOT_RETURN) {
					flagr = ERKStepGetRootInfo(arkode_mem, rootsfound);
					if (check_flag(&flagr, "ERKStepGetRootInfo", 1)) {
						_LogMessage(" Erreur ERKStepGetRootInfo");
				        //strftime(message, 50, "%H:%M:%S", localtime(&current)) ; cout <<  "at " << message << " :  exiting solver" << endl; 
						Update_Output(); return i;
					}
					for (j = 0; j < nrootfns; j++) {
						if (rootsfound[j] == 1) {
							(void)sprintf(message, "Eq.[%d] : Root Found at t = %g :", j + 1, t); _LogMessage(message);
							Update_Output();
						}
					}
					aux(t, y_);
				}
				else {
					if (flag != ARK_SUCCESS) {
						//    if (check_flag(&flag, "CVode", 1)) break ; // routine a  reprendre...
						(void)sprintf(message, "error-flag Arkode = %d", flag); _LogMessage(message);
						//strftime(message, 50, "%H:%M:%S", localtime(&current)) ; cout <<  "at " << message << " :  exiting solver" << endl; 
						Update_Output(true);
						return i;
					}
				}
			}
			if (nbVar_dot > 0) {
				dt = t - t_sav;
				for (j = 0; j < nbVar_dot; j++) {
					Var_dot[j]->set((*(Var_primitive[j]) - *(Var_primitive_sav[j])) / dt);
					Var_primitive_sav[j]->set(*(Var_primitive[j]));
				}
				if (t < tout) 	t_sav = t;
			}
		} // fin boucle  while(t < tout) : les 8 lignes suivantes sont donc executees si (t > tout) :
		if (aux != NULL) {		// calculs auxiliaires, par ex. sauvegarder var. intermediaires
								// a ce stade t_curr est juste superieur a tout : on interpole y a y_out = y(tout) :
			flag = ERKStepGetDky(arkode_mem, tout, 0, yy);
			check_flag(&flag, "ERKStepGetDky", 1);
			aux(T[i], y_);
		}
		t_sav = t;
	} // fin boucle i (T[i])
	if (verbose) {  // Print some final statistics :
		long int nst, nfe, netf;
		flag = ERKStepGetNumSteps(arkode_mem, &nst);
		check_flag(&flag, "ERKStepGetNumSteps", 1);
		flag = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
		check_flag(&flag, "ERKStepGetNumRhsEvals", 1);
		flag = ERKStepGetNumErrTestFails(arkode_mem, &netf);
		check_flag(&flag, "ERKStepGetNumErrTestFails", 1);
		_LogMessage("\nFinal Statistics:");
		(void)sprintf(message, "nst = %-6ld nfe  = %-6ld netf = %-6ld\n", nst, nfe, netf);
		_LogMessage(message);
	}
	/* Free integrator memory : */
	ERKStepFree(&arkode_mem);
	N_VDestroy_Serial(yy); N_VDestroy_Serial(abstol); // heureusement, ne desallouent pas leurs 'double* NV_DATA_S()' car issus de 'N_VMake_Serial'
	if (rootfind != NULL) delete[] rootsfound;
	if (nbVar_dot) {
		for (j = 0; j < nbVar_dot; j++) 	 delete Var_primitive_sav[j];
		delete[] Var_primitive_sav;
	}
//	_strtime_s(message, 100); cout << "at " << message << " :  exiting solver" << endl; Update_Output();
	//strftime(message, 50, "%H:%M:%S", localtime(&current)) ; cout <<  "at " << message << " :  exiting solver" << endl; 
	Update_Output();
	return 0;
}


