/*
* PiafMunch (v.2) -- Implementing and Solving the Munch model of phloem sap flow in higher plants
*
* Copyright (C) 2004-2019 INRA
*
* Author: A. Lacointe, UMR PIAF, Clermont-Ferrand, France
*
* File: PM_KLU.cpp
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

# include <float.h> // pour DBL_EPSILON


#include <PiafMunch/PM_arrays.h>


extern double r ;
char KLU_message[1000] ;

void KLU_LogMessage(const char* message) {
	_LogMessage(message) ;
}

#include "klu.h"
#include <nvector/nvector_serial.h>     /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
#include <sundials/sundials_types.h>    /* defs. of realtype, sunindextype      */
#include <cvode/cvode_impl.h>
#include <sundials/sundials_math.h>

double RCOND_LIM = pow(DBL_EPSILON, 2.0 / 3.0);
double CONDEST_LIM = 1. / RCOND_LIM;

int set_KLU_solve_(double*x, int p, const Sparse_matrix &S, double*b, int** ipiv_ptr, void** TM_ptr, bool full_check = false) ; // resoud   Sx = b    en retournant  x =  inv(S) * b

void KLU_allocate(const Sparse_matrix &S, int*** ipiv_ptr_ptr, void*** TM_ptr_ptr) { // cree les espaces de travail appropries pour la matrice S
	assert((!(*ipiv_ptr_ptr)) && (!(*TM_ptr_ptr)));
	int i, r, n = S.n_, m = S.m_, ntnz = S.ntnz_;
	assert(m == n);
	double* Ax = S.v_; int* ni_ = S.ij_; int* cmi_ = S.ij_ + n + m;
	klu_common* Common_ptr;
	klu_symbolic* Symbolic;
	klu_numeric* Numeric;
	int** ipiv_ptr = *ipiv_ptr_ptr = new int*[2] ;  void** TM_ptr = *TM_ptr_ptr = new void*[3] ;
	int* Ap = ipiv_ptr[0] = new int[n + 1]; Ap[0] = 0; for (i = 0; i < n; i++) Ap[i + 1] = Ap[i] + ni_[i];
	int* Ai = ipiv_ptr[1] = new int[ntnz]; for (i = 0; i < ntnz; i++)  Ai[i] = cmi_[i] - 1;
	TM_ptr[0] = Common_ptr = new klu_common; klu_defaults(Common_ptr);
	TM_ptr[1] = Symbolic = klu_analyze(n, Ap, Ai, Common_ptr);
	TM_ptr[2] = Numeric = klu_factor(Ap, Ai, Ax, Symbolic, Common_ptr);
	r = klu_sort(Symbolic, Numeric, Common_ptr);  assert(r);
}

void KLU_free(int*** ipiv_ptr_ptr, void*** TM_ptr_ptr) { // libere et anNULLe les espaces de travail
	if (*ipiv_ptr_ptr) {
		assert(*TM_ptr_ptr);
		if ((*ipiv_ptr_ptr)[0]){ delete[] (*ipiv_ptr_ptr)[0];} (*ipiv_ptr_ptr)[0] = NULL;
		if ((*ipiv_ptr_ptr)[1]){ delete[](*ipiv_ptr_ptr)[1];} (*ipiv_ptr_ptr)[1] = NULL;
		delete[](*ipiv_ptr_ptr); *ipiv_ptr_ptr = NULL;
		assert((*TM_ptr_ptr)[0]);
		if ((*TM_ptr_ptr)[2]) klu_free_numeric((klu_numeric**)(&(*TM_ptr_ptr)[2]), (klu_common*)(*TM_ptr_ptr)[0]);
		if ((*TM_ptr_ptr)[1]) klu_free_symbolic((klu_symbolic**)(&(*TM_ptr_ptr)[1]), (klu_common*)(*TM_ptr_ptr)[0]);
		delete (*TM_ptr_ptr)[0];  (*TM_ptr_ptr)[0] = (*TM_ptr_ptr)[1] = (*TM_ptr_ptr)[2] = NULL;
		delete[](*TM_ptr_ptr); *TM_ptr_ptr = NULL;
	}
	else {
		assert(! (*TM_ptr_ptr));
		cout << "KLU_free : Warning : les espaces ipiv_ptr (Ap, Ax) et TM_ptr (Common, Symbolic, Numeric) avaient deja ete liberes" << endl;
	}
}

Fortran_vector KLU_solve(const Sparse_matrix &S, const Fortran_vector &y, int** ipiv_ptr, void** TM_ptr, bool full_check) { // resoud   Sx = y    en retournant  x =  inv(S) * y
	// si les espaces de travail (ipiv_ptr et TM_ptr) propres a la matrice S  sont NULLs, ils sont crees provisoirement et detruits a la sortie ; sinon, ils sont MAJ si necessaire.
	Fortran_vector x(y.size());
	r = x.set_KLU_solve(S, y, ipiv_ptr, TM_ptr, full_check);
	if (!r) {
		cout << "KLU_solve()  failed -- zero-dim vector is returned" << endl;
		return(Fortran_vector(0));
	}
	return x;
}

int Fortran_vector::set_KLU_solve(const Sparse_matrix &S, const Fortran_vector &y, int** ipiv_ptr, void** TM_ptr, bool full_check) { // resoud   Sx = y    en retournant  x =  inv(S) * y
	// si les espaces de travail (ipiv_ptr et TM_ptr) propres a la matrice S  sont NULLs, ils sont crees provisoirement et detruits a la sortie ; sinon, ils sont MAJ si necessaire.
	int n = S.n_; assert((int)y.v_[0] == n) ; assert((int)v_[0] == n) ;
	return(set_KLU_solve_(v_ + 1, 1, S, y.v_ + 1, ipiv_ptr, TM_ptr, full_check));
}

int set_KLU_solve_(double*x, int p, const Sparse_matrix &S, double*b, int** ipiv_ptr, void** TM_ptr, bool full_check) { // resoud   Sx = b    en retournant  x =  inv(S) * b
	// p = nb de cols de la matrice B, donnee en 1! array=B.v_+1 -- si p=1 c'est un vecteur. !!! suppose l'array x deja reserve (b l'est par definition du probleme Sx = b) !!!
	// si les espaces de travail (ipiv_ptr et TM_ptr) propres a la matrice S  sont NULLs, ils sont crees provisoirement et detruits a la sortie ; sinon, ils sont MAJ si necessaire.
	int i, r, n = S.n_, m = S.m_, ntnz = S.ntnz_; bool null_input_buffers(false);
	assert(m == n);
	for (i = 0; i < n*p; i++)   x[i] = b[i] ;
	int* Ap; int* Ai; double* Ax = S.v_; int* ni_ = S.ij_; int* cmi_ = S.ij_ + n + m;
	klu_common* Common_ptr;
	klu_symbolic* Symbolic;
	klu_numeric* Numeric;
	if (ipiv_ptr) {
		assert(TM_ptr);
		assert(*ipiv_ptr); assert(*TM_ptr);
		Common_ptr = (klu_common*)TM_ptr[0]; assert(Common_ptr->btf <= 1);
		Symbolic = (klu_symbolic*)TM_ptr[1]; assert(Symbolic); assert(Symbolic->n == n);
		Numeric = (klu_numeric*)TM_ptr[2]; assert(Numeric); assert(Numeric->n == n);
		Ap = ipiv_ptr[0]; Ai = ipiv_ptr[1]; assert(Ap); assert(Ai);
		assert(Ap[0] == 0); assert(Ap[n] == Symbolic->nz);
		if (Ap[n] != ntnz) {
			for (i = 0; i < n; i++) Ap[i + 1] = Ap[i] + ni_[i];
			delete[] Ai;  ipiv_ptr[1] = Ai = new int[ntnz]; for (i = 0; i < ntnz; i++)  Ai[i] = cmi_[i] - 1;
			r = klu_free_symbolic(&Symbolic, Common_ptr); assert(r);
			TM_ptr[1] = Symbolic = klu_analyze(n, Ap, Ai, Common_ptr);
			r = klu_free_numeric(&Numeric, Common_ptr); assert(r);
			TM_ptr[2] = Numeric = klu_factor(Ap, Ai, Ax, Symbolic, Common_ptr);
			r = klu_sort(Symbolic, Numeric, Common_ptr);  assert(r);
		}
		else {
			assert(Ai[0] == cmi_[0] - 1); assert(Ai[ntnz - 1] == cmi_[ntnz-1] - 1); assert(Ap[1] == ni_[0]);
			if (full_check) {
				for (i = 1 ; i < n ; i ++)  assert(Ap[i+1] == Ap[i ] + ni_[i]);
				for (i = 1 ; i < ntnz - 1 ; i ++)  assert(Ai[i] == cmi_[i] - 1);
			}
			r = klu_refactor(Ap, Ai, Ax, Symbolic, Numeric, Common_ptr); assert(r);
			r = klu_rcond(Symbolic, Numeric, Common_ptr); assert(r);
			if (Common_ptr->rcond < RCOND_LIM) {
				r = klu_condest(Ap, Ax, Symbolic, Numeric, Common_ptr); assert(r);
				if (Common_ptr->condest > CONDEST_LIM) {
					r = klu_free_numeric(&Numeric, Common_ptr); assert(r);
					TM_ptr[2] = Numeric = klu_factor(Ap, Ai, Ax, Symbolic, Common_ptr);
					r = klu_sort(Symbolic, Numeric, Common_ptr);  assert(r);
				}
			}
		}
	}
	else {
		null_input_buffers = true; assert(!TM_ptr);
		ipiv_ptr = new int*[2]; TM_ptr = new void*[3];
		ipiv_ptr[0] = Ap = new int[n + 1]; Ap[0] = 0; for (i = 0; i < n; i++) Ap[i + 1] = Ap[i] + ni_[i];
		ipiv_ptr[1] = Ai = new int[ntnz]; for (i = 0; i < ntnz; i++)  Ai[i] = cmi_[i] - 1;
		TM_ptr[0] = Common_ptr = new klu_common; klu_defaults(Common_ptr);
		TM_ptr[1] = Symbolic = klu_analyze(n, Ap, Ai, Common_ptr);
		TM_ptr[2] = Numeric = klu_factor(Ap, Ai, Ax, Symbolic, Common_ptr);
		r = klu_sort(Symbolic, Numeric, Common_ptr);  assert(r);
	}
	r = klu_solve(Symbolic, Numeric, n, p, x, Common_ptr); assert(r);
	if (null_input_buffers) {
		KLU_free(&ipiv_ptr, &TM_ptr);
	}
	return(r);
}


