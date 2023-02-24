/*
* PiafMunch (v.2) -- Implementing and Solving the Munch model of phloem sap flow in higher plants
*
* Copyright (C) 2004-2019 INRA
*
* Author: A. Lacointe, UMR PIAF, Clermont-Ferrand, France
*
* File: sparse_matrix.cpp
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

#include <vector>
#include <list>

#include <PiafMunch/PM_arrays.h>


bool Sparse_add_set_ij_check(int * ij_, int * ij_1, int * ij_2) ; //
void Sparse_add_fill_ij_(int* ij, int*ij1, int* ij2, const int& m, const int& n, const int& nnz, const int& nnz1, const int& nnz2);
void Sparse_add_set_ij_count_update(int * ij_, int * ij_1, int * ij_2) ;
void Sparse_add_set_ij_uncount(int * ij_, int * ij_1, int * ij_2) ;
int Sparse_add_set_nnz_(int * ij1, int * ij2, const int & m, const int & n, const int & nnz1, const int & nnz2) ;
int Sparse_elemult_set_nnz_(int * ij1, int * ij2, const int & m, const int & n, const int & nnz1, const int & nnz2) ;
bool Sparse_elemult_set_ij_check(int * ij_, int * ij_1, int * ij_2) ;
void Sparse_elemult_fill_ij_(int* ij, int*ij1, int* ij2, const int &m, const int &n, const int& nnz, const int& nnz1, const int& nnz2);
void Sparse_elemult_set_ij_count_update(int * ij_, int * ij_1, int * ij_2) ;
void Sparse_elemult_set_ij_uncount(int * ij_, int * ij_1, int * ij_2) ;
int Sparse_matmult_set_nnz_(int * ij1, int * ij2, const int &m, const int &n, const int &p, const int & nnz1, const int & nnz2, double* v) ;
bool Sparse_matmult_set_ij_check(int * ij_, int * ij_1, int * ij_2) ;
void Sparse_matmult_fill_ij_(int* ij, int*ij1, int* ij2, const int &m, const int &n, const int &p, const int &nnz, const int& nnz1, const int& nnz2, double* v);
void Sparse_matmult_set_ij_count_update(int * ij_, int * ij_1, int * ij_2) ;
void Sparse_matmult_set_ij_uncount(int * ij_, int * ij_1, int * ij_2) ;
void Sparse_XXX_set_ij_uncount(int * ij_) ;
int SpUnit_matmult_set_nnz_(int * ij1, int * ij2, const int &m, const int &n, const int &p, const int & nnz1, const int & nnz2, int* v) ;
void SpUnit_matmult_fill_ij_(int* ij, int*ij1, int* ij2, const int &m, const int &n, const int &p, const int &nnz, const int& nnz1, const int& nnz2, int* v);

double r ;

int Sparse_add_set_ij_count = 0 ;
int ** Sparse_add_set_ij_ = NULL ;
int Sparse_add_set_ij_max = 100 ;

int Sparse_elemult_set_ij_count = 0 ;
int ** Sparse_elemult_set_ij_ = NULL ;
int Sparse_elemult_set_ij_max = 100 ;

int Sparse_matmult_set_ij_count = 0 ;
int ** Sparse_matmult_set_ij_ = NULL ;
int Sparse_matmult_set_ij_max = 100 ;

int Sparse_addsubdiag_set_ij_count = 0 ;
int ** Sparse_addsubdiag_set_ij_ = NULL ;
int * Sparse_addsubdiag_i_set_ij_ = NULL ;
int Sparse_addsubdiag_set_ij_max = 100;


Sparse_matrix::Sparse_matrix(): m_(0), n_(0), ntnz_(0), npnz_(0), ij_(NULL), v_(NULL) { // constructeur par defaut, ne cree pas les arrays ij_ et v_
}

Sparse_matrix::Sparse_matrix(const int &m, const int &n, const int &npnz): m_(m), n_(n), ntnz_(0), npnz_(npnz), ij_(NULL), v_(NULL) {
	ij_ = new int[n + m + npnz + npnz] ; v_ = new double[npnz + npnz] ; // marche pour npnz_ = 0 ??
	if (( v_ == NULL ) || (ij_ == NULL)) {
		_LogMessage("ConstrSpMatrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
		assert(false);
	}
	int i, nt = n + m ;
	for (i = 0 ; i < nt ; i ++) ij_[i] = 0 ; // ni_[] et nj_[]
}

Sparse_matrix::Sparse_matrix(const Sparse_matrix &S, const int &npnz): m_(S.m_), n_(S.n_), ntnz_(S.ntnz_), npnz_(npnz), ij_(NULL), v_(NULL) {
	// constructeur par recopie : si npnz = -1 (defaut), garde le S.npnz_ ; si = -2, restreint npnz_ a ntnz_ .
	if (npnz == -1) {
		npnz_ = S.npnz_ ;
	}else {
		if (npnz == -2) {
			npnz_ = ntnz_ ;
		} else {
			if (npnz < ntnz_) {
				cout << "Constructeur par recopie (Sparse_matrix) : expanding npnz_  from " << npnz << " to actual ntnz_ = " << ntnz_ << endl;
				npnz_ = ntnz_ ;
			}
		}
	}
	ij_ = new int[n_ + m_ + npnz_ + npnz_];
	v_ = new double[npnz_ + npnz_];
	if (( v_ == NULL ) || (ij_ == NULL)) {
		_LogMessage("ConstrSpMatrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
		assert(false);
	}
	int k ; int * tempS_ij = S.ij_ ; int* temp_ij = ij_ ;
	int nt = n_+m_+ ntnz_ ; for (k = 0 ; k < nt ; k ++)  temp_ij[k] = tempS_ij[k] ; // copie de ni_, nj, cmi_
	temp_ij += nt + npnz_ - ntnz_ ; tempS_ij += nt + S.npnz_ - ntnz_ ;  // pointent sur rmj_
	for (k = 0 ; k < ntnz_ ; k ++)  temp_ij[k] = tempS_ij[k] ; // copie de rmj_
	double* tempS_v = S.v_ ; double* temp_v = v_ ;
	for (k = 0 ; k < ntnz_ ; k ++)  temp_v[k] = tempS_v[k] ; // copie de cmv_
	temp_v += npnz_ ; tempS_v += S.npnz_ ;  // pointent sur rmv_
	for (k = 0 ; k < ntnz_ ; k ++)  temp_v[k] = tempS_v[k] ; // copie de rmv_
}

Sparse_matrix::Sparse_matrix(const SpUnit_matrix &U): m_(U.m_), n_(U.n_), ntnz_(U.nnz_), npnz_(U.nnz_), ij_(NULL), v_(NULL) {
	assert(ntnz_ != 0) ;
	int k, nt = n_ + m_ + npnz_ + npnz_ ;
	ij_ = new int[nt] ;	v_ = new double[npnz_ + npnz_] ;
	if (( v_ == NULL ) || (ij_ == NULL)) {
		_LogMessage("Sparse_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
		assert(false);
	}
	int * temp_ij = U.ij_ ;
	for (k = 0 ; k < nt ; k ++)  ij_[k] = temp_ij[k] ;
	temp_ij = U.ij_ + nt ; nt = ntnz_ + ntnz_ ;
	for (k = 0 ; k < nt ; k ++) {
		if (temp_ij[k]) v_[k] = 1.;
		else v_[k] = -1.;
	}
}

Sparse_matrix& Sparse_matrix::operator=(const Sparse_matrix &S) { // operateur de recopie inconditionnelle => perte des infos des registres Sparse_XXX_set_ij_ :
	Sparse_XXX_set_ij_uncount(ij_) ;
	int nt = S.n_ + S.m_ + S.npnz_ + S.npnz_ ;
	if (nt != n_+m_+npnz_+npnz_) {
		if (ij_ != NULL) delete [] ij_ ;
		ij_ = new int[nt];
	}
	if (npnz_ != S.npnz_) {
		if (v_ != NULL)  delete [] v_ ;
		v_ = new double[S.npnz_ + S.npnz_];
	}
	if (( v_ == NULL ) || (ij_ == NULL)) {
		_LogMessage("ConstrSpMatrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
		assert(false);
	}
	m_ = S.m_ ; n_ = S.n_ ; npnz_ = S.npnz_ ; ntnz_ = S.ntnz_ ;
	int k ; int * tempS_ij = S.ij_ ; int* temp_ij = ij_ ;
	nt = n_+m_+ ntnz_ ; for (k = 0 ; k < nt ; k ++)  temp_ij[k] = tempS_ij[k] ;
	temp_ij += nt + npnz_ - ntnz_ ; tempS_ij += nt + npnz_ - ntnz_ ;
	for (k = 0 ; k < ntnz_ ; k ++)  temp_ij[k] = tempS_ij[k] ;
	double* tempS_v = S.v_ ; double* temp_v = v_ ;
	for (k = 0 ; k < ntnz_ ; k ++)  temp_v[k] = tempS_v[k] ;
	temp_v += npnz_ ; tempS_v += npnz_ ;
	for (k = 0 ; k < ntnz_ ; k ++)  temp_v[k] = tempS_v[k] ;
	return *this;
}

void Sparse_matrix::reset_npnz(int npnz) {
	assert((ij_ != NULL) && (v_ != NULL)) ;
	if (npnz <= 0) npnz = ntnz_ ;
	else 	assert(npnz >= ntnz_) ;
	if (npnz != npnz_) {
		int* new_ij = new int[n_ + m_ + npnz + npnz];
		double* new_v = new double[npnz + npnz];
		if (( new_v == NULL ) || (new_ij == NULL)) {
			_LogMessage("ConstrSpMatrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
			assert(false);
		}
		int k ; int * tempN_ij = new_ij ; int* temp_ij = ij_ ;
		int nt = n_+m_+ ntnz_ ; for (k = 0 ; k < nt ; k ++)  tempN_ij[k] = temp_ij[k] ;
		temp_ij += nt + npnz_ - ntnz_ ; tempN_ij += nt + npnz - ntnz_ ;
		for (k = 0 ; k < ntnz_ ; k ++)  tempN_ij[k] = temp_ij[k] ;
		double* tempN_v = new_v ; double* temp_v = v_ ;
		for (k = 0 ; k < ntnz_ ; k ++)  tempN_v[k] = temp_v[k] ;
		temp_v += npnz_ ; tempN_v += npnz ;
		for (k = 0 ; k < ntnz_ ; k ++)  tempN_v[k] = temp_v[k] ;
		Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; delete [] v_ ;
		ij_ = new_ij ; v_ = new_v ;
		npnz_ = npnz ;
	}
}

const double& Sparse_matrix::operator()(int i, int j) const { // retourne la valeur de l'element d'indices i,j
	if ((i < 1) || (i > m_) || (j < 1) || (j > n_)) assert(false) ;
	int* ni_ = ij_  ; int* nj_ = ij_ + n_ ; r = 0. ;
	if ((nj_[i-1] == 0) || (ni_[j-1] == 0))   	return r ;
	int k ; int p = 0 ;
	if (i < j) {
		int* rmj_ = nj_ + m_ + npnz_ ;
		for (k = 0 ; k < i-1 ; k ++)   p += nj_[k] ;
		rmj_ += p ;
		for (k = 0 ; k < nj_[i-1] ; k ++) {
			if (rmj_[k] == j)  	return(v_[npnz_ + p + k]) ;
		}
		return r ;
	} else  {
		int* cmi_ = nj_ + m_ ;
		for (k = 0 ; k < j-1 ; k ++)   p += ni_[k] ;
		cmi_ += p ;
		for (k = 0 ; k < ni_[j-1] ; k ++) {
			if (cmi_[k] == i)  return(v_[p + k]) ;
		}
		return r ;
	}
}

void Sparse_matrix::add_(int i, int j, double a, int ad) {
	// affecte (si ad =0, macro 'set_(i,j,a)'), ajoute (si (ad =1) ou soustrait (si ad =-1, macro 'sub_(i,j,a)') la valeur a a l'element d'indices (i,j)
	if ((i < 1) || (i > m_) || (j < 1) || (j > n_)) assert(false) ;
	assert((ad == 0) || (ad == 1) || (ad == -1)) ;
	int* ni_ = ij_  ; int* nj_ = ij_ + n_ ;
	int k, ki, kj ; int p ;
	int* cmi_ = nj_ + m_ ; int* rmj_ = cmi_ + npnz_ ;
	double* cmv_ = v_ ; double* rmv_ = v_ + npnz_ ;
	p = 0 ; for (k = 0 ; k < j-1 ; k ++)   p += ni_[k] ;
	cmi_ += p ;
	ki = -1 ;
	for (k = 0 ; k < ni_[j-1] ; k ++) {
		if (cmi_[k] >= i) {
			if (cmi_[k] == i)  ki = k ;
			break ;
		}
	}
	if (ki >= 0) {
		if (ad == 0)  cmv_[p + k] = a ; // macro 'set_()'
		else {
			if (ad == 1)  cmv_[p+k] += a ; // 'add_()' s.s.
			else  cmv_[p+k] -= a ; // macro 'sub_()'
		}
		p = 0 ; for (k = 0 ; k < i-1 ; k ++)   p += nj_[k] ;
		rmj_ += p ;
		kj = -1 ;
		for (k = 0 ; k < nj_[i-1] ; k ++) {
			if (rmj_[k] >= j) {
				if (rmj_[k] == j)  kj = k ;
				break ;
			}
		}
		assert(kj >= 0) ;
		if (ad == 0)  rmv_[p + k] = a ; // macro 'set_()'
		else {
			if (ad == 1)  rmv_[p+k] += a ; // 'add_()' s.s.
			else  rmv_[p+k] -= a ; // macro 'sub_()'
		}
	} else {
		Sparse_XXX_set_ij_uncount(ij_) ;
		int *T_ij ; int* tempT_ij ; int* temp_ij ; double *T_v ; double *tempT_v ; double* temp_v ; int nt, T_npnz ;
		if (npnz_ == ntnz_) {
			T_npnz = 3 * npnz_ / 2 ; T_ij = new int[m_ + n_ + T_npnz + T_npnz] ; T_v = new double[T_npnz + T_npnz] ;
			nt = n_ + m_ + p + k ; for (ki = 0 ; ki < nt ; ki ++)  T_ij[ki] = ij_[ki] ;
			nt = p + k ; for (ki = 0 ; ki < nt ; ki ++)  T_v[ki] = v_[ki] ;
		} else  {
			T_ij = ij_ ; T_v = v_ ; T_npnz = npnz_ ;
		}
		T_ij[j-1] = ni_[j-1] + 1 ;
		temp_ij = ij_ ; tempT_ij = T_ij ; temp_v = v_ ; tempT_v = T_v ;
		tempT_ij += n_ + m_ + p + k ; temp_ij += n_ + m_ + p + k - 1 ;
		nt = ntnz_ - p - k ; for (ki = nt ; ki > 0 ; ki --)  tempT_ij[ki] = temp_ij[ki] ;
		*tempT_ij = i ;
		tempT_v += p + k ; temp_v += p + k - 1 ;
		for (ki = nt ; ki > 0 ; ki --)  tempT_v[ki] = temp_v[ki] ;
		if (ad >= 0)  *tempT_v = a ;
		else   *tempT_v = -a ;
		p = 0 ; for (k = 0 ; k < i-1 ; k ++)   p += nj_[k] ;
		rmj_ += p ;
		kj = -1 ;
		for (k = 0 ; k < nj_[i-1] ; k ++) {
			if (rmj_[k] >= j) {
				if (rmj_[k] == j)  kj = k ;
				break ;
			}
		}
		assert (kj < 0) ;
		temp_ij = ij_ + n_ + m_ + npnz_ ; tempT_ij = T_ij + n_ + m_ + T_npnz ; temp_v = v_ + npnz_ ; tempT_v = T_v + T_npnz ;
		if (T_ij != ij_) {
			nt = p + k ; for (ki = 0 ; ki < nt ; ki ++)  tempT_ij[ki] = temp_ij[ki] ;
			for (ki = 0 ; ki < nt ; ki ++)  tempT_v[ki] = temp_v[ki] ;
		}
		T_ij[n_ + i-1] = nj_[i-1] + 1 ;
		tempT_ij += p + k ; temp_ij += p + k - 1 ;
		nt = ntnz_ - p - k ; for (kj = nt ; kj > 0 ; kj --)  tempT_ij[kj] = temp_ij[kj] ; // copie de la fin de rmj_, decalee de 1
		*tempT_ij = j ;
		tempT_v += p + k ; temp_v += p + k - 1 ;
		for (ki = nt ; ki > 0 ; ki --)  tempT_v[ki] = temp_v[ki] ;
		if (ad >= 0)  *tempT_v = a ;
		else   *tempT_v = -a ;
		ntnz_ ++ ; if (T_ij != ij_) {delete[] v_ ; delete[] ij_ ; v_ = T_v ; ij_ = T_ij  ; npnz_ = T_npnz ;}
	}
}

Sparse_matrix::~Sparse_matrix() {
	Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; delete [] v_ ;
}

int Sparse_matrix::nblin() const { return m_; }

int Sparse_matrix::nbcol() const { return n_; }

int Sparse_matrix::ntnz() const { return ntnz_; }

int Sparse_matrix::npnz() const { return npnz_; }

void Sparse_matrix::set(const Sparse_matrix &S, bool check) { // n'impose pas l'egalite des npnz_ .  si (check) --> verifie le nombre et la position  des NZ
	if ((S.m_ != m_) || (S.n_ != n_)) assert(false) ;
	if (S.ntnz_ != ntnz_) { operator=(S) ; return ; }
	int k, kk, nt = n_+m_+ ntnz_ , *tempS_ij = S.ij_ , *temp_ij = ij_ ; bool OK = true ;
	if(check) {
		kk = 0 ; for (k = kk ; k < n_ ; k ++)  if(temp_ij[k] != tempS_ij[k])  { OK = false ; kk= k ; break ; }
		if(! OK) { for (k = kk ; k < n_ ; k ++)  temp_ij[k] = tempS_ij[k] ; }
		kk = n_+m_ ; if(OK)  for (k = kk ; k < nt ; k ++)  if(temp_ij[k] != tempS_ij[k])  { OK = false ; kk = k ; break ; }
		if(! OK) {
			for (k = kk ; k < nt ; k ++){  temp_ij[k] = tempS_ij[k] ;} Sparse_XXX_set_ij_uncount(ij_) ;
			for (k = n_ ; k < m_+n_ ; k ++)  temp_ij[k] = tempS_ij[k] ;
			temp_ij += n_ + m_ + npnz_ ; tempS_ij += n_ + m_ + S.npnz_ ;
			for (k = 0 ; k < ntnz_ ; k ++)  temp_ij[k] = tempS_ij[k] ;
		}
	}
	double* tempS_v = S.v_ ; double* temp_v = v_ ;
	for (k = 0 ; k < ntnz_ ; k ++)  temp_v[k] = tempS_v[k] ;
	temp_v += npnz_ ; tempS_v += S.npnz_ ;
	for (k = 0 ; k < ntnz_ ; k ++)  temp_v[k] = tempS_v[k] ;
}

void Sparse_matrix::dump() {
	int i ;
	cout << "Sparse_matrix : m_ = " << m_ << " ; n_ = " << n_ << " ; ntnz_ = " << ntnz_ << " ; npnz_ = " << npnz_ << endl << "ni_[] =" ;
	for (i = 0 ; i < n_ ; i ++)  cout << " " << ij_[i] ;
	cout << endl << "nj_[] =";
	for (i = 0 ; i < m_ ; i ++)  cout << " " << ij_[n_ + i] ;
	cout << endl << "cmi_[] =";
	for (i = 0 ; i < ntnz_ ; i ++)  cout << " " << ij_[n_ + m_ + i] ;
	cout << endl << "rmj_[] =";
	for (i = 0 ; i < ntnz_ ; i ++)  cout << " " << ij_[n_ + m_ + npnz_ + i] ;
	cout << endl << "cmv_[] =";
	for (i = 0 ; i < ntnz_ ; i ++)  cout << " " << v_[i] ;
	cout << endl << "rmv_[] =";
	for (i = 0 ; i < ntnz_ ; i ++)  cout << " " << v_[npnz_ + i] ;
	cout << endl << endl ;
}

void Sparse_matrix::add_subdiag(const int & i1, const int & i2, const double &a, int ad) {
	if ((i1 < 1) || (i1 > i2) || (i2 > m_) || (i2 > n_)) assert(false);
	int mn = i2 - i1 + 1 ;
	if (ad == 0 ) { // fonction 'Sparse_matrix::set_subdiag(i1, i2, a)''
		int k ;
		if (npnz_ < mn) {
			(*this) = Sparse_matrix(m_, n_, mn) ;
			cout << "Sparse_matrix::add_subdiag() : expanding npnz_  from " << npnz_ << " to actual ntnz_ = " << mn << endl;
		}
		(*this).check_addsubdiag_set_ij_(i1,i2) ;
		ntnz_ = mn ;
		int * ij = ij_ ; double * v = v_ ;
		int * ij_max = ij + i1-1 ; while (ij < ij_max)  { *ij = 0 ; ij ++ ; }
		ij_max += mn ; while (ij < ij_max)  { *ij = 1 ; ij ++ ; }
		ij_max += n_ - i2 ; while (ij < ij_max)  { *ij = 0 ; ij ++ ; }
		ij_max += i1-1 ; while (ij < ij_max)  { *ij = 0 ; ij ++ ; }
		ij_max += mn ; while (ij < ij_max)  { *ij = 1 ; ij ++ ; }
		ij_max += m_ - i2 ; while (ij < ij_max)  { *ij = 0 ; ij ++ ; }
		for (k = i1 ; k <= i2 ; k ++)  { *ij = k ; ij ++ ; }
		ij += npnz_ - ntnz_ ; for (k = i1 ; k <= i2 ; k ++)  { *ij = k ; ij ++ ; }
		double* v_max = v_ + mn ; while (v < v_max) { *v = a; v ++ ; }
		v += npnz_ - ntnz_ ; v_max = v_ + npnz_ + mn ; while (v < v_max) { *v = a; v ++ ; }
	}
	else { // fonction  'add_subdiag(i1, i2, a)'  ou  'sub_subdiag(i1, i2, a)' suivant ad = 1 ou -1
		int i, j ;
		(*this).check_addsubdiag_set_ij_(i1,i2) ;
		double * rmv_ = v_ + npnz_ ;
		int * nj_ = ij_ + n_  - 1 ;
		int * rmj_, * rmj_max ;
		int * rmj_0 = ij_ + m_ + n_ + npnz_ ;
		int * rmj_min = rmj_0 ;
		int * ni_ = ij_ - 1 ;
		double * cmv_ = v_ ;
		int * cmi_, * cmi_max ;
		int * cmi_0 = ij_ + m_ + n_ ;
		int * cmi_min = cmi_0 ;
		if (ad == 1) { // fonction 'add_subdiag(i1, i2, a)'
			for (i = 1 ; i < i1 ; i ++)   rmj_min += nj_[i] ;
			for (i = i1 ; i <= i2 ; i ++) {
				rmj_max = rmj_min + nj_[i] ;
				for (rmj_ = rmj_min ; rmj_ < rmj_max ; rmj_ ++) {
					if ((*rmj_) == i) { rmv_[rmj_- rmj_0] += a ; break ; }
				}
				rmj_min = rmj_max ;
			}
			for (j = 1 ; j < i1 ; j ++)   cmi_min += ni_[j] ;
			for (j = i1 ; j <= i2 ; j ++) {
				cmi_max = cmi_min + ni_[j] ;
				for (cmi_ = cmi_min ; cmi_ < cmi_max ; cmi_ ++) {
					if ((*cmi_) == j) { cmv_[cmi_- cmi_0] += a ; break ; }
				}
				cmi_min = cmi_max ;
			}
		} else {  // fonction 'sub_subdiag(i1, i2, a)'
			for (i = 1 ; i < i1 ; i ++)   rmj_min += nj_[i] ;
			for (i = i1 ; i <= i2 ; i ++) {
				rmj_max = rmj_min + nj_[i] ;
				for (rmj_ = rmj_min ; rmj_ < rmj_max ; rmj_ ++) {
					if ((*rmj_) == i) { rmv_[rmj_- rmj_0] -= a ; break ; }
				}
				rmj_min = rmj_max ;
			}
			for (j = 1 ; j < i1 ; j ++)   cmi_min += ni_[j] ;
			for (j = i1 ; j <= i2 ; j ++) {
				cmi_max = cmi_min + ni_[j] ;
				for (cmi_ = cmi_min ; cmi_ < cmi_max ; cmi_ ++) {
					if ((*cmi_) == j) { cmv_[cmi_- cmi_0] -= a ; break ; }
				}
				cmi_min = cmi_max ;
			}
		}
	}
}

void Sparse_matrix::add_subdiag(const int & i1, const int & i2, const Fortran_vector & V, int ad) {
	double * tempV = V.v_ ;
	if ((i1 < 1) || (i1 > i2) || (i2 > m_) || (i2 > n_)) assert(false);
	int mn = i2 - i1 + 1 ;
	if (mn != (int)tempV[0]) assert(false);
	double * cmv_ ; double * rmv_ ;
	if (ad == 0 ) { // fonction 'Sparse_matrix::set_subdiag(i1, i2, V)''
		int k ;
		if (npnz_ < mn) {
			(*this) = Sparse_matrix(m_, n_, mn) ;
			cout << "Sparse_matrix::add_subdiag() : expanding npnz_  from " << npnz_ << " to actual ntnz_ = " << mn << endl;
		}
		(*this).check_addsubdiag_set_ij_(i1,i2) ;
		ntnz_ = mn ;
		int * ij = ij_ ;
		int * ij_max = ij + i1-1 ; while (ij < ij_max)  { *ij = 0 ; ij ++ ; }
		ij_max += mn ; while (ij < ij_max)  { *ij = 1 ; ij ++ ; }
		ij_max += n_ - i2 ; while (ij < ij_max)  { *ij = 0 ; ij ++ ; }
		ij_max += i1-1 ; while (ij < ij_max)  { *ij = 0 ; ij ++ ; }
		ij_max += mn ; while (ij < ij_max)  { *ij = 1 ; ij ++ ; }
		ij_max += m_ - i2 ; while (ij < ij_max)  { *ij = 0 ; ij ++ ; }
		for (k = i1 ; k <= i2 ; k ++)  { *ij = k ; ij ++ ; }
		ij += npnz_ - ntnz_ ; for (k = i1 ; k <= i2 ; k ++)  { *ij = k ; ij ++ ; }
		cmv_ = v_ - 1 ; rmv_ = cmv_ + npnz_ ;
		for (k = 1 ; k <= mn ; k ++)  { cmv_[k] = rmv_[k] = tempV[k] ; }
	}
	else { // fonction  'add_subdiag(i1, i2, V)'  ou  'sub_subdiag(i1, i2, V)' suivant ad = 1 ou -1
		int i, j ;
		(*this).check_addsubdiag_set_ij_(i1,i2) ;
		rmv_ = v_ + npnz_ ;
		int * nj_ = ij_ + n_ - 1 ;
		int * rmj_, * rmj_max ;
		int * rmj_0 = ij_ + m_ + n_ + npnz_ ;
		int * rmj_min = rmj_0 ;
		int * ni_ = ij_ - 1 ;
		cmv_ = v_ ;
		int * cmi_, * cmi_max ;
		int * cmi_0 = ij_ + m_ + n_ ;
		int * cmi_min = cmi_0 ;
		if (ad == 1) { // fonction 'add_subdiag(i1, i2, V'
			for (i = 1 ; i < i1 ; i ++)   rmj_min += nj_[i] ;
			for (i = i1 ; i <= i2 ; i ++) {
				rmj_max = rmj_min + nj_[i] ; tempV ++ ;
				for (rmj_ = rmj_min ; rmj_ < rmj_max ; rmj_ ++) {
					if ((*rmj_) == i) { rmv_[rmj_- rmj_0] += *tempV ; break ; }
				}
				rmj_min = rmj_max ;
			}
			tempV = V.v_ ;
			for (j = 1 ; j < i1 ; j ++)   cmi_min += ni_[j] ;
			for (j = i1 ; j <= i2 ; j ++) {
				cmi_max = cmi_min + ni_[j] ; tempV ++ ;
				for (cmi_ = cmi_min ; cmi_ < cmi_max ; cmi_ ++) {
					if ((*cmi_) == j) { cmv_[cmi_- cmi_0] += *tempV ; break ; }
				}
				cmi_min = cmi_max ;
			}
		} else {  // fonction 'sub_subdiag(i1, i2, V'
			for (i = 1 ; i < i1 ; i ++)   rmj_min += nj_[i] ;
			for (i = i1 ; i <= i2 ; i ++) {
				rmj_max = rmj_min + nj_[i] ; tempV ++ ;
				for (rmj_ = rmj_min ; rmj_ < rmj_max ; rmj_ ++) {
					if ((*rmj_) == i) { rmv_[rmj_- rmj_0] -= *tempV ; break ; }
				}
				rmj_min = rmj_max ;
			}
			tempV = V.v_ ;
			for (j = 1 ; j < i1 ; j ++)   cmi_min += ni_[j] ;
			for (j = i1 ; j <= i2 ; j ++) {
				cmi_max = cmi_min + ni_[j] ; tempV ++ ;
				for (cmi_ = cmi_min ; cmi_ < cmi_max ; cmi_ ++) {
					if ((*cmi_) == j) { cmv_[cmi_- cmi_0] -= *tempV ; break ; }
				}
				cmi_min = cmi_max ;
			}
		}
	}
}

void Sparse_matrix::set_diag(const double &a) {
	(*this).add_subdiag(1, min(m_, n_), a, 0) ;
}

void Sparse_matrix::set_diag(const Fortran_vector &V) {
	(*this).add_subdiag(1, V.size(), V, 0) ;
}

void Sparse_matrix::set_subdiag(const int & i1, const int & i2, const double &a) {
	(*this).add_subdiag(i1, i2, a, 0) ;
}

void Sparse_matrix::set_subdiag(const int & i1, const int & i2, const Fortran_vector &V) {
	(*this).add_subdiag(i1, i2, V, 0) ;
}

void Sparse_matrix::add_diag(const double &a, int ad) {
	(*this).add_subdiag(1, min(m_, n_), a, ad) ;
}

void Sparse_matrix::add_diag(const Fortran_vector &V, int ad) {
	(*this).add_subdiag(1, V.size(), V, ad) ;
}

Sparse_matrix diag(const double &a, const int m, int n) {
	int ntnz ;
	if (n == -1) { n = m ; ntnz = m ;} // (n = -1) : valeur par defaut = convention signalant une matrice carree
	else  ntnz = min(m, n) ;
	Sparse_matrix T(m, n, ntnz) ;
	T.set_diag(a);
	return T ;
}

Sparse_matrix diag(const Fortran_vector & V, const int m, const int n) {
	int ntnz = min(m, n) ;
	Sparse_matrix T(m, n, ntnz) ;
	T.set_diag(V);
	return T ;
}

Sparse_matrix diag(const Fortran_vector & V) {
	double * tempV ; tempV = V.v_ ;
	int m = (int)tempV[0] ;
	Sparse_matrix T(m, m, m) ;
	T.set_diag(V);
	return T ;
}

Sparse_matrix transpose(const Sparse_matrix & S, int npnz) { // si npnz = -1 (defaut), garde le S.npnz_ ; si = -2, restreint npnz_ a ntnz_ .
	int ntnz = S.ntnz_ ;
	if (npnz == -1) {
		npnz = S.npnz_ ;
	}else {
		if (npnz == -2) {
			npnz = ntnz ;
		} else {
			if (npnz < ntnz) {
				cout << "Constructeur par recopie (Sparse_matrix) : expanding npnz_  from " << npnz << " to actual ntnz_ = " << ntnz << endl;
				npnz = ntnz ;
			}
		}
	}
	assert(npnz != 0) ;
	Sparse_matrix T ; // sans arguments ----> matrice nulle (constructeur par defaut)
	T.m_ = S.n_ ; T.n_ = S.m_ ; T.ntnz_ = ntnz ; T.npnz_ = npnz ;
	int k, nt = S.n_ + S.m_ + npnz + npnz ;
	T.ij_ = new int[nt] ;	T.v_ = new double[npnz + npnz] ;
	if (( T.v_ == NULL ) || (T.ij_ == NULL)) {
		_LogMessage("Sparse_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
		assert(false);
	}
	int * temp_S, * temp_T ;
	temp_S = S.ij_ ; temp_T = T.ij_ + T.n_ ; nt = T.m_ ;
	for (k = 0 ; k < nt ; k ++)  temp_T[k] = temp_S[k] ;
	temp_S = S.ij_ + S.n_ ; temp_T = T.ij_ ; nt = T.n_ ;
	for (k = 0 ; k < nt ; k ++)  temp_T[k] = temp_S[k] ;
	temp_S = S.ij_ + S.n_ + S.m_ ; temp_T = T.ij_ + T.n_ + T.m_ + T.npnz_ ;
	for (k = 0 ; k < ntnz ; k ++)  temp_T[k] = temp_S[k] ;
	temp_S = S.ij_ + S.n_ + S.m_ + S.npnz_ ; temp_T = T.ij_ + T.n_ + T.m_ ;
	for (k = 0 ; k < ntnz ; k ++)  temp_T[k] = temp_S[k] ;
	double* temp_S_v = S.v_ ; double* temp_T_v = T.v_ + T.npnz_ ;
	for (k = 0 ; k < ntnz ; k ++)  temp_T_v[k] = temp_S_v[k] ;
	temp_S_v = S.v_ + S.npnz_ ; temp_T_v = T.v_ ;
	for (k = 0 ; k < ntnz ; k ++)  temp_T_v[k] = temp_S_v[k] ;
	return T;
}

/**** arithmetique matricielle ordinaire ******************************************************************/

Sparse_matrix Sparse_matrix::operator-() {
	//   - S (operateur unaire de changement de signe)
	assert(npnz_ != 0) ;
	Sparse_matrix T ;
	int k, nt = n_ + m_ + npnz_ + npnz_ ;
	T.ij_ = new int[nt] ;	T.v_ = new double[npnz_ + npnz_] ;
	if (( T.v_ == NULL ) || (T.ij_ == NULL)) {
		_LogMessage("Sparse_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
		assert(false);
	}
	T.m_ = m_ ; T.n_ = n_ ; T.ntnz_ = ntnz_ ; T.npnz_ = npnz_ ;
	int * tempT_ij = T.ij_ ; int* temp_ij = ij_ ;
	nt = n_+m_+ ntnz_ ; for (k = 0 ; k < nt ; k ++)  tempT_ij[k] = temp_ij[k] ;
	tempT_ij += nt + npnz_ - ntnz_ ; temp_ij += nt + npnz_ - ntnz_ ;
	for (k = 0 ; k < ntnz_ ; k ++)  tempT_ij[k] = temp_ij[k] ;
	double* tempT_v = T.v_ ; double* temp_v = v_ ;
	for (k = 0 ; k < ntnz_ ; k ++)  tempT_v[k] = -temp_v[k] ;
	temp_v += npnz_ ; tempT_v += npnz_ ;
	for (k = 0 ; k < ntnz_ ; k ++)  tempT_v[k] = -temp_v[k] ;
	return T;
}

Sparse_matrix Sparse_matrix::operator+(Sparse_matrix &S) {
	//   M + S
	Sparse_matrix M(m_, n_);
	M.add_add((*this), S, 0, false);
	return M;
}

Sparse_matrix Sparse_matrix::operator+(const SpUnit_matrix &U) {
	//   M + U
	Sparse_matrix M(m_, n_);
	M.add_add((*this), U, 0, false);
	return M;
}

Sparse_matrix Sparse_matrix::operator-(const Sparse_matrix &S) {
	//   M - S
	Sparse_matrix M(m_, n_);
	M.add_sub((*this), S, 0, false);
	return M;
}

Sparse_matrix Sparse_matrix::operator-(const SpUnit_matrix &U) {
	//   M - U
	Sparse_matrix M(m_, n_);
	M.add_sub((*this), U, 0, false);
	return M;
}

Sparse_matrix Sparse_matrix::operator*(const Sparse_matrix &S) {
	//   M * S  elementwise
	Sparse_matrix M(m_, n_);
	M.add_elemult((*this), S, 0, false);
	return M;
}

Sparse_matrix Sparse_matrix::operator*(const SpUnit_matrix &U) {
	//   M * S  elementwise
	Sparse_matrix M(m_, n_);
	M.add_elemult(U, (*this), 0, false);
	return M;
}

Sparse_matrix operator*(const double &a, const Sparse_matrix &S) {
	// a * S : multiplication par un scalaire,
	Sparse_matrix M(S.m_, S.n_, S.ntnz_) ;
	M.add_mult(a, S, 0, false);
	return M;
}

Sparse_matrix Sparse_matrix::operator/(const Fortran_vector &v) {
	//   S / v : division elementwise, colonne par colonne, par le vecteur v  (v.size = m_)
	Sparse_matrix M(m_, n_, ntnz_) ;
	M.add_elediv((*this), v, 0, false);
	return M;
}

Sparse_matrix Sparse_matrix::operator/(const double &a) {
	//   S / a : division par un scalaire
	Sparse_matrix M(m_, n_, ntnz_) ;
	M.add_div((*this), a, 0, false);
	return M;
}

Sparse_matrix row_elemult(const Sparse_matrix &S, const Fortran_vector &v) {
	// S * v : multiplication elementwise PAR LIGNE, par un vecteur de taille n = S.nbcol()
	Sparse_matrix M(S.m_, S.n_, S.ntnz_) ;
	M.add_row_elemult(S, v, 0, false);
	return M;
}

Sparse_matrix row_elediv(const Sparse_matrix &S, const Fortran_vector &v) {
	// S / v : division elementwise PAR LIGNE, par un vecteur de taille n = S.nbcol()
	Sparse_matrix M(S.m_, S.n_, S.ntnz_) ;
	M.add_row_elediv(S, v, 0, false);
	return M;
}

Fortran_vector matmult(const Sparse_matrix &S, const Fortran_vector &v) {
	//   S * v, multiplication matricielle
	Fortran_vector Tv(S.m_) ;
	Tv.set_matmult(S, v) ;
	return Tv;
}

Sparse_matrix matmult(const Sparse_matrix &S1, const Sparse_matrix &S2, Fortran_vector * temp_v) {
	// si temp_v est fourni, il doit pointer sur un vecteur  de taille  n = S1.n_  qui servira de tampon.
	Sparse_matrix TM(S1.m_, S2.n_) ;
	TM.set_matmult(S1, S2, temp_v, false) ;
	return TM;
}

/**** arithmetique matricielle 'in-place' *****************************************************************/

Sparse_matrix & Sparse_matrix::operator+=(Sparse_matrix &S) {
	//  M += S
	if ((S.n_ != n_) || (S.m_ != m_)) assert(false);
	int i, j ;
	bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, ij_, S.ij_) ;
	if (!(Sparse_add_set_ij_done)) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		(*this).add_add(*this, S, 0) ;
	}
	else {
		int jS ;
		int * nj_ = ij_ + n_ ; int * nj_S = S.ij_ + n_  ;
		int * rmj_ = nj_ + m_ + npnz_ ; int * rmj_S = nj_S + m_ + S.npnz_ ; int * rmj_max, * rmjS_max ;
		double * rmv_ = v_ + npnz_ ;	double * rmv_S = S.v_ + S.npnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmj_max = rmj_ + nj_[i] - 1 ; rmjS_max = rmj_S + nj_S[i] - 1 ; jS = (*rmj_S) ;
			while (rmj_S <= rmjS_max) {
				if ((*rmj_) == jS) {
					*rmv_ += (*rmv_S) ; rmv_ ++ ; rmv_S ++ ; rmj_ ++ ; rmj_S ++ ; jS = (*rmj_S) ;
				} else
					{ rmv_ ++ ; rmj_ ++ ; }
			}
			rmv_ += rmj_max - rmj_ + 1 ; rmj_ = rmj_max + 1 ;
		}
		int iS ;
		int * ni_ = ij_ ; int * ni_S = S.ij_  ;
		int * cmi_ = ni_ + m_ + n_ ; int * cmi_S = ni_S + m_ + n_ ;	int * cmi_max, * cmiS_max ;
		double * cmv_ = v_ ; double * cmv_S = S.v_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmi_max = cmi_ + ni_[j] - 1 ; cmiS_max = cmi_S + ni_S[j] - 1 ;
			iS = (*cmi_S) ;
			while (cmi_S <= cmiS_max) {
				if ((*cmi_) == iS) {
					*cmv_ += (*cmv_S) ; cmv_ ++ ; cmv_S ++ ; cmi_ ++ ; cmi_S ++ ; iS = (*cmi_S) ;
				} else
					{ cmv_ ++ ; cmi_ ++ ; }
			}
			cmv_ += cmi_max - cmi_ + 1 ; cmi_ = cmi_max + 1 ;
		}
	}
	return *this;
}

Sparse_matrix & Sparse_matrix::operator+=(const SpUnit_matrix &U) {
	//  M += U
	if ((U.n_ != n_) || (U.m_ != m_)) assert(false);
	int i, j ;
	bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, ij_, U.ij_) ;
	if (!(Sparse_add_set_ij_done)) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		(*this).add_add(*this, U, 0) ;
	}
	else {
		int jU ;
		int * nj_ = ij_ + n_ ; int * nj_U = U.ij_ + n_  ;
		int * rmj_ = nj_ + m_ + npnz_ ; int * rmj_U = nj_U + m_ + U.nnz_ ; int * rmj_max, * rmjU_max ;
		double * rmv_ = v_ + npnz_ ;	int * rmv_U = rmj_U + 2 * U.nnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmj_max = rmj_ + nj_[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ; jU = (*rmj_U) ;
			while (rmj_U <= rmjU_max) {
				if ((*rmj_) == jU) {
					if (*rmv_U)  *rmv_ += 1. ; else  *rmv_ -= 1. ;
					rmv_ ++ ; rmv_U ++ ; rmj_ ++ ; rmj_U ++ ; jU = (*rmj_U) ;
				} else
					{ rmv_ ++ ; rmj_ ++ ; }
			}
			rmv_ += rmj_max - rmj_ + 1 ; rmj_ = rmj_max + 1 ;
		}
		int iU ;
		int * ni_ = ij_ ; int * ni_U = U.ij_ ;
		int * cmi_ = ni_ + m_ + n_ ; int * cmi_U = ni_U + m_ + n_ ;	int * cmi_max, * cmiU_max ;
		double * cmv_ = v_ ; int * cmv_U = cmi_U + 2 * U.nnz_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmi_max = cmi_ + ni_[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ; iU = (*cmi_U) ;
			while (cmi_U <= cmiU_max) {
				if ((*cmi_) == iU) {
					if (*cmv_U)  *cmv_ += 1. ; else  *cmv_ -= 1. ;
					cmv_ ++ ; cmv_U ++ ; cmi_ ++ ; cmi_U ++ ; iU = (*cmi_U) ;
				} else
					{ cmv_ ++ ; cmi_ ++ ; }
			}
			cmv_ += cmi_max - cmi_ + 1 ; cmi_ = cmi_max + 1 ;
		}
	}
	return *this;
}

Sparse_matrix & Sparse_matrix::operator-=(const Sparse_matrix &S) {
	//  M -= S
	if ((S.n_ != n_) || (S.m_ != m_)) assert(false);
	int i, j ;
	bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, ij_, S.ij_) ;
	if (!(Sparse_add_set_ij_done)) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		(*this).add_sub(*this, S, 0) ;
	}
	else {
		int jS ;
		int * nj_ = ij_ + n_ ; int * nj_S = S.ij_ + n_  ;
		int * rmj_ = nj_ + m_ + npnz_ ; int * rmj_S = nj_S + m_ + S.npnz_ ; int * rmj_max, * rmjS_max ;
		double * rmv_ = v_ + npnz_ ;	double * rmv_S = S.v_ + S.npnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmj_max = rmj_ + nj_[i] - 1 ; rmjS_max = rmj_S + nj_S[i] - 1 ; jS = (*rmj_S) ;
			while (rmj_S <= rmjS_max) {
				if ((*rmj_) == jS) {
					*rmv_ -= (*rmv_S) ; rmv_ ++ ; rmv_S ++ ; rmj_ ++ ; rmj_S ++ ; jS = (*rmj_S) ;
				} else
					{ rmv_ ++ ; rmj_ ++ ; }
			}
			rmv_ += rmj_max - rmj_ + 1 ; rmj_ = rmj_max + 1 ;
		}
		int iS ;
		int * ni_ = ij_ ; int * ni_S = S.ij_  ;
		int * cmi_ = ni_ + m_ + n_ ; int * cmi_S = ni_S + m_ + n_ ;	int * cmi_max, * cmiS_max ;
		double * cmv_ = v_ ; double * cmv_S = S.v_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmi_max = cmi_ + ni_[j] - 1 ; cmiS_max = cmi_S + ni_S[j] - 1 ;
			iS = (*cmi_S) ;
			while (cmi_S <= cmiS_max) {
				if ((*cmi_) == iS) {
					*cmv_ -= (*cmv_S) ; cmv_ ++ ; cmv_S ++ ; cmi_ ++ ; cmi_S ++ ; iS = (*cmi_S) ;
				} else
					{ cmv_ ++ ; cmi_ ++ ; }
			}
			cmv_ += cmi_max - cmi_ + 1 ; cmi_ = cmi_max + 1 ;
		}
	}
	return *this;
}

Sparse_matrix & Sparse_matrix::operator-=(const SpUnit_matrix &U) {
	//  M -= U
	if ((U.n_ != n_) || (U.m_ != m_)) assert(false);
	int i, j ;
	bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, ij_, U.ij_) ;
	if (!(Sparse_add_set_ij_done)) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		(*this).add_sub(*this, U, 0) ;
	}
	else {
		int jU ;
		int * nj_ = ij_ + n_ ; int * nj_U = U.ij_ + n_  ;
		int * rmj_ = nj_ + m_ + npnz_ ; int * rmj_U = nj_U + m_ + U.nnz_ ; int * rmj_max, * rmjU_max ;
		double * rmv_ = v_ + npnz_ ;	int * rmv_U = rmj_U + 2 * U.nnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmj_max = rmj_ + nj_[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ; jU = (*rmj_U) ;
			while (rmj_U <= rmjU_max) {
				if ((*rmj_) == jU) {
					if (*rmv_U)  *rmv_ -= 1. ; else  *rmv_ += 1. ;
					rmv_ ++ ; rmv_U ++ ; rmj_ ++ ; rmj_U ++ ; jU = (*rmj_U) ;
				} else
					{ rmv_ ++ ; rmj_ ++ ; }
			}
			rmv_ += rmj_max - rmj_ + 1 ; rmj_ = rmj_max + 1 ;
		}
		int iU ;
		int * ni_ = ij_ ; int * ni_U = U.ij_ ;
		int * cmi_ = ni_ + m_ + n_ ; int * cmi_U = ni_U + m_ + n_ ;	int * cmi_max, * cmiU_max ;
		double * cmv_ = v_ ; int * cmv_U = cmi_U + 2 * U.nnz_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmi_max = cmi_ + ni_[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ; iU = (*cmi_U) ;
			while (cmi_U <= cmiU_max) {
				if ((*cmi_) == iU) {
					if (*cmv_U)  *cmv_ -= 1. ; else  *cmv_ += 1. ;
					cmv_ ++ ; cmv_U ++ ; cmi_ ++ ; cmi_U ++ ; iU = (*cmi_U) ;
				} else
					{ cmv_ ++ ; cmi_ ++ ; }
			}
			cmv_ += cmi_max - cmi_ + 1 ; cmi_ = cmi_max + 1 ;
		}
	}
	return *this;
}

Sparse_matrix & Sparse_matrix::operator*=(const Sparse_matrix &S) {
	//  M *= S , elementwise
	if ((S.n_ != n_) || (S.m_ != m_)) assert(false);
	int i, j ;
	bool Sparse_elemult_set_ij_done = Sparse_elemult_set_ij_check(ij_, ij_, S.ij_) ;
	if (!(Sparse_elemult_set_ij_done)) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		(*this).add_elemult(*this, S, 0) ;
	}
	else {
		int * nj_ = ij_ + n_ ; int * nj_S = S.ij_ + n_ ; int * rmj_ = nj_ + m_ + npnz_ ;
		int * rmj_S = nj_S + m_ + S.npnz_ ; int * rmj_max, * rmjS_max ;
		double * rmv_ = v_ + npnz_ ; double * rmv_S = S.v_ + S.npnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmj_max = rmj_ + nj_[i] - 1 ; rmjS_max = rmj_S + nj_S[i] - 1 ;
			j = (*rmj_) ;
			while (rmj_ <= rmj_max) {
				if (j == (*rmj_S)) {
					*rmv_ *= (*rmv_S) ; rmv_ ++ ; rmv_S ++ ; rmj_ ++ ; j = (*rmj_) ; rmj_S ++ ;
				} else { rmv_S ++ ; rmj_S ++ ; }
			}
			rmv_S += rmjS_max + 1 - rmj_S; rmj_S = 1 + rmjS_max ;
		}
		int * ni_ = ij_ ; int * ni_S = S.ij_ ; int * cmi_ = ni_ + m_ + n_ ;
		int * cmi_S = ni_S + m_ + n_ ; int * cmi_max, * cmiS_max ;
		double * cmv_ = v_ ; double * cmv_S = S.v_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmi_max = cmi_ + ni_[j] - 1 ; cmiS_max = cmi_S + ni_S[j] - 1 ;
			i = (*cmi_) ; ;
			while (cmi_ <= cmi_max) {
				if (i == (*cmi_S)) {
					*cmv_ *= (*cmv_S) ; cmv_ ++ ; cmv_S ++ ; cmi_ ++ ; i = (*cmi_) ; cmi_S ++ ;
				} else  { cmv_S ++ ; cmi_S ++ ; }
			}
			cmv_S += cmiS_max + 1 - cmi_S; cmi_S = 1 + cmiS_max ;
		}
	}
	return *this;
}

Sparse_matrix & Sparse_matrix::operator*=(const SpUnit_matrix &U) {
	//  M *= U , elementwise
	if ((U.n_ != n_) || (U.m_ != m_)) assert(false);
	int i, j ;
	bool Sparse_elemult_set_ij_done = Sparse_elemult_set_ij_check(ij_, ij_, U.ij_) ;
	if (!(Sparse_elemult_set_ij_done)) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		(*this).add_elemult(*this, U, 0) ;
	}
	else {
		int * nj_ = ij_ + n_ ; int * nj_U = U.ij_ + n_ ; int * rmj_ = nj_ + m_ + npnz_ ;
		int * rmj_U = nj_U + m_ + U.nnz_ ; int * rmj_max, * rmjU_max ;
		double * rmv_ = v_ + npnz_ ;	int * rmv_U = rmj_U + 2 * U.nnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmj_max = rmj_ + nj_[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ;
			j = (*rmj_) ;
			while (rmj_ <= rmj_max) {
				if (j == (*rmj_U)) {
					if (*rmv_U == 0)  *rmv_ = -(*rmv_) ;
					rmv_ ++ ; rmv_U ++ ; rmj_ ++ ; j = (*rmj_) ; rmj_U ++ ;
				} else	{ rmv_U ++ ; rmj_U ++ ; }
			}
			rmv_U += rmjU_max + 1 - rmj_U ; rmj_U = 1 + rmjU_max ;
		}
		int * ni_ = ij_ ; int * ni_U = U.ij_ ; int * cmi_ = ni_ + m_ + n_ ;
		int * cmi_U = ni_U + m_ + n_ ; int * cmi_max, * cmiU_max ;
		double * cmv_ = v_ ; int * cmv_U = cmi_U + 2 * U.nnz_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmi_max = cmi_ + ni_[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ;
			i = (*cmi_) ;
			while (cmi_ <= cmi_max) {
				if (i == (*cmi_U)) {
					if (*cmv_U == 0)  *cmv_ = -(*cmv_) ;
					cmv_ ++ ; cmv_U ++ ; cmi_ ++ ; i = (*cmi_) ; cmi_U ++ ;
				} else	{ cmv_U ++ ; cmi_U ++ ; }
			}
			cmv_U += cmiU_max + 1 - cmi_U ; cmi_U = 1 + cmiU_max ;
		}
	}
	return *this;
}

Sparse_matrix & Sparse_matrix::operator*=(const double &a) {
	// S *= a , multiplication par un scalaire in place
	double * temp = v_ - 1 ;
	for (int k = 0 ; k < ntnz_ ; k ++) {
		temp ++ ; *temp *= a ;
	}
	temp += npnz_ - ntnz_ ;
	for (int k = 0 ; k < ntnz_ ; k ++) {
		temp ++ ; *temp *= a ;
	}
	return *this;
}

Sparse_matrix & Sparse_matrix::operator*=(const Fortran_vector &v) {
	// multiplication elementwise par colonne IN PLACE, par un vecteur de taille m = S.nblin()
	double * temp_v = v.v_ ;
	if ((int) temp_v[0] != m_) assert(false);
	int * temp_nj = ij_ + n_ - 1 ;
	double * temp_s = v_ + npnz_ - 1 ;
	int i, k, nj ;
	for (i = 1 ; i <= m_ ; i ++) {
		temp_v ++ ;
		nj = temp_nj[i] ;
		for (k = 1 ; k <= nj ; k ++) {
			temp_s ++ ; (*temp_s) *= *temp_v ;
		}
	}
	int * temp_cmi = temp_nj + m_ ;
	temp_s = v_ - 1 ; temp_v = v.v_ ;
	for (k = 1 ; k <= ntnz_ ; k ++) {
		temp_s ++ ; temp_cmi ++ ;
		(*temp_s) *= temp_v[*temp_cmi] ;
	}
	return *this;
}

Sparse_matrix & Sparse_matrix::operator/=(const double &a) {
	// S /= a , division par un scalaire in place
	double * temp = v_ - 1 ;
	for (int k = 0 ; k < ntnz_ ; k ++) {
		temp ++ ; *temp /= a ;
	}
	temp += npnz_ - ntnz_ ;
	for (int k = 0 ; k < ntnz_ ; k ++) {
		temp ++ ; *temp /= a ;
	}
	return *this;
}

Sparse_matrix & Sparse_matrix::operator/=(const Fortran_vector &v) {
	// division elementwise par colonne IN PLACE, par un vecteur de taille m = S.nblin()
	double * temp_v = v.v_ ;
	if ((int) temp_v[0] != m_) assert(false);
	int * temp_nj = ij_ + n_ - 1 ;
	double * temp_s = v_ + npnz_ - 1 ;
	int i, k, nj ;
	for (i = 1 ; i <= m_ ; i ++) {
		temp_v ++ ;
		nj = temp_nj[i] ;
		for (k = 1 ; k <= nj ; k ++) {
			temp_s ++ ; (*temp_s) /= *temp_v ;
		}
	}
	int * temp_cmi = temp_nj + m_ ;
	temp_s = v_ - 1 ; temp_v = v.v_ ;
	for (k = 1 ; k <= ntnz_ ; k ++) {
		temp_s ++ ; temp_cmi ++ ;
		(*temp_s) /= temp_v[*temp_cmi] ;
	}
	return *this;
}

void Sparse_matrix::row_elemult(const Fortran_vector &v) {
	// multiplication elementwise PAR LIGNE, IN PLACE, par un vecteur de taille n = S.nbcol()
	double * temp_v = v.v_ ;
	if ((int) temp_v[0] != n_) assert(false);
	int * temp_ni = ij_ - 1 ;
	double * temp_s = v_ - 1 ;
	int j, k, ni ;
	for (j = 1 ; j <= n_ ; j ++) {
		temp_v ++ ;
		ni = temp_ni[j] ;
		for (k = 1 ; k <= ni ; k ++) {
			temp_s ++ ; (*temp_s) *= *temp_v ;
		}
	}
	int * temp_rmj = temp_ni + m_ + n_ + npnz_ ;
	temp_s = v_ + npnz_ - 1 ; temp_v = v.v_ ;
	for (k = 1 ; k <= ntnz_ ; k ++) {
		temp_s ++ ; temp_rmj ++ ;
		(*temp_s) *= temp_v[*temp_rmj] ;
	}
}

void Sparse_matrix::row_elediv(const Fortran_vector &v) {
	// division elementwise PAR LIGNE, IN PLACE, par un vecteur de taille n = S.nbcol()
	double * temp_v = v.v_ ;
	if ((int) temp_v[0] != n_) assert(false);
	int * temp_ni = ij_ - 1 ;
	double * temp_s = v_ - 1 ;
	int j, k, ni ;
	for (j = 1 ; j <= n_ ; j ++) {
		temp_v ++ ;
		ni = temp_ni[j] ;
		for (k = 1 ; k <= ni ; k ++) {
			temp_s ++ ; (*temp_s) /= *temp_v ;
		}
	}
	int * temp_rmj = temp_ni + m_ + n_ + npnz_ ;
	temp_s = v_ + npnz_ - 1 ; temp_v = v.v_ ;
	for (k = 1 ; k <= ntnz_ ; k ++) {
		temp_s ++ ; temp_rmj ++ ;
		(*temp_s) /= temp_v[*temp_rmj] ;
	}
}

/**** Arithmetique matricielle 'inplace' composite ********************************************************************/

void Sparse_matrix::add_add(const Sparse_matrix &S1, const Sparse_matrix &S2, int ad, bool save_count) {
	if ((n_ != S1.n_) || (m_ != S1.m_) || (n_ != S2.n_) || (m_ != S2.m_)) assert(false);
	if (ad != 0) assert(false) ;	// seule la fonction  set_add(S1,S2), soit (*this) = S1+S2, est implementee
	bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, S1.ij_, S2.ij_) ;
	if (! Sparse_add_set_ij_done) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		(*this).Sparse_add_set_ij_set(S1, S2);
		(*this).add_add(S1, S2, 0, false) ;
		if (! save_count)  Sparse_add_set_ij_uncount(ij_, S1.ij_, S2.ij_) ;
	}
	else {
		int i, j, i1, i2, j1, j2 ;
//		int * nj_ = ij_ + n_ ;
        int * nj_1 = S1.ij_ + n_ ; int * nj_2 = S2.ij_ + n_  ;
		int * rmj_1 = nj_1 + m_ + S1.npnz_ ;	int * rmj_2 = nj_2 + m_ + S2.npnz_ ;	int * rmj1_max, * rmj2_max ;
		double * rmv_ = v_ + npnz_ ;	double * rmv_1 = S1.v_ + S1.npnz_ ; double * rmv_2 = S2.v_ + S2.npnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmj1_max = rmj_1 + nj_1[i] - 1 ; rmj2_max = rmj_2 + nj_2[i] - 1 ;
			j1 = (*rmj_1) ; j2 = (*rmj_2) ;
			while ((rmj_1 <= rmj1_max) && (rmj_2 <= rmj2_max)) {
				if (j1 == j2) {
					*rmv_ = (*rmv_1) + (*rmv_2) ; rmv_ ++ ; rmv_1 ++ ; rmv_2 ++ ;
					 rmj_1 ++ ; j1 = (*rmj_1) ; rmj_2 ++ ; j2 = (*rmj_2) ;
				} else
					if (j1 < j2) { *rmv_ = (*rmv_1) ; rmv_ ++ ; rmv_1 ++ ; rmj_1 ++ ; j1 = (*rmj_1) ; }
					else { *rmv_ = (*rmv_2) ; rmv_ ++ ; rmv_2 ++ ; rmj_2 ++ ; j2 = (*rmj_2) ; }
			}
			while (rmj_1 <= rmj1_max) { *rmv_ = (*rmv_1) ; rmv_ ++ ; rmv_1 ++ ; rmj_1 ++ ; }
			while (rmj_2 <= rmj2_max) { *rmv_ = (*rmv_2) ; rmv_ ++ ; rmv_2 ++ ; rmj_2 ++ ; }
		}
//		int * ni_ = ij_ ;
        int * ni_1 = S1.ij_ ; int * ni_2 = S2.ij_  ;
		int * cmi_1 = ni_1 + m_ + n_ ; int * cmi_2 = ni_2 + m_ + n_ ; int * cmi1_max, * cmi2_max ;
		double * cmv_ = v_ ; double * cmv_1 = S1.v_ ; double * cmv_2 = S2.v_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmi1_max = cmi_1 + ni_1[j] - 1 ; cmi2_max = cmi_2 + ni_2[j] - 1 ;
			i1 = (*cmi_1) ; i2 = (*cmi_2) ;
			while ((cmi_1 <= cmi1_max) && (cmi_2 <= cmi2_max)) {
				if (i1 == i2) {
					*cmv_ = (*cmv_1) + (*cmv_2) ; cmv_ ++ ; cmv_1 ++ ; cmv_2 ++ ;
					 cmi_1 ++ ; i1 = (*cmi_1) ; cmi_2 ++ ; i2 = (*cmi_2) ;
				} else
					if (i1 < i2) { *cmv_ = (*cmv_1) ; cmv_ ++ ; cmv_1 ++ ; cmi_1 ++ ; i1 = (*cmi_1) ; }
					else { *cmv_ = (*cmv_2) ; cmv_ ++ ; cmv_2 ++ ; cmi_2 ++ ; i2 = (*cmi_2) ; }
			}
			while (cmi_1 <= cmi1_max) { *cmv_ = (*cmv_1) ; cmv_ ++ ; cmv_1 ++ ; cmi_1 ++ ; }
			while (cmi_2 <= cmi2_max) { *cmv_ = (*cmv_2) ; cmv_ ++ ; cmv_2 ++ ; cmi_2 ++ ; }
		}
	}
}

void Sparse_matrix::add_add(const Sparse_matrix &S, const SpUnit_matrix &U, int ad, bool save_count) {
	if ((n_ != S.n_) || (m_ != S.m_) || (n_ != U.n_) || (m_ != U.m_)) assert(false);
	if (ad != 0) assert(false) ;	// seule la fonction  set_add(S,U), soit (*this) = S+U, est implementee
	bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, S.ij_, U.ij_) ;
	if (! Sparse_add_set_ij_done) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		Sparse_matrix T(m_, n_) ; T.ntnz_ = U.nnz_ ; T.npnz_ = U.nnz_ ; T.ij_ = U.ij_ ;
		(*this).Sparse_add_set_ij_set(S, T);
		T.ij_ = NULL ;
		(*this).add_add(S, U, 0, false) ;
		if (! save_count)  Sparse_add_set_ij_uncount(ij_, S.ij_, U.ij_) ;
	}
	else {
		int i, j, jS, jU, iS, iU ;
//		int * nj_ = ij_ + n_ ;
        int * nj_S = S.ij_ + n_ ; int * nj_U = U.ij_ + n_  ;
		int * rmj_S = nj_S + m_ + S.npnz_ ; int * rmj_U = nj_U + m_ + U.nnz_ ; int * rmjS_max, * rmjU_max ;
		double * rmv_ = v_ + npnz_ ; double * rmv_S = S.v_ + S.npnz_ ; int * rmv_U = rmj_U + 2 * U.nnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmjS_max = rmj_S + nj_S[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ;
			jS = (*rmj_S) ; jU = (*rmj_U) ;
			while ((rmj_S <= rmjS_max) && (rmj_U <= rmjU_max)) {
				if (jS == jU) {
					if (*rmv_U)  *rmv_ = (*rmv_S) + 1. ;  else  *rmv_ = (*rmv_S) - 1. ;
					rmv_ ++ ; rmv_S ++ ; rmv_U ++ ; rmj_S ++ ; jS = (*rmj_S) ; rmj_U ++ ; jU = (*rmj_U) ;
				} else
					if (jS < jU) { *rmv_ = (*rmv_S) ; rmv_ ++ ; rmv_S ++ ; rmj_S ++ ; jS = (*rmj_S) ; }
					else { if (*rmv_U)  *rmv_ = 1. ;  else  *rmv_ = -1. ;
						   rmv_ ++ ; rmv_U ++ ; rmj_U ++ ; jU = (*rmj_U) ; }
			}
			while (rmj_S <= rmjS_max) { *rmv_ = (*rmv_S) ; rmv_ ++ ; rmv_S ++ ; rmj_S ++ ; }
			while (rmj_U <= rmjU_max) { if (*rmv_U)  *rmv_ = 1. ;  else  *rmv_ = -1. ;
										rmv_ ++ ; rmv_U ++ ; rmj_U ++ ; }
		}
//		int * ni_ = ij_ ;
        int * ni_S = S.ij_ ; int * ni_U = U.ij_ ;
		int * cmi_S = ni_S + m_ + n_ ; int * cmi_U = ni_U + m_ + n_ ; int * cmiS_max, * cmiU_max ;
		double * cmv_ = v_ ; double * cmv_S = S.v_ ; int * cmv_U = cmi_U + 2 * U.nnz_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmiS_max = cmi_S + ni_S[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ;
			iS = (*cmi_S) ; iU = (*cmi_U) ;
			while ((cmi_S <= cmiS_max) && (cmi_U <= cmiU_max)) {
				if (iS == iU) {
					if (*cmv_U)  *cmv_ = (*cmv_S) + 1. ;  else  *cmv_ = (*cmv_S) - 1. ;
					cmv_ ++ ; cmv_S ++ ; cmv_U ++ ; cmi_S ++ ; iS = (*cmi_S) ; cmi_U ++ ; iU = (*cmi_U) ;
				} else
					if (iS < iU) { *cmv_ = (*cmv_S) ; cmv_ ++ ; cmv_S ++ ; cmi_S ++ ; iS = (*cmi_S) ; }
					else { if (*cmv_U)  *cmv_ = 1. ;  else  *cmv_ = -1. ;
						   cmv_ ++ ; cmv_U ++ ; cmi_U ++ ; iU = (*cmi_U) ; }
			}
			while (cmi_S <= cmiS_max) { *cmv_ = (*cmv_S) ; cmv_ ++ ; cmv_S ++ ; cmi_S ++ ; }
			while (cmi_U <= cmiU_max) { if (*cmv_U)  *cmv_ = 1. ;  else  *cmv_ = -1. ;
										cmv_ ++ ; cmv_U ++ ; cmi_U ++ ; }
		}
	}
}

void Sparse_matrix::add_add(const SpUnit_matrix &U, const Sparse_matrix &S, int ad, bool save_count) {
	(*this).add_add(S, U, ad, save_count) ;
}

void Sparse_matrix::add_add(const SpUnit_matrix &U1, const SpUnit_matrix &U2, int ad, bool save_count) {
	if ((n_ != U1.n_) || (m_ != U1.m_) || (n_ != U2.n_) || (m_ != U2.m_)) assert(false);
	if (ad != 0) assert(false) ;	// seule la fonction  set_add(U1,U2), soit (*this) = U1+U2, est implementee
	bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, U1.ij_, U2.ij_) ;
	if (! Sparse_add_set_ij_done) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		Sparse_matrix T1(m_, n_) ; T1.ij_ = U1.ij_ ; T1.ntnz_ = U1.nnz_ ; T1.npnz_ = U1.nnz_ ;
		Sparse_matrix T2(m_, n_) ; T2.ij_ = U2.ij_ ; T2.ntnz_ = U2.nnz_ ; T2.npnz_ = U2.nnz_ ;
		(*this).Sparse_add_set_ij_set(T1, T2);
		T1.ij_ = T2.ij_ = NULL ;
		(*this).add_add(U1, U2, 0, false) ;
		if (! save_count)  Sparse_add_set_ij_uncount(ij_, U1.ij_, U2.ij_) ;
	}
	else {
		int i, j, jU1, jU2, iU1, iU2 ;
//		int * nj_ = ij_ + n_ ;
        int * nj_U1 = U1.ij_ + n_ ; int * nj_U2 = U2.ij_ + n_  ;
		int * rmj_U1 = nj_U1 + m_ + U1.nnz_ ; int * rmj_U2 = nj_U2 + m_ + U2.nnz_ ;	int * rmjU1_max, * rmjU2_max ;
		double * rmv_ = v_ + npnz_ ;	int * rmv_U1 = rmj_U1 + 2 * U1.nnz_ ; int * rmv_U2 = rmj_U2 + 2 * U2.nnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmjU1_max = rmj_U1 + nj_U1[i] - 1 ;	rmjU2_max = rmj_U2 + nj_U2[i] - 1 ;
			jU1 = (*rmj_U1) ; jU2 = (*rmj_U2) ;
			while ((rmj_U1 <= rmjU1_max) && (rmj_U2 <= rmjU2_max)) {
				if (jU1 == jU2) {
					if (*rmv_U2) {if (*rmv_U1) *rmv_ = 2.; else *rmv_ =  0.;}
					else		 {if (*rmv_U1) *rmv_ = 0.; else *rmv_ = -2.;}
					rmv_ ++ ; rmv_U1 ++ ; rmv_U2 ++ ; rmj_U1 ++ ; jU1 = (*rmj_U1) ; rmj_U2 ++ ; jU2 = (*rmj_U2) ;
				} else
					if (jU1 < jU2) { if (*rmv_U1)  *rmv_ = 1. ;  else  *rmv_ = -1. ;
									 rmv_ ++ ; rmv_U1 ++ ; rmj_U1 ++ ; jU1 = (*rmj_U1) ; }
					else		   { if (*rmv_U2)  *rmv_ = 1. ;  else  *rmv_ = -1. ;
									 rmv_ ++ ; rmv_U2 ++ ; rmj_U2 ++ ; jU2 = (*rmj_U2) ; }
			}
			while (rmj_U1 <= rmjU1_max) { if (*rmv_U1)  *rmv_ = 1. ;  else  *rmv_ = -1. ;
										  rmv_ ++ ; rmv_U1 ++ ; rmj_U1 ++ ; }
			while (rmj_U2 <= rmjU2_max) { if (*rmv_U2)  *rmv_ = 1. ;  else  *rmv_ = -1. ;
										  rmv_ ++ ; rmv_U2 ++ ; rmj_U2 ++ ; }
		}
//		int * ni_ = ij_ ;
        int * ni_U1 = U1.ij_ ; int * ni_U2 = U2.ij_  ;
		int * cmi_U1 = ni_U1 + m_ + n_ ; int * cmi_U2 = ni_U2 + m_ + n_ ; int * cmiU1_max, * cmiU2_max ;
		double * cmv_ = v_ ; int * cmv_U1 = cmi_U1 + 2 * U1.nnz_ ; int * cmv_U2 = cmi_U2 + 2 * U2.nnz_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmiU1_max = cmi_U1 + ni_U1[j] - 1 ; cmiU2_max = cmi_U2 + ni_U2[j] - 1 ;
			iU1 = (*cmi_U1) ; iU2 = (*cmi_U2) ;
			while ((cmi_U1 <= cmiU1_max) && (cmi_U2 <= cmiU2_max)) {
				if (iU1 == iU2) {
					if (*cmv_U2) {if (*cmv_U1) *cmv_ = 2.; else *cmv_ =  0.;}
					else		 {if (*cmv_U1) *cmv_ = 0.; else *cmv_ = -2.;}
					cmv_ ++ ; cmv_U1 ++ ; cmv_U2 ++ ; cmi_U1 ++ ; iU1 = (*cmi_U1) ; cmi_U2 ++ ; iU2 = (*cmi_U2) ;
				} else
					if (iU1 < iU2) { if (*cmv_U1)  *cmv_ = 1. ;  else  *cmv_ = -1. ;
									 cmv_ ++ ; cmv_U1 ++ ; cmi_U1 ++ ; iU1 = (*cmi_U1) ; }
					else		   { if (*cmv_U2)  *cmv_ = 1. ;  else  *cmv_ = -1. ;
									 cmv_ ++ ; cmv_U2 ++ ; cmi_U2 ++ ; iU2 = (*cmi_U2) ; }
			}
			while (cmi_U1 <= cmiU1_max) { if (*cmv_U1)  *cmv_ = 1. ;  else  *cmv_ = -1. ;
										  cmv_ ++ ; cmv_U1 ++ ; cmi_U1 ++ ; }
			while (cmi_U2 <= cmiU2_max) { if (*cmv_U2)  *cmv_ = 1. ;  else  *cmv_ = -1. ;
										  cmv_ ++ ; cmv_U2 ++ ; cmi_U2 ++ ; }
		}
	}
}

void Sparse_matrix::add_sub(const Sparse_matrix &S1, const Sparse_matrix &S2, int ad, bool save_count) {
	if ((n_ != S1.n_) || (m_ != S1.m_) || (n_ != S2.n_) || (m_ != S2.m_)) assert(false);
	if (ad != 0) assert(false) ;	// seule la fonction  set_sub(S1,S2), soit (*this) = S1-S2, est implementee
	bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, S1.ij_, S2.ij_) ;
	if (! Sparse_add_set_ij_done) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		(*this).Sparse_add_set_ij_set(S1, S2);
		(*this).add_sub(S1, S2, 0, false) ;
		if (! save_count)  Sparse_add_set_ij_uncount(ij_, S1.ij_, S2.ij_) ;
	}
	else {
		int i, j, i1, i2, j1, j2 ;
//		int * nj_ = ij_ + n_ ;
        int * nj_1 = S1.ij_ + n_ ; int * nj_2 = S2.ij_ + n_  ;
		int * rmj_1 = nj_1 + m_ + S1.npnz_ ; int * rmj_2 = nj_2 + m_ + S2.npnz_ ;	int * rmj1_max, * rmj2_max ;
		double * rmv_ = v_ + npnz_ ;	double * rmv_1 = S1.v_ + S1.npnz_ ; double * rmv_2 = S2.v_ + S2.npnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmj1_max = rmj_1 + nj_1[i] - 1 ; rmj2_max = rmj_2 + nj_2[i] - 1 ;
			j1 = (*rmj_1) ; j2 = (*rmj_2) ;
			while ((rmj_1 <= rmj1_max) && (rmj_2 <= rmj2_max)) {
				if (j1 == j2) {
					*rmv_ = (*rmv_1) - (*rmv_2) ; rmv_ ++ ; rmv_1 ++ ; rmv_2 ++ ;
					 rmj_1 ++ ; j1 = (*rmj_1) ; rmj_2 ++ ; j2 = (*rmj_2) ;
				} else
					if (j1 < j2) { *rmv_ = (*rmv_1) ; rmv_ ++ ; rmv_1 ++ ; rmj_1 ++ ; j1 = (*rmj_1) ; }
					else { *rmv_ = -(*rmv_2) ; rmv_ ++ ; rmv_2 ++ ; rmj_2 ++ ; j2 = (*rmj_2) ; }
			}
			while (rmj_1 <= rmj1_max) { *rmv_ = (*rmv_1) ; rmv_ ++ ; rmv_1 ++ ; rmj_1 ++ ; }
			while (rmj_2 <= rmj2_max) { *rmv_ = -(*rmv_2) ; rmv_ ++ ; rmv_2 ++ ; rmj_2 ++ ; }
		}
//		int * ni_ = ij_ ;
        int * ni_1 = S1.ij_ ; int * ni_2 = S2.ij_ ;
		int * cmi_1 = ni_1 + m_ + n_ ; int * cmi_2 = ni_2 + m_ + n_ ; int * cmi1_max, * cmi2_max ;
		double * cmv_ = v_ ; double * cmv_1 = S1.v_ ; double * cmv_2 = S2.v_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmi1_max = cmi_1 + ni_1[j] - 1 ; cmi2_max = cmi_2 + ni_2[j] - 1 ;
			i1 = (*cmi_1) ; i2 = (*cmi_2) ;
			while ((cmi_1 <= cmi1_max) && (cmi_2 <= cmi2_max)) {
				if (i1 == i2) {
					*cmv_ = (*cmv_1) - (*cmv_2) ; cmv_ ++ ; cmv_1 ++ ; cmv_2 ++ ;
					 cmi_1 ++ ; i1 = (*cmi_1) ; cmi_2 ++ ; i2 = (*cmi_2) ;
				} else
					if (i1 < i2) { *cmv_ = (*cmv_1) ; cmv_ ++ ; cmv_1 ++ ; cmi_1 ++ ; i1 = (*cmi_1) ; }
					else { *cmv_ = -(*cmv_2) ; cmv_ ++ ; cmv_2 ++ ; cmi_2 ++ ; i2 = (*cmi_2) ; }
			}
			while (cmi_1 <= cmi1_max) { *cmv_ = (*cmv_1) ; cmv_ ++ ; cmv_1 ++ ; cmi_1 ++ ; }
			while (cmi_2 <= cmi2_max) { *cmv_ = -(*cmv_2) ; cmv_ ++ ; cmv_2 ++ ; cmi_2 ++ ; }
		}
	}
}

void Sparse_matrix::add_sub(const SpUnit_matrix &U, const Sparse_matrix &S, int ad, bool save_count) {
	if ((n_ != S.n_) || (m_ != S.m_) || (n_ != U.n_) || (m_ != U.m_)) assert(false);
	if (ad != 0) assert(false) ;	// seule la fonction  set_sub(U,S), soit (*this) = U-S, est implementee
	bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, S.ij_, U.ij_) ;
	if (! Sparse_add_set_ij_done) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		Sparse_matrix T(m_, n_) ; T.ntnz_ = U.nnz_ ; T.npnz_ = U.nnz_ ; T.ij_ = U.ij_ ;
		(*this).Sparse_add_set_ij_set(T, S);
		T.ij_ = NULL ;
		(*this).add_sub(U, S, 0, false) ;
		if (! save_count)  Sparse_add_set_ij_uncount(ij_, S.ij_, U.ij_) ;
	}
	else {
		int i, j, jS, jU, iS, iU ;
//		int * nj_ = ij_ + n_ ;
        int * nj_S = S.ij_ + n_ ; int * nj_U = U.ij_ + n_ ;
		int * rmj_S = nj_S + m_ + S.npnz_ ; int * rmj_U = nj_U + m_ + U.nnz_ ; int * rmjS_max, * rmjU_max ;
		double * rmv_ = v_ + npnz_ ;	double * rmv_S = S.v_ + S.npnz_ ; int * rmv_U = rmj_U + 2 * U.nnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmjS_max = rmj_S + nj_S[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ;
			jS = (*rmj_S) ; jU = (*rmj_U) ;
			while ((rmj_S <= rmjS_max) && (rmj_U <= rmjU_max)) {
				if (jS == jU) {
					if (*rmv_U)  *rmv_ = 1. - (*rmv_S);  else  *rmv_ = -1. - (*rmv_S) ;
					rmv_ ++ ; rmv_S ++ ; rmv_U ++ ; rmj_S ++ ; jS = (*rmj_S) ; rmj_U ++ ; jU = (*rmj_U) ;
				} else
					if (jS < jU) { *rmv_ = -(*rmv_S) ; rmv_ ++ ; rmv_S ++ ; rmj_S ++ ; jS = (*rmj_S) ; }
					else { if (*rmv_U)  *rmv_ = 1. ;  else  *rmv_ = -1. ;
						   rmv_ ++ ; rmv_U ++ ; rmj_U ++ ; jU = (*rmj_U) ; }
			}
			while (rmj_S <= rmjS_max) { *rmv_ = -(*rmv_S) ; rmv_ ++ ; rmv_S ++ ; rmj_S ++ ; }
			while (rmj_U <= rmjU_max) { if (*rmv_U)  *rmv_ = 1. ;  else  *rmv_ = -1. ;
										rmv_ ++ ; rmv_U ++ ; rmj_U ++ ; }
		}
//		int * ni_ = ij_ ;
        int * ni_S = S.ij_ ; int * ni_U = U.ij_ ;
		int * cmi_S = ni_S + m_ + n_ ; int * cmi_U = ni_U + m_ + n_ ; int * cmiS_max, * cmiU_max ;
		double * cmv_ = v_ ; double * cmv_S = S.v_ ; int * cmv_U = cmi_U + 2 * U.nnz_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmiS_max = cmi_S + ni_S[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ;
			iS = (*cmi_S) ; iU = (*cmi_U) ;
			while ((cmi_S <= cmiS_max) && (cmi_U <= cmiU_max)) {
				if (iS == iU) {
					if (*cmv_U)  *cmv_ = 1. - (*cmv_S) ;  else  *cmv_ = -1. - (*cmv_S) ;
					cmv_ ++ ; cmv_S ++ ; cmv_U ++ ; cmi_S ++ ; iS = (*cmi_S) ; cmi_U ++ ; iU = (*cmi_U) ;
				} else
					if (iS < iU) { *cmv_ = -(*cmv_S) ; cmv_ ++ ; cmv_S ++ ; cmi_S ++ ; iS = (*cmi_S) ; }
					else { if (*cmv_U)  *cmv_ = 1. ;  else  *cmv_ = -1. ;
						   cmv_ ++ ; cmv_U ++ ; cmi_U ++ ; iU = (*cmi_U) ; }
			}
			while (cmi_S <= cmiS_max) { *cmv_ = -(*cmv_S) ; cmv_ ++ ; cmv_S ++ ; cmi_S ++ ; }
			while (cmi_U <= cmiU_max) { if (*cmv_U)  *cmv_ = 1. ;  else  *cmv_ = -1. ;
										cmv_ ++ ; cmv_U ++ ; cmi_U ++ ; }
		}
	}
}

void Sparse_matrix::add_sub(const Sparse_matrix &S, const SpUnit_matrix &U, int ad, bool save_count) {
	if ((n_ != S.n_) || (m_ != S.m_) || (n_ != U.n_) || (m_ != U.m_)) assert(false);
	if (ad != 0) assert(false) ;	// seule la fonction  set_sub(S,U), soit (*this) = S-U, est implementee
	bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, S.ij_, U.ij_) ;
	if (! Sparse_add_set_ij_done) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		Sparse_matrix T(m_, n_) ; T.ntnz_ = U.nnz_ ; T.npnz_ = U.nnz_ ; T.ij_ = U.ij_ ;
		(*this).Sparse_add_set_ij_set(S, T);
		T.ij_ = NULL ;
		(*this).add_sub(S, U, 0, false) ;
		if (! save_count)  Sparse_add_set_ij_uncount(ij_, S.ij_, U.ij_) ;
	}
	else {
		int i, j, jS, jU, iS, iU ;
//		int * nj_ = ij_ + n_ ;
        int * nj_S = S.ij_ + n_ ; int * nj_U = U.ij_ + n_  ;
		int * rmj_S = nj_S + m_ + S.npnz_ ; int * rmj_U = nj_U + m_ + U.nnz_ ; int * rmjS_max, * rmjU_max ;
		double * rmv_ = v_ + npnz_ ;	double * rmv_S = S.v_ + S.npnz_ ; int * rmv_U = rmj_U + 2 * U.nnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmjS_max = rmj_S + nj_S[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ;
			jS = (*rmj_S) ; jU = (*rmj_U) ;
			while ((rmj_S <= rmjS_max) && (rmj_U <= rmjU_max)) {
				if (jS == jU) {
					if (*rmv_U)  *rmv_ = (*rmv_S) - 1. ;  else  *rmv_ = (*rmv_S) + 1. ;
					rmv_ ++ ; rmv_S ++ ; rmv_U ++ ; rmj_S ++ ; jS = (*rmj_S) ; rmj_U ++ ; jU = (*rmj_U) ;
				} else
					if (jS < jU) { *rmv_ = (*rmv_S) ; rmv_ ++ ; rmv_S ++ ; rmj_S ++ ; jS = (*rmj_S) ; }
					else { if (*rmv_U)  *rmv_ = -1. ;  else  *rmv_ = 1. ;
						   rmv_ ++ ; rmv_U ++ ; rmj_U ++ ; jU = (*rmj_U) ; }
			}
			while (rmj_S <= rmjS_max) { *rmv_ = (*rmv_S) ; rmv_ ++ ; rmv_S ++ ; rmj_S ++ ; }
			while (rmj_U <= rmjU_max) { if (*rmv_U)  *rmv_ = -1. ;  else  *rmv_ = 1. ;
										rmv_ ++ ; rmv_U ++ ; rmj_U ++ ; }
		}
//		int * ni_ = ij_ ;
        int * ni_S = S.ij_ ; int * ni_U = U.ij_  ;
		int * cmi_S = ni_S + m_ + n_ ; int * cmi_U = ni_U + m_ + n_ ; int * cmiS_max, * cmiU_max ;
		double * cmv_ = v_ ; double * cmv_S = S.v_ ; int * cmv_U = cmi_U + 2 * U.nnz_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmiS_max = cmi_S + ni_S[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ;
			iS = (*cmi_S) ; iU = (*cmi_U) ;
			while ((cmi_S <= cmiS_max) && (cmi_U <= cmiU_max)) {
				if (iS == iU) {
					if (*cmv_U)  *cmv_ = (*cmv_S) - 1. ;  else  *cmv_ = (*cmv_S) + 1. ;
					cmv_ ++ ; cmv_S ++ ; cmv_U ++ ; cmi_S ++ ; iS = (*cmi_S) ; cmi_U ++ ; iU = (*cmi_U) ;
				} else
					if (iS < iU) { *cmv_ = (*cmv_S) ; cmv_ ++ ; cmv_S ++ ; cmi_S ++ ; iS = (*cmi_S) ; }
					else { if (*cmv_U)  *cmv_ = -1. ;  else  *cmv_ = 1. ;
						   cmv_ ++ ; cmv_U ++ ; cmi_U ++ ; iU = (*cmi_U) ; }
			}
			while (cmi_S <= cmiS_max) { *cmv_ = (*cmv_S) ; cmv_ ++ ; cmv_S ++ ; cmi_S ++ ; }
			while (cmi_U <= cmiU_max) { if (*cmv_U)  *cmv_ = -1. ;  else  *cmv_ = 1. ;
										cmv_ ++ ; cmv_U ++ ; cmi_U ++ ; }
		}
	}
}

void Sparse_matrix::add_sub(const SpUnit_matrix &U1, const SpUnit_matrix &U2, int ad, bool save_count) {
	if ((n_ != U1.n_) || (m_ != U1.m_) || (n_ != U2.n_) || (m_ != U2.m_)) assert(false);
	if (ad != 0) assert(false) ;	// seule la fonction  set_sub(U1,U2), soit (*this) = U1-U2, est implementee
	bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, U1.ij_, U2.ij_) ;
	if (! Sparse_add_set_ij_done) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		Sparse_matrix T1(m_, n_) ; T1.ij_ = U1.ij_ ; T1.ntnz_ = U1.nnz_ ; T1.npnz_ = U1.nnz_ ;
		Sparse_matrix T2(m_, n_) ; T2.ij_ = U2.ij_ ; T2.ntnz_ = U2.nnz_ ; T2.npnz_ = U2.nnz_ ;
		(*this).Sparse_add_set_ij_set(T1, T2);
		T1.ij_ = T2.ij_ = NULL ;
		(*this).add_sub(U1, U2, 0, false) ;
		if (! save_count)  Sparse_add_set_ij_uncount(ij_, U1.ij_, U2.ij_) ;
	}
	else {
		int i, j, jU1, jU2, iU1, iU2 ;
        int * nj_U1 = U1.ij_ + n_ ; int * nj_U2 = U2.ij_ + n_ ;
		int * rmj_U1 = nj_U1 + m_ + U1.nnz_ ; int * rmj_U2 = nj_U2 + m_ + U2.nnz_ ;	int * rmjU1_max, * rmjU2_max ;
		double * rmv_ = v_ + npnz_ ;	int * rmv_U1 = rmj_U1 + 2 * U1.nnz_ ; int * rmv_U2 = rmj_U2 + 2 * U2.nnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmjU1_max = rmj_U1 + nj_U1[i] - 1 ;	rmjU2_max = rmj_U2 + nj_U2[i] - 1 ;
			jU1 = (*rmj_U1) ; jU2 = (*rmj_U2) ;
			while ((rmj_U1 <= rmjU1_max) && (rmj_U2 <= rmjU2_max)) {
				if (jU1 == jU2) {
					if (*rmv_U2) {if (*rmv_U1) *rmv_ = 0.; else *rmv_ = -2.;}
					else		 {if (*rmv_U1) *rmv_ = 2.; else *rmv_ =  0.;}
					rmv_ ++ ; rmv_U1 ++ ; rmv_U2 ++ ; rmj_U1 ++ ; jU1 = (*rmj_U1) ; rmj_U2 ++ ; jU2 = (*rmj_U2) ;
				} else
					if (jU1 < jU2) { if (*rmv_U1)  *rmv_ = 1. ;  else  *rmv_ = -1. ;
									 rmv_ ++ ; rmv_U1 ++ ; rmj_U1 ++ ; jU1 = (*rmj_U1) ; }
					else		   { if (*rmv_U2)  *rmv_ = -1. ;  else  *rmv_ = 1. ;
									 rmv_ ++ ; rmv_U2 ++ ; rmj_U2 ++ ; jU2 = (*rmj_U2) ; }
			}
			while (rmj_U1 <= rmjU1_max) { if (*rmv_U1)  *rmv_ = 1. ;  else  *rmv_ = -1. ;
										  rmv_ ++ ; rmv_U1 ++ ; rmj_U1 ++ ; }
			while (rmj_U2 <= rmjU2_max) { if (*rmv_U2)  *rmv_ = -1. ;  else  *rmv_ = 1. ;
										  rmv_ ++ ; rmv_U2 ++ ; rmj_U2 ++ ; }
		}
        int * ni_U1 = U1.ij_ ; int * ni_U2 = U2.ij_ ;
		int * cmi_U1 = ni_U1 + m_ + n_ ; int * cmi_U2 = ni_U2 + m_ + n_ ; int * cmiU1_max, * cmiU2_max ;
		double * cmv_ = v_ ; int * cmv_U1 = cmi_U1 + 2 * U1.nnz_ ; int * cmv_U2 = cmi_U2 + 2 * U2.nnz_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmiU1_max = cmi_U1 + ni_U1[j] - 1 ; cmiU2_max = cmi_U2 + ni_U2[j] - 1 ;
			iU1 = (*cmi_U1) ; iU2 = (*cmi_U2) ;
			while ((cmi_U1 <= cmiU1_max) && (cmi_U2 <= cmiU2_max)) {
				if (iU1 == iU2) {
					if (*cmv_U2) {if (*cmv_U1) *cmv_ = 0.; else *cmv_ = -2.;}
					else		 {if (*cmv_U1) *cmv_ = 2.; else *cmv_ =  0.;}
					cmv_ ++ ; cmv_U1 ++ ; cmv_U2 ++ ; cmi_U1 ++ ; iU1 = (*cmi_U1) ; cmi_U2 ++ ; iU2 = (*cmi_U2) ;
				} else
					if (iU1 < iU2) { if (*cmv_U1)  *cmv_ = 1. ;  else  *cmv_ = -1. ;
									 cmv_ ++ ; cmv_U1 ++ ; cmi_U1 ++ ; iU1 = (*cmi_U1) ; }
					else		   { if (*cmv_U2)  *cmv_ = -1. ;  else  *cmv_ = 1. ;
									 cmv_ ++ ; cmv_U2 ++ ; cmi_U2 ++ ; iU2 = (*cmi_U2) ; }
			}
			while (cmi_U1 <= cmiU1_max) { if (*cmv_U1)  *cmv_ = 1. ;  else  *cmv_ = -1. ;
										  cmv_ ++ ; cmv_U1 ++ ; cmi_U1 ++ ; }
			while (cmi_U2 <= cmiU2_max) { if (*cmv_U2)  *cmv_ = -1. ;  else  *cmv_ = 1. ;
										  cmv_ ++ ; cmv_U2 ++ ; cmi_U2 ++ ; }
		}
	}
}

void Sparse_matrix::set_matmult(const Sparse_matrix &S1, const Sparse_matrix &S2, Fortran_vector* temp_v, bool save_count) {
	// n'autorise (*this = S1) ou (*this = S2) que si nnz_ reste inchange.
	//   si temp_v est fourni, il doit pointer sur un vecteur  de taille  n = S1.n_  qui servira de tampon.
	int m = S1.m_, n = S1.n_, p = S2.n_ ;
	if ((S2.m_ != n) || (m_ != m) || (n_ != p)) assert(false);
	double * v ; if (temp_v == NULL)   v = new double[n+1] ;
				 else { v = (*temp_v).v_ ; if ((int)v[0] != n) assert(false); }
	bool Sparse_matmult_set_ij_done = Sparse_matmult_set_ij_check(ij_, S1.ij_, S2.ij_) ;
	if (! Sparse_matmult_set_ij_done) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		(*this).Sparse_matmult_set_ij_set(S1, S2, v);
		Fortran_vector temp_v2 ; temp_v2.v_ = v ; v[0] = n ;
		(*this).set_matmult(S1, S2, &temp_v2, false) ;
		temp_v2.v_ = NULL ;
		if (! save_count)  Sparse_matmult_set_ij_uncount(ij_, S1.ij_, S2.ij_) ;
	}
	else {
		int i, j, i1, j2, k ;
		int * cmi_2 = S2.ij_ + p + n ; double * cmv_2 = S2.v_ ; int * rmj_1 ; double * rmv_1 ;
		int * ni_ = ij_ - 1 ; int * ni_2 = S2.ij_ - 1 ; int * nj_1 = S1.ij_ + n - 1 ;
		int * cmi_ = ij_ + m_ + n_ ; double * cmv_ = v_ ;
		int * rmj1_max, * cmi2_max, * rmj_max, * cmi_max ;
		for (j = 1 ; j <= n_ ; j ++) {
			for (k = 1 ; k <= n ; k++)	 v[k] = 0.;
			cmi2_max = cmi_2 + ni_2[j] - 1 ;
			while (cmi_2 <= cmi2_max)	{ v[*cmi_2] = *cmv_2 ; cmi_2 ++ ; cmv_2 ++ ; }
			rmj_1 = S1.ij_ + m + n + S1.npnz_ ; rmv_1 = S1.v_ + S1.npnz_ ;
			i = 0 ; cmi_max = cmi_ + ni_[j] - 1 ;
			while (cmi_ <= cmi_max) {
				i1 = i+1 ; i = *cmi_ ; *cmv_ = 0 ;
				for (k = i1 ; k < i ; k ++) { rmj_1 += nj_1[k] ; rmv_1 += nj_1[k] ; }
				rmj1_max = rmj_1 + nj_1[i] - 1 ;
				while (rmj_1 <= rmj1_max) {	(*cmv_) += (*rmv_1) * v[*rmj_1] ; rmv_1 ++ ; rmj_1 ++ ;	}
				cmv_ ++ ; cmi_ ++ ;
			}
		}
		rmj_1 = S1.ij_ + m + n + S1.npnz_ ; rmv_1 = S1.v_ + S1.npnz_ ;
		int * nj_ = ij_ - 1 + n_ ; int * rmj_ = ij_ + m_ + n_ + npnz_ ; double * rmv_ = v_ + npnz_ ;
		for (i = 1 ; i <= m_ ; i ++) {
			for (k = 1 ; k <= n ; k++)	 v[k] = 0.;
			rmj1_max = rmj_1 + nj_1[i] - 1 ;
			while (rmj_1 <= rmj1_max)	{ v[*rmj_1] = *rmv_1 ; rmj_1 ++ ; rmv_1 ++ ; }
			cmi_2 = S2.ij_ + p + n ; cmv_2 = S2.v_ ;
			j = 0 ; rmj_max = rmj_ + nj_[i] - 1 ;
			while (rmj_ <= rmj_max) {
				j2 = j+1 ; j = *rmj_ ; *rmv_ = 0 ;
				for (k = j2 ; k < j ; k ++) { cmi_2 += ni_2[k] ; cmv_2 += ni_2[k] ; }
				cmi2_max = cmi_2 + ni_2[j] - 1 ;
				while (cmi_2 <= cmi2_max) {	(*rmv_) += (*cmv_2) * v[*cmi_2] ; cmv_2 ++ ; cmi_2 ++ ;	}
				rmv_ ++ ; rmj_ ++ ;
			}
		}
	}
	if (temp_v == NULL) delete [] v ;
}

void Sparse_matrix::set_matmult(const SpUnit_matrix &U1, const Sparse_matrix &S2, Fortran_vector* temp_v, bool save_count) {
	// n'autorise (*this = U1) ou (*this = S2) que si nnz_ reste inchange.
	//   si temp_v est fourni, il doit pointer sur un vecteur  de taille  n = U1.n_  qui servira de tampon.
	int m = U1.m_, n = U1.n_, p = S2.n_ ;
	if ((S2.m_ != n) || (m_ != m) || (n_ != p)) assert(false);
	double * v ; if (temp_v == NULL)   v = new double[n+1] ;
				 else { v = (*temp_v).v_ ; if ((int)v[0] != n) assert(false); }
	bool Sparse_matmult_set_ij_done = Sparse_matmult_set_ij_check(ij_, U1.ij_, S2.ij_) ;
	if (! Sparse_matmult_set_ij_done) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		Sparse_matrix T1(m, n) ; T1.ij_ = U1.ij_ ; T1.ntnz_ = U1.nnz_ ; T1.npnz_ = U1.nnz_ ;
		(*this).Sparse_matmult_set_ij_set(T1, S2, v);
		T1.ij_ = NULL ;
		Fortran_vector temp_v2 ; temp_v2.v_ = v ; v[0] = n ;
		(*this).set_matmult(U1, S2, &temp_v2, false) ;
		temp_v2.v_ = NULL ;
		if (! save_count)  Sparse_matmult_set_ij_uncount(ij_, U1.ij_, S2.ij_) ;
	}
	else {
		int i, j, i1, j2, k ;
		int * cmi_2 = S2.ij_ + p + n ; double * cmv_2 = S2.v_ ; int * rmj_1 ; int * rmv_1 ;
		int * ni_ = ij_ - 1 ; int * ni_2 = S2.ij_ - 1 ; int * nj_1 = U1.ij_ + n - 1 ;
		int * cmi_ = ij_ + m_ + n_ ; double * cmv_ = v_ ;
		int * rmj1_max, * cmi2_max, * rmj_max, * cmi_max ;
		for (j = 1 ; j <= n_ ; j ++) {
			for (k = 1 ; k <= n ; k++)	 v[k] = 0.;
			cmi2_max = cmi_2 + ni_2[j] - 1 ;
			while (cmi_2 <= cmi2_max)	{ v[*cmi_2] = *cmv_2 ; cmi_2 ++ ; cmv_2 ++ ; }
			rmj_1 = U1.ij_ + m + n + U1.nnz_ ; rmv_1 = rmj_1 + U1.nnz_ + U1.nnz_ ;
			i = 0 ; cmi_max = cmi_ + ni_[j] - 1 ;
			while (cmi_ <= cmi_max) {
				i1 = i+1 ; i = *cmi_ ; *cmv_ = 0 ;
				for (k = i1 ; k < i ; k ++) { rmj_1 += nj_1[k] ; rmv_1 += nj_1[k] ; }
				rmj1_max = rmj_1 + nj_1[i] - 1 ;
				while (rmj_1 <= rmj1_max)
					{ if (*rmv_1)  (*cmv_) += v[*rmj_1] ;  else  (*cmv_) -= v[*rmj_1] ;
					  rmv_1 ++ ; rmj_1 ++ ;	}
				cmv_ ++ ; cmi_ ++ ;
			}
		}
		rmj_1 = U1.ij_ + m + n + U1.nnz_ ; rmv_1 = rmj_1 + U1.nnz_ + U1.nnz_ ;
		int * nj_ = ij_ - 1 + n_ ; int * rmj_ = ij_ + m_ + n_ + npnz_ ; double * rmv_ = v_ + npnz_ ;
		for (i = 1 ; i <= m_ ; i ++) {
			for (k = 1 ; k <= n ; k++)	 v[k] = 0.;
			rmj1_max = rmj_1 + nj_1[i] - 1 ;
			while (rmj_1 <= rmj1_max)	{ if (*rmv_1) v[*rmj_1] = 1 ; else  v[*rmj_1] = -1 ;
											rmj_1 ++ ; rmv_1 ++ ; }
			cmi_2 = S2.ij_ + p + n ; cmv_2 = S2.v_ ;
			j = 0 ; rmj_max = rmj_ + nj_[i] - 1 ;
			while (rmj_ <= rmj_max) {
				j2 = j+1 ; j = *rmj_ ; *rmv_ = 0 ;
				for (k = j2 ; k < j ; k ++) { cmi_2 += ni_2[k] ; cmv_2 += ni_2[k] ; }
				cmi2_max = cmi_2 + ni_2[j] - 1 ;
				while (cmi_2 <= cmi2_max) {
					if (v[*cmi_2])  { if (v[*cmi_2] == 1)  (*rmv_) += (*cmv_2) ;  else  (*rmv_) -= (*cmv_2) ; }
					cmv_2 ++ ; cmi_2 ++ ;
				}
				rmv_ ++ ; rmj_ ++ ;
			}
		}
	}
	if (temp_v == NULL) delete [] v ;
}

void Sparse_matrix::set_matmult(const Sparse_matrix &S1, const SpUnit_matrix &U2, Fortran_vector* temp_v, bool save_count) {
	// n'autorise (*this = S1) ou (*this = U2) que si nnz_ reste inchange.
	//   si temp_v est fourni, il doit pointer sur un vecteur  de taille  n = S1.n_  qui servira de tampon.
	int m = S1.m_, n = S1.n_, p = U2.n_ ;
	if ((U2.m_ != n) || (m_ != m) || (n_ != p)) assert(false);
	double * v ; if (temp_v == NULL)   v = new double[n+1] ;
				 else { v = (*temp_v).v_ ; if ((int)v[0] != n) assert(false); }
	bool Sparse_matmult_set_ij_done = Sparse_matmult_set_ij_check(ij_, S1.ij_, U2.ij_) ;
	if (! Sparse_matmult_set_ij_done) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		Sparse_matrix T2(n, p) ; T2.ij_ = U2.ij_ ; T2.ntnz_ = U2.nnz_ ; T2.npnz_ = U2.nnz_ ;
		(*this).Sparse_matmult_set_ij_set(S1, T2, v);
		T2.ij_ = NULL ;
		Fortran_vector temp_v2 ; temp_v2.v_ = v ; v[0] = n ;
		(*this).set_matmult(S1, U2, &temp_v2, false) ;
		temp_v2.v_ = NULL ;
		if (! save_count)  Sparse_matmult_set_ij_uncount(ij_, S1.ij_, U2.ij_) ;
	}
	else {
		int i, j, i1, j2, k ;
		int * cmi_2 = U2.ij_ + p + n ; int * cmv_2 = cmi_2 + U2.nnz_ + U2.nnz_ ; int * rmj_1 ; double * rmv_1 ;
		int * ni_ = ij_ - 1 ; int * ni_2 = U2.ij_ - 1 ; int * nj_1 = S1.ij_ + n - 1 ;
		int * cmi_ = ij_ + m_ + n_ ; double * cmv_ = v_ ;
		int * rmj1_max, * cmi2_max, * rmj_max, * cmi_max ;
		for (j = 1 ; j <= n_ ; j ++) {
			for (k = 1 ; k <= n ; k++)	 v[k] = 0.;
			cmi2_max = cmi_2 + ni_2[j] - 1 ;
			while (cmi_2 <= cmi2_max)	{ if (*cmv_2)  v[*cmi_2] = 1 ; else v[*cmi_2] = -1 ;
											cmi_2 ++ ; cmv_2 ++ ; }
			rmj_1 = S1.ij_ + m + n + S1.npnz_ ; rmv_1 = S1.v_ + S1.npnz_ ;
			i = 0 ; cmi_max = cmi_ + ni_[j] - 1 ;
			while (cmi_ <= cmi_max) {
				i1 = i+1 ; i = *cmi_ ; *cmv_ = 0 ;
				for (k = i1 ; k < i ; k ++) { rmj_1 += nj_1[k] ; rmv_1 += nj_1[k] ; }
				rmj1_max = rmj_1 + nj_1[i] - 1 ;
				while (rmj_1 <= rmj1_max) {
					if (v[*rmj_1])  { if (v[*rmj_1] == 1)  (*cmv_) += (*rmv_1) ; else  (*cmv_) -= (*rmv_1) ; }
					rmv_1 ++ ; rmj_1 ++ ;
				}
				cmv_ ++ ; cmi_ ++ ;
			}
		}
		rmj_1 = S1.ij_ + m + n + S1.npnz_ ; rmv_1 = S1.v_ + S1.npnz_ ;
		int * nj_ = ij_ - 1 + n_ ; int * rmj_ = ij_ + m_ + n_ + npnz_ ; double * rmv_ = v_ + npnz_ ;
		for (i = 1 ; i <= m_ ; i ++) {
			for (k = 1 ; k <= n ; k++)	 v[k] = 0.;
			rmj1_max = rmj_1 + nj_1[i] - 1 ;
			while (rmj_1 <= rmj1_max)	{ v[*rmj_1] = *rmv_1 ; rmj_1 ++ ; rmv_1 ++ ; }
			cmi_2 = U2.ij_ + p + n ; cmv_2 = cmi_2 + U2.nnz_ + U2.nnz_ ;
			j = 0 ; rmj_max = rmj_ + nj_[i] - 1 ;
			while (rmj_ <= rmj_max) {
				j2 = j+1 ; j = *rmj_ ; *rmv_ = 0 ;
				for (k = j2 ; k < j ; k ++) { cmi_2 += ni_2[k] ; cmv_2 += ni_2[k] ; }
				cmi2_max = cmi_2 + ni_2[j] - 1 ;
				while (cmi_2 <= cmi2_max) {	if (*cmv_2) (*rmv_) += v[*cmi_2] ; else (*rmv_) -= v[*cmi_2] ;
											cmv_2 ++ ; cmi_2 ++ ;	}
				rmv_ ++ ; rmj_ ++ ;
			}
		}
	}
	if (temp_v == NULL) delete [] v ;
}

void Sparse_matrix::set_matmult(const SpUnit_matrix &U1, const SpUnit_matrix &U2, Index_vector* temp_v, bool save_count) {
	// n'autorise (*this = U1) ou (*this = U2) que si nnz_ reste inchange.
	//   si temp_v est fourni, il doit pointer sur un vecteur d'entiers  de taille  n = U1.n_  qui servira de tampon.
	int m = U1.m_, n = U1.n_, p = U2.n_ ;
	if ((U2.m_ != n) || (m_ != m) || (n_ != p)) assert(false);
	int * v ; if (temp_v == NULL)   v = new int[n+1] ;
				 else { v = (*temp_v).v_ ; if (v[0] != n) assert(false); }
	bool Sparse_matmult_set_ij_done = Sparse_matmult_set_ij_check(ij_, U1.ij_, U2.ij_) ;
	if (! Sparse_matmult_set_ij_done) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		(*this).SpUnit_matmult_set_ij_set(U1, U2, v);
		Index_vector temp_v2 ; temp_v2.v_ = v ; v[0] = n ;
		(*this).set_matmult(U1, U2, &temp_v2, false) ;
		temp_v2.v_ = NULL ;
		if (! save_count)  Sparse_matmult_set_ij_uncount(ij_, U1.ij_, U2.ij_) ;
	}
	else {
		int i, j, i1, j2, k ;
		int * cmi_2 = U2.ij_ + p + n ; int * cmv_2 = cmi_2 + U2.nnz_ + U2.nnz_ ; int * rmj_1 ; int * rmv_1 ;
		int * ni_ = ij_ - 1 ; int * ni_2 = U2.ij_ - 1 ; int * nj_1 = U1.ij_ + n - 1 ;
		int * cmi_ = ij_ + m_ + n_ ; double * cmv_ = v_ ;
		int * rmj1_max, * cmi2_max, * rmj_max, * cmi_max ;
		for (j = 1 ; j <= n_ ; j ++) {
			for (k = 1 ; k <= n ; k++)	 v[k] = 0;
			cmi2_max = cmi_2 + ni_2[j] - 1 ;
			while (cmi_2 <= cmi2_max)	{ if (*cmv_2)  v[*cmi_2] = 1 ; else v[*cmi_2] = -1 ;
											cmi_2 ++ ; cmv_2 ++ ; }
			rmj_1 = U1.ij_ + m + n + U1.nnz_ ; rmv_1 = rmj_1 + U1.nnz_ + U1.nnz_ ;
			i = 0 ; cmi_max = cmi_ + ni_[j] - 1 ;
			while (cmi_ <= cmi_max) {
				i1 = i+1 ; i = *cmi_ ; *cmv_ = 0 ;
				for (k = i1 ; k < i ; k ++) { rmj_1 += nj_1[k] ; rmv_1 += nj_1[k] ; }
				rmj1_max = rmj_1 + nj_1[i] - 1 ;
				while (rmj_1 <= rmj1_max) {
					if (v[*rmj_1])  { if (*rmv_1)  (*cmv_) += v[*rmj_1] ; else  (*cmv_) -= v[*rmj_1] ; }
					rmv_1 ++ ; rmj_1 ++ ;
				}
				cmv_ ++ ; cmi_ ++ ;
			}
		}
		rmj_1 = U1.ij_ + m + n + U1.nnz_ ; rmv_1 = rmj_1 + U1.nnz_ + U1.nnz_ ;
		int * nj_ = ij_ - 1 + n_ ; int * rmj_ = ij_ + m_ + n_ + npnz_ ; double * rmv_ = v_ + npnz_ ;
		for (i = 1 ; i <= m_ ; i ++) {
			for (k = 1 ; k <= n ; k++)	 v[k] = 0;
			rmj1_max = rmj_1 + nj_1[i] - 1 ;
			while (rmj_1 <= rmj1_max)	{ if (*rmv_1) v[*rmj_1] = 1 ; else  v[*rmj_1] = -1 ;
											rmj_1 ++ ; rmv_1 ++ ; }
			cmi_2 = U2.ij_ + p + n ; cmv_2 = cmi_2 + U2.nnz_ + U2.nnz_ ;
			j = 0 ; rmj_max = rmj_ + nj_[i] - 1 ;
			while (rmj_ <= rmj_max) {
				j2 = j+1 ; j = *rmj_ ; *rmv_ = 0 ;
				for (k = j2 ; k < j ; k ++) { cmi_2 += ni_2[k] ; cmv_2 += ni_2[k] ; }
				cmi2_max = cmi_2 + ni_2[j] - 1 ;
				while (cmi_2 <= cmi2_max) {
					if (v[*cmi_2])  { if (*cmv_2)  (*rmv_) += v[*cmi_2] ;  else  (*rmv_) -= v[*cmi_2] ; }
					cmv_2 ++ ; cmi_2 ++ ;
				}
				rmv_ ++ ; rmj_ ++ ;
			}
		}
	}
	if (temp_v == NULL) delete [] v ;
}

void Sparse_matrix::add_elemult(const Sparse_matrix &S1, const Sparse_matrix &S2, int ad, bool save_count) {
	if ((n_ != S1.n_) || (m_ != S1.m_) || (n_ != S2.n_) || (m_ != S2.m_)) assert(false);
	if (ad == 0) {	// seule la fonction  set_elemult(S1,S2), soit (*this) = S1*S2, est implementee... pour le moment
		bool Sparse_elemult_set_ij_done = Sparse_elemult_set_ij_check(ij_, S1.ij_, S2.ij_) ;
		if (! Sparse_elemult_set_ij_done) {
			Sparse_XXX_set_ij_uncount(ij_) ;
			(*this).Sparse_elemult_set_ij_set(S1, S2);
			(*this).add_elemult(S1, S2, 0, false) ;
			if (! save_count)  Sparse_elemult_set_ij_uncount(ij_, S1.ij_, S2.ij_) ;
		}
		else {
			int i, j, i1, i2, j1, j2 ;
//			int * nj_ = ij_ + n_ ;
            int * nj_1 = S1.ij_ + n_  ; int * nj_2 = S2.ij_ + n_  ;
			int * rmj_1 = nj_1 + m_ + S1.npnz_ ;	int * rmj_2 = nj_2 + m_ + S2.npnz_ ;	int * rmj1_max, * rmj2_max ;
			double * rmv_ = v_ + npnz_ ; double * rmv_1 = S1.v_ + S1.npnz_ ; double * rmv_2 = S2.v_ + S2.npnz_ ;
			for (i = 0 ; i < m_ ; i ++) {
				rmj1_max = rmj_1 + nj_1[i] - 1 ; rmj2_max = rmj_2 + nj_2[i] - 1 ;
				j1 = (*rmj_1) ; j2 = (*rmj_2) ;
				while ((rmj_1 <= rmj1_max) && (rmj_2 <= rmj2_max)) {
					if (j1 == j2) {
						*rmv_ = (*rmv_1) * (*rmv_2) ; rmv_ ++ ; rmv_1 ++ ; rmv_2 ++ ;
						 rmj_1 ++ ; j1 = (*rmj_1) ; rmj_2 ++ ; j2 = (*rmj_2) ;
					} else
						if (j1 < j2) { rmv_1 ++ ; rmj_1 ++ ; j1 = (*rmj_1) ; }
						else		 { rmv_2 ++ ; rmj_2 ++ ; j2 = (*rmj_2) ; }
				}
				rmv_1 += rmj1_max + 1 - rmj_1 ; rmj_1 = 1 + rmj1_max ;
				rmv_2 += rmj2_max + 1 - rmj_2 ; rmj_2 = 1 + rmj2_max ;
			}
//			int * ni_ = ij_ ;
            int * ni_1 = S1.ij_ ;	int * ni_2 = S2.ij_  ;
			int * cmi_1 = ni_1 + m_ + n_ ; int * cmi_2 = ni_2 + m_ + n_ ; int * cmi1_max, * cmi2_max ;
			double * cmv_ = v_ ; double * cmv_1 = S1.v_ ; double * cmv_2 = S2.v_ ;
			for (j = 0 ; j < n_ ; j ++) {
				cmi1_max = cmi_1 + ni_1[j] - 1 ; cmi2_max = cmi_2 + ni_2[j] - 1 ;
				i1 = (*cmi_1) ; i2 = (*cmi_2) ;
				while ((cmi_1 <= cmi1_max) && (cmi_2 <= cmi2_max)) {
					if (i1 == i2) {
						*cmv_ = (*cmv_1) * (*cmv_2) ; cmv_ ++ ; cmv_1 ++ ; cmv_2 ++ ;
						 cmi_1 ++ ; i1 = (*cmi_1) ; cmi_2 ++ ; i2 = (*cmi_2) ;
					} else
						if (i1 < i2) { cmv_1 ++ ; cmi_1 ++ ; i1 = (*cmi_1) ; }
						else		 { cmv_2 ++ ; cmi_2 ++ ; i2 = (*cmi_2) ; }
				}
				cmv_1 += cmi1_max + 1 - cmi_1 ; cmi_1 = 1 + cmi1_max ;
				cmv_2 += cmi2_max + 1 - cmi_2 ; cmi_2 = 1 + cmi2_max ;
			}
		}
	}
	else assert(false) ;// seule la fonction  set_elemult(S1,S2), soit (*this) = S1*S2, est implementee
}

void Sparse_matrix::add_elemult(const SpUnit_matrix &U1, const Sparse_matrix &S2, int ad, bool save_count) {
	if ((n_ != U1.n_) || (m_ != U1.m_) || (n_ != S2.n_) || (m_ != S2.m_)) assert(false);
	if (ad == 0) {	// seule la fonction  set_elemult(U1,S2), soit (*this) = U1*S2, est implementee
		bool Sparse_elemult_set_ij_done = Sparse_elemult_set_ij_check(ij_, U1.ij_, S2.ij_) ;
		if (! Sparse_elemult_set_ij_done) {
			Sparse_XXX_set_ij_uncount(ij_) ;
			Sparse_matrix T1(m_, n_) ; T1.ntnz_ = U1.nnz_ ; T1.npnz_ = U1.nnz_ ; T1.ij_ = U1.ij_ ;
			(*this).Sparse_elemult_set_ij_set(T1, S2);
			T1.ij_ = NULL ;
			(*this).add_elemult(U1, S2, 0, false) ;
			if (! save_count)  	Sparse_elemult_set_ij_uncount(ij_, U1.ij_, S2.ij_) ;
		}
		else {
			int i, j, i1, i2, j1, j2 ;
            int * nj_1 = U1.ij_ + n_ ; int * nj_2 = S2.ij_ + n_  ;
			int * rmj_1 = nj_1 + m_ + U1.nnz_ ; int * rmj_2 = nj_2 + m_ + S2.npnz_ ; int * rmj1_max, * rmj2_max ;
			double * rmv_ = v_ + npnz_ ; int * rmv_1 = rmj_1 + 2 * U1.nnz_ ;	double * rmv_2 = S2.v_ + S2.npnz_ ;
			for (i = 0 ; i < m_ ; i ++) {
				rmj1_max = rmj_1 + nj_1[i] - 1 ; rmj2_max = rmj_2 + nj_2[i] - 1 ;
				j1 = (*rmj_1) ; j2 = (*rmj_2) ;
				while ((rmj_1 <= rmj1_max) && (rmj_2 <= rmj2_max)) {
					if (j1 == j2) {
						if (*rmv_1)  *rmv_ = (*rmv_2) ;  else  *rmv_ = -(*rmv_2) ;
						rmv_ ++ ; rmv_1 ++ ; rmv_2 ++ ; rmj_1 ++ ; j1 = (*rmj_1) ; rmj_2 ++ ; j2 = (*rmj_2) ;
					} else
						if (j1 < j2) { rmv_1 ++ ; rmj_1 ++ ; j1 = (*rmj_1) ; }
						else		 { rmv_2 ++ ; rmj_2 ++ ; j2 = (*rmj_2) ; }
				}
				rmv_1 += rmj1_max + 1 - rmj_1 ; rmj_1 = 1 + rmj1_max ;
				rmv_2 += rmj2_max + 1 - rmj_2 ; rmj_2 = 1 + rmj2_max ;
			}
            int * ni_1 = U1.ij_ ; int * ni_2 = S2.ij_ ;
			int * cmi_1 = ni_1 + m_ + n_ ; int * cmi_2 = ni_2 + m_ + n_ ; int * cmi1_max, * cmi2_max ;
			double * cmv_ = v_ ; int * cmv_1 = cmi_1 + 2 * U1.nnz_ ; double * cmv_2 = S2.v_ ;
			for (j = 0 ; j < n_ ; j ++) {
				cmi1_max = cmi_1 + ni_1[j] - 1 ; cmi2_max = cmi_2 + ni_2[j] - 1 ;
				i1 = (*cmi_1) ; i2 = (*cmi_2) ;
				while ((cmi_1 <= cmi1_max) && (cmi_2 <= cmi2_max)) {
					if (i1 == i2) {
						if (*cmv_1)  *cmv_ = (*cmv_2) ;  else  *cmv_ = -(*cmv_2) ;
						cmv_ ++ ; cmv_1 ++ ; cmv_2 ++ ; cmi_1 ++ ; i1 = (*cmi_1) ; cmi_2 ++ ; i2 = (*cmi_2) ;
					} else
						if (i1 < i2) { cmv_1 ++ ; cmi_1 ++ ; i1 = (*cmi_1) ; }
						else		 { cmv_2 ++ ; cmi_2 ++ ; i2 = (*cmi_2) ; }
				}
				cmv_1 += cmi1_max + 1 - cmi_1 ; cmi_1 = 1 + cmi1_max ;
				cmv_2 += cmi2_max + 1 - cmi_2 ; cmi_2 = 1 + cmi2_max ;
			}
		}
	}
	else assert(false) ;// seule la fonction  set_elemult(U1,S2), soit (*this) = U1*S2, est implementee
}

void Sparse_matrix::add_elemult(const Sparse_matrix &S, const SpUnit_matrix &U, int ad, bool save_count) {
	(*this).add_elemult(U, S, ad) ;
}

void Sparse_matrix::add_elemult(const Fortran_vector &v, const Sparse_matrix &S, int ad, bool save_count, bool transpose) {
	// v * S : multiplication elementwise par colonne, par un vecteur de taille m = S.nblin()
	double * temp_v = v.v_ ;
	if ((n_ != S.n_) || (m_ != S.m_)) assert(false);
	int m, n ; double * cmv_ ; double * rmv_ ; double * cmv_S ; double * rmv_S ; int * nj_S ; int * ni_S ;
	if (transpose) {cmv_S = S.v_ + S.npnz_ ; rmv_S = S.v_ ; nj_S = S.ij_ ; ni_S = S.ij_ + n_ ; m = n_ ; n = m_ ; }
	else {cmv_S = S.v_ ; rmv_S = S.v_ + S.npnz_ ; nj_S = S.ij_ + n_ ; ni_S = S.ij_ ; m = m_ ; n = n_ ; }
	if ((int)temp_v[0] != m) assert(false);
	if (ad == 0) { // fonction 'set_elemult()'
		if (ntnz_ != S.ntnz_) {
			Sparse_matrix T(m_, n_, S.ntnz_) ; T.ntnz_ = S.ntnz_ ;
			T.add_elemult(v, S, 0, false, transpose) ;
			Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; delete [] v_ ; ij_ = T.ij_ ; v_ = T.v_ ; ntnz_ = T.ntnz_ ; npnz_ = T.npnz_ ;
			T.ij_ = NULL ; T.v_ = NULL ;
		}
		else {
			int * nj_ ; int * cmi_ ;
			if (transpose)	{nj_ = ij_ ; cmi_ = ij_ + m + n + npnz_ ; rmv_ = v_ ; cmv_ = v_ + npnz_ ;}
			else			{nj_ = ij_ + n_ ; cmi_ = ij_ + m + n ; cmv_ = v_ ; rmv_ = v_ + npnz_ ; }
			int i, k, nj, nt = n_ + m_ + ntnz_ ; int * S_ij = S.ij_ ; for (k = 0 ; k < nt ; k ++)  ij_[k] = S_ij[k] ;
			S_ij += S.npnz_ - npnz_ ; nt = n_ + m_ + npnz_ + ntnz_ ; for (k = n_ + m_ + npnz_ ; k < nt ; k ++)   ij_[k] = S_ij[k] ;
			for (i = 0 ; i < m ; i ++) {
				temp_v ++ ; nj = nj_[i] ;
				for (k = 1 ; k <= nj ; k ++) {
					(*rmv_) = (*rmv_S) * (*temp_v) ; rmv_S ++ ; rmv_ ++ ;
				}
			}
			temp_v = v.v_ ;
			for (k = 1 ; k <= ntnz_ ; k ++) {
				(*cmv_) = (*cmv_S) * temp_v[*cmi_] ;
				cmv_S ++ ; cmv_ ++ ; cmi_ ++ ;
			}
		}
	}
	else { // ad = 1  ou  ad = -1
		bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, ij_, S.ij_) ;
		if (!(Sparse_add_set_ij_done)) {
			Sparse_XXX_set_ij_uncount(ij_) ;
			Sparse_matrix T(m_, n_, S.ntnz_) ; T.ntnz_ = S.ntnz_ ; T.add_elemult(v, S, 0, false, transpose) ;
			if (ad == 1)  (*this).add_add(*this, T, 0, false) ; else  (*this).add_sub(*this, T, 0, false) ;
			if (save_count)  Sparse_add_set_ij_count_update(ij_, ij_, S.ij_) ;
		}
		else {
			int i, j, iS, jS ;
			int * rmj_S ; int * cmi_S ; int * rmjS_max ; int * cmiS_max ;
			if (transpose) { rmv_ = v_ ; cmv_ = v_ + npnz_ ; rmj_S = S.ij_ + n + m ; cmi_S = S.ij_ + n + m + S.npnz_ ; }
			else		   { cmv_ = v_ ; rmv_ = v_ + npnz_ ; rmj_S = S.ij_ + n + m + S.npnz_ ; cmi_S = S.ij_ + n + m ; }
			int * nj_  ; int * ni_  ; int * cmi_ ; int * rmj_  ; int * rmj_max ; int * cmi_max ;
			if (transpose)  { nj_ = ij_ ; ni_ = ij_ + n_ ; cmi_ = ij_ + m + n + npnz_ ; rmj_ = ij_ + m + n ; }
			else			{ nj_ = ij_ + n_ ; ni_ = ij_ ; cmi_ = ij_ + m + n ; rmj_ = ij_ + m + n + npnz_ ; }
			if (ad == 1) {
				for (i = 0 ; i < m ; i ++) {
					temp_v ++ ;	rmj_max = rmj_ + nj_[i] - 1 ; rmjS_max = rmj_S + nj_S[i] - 1 ;
					j = (*rmj_) ; jS = (*rmj_S) ;
					while (rmj_S <= rmjS_max) {
						if (j == jS) {
							*rmv_ += (*temp_v) * (*rmv_S) ; rmv_ ++ ; rmv_S ++ ;
							 rmj_ ++ ; j = (*rmj_) ; rmj_S ++ ; jS = (*rmj_S) ;
						} else
							{ rmv_ ++ ; rmj_ ++ ; j = (*rmj_) ; }
					}
					j = rmj_max - rmj_ + 1 ; rmj_ += j ; rmv_ += j ;
				}
				temp_v = v.v_ ;
				for (j = 0 ; j < n ; j ++) {
					cmi_max = cmi_ + ni_[j] - 1 ; cmiS_max = cmi_S + ni_S[j] - 1 ;
					i = (*cmi_) ; iS = (*cmi_S) ;
					while (cmi_S <= cmiS_max) {
						if (i == iS) {
							*cmv_ += temp_v[iS] * (*cmv_S) ; cmv_ ++ ; cmv_S ++ ;
							 cmi_ ++ ; i = (*cmi_) ; cmi_S ++ ; iS = (*cmi_S) ;
						} else
							{ cmv_ ++ ; cmi_ ++ ; i = (*cmi_) ; }
					}
					i = cmi_max - cmi_ + 1 ; cmi_ += i ; cmv_ += i ;
				}
			}
			else { // ad = -1 : fonction 'sub_elemult()'
				for (i = 0 ; i < m ; i ++) {
					temp_v ++ ;	rmj_max = rmj_ + nj_[i] - 1 ; rmjS_max = rmj_S + nj_S[i] - 1 ;
					j = (*rmj_) ; jS = (*rmj_S) ;
					while (rmj_S <= rmjS_max) {
						if (j == jS) {
							*rmv_ -= (*temp_v) * (*rmv_S) ; rmv_ ++ ; rmv_S ++ ;
							 rmj_ ++ ; j = (*rmj_) ; rmj_S ++ ; jS = (*rmj_S) ;
						} else
							{ rmv_ ++ ; rmj_ ++ ; j = (*rmj_) ; }
					}
					j = rmj_max - rmj_ + 1 ; rmj_ += j ; rmv_ += j ;
				}
				temp_v = v.v_ ;
				for (j = 0 ; j < n ; j ++) {
					cmi_max = cmi_ + ni_[j] - 1 ; cmiS_max = cmi_S + ni_S[j] - 1 ;
					i = (*cmi_) ; iS = (*cmi_S) ;
					while (cmi_S <= cmiS_max) {
						if (i == iS) {
							*cmv_ -= temp_v[iS] * (*cmv_S) ; cmv_ ++ ; cmv_S ++ ;
							 cmi_ ++ ; i = (*cmi_) ; cmi_S ++ ; iS = (*cmi_S) ;
						} else
							{ cmv_ ++ ; cmi_ ++ ; i = (*cmi_) ; }
					}
					i = cmi_max - cmi_ + 1 ; cmi_ += i ; cmv_ += i ;
				}
			}
		}
	}
}

void Sparse_matrix::add_elemult(const Fortran_vector &v, const SpUnit_matrix &U, int ad, bool save_count, bool transpose) {
	// v * U : multiplication elementwise par colonne, par un vecteur de taille m = S.nblin()
	double * temp_v = v.v_ ;
	if ((n_ != U.n_) || (m_ != U.m_)) assert(false);
	int m, n ; double * cmv_ ; double * rmv_ ; int * cmv_U ; int * rmv_U ; int * nj_U ; int * ni_U ;
	if (transpose) {nj_U = U.ij_; ni_U = U.ij_+ n_; m = n_; n = m_; rmv_U = ni_U + m_+ 2*U.nnz_; cmv_U = rmv_U + U.nnz_;}
	else {nj_U = U.ij_+ n_; ni_U = U.ij_; m = m_; n = n_; cmv_U = nj_U + m_+ 2*U.nnz_; rmv_U = cmv_U + U.nnz_; }
	if ((int)temp_v[0] != m) assert(false);
	if (ad == 0) { // fonction 'set_elemult()'
		if (ntnz_ != U.nnz_) {
			Sparse_matrix T(m_, n_, U.nnz_) ; T.ntnz_ = U.nnz_ ;
			T.add_elemult(v, U, 0, false, transpose) ;
			Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; delete [] v_ ; ij_ = T.ij_ ; v_ = T.v_ ; ntnz_ = T.ntnz_ ; npnz_ = T.npnz_ ;
			T.ij_ = NULL ; T.v_ = NULL ;
		}
		else {
			int * nj_ ; int * cmi_ ;
			if (transpose) {nj_ = ij_ ; cmi_ = ij_ + m + n + npnz_ ; rmv_ = v_ ; cmv_ = v_ + npnz_ ;}
			else {nj_ = ij_ + n_ ; cmi_ = ij_ + m + n ; cmv_ = v_ ; rmv_ = v_ + npnz_ ; }
			int i, k, nj, nt = n_ + m_ + ntnz_ ; int * U_ij = U.ij_ ; for (k = 0 ; k < nt ; k ++)  ij_[k] = U_ij[k] ;
			U_ij += U.nnz_ - npnz_ ; nt = n_ + m_ + npnz_ + ntnz_ ; for (k = n_ + m_ + npnz_ ; k < nt ; k ++)  ij_[k] = U_ij[k] ;
			for (i = 0 ; i < m ; i ++) {
				temp_v ++ ; nj = nj_[i] ;
				for (k = 1 ; k <= nj ; k ++) {
					if (*rmv_U)  *rmv_ = (*temp_v) ;  else *rmv_ = -(*temp_v) ;
					rmv_U ++ ; rmv_ ++ ;
				}
			}
			temp_v = v.v_ ;
			for (k = 1 ; k <= ntnz_ ; k ++) {
				if (*cmv_U)  *cmv_ = temp_v[*cmi_] ;  else *cmv_ = -temp_v[*cmi_] ;
				cmv_U ++ ; cmv_ ++ ; cmi_ ++ ;
			}
		}
	}
	else { // ad = 1  ou  ad = -1
		bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, ij_, U.ij_) ;
		if (!(Sparse_add_set_ij_done)) {
			Sparse_XXX_set_ij_uncount(ij_) ;
			Sparse_matrix T(m_, n_, U.nnz_) ; T.add_elemult(v, U, 0, false, transpose) ;
			if (ad == 1)  (*this).add_add(*this, T, 0, false) ; else  (*this).add_sub(*this, T, 0, false) ;
			if (save_count)  Sparse_add_set_ij_count_update(ij_, ij_, U.ij_) ;
		}
		else { // done...
			int i, j, iU, jU ;
			int * rmj_U ; int * cmi_U ; int * rmjU_max ; int * cmiU_max ;
			if (transpose) { rmv_ = v_ ; cmv_ = v_ + npnz_ ; rmj_U = U.ij_ + n + m ; cmi_U = U.ij_ + n + m + U.nnz_ ; }
			else		   { cmv_ = v_ ; rmv_ = v_ + npnz_ ; rmj_U = U.ij_ + n + m + U.nnz_ ; cmi_U = U.ij_ + n + m ; }
			int * nj_  ; int * ni_  ; int * cmi_ ; int * rmj_  ; int * rmj_max ; int * cmi_max ;
			if (transpose)  { nj_ = ij_ ; ni_ = ij_ + n_ ; cmi_ = ij_ + m + n + npnz_ ; rmj_ = ij_ + m + n ; }
			else			{ nj_ = ij_ + n_ ; ni_ = ij_ ; cmi_ = ij_ + m + n ; rmj_ = ij_ + m + n + npnz_ ; }
			if (ad == 1) {
				for (i = 0 ; i < m ; i ++) {
					temp_v ++ ;	rmj_max = rmj_ + nj_[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ;
					j = (*rmj_) ; jU = (*rmj_U) ;
					while (rmj_U <= rmjU_max) {
						if (j == jU) {
							if (*rmv_U)  *rmv_ += *temp_v ; else  *rmv_ -= *temp_v ;
							rmv_ ++ ; rmv_U ++ ; rmj_ ++ ; j = (*rmj_) ; rmj_U ++ ; jU = (*rmj_U) ;
						} else
							{ rmv_ ++ ; rmj_ ++ ; j = (*rmj_) ; }
					}
					j = rmj_max - rmj_ + 1 ; rmj_ += j ; rmv_ += j ;
				}
				temp_v = v.v_ ;
				for (j = 0 ; j < n ; j ++) {
					cmi_max = cmi_ + ni_[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ;
					i = (*cmi_) ; iU = (*cmi_U) ;
					while (cmi_U <= cmiU_max) {
						if (i == iU) {
							if (*cmv_U)  *cmv_ += temp_v[iU] ; else  *cmv_ -= temp_v[iU] ;
							cmv_ ++ ; cmv_U ++ ; cmi_ ++ ; i = (*cmi_) ; cmi_U ++ ; iU = (*cmi_U) ;
						} else
							{ cmv_ ++ ; cmi_ ++ ; i = (*cmi_) ; }
					}
					i = cmi_max - cmi_ + 1 ; cmi_ += i ; cmv_ += i ;
				}
			}
			else { // ad = -1 : fonction 'sub_elemult()'
				for (i = 0 ; i < m ; i ++) {
					temp_v ++ ;	rmj_max = rmj_ + nj_[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ;
					j = (*rmj_) ; jU = (*rmj_U) ;
					while (rmj_U <= rmjU_max) {
						if (j == jU) {
							if (*rmv_U)  *rmv_ -= *temp_v ; else  *rmv_ += *temp_v ;
							rmv_ ++ ; rmv_U ++ ; rmj_ ++ ; j = (*rmj_) ; rmj_U ++ ; jU = (*rmj_U) ;
						} else
							{ rmv_ ++ ; rmj_ ++ ; j = (*rmj_) ; }
					}
					j = rmj_max - rmj_ + 1 ; rmj_ += j ; rmv_ += j ;
				}
				temp_v = v.v_ ;
				for (j = 0 ; j < n ; j ++) {
					cmi_max = cmi_ + ni_[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ;
					i = (*cmi_) ; iU = (*cmi_U) ;
					while (cmi_U <= cmiU_max) {
						if (i == iU) {
							if (*cmv_U)  *cmv_ -= temp_v[iU] ; else  *cmv_ += temp_v[iU] ;
							cmv_ ++ ; cmv_U ++ ; cmi_ ++ ; i = (*cmi_) ; cmi_U ++ ; iU = (*cmi_U) ;
						} else
							{ cmv_ ++ ; cmi_ ++ ; i = (*cmi_) ; }
					}
					i = cmi_max - cmi_ + 1 ; cmi_ += i ; cmv_ += i ;
				}
			}
		}
	}
}

void Sparse_matrix::add_mult(const double &a, const Sparse_matrix &S, int ad, bool save_count) {
	// a * S : multiplication par un scalaire,
	if ((n_ != S.n_) || (m_ != S.m_)) assert(false);
	if (ad == 0) { // fonction 'set_mult()'
		if (ntnz_ != S.ntnz_) {
			Sparse_matrix T(m_, n_, S.ntnz_) ; T.ntnz_ = S.ntnz_ ;
			T.set_mult(a, S) ;
			Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; delete [] v_ ; ij_ = T.ij_ ; v_ = T.v_ ; ntnz_ = T.ntnz_ ; npnz_ = T.npnz_ ;
			T.ij_ = NULL ; T.v_ = NULL ;
		}
		else {
			int k, nt = n_ + m_ + ntnz_ ; int * temp_ij = S.ij_ ; for (k = 0 ; k < nt ; k ++)  ij_[k] = temp_ij[k] ;
			temp_ij += S.npnz_ - npnz_ ; nt = n_ + m_ + npnz_ + ntnz_ ; for (k = n_ + m_ + npnz_ ; k < nt ; k ++)   ij_[k] = temp_ij[k] ;
			double* temp_v = v_ ; double* temp_S = S.v_ ; for (k = 0 ; k < ntnz_ ; k ++)  {temp_v[k] = a * temp_S[k] ; }
			temp_v += npnz_ ; temp_S += S.npnz_ ; for (k = 0 ; k < ntnz_ ; k ++)  {temp_v[k] = a * temp_S[k] ; }
		}
	}
	else { // ad = 1  ou  ad = -1
		double aa ; if (ad == 1) aa = a  ; else aa = -a ;
		bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, ij_, S.ij_) ;
		if (!(Sparse_add_set_ij_done)) {
			Sparse_XXX_set_ij_uncount(ij_) ;
			Sparse_matrix T(m_, n_, S.ntnz_) ; T.ntnz_ = S.ntnz_ ; T.set_mult(aa, S) ;
			(*this).add_add(*this, T, 0, false) ;
			if (save_count)  Sparse_add_set_ij_count_update(ij_, ij_, S.ij_) ;
		}
		else {
			int i, j, iS, jS ;
			int * nj_S = S.ij_ + n_  ; int * rmj_S = nj_S + m_ + S.npnz_ ; int * rmjS_max ;
			double * cmv_S = S.v_ ;	double * rmv_S = S.v_ + S.npnz_ ;
			int * ni_S = S.ij_ ; int * cmi_S = ni_S + m_ + n_ ;	int * cmiS_max ;
			double * cmv_ = v_ ; double * rmv_ = v_ + npnz_ ;
			int * nj_ = ij_ + n_ ; int * rmj_ = nj_ + m_ + npnz_ ; int * rmj_max ;
			for (i = 0 ; i < m_ ; i ++) {
				rmj_max = rmj_ + nj_[i] - 1 ; rmjS_max = rmj_S + nj_S[i] - 1 ;
				j = (*rmj_) ; jS = (*rmj_S) ;
				while (rmj_S <= rmjS_max) {
					if (j == jS) {
						*rmv_ += aa * (*rmv_S) ; rmv_ ++ ; rmv_S ++ ;
						 rmj_ ++ ; j = (*rmj_) ; rmj_S ++ ; jS = (*rmj_S) ;
					} else
						{ rmv_ ++ ; rmj_ ++ ; j = (*rmj_) ; }
				}
				j = rmj_max - rmj_ + 1 ; rmj_ += j ; rmv_ += j ;
			}
			int * ni_ = ij_ ; int * cmi_ = ni_ + m_ + n_ ; int * cmi_max ;
			for (j = 0 ; j < n_ ; j ++) {
				cmi_max = cmi_ + ni_[j] - 1 ; cmiS_max = cmi_S + ni_S[j] - 1 ;
				i = (*cmi_) ; iS = (*cmi_S) ;
				while (cmi_S <= cmiS_max) {
					if (i == iS) {
						*cmv_ += aa * (*cmv_S) ; cmv_ ++ ; cmv_S ++ ;
						 cmi_ ++ ; i = (*cmi_) ; cmi_S ++ ; iS = (*cmi_S) ;
					} else
						{ cmv_ ++ ; cmi_ ++ ; i = (*cmi_) ; }
				}
				i = cmi_max - cmi_ + 1 ; cmi_ += i ; cmv_ += i ;
			}
		}
	}
}

void Sparse_matrix::add_mult(const double &a, const SpUnit_matrix &U, int ad, bool save_count) {
	// a * U : multiplication par un scalaire
	if ((n_ != U.n_) || (m_ != U.m_)) assert(false);
	if (ad == 0) {
		if (ntnz_ != U.nnz_) {
			Sparse_matrix T(m_, n_, U.nnz_) ; T.ntnz_ = U.nnz_ ;
			T.set_mult(a, U) ;
			Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; delete [] v_ ; ij_ = T.ij_ ; v_ = T.v_ ; ntnz_ = T.ntnz_ ; npnz_ = T.npnz_ ;
			T.ij_ = NULL ; T.v_ = NULL ;
		}
		else {
			int k, nt = n_ + m_ + ntnz_ ; int * U_ij = U.ij_ ; for (k = 0 ; k < nt ; k ++)  ij_[k] = U_ij[k] ;
			U_ij += U.nnz_ - npnz_ ; nt = n_ + m_ + npnz_ + ntnz_ ; for (k = n_ + m_ + npnz_ ; k < nt ; k ++)   ij_[k] = U_ij[k] ;
			double* temp_v = v_ ; U_ij += nt ; for (k = 0 ; k < ntnz_ ; k ++) {	if (U_ij[k]) v_[k] = a ;	else  v_[k] = -a ;}
			temp_v += npnz_ ; U_ij += ntnz_ ; for (k = 0 ; k < ntnz_ ; k ++) {	if (U_ij[k]) v_[k] = a ;	else  v_[k] = -a ;}
		}
	}
	else { // ad = 1  ou  ad = -1
		double aa ; if (ad == 1) aa = a  ; else aa = -a ;
		bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, ij_, U.ij_) ;
		if (!(Sparse_add_set_ij_done)) {
			Sparse_XXX_set_ij_uncount(ij_) ;
			Sparse_matrix T(m_, n_, U.nnz_) ; T.ntnz_ = U.nnz_ ; T.set_mult(aa, U) ;
			(*this).add_add(*this, T, 0, false) ;
			if (save_count)  Sparse_add_set_ij_count_update(ij_, ij_, U.ij_) ;
		}
		else {
			int i, j, iU, jU ;
			int * nj_U = U.ij_ + n_  ; int * rmj_U = nj_U + m_ + U.nnz_ ; int * rmjU_max ;
			int * cmv_U = rmj_U + U.nnz_ ;	int * rmv_U = cmv_U + U.nnz_ ;
			int * ni_U = U.ij_ ; int * cmi_U = ni_U + m_ + n_ ;	int * cmiU_max ;
			double * cmv_ = v_ ; double * rmv_ = v_ + npnz_ ;
			int * nj_ = ij_ + n_ ; int * rmj_ = nj_ + m_ + npnz_ ; int * rmj_max ;
			for (i = 0 ; i < m_ ; i ++) {
				rmj_max = rmj_ + nj_[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ;
				j = (*rmj_) ; jU = (*rmj_U) ;
				while (rmj_U <= rmjU_max) {
					if (j == jU) {
						if (*rmv_U)  *rmv_ += aa ; else *rmv_ -= aa ;
						rmv_ ++ ; rmv_U ++ ; rmj_ ++ ; j = (*rmj_) ; rmj_U ++ ; jU = (*rmj_U) ;
					} else
						{ rmv_ ++ ; rmj_ ++ ; j = (*rmj_) ; }
				}
				j = rmj_max - rmj_ + 1 ; rmj_ += j ; rmv_ += j ;
			}
			int * ni_ = ij_ ; int * cmi_ = ni_ + m_ + n_ ; int * cmi_max ;
			for (j = 0 ; j < n_ ; j ++) {
				cmi_max = cmi_ + ni_[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ;
				i = (*cmi_) ; iU = (*cmi_U) ;
				while (cmi_U <= cmiU_max) {
					if (i == iU) {
						if (*cmv_U)  *cmv_ += aa ; else  *cmv_ -= aa ;
						cmv_ ++ ; cmv_U ++ ; cmi_ ++ ; i = (*cmi_) ; cmi_U ++ ; iU = (*cmi_U) ;
					} else
						{ cmv_ ++ ; cmi_ ++ ; i = (*cmi_) ; }
				}
				i = cmi_max - cmi_ + 1 ; cmi_ += i ; cmv_ += i ;
			}
		}
	}
}

void Sparse_matrix::add_elediv(const Sparse_matrix &S, const Fortran_vector &v, int ad, bool save_count, bool transpose) {
	// S / v : division elementwise par colonne, par un vecteur de taille m = S.nblin()
	double * temp_v = v.v_ ;
	if ((n_ != S.n_) || (m_ != S.m_)) assert(false);
	int m, n ; double * cmv_ ; double * rmv_ ; double * cmv_S ; double * rmv_S ; int * nj_S ; int * ni_S ;
	if (transpose) {cmv_S = S.v_ + S.npnz_ ; rmv_S = S.v_ ; nj_S = S.ij_ ; ni_S = S.ij_ + n_ ; m = n_ ; n = m_ ; }
	else {cmv_S = S.v_ ; rmv_S = S.v_ + S.npnz_ ; nj_S = S.ij_ + n_ ; ni_S = S.ij_ ; m = m_ ; n = n_ ; }
	if ((int)temp_v[0] != m) assert(false);
	if (ad == 0) { // fonction 'set_elemult()'
		if (ntnz_ != S.ntnz_) {
			Sparse_matrix T(m_, n_, S.ntnz_) ; T.ntnz_ = S.ntnz_ ;
			T.add_elediv(S, v, 0, false, transpose) ;
			Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; delete [] v_ ; ij_ = T.ij_ ; v_ = T.v_ ; ntnz_ = T.ntnz_ ; npnz_ = T.npnz_ ;
			T.ij_ = NULL ; T.v_ = NULL ;
		}
		else {
			int * nj_ ; int * cmi_ ;
			if (transpose)	{nj_ = ij_ ; cmi_ = ij_ + m + n + npnz_ ; rmv_ = v_ ; cmv_ = v_ + npnz_ ;}
			else			{nj_ = ij_ + n_ ; cmi_ = ij_ + m + n ; cmv_ = v_ ; rmv_ = v_ + npnz_ ; }
			int i, k, nj, nt = n_ + m_ + ntnz_ ; int * S_ij = S.ij_ ; for (k = 0 ; k < nt ; k ++)  ij_[k] = S_ij[k] ;
			S_ij += S.npnz_ - npnz_ ; nt = n_ + m_ + npnz_ + ntnz_ ; for (k = n_ + m_ + npnz_ ; k < nt ; k ++)   ij_[k] = S_ij[k] ;
			for (i = 0 ; i < m ; i ++) {
				temp_v ++ ; nj = nj_[i] ;
				for (k = 1 ; k <= nj ; k ++) {
					(*rmv_) = (*rmv_S) / (*temp_v) ; rmv_S ++ ; rmv_ ++ ;
				}
			}
			temp_v = v.v_ ;
			for (k = 1 ; k <= ntnz_ ; k ++) {
				(*cmv_) = (*cmv_S) / temp_v[*cmi_] ;
				cmv_S ++ ; cmv_ ++ ; cmi_ ++ ;
			}
		}
	}
	else { // ad = 1  ou  ad = -1
		bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, ij_, S.ij_) ;
		if (!(Sparse_add_set_ij_done)) {
			Sparse_XXX_set_ij_uncount(ij_) ;
			Sparse_matrix T(m_, n_, S.ntnz_) ; T.ntnz_ = S.ntnz_ ; T.add_elediv(S, v, 0, false, transpose) ;
			if (ad == 1)  (*this).add_add(*this, T, 0, false) ; else  (*this).add_sub(*this, T, 0, false) ;
			if (save_count)  Sparse_add_set_ij_count_update(ij_, ij_, S.ij_) ;
		}
		else {
			int i, j, iS, jS ;
			int * rmj_S ; int * cmi_S ; int * rmjS_max ; int * cmiS_max ;
			if (transpose) { rmv_ = v_ ; cmv_ = v_ + npnz_ ; rmj_S = S.ij_ + n + m ; cmi_S = S.ij_ + n + m + S.npnz_ ; }
			else		   { cmv_ = v_ ; rmv_ = v_ + npnz_ ; rmj_S = S.ij_ + n + m + S.npnz_ ; cmi_S = S.ij_ + n + m ; }
			int * nj_  ; int * ni_  ; int * cmi_ ; int * rmj_  ; int * rmj_max ; int * cmi_max ;
			if (transpose)  { nj_ = ij_ ; ni_ = ij_ + n_ ; cmi_ = ij_ + m + n + npnz_ ; rmj_ = ij_ + m + n ; }
			else			{ nj_ = ij_ + n_ ; ni_ = ij_ ; cmi_ = ij_ + m + n ; rmj_ = ij_ + m + n + npnz_ ; }
			if (ad == 1) {
				for (i = 0 ; i < m ; i ++) {
					temp_v ++ ;	rmj_max = rmj_ + nj_[i] - 1 ; rmjS_max = rmj_S + nj_S[i] - 1 ;
					j = (*rmj_) ; jS = (*rmj_S) ;
					while (rmj_S <= rmjS_max) {
						if (j == jS) {
							*rmv_ += (*rmv_S) / (*temp_v) ; rmv_ ++ ; rmv_S ++ ;
							 rmj_ ++ ; j = (*rmj_) ; rmj_S ++ ; jS = (*rmj_S) ;
						} else
							{ rmv_ ++ ; rmj_ ++ ; j = (*rmj_) ; }
					}
					j = rmj_max - rmj_ + 1 ; rmj_ += j ; rmv_ += j ;
				}
				temp_v = v.v_ ;
				for (j = 0 ; j < n ; j ++) {
					cmi_max = cmi_ + ni_[j] - 1 ; cmiS_max = cmi_S + ni_S[j] - 1 ;
					i = (*cmi_) ; iS = (*cmi_S) ;
					while (cmi_S <= cmiS_max) {
						if (i == iS) {
							*cmv_ += (*cmv_S) / temp_v[iS] ; cmv_ ++ ; cmv_S ++ ;
							 cmi_ ++ ; i = (*cmi_) ; cmi_S ++ ; iS = (*cmi_S) ;
						} else
							{ cmv_ ++ ; cmi_ ++ ; i = (*cmi_) ; }
					}
					i = cmi_max - cmi_ + 1 ; cmi_ += i ; cmv_ += i ;
				}
			}
			else { // ad = -1 : fonction 'sub_elemult()'
				for (i = 0 ; i < m ; i ++) {
					temp_v ++ ;	rmj_max = rmj_ + nj_[i] - 1 ; rmjS_max = rmj_S + nj_S[i] - 1 ;
					j = (*rmj_) ; jS = (*rmj_S) ;
					while (rmj_S <= rmjS_max) {
						if (j == jS) {
							*rmv_ -= (*rmv_S) / (*temp_v) ; rmv_ ++ ; rmv_S ++ ;
							 rmj_ ++ ; j = (*rmj_) ; rmj_S ++ ; jS = (*rmj_S) ;
						} else
							{ rmv_ ++ ; rmj_ ++ ; j = (*rmj_) ; }
					}
					j = rmj_max - rmj_ + 1 ; rmj_ += j ; rmv_ += j ;
				}
				temp_v = v.v_ ;
				for (j = 0 ; j < n ; j ++) {
					cmi_max = cmi_ + ni_[j] - 1 ; cmiS_max = cmi_S + ni_S[j] - 1 ;
					i = (*cmi_) ; iS = (*cmi_S) ;
					while (cmi_S <= cmiS_max) {
						if (i == iS) {
							*cmv_ -= (*cmv_S) / temp_v[iS] ; cmv_ ++ ; cmv_S ++ ;
							 cmi_ ++ ; i = (*cmi_) ; cmi_S ++ ; iS = (*cmi_S) ;
						} else
							{ cmv_ ++ ; cmi_ ++ ; i = (*cmi_) ; }
					}
					i = cmi_max - cmi_ + 1 ; cmi_ += i ; cmv_ += i ;
				}
			}
		}
	}
}

void Sparse_matrix::add_elediv(const SpUnit_matrix &U, const Fortran_vector &v, int ad, bool save_count, bool transpose) {
	// U / v : division elementwise par colonne, par un vecteur de taille m = S.nblin()
	double * temp_v = v.v_ ;
	if ((n_ != U.n_) || (m_ != U.m_)) assert(false);
	int m, n ; double * cmv_ ; double * rmv_ ; int * cmv_U ; int * rmv_U ; int * nj_U ; int * ni_U ;
	if (transpose) {nj_U = U.ij_; ni_U = U.ij_+ n_; m = n_; n = m_; rmv_U = ni_U + m_+ 2*U.nnz_; cmv_U = rmv_U + U.nnz_;}
	else {nj_U = U.ij_+ n_; ni_U = U.ij_; m = m_; n = n_; cmv_U = nj_U + m_+ 2*U.nnz_; rmv_U = cmv_U + U.nnz_; }
	if ((int)temp_v[0] != m) assert(false);
	if (ad == 0) { // fonction 'set_elemult()'
		if (ntnz_ != U.nnz_) {
			Sparse_matrix T(m_, n_, U.nnz_) ; T.ntnz_ = U.nnz_ ;
			T.add_elediv(U, v, 0, false, transpose) ;
			Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; delete [] v_ ; ij_ = T.ij_ ; v_ = T.v_ ; ntnz_ = T.ntnz_ ; npnz_ = T.npnz_ ;
			T.ij_ = NULL ; T.v_ = NULL ;
		}
		else {
			int * nj_ ; int * cmi_ ;
			if (transpose) {nj_ = ij_ ; cmi_ = ij_ + m + n + npnz_ ; rmv_ = v_ ; cmv_ = v_ + npnz_ ;}
			else {nj_ = ij_ + n_ ; cmi_ = ij_ + m + n ; cmv_ = v_ ; rmv_ = v_ + npnz_ ; }
			int i, k, nj, nt = n_ + m_ + ntnz_ ; int * U_ij = U.ij_ ; for (k = 0 ; k < nt ; k ++)  ij_[k] = U_ij[k] ;
			U_ij += U.nnz_ - npnz_ ; nt = n_ + m_ + npnz_ + ntnz_ ; for (k = n_ + m_ + npnz_ ; k < nt ; k ++)   ij_[k] = U_ij[k] ;
			for (i = 0 ; i < m ; i ++) {
				temp_v ++ ; nj = nj_[i] ;
				for (k = 1 ; k <= nj ; k ++) {
					if (*rmv_U)  *rmv_ = 1./ (*temp_v) ;  else *rmv_ = -1./ (*temp_v) ;
					rmv_U ++ ; rmv_ ++ ;
				}
			}
			temp_v = v.v_ ;
			for (k = 1 ; k <= ntnz_ ; k ++) {
				if (*cmv_U)  *cmv_ = 1./ temp_v[*cmi_] ;  else *cmv_ = -1./ temp_v[*cmi_] ;
				cmv_U ++ ; cmv_ ++ ; cmi_ ++ ;
			}
		}
	}
	else { // ad = 1  ou  ad = -1
		bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, ij_, U.ij_) ;
		if (!(Sparse_add_set_ij_done)) {
			Sparse_XXX_set_ij_uncount(ij_) ;
			Sparse_matrix T(m_, n_, U.nnz_) ; T.ntnz_ = U.nnz_ ; T.add_elediv(U, v, 0, false, transpose) ;
			if (ad == 1)  (*this).add_add(*this, T, 0, false) ; else  (*this).add_sub(*this, T, 0, false) ;
			if (save_count)  Sparse_add_set_ij_count_update(ij_, ij_, U.ij_) ;
		}
		else { // done...
			int i, j, iU, jU ;
			int * rmj_U ; int * cmi_U ; int * rmjU_max ; int * cmiU_max ;
			if (transpose) { rmv_ = v_ ; cmv_ = v_ + npnz_ ; rmj_U = U.ij_ + n + m ; cmi_U = U.ij_ + n + m + U.nnz_ ; }
			else		   { cmv_ = v_ ; rmv_ = v_ + npnz_ ; rmj_U = U.ij_ + n + m + U.nnz_ ; cmi_U = U.ij_ + n + m ; }
			int * nj_  ; int * ni_  ; int * cmi_ ; int * rmj_  ; int * rmj_max ; int * cmi_max ;
			if (transpose)  { nj_ = ij_ ; ni_ = ij_ + n_ ; cmi_ = ij_ + m + n + npnz_ ; rmj_ = ij_ + m + n ; }
			else			{ nj_ = ij_ + n_ ; ni_ = ij_ ; cmi_ = ij_ + m + n ; rmj_ = ij_ + m + n + npnz_ ; }
			if (ad == 1) {
				for (i = 0 ; i < m ; i ++) {
					temp_v ++ ;	rmj_max = rmj_ + nj_[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ;
					j = (*rmj_) ; jU = (*rmj_U) ;
					while (rmj_U <= rmjU_max) {
						if (j == jU) {
							if (*rmv_U)  *rmv_ += 1./(*temp_v) ; else  *rmv_ -= 1./(*temp_v) ;
							rmv_ ++ ; rmv_U ++ ; rmj_ ++ ; j = (*rmj_) ; rmj_U ++ ; jU = (*rmj_U) ;
						} else
							{ rmv_ ++ ; rmj_ ++ ; j = (*rmj_) ; }
					}
					j = rmj_max - rmj_ + 1 ; rmj_ += j ; rmv_ += j ;
				}
				temp_v = v.v_ ;
				for (j = 0 ; j < n ; j ++) {
					cmi_max = cmi_ + ni_[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ;
					i = (*cmi_) ; iU = (*cmi_U) ;
					while (cmi_U <= cmiU_max) {
						if (i == iU) {
							if (*cmv_U)  *cmv_ += 1./ temp_v[iU] ; else  *cmv_ -= 1./ temp_v[iU] ;
							cmv_ ++ ; cmv_U ++ ; cmi_ ++ ; i = (*cmi_) ; cmi_U ++ ; iU = (*cmi_U) ;
						} else
							{ cmv_ ++ ; cmi_ ++ ; i = (*cmi_) ; }
					}
					i = cmi_max - cmi_ + 1 ; cmi_ += i ; cmv_ += i ;
				}
			}
			else { // ad = -1 : fonction 'sub_elemult()'
				for (i = 0 ; i < m ; i ++) {
					temp_v ++ ;	rmj_max = rmj_ + nj_[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ;
					j = (*rmj_) ; jU = (*rmj_U) ;
					while (rmj_U <= rmjU_max) {
						if (j == jU) {
							if (*rmv_U)  *rmv_ -= 1./(*temp_v) ; else  *rmv_ += 1./(*temp_v) ;
							rmv_ ++ ; rmv_U ++ ; rmj_ ++ ; j = (*rmj_) ; rmj_U ++ ; jU = (*rmj_U) ;
						} else
							{ rmv_ ++ ; rmj_ ++ ; j = (*rmj_) ; }
					}
					j = rmj_max - rmj_ + 1 ; rmj_ += j ; rmv_ += j ;
				}
				temp_v = v.v_ ;
				for (j = 0 ; j < n ; j ++) {
					cmi_max = cmi_ + ni_[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ;
					i = (*cmi_) ; iU = (*cmi_U) ;
					while (cmi_U <= cmiU_max) {
						if (i == iU) {
							if (*cmv_U)  *cmv_ -= 1./ temp_v[iU] ; else  *cmv_ += 1./ temp_v[iU] ;
							cmv_ ++ ; cmv_U ++ ; cmi_ ++ ; i = (*cmi_) ; cmi_U ++ ; iU = (*cmi_U) ;
						} else
							{ cmv_ ++ ; cmi_ ++ ; i = (*cmi_) ; }
					}
					i = cmi_max - cmi_ + 1 ; cmi_ += i ; cmv_ += i ;
				}
			}
		}
	}
}

void Sparse_matrix::add_div(const Sparse_matrix &S, const double &a, int ad, bool save_count) {
	//   S / a : division par un scalaire
	if ((n_ != S.n_) || (m_ != S.m_)) assert(false);
	if (ad == 0) { // fonction 'set_elemult()'
		if (ntnz_ != S.ntnz_) {
			Sparse_matrix T(m_, n_, S.ntnz_) ; T.ntnz_ = S.ntnz_ ;
			T.set_div(S, a) ;
			Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; delete [] v_ ; ij_ = T.ij_ ; v_ = T.v_ ; ntnz_ = T.ntnz_ ; npnz_ = T.npnz_ ;
			T.ij_ = NULL ; T.v_ = NULL ;
		}
		else {
			int k, nt = n_ + m_ + ntnz_ ; int * temp_ij = S.ij_ ; for (k = 0 ; k < nt ; k ++)  ij_[k] = temp_ij[k] ;
			temp_ij += S.npnz_ - npnz_ ; nt = n_ + m_ + npnz_ + ntnz_ ; for (k = n_ + m_ + npnz_ ; k < nt ; k ++)   ij_[k] = temp_ij[k] ;
			double* temp_v = v_ ; double* temp_S = S.v_ ; for (k = 0 ; k < ntnz_ ; k ++)  {temp_v[k] = temp_S[k] / a ; }
			temp_v += npnz_ ; temp_S += S.npnz_ ; for (k = 0 ; k < ntnz_ ; k ++)  {temp_v[k] = temp_S[k] / a ; }
		}
	}
	else { // ad = 1  ou  ad = -1
		double aa ; if (ad == 1) aa = a  ; else aa = -a ;
		bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, ij_, S.ij_) ;
		if (!(Sparse_add_set_ij_done)) {
			Sparse_XXX_set_ij_uncount(ij_) ;
			Sparse_matrix T(m_, n_, S.ntnz_) ; T.ntnz_ = S.ntnz_ ; T.set_div(S, aa) ;
			(*this).add_add(*this, T, 0, false) ;
			if (save_count)  Sparse_add_set_ij_count_update(ij_, ij_, S.ij_) ;
		}
		else {
			int i, j, iS, jS ;
			int * nj_S = S.ij_ + n_  ; int * rmj_S = nj_S + m_ + S.npnz_ ; int * rmjS_max ;
			double * cmv_S = S.v_ ;	double * rmv_S = S.v_ + S.npnz_ ;
			int * ni_S = S.ij_ ; int * cmi_S = ni_S + m_ + n_ ;	int * cmiS_max ;
			double * cmv_ = v_ ; double * rmv_ = v_ + npnz_ ;
			int * nj_ = ij_ + n_ ; int * rmj_ = nj_ + m_ + npnz_ ; int * rmj_max ;
			for (i = 0 ; i < m_ ; i ++) {
				rmj_max = rmj_ + nj_[i] - 1 ; rmjS_max = rmj_S + nj_S[i] - 1 ;
				j = (*rmj_) ; jS = (*rmj_S) ;
				while (rmj_S <= rmjS_max) {
					if (j == jS) {
						*rmv_ += (*rmv_S) / aa ; rmv_ ++ ; rmv_S ++ ;
						 rmj_ ++ ; j = (*rmj_) ; rmj_S ++ ; jS = (*rmj_S) ;
					} else
						{ rmv_ ++ ; rmj_ ++ ; j = (*rmj_) ; }
				}
				j = rmj_max - rmj_ + 1 ; rmj_ += j ; rmv_ += j ;
			}
			int * ni_ = ij_ ; int * cmi_ = ni_ + m_ + n_ ; int * cmi_max ;
			for (j = 0 ; j < n_ ; j ++) {
				cmi_max = cmi_ + ni_[j] - 1 ; cmiS_max = cmi_S + ni_S[j] - 1 ;
				i = (*cmi_) ; iS = (*cmi_S) ;
				while (cmi_S <= cmiS_max) {
					if (i == iS) {
						*cmv_ += (*cmv_S) / aa ; cmv_ ++ ; cmv_S ++ ;
						 cmi_ ++ ; i = (*cmi_) ; cmi_S ++ ; iS = (*cmi_S) ;
					} else
						{ cmv_ ++ ; cmi_ ++ ; i = (*cmi_) ; }
				}
				i = cmi_max - cmi_ + 1 ; cmi_ += i ; cmv_ += i ;
			}
		}
	}
}

void Sparse_matrix::add_div(const SpUnit_matrix &U, const double &a, int ad, bool save_count) {
	//   U / a : division par un scalaire
	if ((n_ != U.n_) || (m_ != U.m_)) assert(false);
	if (ad == 0) {

		if (ntnz_ != U.nnz_) {
			Sparse_matrix T(m_, n_, U.nnz_) ; T.ntnz_ = U.nnz_ ;
			T.set_div(U, a) ;
			Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; delete [] v_ ; ij_ = T.ij_ ; v_ = T.v_ ; ntnz_ = T.ntnz_ ; npnz_ = T.npnz_ ;
			T.ij_ = NULL ; T.v_ = NULL ;
		}

		else {
			int k, nt = n_ + m_ + ntnz_ ; int * U_ij = U.ij_ ; for (k = 0 ; k < nt ; k ++)  ij_[k] = U_ij[k] ;
			U_ij += U.nnz_ - npnz_ ; nt = n_ + m_ + npnz_ + ntnz_ ; for (k = n_ + m_ + npnz_ ; k < nt ; k ++)   ij_[k] = U_ij[k] ;
			double* temp_v = v_ ; U_ij += nt ; for (k = 0 ; k < ntnz_ ; k ++) {	if (U_ij[k]) v_[k] = 1./ a ;	else  v_[k] = -1./ a ;}
			temp_v += npnz_ ; U_ij += ntnz_ ; for (k = 0 ; k < ntnz_ ; k ++) {	if (U_ij[k]) v_[k] = 1./ a ;	else  v_[k] = -1./ a ;}
		}
	}
	else { // ad = 1  ou  ad = -1
		double aa ; if (ad == 1) aa = a  ; else aa = -a ;
		bool Sparse_add_set_ij_done = Sparse_add_set_ij_check(ij_, ij_, U.ij_) ;
		if (!(Sparse_add_set_ij_done)) {
			Sparse_XXX_set_ij_uncount(ij_) ;
			Sparse_matrix T(m_, n_, U.nnz_) ; T.ntnz_ = U.nnz_ ; T.set_div(U, aa) ;
			(*this).add_add(*this, T, 0, false) ;
			if (save_count)  Sparse_add_set_ij_count_update(ij_, ij_, U.ij_) ;
		}
		else {
			int i, j, iU, jU ;
			int * nj_U = U.ij_ + n_  ; int * rmj_U = nj_U + m_ + U.nnz_ ; int * rmjU_max ;
			int * cmv_U = rmj_U + U.nnz_ ;	int * rmv_U = cmv_U + U.nnz_ ;
			int * ni_U = U.ij_ ; int * cmi_U = ni_U + m_ + n_ ;	int * cmiU_max ;
			double * cmv_ = v_ ; double * rmv_ = v_ + npnz_ ;
			int * nj_ = ij_ + n_ ; int * rmj_ = nj_ + m_ + npnz_ ; int * rmj_max ;
			for (i = 0 ; i < m_ ; i ++) {
				rmj_max = rmj_ + nj_[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ;
				j = (*rmj_) ; jU = (*rmj_U) ;
				while (rmj_U <= rmjU_max) {
					if (j == jU) {
						if (*rmv_U)  *rmv_ += 1./ aa ; else *rmv_ -= 1./ aa ;
						rmv_ ++ ; rmv_U ++ ; rmj_ ++ ; j = (*rmj_) ; rmj_U ++ ; jU = (*rmj_U) ;
					} else
						{ rmv_ ++ ; rmj_ ++ ; j = (*rmj_) ; }
				}
				j = rmj_max - rmj_ + 1 ; rmj_ += j ; rmv_ += j ;
			}
			int * ni_ = ij_ ; int * cmi_ = ni_ + m_ + n_ ; int * cmi_max ;
			for (j = 0 ; j < n_ ; j ++) {
				cmi_max = cmi_ + ni_[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ;
				i = (*cmi_) ; iU = (*cmi_U) ;
				while (cmi_U <= cmiU_max) {
					if (i == iU) {
						if (*cmv_U)  *cmv_ += 1./ aa ; else  *cmv_ -= 1./ aa ;
						cmv_ ++ ; cmv_U ++ ; cmi_ ++ ; i = (*cmi_) ; cmi_U ++ ; iU = (*cmi_U) ;
					} else
						{ cmv_ ++ ; cmi_ ++ ; i = (*cmi_) ; }
				}
				i = cmi_max - cmi_ + 1 ; cmi_ += i ; cmv_ += i ;
			}
		}
	}
}

void Sparse_matrix::add_row_elemult(const Sparse_matrix &S, const Fortran_vector &v, int ad, bool save_count) {
	// S * v : multiplication elementwise PAR LIGNE, par un vecteur de taille n = S.nbcol()
	(*this).add_elemult(v, S, ad, save_count, true) ;
}

void Sparse_matrix::add_row_elemult(const SpUnit_matrix &U, const Fortran_vector &v, int ad, bool save_count) {
	// U * v : multiplication elementwise PAR LIGNE, par un vecteur de taille n = S.nbcol()
	(*this).add_elemult(v, U, ad, save_count, true) ;
}

void Sparse_matrix::add_row_elediv(const Sparse_matrix &S, const Fortran_vector &v, int ad, bool save_count) {
	// S / v : division elementwise PAR LIGNE, par un vecteur de taille n = S.nbcol()
	(*this).add_elediv(S, v, ad, save_count, true) ;
}

void Sparse_matrix::add_row_elediv(const SpUnit_matrix &U, const Fortran_vector &v, int ad, bool save_count) {
	// U / v : division elementwise PAR LIGNE, par un vecteur de taille n = S.nbcol()
	(*this).add_elediv(U, v, ad, save_count, true) ;
}

/***** Classe SpUnit_matrix *********************************************************************/

SpUnit_matrix::SpUnit_matrix(): m_(0), n_(0), nnz_(0), ij_(NULL)
{	// constructeur par defaut, donne une matrice nulle
}

SpUnit_matrix::SpUnit_matrix(const int &m, const int &n, const int &nnz): m_(m), n_(n), nnz_(nnz), ij_(NULL) {
	ij_ = new int[n + m + 4 * nnz] ;
	if (ij_ == NULL) {
		_LogMessage("Sparse_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
		assert(false);
	}
	int i ;
	for (i = 0 ; i < n + m ; i ++) ij_[i] = 0 ; // ni_[] et nj_[]
}

SpUnit_matrix::SpUnit_matrix(const Sparse_matrix &S): m_(S.m_), n_(S.n_), nnz_(S.ntnz_), ij_(NULL) {
	assert(nnz_ != 0) ;
	int k, nt = n_ + m_ + 4 * nnz_ ;
	ij_ = new int[nt] ;
	if (ij_ == NULL) {
		_LogMessage("SpUnit_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
		assert(false);
	}
	int * S_ij = S.ij_ ; nt = n_ + m_ + nnz_ ; int * U_ij = ij_ ; for (k = 0 ; k < nt ; k ++)  U_ij[k] = S_ij[k] ; // copie ni_[], nj_[] et cmi[]
	U_ij += nt ; S_ij += nt + S.npnz_ - nnz_ ; for (k = 0 ; k < nnz_ ; k ++)   U_ij[k] = S_ij[k] ; // copie rmj_[]
	double * S_v = S.v_ ; U_ij += nnz_ ;
	for (k = 0 ; k < nnz_ ; k ++) {
		if (S_v[k] == 1.) U_ij[k] = true ;
		else {
			if (S_v[k] == -1.) U_ij[k] = false ;
			else {
				(void)sprintf(message, "SpUnit_matrix(Sparse_matrix &S): constructor failed: S.cmv[%d] is not 0, 1 or -1", k+1) ; _LogMessage(message) ;
				assert(false) ;
			}
		}
	}
	S_v += S.npnz_ ; U_ij += nnz_ ;
	for (k = 0 ; k < nnz_ ; k ++) {
		if (S_v[k] == 1.) U_ij[k] = true ;
		else {
			if (S_v[k] == -1.) U_ij[k] = false ;
			else {
				(void)sprintf(message, "SpUnit_matrix(Sparse_matrix &S): constructor failed: S.rmv[%d] is not 0, 1 or -1", k+1) ; _LogMessage(message) ;
				assert(false) ;
			}
		}
	}
}

SpUnit_matrix::SpUnit_matrix(const SpUnit_matrix &U): m_(U.m_), n_(U.n_), nnz_(U.nnz_), ij_(NULL) {
	assert(nnz_ != 0) ;
	int k, nt = n_ + m_ + 4 * nnz_ ;
	ij_ = new int[nt] ;
	if (ij_ == NULL) {
		_LogMessage("SpUnit_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
		assert(false);
	}
	int * temp_ij = U.ij_ ;
	for (k = 0 ; k < nt ; k ++)  ij_[k] = temp_ij[k] ;
}

SpUnit_matrix & SpUnit_matrix::operator=(const SpUnit_matrix &U) {
	assert(U.nnz_ != 0) ;
	int nt = U.n_ + U.m_ + 4 * U.nnz_ ;
	if (nt != n_ + m_ + 4 * nnz_) {
		if (ij_ != NULL) { Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; }
		ij_ = new int[nt];
	}
	if (ij_ == NULL) {
		_LogMessage("SpUnit_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
		assert(false);
	}
	m_ = U.m_ ; n_ = U.n_ ; nnz_ = U.nnz_ ;
	int k ; int * temp_ij = U.ij_ ;
	for (k = 0 ; k < nt ; k ++)  ij_[k] = temp_ij[k] ;
	return *this;
}

SpUnit_matrix::~SpUnit_matrix() {
	Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ;
}

int SpUnit_matrix::nblin() const { return m_; }

int SpUnit_matrix::nbcol() const { return n_; }

int SpUnit_matrix::nbnz() const { return nnz_; }

void SpUnit_matrix::set(const SpUnit_matrix &U) {
	if ((U.m_ != m_) || (U.n_ != n_) || (U.nnz_ != nnz_)) assert(false);
	int k, nt = n_ + m_ + 4 * nnz_ ; int * temp_ij = U.ij_ ;
	for (k = 0 ; k < nt ; k ++)  ij_[k] = temp_ij[k] ;
}

SpUnit_matrix transpose(const SpUnit_matrix &U) {
	int nnz = U.nnz_ ;
	assert(nnz != 0) ;
	SpUnit_matrix T ;
	int k, nt = U.n_ + U.m_ + 4 * nnz ;
	T.ij_ = new int[nt] ;
	if (T.ij_ == NULL) {
		_LogMessage("SpUnit_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
		assert(false);
	}
	T.m_ = U.n_ ; T.n_ = U.m_ ; T.nnz_ = nnz ;
	int * temp_U, * temp_T ;
	temp_U = U.ij_ ; temp_T = T.ij_ + T.n_ ; nt = T.m_ ;
	for (k = 0 ; k < nt ; k ++)  temp_T[k] = temp_U[k] ;
	temp_U += U.n_ ; temp_T = T.ij_ ; nt = T.n_ ;
	for (k = 0 ; k < nt ; k ++)  temp_T[k] = temp_U[k] ;
	temp_U += U.m_ ; temp_T = T.ij_ + T.n_ + T.m_ + nnz ;
	for (k = 0 ; k < nnz ; k ++)  temp_T[k] = temp_U[k] ;
	temp_U += nnz ; temp_T = T.ij_ + T.n_ + T.m_ ;
	for (k = 0 ; k < nnz ; k ++)  temp_T[k] = temp_U[k] ;
	temp_U += nnz ; temp_T = T.ij_ + T.n_ + T.m_ + 3 * nnz ;
	for (k = 0 ; k < nnz ; k ++)  temp_T[k] = temp_U[k] ;
	temp_U += nnz ; temp_T = T.ij_ + T.n_ + T.m_ + 2 * nnz ;
	for (k = 0 ; k < nnz ; k ++)  temp_T[k] = temp_U[k] ;
	return T;
}

SpUnit_matrix SpUnit_matrix::operator-() {
	//   - S (operateur unaire de changement de signe)
	assert(nnz_ != 0) ;
	SpUnit_matrix T ;
	int k, nt = n_ + m_ + 4 * nnz_ ;
	T.ij_ = new int[nt] ;
	if (T.ij_ == NULL) {
		_LogMessage("SpUnit_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
		assert(false);
	}
	T.m_ = m_ ; T.n_ = n_ ; T.nnz_ = nnz_ ;
	int * temp_ij = T.ij_ ; nt = n_ + m_ + nnz_ + nnz_;
	for (k = 0 ; k < nt ; k ++)  temp_ij[k] = ij_[k] ;
	nt = n_ + m_ + 4 * nnz_ ;
	for (k = nt-nnz_-nnz_ ; k < nt ; k ++) temp_ij[k] = !(ij_[k]) ;
	return T;
}

Sparse_matrix SpUnit_matrix::operator+(const Sparse_matrix &S) {
	//   M + S
	Sparse_matrix M(m_, n_);
	M.add_add(S, (*this), 0, false);
	return M;
}

Sparse_matrix SpUnit_matrix::operator+(const SpUnit_matrix &U) {
	//   M + U
	Sparse_matrix M(m_, n_);
	M.add_add((*this), U, 0, false);
	return M;
}

Sparse_matrix SpUnit_matrix::operator-(const Sparse_matrix &S) {
	//   M - S
	Sparse_matrix M(m_, n_);
	M.add_sub((*this), S, 0, false);
	return M;
}

Sparse_matrix SpUnit_matrix::operator-(const SpUnit_matrix &U) {
	//   M - U
	Sparse_matrix M(m_, n_);
	M.add_sub((*this), U, 0, false);
	return M;
}

Sparse_matrix SpUnit_matrix::operator*(const Sparse_matrix &S) {
	//   M * S  elementwise
	Sparse_matrix M(m_, n_);
	M.add_elemult((*this), S, 0, false);
	return M;
}

SpUnit_matrix SpUnit_matrix::operator*(const SpUnit_matrix &U) {
	//   M * U  elementwise
	SpUnit_matrix M(m_, n_);
	M.add_elemult(U, (*this), 0, false);
	return M;
}

Sparse_matrix operator*(const double &a, const SpUnit_matrix &U) {
	// a * U : multiplication par un scalaire,
	Sparse_matrix M(U.m_, U.n_, U.nnz_) ;
	M.add_mult(a, U, 0, false);
	return M;
}

Sparse_matrix SpUnit_matrix::operator/(const double &a) {
	//   U / a : division par un scalaire
	Sparse_matrix M(m_, n_, nnz_) ;
	M.add_div((*this), a, 0, false);
	return M;
}

Sparse_matrix SpUnit_matrix::operator/(const Fortran_vector &v) {
	// U / v : division elementwise par colonne, par un vecteur de taille m = U.nblin()
	Sparse_matrix M(m_, n_, nnz_) ;
	M.add_elediv((*this), v, 0, false);
	return M;
}

Sparse_matrix row_elemult(const SpUnit_matrix &U, const Fortran_vector &v) {
	// U * v : multiplication elementwise PAR LIGNE, par un vecteur de taille n = U.nbcol()
	Sparse_matrix M(U.m_, U.n_, U.nnz_) ;
	M.add_row_elemult(U, v, 0, false);
	return M;
}

Sparse_matrix row_elediv(const SpUnit_matrix &U, const Fortran_vector &v) {
	// U / v : division elementwise PAR LIGNE, par un vecteur de taille n = S.nbcol()
	Sparse_matrix M(U.m_, U.n_, U.nnz_) ;
	M.add_row_elediv(U, v, 0, false);
	return M;
}

Fortran_vector matmult(const SpUnit_matrix &U, const Fortran_vector &v) {
	//   S * v, multiplication matricielle
	Fortran_vector y(U.m_) ;
	y.set_matmult(U, v) ;
	return y;
}

Sparse_matrix matmult(const SpUnit_matrix &U, const Sparse_matrix &S, Fortran_vector * temp_v) {
	// si temp_v est fourni, il doit pointer sur un vecteur  de taille  n = U.n_  qui servira de tampon.
	Sparse_matrix TM(U.m_, S.n_) ;
	TM.set_matmult(U, S, temp_v, false) ;
	return TM;
}
Sparse_matrix matmult(const Sparse_matrix &S, const SpUnit_matrix &U, Fortran_vector * temp_v) {
	// si temp_v est fourni, il doit pointer sur un vecteur  de taille  n = S.n_  qui servira de tampon.
	Sparse_matrix TM(S.m_, U.n_) ;
	TM.set_matmult(S, U, temp_v, false) ;
	return TM;
}

Sparse_matrix matmult(const SpUnit_matrix &U1, const SpUnit_matrix &U2, Index_vector * temp_v) {
	// si temp_v est fourni, il doit pointer sur un vecteur d'entiers de taille  n = U1.n_  qui servira de tampon.
	Sparse_matrix TM(U1.m_, U2.n_) ;
	TM.set_matmult(U1, U2, temp_v, false) ;
	return TM;
}

SpUnit_matrix & SpUnit_matrix::operator*=(const SpUnit_matrix &U) {
	//  M *= U ,  elementwise ;
	if ((U.n_ != n_) || (U.m_ != m_)) assert(false);
	int i, j ;
	bool Sparse_elemult_set_ij_done = Sparse_elemult_set_ij_check(ij_, ij_, U.ij_) ;
	if (!(Sparse_elemult_set_ij_done)) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		(*this).add_elemult(*this, U, 0) ;
	}
	else {
		int * nj_ = ij_ + n_ ; int * nj_U = U.ij_ + n_  ;
		int * rmj_ = nj_ + m_ + nnz_ ; int * rmj_U = nj_U + m_ + U.nnz_ ; int * rmj_max, * rmjU_max ;
		int * rmv_ = nj_ + m_ + 3 * nnz_ ; int * rmv_U = rmj_U + 2 * U.nnz_ ;
		for (i = 0 ; i < m_ ; i ++) {
			rmj_max = rmj_ + nj_[i] - 1 ; rmjU_max = rmj_U + nj_U[i] - 1 ; j = (*rmj_) ;
			while (rmj_ <= rmj_max) {
				if (j == (*rmj_U)) {
					if (*rmv_)  *rmv_ = (*rmv_U != 0) ;	else  *rmv_ = (*rmv_U == 0) ;
					rmv_ ++ ; rmv_U ++ ; rmj_ ++ ; j = (*rmj_) ; rmj_U ++ ;
				}
				else { rmv_U ++ ; rmj_U ++ ; }
			}
			rmv_U += rmjU_max + 1 - rmj_U ; rmj_U = 1 + rmjU_max ;
		}
		int * ni_ = ij_ ; int * ni_U = U.ij_  ;
		int * cmi_ = ni_ + m_ + n_ ; int * cmi_U = ni_U + m_ + n_ ;	int * cmi_max, * cmiU_max ;
		int * cmv_ = nj_ + m_ + 2 * nnz_ ; int * cmv_U = cmi_U + 2 * U.nnz_ ;
		for (j = 0 ; j < n_ ; j ++) {
			cmi_max = cmi_ + ni_[j] - 1 ; cmiU_max = cmi_U + ni_U[j] - 1 ; i = (*cmi_) ;
			while (cmi_ <= cmi_max) {
				if (i == (*cmi_U)) {
					if (*cmv_)  *cmv_ = (*cmv_U != 0) ;	else  *cmv_ = (*cmv_U == 0) ;
					cmv_ ++ ; cmv_U ++ ; cmi_ ++ ; i = (*cmi_) ; cmi_U ++ ;
				}
				else { cmv_U ++ ; cmi_U ++ ; }
			}
			cmv_U += cmiU_max + 1 - cmi_U ; cmi_U = 1 + cmiU_max ;
		}
	}
	return *this;
}

void SpUnit_matrix::add_elemult(const SpUnit_matrix &U1, const SpUnit_matrix &U2, int ad, bool save_count) {
	if ((n_ != U1.n_) || (m_ != U1.m_) || (n_ != U2.n_) || (m_ != U2.m_)) assert(false);
	if (ad == 0) {	// seule la fonction  set_elemult(U1,U2), soit (*this) = U1*U2, est implementee
		bool Sparse_elemult_set_ij_done = Sparse_elemult_set_ij_check(ij_, U1.ij_, U2.ij_) ;
		if (! Sparse_elemult_set_ij_done) {
			Sparse_XXX_set_ij_uncount(ij_) ;
			(*this).SpUnit_elemult_set_ij_set(U1, U2);
			(*this).add_elemult(U1, U2, 0, false) ;
			if (! save_count)  Sparse_elemult_set_ij_uncount(ij_, U1.ij_, U2.ij_) ;
		}
		else {
			int i, j, i1, i2, j1, j2 ;
			int * nj_ = ij_ + n_ ; int * nj_1 = U1.ij_ + n_  ; int * nj_2 = U2.ij_ + n_  ;
			int * rmj_1 = nj_1 + m_ + U1.nnz_ ;	int * rmj_2 = nj_2 + m_ + U2.nnz_ ;	int * rmj1_max, * rmj2_max ;
			int * rmv_ = nj_ + m_ + 3 * nnz_ ; int * rmv_1 = rmj_1 + 2 * U1.nnz_ ; int * rmv_2 = rmj_2 + 2 * U2.nnz_ ;
			for (i = 0 ; i < m_ ; i ++) {
				rmj1_max = rmj_1 + nj_1[i] - 1 ; rmj2_max = rmj_2 + nj_2[i] - 1 ;
				j1 = (*rmj_1) ; j2 = (*rmj_2) ;
				while ((rmj_1 <= rmj1_max) && (rmj_2 <= rmj2_max)) {
					if (j1 == j2) {
						if (*rmv_1)  *rmv_ = (*rmv_2 != 0) ;  else  *rmv_ = (*rmv_2 == 0) ;
						rmv_ ++ ; rmv_1 ++ ; rmv_2 ++ ; rmj_1 ++ ; j1 = (*rmj_1) ; rmj_2 ++ ; j2 = (*rmj_2) ;
					} else
						if (j1 < j2) { rmv_1 ++ ; rmj_1 ++ ; j1 = (*rmj_1) ; }
						else		 { rmv_2 ++ ; rmj_2 ++ ; j2 = (*rmj_2) ; }
				}
				rmv_1 += rmj1_max + 1 - rmj_1 ; rmj_1 = 1 + rmj1_max ;
				rmv_2 += rmj2_max + 1 - rmj_2 ; rmj_2 = 1 + rmj2_max ;
			}
//			int * ni_ = ij_ ;
            int * ni_1 = U1.ij_ ;	int * ni_2 = U2.ij_  ;
			int * cmi_1 = ni_1 + m_ + n_ ; int * cmi_2 = ni_2 + m_ + n_ ; int * cmi1_max, * cmi2_max ;
			int * cmv_ = nj_ + m_ + 2 * nnz_ ; int * cmv_1 = cmi_1 + 2 * U1.nnz_ ; int * cmv_2 = cmi_2 + 2 * U2.nnz_ ;
			for (j = 0 ; j < n_ ; j ++) {
				cmi1_max = cmi_1 + ni_1[j] - 1 ; cmi2_max = cmi_2 + ni_2[j] - 1 ;
				i1 = (*cmi_1) ; i2 = (*cmi_2) ;
				while ((cmi_1 <= cmi1_max) && (cmi_2 <= cmi2_max)) {
					if (i1 == i2) {
						if (*cmv_1)  *cmv_ = (*cmv_2 != 0) ;  else  *cmv_ = (*cmv_2 == 0) ;
						cmv_ ++ ; cmv_1 ++ ; cmv_2 ++ ; cmi_1 ++ ; i1 = (*cmi_1) ; cmi_2 ++ ; i2 = (*cmi_2) ;
					} else
						if (i1 < i2) { cmv_1 ++ ; cmi_1 ++ ; i1 = (*cmi_1) ; }
						else		 { cmv_2 ++ ; cmi_2 ++ ; i2 = (*cmi_2) ; }
				}
				cmv_1 += cmi1_max + 1 - cmi_1 ; cmi_1 = 1 + cmi1_max ;
				cmv_2 += cmi2_max + 1 - cmi_2 ; cmi_2 = 1 + cmi2_max ;
			}
		}
	}
	else assert(false) ;// seule la fonction  set_elemult(U1,U2), soit (*this) = U1*U2, est implementee
}


/**** les fonctions suivantes sont appelees par les fonctions 'add_XXX()' ************************************************/

bool Sparse_add_set_ij_check(int * ij_, int * ij_1, int * ij_2) { // indique si la matrice(ij_) a deja ete calculee comme somme/diff. des matrices (ij_1, ij_2)
	int i, k ;
	if (Sparse_add_set_ij_count == 0)
		Sparse_add_set_ij_ = new int* [30] ;
	bool Sparse_add_set_ij_done = false ;
	k = 3 * Sparse_add_set_ij_count ;
	for (i = 0 ; i < k ; i += 3) {
		if (Sparse_add_set_ij_[i] == ij_) {
			if (((Sparse_add_set_ij_[i+1] == ij_1) && (Sparse_add_set_ij_[i+2] == ij_2)) || ((Sparse_add_set_ij_[i+2] == ij_1) && (Sparse_add_set_ij_[i+1] == ij_2))) {
				Sparse_add_set_ij_done = true ; break ;
			}
		}
	}
	return Sparse_add_set_ij_done ;
}

void Sparse_matrix::Sparse_add_set_ij_set(const Sparse_matrix & S1, const Sparse_matrix & S2) {
	int * S1_ij = S1.ij_ ; int S1_npnz = S1.npnz_ ; int * S2_ij = S2.ij_ ; int S2_npnz = S2.npnz_ ;
	int ntnz = Sparse_add_set_nnz_(S1_ij, S2_ij, m_, n_, S1_npnz, S2_npnz) ;
	if (npnz_ < ntnz) {
		npnz_ = ntnz ;
		int * save_ij = ij_ ; double * save_v = v_ ;
		ij_ = new int[n_ + m_ + npnz_ + npnz_] ;	v_ = new double[npnz_ + npnz_] ;
		if (( v_ == NULL ) || (ij_ == NULL)) {
		_LogMessage("Sparse_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
			assert(false) ;
		}
		Sparse_XXX_set_ij_uncount(save_ij) ; delete [] save_ij ; delete [] save_v ;
	}
	ntnz_ = ntnz ;
	Sparse_add_fill_ij_(ij_, S1_ij, S2_ij, m_, n_, npnz_, S1_npnz, S2_npnz) ;
	Sparse_add_set_ij_count_update(ij_, S1_ij, S2_ij) ;
}

int Sparse_add_set_nnz_(int * ij1, int * ij2, const int & m, const int & n, const int & npnz1, const int & npnz2) { // donne ntnz de la matrice resultante
	int i, j1, j2, ntnz = 0 ;
	int * nj_1 = ij1 + n  ;	int * nj_2 = ij2 + n  ;
	int * rmj_1 = nj_1 + m + npnz1 ;	int * rmj_2 = nj_2 + m + npnz2 ;	int * rmj1_max, * rmj2_max ;
	for (i = 0 ; i < m ; i ++) {
		rmj1_max = rmj_1 + nj_1[i] - 1 ; rmj2_max = rmj_2 + nj_2[i] - 1 ;
		j1 = (*rmj_1) ; j2 = (*rmj_2) ;
		while ((rmj_1 <= rmj1_max) && (rmj_2 <= rmj2_max)) {
			if (j1 == j2) {
				ntnz++ ; rmj_2 ++ ; j2 = (*rmj_2) ; rmj_1 ++ ; j1 = (*rmj_1) ;
			} else
				if (j1 < j2) { ntnz++ ; rmj_1 ++ ; j1 = (*rmj_1) ; }
				else { ntnz++ ; rmj_2 ++ ; j2 = (*rmj_2) ; }
		}
		ntnz += rmj1_max - rmj_1 + 1 ; rmj_1 = 1 + rmj1_max ;
		ntnz += rmj2_max - rmj_2 + 1 ; rmj_2 = 1 + rmj2_max ;
	}
	return ntnz ;
}

void Sparse_add_fill_ij_(int* ij, int*ij1, int* ij2, const int& m, const int& n, const int& npnz, const int& npnz1, const int& npnz2) {
	int i, j, j1, j2, i1, i2, ni, nj ;
	int * nj_ = ij + n ; int * nj_1 = ij1 + n ; int * nj_2 = ij2 + n ; int * rmj1_max, * rmj2_max ;
	int * rmj_ = nj_ + m + npnz ; int * rmj_1 = nj_1 + m + npnz1 ; int * rmj_2 = nj_2 + m + npnz2 ;
	for (i = 0 ; i < m ; i ++) {
		rmj1_max = rmj_1 + nj_1[i] - 1 ; rmj2_max = rmj_2 + nj_2[i] - 1 ;
		nj = 0 ; j1 = (*rmj_1) ; j2 = (*rmj_2) ;
		while ((rmj_1 <= rmj1_max) && (rmj_2 <= rmj2_max)) {
			if (j1 == j2) {
				*rmj_ = j2 ; rmj_ ++ ; nj ++ ; rmj_2 ++ ; j2 = (*rmj_2) ; rmj_1 ++ ; j1 = (*rmj_1) ;
			} else
				if (j1 < j2) { *rmj_ = j1 ; rmj_ ++ ; nj ++ ; rmj_1 ++ ; j1 = (*rmj_1) ; }
				else { *rmj_ = j2 ; rmj_ ++ ; nj ++ ; rmj_2 ++ ; j2 = (*rmj_2) ; }
		}
		while (rmj_1 <= rmj1_max) { *rmj_ = (*rmj_1) ; rmj_ ++ ; nj ++ ; rmj_1 ++ ;	}
		while (rmj_2 <= rmj2_max) {	*rmj_ = (*rmj_2) ; rmj_ ++ ; nj ++ ; rmj_2 ++ ; }
		(*nj_) = nj ; nj_ ++ ;
	}
	int * ni_ = ij ; int * ni_1 = ij1 ; int * ni_2 = ij2 ; int * cmi1_max, * cmi2_max ;
	int * cmi_ = ni_ + m + n ; int * cmi_1 = ni_1 + m + n ; int * cmi_2 = ni_2 + m + n ;
	for (j = 0 ; j < n ; j ++) {
		cmi1_max = cmi_1 + ni_1[j] - 1 ; cmi2_max = cmi_2 + ni_2[j] - 1 ;
		ni = 0 ; i1 = (*cmi_1) ; i2 = (*cmi_2) ;
		while ((cmi_1 <= cmi1_max) && (cmi_2 <= cmi2_max)) {
			if (i1 == i2) {
				*cmi_ = i2 ; cmi_ ++ ; ni ++ ; cmi_2 ++ ; i2 = (*cmi_2) ; cmi_1 ++ ; i1 = (*cmi_1) ;
			} else
				if (i1 < i2) { *cmi_ = i1 ; cmi_ ++ ; ni ++ ; cmi_1 ++ ; i1 = (*cmi_1) ; }
				else { *cmi_ = i2 ; cmi_ ++ ; ni ++ ; cmi_2 ++ ; i2 = (*cmi_2) ; }
		}
		while (cmi_1 <= cmi1_max) { *cmi_ = (*cmi_1) ; cmi_ ++ ; ni ++ ; cmi_1 ++ ; }
		while (cmi_2 <= cmi2_max) { *cmi_ = (*cmi_2) ; cmi_ ++ ; ni ++ ; cmi_2 ++ ; }
		(*ni_) = ni ; ni_ ++ ;
	}
}

void Sparse_add_set_ij_count_update(int * ij_, int * ij_1, int * ij_2) {
	int k = 3 * Sparse_add_set_ij_count ;
	if (Sparse_add_set_ij_count == Sparse_add_set_ij_max) {
		cout << "on a atteint Sparse_add_set_ij_max : extension de la zone reservee correspondante" << endl ;
		Sparse_add_set_ij_max += 10 ;
		int ** temp = new int* [3 * Sparse_add_set_ij_max] ;
		for (int i = 0 ; i < k ; i ++)  temp[i] = Sparse_add_set_ij_[i] ;
		delete [] Sparse_add_set_ij_ ;
		Sparse_add_set_ij_ = temp ;
	}
	Sparse_add_set_ij_count ++ ;
	if (Sparse_add_set_ij_count > 100)
		{ _LogMessage("ATTENTION : Sparse_add_set_ij_count > 100 !!!") ; assert(false) ; }
	Sparse_add_set_ij_[k] = ij_ ; Sparse_add_set_ij_[k+1] = ij_1 ; Sparse_add_set_ij_[k+2] = ij_2 ;
}

void Sparse_XXX_set_ij_uncount(int * ij_) { // dereferencie ij_ de tous les registres  Sparse_XXX_set_ij : a appeler  juste avant  tout  'delete[ ] ij_'
	int i, j, *r_count, **reg, **reg_count = new int*[3], *** registres = new int**[3] ;
	reg_count[0] = &Sparse_add_set_ij_count ; reg_count[1] = &Sparse_elemult_set_ij_count ; reg_count[2] = &Sparse_matmult_set_ij_count ;
	registres[0] = Sparse_add_set_ij_ ; registres[1] = Sparse_elemult_set_ij_ ; registres[2] = Sparse_matmult_set_ij_ ;
	for(j = 0 ; j < 3 ; j ++) {
		reg = registres[j] ; r_count = reg_count[j] ;
		for (i = 3 * (*r_count) - 3 ; i >= 0 ; i -= 3) {
			if ((reg[i] == ij_) || (reg[i+1] == ij_) || (reg[i+2] == ij_)) {
				reg[i] = reg[3 * (*r_count) - 3] ;
				reg[i + 1] = reg[3 * (*r_count) - 2] ;
				reg[i + 2] = reg[3 * (*r_count) - 1] ;
				(*r_count) -- ;
			}
		}
	}
	delete[ ] reg_count ; delete[ ] registres ;
	reg = Sparse_addsubdiag_set_ij_ ; r_count = &Sparse_addsubdiag_set_ij_count ;// int*reg_i = Sparse_addsubdiag_i_set_ij_ ;
	for (i = (*r_count) - 1 ; i >= 0 ; i --) {
		if (reg[i] == ij_) {
			reg[i] = reg[(*r_count) - 1] ; Sparse_addsubdiag_i_set_ij_[i+i] = Sparse_addsubdiag_i_set_ij_[2 * (*r_count) - 2] ; Sparse_addsubdiag_i_set_ij_[i+i+1] = Sparse_addsubdiag_i_set_ij_[2 * (*r_count) - 1] ;
			(*r_count) -- ;
		}
	}
}

void Sparse_add_set_ij_uncount(int * ij_, int * ij_1, int * ij_2) {
	int i, ii = 0 ;
	bool Sparse_combination_found = false ;
	for (i = 3 * Sparse_add_set_ij_count - 3 ; i >= 0 ; i -= 3) {
		if (Sparse_add_set_ij_[i] == ij_) {
			if (((Sparse_add_set_ij_[i+1] == ij_1) && (Sparse_add_set_ij_[i+2] == ij_2)) || ((Sparse_add_set_ij_[i+2] == ij_1) && (Sparse_add_set_ij_[i+1] == ij_2))) {
				Sparse_combination_found = true ; ii = i ;break ;
			}
		}
	}
	if (! Sparse_combination_found) assert(false) ;
	if (ii != 3 * Sparse_add_set_ij_count - 3) {
		Sparse_add_set_ij_[ii] = Sparse_add_set_ij_[3 * Sparse_add_set_ij_count - 3] ;
		Sparse_add_set_ij_[ii + 1] = Sparse_add_set_ij_[3 * Sparse_add_set_ij_count - 2] ;
		Sparse_add_set_ij_[ii + 2] = Sparse_add_set_ij_[3 * Sparse_add_set_ij_count - 1] ;
	}
	Sparse_add_set_ij_count -- ;
}

bool Sparse_elemult_set_ij_check(int * ij_, int * ij_1, int * ij_2) {
	int i, k ;
	if (Sparse_elemult_set_ij_count == 0)
		Sparse_elemult_set_ij_ = new int* [30] ;
	bool Sparse_elemult_set_ij_done = false ;
	k = 3 * Sparse_elemult_set_ij_count ;
	for (i = 0 ; i < k ; i += 3) {
		if (Sparse_elemult_set_ij_[i] == ij_) {
			if (((Sparse_elemult_set_ij_[i+1] == ij_1) && (Sparse_elemult_set_ij_[i+2] == ij_2)) || ((Sparse_elemult_set_ij_[i+2] == ij_1) && (Sparse_elemult_set_ij_[i+1] == ij_2))) {
				Sparse_elemult_set_ij_done = true ; break ;
			}
		}
	}
	return Sparse_elemult_set_ij_done ;
}

void SpUnit_matrix::SpUnit_elemult_set_ij_set(const SpUnit_matrix & U1, const SpUnit_matrix & U2) {
	int * U1_ij = U1.ij_ ; int U1_nnz = U1.nnz_ ; int * U2_ij = U2.ij_ ; int U2_nnz = U2.nnz_ ;
	int nnz = Sparse_elemult_set_nnz_(U1_ij, U2_ij, m_, n_, U1_nnz, U2_nnz) ;
	if (nnz != nnz_) {
		int * save_ij = ij_ ; nnz_ = nnz ;
		ij_ = new int[n_ + m_ + 4 * nnz_] ;
		if (ij_ == NULL) {
			_LogMessage("SpUnit_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
			assert(false) ;
		}
		Sparse_XXX_set_ij_uncount(save_ij) ; delete [] save_ij ;
	}
	Sparse_elemult_fill_ij_(ij_, U1_ij, U2_ij, m_, n_, nnz_, U1_nnz, U2_nnz) ;
	Sparse_elemult_set_ij_count_update(ij_, U1_ij, U2_ij) ;
}

void Sparse_matrix::Sparse_elemult_set_ij_set(const Sparse_matrix & S1, const Sparse_matrix & S2) {
	int * S1_ij = S1.ij_ ; int S1_npnz = S1.npnz_ ; int * S2_ij = S2.ij_ ; int S2_npnz = S2.npnz_ ;
	int ntnz = Sparse_elemult_set_nnz_(S1_ij, S2_ij, m_, n_, S1_npnz, S2_npnz) ;
	if (npnz_ < ntnz) {
		int * save_ij = ij_ ; double * save_v = v_ ;
		npnz_ = ntnz ;
		ij_ = new int[n_ + m_ + npnz_ + npnz_] ;	v_ = new double[npnz_ + npnz_] ;
		if (( v_ == NULL ) || (ij_ == NULL)) {
			_LogMessage("SpUnit_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
			assert(false) ;
		}
		Sparse_XXX_set_ij_uncount(save_ij) ; delete [] save_ij ; delete [] save_v ;
	}
	ntnz_ = ntnz ;
	Sparse_elemult_fill_ij_(ij_, S1_ij, S2_ij, m_, n_, npnz_, S1_npnz, S2_npnz) ;
	Sparse_elemult_set_ij_count_update(ij_, S1_ij, S2_ij) ;
}

int Sparse_elemult_set_nnz_(int * ij1, int * ij2, const int & m, const int & n, const int & npnz1, const int & npnz2) {
	int i, j1, j2, ntnz = 0 ;
	int * nj_1 = ij1 + n  ;	int * nj_2 = ij2 + n  ;
	int * rmj_1 = nj_1 + m + npnz1 ;	int * rmj_2 = nj_2 + m + npnz2 ;	int * rmj1_max, * rmj2_max ;
	for (i = 0 ; i < m ; i ++) {
		rmj1_max = rmj_1 + nj_1[i] - 1 ; rmj2_max = rmj_2 + nj_2[i] - 1 ;
		j1 = (*rmj_1) ; j2 = (*rmj_2) ;
		while ((rmj_1 <= rmj1_max) && (rmj_2 <= rmj2_max)) {
			if (j1 == j2) {
				ntnz++ ; rmj_2 ++ ; j2 = (*rmj_2) ; rmj_1 ++ ; j1 = (*rmj_1) ;
			} else
				if (j1 < j2) { rmj_1 ++ ; j1 = (*rmj_1) ; }
				else		 { rmj_2 ++ ; j2 = (*rmj_2) ; }
		}
		rmj_1 = 1 + rmj1_max ;
		rmj_2 = 1 + rmj2_max ;
	}
	return ntnz ;
}

void Sparse_elemult_fill_ij_(int* ij, int*ij1, int* ij2, const int& m, const int& n, const int& npnz, const int& npnz1, const int& npnz2) {
	int i, j, j1, j2, i1, i2, ni, nj ;
	int * nj_ = ij + n ; int * nj_1 = ij1 + n ; int * nj_2 = ij2 + n ; int * rmj1_max, * rmj2_max ;
	int * rmj_ = nj_ + m + npnz ; int * rmj_1 = nj_1 + m + npnz1 ; int * rmj_2 = nj_2 + m + npnz2 ;
	for (i = 0 ; i < m ; i ++) {
		rmj1_max = rmj_1 + nj_1[i] - 1 ; rmj2_max = rmj_2 + nj_2[i] - 1 ;
		nj = 0 ; j1 = (*rmj_1) ; j2 = (*rmj_2) ;
		while ((rmj_1 <= rmj1_max) && (rmj_2 <= rmj2_max)) {
			if (j1 == j2) {
				*rmj_ = j2 ; rmj_ ++ ; nj ++ ; rmj_2 ++ ; j2 = (*rmj_2) ; rmj_1 ++ ; j1 = (*rmj_1) ;
			} else
				if (j1 < j2) { rmj_1 ++ ; j1 = (*rmj_1) ; }
				else		 { rmj_2 ++ ; j2 = (*rmj_2) ; }
		}
		rmj_1 = 1 + rmj1_max ;
		rmj_2 = 1 + rmj2_max ;
		(*nj_) = nj ; nj_ ++ ;
	}
	int * ni_ = ij ; int * ni_1 = ij1 ; int * ni_2 = ij2 ; int * cmi1_max, * cmi2_max ;
	int * cmi_ = ni_ + m + n ; int * cmi_1 = ni_1 + m + n ; int * cmi_2 = ni_2 + m + n ;
	for (j = 0 ; j < n ; j ++) {
		cmi1_max = cmi_1 + ni_1[j] - 1 ; cmi2_max = cmi_2 + ni_2[j] - 1 ;
		ni = 0 ; i1 = (*cmi_1) ; i2 = (*cmi_2) ;
		while ((cmi_1 <= cmi1_max) && (cmi_2 <= cmi2_max)) {
			if (i1 == i2) {
				*cmi_ = i2 ; cmi_ ++ ; ni ++ ; cmi_2 ++ ; i2 = (*cmi_2) ; cmi_1 ++ ; i1 = (*cmi_1) ;
			} else
				if (i1 < i2) { cmi_1 ++ ; i1 = (*cmi_1) ; }
				else		 { cmi_2 ++ ; i2 = (*cmi_2) ; }
		}
		cmi_1 = 1 + cmi1_max ;
		cmi_2 = 1 + cmi2_max ;
		(*ni_) = ni ; ni_ ++ ;
	}
}

void Sparse_elemult_set_ij_count_update(int * ij_, int * ij_1, int * ij_2) {
	int k = 3 * Sparse_elemult_set_ij_count ;
	if (Sparse_elemult_set_ij_count == Sparse_elemult_set_ij_max) {
		Sparse_elemult_set_ij_max += 10 ;
		int ** temp = new int* [3 * Sparse_elemult_set_ij_max] ;
		for (int i = 0 ; i < k ; i ++)  temp[i] = Sparse_elemult_set_ij_[i] ;
		delete [] Sparse_elemult_set_ij_ ;
		Sparse_elemult_set_ij_ = temp ;
	}
	Sparse_elemult_set_ij_count ++ ;
	if (Sparse_elemult_set_ij_count > 100)
		{ _LogMessage("ATTENTION : Sparse_elemult_set_ij_count > 100 !!!") ; assert(false) ; }
	Sparse_elemult_set_ij_[k] = ij_ ; Sparse_elemult_set_ij_[k+1] = ij_1 ; Sparse_elemult_set_ij_[k+2] = ij_2 ;
}

void Sparse_elemult_set_ij_uncount(int * ij_, int * ij_1, int * ij_2) {
	int i, ii = 0 ;
	bool Sparse_combination_found = false ;
	for (i = 3 * Sparse_elemult_set_ij_count - 3 ; i >= 0 ; i -= 3) {
		if (Sparse_elemult_set_ij_[i] == ij_) {
			if (((Sparse_elemult_set_ij_[i+1] == ij_1) && (Sparse_elemult_set_ij_[i+2] == ij_2)) || ((Sparse_elemult_set_ij_[i+2] == ij_1) && (Sparse_elemult_set_ij_[i+1] == ij_2))) {
				Sparse_combination_found = true ; ii = i ;break ;
			}
		}
	}
	if (! Sparse_combination_found) assert(false) ;
	if (ii != 3 * Sparse_elemult_set_ij_count - 3) {
		Sparse_elemult_set_ij_[ii] = Sparse_elemult_set_ij_[3 * Sparse_elemult_set_ij_count - 3] ;
		Sparse_elemult_set_ij_[ii + 1] = Sparse_elemult_set_ij_[3 * Sparse_elemult_set_ij_count - 2] ;
		Sparse_elemult_set_ij_[ii + 2] = Sparse_elemult_set_ij_[3 * Sparse_elemult_set_ij_count - 1] ;
	}
	Sparse_elemult_set_ij_count -- ;
}

bool Sparse_matmult_set_ij_check(int * ij_, int * ij_1, int * ij_2) {
	int i, k ;
	if (Sparse_matmult_set_ij_count == 0)
		Sparse_matmult_set_ij_ = new int* [30] ;
	bool Sparse_matmult_set_ij_done = false ;
	k = 3 * Sparse_matmult_set_ij_count ;
	for (i = 0 ; i < k ; i += 3) {
		if (Sparse_matmult_set_ij_[i] == ij_) {
			if (((Sparse_matmult_set_ij_[i+1] == ij_1) && (Sparse_matmult_set_ij_[i+2] == ij_2)) || ((Sparse_matmult_set_ij_[i+2] == ij_1) && (Sparse_matmult_set_ij_[i+1] == ij_2))) {
				Sparse_matmult_set_ij_done = true ; break ;
			}
		}
	}
	return Sparse_matmult_set_ij_done ;
}

void Sparse_matrix::Sparse_matmult_set_ij_set(const Sparse_matrix & S1, const Sparse_matrix & S2, double * v) {
	int * S1_ij = S1.ij_ ; int S1_npnz = S1.npnz_ ; int * S2_ij = S2.ij_ ; int S2_npnz = S2.npnz_ ;
	int ntnz = Sparse_matmult_set_nnz_(S1_ij, S2_ij, m_, S1.n_, n_, S1_npnz, S2_npnz, v) ;
	if (npnz_ < ntnz) {
		if ((ij_ == S1_ij) || (ij_ == S2_ij)) assert(false) ;
		Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; delete [] v_ ;
		npnz_ = ntnz ;
		ij_ = new int[n_ + m_ + npnz_ + npnz_] ;	v_ = new double[npnz_ + npnz_] ;
		if (( v_ == NULL ) || (ij_ == NULL)) {
			_LogMessage("Sparse_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
			assert(false) ;
		}
	}
	ntnz_ = ntnz ;
	Sparse_matmult_fill_ij_(ij_, S1_ij, S2_ij, m_, S1.n_, n_, npnz_, S1_npnz, S2_npnz, v) ;
	Sparse_matmult_set_ij_count_update(ij_, S1_ij, S2_ij) ;
}

int Sparse_matmult_set_nnz_(int* ij1, int* ij2, const int &m, const int &n, const int &p, const int &npnz1, const int &npnz2, double * v) {
	int i, j, k, ntnz = 0 ;
	int * cmi_2 = ij2 + p + n ; int * rmj_1 ; int * rmj1_max, * cmi2_max ;
	int * ni_2 = ij2 - 1 ; int * nj_1 = ij1 + n - 1 ;
	for (j = 1 ; j <= p ; j ++) {
		for (k = 1 ; k <= n ; k++)	 v[k] = 0.;
		cmi2_max = cmi_2 + ni_2[j] - 1 ;
		while (cmi_2 <= cmi2_max)	{ v[*cmi_2] = 1. ; cmi_2 ++ ; }
		rmj_1 = ij1 + m + n + npnz1 ;
		for (i = 1 ; i <= m ; i ++) {
			rmj1_max = rmj_1 + nj_1[i] - 1 ;
			while (rmj_1 <= rmj1_max) {
				if (v[*rmj_1]) { ntnz ++ ; break ; }
				rmj_1 ++ ;
			}
			rmj_1 = 1 + rmj1_max ;
		}
	}
	return ntnz ;
}

void Sparse_matmult_fill_ij_(int* ij, int*ij1, int* ij2, const int &m, const int &n, const int &p, const int &npnz, const int &npnz1, const int &npnz2, double* v) {
	// v (requis) = un double* dimensionne a  (n = S1.n_) + 1
	int i, j, k ;
	int * ni_2 = ij2 - 1 ; int * nj_1 = ij1 + n - 1 ;
	int * ni_ = ij - 1 ; int * nj_ = ij + p - 1 ; int * cmi_ = ij + m + p ; int * rmj_ = cmi_ + npnz ;
	int * cmi_2 = ij2 + p + n ; int * rmj_1 ; int * rmj1_max, * cmi2_max ;
	for (j = 1 ; j <= p ; j ++) {
		ni_ ++ ; *ni_ = 0 ;
		for (k = 1 ; k <= n ; k++)	 v[k] = 0.;
		cmi2_max = cmi_2 + ni_2[j] - 1 ;
		while (cmi_2 <= cmi2_max)	{ v[*cmi_2] = 1. ; cmi_2 ++ ; }
		rmj_1 = ij1 + m + n + npnz1 ;
		for (i = 1 ; i <= m ; i ++) {
			rmj1_max = rmj_1 + nj_1[i] - 1 ;
			while (rmj_1 <= rmj1_max) {
				if (v[*rmj_1]) { (*ni_) ++ ; *cmi_ = i ; cmi_ ++ ; break ; }
				rmj_1 ++ ;
			}
			rmj_1 = 1 + rmj1_max ;
		}
	}
	rmj_1 = ij1 + m + n + npnz1 ;
	for (i = 1 ; i <= m ; i ++) {
		nj_ ++ ; *nj_ = 0 ;
		for (k = 1 ; k <= n ; k++)	 v[k] = 0.;
		rmj1_max = rmj_1 + nj_1[i] - 1 ;
		while (rmj_1 <= rmj1_max)	{ v[*rmj_1] = 1. ; rmj_1 ++ ; }
		cmi_2 = ij2 + p + n ;
		for (j = 1 ; j <= p ; j ++) {
			cmi2_max = cmi_2 + ni_2[j] - 1 ;
			while (cmi_2 <= cmi2_max) {
				if (v[*cmi_2]) { (*nj_) ++ ; *rmj_ = j ; rmj_ ++ ; break ; }
				cmi_2 ++ ;
			}
			cmi_2 = 1 + cmi2_max ;
		}
	}
}

void Sparse_matrix::SpUnit_matmult_set_ij_set(const SpUnit_matrix & U1, const SpUnit_matrix & U2, int * v) {
	// v (requis) = un double* dimensionne a  (n = U1.n_) + 1
	int * U1_ij = U1.ij_ ; int U1_nnz = U1.nnz_ ; int * U2_ij = U2.ij_ ; int U2_nnz = U2.nnz_ ;
	int ntnz = SpUnit_matmult_set_nnz_(U1_ij, U2_ij, m_, U1.n_, n_, U1_nnz, U2_nnz, v) ;
	if (npnz_ < ntnz) {
		if ((ij_ == U1_ij) || (ij_ == U2_ij)) assert(false) ;
		Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ; delete []  v_ ;
		npnz_ = ntnz ;
		ij_ = new int[n_ + m_ + npnz_ + npnz_] ;	v_ = new double[npnz_ + npnz_] ;
		if (( v_ == NULL ) || (ij_ == NULL)) {
			_LogMessage("Sparse_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
			assert(false) ;
		}
	}
	ntnz_ = ntnz ;
	SpUnit_matmult_fill_ij_(ij_, U1_ij, U2_ij, m_, U1.n_, n_, npnz_, U1_nnz, U2_nnz, v) ;
	Sparse_matmult_set_ij_count_update(ij_, U1_ij, U2_ij) ;
}

int SpUnit_matmult_set_nnz_(int* ij1, int* ij2, const int &m, const int &n, const int &p, const int &npnz1, const int &npnz2, int * v) {
	int i, j, k, ntnz = 0 ;
	int * cmi_2 = ij2 + p + n ; int * rmj_1 ; int * rmj1_max, * cmi2_max ;
	int * ni_2 = ij2 - 1 ; int * nj_1 = ij1 + n - 1 ;
	for (j = 1 ; j <= p ; j ++) {
		for (k = 1 ; k <= n ; k++)	 v[k] = 0;
		cmi2_max = cmi_2 + ni_2[j] - 1 ;
		while (cmi_2 <= cmi2_max)	{ v[*cmi_2] = 1 ; cmi_2 ++ ; }
		rmj_1 = ij1 + m + n + npnz1 ;
		for (i = 1 ; i <= m ; i ++) {
			rmj1_max = rmj_1 + nj_1[i] - 1 ;
			while (rmj_1 <= rmj1_max) {
				if (v[*rmj_1]) { ntnz ++ ; break ; }
				rmj_1 ++ ;
			}
			rmj_1 = 1 + rmj1_max ;
		}
	}
	return ntnz ;
}

void SpUnit_matmult_fill_ij_(int* ij, int*ij1, int* ij2, const int &m, const int &n, const int &p, const int &npnz, const int &npnz1, const int &npnz2, int* v) {
	int i, j, k ;
	int * ni_2 = ij2 - 1 ; int * nj_1 = ij1 + n - 1 ;
	int * ni_ = ij - 1 ; int * nj_ = ij + p - 1 ; int * cmi_ = ij + m + p ; int * rmj_ = cmi_ + npnz ;
	int * cmi_2 = ij2 + p + n ; int * rmj_1 ; int * rmj1_max, * cmi2_max ;
	for (j = 1 ; j <= p ; j ++) {
		ni_ ++ ; *ni_ = 0 ;
		for (k = 1 ; k <= n ; k++)	 v[k] = 0;
		cmi2_max = cmi_2 + ni_2[j] - 1 ;
		while (cmi_2 <= cmi2_max)	{ v[*cmi_2] = 1 ; cmi_2 ++ ; }
		rmj_1 = ij1 + m + n + npnz1 ;
		for (i = 1 ; i <= m ; i ++) {
			rmj1_max = rmj_1 + nj_1[i] - 1 ;
			while (rmj_1 <= rmj1_max) {
				if (v[*rmj_1]) { (*ni_) ++ ; *cmi_ = i ; cmi_ ++ ; break ; }
				rmj_1 ++ ;
			}
			rmj_1 = 1 + rmj1_max ;
		}
	}
	rmj_1 = ij1 + m + n + npnz1 ;
	for (i = 1 ; i <= m ; i ++) {
		nj_ ++ ; *nj_ = 0 ;
		for (k = 1 ; k <= n ; k++)	 v[k] = 0;
		rmj1_max = rmj_1 + nj_1[i] - 1 ;
		while (rmj_1 <= rmj1_max)	{ v[*rmj_1] = 1 ; rmj_1 ++ ; }
		cmi_2 = ij2 + p + n ;
		for (j = 1 ; j <= p ; j ++) {
			cmi2_max = cmi_2 + ni_2[j] - 1 ;
			while (cmi_2 <= cmi2_max) {
				if (v[*cmi_2]) { (*nj_) ++ ; *rmj_ = j ; rmj_ ++ ; break ; }
				cmi_2 ++ ;
			}
			cmi_2 = 1 + cmi2_max ;
		}
	}
}

void Sparse_matmult_set_ij_count_update(int * ij_, int * ij_1, int * ij_2) {
	int k = 3 * Sparse_matmult_set_ij_count ;
	if (Sparse_matmult_set_ij_count == Sparse_matmult_set_ij_max) {
		Sparse_matmult_set_ij_max += 10 ;
		int ** temp = new int* [3 * Sparse_matmult_set_ij_max] ;
		for (int i = 0 ; i < k ; i ++)  temp[i] = Sparse_matmult_set_ij_[i] ;
		delete [] Sparse_matmult_set_ij_ ;
		Sparse_matmult_set_ij_ = temp ;
	}
	Sparse_matmult_set_ij_count ++ ;
	if (Sparse_matmult_set_ij_count > 100)
		{ _LogMessage("ATTENTION : Sparse_matmult_set_ij_count > 100 !!!") ; assert(false) ; }
	Sparse_matmult_set_ij_[k] = ij_ ; Sparse_matmult_set_ij_[k+1] = ij_1 ; Sparse_matmult_set_ij_[k+2] = ij_2 ;
}

void Sparse_matmult_set_ij_uncount(int * ij_, int * ij_1, int * ij_2) {
	int i, ii = 0 ;
	bool Sparse_combination_found = false ;
	for (i = 3 * Sparse_matmult_set_ij_count - 3 ; i >= 0 ; i -= 3) {
		if (Sparse_matmult_set_ij_[i] == ij_) {
			if (((Sparse_matmult_set_ij_[i+1] == ij_1) && (Sparse_matmult_set_ij_[i+2] == ij_2)) || ((Sparse_matmult_set_ij_[i+2] == ij_1) && (Sparse_matmult_set_ij_[i+1] == ij_2))) {
				Sparse_combination_found = true ; ii = i ;break ;
			}
		}
	}
	if (! Sparse_combination_found) assert(false) ;
	if (ii != 3 * Sparse_matmult_set_ij_count - 3) {
		Sparse_matmult_set_ij_[ii] = Sparse_matmult_set_ij_[3 * Sparse_matmult_set_ij_count - 3] ;
		Sparse_matmult_set_ij_[ii + 1] = Sparse_matmult_set_ij_[3 * Sparse_matmult_set_ij_count - 2] ;
		Sparse_matmult_set_ij_[ii + 2] = Sparse_matmult_set_ij_[3 * Sparse_matmult_set_ij_count - 1] ;
	}
	Sparse_matmult_set_ij_count -- ;
}

void Sparse_matrix::check_addsubdiag_set_ij_(const int & i1, const int & i2) {
	int i, j, k ;
	if (Sparse_addsubdiag_set_ij_count == 0) {
		Sparse_addsubdiag_set_ij_ = new int* [10] ;
		Sparse_addsubdiag_i_set_ij_ = new int [20] ;
	}
	bool Sparse_addsubdiag_set_ij_done = false ;
	for (i = 0 ; i < Sparse_addsubdiag_set_ij_count ; i ++) {
		if (Sparse_addsubdiag_set_ij_[i] == ij_) {
			if ((Sparse_addsubdiag_i_set_ij_[i+i] == i1) && (Sparse_addsubdiag_i_set_ij_[i+i+1] == i2)) {
				Sparse_addsubdiag_set_ij_done = true ; break ;
			}
		}
	}
	if (!(Sparse_addsubdiag_set_ij_done)) {
		Sparse_XXX_set_ij_uncount(ij_) ;
		int ntnz = ntnz_ ;
		int * nj_ = ij_ + n_ - 1 ;
		int * rmj_, * rmj_max ;
		int * rmj_min = ij_ + n_ + m_ + npnz_ ;
		for (i = 1 ; i < i1 ; i ++) rmj_min += nj_[i] ;
		for (i = i1 ; i <= i2 ; i ++) {
			rmj_max = rmj_min + nj_[i] ;
			for (rmj_ = rmj_min ; rmj_ < rmj_max ; rmj_ ++) {
				if ((*rmj_) == i) { ntnz -- ; break ; }
			}
			ntnz ++ ; rmj_min = rmj_max ;
		}
		if (ntnz != ntnz_) {
			npnz_ = ntnz ; // a revoir ?
			int * temp_ij = new int[n_ + m_ + npnz_ + npnz_] ;
			double * temp_v = new double[npnz_ + npnz_] ;
			if (( temp_v == NULL ) || (temp_ij == NULL)) {
				_LogMessage("Sparse_matrix constructor failed -- insufficient memory. Returned matrix is of size 0 x 0") ;
				assert(false) ;
			}
			int * temp_nj = temp_ij + n_ - 1 ;
			int * temp_rmj = temp_ij + m_ + n_ + npnz_ ;
			double * temp_rmv = temp_v + npnz_ ;
			double * rmv_ = v_ + npnz_ ;
			bool diag_found ;
			rmj_ = ij_ + m_ + n_ + npnz_ ;
			j = (*rmj_) ;
			for (i = 1 ; i < i1 ; i ++) {
				rmj_max = rmj_ + nj_[i] - 1 ;
				while (rmj_ <= rmj_max)
					{*temp_rmj = j ; temp_rmj ++ ; rmj_ ++ ; j = *rmj_ ; *temp_rmv = *rmv_ ; rmv_ ++ ; temp_rmv ++ ;}
			}
			for (i = i1 ; i <= i2 ; i ++) {
				diag_found = false ;
				rmj_max = rmj_ + nj_[i] - 1 ;
				while ((rmj_ <= rmj_max) && (j <= i)) {
					*temp_rmj = j ; *temp_rmv = *rmv_ ;
					if (j == i)  diag_found = true ;
					temp_rmj ++ ; rmj_ ++ ; temp_rmv ++ ; rmv_ ++ ; j = (*rmj_) ;
				}
				if (diag_found)  temp_nj[i] = nj_[i] ;
				else  { temp_nj[i] = nj_[i] + 1 ; *temp_rmj = i ; temp_rmj ++ ; *temp_rmv = 0 ; temp_rmv ++ ; }
				while (rmj_ <= rmj_max)
					{*temp_rmj = j ; temp_rmj ++ ; rmj_ ++ ; j = *rmj_ ; *temp_rmv = *rmv_ ; rmv_ ++ ; temp_rmv ++ ;}
			}
			for (i = i2 + 1 ; i <= m_ ; i ++) {
				rmj_max = rmj_ + nj_[i] - 1 ;
				while (rmj_ <= rmj_max)
					{*temp_rmj = j ; temp_rmj ++ ; rmj_ ++ ; j = *rmj_ ; *temp_rmv = *rmv_ ; rmv_ ++ ; temp_rmv ++ ;}
			}
			int * ni_ = ij_ - 1 ; int * cmi_max ;
			int * cmi_ = ij_ + m_ + n_ ;
			int * temp_ni = temp_ij - 1 ;
			int * temp_cmi = temp_ij + m_ + n_ ;
			double * temp_cmv = temp_v ;
			double * cmv_ = v_ ;
			i = (*cmi_) ;
			for (j = 1 ; j < i1 ; j ++) {
				cmi_max = cmi_ + ni_[j] - 1 ;
				while (cmi_ <= cmi_max)
					{*temp_cmi = i ; temp_cmi ++ ; cmi_ ++ ; i = *cmi_ ; *temp_cmv = *cmv_ ; cmv_ ++ ; temp_cmv ++ ;}
			}
			for (j = i1 ; j <= i2 ; j ++) {
				diag_found = false ;
				cmi_max = cmi_ + ni_[j] - 1 ;
				while ((cmi_ <= cmi_max) && (i <= j)) {
					*temp_cmi = i ; *temp_cmv = *cmv_ ;
					if (i == j)  diag_found = true ;
					temp_cmi ++ ; cmi_ ++ ; temp_cmv ++ ; cmv_ ++ ; i = (*cmi_) ;
				}
				if (diag_found)  temp_ni[j] = ni_[j] ;
				else  { temp_ni[j] = ni_[j] + 1 ; *temp_cmi = j ; temp_cmi ++ ; *temp_cmv = 0 ; temp_cmv ++ ; }
				while (cmi_ <= cmi_max)
					{*temp_cmi = i ; temp_cmi ++ ; cmi_ ++ ; i = *cmi_ ; *temp_cmv = *cmv_ ; cmv_ ++ ; temp_cmv ++ ;}
			}
			for (j = i2 + 1 ; j <= n_ ; j ++) {
				cmi_max = cmi_ + ni_[j] - 1 ;
				while (cmi_ <= cmi_max)
					{*temp_cmi = i ; temp_cmi ++ ; cmi_ ++ ; i = *cmi_ ; *temp_cmv = *cmv_ ; cmv_ ++ ; temp_cmv ++ ;}
			}
			delete [] v_ ; Sparse_XXX_set_ij_uncount(ij_) ; delete [] ij_ ;
			v_ = temp_v ; ij_ = temp_ij ; ntnz_ = ntnz ;
		}
		if (Sparse_addsubdiag_set_ij_count == Sparse_addsubdiag_set_ij_max) {
			Sparse_addsubdiag_set_ij_max += 10 ;
			int ** temp = new int* [Sparse_addsubdiag_set_ij_max] ;
			int * temp_i = new int [2 * Sparse_addsubdiag_set_ij_max] ;
			for (i = 0 ; i < Sparse_addsubdiag_set_ij_count ; i ++)  temp[i] = Sparse_addsubdiag_set_ij_[i] ;
			for (i = 0 ; i < 2 ; i ++)  temp_i[i] = Sparse_addsubdiag_i_set_ij_[i] ;
			delete [] Sparse_addsubdiag_set_ij_ ; delete [] Sparse_addsubdiag_i_set_ij_ ;
			Sparse_addsubdiag_set_ij_ = temp ;	Sparse_addsubdiag_i_set_ij_ = temp_i ;
		}
		Sparse_addsubdiag_set_ij_[Sparse_addsubdiag_set_ij_count] = ij_ ;
		k =  2 * Sparse_addsubdiag_set_ij_count ;
		Sparse_addsubdiag_i_set_ij_[k] = i1 ; Sparse_addsubdiag_i_set_ij_[k + 1] = i2 ;
		Sparse_addsubdiag_set_ij_count ++ ;
		if (Sparse_addsubdiag_set_ij_count > 100)
			{ _LogMessage("ATTENTION : Sparse_addsubdiag_set_ij_count > 100 !!!") ; assert(false) ; }
	}
}
