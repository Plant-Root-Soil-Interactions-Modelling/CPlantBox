/*
* PiafMunch (v.2) -- Implementing and Solving the Munch model of phloem sap flow in higher plants
*
* Copyright (C) 2004-2019 INRA
*
* Author: A. Lacointe, UMR PIAF, Clermont-Ferrand, France
*
* File: PM_matrix.cpp
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

#include <string.h>
#ifdef QT_GUI_LIB
#include <QtGui>
#endif // QT_GUI_LIB

#include <PiafMunch/PM_arrays.h>


char message[10000], name_[300] ;

void Update_Output(bool unconditional) {}

void _LogMessage(const char* message) {
		cout << message << endl ;
}

int MsgBox(const char* message, const char* titre, int button0, int button1, int button2) {
	#ifdef QT_GUI_LIB
		 return QMessageBox::warning(0, QString(titre), QString(message),  button0,  button1, button2) ;
	#else
		cout << message << endl ;
		return 1024 ; // simule 'OK'
	#endif
}

Fortran_matrix::Fortran_matrix(int m) : m_(m), n_(m), v_(new double [1 + m*m]) {
	if (m < 0)	assert(false);
	if ( v_ == NULL ) {
		_LogMessage("Fortran_matrix constructor failed -- insufficient memory.\n Returned matrix is of size 0 x 0");
		assert(false);
	}
	v_[0] = m_ * m_;
}

Fortran_matrix::Fortran_matrix(int m, int n) : m_(m), n_(n), v_(new double [1 + m*n]) {
	if ((m < 0) || (n < 0))	assert(false);
	if ( v_ == NULL ) {
		_LogMessage("Fortran_matrix constructor failed -- insufficient memory.\n Returned matrix is of size 0 x 0");
		assert(false);
	}
	v_[0] = m_ * n_;
}

Fortran_matrix::Fortran_matrix(int m, int n, const double &a) : m_(m), n_(n), v_(new double [1 + m*n]) {
	if ((m < 0) || (n < 0))	assert(false);
	if ( v_ == NULL ) {
		_LogMessage("Fortran_matrix constructor failed -- insufficient memory.\n Returned matrix is of size 0 x 0");
		assert(false);
	}
	int mn = m * n ; double* temp = v_  ; *temp = mn ;
	for (int i = 0 ; i < mn ; i++) { temp ++ ; *temp = a ; }
}

Fortran_matrix::Fortran_matrix(const Sparse_matrix &S) : m_(S.m_), n_(S.n_), v_(NULL) {
	int mn = m_ * n_;
	v_ = new double [mn + 1];
	if ( v_ == NULL ) {
		_LogMessage("Fortran_matrix constructor failed -- insufficient memory.\n Returned matrix is of size 0 x 0");
		assert(false);
	}
	v_[0] = mn ; double * temp ;
	int * temp_ni = S.ij_ - 1 ; int * temp_cmi = temp_ni + m_ + n_ ;
	double * temp_cmv = S.v_ - 1 ;
	int j, k, ni, ks = 0 ;
	for (k = 1 ; k <= mn ; k++)	 v_[k] = 0.;
	for (j = 1 ; j <= n_ ; j++)	{
		temp = v_ + (j-1) * m_ ; ni = temp_ni[j] ;
		for (k = 1 ; k <= ni ; k++)	{ ks ++ ; temp[temp_cmi[ks]] = temp_cmv[ks] ; }
	}
}

Fortran_matrix::Fortran_matrix(const SpUnit_matrix &U) : m_(U.m_), n_(U.n_), v_(NULL) {
	int mn = m_ * n_;
	v_ = new double [mn + 1];
	if ( v_ == NULL ) {
		_LogMessage("Fortran_matrix constructor failed -- insufficient memory.\n Returned matrix is of size 0 x 0");
		assert(false);
	}
	v_[0] = mn ; int nnz = U.nnz_ ; double * temp ;
	int * temp_ni = U.ij_ - 1 ; int * temp_cmi = temp_ni + m_ + n_ ;
	int * temp_cmv = temp_cmi + nnz + nnz ;
	int j, k, ni, ks = 0 ;
	for (k = 1 ; k <= mn ; k++)	 v_[k] = 0.;
	for (j = 1 ; j <= n_ ; j++)	{
		temp = v_ + (j-1) * m_ ; ni = temp_ni[j] ;
		for (k = 1 ; k <= ni ; k++)	{
			ks ++ ;
			if (temp_cmv[ks]) temp[temp_cmi[ks]] = 1. ;
			else temp[temp_cmi[ks]] = -1. ;
		}
	}
}

Fortran_matrix & Fortran_matrix::operator=(const Fortran_matrix &A) {
	int mnp1 = 1 + A.m_ * A.n_  ; //, incx = 1;
	if ((m_ != A.m_) || (n_ != A.n_)) {
		m_ = A.m_ ; n_ = A.n_;
		delete [ ] v_;
		v_ = new double [mnp1];
		if ( v_ == NULL ) {
			_LogMessage("Fortran_matrix constructor failed -- insufficient memory.\n Returned matrix is of size 0 x 0");
			assert(false);
		}
	}
	double* temp(v_), *temp0(A.v_) ;
	for(int i = 0 ; i < mnp1 ; i ++) { *temp = *temp0 ; temp ++ ; temp0 ++ ; }
	return *this;
}

Fortran_matrix::~Fortran_matrix() {
	delete [ ] v_;
}

int Fortran_matrix::nblin() const { return m_; }

int Fortran_matrix::nbcol() const { return n_; }

void Fortran_matrix::display(int i1, int i2, int j1, int j2) {
	if (i2 == 0)  i2 = m_ ;
	if (j2 == 0)  j2 = n_ ;
	if ((i1 < 1) || (i2 > m_) || (j1 < 1) || (j2 > n_)) {
		assert(false) ;
	}
	int i, j, l = sprintf( message, "   "), ll = 2 ;
	for (j=j1 ; j <= j2 ; j++) {
		l = sprintf( message + ll, "             ") ; ll = 14 * (j-j1+1) - 1 ;
		l = sprintf( message + ll, "[%4d]", j) ; ll += l ;
	}
	_LogMessage(message) ;
	for (i = i1 ; i <= i2 ; i++) {
		l = sprintf( message, "[%4d]", i); ll = l ;
		for (j=j1 ; j <= j2 ; j++) {
			l = sprintf( message + ll, "     ") ; ll = 14 * (j-j1+1) - 7 ;
			l = sprintf( message + ll, "%#13g", v_[ i + (j-1)*m_]) ; ll += l ;
		}
		_LogMessage(message) ;
	}
}

Fortran_matrix submatrix(const Fortran_matrix &A, const int &i1, const int &i2, const int &j1, const int &j2) {
	// M = A[i1:i2, j1:j2]
	int i, j, Am = A.m_, An = A.n_, m = i2-i1+1, n = j2-j1+1;
	if ((i1 < 1) || (i1 > i2) || (i2 > Am)) assert(false);
	if ((j1 < 1) || (j1 > j2) || (j2 > An)) assert(false);
	Fortran_matrix T(m,n);
	double *temp, *tempA ;
	for (j = 1 ; j <= n ; j++) {
		tempA = A.v_ + Am * (j + j1 - 2) + i1 - 1 ;
		temp = T.v_ + m*(j-1);
		for (i = 1 ; i <= m ; i++)
			temp[i] = tempA[i] ;
	}
	return T;
}
