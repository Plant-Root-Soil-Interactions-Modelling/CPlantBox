/*
* PiafMunch (v.2) -- Implementing and Solving the Munch model of phloem sap flow in higher plants
*
* Copyright (C) 2004-2019 INRA
*
* Author: A. Lacointe, UMR PIAF, Clermont-Ferrand, France
*
* File: PM_vector.cpp
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


#include <PiafMunch/PM_arrays.h>

vector<double> Fortran_vector::toCppVector() {
	int n = (int)v_[0];
	vector<double> vecOut(n);
	for (int i = 1 ; i <= n ; i++){vecOut[i-1] = v_[i];}
	return vecOut;
}

Fortran_vector::Fortran_vector(vector<double> vecIn, double factor ) : v_(new double [1 + vecIn.size()]) {
	if (vecIn.size() < 0)	assert(false);
	if ( v_ == NULL ) {
		_LogMessage("Fortran_vector constructor failed -- insufficient memory.\n Returned vector is of size 0");
		assert(false);
	}
	v_[0] = vecIn.size();
	for (int i = 1 ; i <= vecIn.size() ; i++)  v_[i] = vecIn[i-1] * factor;
}

Fortran_vector::Fortran_vector(int size) : v_(new double [1 + size]) {
	if (size < 0)	assert(false);
	if ( v_ == NULL ) {
		_LogMessage("Fortran_vector constructor failed -- insufficient memory.\n Returned vector is of size 0");
		assert(false);
	}
	v_[0] = size;
}

Fortran_vector::Fortran_vector(int size, const double &a) : v_(new double [1 + size]) {
	if (size < 0)	assert(false);
	if ( v_ == NULL ) {
		_LogMessage("Fortran_vector constructor failed -- insufficient memory. Returned vector is of size 0");
		assert(false);
	}
	v_[0] = size;
	for (int i = 1 ; i <= size ; i++)  v_[i] = a;
}

Fortran_vector::Fortran_vector(int size, const double *b, bool ignore_zeroindex) : v_(new double [1 + size]) {
	if (size < 1)	assert(false);
	if ( v_ == NULL ) {
		_LogMessage("Fortran_vector constructor failed -- insufficient memory. Returned vector is of size 0");
		assert(false);
	}
	v_[0] = size;
	if (ignore_zeroindex)
		for (int i = 1 ; i <= size ; i++)  v_[i] = b[i] ;
	else {
		const double *c ;
		c = b-1 ;
		for (int i = 1 ; i <= size ; i++)  v_[i] = c[i] ;
	}
}

Fortran_vector::Fortran_vector(const Fortran_vector &v) : v_(NULL) {
	int i, np1 = 1 + (int)v.v_[0];//, incx = 1;
	v_ = new double [np1];
	if ( v_ == NULL ) {
		_LogMessage("Fortran_vector constructor failed -- insufficient memory. Returned vector is of size 0");
		assert(false);
	}
	double* temp(v_), *temp0(v.v_) ;
	for(i=0 ; i < np1 ; i++) {*temp = *temp0 ; temp ++ ; temp0 ++ ; }
}

Fortran_vector & Fortran_vector::operator=(const Fortran_vector &v) {
	int i, n = (int)v.v_[0] ;
	if (n != (int)v_[0]) {
		delete [] v_ ;
		v_ = new double [n + 1];
	}
	if ( v_ == NULL ) {
		_LogMessage("Fortran_vector constructor failed -- insufficient memory. Returned vector is of size 0");
		assert(false);
	}
	double* temp(v_), *temp0(v.v_) ;
	for(i=0 ; i <= n ; i++) { *temp = *temp0 ; temp ++ ; temp0 ++ ; }
	return *this;
}

Fortran_vector::~Fortran_vector() {
	delete [ ] v_;
}

int Fortran_vector::size() const { return (int)v_[0]; }

void Fortran_vector::set(const double &a) {
	int n = (int)v_[0];
	for (int i = 1 ; i <= n ; i++)
		v_[i] = a;
}

void Fortran_vector::set(const Fortran_vector &v) {
	// identique a l'operateur=(), mais impose l'egalite des tailles des 2 vecteurs arguments.
	int i, n = (int)v_[0] ;
	
	if (n != (int)v.v_[0])	{
		//std::cout<<"in set "<<n<<" "<<((int)v.v_[0])<<std::endl;
		for(i=0 ; i < n ; i++) {std::cout<<v_[i]<<" ";}
		//std::cout<<std::endl;
		for(i=0 ; i < ((int)v.v_[0]) ; i++) {std::cout<<v.v_[i]<<" ";}
		//std::cout<<std::endl;
		//assert(false);
	}
	double* temp(v_), *temp0(v.v_) ;
	for(i=0 ; i < n ; i++) { temp ++ ; temp0 ++ ; *temp = *temp0 ; }
}

void Fortran_vector::add_subvectorx(const Index_vector &index, const Fortran_vector &v, int ad) {
	// remplace (ou ajoute ou soustrait a, si ad != 0) les valeurs d'indices specifies par index, (par) celles du vecteur v de meme taille que index
	int i, *temp = index.v_ ;
	int n = temp[0], n_ = (int)v_[0];
	double *temp1 = v.v_ ;
	if ((int)temp1[0] != n) assert(false);
	if (ad == 0) {
		for (i = 0 ; i < n ; i++) { temp ++ ; temp1 ++ ; if ((*temp < 1) || (*temp > n_)) assert(false) ; else v_[*temp] = *temp1 ; }
	} else if (ad == 1) {
		for (i = 0 ; i < n ; i++) { temp ++ ; temp1 ++ ; if ((*temp < 1) || (*temp > n_)) assert(false) ; else v_[*temp] += *temp1 ; }
	} else { // (ad = -1)
		for (i = 0 ; i < n ; i++) { temp ++ ; temp1 ++ ; if ((*temp < 1) || (*temp > n_)) assert(false) ; else v_[*temp] -= *temp1 ; }
	}
}

void Fortran_vector::add_subvector(const int &i1, const int &i2, const Fortran_vector &v, int ad) {
	// remplace (ou additionne, ou soustrait a...) les elements d'indices i1 a i2 (par) les elements de v, dont la taille doit correspondre :
	int i, nx = i2 - i1 + 1 ; //, incx = 1;
	int n_ = (int)v_[0];
	if ((i1 < 1) || (i1 > i2) || (i2 > n_)) assert(false);
	if (nx != (int)v.v_[0]) assert(false);
	double * temp(v_ + i1), *temp0(v.v_ + 1) ;
	if (ad == 0)
	//	dcopy_(&nx, v.v_ + 1, &incx, v_ + i1, &incx);
		for(i=1 ; i <= nx ; i++) {*temp = *temp0 ; temp ++ ; temp0 ++ ; }
	else if (ad == 1)
		for(i=1 ; i <= nx ; i++) {*temp += *temp0 ; temp ++ ; temp0 ++ ; }
	else // (ad = -1)
		for(i=1 ; i <= nx ; i++) {*temp -= *temp0 ; temp ++ ; temp0 ++ ; }
}

void Fortran_vector::append(const Fortran_vector & v) {
	double * temp2 = v.v_; double* t_v = v_ ;
	int i, n1 = (int)v_[0] ; int n2 = (int)temp2[0];
	double * temp = new double [1 + n1 + n2];
	temp[0] = n1 + n2 ;
	for (i = 1 ; i <= n1 ; i++)
		temp[i] = v_[i];
	v_ = temp;
	temp = v_ + n1 ;
	for (i = 1 ; i <= n2 ; i++)
		temp[i] = temp2[i];
	delete [] t_v ;
}





void Fortran_vector::sequentialFill(std::vector<double> vecd, int smallVal, int bigVal) {
	int n1 = (int)v_[0] ; int n2 = (int)vecd.size();
	assert((smallVal <= bigVal)&&"Fortran_vector::sequentiallFill: smallVal > bigVal");
	assert(( smallVal <= n2 )&&"Fortran_vector::sequentiallFill: vecd.size() < smallVal");
	int numSeq = n1/bigVal;
	for(int z = 0; z< numSeq; z++)
	{
		for(int zz = 0; zz< smallVal; zz++)
		{
			v_[z*bigVal + zz + 1 ] = vecd[z*smallVal + zz];
		}
	}
	
}

void Fortran_vector::zero(const Fortran_vector& z) {
	double * temp = v_, *tempz = z.v_ ;
	int n = (int)(*temp) ;
	assert(n == (int)(*tempz)) ;
	for (int i = 0 ; i < n ; i++) {
		tempz ++ ; temp ++ ;
		if(abs_(*temp) < abs_(*tempz)) *temp = 0. ;
	}
}
//display(int i1 = 1, int i2 = size(), int n_per_line = 4)
void Fortran_vector::display(int i1, int i2, int n_per_line) {
	int n = size() ;
	if (i2 == 0)  i2 = n ;
	if ((i1 < 1) || (i2 > n)) {
		assert(false) ;
	}
	if (n_per_line == 0) {
			n_per_line = 4 ;
	}
	int i, j, k = (i2-i1+1)/n_per_line, l = 0, ll ;
	//std::cout<<"diplay "<<n<<" "<<i1<<" "<<i2<<" "<<k<<" "<<n_per_line<<std::endl;
	for (j = 1 ; j <= k ; j++) {
		//std::cout<<"		"<<j<<std::endl;
		l = sprintf(message, "[%4d:%4d]", i1 + (j-1) * n_per_line, i1 - 1 + j * n_per_line); ll = l ;
		for (i=1 ; i <= n_per_line ; i++) {
			l = sprintf( message + ll, "     ") ; ll = i * 14 - 1 ;
			l = sprintf( message + ll, "%#13g", v_[i + i1 - 1 + (j-1) * n_per_line]) ;ll += l ;
		}
		_LogMessage(message) ;
	}
	//std::cout<<"	"<<j<<" "<<(n - j*n_per_line)<<" "<<double(n - j*n_per_line)<<std::endl;
	if(n - k*n_per_line > 0){
		l = sprintf(message, "[%4d:%4d]", i1 + (j-1) * n_per_line, i2); ll = l ;
		for (i=1 ; i <= (i2-i1+1) % n_per_line ; i++) {
			l = sprintf( message + ll, "     ") ; ll = i * 14 - 1 ;
			l = sprintf( message + ll, "%#13g", v_[i + i1 - 1 + (j-1) * n_per_line]) ; ll += l ;
		}
		_LogMessage(message) ;
	}
}

/**** sous-vecteur, fusion de 2 vecteurs, diff, diagonale d'une matrice (extraction et affectation) *******/

Fortran_vector subvectorx(const Fortran_vector & V, const Index_vector &index) {
	int i, *temp ; double *tempV ;
	temp = index.v_ ; tempV = V.v_;
	int n = temp[0], n_ = (int)tempV[0];
	Fortran_vector T(n);
	double *temp1 = T.v_ ;
	for (i = 1 ; i <= n ; i++) {
		temp ++ ;
		if ((*temp < 1) || (*temp > n_)) assert(false);
		temp1[i] = tempV[*temp];
	}
	return T;
}

Fortran_vector subvector(const Fortran_vector & V, const int &i1, const int &i2) {
	int i, nx = i2 - i1 + 1 ;
	int n_ = (int)V.v_[0];
	if ((i1 < 1) || (i1 > i2) || (i2 > n_)) assert(false);
	Fortran_vector T(nx);
	double * temp(T.v_ + 1), *temp0(V.v_ + i1) ;
	for(i=1 ; i <= nx ; i++) {*temp = *temp0 ; temp ++ ; temp0 ++ ; }
	return T;
}

/**** arithmetique vectorielle ordinaire ******************************************************************/

Fortran_vector Fortran_vector::operator-() {
	//   - v (operateur unaire de changement de signe)
	int n = (int)v_[0];
	Fortran_vector v(n);
	double *temp = v.v_ ;
	for (int i = 1 ; i <= n ; i++)
		temp[i] = -v_[i];
	return v;
}

Fortran_vector Fortran_vector::operator+(const Fortran_vector &v2) {
	//   v + v2
	Fortran_vector v((int)v_[0]);
	v.add_add((*this), v2, 0);
	return v;
}

Fortran_vector Fortran_vector::operator+(const double &a) {
	//   v + Fortran_vector(v.size, a)
	Fortran_vector v((int)v_[0]);
	v.set_add((*this), a);
	return v;
}

Fortran_vector Fortran_vector::operator-(const Fortran_vector &v2) {
	//   v - v2
	Fortran_vector v((int)v_[0]);
	v.set_sub((*this), v2);
	return v;
}

Fortran_vector Fortran_vector::operator-(const double &a) {
	//   v - Fortran_vector(v.size, a)
	Fortran_vector v((int)v_[0]);
	v.set_sub((*this), a);
	return v;
}

Fortran_vector operator-(const double &a, const Fortran_vector &v) {
	//   Fortran_vector(v.size, a) - v
	Fortran_vector Tv((int)v.v_[0]);
	Tv.set_sub(a, v);
	return Tv;
}

Fortran_vector Fortran_vector::operator*(const Fortran_vector &v2) {
	//   v * v2 , multiplication elementwise (v et v2 doivent avoir la meme taille)
	Fortran_vector v((int)v_[0]);
	v.set_elemult((*this), v2);
	return v;
}

Fortran_vector operator*(double a, const Fortran_vector &v) {
	// a * v : multiplication par un scalaire
	Fortran_vector Tv((int)v.v_[0]);
	Tv.set_mult(a, v);
	return Tv;
}

Sparse_matrix Fortran_vector::operator*(const Sparse_matrix &S) {
	// v * S : multiplication elementwise par colonne, par un vecteur de taille m = S.nblin()
	Sparse_matrix M(S.m_, S.n_, S.ntnz_) ;
	M.set_elemult((*this), S);
	return M;
}

Sparse_matrix Fortran_vector::operator*(const SpUnit_matrix &U) {
	// v * U : multiplication elementwise par colonne, par un vecteur de taille m = U.nblin()
	Sparse_matrix M(U.m_, U.n_, U.nnz_) ;
	M.set_elemult((*this), U);
	return M;
}

Fortran_vector Fortran_vector::operator/(const Fortran_vector &v2) {
	//   v / v2 , division elementwise (v et v2 doivent avoir la meme taille)
	Fortran_vector v((int)v_[0]);
	v.set_elediv((*this), v2);
	return v;
}

Fortran_vector Fortran_vector::operator/(const double &a) {
	//   v/a : division par un scalaire
	Fortran_vector v((int)v_[0]);
	v.set_div((*this), a);
	return v;
}

Fortran_vector operator/(const double &a, const Fortran_vector &v) {
	//   Fortran_vector(v.size, a/v[i])
	Fortran_vector Tv((int)v.v_[0]);
	Tv.set_div(a, v);
	return Tv;
}

/**** arithmetique vectorielle 'in-place' *****************************************************************/

Fortran_vector & Fortran_vector::operator+=(const Fortran_vector &v2) {
	//  v += v2
	int i, n = (int)v_[0];
	double*temp(v_), *temp2 = v2.v_ ;
	if ((int)temp2[0] != n) assert(false);
	for (i = 0 ; i < n ; i++) { temp ++ ; temp2 ++ ; *temp += *temp2 ; }
	return *this;
}

Fortran_vector & Fortran_vector::operator+=(const double &a) {
	//  v += Fortran_vector(v.size, a)
	int i, n = (int)v_[0];
	double*temp(v_) ;
	for (i = 0 ; i < n ; i++) { temp ++ ; *temp += a ; }
	return *this;
}

Fortran_vector & Fortran_vector::operator-=(const Fortran_vector &v2) {
	//  v -= v2
	int i, n = (int)v_[0];
	double*temp(v_), *temp2 = v2.v_ ;
	if ((int)temp2[0] != n) assert(false);
	for (i = 0 ; i < n ; i++) { temp ++ ; temp2 ++ ; *temp -= *temp2 ; }
	return *this;
}

Fortran_vector & Fortran_vector::operator-=(const double &a) {
	//  v -= Fortran_vector(v.size, a)
	int i, n = (int)v_[0];
	double*temp(v_) ;
	for (i = 0 ; i < n ; i++) { temp ++ ; *temp -= a ; }
	return *this;
}

Fortran_vector & Fortran_vector::operator*=(const Fortran_vector &v2) {
	//   v = v * v2 , multiplication elementwise (v et v2 doivent avoir la meme taille)
	int i, n = (int)v_[0];
	double*temp(v_), *temp2 = v2.v_ ;
	if ((int)temp2[0] != n) assert(false);
	for (i = 0 ; i < n ; i++) { temp ++ ; temp2 ++ ; *temp *= *temp2 ; }
	return *this;
}

Fortran_vector & Fortran_vector::operator*=(double a) {
	int i, n = (int)v_[0];
	double*temp(v_) ;
	for (i = 0 ; i < n ; i++) { temp ++ ; *temp *= a ; }
	return *this;
}

Fortran_vector & Fortran_vector::operator/=(const Fortran_vector &v2) {
	//   v = v/v2 , division elementwise (v et v2 doivent avoir la meme taille)
	int i, n = (int)v_[0];
	double*temp(v_), *temp2 = v2.v_ ;
	if ((int)temp2[0] != n) assert(false);
	for (i = 0 ; i < n ; i++) { temp ++ ; temp2 ++ ; *temp /= *temp2 ; }
	return *this;
}

Fortran_vector & Fortran_vector::operator/=(const double &a) {
	int i, n = (int)v_[0];
	double*temp(v_) ;
	for (i = 0 ; i < n ; i++) { temp ++ ; *temp /= a ; }
	return *this;
}

/**** Arithmetique vectorielle 'inplace' composite (fonctions membres de la classe Fortran_vector) ************/

void Fortran_vector::set_matmult(const Sparse_matrix &S, const Fortran_vector &v1) {
	// v = S*v1 (mult.matricielle); S.m_= v.size ; S.n_= v1.size
	
	//std::cout<<"Fortran_vector::set_matmult "<<std::endl;
	(*this).add_matmult(S, v1, 0);
}

void Fortran_vector::set_matmult(const SpUnit_matrix &U, const Fortran_vector &v1) {
	// v = U*v1 (mult.matricielle); U.m_= v.size ; U.n_= v1.size
	(*this).add_matmult(U, v1, 0);
}

void Fortran_vector::add_add(const Fortran_vector &v1, const Fortran_vector &v2, int ad) {
	// v += v1 + v2  (ou :  v = v1 + v2 , ou   v -= v1 + v2  suivant la valeur de ad)
	int n = (int)v_[0];
	double *temp2 = v2.v_ ;
	double *temp1 = v1.v_ ;
	if ((temp1[0] != n) || (temp2[0] != n)) assert(false);
	if (ad == 0)
		for (int i = 1 ; i <= n ; i++) v_[i] = temp1[i] + temp2[i] ;
	else
		if (ad == 1)
			for (int i = 1 ; i <= n ; i++) v_[i] += temp1[i] + temp2[i] ;
		else {
			if (ad != -1) assert(false) ;
			for (int i = 1 ; i <= n ; i++) v_[i] -= temp1[i] + temp2[i] ;
		}
}

void Fortran_vector::add_add(const Fortran_vector &v1, const double &a, int ad) {
	// v[i] += v1[i] + a  (ou :  v[i] = v1[i] + a , ou   v[i] -= v1[i] + a  suivant la valeur de ad)
	int n = (int)v_[0];
	double *temp1 = v1.v_ ;
	if (temp1[0] != n) assert(false);
	if (ad == 0)
		for (int i = 1 ; i <= n ; i++) v_[i] = a + temp1[i] ;
	else
		if (ad == 1)
			for (int i = 1 ; i <= n ; i++) v_[i] += a + temp1[i] ;
		else {
			if (ad != -1) assert(false) ;
			for (int i = 1 ; i <= n ; i++) v_[i] -= a + temp1[i] ;
		}
}

void Fortran_vector::add_sub(const Fortran_vector &v1, const Fortran_vector &v2, int ad) {
	// v += v1 - v2  (ou :  v = v1 - v2 , ou   v -= v1 - v2  suivant la valeur de ad)
	int n = (int)v_[0];
	double *temp2 = v2.v_ ;
	double *temp1 = v1.v_ ;
	if ((temp1[0] != n) || (temp2[0] != n)) assert(false);
	if (ad == 0)
		for (int i = 1 ; i <= n ; i++) v_[i] = temp1[i] - temp2[i] ;
	else
		if (ad == 1)
			for (int i = 1 ; i <= n ; i++) v_[i] += temp1[i] - temp2[i] ;
		else {
			if (ad != -1) assert(false) ;
			for (int i = 1 ; i <= n ; i++) v_[i] -= temp1[i] - temp2[i] ;
		}
}

void Fortran_vector::add_sub(const Fortran_vector &v1, const double &a, int ad) {
	// v[i] += v1[i] + a  (ou :  v[i] = v1[i] + a , ou   v[i] -= v1[i] + a  suivant la valeur de ad)
	int n = (int)v_[0];
	double *temp1 = v1.v_ ;
	if (temp1[0] != n) assert(false);
	if (ad == 0)
		for (int i = 1 ; i <= n ; i++) v_[i] = temp1[i] - a ;
	else
		if (ad == 1)
			for (int i = 1 ; i <= n ; i++) v_[i] += temp1[i] - a ;
		else {
			if (ad != -1) assert(false) ;
			for (int i = 1 ; i <= n ; i++) v_[i] -= temp1[i] - a ;
		}
}

void Fortran_vector::add_sub(const double &a, const Fortran_vector &v1, int ad) {
	// v[i] += a - v1[i]  (ou :  v[i] = a - v1[i] , ou   v[i] -= a - v1[i]  suivant la valeur de ad)
	int n = (int)v_[0];
	double *temp1 = v1.v_ ;
	if (temp1[0] != n) assert(false);
	if (ad == 0)
		for (int i = 1 ; i <= n ; i++) v_[i] = a - temp1[i] ;
	else
		if (ad == 1)
			for (int i = 1 ; i <= n ; i++) v_[i] += a - temp1[i] ;
		else {
			if (ad != -1) assert(false) ;
			for (int i = 1 ; i <= n ; i++) v_[i] -= a - temp1[i] ;
		}
}

void Fortran_vector::add_mult(const double &a, const Fortran_vector &v1, int ad) {
	// v += a * v1  (ou :  v = a * v1 , ou   v -= a * v1  suivant la valeur de ad)
	int n = (int)v_[0];
	double *temp1 = v1.v_ ;
	if (temp1[0] != n) assert(false);
	if (ad == 0)
		for (int i = 1 ; i <= n ; i++) v_[i] = a * temp1[i] ;
	else {
		if (ad == 1)
			for (int i = 1; i <= n; i++) v_[i] += a * temp1[i];
		else 	// (ad == -1)
			for (int i = 1; i <= n; i++) v_[i] -= a * temp1[i];
	}
}

void Fortran_vector::add_elemult(const Fortran_vector &v1, const Fortran_vector &v2, int ad) {
	// v[i] += v1[i] * v2[i]  (ou :  v[i] = v1[i] * v2[i] , ou   v[i] -= v1[i] * v2[i]  suivant la valeur de ad)
	int n = (int)v_[0];
	double *temp2 = v2.v_ ;
	double *temp1 = v1.v_ ;
	if ((temp1[0] != n) || (temp2[0] != n)) assert(false);
	if (ad == 0)
		for (int i = 1 ; i <= n ; i++) v_[i] = temp1[i] * temp2[i] ;
	else
		if (ad == 1)
			for (int i = 1 ; i <= n ; i++) v_[i] += temp1[i] * temp2[i] ;
		else {
			if (ad != -1) assert(false) ;
			for (int i = 1 ; i <= n ; i++) v_[i] -= temp1[i] * temp2[i] ;
		}
}

void Fortran_vector::add_div(const Fortran_vector &v1, const double &a, int ad) {
	// v[i] += v1[i] / a  (ou :  v[i] = v1[i] / a , ou   v[i] -= v1[i] / a  suivant la valeur de ad)
	int n = (int)v_[0];
	double *temp1 = v1.v_ ;
	if (temp1[0] != n) assert(false);
	if (ad == 0)
		for (int i = 1 ; i <= n ; i++) v_[i] = temp1[i] / a ;
	else
		if (ad == 1)
			for (int i = 1 ; i <= n ; i++) v_[i] += temp1[i] / a ;
		else {
			if (ad != -1) assert(false) ;
			for (int i = 1 ; i <= n ; i++) v_[i] -= temp1[i] / a ;
		}
}

void Fortran_vector::add_div(const double &a, const Fortran_vector &v1, int ad) {
	// v[i] += a / v1[i]  (ou :  v[i] = a / v1[i] , ou   v[i] -= a / v1[i]  suivant la valeur de ad)
	int n = (int)v_[0];
	double *temp1 = v1.v_ ;
	if (temp1[0] != n) assert(false);
	if (ad == 0)
		for (int i = 1 ; i <= n ; i++) v_[i] = a / temp1[i] ;
	else
		if (ad == 1)
			for (int i = 1 ; i <= n ; i++) v_[i] += a / temp1[i] ;
		else {
			if (ad != -1) assert(false) ;
			for (int i = 1 ; i <= n ; i++) v_[i] -= a / temp1[i] ;
		}
}

void Fortran_vector::add_elediv(const Fortran_vector &v1, const Fortran_vector &v2, int ad) {
	// v[i] += v1[i] / v2[i]  (ou :  v[i] = v1[i] / v2[i] , ou   v[i] -= v1[i] / v2[i]  suivant la valeur de ad)
	int n = (int)v_[0];
	double *temp2 = v2.v_ ;
	double *temp1 = v1.v_ ;
	if ((temp1[0] != n) || (temp2[0] != n)) assert(false);
	if (ad == 0)
		for (int i = 1 ; i <= n ; i++) v_[i] = temp1[i] / temp2[i] ;
	else
		if (ad == 1)
			for (int i = 1 ; i <= n ; i++) v_[i] += temp1[i] / temp2[i] ;
		else {
			if (ad != -1) assert(false) ;
			for (int i = 1 ; i <= n ; i++) v_[i] -= temp1[i] / temp2[i] ;
		}
}

void Fortran_vector::add_matmult(const Sparse_matrix &S, const Fortran_vector &v1, int ad) {
	//  v += S*v1 : multiplication par une matrice creuse (ou :  v = S*v1, ou  v -= S*v1, suivant la valeur de ad)
	//std::cout<<"Fortran_vector::add_matmult "<<ad<<std::endl;
	int m = S.m_, n = S.n_, npnz = S.npnz_ ;
	double * temp_v = v1.v_ ;
	if (((int)v_[0] != m) || ((int)temp_v[0]!= n)) assert(false);
	double * temp_y = v_ ;
	int * temp_nj = S.ij_ + n - 1 ;
	double * temp_s = S.v_ + npnz - 1 ;					// temp_s = rmv_
	int * temp_rmj = temp_nj + m + npnz ;
	int i, k, nj ;
	if ((ad == 0) || (ad == 1)) {
		if (ad == 0) {for (i = 1 ; i <= m ; i ++)  v_[i] = 0 ;}
		for (i = 1 ; i <= m ; i ++) {
			temp_y ++ ; nj = temp_nj[i] ;
			for (k = 1 ; k <= nj ; k ++) {
				temp_s ++ ; temp_rmj ++ ;
				//std::cout<<"before "<<(*temp_y)<<" "<<(*temp_s)<<" "<<(*temp_rmj)<<" "<<temp_v[*temp_rmj];
				//std::cout<<" "<<((*temp_s) * temp_v[*temp_rmj])<<" "<<((*temp_y) +(*temp_s) * temp_v[*temp_rmj]);
				(*temp_y) += (*temp_s) * temp_v[*temp_rmj] ;
				//std::cout<<" after "<<(*temp_y);
			}
		}
	}
	else {
		if (ad != -1) assert(false) ;
		for (i = 1 ; i <= m ; i ++) {
			temp_y ++ ; nj = temp_nj[i] ;
			for (k = 1 ; k <= nj ; k ++) {
				temp_s ++ ; temp_rmj ++ ;
				(*temp_y) -= (*temp_s) * temp_v[*temp_rmj] ;
			}
		}
	}
	//std::cout<<"done matmult "<<std::endl;
}

void Fortran_vector::add_matmult(const SpUnit_matrix &U, const Fortran_vector &v1, int ad) {
	//  v += U*v1 : multiplication par une matrice a +-1 (ou :  v = U*v1, ou  v -= U*v1, suivant la valeur de ad)
	int m = U.m_, n = U.n_, nnz = U.nnz_ ;
	double * temp_v = v1.v_ ;
	if (((int)v_[0] != m) || ((int)temp_v[0]!= n)) assert(false);
	double * temp_y = v_ ;
	int * temp_nj = U.ij_ + n - 1 ;
	int * temp_rmj = temp_nj + m + nnz ;
	int * temp_s = temp_rmj + nnz + nnz ;			// temp_s = rmv_
	int i, k, nj ;
	if ((ad == 0) || (ad == 1)) {
		if (ad == 0) {for (i = 1 ; i <= m ; i ++)  v_[i] = 0 ;}
		for (i = 1 ; i <= m ; i ++) {
			temp_y ++ ; nj = temp_nj[i] ;
			for (k = 1 ; k <= nj ; k ++) {
				temp_s ++ ; temp_rmj ++ ;
				if (*temp_s)  (*temp_y) += temp_v[*temp_rmj] ;
				else		  (*temp_y) -= temp_v[*temp_rmj] ;
			}
		}
	}
	else {
		if (ad != -1) assert(false) ;
		for (i = 1 ; i <= m ; i ++) {
			temp_y ++ ; nj = temp_nj[i] ;
			for (k = 1 ; k <= nj ; k ++) {
				temp_s ++ ; temp_rmj ++ ;
				if (*temp_s)  (*temp_y) -= temp_v[*temp_rmj] ;
				else		  (*temp_y) += temp_v[*temp_rmj] ;
			}
		}
	}
}

void Fortran_vector::sub_matmult(const Sparse_matrix &S, const Fortran_vector &v1) {
	// v -= S*v1 (mult.matricielle); S.m_= v.size ; S.n_= v1.size
	(*this).add_matmult(S, v1, -1);
}

void Fortran_vector::sub_matmult(const SpUnit_matrix &U, const Fortran_vector &v1) {
	// v -= U*v1 (mult.matricielle); U.m_= v.size ; U.n_= v1.size
	(*this).add_matmult(U, v1, -1);
}
