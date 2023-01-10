/*
* PiafMunch (v.2) -- Implementing and Solving the Munch model of phloem sap flow in higher plants
*
* Copyright (C) 2004-2019 INRA
*
* Author: A. Lacointe, UMR PIAF, Clermont-Ferrand, France
*
* File: PM_arrays.h
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

#ifndef PM_ARRAYS_H
#define PM_ARRAYS_H

#include <assert.h>
#include <iostream>
using namespace std;
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <locale>
#include <sstream>
#include <vector>

extern char message[], name_[] ;
void _LogMessage(const char* message)  ;
void Update_Output(bool unconditional = false) ;
int MsgBox(const char* message, const char* titre = "Message", int button0 = 1024, int button1 = 0, int button2 = 0) ;

// -- NOTATION DES MULTIPLICATIONS ENTRE DEUX ARRAYS (matrice x matrice, matrice x vecteur, vecteur x vecteur) ----------------
//
//    !!!     Le symbole '*' represente la multiplication ELEMENTWISE, COLONNE PAR COLONNE     !!!
//    La multiplication matricielle est notee  matmult( ), le produit scalaire est note  dotprod( ),
//    et la multiplication elementwise ligne par ligne (in-place ou non)  est notee  row_elemult( )
//
// ----------------------------------------------------------------------------------------------------------------------------

class Index_vector ;
class Fortran_vector ;
class Index_matrix ;
class Sparse_matrix ;
class SpUnit_matrix ;

string Upper(string &S) ;

/******** Arithmetique composite 'in-place' commune a (presque) toutes les classes, en plus de celles propres a chacune *****/
#define set_(i, j, a)			add_(i, j, a, 0)			// X(i,j) = a  (X: Sparse_matrix)
#define sub_(i, j, a)			add_(i, j, a, -1)			// X(i,j) -= a  (X: Sparse_matrix)
#define set_add(a, b)			add_add(a, b, 0)			// X = a + b  (X: Fortran_matrix ou Fortran_vector)
#define set_sub(a, b)			add_sub(a, b, 0)			// X = a - b  (X: Fortran_matrix ou Fortran_vector)
#define set_mult(a, b)			add_mult(a, b, 0)			// X = a * b  (X: Fortran_matrix ou Fortran_vector)
#define set_div(a, b)			add_div(a, b, 0)			// X = a / b  (X: Fortran_matrix ou Fortran_vector)
#define set_elemult(a, b)		add_elemult(a, b, 0)		// X = elemult(a,b) (X: Fortran_matrix ou Fortran_vector)
#define set_elediv(a, b)		add_elediv(a, b, 0)			// X = elediv(a,b) (X: Fortran_matrix ou Fortran_vector)
#define set_row_elemult(a, b)	add_row_elemult(a, b, 0)	// X = row.elemult(a,b) (X: Fortran_matrix ou F_vector)
#define set_row_elediv(a, b)	add_row_elediv(a, b, 0)		// X = row.elediv(a,b) (X: Fortran_matrix ou F_vector)
#define set_subvector(i1,i2,a)	add_subvector(i1,i2,a,0)	// X[i1:12] = a     (X: Fortran_vector)
#define set_subvectorx(iv, a)	add_subvectorx(iv, a, 0)	// X[indices specified by iv] = a     (X: Fortran_vector)
#define set_submatrix(i1,i2,j1,j2,a)  add_submatrix(i1,i2,j1,j2,a,0)	// X[i1:12, j1:j2] = a     (X: Fortran_matrix)
#define sub_add(a, b)			add_add(a, b, -1)			// X -= a + b  (X: Fortran_matrix ou Fortran_vector)
#define sub_sub(a, b)			add_sub(a, b, -1)			// X -= a - b  (X: Fortran_matrix ou Fortran_vector)
#define sub_mult(a, b)			add_mult(a, b, -1)			// X -= a * b  (X: Fortran_matrix ou Fortran_vector)
#define sub_div(a, b)			add_div(a, b, -1)			// X -= a / b  (X: Fortran_matrix ou Fortran_vector)
#define sub_elemult(a, b)		add_elemult(a, b, -1)		// X -= elemult(a,b) (X: Fortran_matrix ou Fortran_vector)
#define sub_elediv(a, b)		add_elediv(a, b, -1)		// X -= elediv(a,b) (X: Fortran_matrix ou Fortran_vector)
#define sub_row_elemult(a, b)	add_row_elemult(a, b, -1)	// X -= row.elemult(a,b) (X: Fortran_matrix ou F_vector)
#define sub_row_elediv(a, b)	add_row_elediv(a, b, -1)	// X -= row.elediv(a,b) (X: Fortran_matrix ou F_vector)
#define sub_diag(a)				add_diag(a, -1)				// (X: Fortran_matrix ou Sparse_matrix)
#define sub_subdiag(i1, i2, a)	add_subdiag(i1, i2, a, -1)	// (X: Fortran_matrix ou Sparse_matrix)
#define sub_subvector(i1,i2,a)	add_subvector(i1,i2,a,-1)	// X[i1:12] -= a     (X: Fortran_vector)
#define sub_subvectorx(iv, a)	add_subvectorx(iv, a, -1)	// X[indices specified by iv] -= a     (X: Fortran_vector)
#define sub_submatrix(i1,i2,j1,j2,a)  add_submatrix(i1,i2,j1,j2,a,-1)	// X[i1:12, j1:j2] -= a     (X: Fortran_matrix)
#define abs_(a) ( ((a) < (0.)) ? (-a) : (a) )

//#ifndef min
//	#define min(a,b) ((a) <= (b) ? (a) : (b))
//	#define max(a,b) ((a) >= (b) ? (a) : (b))
//#endif // min


class Fortran_matrix
{
  private:
	int m_;							// nombre de lignes
	int n_;							// nombre de colonnes
	double* v_;						// 'vecteur'[1 + m_*n_] contenant, ranges par colonnes successives, ...
									// ...les elements de la matrice dont le nombre est stocke en v_[0] = m_*n_
  public:
	Fortran_matrix(int m = 0);							// matrice carree d'ordre m, NON INITIALISEE
	Fortran_matrix(int m, int n);						// matrice a m lignes x n colonnes, NON INITIALISEE
	Fortran_matrix(int m, int n, const double &a);		// matrice a m lignes x n colonnes, initialisee a 'a'
	Fortran_matrix(const Sparse_matrix &S);				// cree une matrice pleine a partir d'une creuse existante
	Fortran_matrix(const SpUnit_matrix &U);				// cree une matrice pleine a partir d'une creuse a +- 1.
	Fortran_matrix & operator=(const Fortran_matrix &A);// operateur de recopie
	~Fortran_matrix();
	int nblin() const;									// retourne le nombre de lignes de la matrice
	int nbcol() const;									// retourne le nombre de colonnes de la matrice
	void display(int i1 = 1, int i2 = 0, int j1 = 1, int j2 = 0);// display M[i1:i2, j1:j2] (default: i2 = 0 -> = m_ ; j2 = 0 -> = n_

	inline const double& operator()(int i, int j) const { // retourne la valeur de l'element d'indices i,j
		if ((i < 1) || (i > m_) || (j < 1) || (j > n_)) assert(false) ;
		return v_[ i + (j-1)*m_];
	}
	inline double & operator()(int i, int j) {			// reaffecte l'element d'indices (i,j)
		if ((i < 1) || (i > m_) || (j < 1) || (j > n_)) assert(false) ;
		return v_[ i + (j-1)*m_];
	}
	/**** Fonctions externes declarees 'friend' de la classe Fortran_matrix ***********************************/
	friend class Fortran_vector;							// impose par un croisement entre les 2 classes / arithmetique composite
	friend class Sparse_matrix;
	friend class SpUnit_matrix;
	friend Fortran_matrix submatrix(const Fortran_matrix &A, const int &i1, const int &i2, const int &j1, const int &j2) ;
};
	/**** sous-matrices, fusion, transposition, diff, diagonales (extraction et affectation) ******************/
	Fortran_matrix submatrix(const Fortran_matrix &A, const int &i1, const int &i2, const int &j1, const int &j2) ;


class Fortran_vector
{
  private:
	double* v_;						// 'vecteur'[n+1]: n elements stockes en v[1]..v[n] , et v[0] = n ;

  public:
	Fortran_vector(int size = 0);							// 'size' elements, non initialises
	Fortran_vector(vector<double> vecIn, double factor = 1.);// transform from cppvec to fortran vec
	Fortran_vector(int size, const double &a);				// 'size' elements, initialises a la valeur a
	Fortran_vector(int size, const double *b, bool ignore_zeroindex = true); // v.v_ = recopie de b (donc independant)
    Fortran_vector(const Fortran_vector &v);				// cree une copie du vecteur v
	Fortran_vector & operator=(const Fortran_vector &v);	// operateur de recopie
	~Fortran_vector();
	vector<double> toCppVector();							//fortran to cpp vector
	int size() const;										// retourne la taille du vecteur
	void display(int i1 = 1, int i2 = 0, int n_per_line = 0);// display v[i1:i2] (default: i2 = 0 -> = v.n_ ; n_p_l = 0 -> = 5 (DOS) or 10 (GUI)
	void set(const double &a);								// affecte la valeur indiquee a tous les elements
	void set(const Fortran_vector &v);						// vecteur argument = v, qui doit etre de meme taille
	void append(const Fortran_vector & v) ;					// fusion avec v en queue, eq. a  X = bind(X, v)
	void sequentialFill(std::vector<double> vecd, int smallVal, int bigVal);
	void zero(const Fortran_vector &z) ; // annule toute composante [i] inferieure en val.abs. a z[i]

	inline const double& operator()(int i) const {			// retourne la valeur de l'element d'indice i
		if ((i < 1) || (i > (int)v_[0])) assert(false) ;
		return v_[i];
	}
	inline double & operator()(int i) {						// reaffecte l'element d'indice i
		if ((i < 1) || (i > (int)v_[0])) assert(false) ;
		return v_[i];
	}
	inline const double& operator[](int i) const {			// retourne la valeur de l'element d'indice i
		if ((i < 1) || (i > (int)v_[0])) assert(false) ;
		return v_[i];
	}
	inline double & operator[](int i) {						// reaffecte l'element d'indice i
		if ((i < 1) || (i > (int)v_[0])) assert(false) ;
		return v_[i];
	}
	Fortran_vector operator-() ;							// operateur unaire de changement de signe
	Fortran_vector operator+(const Fortran_vector &v2);		//  v + v2
	Fortran_vector operator+(const double &a) ;				//  v + Fortran_vector(v.size, a)
	Fortran_vector operator-(const Fortran_vector &v2);		//  v - v2
	Fortran_vector operator-(const double &a) ;				//  v - Fortran_vector(v.size, a)
	Fortran_vector operator*(const Fortran_vector &v2);		//  v * v2 (mult. elementwise)
	Sparse_matrix operator*(const Sparse_matrix &S);		// matrice V*S (v.size = S.m_)
	Sparse_matrix operator*(const SpUnit_matrix &U);		// matrice V*S (v.size = U.m_)
	Fortran_vector operator/(const double &a);				//  v/a (div. par un scalaire);
	Fortran_vector operator/(const Fortran_vector &v2);		//  v/v2 (division elementwise)

	Fortran_vector & operator+=(const Fortran_vector &v2);	// v += v2
	Fortran_vector & operator+=(const double &a);			// v += Fortran_vector(v.size, a)
	Fortran_vector & operator-=(const Fortran_vector &v2);	// v -= v2
	Fortran_vector & operator-=(const double &a);			// v -= Fortran_vector(v.size, a)
	Fortran_vector & operator*=(double a);					// v *= a :  mult. par un scalaire
	Fortran_vector & operator*=(const Fortran_vector &v2);	// v = v * v2 (mult. elementwise)
	Fortran_vector & operator/=(const double &a);			// v /= a : division par un scalaire
	Fortran_vector & operator/=(const Fortran_vector &v2);	// v = v/v2 (division elementwise)

	/******************************* Arithmetique vectorielle 'inplace' composite  ******************************************/
	//  Dans la serie 'void add_XXX()', le dernier argument 'ad = 1', optionnel et invisible a l'utilisateur final, oriente vers
	// 'set_xxx()' ou 'sub_xxx()', soit par une macro #define (cf. en-tete de ce fichier), soit par une fonction vraie (X_matmult())
	void add_subvector(const int &i1, const int &i2, const Fortran_vector &v1, int ad = 1); // v[i1, i2] += v1
	void add_subvectorx(const Index_vector &index, const Fortran_vector &v, int ad = 1); //...specifies par 'index'
	void add_add(const Fortran_vector &v1, const Fortran_vector &v2, int ad = 1);	// v += v1 + v2
	void add_add(const Fortran_vector &v1, const double &a, int ad = 1);			// v[i] += v1[i] + a
	void add_sub(const Fortran_vector &v1, const Fortran_vector &v2, int ad = 1);	// v += v1 - v2
	void add_sub(const Fortran_vector &v1, const double &a, int ad = 1);			// v[i] += v[i]1 - a
	void add_sub(const double &a, const Fortran_vector &v1, int ad = 1);			// v[i] += a - v1[i]
	void add_mult(const double &a, const Fortran_vector &v1, int ad = 1);			// v += a * v1
	void add_elemult(const Fortran_vector &v1, const Fortran_vector &v2, int ad = 1);// v[i] += v1[i] * v2[i]
	void add_div(const Fortran_vector &v1, const double &a, int ad = 1);			// v[i] += v1[i] / a
	void add_div(const double &a, const Fortran_vector &v1, int ad = 1);			// v[i] += a / v1[i]
	void add_elediv(const Fortran_vector &v1, const Fortran_vector &v2, int ad = 1);// v[i] += v1[i] / v2[i]
	void add_matmult(const Sparse_matrix &S, const Fortran_vector &v1, int ad = 1);	// v += S*v1 (mult.matricielle)
	void add_matmult(const SpUnit_matrix &U, const Fortran_vector &v1, int ad = 1);	// v += U*v1 (mult.matricielle)
	void set_matmult(const Sparse_matrix &S, const Fortran_vector &v1);				// v = S*v1 (mult.matricielle)
	void set_matmult(const SpUnit_matrix &U, const Fortran_vector &v1);				// v = U*v1 (mult.matricielle)
	void sub_matmult(const Sparse_matrix &S, const Fortran_vector &v1);				// v -= S*v1 (mult.matricielle)
	void sub_matmult(const SpUnit_matrix &U, const Fortran_vector &v1);				// v -= U*v1 (mult.matricielle)
	int set_KLU_solve(const Sparse_matrix &S, const Fortran_vector &y, int** ipiv_ptr = NULL, void** TM_ptr = NULL, bool full_check = false);// resoud   S*x = y    en retournant  x =  inv(S) * y

	/**** Fonctions externes declarees 'friend' de la classe Fortran_vector ***********************************/
	friend class Fortran_matrix;
	friend class Sparse_matrix;
	friend class SpUnit_matrix;
	friend double* InPlace_Array(Fortran_vector &v) ;
	friend Fortran_vector subvectorx(const Fortran_vector &V, const Index_vector &index);
	friend Fortran_vector subvector(const Fortran_vector &V, const int &i1, const int &i2);
	friend Sparse_matrix diag(const Fortran_vector & V);
	friend Fortran_vector operator-(const double &a, const Fortran_vector &v) ;
	friend Fortran_vector operator*(double a, const Fortran_vector &v);
	friend Fortran_vector operator/(const double &a, const Fortran_vector &v);
	friend Sparse_matrix row_elemult(const Sparse_matrix &S, const Fortran_vector &v);
	friend Sparse_matrix row_elediv(const Sparse_matrix &S, const Fortran_vector &v);
	friend Sparse_matrix row_elemult(const SpUnit_matrix &U, const Fortran_vector &v);
	friend Sparse_matrix row_elediv(const SpUnit_matrix &U, const Fortran_vector &v);
	friend Fortran_vector matmult(const Sparse_matrix &S, const Fortran_vector &v);
};
    inline double* InPlace_Array(Fortran_vector &v) {return v.v_;} // attention, fonction dangereuse a utiliser avec precaution !!

	/**** sous-vecteur, fusion de 2 vecteurs, diff, diagonale d'une matrice (extraction et affectation) *******/
	Fortran_vector subvectorx(const Fortran_vector &V, const Index_vector &index) ;	// extraction d'un sous-vecteur
	Fortran_vector subvector(const Fortran_vector &V, const int &i1, const int &i2);// id./ ind. entre i1 et i2

	/**** arithmetique vectorielle ordinaire ******************************************************************/
	inline Fortran_vector operator+(const double &a, Fortran_vector &v) {return (v + a) ;} // Fortran_vector(v.size, a) + v
	Fortran_vector operator-(const double &a, const Fortran_vector &v) ;		// Fortran_vector(v.size, a) - v
	Fortran_vector operator*(double a, const Fortran_vector &v);				// a * v (mult. par un scalaire)
	Fortran_vector operator/(const double &a, const Fortran_vector &v);			// vecteur [a / v(i)]

class Index_vector
{
  private:
	int* v_; // 'vecteur'[n+1]: n elements (indices) stockes en v[1]..v[n] , et v[0] = n ;

  public:
	Index_vector(int size = 0); // cree un tableau de 'size' valeurs d'indices, initialisees a 0
	Index_vector(int size, const int start, const int incr = 1); // sequence (start=debut, incr=increment)
    Index_vector(const Index_vector &iv);	// cree une copie d'un index_vector existant
	Index_vector & operator=(const Index_vector &iv); // operateur de recopie
	int & operator()(int i);				// affectation d'un element du tableau d'indices
	const int & operator()(int i) const ; // retourne un element du tableau d'indices
	int & operator[](int i);				// affectation d'un element du tableau d'indices
	const int & operator[](int i) const ; // retourne un element du tableau d'indices
	int size() const;				// nombre d'elements (valeurs d'indices) du vecteur considere
	void append(const Index_vector & iv) ; // ajoute iv.size() valeurs d'indices supplementaires
	~Index_vector();

	friend class Fortran_vector;
	friend class Fortran_matrix;
	friend class Sparse_matrix;
	friend Fortran_vector subvectorx(const Fortran_vector &V, const Index_vector &index);
	friend Index_vector bind(const Index_vector & iv1, const Index_vector & iv2) ;
	friend Index_vector subvectorx(const Index_vector &IV, const Index_vector &index) ;
};
	Index_vector subvectorx(const Index_vector &IV, const Index_vector &index) ; // copie partielle = un sous-ensemble des indices
	Index_vector bind(const Index_vector & iv1, const Index_vector & iv2); // fusion de 2 index_vectors

class Sparse_matrix
{
  private:
	int m_, n_, ntnz_ ,npnz_;  // resp.: nb de lignes, de colonnes et d'elements POTENTIELLEMENT non-nuls de la matrice , et d'elements REELLEMENT non-nuls ('number of true NZ').
	int * ij_ ;			// array a (n_ + m_ + 2 * npnz_) elements, successivement :
						//		ni_[j = 1..n_] et nj_[i = 1..m_] = nb d'elements reellement non nuls, resp. de la colonne #j et de la ligne #i  -- sigma(ni_[j-1]) = sigma(nj_[i-1]) = ntnz_
						//		cmi_[1..ntnz_]  = indices i des ntnz_ elements reellement non-nuls (col-major),   SUIVI DE (npnz-ntnz) ZeROS  =>  LONGUEUR TOTALE DU SUB-ARRAY cmi_ = npnz
						//		rmj_[1..ntnz_]  = indices j des ntnz_ elements reellement non-nuls (row-major),   SUIVI DE (npnz-ntnz) ZeROS  =>  LONGUEUR TOTALE DU SUB-ARRAY rmj_ = npnz
	double * v_ ;		// (2 * npnz elts) : cmv_[1..ntnz_] suivi de (npnz-ntnz) zeros, puis rmv_[1..ntnz_] suivi de (npnz-ntnz) zeros = valeurs des npnz elts 'non nuls', resp. col_maj. et row_maj.

  public:
	Sparse_matrix() ;							// constructeur par defaut (matrice nulle, 0 x 0)
	Sparse_matrix(const int &m, const int &n, const int &npnz = 0) ; // initialize les dimensions, y compris npnz_
	Sparse_matrix(const Sparse_matrix &S, const int &npnz = -1) ;		// constructeur par recopie ; si npnz = -1 (defaut), garde le S.npnz_ ; si = -2, restreint npnz_ a ntnz_ .
	Sparse_matrix(const SpUnit_matrix &U) ;		// remplace 'true' par 1. et 'false' par -1.
	Sparse_matrix & operator=(const Sparse_matrix &S) ;// operateur de recopie : garde le S.npnz_
	~Sparse_matrix();
	int nblin() const;									// retourne le nombre de lignes de la matrice
	int nbcol() const;									// retourne le nombre de colonnes de la matrice
	int ntnz() const ;									// retourne le nombre total d'elements non nuls
	int npnz() const ;									// retourne le nombre maximal d'elements non nuls
	void reset_npnz(int npnz = 0) ;	// si npnz <= 0, npnz_ = ntnz_ ; sinon, npnz_ doit etre >= ntnz_.
	void display(int i1 = 1, int i2 = 0, int j1 = 1, int j2 = 0) // display M[i1:i2, j1:j2] (default: i2 = 0 -> = m_ ; j2 = 0 -> = n_
		{ (Fortran_matrix(*this)).display(i1, i2, j1, j2) ; }
	void dump() ; // affiche tous les attributs et les 4 arrays.
	void set(const Sparse_matrix &S, bool check = true) ;	// S_ = S ; les 2 matrices doivent etre de memes dimensions et ntnz.  (n'impose pas l'egalite des npnz_) -- check : verifie nb et position des nz
	const double& operator()(int i, int j) const ;  // retourne la valeur de l'element d'indices i,j
	void add_(int i, int j, double a, int ad = 1) ;	// affecte (si ad =0 : macro 'set_'), ajoute (si (ad =1 : 'add_' s.s.) ou soustrait (si ad =-1 : macro 'sub_') la valeur a a l'element d'indices (i,j)
	Sparse_matrix operator-() ;							// operateur unaire de changement de signe
	Sparse_matrix operator+(Sparse_matrix &S) ;
	Sparse_matrix operator+(const SpUnit_matrix &U) ;
	Sparse_matrix operator-(const Sparse_matrix &S) ;
	Sparse_matrix operator-(const SpUnit_matrix &U) ;
	Sparse_matrix operator*(const Sparse_matrix &S) ;	// elementwise
	Sparse_matrix operator*(const SpUnit_matrix &U) ;	// elementwise
	Sparse_matrix operator/(const double &a) ;			// matrice [S(i,j) / a]
	Sparse_matrix operator/(const Fortran_vector &v);	// matrice [S(i,j) / v('i')] ; si v.size() = nbcol()

	Sparse_matrix & operator+=(Sparse_matrix &S) ;
	Sparse_matrix & operator+=(const SpUnit_matrix &U) ;
	Sparse_matrix & operator-=(const Sparse_matrix &S) ;
	Sparse_matrix & operator-=(const SpUnit_matrix &U) ;
	Sparse_matrix & operator*=(const Sparse_matrix &S) ;	// elementwise
	Sparse_matrix & operator*=(const SpUnit_matrix &U) ;	// elementwise
	Sparse_matrix & operator*=(const double &a);			// S *= a (mult. par un scalaire)
	Sparse_matrix & operator*=(const Fortran_vector &v);	// S = v * S elemtw.; v.size = S.nblin()
	Sparse_matrix & operator/=(const double &a);			// S = S / a (division par un scalaire)
	Sparse_matrix & operator/=(const Fortran_vector &v);	// S = S / v elemtw.; v.size = S.nblin()
	void row_elemult(const Fortran_vector &v);				// S = v*S, elemtw. PAR LIGNE/ v.size= S.nbcol()
	void row_elediv(const Fortran_vector &v);				// S = S/v, elemtw. PAR LIGNE/ v.size= S.nbcol()

	/********************************** Arithmetique matricielle 'inplace' composite  *********************************/
	//  Dans la serie 'void add_XXX()', le dernier argument 'ad = 1', optionnel et invisible a l'utilisateur final, oriente vers
	// 'set_xxx()' ou 'sub_xxx()', soit par une macro #define (cf. en-tete de ce fichier), soit par une fonction vraie (x_matmult())
	void add_add(const Sparse_matrix &S1, const Sparse_matrix &S2, int ad = 1, bool save_count = true) ; // seul (ad = 0) est accepte
	void add_add(const Sparse_matrix &S, const SpUnit_matrix &U, int ad = 1, bool save_count = true) ; // seul (ad = 0) est accepte
	void add_add(const SpUnit_matrix &U, const Sparse_matrix &S, int ad = 1, bool save_count = true) ; // seul (ad = 0) est accepte
	void add_add(const SpUnit_matrix &U1, const SpUnit_matrix &U2, int ad = 1, bool save_count = true) ; // seul (ad = 0) est accepte
	void add_sub(const Sparse_matrix &S1, const Sparse_matrix &S2, int ad = 1, bool save_count = true) ; // seul (ad = 0) est accepte
	void add_sub(const Sparse_matrix &S, const SpUnit_matrix &U, int ad = 1, bool save_count = true) ; // seul (ad = 0) est accepte
	void add_sub(const SpUnit_matrix &U, const Sparse_matrix &S, int ad = 1, bool save_count = true) ; // seul (ad = 0) est accepte
	void add_sub(const SpUnit_matrix &U1, const SpUnit_matrix &U2, int ad = 1, bool save_count = true) ; // seul (ad = 0) est accepte
	void add_elemult(const Sparse_matrix &S1, const Sparse_matrix &S2, int ad = 1, bool save_count = true); // seul (ad = 0) est accepte
	void add_elemult(const Sparse_matrix &S, const SpUnit_matrix &U, int ad = 1, bool save_count = true); // seul (ad = 0) est accepte
	void add_elemult(const SpUnit_matrix &U, const Sparse_matrix &S, int ad = 1, bool save_count = true); // seul (ad = 0) est accepte
	void add_elemult(const Fortran_vector &v, const Sparse_matrix &S, int ad = 1, bool save_count = true, bool transpose = false);// v*S: multiplication elementwise par colonne, par un vecteur de taille m = S.nblin(, int ad = 1)
	void add_elemult(const Fortran_vector &v, const SpUnit_matrix &U, int ad = 1, bool save_count = true, bool transpose = false);// v*U: multiplication elementwise par colonne, par un vecteur de taille m = S.nblin(, int ad = 1)
	void add_mult(const double &a, const Sparse_matrix &S, int ad = 1, bool save_count = true) ; // a * S : multiplication par un scalaire,
	void add_mult(const double &a, const SpUnit_matrix &U, int ad = 1, bool save_count = true) ; // a * U : multiplication par un scalaire,
	void add_elediv(const Sparse_matrix &S, const Fortran_vector &v, int ad=1, bool save_count=true, bool transpose=false);// S/v: division elementwise par colonne, par un vecteur de taille m = S.nblin(, int ad = 1)
	void add_elediv(const SpUnit_matrix &U, const Fortran_vector &v, int ad=1, bool save_count=true, bool transpose=false);// U/v: division elementwise par colonne, par un vecteur de taille m = S.nblin(, int ad = 1)
	void add_div(const Sparse_matrix &S, const double &a, int ad = 1, bool save_count = true);	//   S / a : division par un scalaire
	void add_div(const SpUnit_matrix &U, const double &a, int ad = 1, bool save_count = true);	//   U / a : division par un scalaire
	void add_row_elemult(const Sparse_matrix &S, const Fortran_vector &v, int ad = 1, bool save_count = true);// S*v: multiplication elementwise PAR LIGNE, par un vecteur de taille n = S.nbcol(, int ad = 1)
	void add_row_elemult(const SpUnit_matrix &U, const Fortran_vector &v, int ad = 1, bool save_count = true);// U*v: multiplication elementwise PAR LIGNE, par un vecteur de taille n = S.nbcol(, int ad = 1)
	void add_row_elediv(const Sparse_matrix &S, const Fortran_vector &v, int ad = 1, bool save_count = true);// S/v: division elementwise PAR LIGNE, par un vecteur de taille n = S.nbcol(, int ad = 1)
	void add_row_elediv(const SpUnit_matrix &U, const Fortran_vector &v, int ad = 1, bool save_count = true);// U/v: division elementwise PAR LIGNE, par un vecteur de taille n = S.nbcol()
	void set_diag(const double &a);
	void set_diag(const Fortran_vector &V);
	void set_subdiag(const int & i1, const int & i2, const double &a); // id.; S(i,i) = a si i1 =< i = < i2, 0 sinon.
	void set_subdiag(const int & i1, const int & i2, const Fortran_vector &V);// id. avec v[i] au lieu de a.
	void add_diag(const double &a, int ad = 1);
	void add_diag(const Fortran_vector &V, int ad = 1);
	void add_subdiag(const int & i1, const int & i2, const double &a, int ad = 1); // id.; S(i,i) = a si i1 =< i = < i2, 0 sinon.
	void add_subdiag(const int & i1, const int & i2, const Fortran_vector &V, int ad = 1);// id. avec v[i] au lieu de a.
	void check_addsubdiag_set_ij_(const int & i1, const int & i2) ;
	//  dans les fonctions suivantes, si temp_v est fourni il doit pointer sur un vecteur  de taille  n = S1.n_  (ou U.n_, S.n_, U1.n_) qui servira de tampon :
	void set_matmult(const Sparse_matrix &S1, const Sparse_matrix &S2, Fortran_vector * temp_v = NULL, bool save_count = true);// M = S1*S2...
	void set_matmult(const SpUnit_matrix &U, const Sparse_matrix &S, Fortran_vector * temp_v = NULL, bool save_count = true);// M = U*S...
	void set_matmult(const Sparse_matrix &S, const SpUnit_matrix &U, Fortran_vector * temp_v = NULL, bool save_count = true);// M = S*U...
	void set_matmult(const SpUnit_matrix &U1, const SpUnit_matrix &U2, Index_vector * temp_v = NULL, bool save_count = true);// M = U1*U2...
	void Sparse_add_set_ij_set(const Sparse_matrix & S1, const Sparse_matrix & S2) ; //
	void Sparse_elemult_set_ij_set(const Sparse_matrix & S1, const Sparse_matrix & S2) ; //
	void Sparse_matmult_set_ij_set(const Sparse_matrix & S1, const Sparse_matrix & S2, double* v) ; // v : tampon facultatif de taille
	void SpUnit_matmult_set_ij_set(const SpUnit_matrix & U1, const SpUnit_matrix & U2, int* v) ; //

	friend class Fortran_matrix;
	friend class Fortran_vector;
	friend class SpUnit_matrix;
	friend Sparse_matrix transpose(const Sparse_matrix & S, int npnz); // si npnz = -1 (defaut), garde le S.npnz_ ; si = -2, restreint npnz_ a ntnz_ .
	friend Sparse_matrix operator*(const double &a, const Sparse_matrix &S);
	friend Sparse_matrix operator*(const double &a, const SpUnit_matrix &S);
	friend Sparse_matrix row_elemult(const SpUnit_matrix &S, const Fortran_vector &v);
	friend Sparse_matrix row_elediv(const SpUnit_matrix &S, const Fortran_vector &v);
	friend Sparse_matrix row_elemult(const Sparse_matrix &S, const Fortran_vector &v);
	friend Sparse_matrix row_elediv(const Sparse_matrix &S, const Fortran_vector &v);
	friend Sparse_matrix matmult(const Sparse_matrix &S1, const Sparse_matrix &S2, Fortran_vector * temp_v);
	friend Sparse_matrix matmult(const SpUnit_matrix &U, const Sparse_matrix &S, Fortran_vector * temp_v);
	friend Sparse_matrix matmult(const Sparse_matrix &S, const SpUnit_matrix &U, Fortran_vector * temp_v);
	friend Fortran_vector matmult(const Sparse_matrix &S, const Fortran_vector &v);
	friend void Fortran_vector::set_matmult(const Sparse_matrix &S, const Fortran_vector &v1);
	friend void KLU_allocate(const Sparse_matrix &S, int*** ipiv_ptr_ptr, void*** TM_ptr_ptr); // cree les espaces de travail appropries pour la matrice S
	friend int set_KLU_solve_(double*x, int p, const Sparse_matrix &S, double*b, int** ipiv_ptr, void** TM_ptr, bool full_check) ; // internal use, not for user !! use above function instead
};
	Sparse_matrix diag(const double &a, const int m, int n = -1);		// diag.: a * Id(m,n); carree (m,m) si n = -1
	Sparse_matrix diag(const Fortran_vector & V);			// construit la matrice carree diagonale <- vect. V
	Sparse_matrix diag(const Fortran_vector & V, const int m, const int n);//id, rectang.; V.size() = min(m,n)
	Sparse_matrix transpose(const Sparse_matrix & S, int npnz = -1); // si npnz = -1 (defaut), restreint npnz_ a ntnz_ ; si = -2, garde le S.npnz_.
	Sparse_matrix operator*(const double &a, const Sparse_matrix &S) ;			// a * S (multipl. par un scalaire)
	Sparse_matrix row_elemult(const Sparse_matrix &S, const Fortran_vector &v);//elemtw. PAR LIGNE/ v.size= n_)
	Sparse_matrix row_elediv(const Sparse_matrix &S, const Fortran_vector &v);// S/v elemtw. PAR LIGNE/ v.size= n_)
	Sparse_matrix matmult(const Sparse_matrix &S1, const Sparse_matrix &S2, Fortran_vector * temp_v = NULL);// S1*S2 (mult.mat...)
	Fortran_vector matmult(const Sparse_matrix &S, const Fortran_vector &v);
	Fortran_vector KLU_solve(const Sparse_matrix &S, const Fortran_vector &y, int** ipiv_ptr = NULL, void** TM_ptr = NULL, bool full_check = false); // resoud   S*x = y    en retournant  x =  inv(S) * y
	void KLU_allocate(const Sparse_matrix &S, int*** ipiv_ptr_ptr, void*** TM_ptr_ptr) ; // cree les espaces de travail appropries pour la matrice S
	void KLU_free(int*** ipiv_ptr_ptr, void*** TM_ptr_ptr); // libere et anNULLe les espaces de travail

class SpUnit_matrix
{
  private:
	int m_, n_, nnz_ ;  // resp.: nb de lignes, de colonnes et d'elements non-nuls de la matrice
	int * ij_ ;			// array a (m_ + n_ + 4*nnz_) elements, successivement :
						//		ni_[j = 1..n_] et nj_[i = 1..m_] = nb d'elements non nuls, resp. de la colonne #j et de la ligne #i,
						//		cmi_[1..nnz_] et rmj_[1..nnz_] = indices des elements non-nuls, resp. i (col-major) et j (row-maj.);
						//		cmv_[1..nnz_] et rmv[1..nnz_] = valeurs des nnz elts, resp. col_maj. et row_maj, codes: +1 = true, -1 = false
  public:
	SpUnit_matrix() ;							// constructeur par defaut (matrice nulle, 0 x 0)
	SpUnit_matrix(const int &m, const int &n, const int &nnz = 0) ; // initialize les dimensions, y compris nnz_
	SpUnit_matrix(const Sparse_matrix &S) ;		// constructeur par recopie
	SpUnit_matrix(const SpUnit_matrix &U) ;		// remplace 'true' par 1. et 'false' par -1.
	SpUnit_matrix & operator=(const SpUnit_matrix &U);// operateur de recopie
	~SpUnit_matrix();
	int nblin() const;									// retourne le nombre de lignes de la matrice
	int nbcol() const;									// retourne le nombre de colonnes de la matrice
	int nbnz() const ;									// retourne le nombre total d'elements marques comme non nuls
	void display(int i1 = 1, int i2 = 0, int j1 = 1, int j2 = 0)// display M[i1:i2, j1:j2] (default: i2 = 0 -> = m_ ; j2 = 0 -> = n_
		{ (Fortran_matrix(*this)).display(i1, i2, j1, j2) ; }

	void set(const SpUnit_matrix &U);					// U_ = U ; LES NZ (vrais) DOIVENT ETRE SYNCHRO dans les 2 matrices
	SpUnit_matrix operator-() ;							// operateur unaire de changement de signe
	Sparse_matrix operator+(const Sparse_matrix &S) ;
	Sparse_matrix operator+(const SpUnit_matrix &U) ;
	Sparse_matrix operator-(const Sparse_matrix &S) ;
	Sparse_matrix operator-(const SpUnit_matrix &U) ;
	Sparse_matrix operator*(const Sparse_matrix &S) ;	// elementwise
	SpUnit_matrix operator*(const SpUnit_matrix &U) ;	// elementwise
	Sparse_matrix operator/(const double &a) ;			// matrice [U(i,j) / a]
	Sparse_matrix operator/(const Fortran_vector &v);	// matrice [U(i,j) / v('i')] ; si v.size() = nbcol()

	SpUnit_matrix & operator*=(const SpUnit_matrix &U) ;// elementwise ; LES NZ (vrais) DOIVENT ETRE SYNCHRO dans les 2 matrices

	void add_elemult(const SpUnit_matrix &U1, const SpUnit_matrix &U2, int ad = 1, bool save_count = true); // seul (ad = 0) est accepte
	void SpUnit_elemult_set_ij_set(const SpUnit_matrix &U1, const SpUnit_matrix &U2) ; //

	friend class Fortran_matrix;
	friend class Fortran_vector;
	friend class Sparse_matrix;
	friend SpUnit_matrix transpose(const SpUnit_matrix &U);
	friend Sparse_matrix operator*(const double &a, const SpUnit_matrix &U);
	friend Sparse_matrix row_elemult(const SpUnit_matrix &U, const Fortran_vector &v);
	friend Sparse_matrix row_elediv(const SpUnit_matrix &U, const Fortran_vector &v);
	friend Fortran_vector matmult(const SpUnit_matrix &U, const Fortran_vector &v);
	friend Sparse_matrix matmult(const SpUnit_matrix &U, const Sparse_matrix &S, Fortran_vector * temp_v);
	friend Sparse_matrix matmult(const Sparse_matrix &S, const SpUnit_matrix &U, Fortran_vector * temp_v);
	friend Sparse_matrix matmult(const SpUnit_matrix &U1, const SpUnit_matrix &U2, Index_vector * temp_v);
	friend void Fortran_vector::set_matmult(const SpUnit_matrix &U, const Fortran_vector &v1);
};
	SpUnit_matrix transpose(const SpUnit_matrix &U);
	Sparse_matrix operator*(const double &a, const SpUnit_matrix &U) ;			// a * U (multipl. par un scalaire)
	Sparse_matrix row_elemult(const SpUnit_matrix &U, const Fortran_vector &v);	//elemtw. PAR LIGNE/ v.size= n_)
	Sparse_matrix row_elediv(const SpUnit_matrix &U, const Fortran_vector &v);	// U/v elemtw. PAR LIGNE/ v.size= n_)
	Fortran_vector matmult(const SpUnit_matrix &U, const Fortran_vector &v);
	Sparse_matrix matmult(const SpUnit_matrix &U, const Sparse_matrix &S, Fortran_vector * temp_v = NULL);// U*S (mult.mat...)
	Sparse_matrix matmult(const Sparse_matrix &S, const SpUnit_matrix &U, Fortran_vector * temp_v = NULL);// S*U (mult.mat...)
	Sparse_matrix matmult(const SpUnit_matrix &U1, const SpUnit_matrix &U2, Index_vector * temp_v = NULL);// U1*U2...

#endif
