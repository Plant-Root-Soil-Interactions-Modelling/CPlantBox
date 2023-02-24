/*
* PiafMunch (v.2) -- Implementing and Solving the Munch model of phloem sap flow in higher plants
*
* Copyright (C) 2004-2019 INRA
*
* Author: A. Lacointe, UMR PIAF, Clermont-Ferrand, France
*
* File: PiafMunch2.h (Light PM_arrays.h + fix _MSVC2015)
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
#ifndef PIAFMUNCH2_H
#define PIAFMUNCH2_H

#include <assert.h>
#include <iostream>  // 'allows for 'system("pause")', which antiquated <iostream.h> did not...
#include <vector>
using namespace std;
#include <math.h>
#include <sstream>
#include <string>

/******** Arithmetique composite 'in-place' commune a (presque) toutes les classes, en plus de celles propores a chacune *****/
#define set_add(a, b)			add_add(a, b, 0)			// X = a + b
#define set_sub(a, b)			add_sub(a, b, 0)			// X = a - b
#define set_mult(a, b)			add_mult(a, b, 0)			// X = a * b
#define set_div(a, b)			add_div(a, b, 0)			// X = a / b
#define set_elemult(a, b)		add_elemult(a, b, 0)		// X = elemult(a,b)
#define set_elediv(a, b)		add_elediv(a, b, 0)			// X = elediv(a,b)
#define set_subvector(i1,i2,a)	add_subvector(i1,i2,a,0)	// X[i1:12] = a     (X: Fortran_vector)
#define set_subvectorx(iv, a)	add_subvectorx(iv, a, 0)	// X[indices specified by iv] = a     (X: Fortran_vector)
#define sub_add(a, b)			add_add(a, b, -1)			// X -= a + b
#define sub_sub(a, b)			add_sub(a, b, -1)			// X -= a - b
#define sub_mult(a, b)			add_mult(a, b, -1)			// X -= a * b
#define sub_div(a, b)			add_div(a, b, -1)			// X -= a / b
#define sub_elemult(a, b)		add_elemult(a, b, -1)		// X -= elemult(a,b)
#define sub_elediv(a, b)		add_elediv(a, b, -1)		// X -= elediv(a,b)
#define sub_row_elemult(a, b)	add_row_elemult(a, b, -1)	// X -= row.elemult(a,b)
#define sub_row_elediv(a, b)	add_row_elediv(a, b, -1)	// X -= row.elediv(a,b)
#define sub_subvector(i1,i2,a)	add_subvector(i1,i2,a,-1)	// X[i1:12] -= a     (X: Fortran_vector)
#define sub_subvectorx(iv, a)	add_subvectorx(iv, a, -1)	// X[indices specified by iv] -= a     (X: Fortran_vector)
#define abs_(a) ( ((a) < (0.)) ? (-a) : (a) )

#ifndef min
	#define min(a,b) ((a) <= (b) ? (a) : (b)) 
	#define max(a,b) ((a) >= (b) ? (a) : (b))
#endif // min

class Index_vector ;

class Fortran_vector
{
  private:
	double* v_;						// 'vecteur'[n+1]: n elements stockes en v[1]..v[n] , et v[0] = n ;

  public:
	Fortran_vector(int size = 0);							// 'size' elements, non initialises
	Fortran_vector(int size, const double &a);				// 'size' elements, initialises a la valeur a
    Fortran_vector(const Fortran_vector &v);				// cree une copie du vecteur v
	Fortran_vector & operator=(const Fortran_vector &v);	// operateur de recopie
	~Fortran_vector();
	void set(const double &a);								// affecte la valeur indiquee a tous les elements
	void set(const Fortran_vector &v);						// vecteur argument = v, qui doit etre de meme taille
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

	/**** Fonctions externes declarees 'friend' de la classe Fortran_vector ***********************************/
	friend Fortran_vector subvectorx(const Fortran_vector &V, const Index_vector &index);
	friend Fortran_vector subvector(const Fortran_vector &V, const int &i1, const int &i2);
	friend Fortran_vector operator-(const double &a, const Fortran_vector &v) ;
	friend Fortran_vector operator*(double a, const Fortran_vector &v);
	friend Fortran_vector operator/(const double &a, const Fortran_vector &v);
};
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
	void append(const Index_vector & iv) ; // ajoute iv.size() valeurs d'indices supplementaires
	~Index_vector();

	friend class Fortran_vector;
	friend Fortran_vector subvectorx(const Fortran_vector &V, const Index_vector &index);
	friend Index_vector bind(const Index_vector & iv1, const Index_vector & iv2) ;
	friend Index_vector subvectorx(const Index_vector &IV, const Index_vector &index) ;
};
//	ostream & operator<<(ostream & sortie, const Index_vector &iv);
	Index_vector subvectorx(const Index_vector &IV, const Index_vector &index) ; // copie partielle = un sous-ensemble des indices
	Index_vector bind(const Index_vector & iv1, const Index_vector & iv2); // fusion de 2 index_vectors

#ifdef _MSVC2015 // fix 'unresolved external symbol ___iob_func' link error with older MSVC compiled libs
	FILE _iob[] = { *stdin, *stdout, *stderr };
	extern "C" FILE * __cdecl __iob_func(void) { return _iob; }
#endif

#endif