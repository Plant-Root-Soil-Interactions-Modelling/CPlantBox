/*
* PiafMunch (v.2) -- Implementing and Solving the Munch model of phloem sap flow in higher plants
*
* Copyright (C) 2004-2019 INRA
*
* Author: A. Lacointe, UMR PIAF, Clermont-Ferrand, France
*
* File: index_vector.cpp
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


#include <PiafMunch/PM_arrays.h>

#define PM_BOUNDS_CHECK // ne concerne que le fichier index_vector.cpp

Index_vector::Index_vector(int size) : v_(new int [1 + size]) {
	#ifdef PM_BOUNDS_CHECK
		assert(size >= 0);
	#endif
	v_[0] = size;
	for (int i = 1 ; i <= size ; i++)  v_[i] = 0 ;
}

Index_vector::Index_vector(int size, const int start, const int incr) : v_(new int [1 + size]) {
	#ifdef PM_BOUNDS_CHECK
		assert(size >= 0);
		assert(start >= 1);
		assert(incr >= 1);
	#endif
	v_[0] = size;
	int r = start ;
	for (int i = 1 ; i <= size ; i++) {
		v_[i] = r;
		r+= incr ;
	}
}

Index_vector::Index_vector(const Index_vector & iv) : v_(NULL) {
	int * temp;
	temp = iv.v_;
	int n = temp[0];
	v_ = new int [1 + n];
	for (int i = 0 ; i <= n ; i++)
		v_[i] = temp[i];
}

Index_vector & Index_vector::operator=(const Index_vector & iv) {
	int * temp;
	temp = iv.v_;
	int n = temp[0];
	delete [ ] v_;
	v_ = new int [1 + n];
	for (int i = 0 ; i <= n ; i++)
		v_[i] = temp[i];
	return *this;
}

int & Index_vector::operator()(int i) {
	#ifdef TNT_BOUNDS_CHECK
		assert(i >= 1);
		assert(i <= v_[0]);
	#endif
	return v_[i];
}

const int& Index_vector::operator()(int i) const {
	#ifdef TNT_BOUNDS_CHECK
		assert(i >= 1);
		assert(i <= v_[0]);
	#endif
	return v_[i];
}

int & Index_vector::operator[](int i) {
	#ifdef PM_BOUNDS_CHECK
		assert(i >= 1);
		assert(i <= v_[0]);
	#endif
	return v_[i];
}

const int& Index_vector::operator[](int i) const {
	#ifdef PM_BOUNDS_CHECK
		assert(i >= 1);
		assert(i <= v_[0]);
	#endif
	return v_[i];
}

Index_vector subvectorx(const Index_vector & IV, const Index_vector &index) {
	int i, *temp, *tempIV ;
	temp = index.v_ ; tempIV = IV.v_;
	int n = temp[0];
	#ifdef PM_BOUNDS_CHECK
		int n_ = tempIV[0];
		for (i = 1 ; i <= n ; i++) {
			assert(temp[i] >= 1);
			assert(temp[i] <= n_);
		}
	#endif
	Index_vector T(n);
	int *temp1 ;
	temp1 = T.v_;
	for (i = 1 ; i <= n ; i++)
		temp1[i] = tempIV[temp[i]];
	return T;
}

int Index_vector::size() const { return v_[0]; }

Index_vector::~Index_vector() {
	delete [ ] v_ ;
}

void Index_vector::append(const Index_vector & iv) {
	int * temp2 = iv.v_;
	int i, n1 = v_[0] ; int n2 = temp2[0];
	int * temp = new int [1 + n1 + n2];
	temp[0] = n1 + n2 ;
	for (i = 1 ; i <= n1 ; i++)
		temp[i] = v_[i];
	delete [] v_;
	v_ = temp;
	temp = v_ + n1 ;
	for (i = 1 ; i <= n2 ; i++)
		temp[i] = temp2[i];
}

Index_vector bind(const Index_vector & iv1, const Index_vector & iv2) {
	int *temp, *temp1, *temp2;
	temp1 = iv1.v_ ; temp2 = iv2.v_ ;
	int i, n1 = temp1[0] ; int n2 = temp2[0];
	Index_vector T(n1 + n2) ;
	temp = T.v_ ;
	for (i = 1 ; i <= n1 ; i++)
		temp[i] = temp1[i];
	temp = T.v_ + n1 ;
	for (i = 1 ; i <= n2 ; i++)
		temp[i] = temp2[i];
	return T ;
}
