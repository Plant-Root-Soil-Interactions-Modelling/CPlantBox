// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef MYMATH_H_
#define MYMATH_H_

/**
 * My own minimalistic non-generic fixed dimension
 * not operator overloading vector matrix classes
 */

#include <cmath>
#include <sstream>
#include <assert.h>
#include <vector>
#include <functional>


namespace CPlantBox {

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/**
 * Vector2i stores two int values
 */
class Vector2i
{
public:

	Vector2i(): x(0), y(0) { } ///< Default constructor
	Vector2i(int x_, int y_): x(x_),y(y_) { } ///< Constructor passing two ints
	Vector2i(const Vector2i& v): x(v.x), y(v.y) { } ///< Copy constructor
	Vector2i(const std::vector<int>& xx): x(xx.at(0)), y(xx.at(1)) { } ///< for Python coupling

	std::string toString() const {
		std::ostringstream strs;
		strs << "( "<<x<<", "<<y<<" )";
		return strs.str();
	} ///< creates a string representing the two doubles

	int x; ///< number 1
	int y; ///< number 2

};



/**
 * Vector2d stores two double values
 */
class Vector2d
{
public:

	Vector2d(): x(0), y(0) { } ///< Default constructor
	Vector2d(double x_, double y_): x(x_),y(y_) { } ///< Constructor passing two doubles
	Vector2d(const Vector2d& v): x(v.x), y(v.y) { } ///< Copy Constructor
	Vector2d(const std::vector<double>& xx): x(xx.at(0)), y(xx.at(1)) { } ///< for Python coupling

	std::string toString() const {
		std::ostringstream strs;
		strs << "( "<<x<<", "<<y<<" )";
		return strs.str();
	} ///< creates a string representing the two doubles

	double x; ///< number 1
	double y; ///< number 2

};



/**
 * Vector3d stores three double values
 */
class Vector3d
{
public:

	Vector3d(): x(0),y(0),z(0) { } ///< Default constructor
	Vector3d(double x_, double y_, double z_): x(x_), y(y_), z(z_) { } ///< Constructor passing three doubles
	Vector3d(const Vector3d& v): x(v.x), y(v.y), z(v.z) { } ///< Copy Constructor
	Vector3d(const std::vector<double>& xx): x(xx.at(0)), y(xx.at(1)), z(xx.at(2)) { } ///< for Python coupling

	static Vector3d rotAB(double a, double b) { ///< first column of Rx(b)*Rz(a)
		double sa = sin(a);
		return Vector3d(cos(a), sa*cos(b), sa*sin(b) );
	};

	void normalize() { double l=length(); x/=l; y/=l; z/=l; } ///< normalizes the vector
	Vector3d normalized() const { const double l=length(); return Vector3d(x/l,y/l,z/l); } ///< returns a normalized copy of the vector

	double times(const Vector3d& v) const { return v.x*x+v.y*y+v.z*z; } ///< inner product
	double length() const { return sqrt(x*x+y*y+z*z); } ///< returns the Euclidian length

	Vector3d times(const double s) const { return Vector3d(s*x,s*y,s*z); } ///< returns the vector multiplied by a scalar value
	Vector3d plus(const Vector3d& v) const { return Vector3d(x+v.x,y+v.y,z+v.z); } ///< adds a vector and returns the result
	Vector3d minus(const Vector3d& v) const { return Vector3d(x-v.x,y-v.y,z-v.z); } ///< subtracts a vector and returns the result
	Vector3d cross(const Vector3d& v) const { return Vector3d(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); } ///< takes the cross product

	std::string toString() const {
		std::ostringstream strs;
		strs << "( "<<x<<", "<<y<<", "<<z<<" )";
		return strs.str();
	} ///< creates a string representing the three doubles

	double x; ///< double number 1
	double y; ///< double number 2
	double z; ///< double number 3

};

inline bool operator==(const Vector3d& lhs, const Vector3d& rhs){ return ((lhs.x==rhs.x) && (lhs.y==rhs.y) && (lhs.z==rhs.z)); } // needed for boost python indexing suite
inline bool operator!=(const Vector3d& lhs, const Vector3d& rhs){ return !(lhs == rhs); }



/**
 * 3x3 Matrix class, compatible with Vector3d for basic linear algebra
 * (i.e. exactly the operations needed for CPlantBox)
 */
class Matrix3d
{
public:

	Matrix3d(): r0(1,0,0), r1(0,1,0), r2(0,0,1) { } ///< Default constructor (identity matrix)
	Matrix3d(double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33):
		r0(m11,m12,m13), r1(m21,m22,m23), r2(m31,m32,m33) { } ///< Constructor passing nine doubles as 3x3 matrix entries
	Matrix3d(const Vector3d& c1, const Vector3d& c2, const Vector3d& c3) :
		r0(c1.x,c2.x,c3.x), r1(c1.y,c2.y,c3.y), r2(c1.z,c2.z,c3.z) { } ///< Constructs the matrix from three column vectors
	Matrix3d(const Matrix3d& m): r0(m.r0), r1(m.r1), r2(m.r2) { } ///< Copy constructor

	static Matrix3d rotX(double a) {
		double ca = cos(a);
		double sa = sin(a);
		return Matrix3d(1,0,0,0,ca,-sa,0,sa,ca);
	} ///< Creates a rotation matrix around the X-axis
	static Matrix3d rotY(double a) {
		double ca = cos(a);
		double sa = sin(a);
		return Matrix3d(ca,0,sa,0,1,0,-sa,0,ca);
	} ///< Creates a rotation matrix around the Y-axis
	static Matrix3d rotZ(double a) {
		double ca = cos(a);
		double sa = sin(a);
		return Matrix3d(ca,-sa,0,sa,ca,0,0,0,1);
	} ///< Creates a rotation matrix around the Z-axis
	static Matrix3d rotAB(double a, double b) { ///< Rx(b)*Rz(a)
		auto rxb = rotX(b);
		rxb.times(rotZ(a));
		return rxb;
	};

	/**
	 * Creates an orthonormal system (ONS) around the vector v
	 *
	 * Remark: This is not unique, in Rootbox the ONS is rotated around a random angle to make it truly random.
	 * This is likely to be not optimal, but I guess it is fast enough anyway.
	 *
	 * @param v           vector, that will be normalized and then be the first column of the resulting ONS
	 *
	 * \return            three orthonormal column vectors
	 */
	static Matrix3d ons(Vector3d& v) {
		Vector3d v2;
		Vector3d v3;
		if ((std::abs(v.x)>=std::abs(v.y)) && (std::abs(v.x)>=std::abs(v.z))) { // choose x and z
			v2 = Vector3d(-v.z, 0, v.x);
			v3 = v.cross(v2);
		} else if ((std::abs(v.y)>=std::abs(v.x)) && (std::abs(v.y)>=std::abs(v.z))) { // choose y and z
			v2 = Vector3d(0,-v.z, v.y);
			v3 = v.cross(v2);
		} else if ((std::abs(v.z)>=std::abs(v.x)) && (std::abs(v.z)>=std::abs(v.y))) { // choose x and z
			v2 = Vector3d(-v.z, 0, v.x);
			v3 = v.cross(v2);
		}
		v.normalize();
		v2.normalize();
		v3.normalize();
		return Matrix3d(v,v2,v3);
	} ///< Creates an orthonormal system (ONS) around the vector v

	double det() const {
		return  r0.x*(r1.y*r2.z-r2.y*r1.z)-r0.y*(r1.x*r2.z-r1.z*r2.x)+r0.z*(r1.x*r2.y-r1.y*r2.x);
	} ///< determinant of the matrix

	Matrix3d inverse() const {
		double d = det();
		double idet = 1. / d;
		Matrix3d A;
		A.r0.x = (r1.y*r2.z-r2.y*r1.z)*idet;
		A.r0.y = (r0.z*r2.y-r0.y*r2.z)*idet;
		A.r0.z = (r0.y*r1.z-r0.z*r1.y)*idet;
		A.r1.x = (r1.z*r2.x-r1.x*r2.z)*idet;
		A.r1.y = (r0.x*r2.z-r0.z*r2.x)*idet;
		A.r1.z = (r1.x*r0.z-r0.x*r1.z)*idet;
		A.r2.x = (r1.x*r2.y-r2.x*r1.y)*idet;
		A.r2.y = (r2.x*r0.y-r0.x*r2.y)*idet;
		A.r2.z = (r0.x*r1.y-r1.x*r0.y)*idet;
		return A;
	} ///< calculates the inverse of the matrix

	Vector3d column(int i) const {
		assert((i>=0) && (i<3));
		switch(i) {
		case 0: return Vector3d(r0.x,r1.x,r2.x);
		case 1: return Vector3d(r0.y,r1.y,r2.y);
		case 2: return Vector3d(r0.z,r1.z,r2.z);
		}
		throw 0; // just to not produce a warning
	} ///< returns the i-th column of the matrix (i=0..2)

	Vector3d row(int i) const {
		assert((i>=0) && (i<3));
		switch(i) {
		case 0: return r0;
		case 1: return r1;
		case 2: return r2;
		}
		throw 0; // just to not produce a warning
	} ///< returns the i-th row of the matrix (i=0..2)

	void times(const Matrix3d& m) {
		r0 = Vector3d(r0.times(m.column(0)), r0.times(m.column(1)), r0.times(m.column(2)) );
		r1 = Vector3d(r1.times(m.column(0)), r1.times(m.column(1)), r1.times(m.column(2)) );
		r2 = Vector3d(r2.times(m.column(0)), r2.times(m.column(1)), r2.times(m.column(2)) );
	} ///< Multiplies matrix m from right

	Vector3d times(const Vector3d& v) const {
		return Vector3d(r0.times(v), r1.times(v), r2.times(v));
	} ///<  Multiplies vector v from right

	std::string toString() const {
		std::ostringstream strs;
		strs << r0.toString() << "\n" << r1.toString() << "\n" << r2.toString();
		return strs.str();
	} ///< creates a string representing the 3x3 matrix

	Vector3d r0; ///< row 1
	Vector3d r1; ///< row 2
	Vector3d r2; ///< row 3

};



/**
 * usefull
 */
class Function {
public:

	/**
	 * trapezoidal rule  (needed in dumux-rosi schröder)
	 */
	static double quad(std::function<double(double)> f, double a, double b, int n) {
		double h = (b-a)/n;
		double s = h*(0.5*f(a) + 0.5*f(b));
		for (int i=1; i<n; i++) {
			s+=h*f(a+i*h);
		}
		return s;
	};

	/**
	 *  mimics numpy.linspace, from stack overflow
	 *  not in the binding, use numpy.linspace
	 */
	static std::vector<double> linspace(double start, double end, int num) {
		std::vector<double> linspaced;
		if (num == 0) {
			return linspaced;
		}
		if (num == 1) {
			linspaced.push_back(start);
			return linspaced;
		}
		double delta = (end - start) / (num - 1);
		for(int i=0; i < num-1; ++i) {
			linspaced.push_back(start + delta * i);
		}
		linspaced.push_back(end); // I want to ensure that start and end are exactly the same as the input
		return linspaced;
	}

	/**
	 * linearly interpolates a single value @param x
	 * not in the binding, use scipy instead
	 */
	static double interp1(double x, std::vector<double> x_,std::vector<double> y_) {
		if (x > x_.back()) { // check bounds
			return y_.back();
		}
		if (x < x_[0]) {
			return y_[0];
		}

		double lookUpIndex = std::distance(x_.begin(), std::lower_bound(x_.begin(), x_.end(), x));     // if we are within bounds find the index of the lower bound
		if (lookUpIndex == 0) {
			return y_[0];
		}
		double ipLinear = (x - x_[lookUpIndex-1])/(x_[lookUpIndex] - x_[lookUpIndex-1]);
		return (1.-ipLinear)*y_[lookUpIndex-1] + (ipLinear)*y_[lookUpIndex];
	}

};





} // end namespace CPlantBox

#endif
