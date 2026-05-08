// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef MYMATH_H_
#define MYMATH_H_

/**
 * My own minimalistic non-generic fixed dimension
 * not operator overloading vector matrix classes
 */

#include <cmath>
#include <deque>
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
	 * linearly interpolates a single value @param x, extrapolates with a constant value
	 * not in the Python binding, use scipy instead
	 */
	static double interp1(double x, std::vector<double> x_,std::vector<double> y_) {
	    if (x > x_.back()) { // check bounds
			return y_.back();
		}
		if (x < x_[0]) {
			return y_[0];
		}
		double lookUpIndex = std::distance(x_.begin(), std::lower_bound(x_.begin(), x_.end(), x)); // if we are within bounds find the index of the lower bound
		if (lookUpIndex == 0) {
			return y_[0];
		}
		double ipLinear = (x - x_[lookUpIndex-1])/(x_[lookUpIndex] - x_[lookUpIndex-1]);
		return (1.-ipLinear)*y_[lookUpIndex-1] + (ipLinear)*y_[lookUpIndex];
	}

};

/**
 * 3-D turtle graphics.
 *
 * Maintains a local coordinate frame (heading H, left L, up U) stored as
 * the columns of a Matrix3d, and a current position.  All rotations are
 * intrinsic (applied in the turtle's own frame).
 *
 * Convention:
 *   column 0 (H) – heading  (forward direction)
 *   column 1 (L) – left
 *   column 2 (U) – up
 *
 * Usage:
 * @code
 *   Turtle3D t;
 *   t.forward(5.0);
 *   t.pitch(M_PI / 6);   // nose up 30°
 *   t.forward(3.0);
 * @endcode
 */
class Turtle3D
{
public:

    /// Initialises at the origin with H = +x, L = +y, U = +z.
    Turtle3D()
        : pos(0., 0., 0.),
          frame(Matrix3d(1,0,0, 0,1,0, 0,0,1)) { }

    /// Initialises with a given position and frame.
    Turtle3D(const Vector3d& position, const Matrix3d& frame)
        : pos(position), frame(frame) { }

    /// Moves forward by @p dist along the current heading (column 0).
    void forward(double dist) {
        Vector3d h = heading();
        pos = pos.plus(h.times(dist));
    }

    /**
     * Yaw: rotate around the up axis (U, column 2) by @p angle [rad].
     * Positive angle turns left when viewed from above.
     */
    void turnLeft(double angle) {
        rotateAroundColumn(2, angle);
    }

    /**
     * Yaw: rotate around the up axis (U, column 2) by @p angle [rad].
     * Positive angle turns right (equivalent to turnLeft(-angle)).
     */
    void turnRight(double angle) {
        rotateAroundColumn(2, -angle);
    }

    /**
     * Pitch: rotate around the left axis (L, column 1) by @p angle [rad].
     * Positive angle pitches the heading upward.
     */
    void pitchUp(double angle) {
        rotateAroundColumn(1, angle);
    }

    /**
     * Pitch: rotate around the left axis (L, column 1) by @p angle [rad].
     * Positive angle pitches the heading downward.
     */
    void pitchDown(double angle) {
        rotateAroundColumn(1, -angle);
    }

    /**
     * Roll: rotate around the heading axis (H, column 0) by @p angle [rad].
     * Positive angle rolls the left side upward.
     */
    void rollLeft(double angle) {
        rotateAroundColumn(0, angle);
    }

    /**
     * Roll: rotate around the heading axis (H, column 0) by @p angle [rad].
     * Positive angle rolls the left side downward.
     */
    void rollRight(double angle) {
        rotateAroundColumn(0, -angle);
    }

    Vector3d getPosition() const { return pos; }         ///< Current turtle position
    Vector3d heading()     const { return frame.column(0); } ///< Current heading direction (unit vector)
    Vector3d left()        const { return frame.column(1); } ///< Current left direction (unit vector)
    Vector3d up()          const { return frame.column(2); } ///< Current up direction (unit vector)
    Matrix3d getFrame()    const { return frame; }        ///< Full local coordinate frame

    /// Teleports to @p p without changing orientation.
    void setPosition(const Vector3d& p) { pos = p; }

    /// Replaces the entire coordinate frame (columns: H, L, U); columns should be orthonormal.
    void setFrame(const Matrix3d& f) { frame = f; }

    std::string toString() const {
        std::ostringstream s;
        s << "pos=" << pos.toString()
          << " H=" << heading().toString()
          << " L=" << left().toString()
          << " U=" << up().toString();
        return s.str();
    }

private:

    /**
     * Rodrigues rotation of the frame around its own column @p axis by @p angle [rad].
     *
     * All three columns are rotated so the frame stays orthonormal.
     */
    void rotateAroundColumn(int axis, double angle) {
        assert(axis >= 0 && axis < 3);
        const double ca = std::cos(angle);
        const double sa = std::sin(angle);
        Vector3d k = frame.column(axis); // unit rotation axis

        for (int c = 0; c < 3; c++) {
            if (c == axis) continue; // axis column is unchanged
            Vector3d v = frame.column(c);
            // Rodrigues: v' = v*cos(a) + (k x v)*sin(a) + k*(k.v)*(1-cos(a))
            // Since k and v are already orthonormal frame columns, k.v == 0,
            // so the formula simplifies to: v' = v*cos(a) + (k x v)*sin(a)
            Vector3d vr = v.times(ca).plus(k.cross(v).times(sa));
            // Write back the rotated column
            switch(c) {
            case 0: frame.r0.x = vr.x; frame.r1.x = vr.y; frame.r2.x = vr.z; break;
            case 1: frame.r0.y = vr.x; frame.r1.y = vr.y; frame.r2.y = vr.z; break;
            case 2: frame.r0.z = vr.x; frame.r1.z = vr.y; frame.r2.z = vr.z; break;
            }
        }
    }

    Vector3d pos;   ///< Current position
    Matrix3d frame; ///< Local frame: columns are H (heading), L (left), U (up)
};


/**
 * A polyline defined in relative turtle-graphics coordinates.
 *
 * Each segment is stored as a TurtleNode: three intrinsic rotations (yaw, pitch,
 * roll) applied to the local frame arriving at that node, followed by a forward
 * step of @c dist.  Node 0 corresponds to deque entry 0 (the initial node);
 * subsequent nodes extend the polyline.
 *
 * Cartesian coordinates are computed on demand by replaying the turtle commands
 * from the front of the deque.  The anchor is the turtle's starting position for
 * playback and is not part of the public node indexing.
 *
 * Nodes can be appended at either end:
 *  - @c addNodeBack()  – grow the tip 
 *  - @c addNodeFront() – insert a at base 
 */
class Meristem
{
public:

    /// One segment of the polyline in turtle-relative coordinates.
    struct TurtleNode {
        double yaw   = 0.; ///< Rotation around U (up) before the forward step [rad]
        double pitch = 0.; ///< Rotation around L (left) before the forward step [rad]
        double roll  = 0.; ///< Rotation around H (heading) before the forward step [rad]
        double dist  = 0.; ///< Forward distance of this segment [cm]
    };

    /// Creates an empty meristem anchored at the origin with the default frame (H=+x, L=+y, U=+z).
    Meristem() : anchor(Vector3d(0., 0., 0.)), anchorFrame(Matrix3d(1,0,0, 0,1,0, 0,0,1))
        { nodes.push_back({0., 0., 0., 0.}); }

    /// Creates an empty meristem with an explicit anchor position and frame.
    Meristem(const Vector3d& anchorPos, const Matrix3d& frame)
        : anchor(anchorPos), anchorFrame(frame)
        { nodes.push_back({0., 0., 0., 0.}); }

    /**
     * Appends a new node at the back (tip) of the polyline.
     *
     * @param dist   forward distance of the new segment [cm]
     * @param yaw    rotation around U before moving forward [rad]
     * @param pitch  rotation around L before moving forward [rad]
     * @param roll   rotation around H before moving forward [rad]
     */
    void addNodeBack(double dist, double yaw = 0., double pitch = 0., double roll = 0.) {
        nodes.push_back({yaw, pitch, roll, dist});
    }

    /**
     * Prepends a new node at the front (base) of the polyline.
     * All existing node indices shift by one.
     *
     * @param dist   forward distance of the new segment [cm]
     * @param yaw    rotation around U before moving forward [rad]
     * @param pitch  rotation around L before moving forward [rad]
     * @param roll   rotation around H before moving forward [rad]
     */
    void addNodeFront(double dist, double yaw = 0., double pitch = 0., double roll = 0.) {
        nodes.push_front({yaw, pitch, roll, dist});
        ++initialNodeIdx;
    }

    /// Returns the number of nodes (equals the deque size).
    int size() const { return static_cast<int>(nodes.size()); }

    /**
     * Returns the Cartesian position of node @p i.
     *
     * Node i is reached by replaying deque entries 0..i starting from the anchor.
     *
     * @param i  node index in [0, size()-1]
     */
    Vector3d getNode(int i) const {
        assert(i >= 0 && i < size());
        Turtle3D t(anchor, anchorFrame);
        int idx = 0;
        for (const auto& n : nodes) {
            t.turnLeft(n.yaw);
            t.pitchUp(n.pitch);
            t.rollLeft(n.roll);
            t.forward(n.dist);
            if (idx == i) break;
            ++idx;
        }
        return t.getPosition();
    }

    /// Returns the TurtleNode at deque index @p i.
    const TurtleNode& getTurtleNode(int i) const { return nodes.at(i); }

    /// Returns all nodes as a vector of Cartesian positions (one per deque entry).
    std::vector<Vector3d> getPolyline() const {
        std::vector<Vector3d> pts;
        pts.reserve(size());
        Turtle3D t(anchor, anchorFrame);
        for (const auto& n : nodes) {
            t.turnLeft(n.yaw);
            t.pitchUp(n.pitch);
            t.rollLeft(n.roll);
            t.forward(n.dist);
            pts.push_back(t.getPosition());
        }
        return pts;
    }

    Vector3d getAnchor()      const { return anchor; }      ///< Anchor position (turtle start; not a numbered node)
    Matrix3d getAnchorFrame() const { return anchorFrame; } ///< Coordinate frame at the anchor
    void setAnchor(const Vector3d& p)      { anchor = p; }
    void setAnchorFrame(const Matrix3d& f) { anchorFrame = f; }

    /// Returns the deque index of the initial node (0 at construction; incremented by addNodeFront()).
    int getInitialNodeIndex() const { return initialNodeIdx; }

    /// Direct read access to the raw turtle-node deque.
    const std::deque<TurtleNode>& getNodes() const { return nodes; }

    std::string toString() const {
        std::ostringstream s;
        s << "Meristem [" << size() << " nodes]\n";
        auto pts = getPolyline();
        for (int i = 0; i < static_cast<int>(pts.size()); ++i) {
            s << "  " << i << ": " << pts[i].toString() << "\n";
        }
        return s.str();
    }

private:

    Vector3d anchor;              ///< Cartesian position of node 0
    Matrix3d anchorFrame;         ///< Local frame at the anchor
    std::deque<TurtleNode> nodes; ///< Turtle commands from base to tip
    int initialNodeIdx = 0;       ///< Deque index of the initial {0,0,0,0} node; incremented by addNodeFront()
};


} // end namespace CPlantBox

#endif
