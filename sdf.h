#ifndef SDF_H
#define SDF_H

#include <iostream>
#include <vector>
#include <stdexcept>

#include "mymath.h"



/**
 * Signed Distance Function (minus is inside, plus is outside)
 *
 * SignedDistanceFunction is the base class of all geometries (e.g. PlantContainer).
 * Furthermore, the class is used when no geometry is set, and describes an unconstrained setting.
 */
class SignedDistanceFunction
{

public:
    virtual ~SignedDistanceFunction() {};

    /**
     * Returns the signed distance to the next boundary
     *
     * @param v     spatial position [cm]
     * \return      signed distance [cm], a minus sign means inside, plus outside
     */
    virtual double getDist(const Vector3d& v) const { return -1e100; } ///< Returns the signed distance to the next boundary

    /**
     * Returns a string representation of the object (for debugging)
     */
    virtual std::string toString()  { return "SignedDistanceFunction"; }

    /**
     * Writes a ParaView Phython script explicitly representing the implicit geometry
     *
     * @param cout      e.g. a file output stream
     * @param c         python object counter for the script (to avoid duplicate names)
     * \return          object counter
     */
    virtual int writePVPScript(std::ostream & cout, int c=1) const { return c; };
    ///< Writes a ParaView Phython script explicitly representing the implicit geometry

    virtual std::string writePVPScript() const; ///< Writes the ParaView Phython script into a string

};



/**
 * PlantBox describes a rectangular box
 */
class SDF_PlantBox : public SignedDistanceFunction
{

public:

    /**
     * Creates a rectangular box: [-x/2,-y/2,0] - [x/2,y/2,-z]
     *
     * @param x  depth
     * @param y  width
     * @param z  height
     */
    SDF_PlantBox(double x, double y, double z) { dim = Vector3d(x/2.,y/2.,z/2.); } ///< creates a rectangular box

    virtual double getDist(const Vector3d& v) const; ///< @see SignedDistanceFunction::getDist

    virtual std::string toString()  { return "SDF_PlantBox"; } ///< @see SignedDistanceFunction::toString

    virtual int writePVPScript(std::ostream & cout, int c=1) const;  ///< @see SignedDistanceFunction::writePVPScript

private:
    Vector3d dim; // dimensions of the box
};



/**
 * Cylindrical or square container
 */
class SDF_PlantContainer : public SignedDistanceFunction
{

public:

    SDF_PlantContainer() { r1=5; r2=5; h=100; square = false; } ///< Default is a cylindrical rhizotron with radius 10 cm and 100 cm depth
    SDF_PlantContainer(double r1_, double r2_, double h_, double sq=false); ///< Creates a cylindrical or square container

    virtual double getDist(const Vector3d& v) const; ///< @see SignedDistanceFunction::getDist

    virtual std::string toString()  { return "SDF_PlantContainer"; } ///< @see SignedDistanceFunction::toString

    virtual int writePVPScript(std::ostream & cout, int c=1) const; ///< @see SignedDistanceFunction::writePVPScript

private:
    double r1;
    double r2;
    double h;
    bool square;
};



/**
 * SDF_RotateTranslate first rotates, and then translates a base geometry
 */
class SDF_RotateTranslate :public SignedDistanceFunction
{

public:

    enum SDF_Axes { xaxis=0, yaxis=1, zaxis=2 };

    SDF_RotateTranslate(SignedDistanceFunction* sdf, double angle=0, int axis=xaxis, const Vector3d& pos = Vector3d(0,0,0)); ///< Constructor
    SDF_RotateTranslate(SignedDistanceFunction* sdf, Vector3d pos): SDF_RotateTranslate(sdf, 0., xaxis, pos) { } ///< Translate only

    virtual double getDist(const Vector3d& v) const; ///< @see SignedDistanceFunction::getDist

    virtual std::string toString()  { return "SDF_RotateTranslate"; } ///< @see SignedDistanceFunction::toString

    virtual int writePVPScript(std::ostream & cout, int c=1)  const; ///< @see SignedDistanceFunction::writePVPScript

private:
    SignedDistanceFunction* sdf; // base geometry
    Vector3d pos; // translate origin to this position
    Matrix3d A; // rotation matrix
    int axis; // axis and angle (in degree) is needed for the python script
    double angle;
};



/**
 * SDF_Intersection computes the geometric intersection between several signed distance functions
 */
class SDF_Intersection : public SignedDistanceFunction
{

public:
    SDF_Intersection(std::vector<SignedDistanceFunction*> sdfs_) { sdfs=sdfs_; }
    ///< Constructs (sdfs_[0] ∩ sdfs_[1] ∩ sdfs_[2] ∩ ... )
    SDF_Intersection(SignedDistanceFunction* sdf1, SignedDistanceFunction* sdf2) { sdfs.push_back(sdf1); sdfs.push_back(sdf2); }
    ///< Constructs (sdf1 ∩ sdf2)

    virtual double getDist(const Vector3d& v) const;  ///< @see SignedDistanceFunction::getDist

    virtual std::string toString()  { return "SDF_Intersection"; } ///< @see SignedDistanceFunction::toString

    virtual int writePVPScript(std::ostream & cout, int c=1) const; ///< @see SignedDistanceFunction::writePVPScript

protected:
    std::vector<SignedDistanceFunction*> sdfs; ///< the set of signed distance functions
};



/**
 * SDF_Union computes the geometric union between several signed distance functions
 */
class SDF_Union : public SDF_Intersection
{

public:
    SDF_Union(std::vector<SignedDistanceFunction*> sdfs): SDF_Intersection(sdfs) { } ///< Constructs (sdfs_[0] U sdfs_[1] U sdfs_[2] U ... )
    SDF_Union(SignedDistanceFunction* sdf1, SignedDistanceFunction* sdf2): SDF_Intersection(sdf1,sdf2) { } ///< Constructs sdf1 U sdf2

    virtual double getDist(const Vector3d& v) const;  ///< @see SignedDistanceFunction::getDist

    virtual std::string toString()  { return "SDF_Union"; } ///< @see SignedDistanceFunction::toString
};



/**
 * SDF_Difference computes the difference between the first and several othter signed distance functions
 */
class SDF_Difference : public SDF_Intersection
{

public:
    SDF_Difference(std::vector<SignedDistanceFunction*> sdfs_) :SDF_Intersection(sdfs) { } ///< Constructs (...((sdfs_[0] \ sdfs_[1]) \ sdfs_[2])...)
    SDF_Difference(SignedDistanceFunction* sdf1, SignedDistanceFunction* sdf2) :SDF_Intersection(sdf1,sdf2) { } ///< Constructs sdf1 \ sdf2

    virtual double getDist(const Vector3d& v) const;  ///< @see SignedDistanceFunction::getDist

    virtual std::string toString()  { return "SDF_Difference"; } ///< @see SignedDistanceFunction::toString
};



/**
 * SDF_Complement computes the geometric complement (i.e. turns the geometry inside out)
 */
class SDF_Complement : public SignedDistanceFunction
{

public:
    SDF_Complement(SignedDistanceFunction* sdf_) { sdf=sdf_; } ///< Constructs the complement (sdf_)^c

    virtual double getDist(const Vector3d& v) const { return -sdf->getDist(v); } ///< @see SignedDistanceFunction::getDist

    virtual int writePVPScript(std::ostream & cout, int c=1) const { return sdf->writePVPScript(cout,c); } ///< same as original geometry

    virtual std::string toString()  { return "SDF_Complement"; } ///< @see SignedDistanceFunction::toString

private:
    SignedDistanceFunction* sdf;
};



/**
 * SDF_HalfPlane defines the signed distance function of a half plane
 */
class SDF_HalfPlane : public SignedDistanceFunction
{

public:
	SDF_HalfPlane(const Vector3d& o, const Vector3d& n_); ///< half plane by origin and normal vector
    SDF_HalfPlane(const Vector3d& o, const Vector3d& p1, const Vector3d& p2);  ///< half plane by origin and two linear independent vectors

    virtual double getDist(const Vector3d& v) const { return n.times(v.minus(o)); } ///< @see SignedDistanceFunction::getDist

    virtual int writePVPScript(std::ostream & cout, int c=1) const; ///< @see SignedDistanceFunction::writePVPScript

    virtual std::string toString() { return "SDF_HalfPlane"; } ///< @see SignedDistanceFunction::toString

    Vector3d o; ///< origin of the plane
    Vector3d n; ///< normal of the plane p1xp2
    Vector3d p1; ///< for visualisation (plane is a simplex: o,p1,p2)
    Vector3d p2; ///< for visualisation (plane is a simplex: o,p1,p2)

};


#endif
