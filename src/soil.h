#ifndef SOIL_H
#define SOIL_H

#include "mymath.h"
#include <cmath>
#include "sdf.h"

namespace CPlantBox {


class Organ;



/**
 * Look up method for a scalar soil property
 */
class SoilLookUp
{
public:
    virtual ~SoilLookUp() {};


    /**
     * Returns a scalar property of the soil scaled from 0..1
     *
     * @param pos       position [cm], (normally, root->getNode(root->getNumberOfNodes()-1))
     * @param root      the root that wants to know the scalar property
     *                  in some situation this might be usefull (e.g. could increase look up speed from a unstructured mesh)
     * \return          scalar soil property
     */
    virtual double getValue(const Vector3d& pos, const Organ* root = nullptr) const { return 1.; } ///< Returns a scalar property of the soil, 1. per default

    virtual std::string toString() const { return "SoilLookUp base class"; } ///< Quick info about the object for debugging

};



/**
 * A static soil property that is defined by a signed distance function
 */
class SoilLookUpSDF : public SoilLookUp
{
public:
    SoilLookUpSDF(): SoilLookUpSDF(nullptr) { } ///< Default constructor

    /**
     * Creaets the soil property from a signed distance function,
     * inside the geometry the value is largest
     *
     * @param sdf_      the signed distance function representing the geometry
     * @param max_      the maximal value of the soil property
     * @param min_      the minimal value of the soil property
     * @param slope_    scales the linear gradient of the sdf (note that |grad(sdf)|= 1)
     */
    SoilLookUpSDF(SignedDistanceFunction* sdf_, double max_=1, double min_=0, double slope_=1) {
        this->sdf=sdf_;
        fmax = max_;
        fmin = min_;
        slope = slope_;
    } ///< Creates the soil property from a signed distance function

    virtual double getValue(const Vector3d& pos, const Organ* root = nullptr) const override {
        double c = -sdf->getDist(pos)/slope*2.; ///< *(-1), because inside the geometry the value is largest
        c += (fmax-fmin)/2.; // thats the value at the boundary
        return std::max(std::min(c,fmax),fmin);
    } ///< returns fmin outside of the domain and fmax inside, and a linear ascend according slope

    virtual std::string toString() const { return "SoilLookUpSDF"; } ///< Quick info about the object for debugging

    SignedDistanceFunction* sdf; ///< signed distance function representing the geometry
    double fmax; ///< maximum is reached within the geometry at the distance slope
    double fmin; ///< minimum is reached outside of the geometry at the distance slope
    double slope; ///< half length of linear interpolation between fmax and fmin
};



/**
 * SoilLookUpSDF scaled from 0..1
 */
class ScaledSoilLookUpSDF : public SoilLookUpSDF
{
public:
	virtual double getValue(const Vector3d& pos, const Organ* root = nullptr) const override {
		double v = SoilLookUpSDF::getValue(pos, root);
		return (v-fmin)/(fmax-fmin);
	} ///< SoilLookUpSDF::getValue() but scaled from 0 to 1

	virtual std::string toString() const { return "ScaledSoilLookUpSDF"; } ///< Quick info about the object for debugging

};


/**
 *  1D linear look up
 */
class SoilLookUp1Dlinear : public SoilLookUp
{
public:

	SoilLookUp1Dlinear(double a, double b, size_t n): a(a), b(b), n(n), data(n) { };
	void setData(std::vector<double>& data_) { assert(data.size()==n); data=data_; }

	double mapIZ(size_t i) const { return a + (i+0.5)*(b-a)/(n-1); }
	size_t mapZI(double z) const { return std::round((z-a)*(n-1)/(b-a)-0.5); }

	virtual double getValue(const Vector3d& pos, const Organ* root = nullptr) const override {
		size_t i =mapZI(pos.z);
		assert(i>=0);
		assert(i<n);
		return data[i];
	}

    virtual std::string toString() const { return "SoilLookUp1Dlinear"; } ///< Quick info about the object for debugging

    double a;
	double b;
	size_t n;
	std::vector<double> data; ///< look up data

};



} // namespace CPlantBox

#endif
