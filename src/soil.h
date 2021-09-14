// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef SOIL_H
#define SOIL_H

#include "mymath.h"
#include "sdf.h"

#include <cmath>
#include <limits>
#include <memory>

namespace CPlantBox  {

class Organ;

/**
 * Base class to look up for a scalar soil property
 */
class SoilLookUp
{
public:

    const double inf = std::numeric_limits<double>::infinity();

    SoilLookUp() { }
    virtual ~SoilLookUp() { }

    virtual std::shared_ptr<SoilLookUp> copy() { return std::make_shared<SoilLookUp>(*this); }

    /**
     * Returns a scalar property of the soil scaled from 0..1.
     *
     * @param pos       position [cm], (normally, root->getNode(root->getNumberOfNodes()-1))
     * @param root      the root that wants to know the scalar property
     *                  in some situation this might be usefull (e.g. could increase look up speed from a unstructured mesh)
     * \return          scalar soil property
     */
    virtual double getValue(const Vector3d& pos, const std::shared_ptr<Organ> organ = nullptr) const { return 1.; } ///< Returns a scalar property of the soil, 1. per default

    virtual std::string toString() const { return "SoilLookUp base class"; } ///< Quick info about the object for debugging

    /**
     * sets the periodic boundaries, periodicity is used if bounds are set,
     * and if it is supportet by the getValue method of the derived classes
     *
     * periodic in x, y, and z
     */
    void setPeriodicDomain(double minx_, double maxx_, double miny_ , double maxy_, double minz_, double maxz_) {
        periodic_ = true;
        minx = minx_;
        xx = maxx_-minx;
        miny = miny_;
        yy = maxy_-miny;
        minz = minz_;
        zz = maxz_-minz;
    }

    /**
     * periodic in x and y
     */
    void setPeriodicDomain(double minx, double maxx, double miny , double maxy) {
        this->setPeriodicDomain(minx, maxx, miny, maxy, 0, inf );
    }

    /**
     * periodic in x
     */
    void setPeriodicDomain(double minx, double maxx) {
        this->setPeriodicDomain(minx, maxx, 0, inf, 0, inf );
    }

    /**
     * maps the point into the periodic domain
     */
    Vector3d periodic(const Vector3d& pos) const {  //< maps point into periodic domain
        if (periodic_) {
            Vector3d p = pos;
            if (!std::isinf(xx)) { // periodic in x
                p.x -= minx;
                p.x = (p.x/xx - (int)(p.x/xx))*xx;
                p.x += minx;
            }
            if (!std::isinf(yy)) { // periodic in y
                p.y -= miny;
                p.y = (p.y/yy - (int)(p.y/yy))*yy;
                p.y += miny;
            }
            if (!std::isinf(zz)) { // periodic in z
                p.z -= minz;
                p.z = (p.z/zz - (int)(p.z/zz))*zz;
                p.z += minz;
            }
            return p;
        } else {
            return pos;
        }
    }

private:

    bool periodic_ = false;
    double minx=0., xx=0., miny=0., yy=0, minz=0., zz=0.;

};



/**
 * Looks up a value based on a signed distance function
 */
class SoilLookUpSDF : public SoilLookUp
{
public:

    SoilLookUpSDF(): SoilLookUpSDF(nullptr) { } ///< Default constructor

    std::shared_ptr<SoilLookUp> copy() override { return std::make_shared<SoilLookUpSDF>(*this); }  // todo?: now its a shallow copy

    /**
     * Creates the soil property from a signed distance function,
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
    }

    /**
     * returns fmin outside of the domain and fmax inside, and a linear ascend according slope
     */
    double getValue(const Vector3d& pos, const std::shared_ptr<Organ> o  = nullptr) const override {
        Vector3d p = periodic(pos);
        double c = -sdf->getDist(p)/slope*2.; ///< *(-1), because inside the geometry the value is largest
        c += (fmax-fmin)/2.; // thats the value at the boundary
        return std::max(std::min(c,fmax),fmin);
    }

    std::string toString() const override { return "SoilLookUpSDF"; } ///< Quick info about the object for debugging

    SignedDistanceFunction* sdf; ///< signed distance function representing the geometry
    double fmax; ///< maximum is reached within the geometry at the distance slope
    double fmin; ///< minimum is reached outside of the geometry at the distance slope
    double slope; ///< half length of linear interpolation between fmax and fmin
};



/**
 * Product of multiple SoilLookUp::getValue
 */
class MultiplySoilLookUps : public SoilLookUp
{
public:

    MultiplySoilLookUps(SoilLookUp* s1, SoilLookUp* s2) { soils = { s1, s2 }; }

    MultiplySoilLookUps(std::vector<SoilLookUp*> soils) :soils(soils) { }

    std::shared_ptr<SoilLookUp> copy() override { return std::make_shared<MultiplySoilLookUps>(*this); } // todo? now its a shallow copy

    double getValue(const Vector3d& pos, const std::shared_ptr<Organ> o = nullptr) const override {
        double v = 1.;
        for (size_t i=0; i<soils.size(); i++) {
            v *= soils[i]->getValue(pos, o);
        }
        return v;
    }

    std::string toString() const override {
        std::string str = "";
        for (size_t i=0; i<soils.size(); i++) {
            str += soils[i]->toString();
            str += "; ";
        }
        return "MultiplySoilLookUps: "+ str;
    } ///< Quick info about the object for debugging

protected:

    std::vector<SoilLookUp*> soils;

};



/**
 * Scales the root elongation with fixed value same for each root
 */
class ProportionalElongation : public SoilLookUp
{
public:

    ProportionalElongation() {
    }

    ProportionalElongation(double scale)
    :scale(scale) {
    }

    ProportionalElongation(double scale, SoilLookUp* baseLookUp)
    :scale(scale), baseLookUp(baseLookUp) {
    }

    std::shared_ptr<SoilLookUp> copy() override { return std::make_shared<ProportionalElongation>(*this); } // todo? now its a shallow copy

    void setScale(double s) { scale = s; }

    void setBaseLookUp(SoilLookUp* baseLookUp) { this->baseLookUp=baseLookUp; } ///< proportionally scales a base soil look up

    double getValue(const Vector3d& pos, const std::shared_ptr<Organ> o = nullptr) const override {
        if (baseLookUp==nullptr) {
            return scale;
        } else {
            Vector3d p = this->periodic(pos);
            return baseLookUp->getValue(p, o)*scale;  // superimpose scaling on a base soil look up function
        }
    }

    std::string toString() const override { return "ProportionalElongation"; } ///< Quick info about the object for debugging

protected:

    double scale = 1.;
    SoilLookUp* baseLookUp = nullptr;

};



/**
 * 1D look up table, where data is located between the grid points
 */
class Grid1D  : public SoilLookUp
{
public:

    Grid1D() {
        n=0;
        data = std::vector<double>(0);
        grid = std::vector<double>(0);
    }

    Grid1D(size_t n, std::vector<double> grid, std::vector<double> data): n(n), grid(grid), data(data) {
        assert(grid.size()==n);
        assert(data.size()==n);
    }


    std::shared_ptr<SoilLookUp> copy()  override { return std::make_shared<Grid1D>(*this); }

    virtual size_t map(double x) const {
        unsigned int jr,jm,jl;
        jl = 0;
        jr = n - 1;
        while (jr-jl > 1) {
            jm=(jr+jl) >> 1; // thats a divided by two
            if (x >= grid[jm])
                jl=jm;
            else
                jr=jm;
        }
        return jl;
    } ///< Generic way to perform look up in an ordered table, overwrite by faster method if appropriate, todo currently floor

    double getValue(const Vector3d& pos, const std::shared_ptr<Organ> o = nullptr) const override {
        Vector3d p = this->periodic(pos);
        return data[map(p.z)];
    } ///< Returns the data of the 1d table, repeats first or last entry if out of bound

    std::string toString() const override { return "RectilinearGrid1D"; } ///< Quick info about the object for debugging

    size_t n;
    std::vector<double> grid;
    std::vector<double> data;
};



/**
 *  1D look up table with equidistant spacing, data is located between the grid points
 */
class EquidistantGrid1D : public Grid1D
{
public:

    EquidistantGrid1D(): EquidistantGrid1D(0,1,0) {
    }

    EquidistantGrid1D(double a, double b, size_t n): a(a), b(b) {
        this->n =n;
        makeGrid(a,b,n);
        this->data = std::vector<double>(n);
    }

    EquidistantGrid1D(double a, double b, const std::vector<double>& data): a(a), b(b) {
        this->n = data.size();
        makeGrid(a,b,n);
        this->data = data;
    }

    std::shared_ptr<SoilLookUp> copy()  override { return std::make_shared<EquidistantGrid1D>(*this); }

    void makeGrid(double a, double b, size_t n) {
        this->grid = std::vector<double>(n);
        for (size_t i=0; i<n; i++) {
            grid[i] = a + (b-a)/double(n-1)*i;
        }
    }

    size_t map(double x) const override {
        return std::floor((x-a)/(b-a)*(n-1));
    } ///< faster than general look up

    std::string toString() const  override{ return "LinearGrid1D"; } ///< Quick info about the object for debugging

    double a;
    double b;

};


/**
 * RectilinearGrid, called tensor product grid (in Dune), data is located between the grid points
 */
class RectilinearGrid3D  : public SoilLookUp
{
public:

    RectilinearGrid3D(Grid1D* xgrid, Grid1D* ygrid, Grid1D* zgrid) :xgrid(xgrid), ygrid(ygrid), zgrid(zgrid) {
        nx = xgrid->n;
        ny = ygrid->n;
        nz = zgrid->n;
        data = std::vector<double>(nx*ny*nz);
    }

    virtual ~RectilinearGrid3D() { };

    std::shared_ptr<SoilLookUp> copy()  override { return std::make_shared<RectilinearGrid3D>(*this); }

    virtual size_t map(double x, double y, double z) const {
        size_t i = xgrid->map(x);
        size_t j = ygrid->map(y);
        size_t k = zgrid->map(z);
//        std::cout << "RectilinearGrid3D::map: " << i << ", " << ", " << j << ", " << k <<
//        		" = " << k*(nx*ny)+j*nx+i << "; " << nx << ", " << ny << "\n";
        return k*(nx*ny)+j*nx+i; // or whatever
    } ///< point to linear data index

    double getData(size_t i, size_t j, size_t k) {
        return data.at(map(i,j,k));
    }

    double getValue(const Vector3d& pos, const std::shared_ptr<Organ> o = nullptr) const override {
        Vector3d p = periodic(pos);
        return data[map(p.x,p.y,p.z)];
    } ///< Returns the data of the 1d table, repeats first or last entry if out of bound

    void setData(size_t i, size_t j, size_t k, double d) {
        data.at(map(i,j,k)) = d;
    }

    Vector3d getGridPoint(size_t i, size_t j, size_t k) {
        return Vector3d(xgrid->grid[i], ygrid->grid[j], zgrid->grid[k]);
    } ///< grid point at indices

    Grid1D* xgrid;
    Grid1D* ygrid;
    Grid1D* zgrid;

    size_t nx,ny,nz;
    std::vector<double> data;
};



/**
 *  Even more fantastic
 */
class EquidistantGrid3D : public RectilinearGrid3D
{
public:

    EquidistantGrid3D() :EquidistantGrid3D(1.,1.,1.,0,0,0) {
    }

    EquidistantGrid3D(double length, double width, double depth, int nx, int ny, int nz)
    :RectilinearGrid3D(new EquidistantGrid1D(-length/2,length/2,nx),new EquidistantGrid1D(-width/2,width/2,ny),new EquidistantGrid1D(-depth,0.,nz)) {
    }

    EquidistantGrid3D(double x0, double xe, int nx, double y0, double ye, int ny, double z0, double ze, int nz)
    :RectilinearGrid3D(new EquidistantGrid1D(x0,xe,nx),new EquidistantGrid1D(y0,ye,ny),new EquidistantGrid1D(z0,ze,nz)) {
    }

    virtual ~EquidistantGrid3D() {
        delete xgrid;
        delete ygrid;
        delete zgrid;
    }

    std::shared_ptr<SoilLookUp> copy()  override { return std::make_shared<EquidistantGrid3D>(*this); }

};

} // end namespace CPlantBox

#endif
