// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef TROPISM_H
#define TROPISM_H

#include "mymath.h"

#include <chrono>
#include <iostream>
#include <vector>

namespace CRootBox {

class SoilLookUp;
class SignedDistanceFunction;
class Organ;
class Organism;

/**
 * Base class for all tropism functions, e.g. Gravitropism, Plagiotropism, Exotropism...
 */
class Tropism
{
public:

    Tropism(Organism* plant):Tropism(plant, 0,0) { } ///< Default constructor is TropismFunction(0,0)

    /**
     * Tropism with n_ number of trials and standard deviation of sigma_
     *
     * Always call the constructor, when overwriting the class!
     * Otherwise it will not work, and the mistake is hard to find.
     *
     * @param n_            number of tries
     * @param sigma_        standard deviation of angular change [1/cm]
     */
    Tropism(Organism* plant, double n_,double sigma_): plant(plant), n(n_), sigma(sigma_), geometry(nullptr) { }
    virtual ~Tropism() {};
    virtual Tropism* copy(Organism* plant); ///< copy object, factory method

    /* parameters */
    void setGeometry(SignedDistanceFunction* geom) { geometry = geom; } ///< sets a confining geometry
    void setTropismParameter(double n_,double sigma_) { n=n_; sigma=sigma_; } ///< sets the tropism parameters

    virtual Vector2d getHeading(const Vector3d& pos, Matrix3d old,  double dx, const Organ* o = nullptr);
    ///< constrained heading, dices n times and takes the best shot (according to the objective function)
    virtual Vector2d getUCHeading(const Vector3d& pos, Matrix3d old, double dx, const Organ* o);
    ///< Get unconfined heading (called by getHeading)

    /**
     * The objective function of the random optimization of getHeading(). Overwrite this function to implement a tropism.
     *
     * @param pos      current root tip position
     * @param old      rotation matrix, old(:,1) is the root tip heading
     * @param a        rotation angle alpha (angular change)
     * @param b        rotation angle beta (radial change)
     * @param dx       small distance to look ahead
     * @param root     points to the root that called getHeading, just in case something else is needed (i.e. iheading for exotropism)
     *
     * \return         the value minimized by getHeading(), it should be in [0,1], in this way combination of various tropisms will be easier
     */
    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* o = nullptr) { std::cout << "TropismFunction::tropismObjective() not overwritten\n"; return 0; }
    ///< The objective function of the random optimization of getHeading().

    static Vector3d getPosition(const Vector3d& pos, Matrix3d old, double a, double b, double dx);
    ///< Auxiliary function: Applies angles a and b and goes dx [cm] into the new direction

protected:

    Organism* plant;

    double n; ///< Number of trials
    double sigma; ///< Standard deviation

    SignedDistanceFunction* geometry; ///< confining geometry
    const int alphaN = 20;
    const int betaN = 5;

};



/**
 * Gravitropism: the tendency to grow downwards
 */
class Gravitropism : public Tropism
{

public:

    Gravitropism(Organism* plant, double n, double sigma) : Tropism(plant, n,sigma) { } ///< @see TropismFunction

    virtual Tropism* copy(Organism* plant) override {
        Gravitropism* nt = new Gravitropism(*this); // default copy constructor
        nt->plant = plant;
        return nt;
    } ///< copy constructor

    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* o = nullptr) override {
        return 0.5*(old.times(Vector3d::rotAB(a,b)).z+1.); // negative values point downwards, transformed to 0..1
    }
    ///< TropismFunction::getHeading minimizes this function, @see TropismFunction::getHeading and @see TropismFunction::tropismObjective



};



/**
 * Plagiotropism: the tendency to stay in a horicontal layer
 */
class Plagiotropism : public Tropism
{

public:

    Plagiotropism(Organism* plant,double n, double sigma) : Tropism(plant, n,sigma) { } ///< @see TropismFunction

    virtual Tropism* copy(Organism* plant) override {
        Plagiotropism* nt = new Plagiotropism(*this); // default copy constructor
        nt->plant = plant;
        return nt;
    } ///< copy constructor

    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* o = nullptr) override {
        return std::abs(old.times(Vector3d::rotAB(a,b)).z); // 0..1
    }
    ///< getHeading() minimizes this function, @see TropismFunction

};



/**
 * Exotropism: the tendency to keep the initial heading
 */
class Exotropism : public Tropism
{

public:

    Exotropism(Organism* plant, double n, double sigma) : Tropism(plant, n,sigma) { } ///< @see TropismFunction

    virtual Tropism* copy(Organism* plant) override {
        Exotropism* nt = new Exotropism(*this); // default copy constructor
        nt->plant = plant;
        return nt;
    } ///< copy constructor

    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* o = nullptr) override;
    ///< getHeading() minimizes this function, @see TropismFunction

};



/**
 * Hydrotropism (or Chemotropism, ...): the tendency to grow towards a higher saturation (or concentration, ...)
 */
class Hydrotropism : public Tropism
{

public:

    Hydrotropism(Organism* plant, double n, double sigma, SoilLookUp* soil) : Tropism(plant, n,sigma), soil(soil) { } ///< @see TropismFunction

    virtual Tropism* copy(Organism* plant) override {
        Hydrotropism* nt = new Hydrotropism(*this); // default copy constructor
        nt->plant = plant;
        return nt;
    } ///< copy constructor

    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* o = nullptr) override;
    ///< getHeading() minimizes this function, @see TropismFunction

private:
    SoilLookUp* soil;
};



/**
 * Combined tropisms, creates a linear combination of the respective objective functions
 */
class CombinedTropism : public Tropism
{

public:

    CombinedTropism(Organism* plant, double n, double sigma, std::vector<Tropism*> tropisms_, std::vector<double> weights_): Tropism(plant,n,sigma), tropisms(tropisms_), weights(weights_) {
        assert(tropisms.size()>0);
        assert(weights.size()>0);
        assert(tropisms.size()==weights.size());
    } ///< linearly combines the objective functions of multiple tropisms

    CombinedTropism(Organism* plant, double n, double sigma, Tropism* t1, double w1, Tropism* t2, double w2);
    ///< linearly combines the objective functions of two tropism functions

    virtual Tropism* copy(Organism* plant) override {
        CombinedTropism* nt = new CombinedTropism(*this); // default copy constructor
        for (size_t i=0; i<tropisms.size(); i++) {
            nt->tropisms[i] = tropisms[i]->copy(plant);
        }
        nt->plant = plant;
        return nt;
    } ///< copy constructor

    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* o = nullptr) override;
    ///< getHeading() minimizes this function, @see TropismFunction

private:
    std::vector<Tropism*> tropisms;
    std::vector<double> weights;
};

} // end namespace CRootBox

#endif
