// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef TROPISM_H
#define TROPISM_H

#include "mymath.h"
#include "Organ.h"
#include "Organism.h"

#include <memory>
#include <chrono>
#include <iostream>
#include <vector>
#include <random>

namespace CPlantBox {

class SoilLookUp;
class SignedDistanceFunction;

/**
 * Base class for all tropism functions, e.g. Gravitropism, Plagiotropism, Exotropism...
 */
class Tropism
{
public:

	Tropism(std::shared_ptr<Organism> plant):Tropism(plant, 0,0) { } ///< Default constructor is TropismFunction(0,0)

	/**
	 * Tropism with n_ number of trials and standard deviation of sigma_
	 *
	 * Always call the constructor, when overwriting the class!
	 * Otherwise it will not work, and the mistake is hard to find.
	 *
	 * @param n_            number of tries
	 * @param sigma_        standard deviation of angular change [1/cm]
	 */
	Tropism(std::shared_ptr<Organism> plant, double n_,double sigma_, double ageSwitch_ = 0):
	    ageSwitch(ageSwitch_),plant(plant), n(n_), sigma(sigma_) { }
	virtual ~Tropism() { }

	virtual std::shared_ptr<Tropism> copy(std::shared_ptr<Organism> plant); ///< copy object, factory method

	/* parameters */
	void setGeometry(std::shared_ptr<SignedDistanceFunction> geom) { geometry = geom; } ///< sets a confining geometry
	void setTropismParameter(double n_,double sigma_) { n=n_; sigma=sigma_; } ///< sets the tropism parameters
	void setSigma(double newSigma) { sigma = newSigma; }

	virtual Vector2d getHeading(const Vector3d& pos, const Matrix3d& old,  double dx, const std::shared_ptr<Organ> o = nullptr, int nodeIdx = -1);
	///< constrained heading, dices n times and takes the best shot (according to the objective function)
	virtual Vector2d getUCHeading(const Vector3d& pos, const Matrix3d& old, double dx, const std::shared_ptr<Organ> o, int nodeIdx );
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
	virtual double tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> o = nullptr)
	    { std::cout << "TropismFunction::tropismObjective() not overwritten\n"; return 0; } ///< The objective function of the random optimization of getHeading().

	static Vector3d getPosition(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx);
	///< Auxiliary function: Applies angles a and b and goes dx [cm] into the new direction

	int alphaN = 20; //stop protecting in case want to increase number of trials => very important to respect soil boundaries when using photosynthesis
	int betaN = 5; //stop protecting in case want to increase number of trials

	bool isExpired() { return (!plant.lock()); } // for debugging
	std::shared_ptr<Organism> getPlant() { return plant.lock(); }

    double ageSwitch; // DL: why not put this in the specialisation (i.g. AntiGravi2Gravitropism)?


protected:

	std::weak_ptr<Organism> plant;

	double n; ///< Number of trials
	double sigma; ///< Standard deviation

	std::weak_ptr<SignedDistanceFunction> geometry; ///< confining geometry
	double randn(int nNode) { if((nNode > 0) && (plant.lock()->getStochastic())) { return ND(gen);} else { return plant.lock()->randn();} } ///< normally distributed random number (0,1)
    double rand(int nNode) { if((nNode > 0) && (plant.lock()->getStochastic())) { return UD(gen);} else { return plant.lock()->rand();} } ///< uniformly distributed random number (0,1)

    std::normal_distribution<double> ND;
    std::uniform_real_distribution<double> UD;
	std::mt19937 gen; ///< random number generator

};



/**
 * Gravitropism: the tendency to grow downwards
 */
class Gravitropism : public Tropism
{
public:

	Gravitropism(std::shared_ptr<Organism> plant, double n, double sigma) : Tropism(plant, n,sigma) { } ///< @see TropismFunction

	std::shared_ptr<Tropism> copy(std::shared_ptr<Organism> plant) override {
		auto nt = std::make_shared<Gravitropism>(*this); // default copy constructor
        nt->plant  =  std::weak_ptr<Organism>(); // necessary?
		nt->plant = plant;
		//std::cout << "Copy tropism: from " << this->plant.lock()->plantId << " to " << nt->plant.lock()->plantId  << "\n";
		return nt;
	} ///< copy constructor

	double tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> o = nullptr) override {
		return 0.5*(old.times(Vector3d::rotAB(a,b)).z+1.); // negative values point downwards, transformed to 0..1
	} ///< TropismFunction::getHeading minimizes this function, @see TropismFunction::getHeading and @see TropismFunction::tropismObjective

};



/**
 * Plagiotropism: the tendency to stay in a horicontal layer
 */
class Plagiotropism : public Tropism
{
public:

	Plagiotropism(std::shared_ptr<Organism> plant,double n, double sigma) : Tropism(plant, n,sigma) { } ///< @see TropismFunction

	std::shared_ptr<Tropism>  copy(std::shared_ptr<Organism> plant) override {
	    auto nt = std::make_shared<Plagiotropism>(*this); // default copy constructor
        nt->plant  =  std::weak_ptr<Organism>(); // necessary?
		nt->plant = plant;
		// std::cout << "Copy tropism: from " << this->plant.lock()->plantId << " to " << nt->plant.lock()->plantId  << "\n";
		return nt;
	} ///< copy constructor

	double tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> o = nullptr) override {
		return std::abs(old.times(Vector3d::rotAB(a,b)).z); // 0..1
	} ///< getHeading() minimizes this function, @see TropismFunction

};



/**
 * Exotropism: the tendency to keep the initial heading
 */
class Exotropism : public Tropism
{
public:

	Exotropism(std::shared_ptr<Organism> plant, double n, double sigma) : Tropism(plant, n,sigma) { } ///< @see TropismFunction

    std::shared_ptr<Tropism> copy(std::shared_ptr<Organism> plant) override {
        auto nt = std::make_shared<Exotropism>(*this); // default copy constructor
        nt->plant  =  std::weak_ptr<Organism>(); // necessary?
        nt->plant = plant;
        // std::cout << "Copy tropism: from " << this->plant.lock()->plantId << " to " << nt->plant.lock()->plantId << "\n";
        return nt;
    } ///< copy constructor

	double tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> o = nullptr) override;
	///< getHeading() minimizes this function, @see TropismFunction

};



/**
 * Hydrotropism (or Chemotropism, ...): the tendency to grow towards a higher saturation (or concentration, ...)
 */
class Hydrotropism : public Tropism
{
public:

	Hydrotropism(std::shared_ptr<Organism> plant, double n, double sigma, std::shared_ptr<SoilLookUp> soil) : Tropism(plant, n,sigma), soil(soil) { } ///< @see TropismFunction

    std::shared_ptr<Tropism>  copy(std::shared_ptr<Organism> plant) override {
        auto nt = std::make_shared<Hydrotropism>(*this); // default copy constructor
        nt->plant = plant;
        return nt;
    } ///< copy constructor

	double tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> o = nullptr) override;
	///< getHeading() minimizes this function, @see TropismFunction

private:

	std::weak_ptr<SoilLookUp> soil;

};



/**
 * Combined tropisms, creates a linear combination of the respective objective functions
 */
class CombinedTropism : public Tropism
{
public:

	CombinedTropism(std::shared_ptr<Organism> plant, double n, double sigma, std::vector<std::shared_ptr<Tropism>> tropisms_, std::vector<double> weights_):
	        Tropism(plant,n,sigma), tropisms(tropisms_), weights(weights_) {
		assert(tropisms.size()>0);
		assert(weights.size()>0);
		assert(tropisms.size()==weights.size());
	} ///< linearly combines the objective functions of multiple tropisms

	CombinedTropism(std::shared_ptr<Organism> plant, double n, double sigma, std::shared_ptr<Tropism> t1, double w1, std::shared_ptr<Tropism> t2, double w2);
	///< linearly combines the objective functions of two tropism functions

	std::shared_ptr<Tropism>  copy(std::shared_ptr<Organism> plant) override {
        auto nt = std::make_shared<CombinedTropism>(*this); // default copy constructor
		for (size_t i=0; i<tropisms.size(); i++) {
			nt->tropisms[i] = tropisms[i]->copy(plant);
		}
		nt->plant = plant; //todo
		return nt;
	} ///< copy constructor, deep copies tropisms

	double tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> o = nullptr) override;
	///< getHeading() minimizes this function, @see TropismFunction

private:

	std::vector<std::shared_ptr<Tropism>> tropisms;
	std::vector<double> weights;

};

/**
 * todo doc
 */
class TwistTropism: public Tropism
{

public:

	TwistTropism(std::shared_ptr<Organism> plant, double n, double sigma) : Tropism(plant, n,sigma) { } ///< @see TropismFunction

    std::shared_ptr<Tropism>  copy(std::shared_ptr<Organism> plant) override {
        auto nt = std::make_shared<TwistTropism>(*this); // default copy constructor
        nt->plant = plant;
        return nt;
    } ///< copy constructor


	double tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> stem) override {
		//    old.times(Matrix3d::rotX(b));
		//    old.times(Matrix3d::rotZ(a));
		return  -0.9*(old.times(Vector3d::rotAB(a+0.5*rand(-1),b+0.5*rand(-1))).z+1.); // negative values point downwards, tranformed to 0..1
	}
	///< TropismFunction::getHeading minimizes this function, @see TropismFunction::getHeading and @see TropismFunction::tropismObjective

};



/**
 * todo doc
 */
class AntiGravitropism : public Tropism
{
public:

	AntiGravitropism(std::shared_ptr<Organism> plant, double n, double sigma) : Tropism(plant, n,sigma) { } ///< @see TropismFunction

    std::shared_ptr<Tropism>  copy(std::shared_ptr<Organism> plant) override {
        auto nt = std::make_shared<AntiGravitropism>(*this); // default copy constructor
        nt->plant = plant;
        return nt;
    } ///< copy constructor


	virtual double tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> stem) override {
		return  -0.5*(old.times(Vector3d::rotAB(a,b)).z+1.); // negative values point downwards, tranformed to 0..1
	}
	///< TropismFunction::getHeading minimizes this function, @see TropismFunction::getHeading and @see TropismFunction::tropismObjective

};


/**
 * for wheat leaves. start by growing straight and then fall to the side.
 */
class AntiGravi2Gravitropism : public Tropism
{
public:

	AntiGravi2Gravitropism(std::shared_ptr<Organism> plant, double n, double sigma, double ageSwitch = 0) :
	Tropism(plant, n,sigma,ageSwitch) { } ///< @see TropismFunction


    std::shared_ptr<Tropism>  copy(std::shared_ptr<Organism> plant) override {
        auto nt = std::make_shared<AntiGravi2Gravitropism>(*this); // default copy constructor
        nt->plant = plant;
        return nt;
    } ///< copy constructor


	virtual double tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> o) override {
		return  -0.5*(old.times(Vector3d::rotAB(a,b)).z+1.)*(o->getAge()<ageSwitch)
		+ 0.5*(old.times(Vector3d::rotAB(a,b)).z+1.)*(o->getAge()>=ageSwitch); // negative values point downwards, tranformed to 0..1
	}
	///< TropismFunction::getHeading minimizes this function, @see TropismFunction::getHeading and @see TropismFunction::tropismObjective

};


} // end namespace CPlantBox

#endif
