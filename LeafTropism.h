#ifndef LEAFTROPISM_H
#define LEAFTROPISM_H

#include <chrono>
#include <random>

#include "soil.h"

class SoilProperty;
class Organ;

/**
* Base class for all tropism functions, e.g. Gravitropism, Plagiotropism, Exotropism...
*/
class LeafTropismFunction
{
public:
  virtual ~LeafTropismFunction() {};

  /**
  * Constructor,
  * Always call the constructor, when overwriting the class!
  * Otherwise it will not work, and the mistake is hard to find.
  *
  * @param n_            number of tries
  * @param sigma_        standard deviation of angular change [1/cm]
  */
  LeafTropismFunction(double n_,double sigma_): n(n_), sigma(sigma_) { }; ///< Tropism with n_ number of trials and standard deviation of sigma_

  LeafTropismFunction():LeafTropismFunction(0,0) { } ///< Default constructor is TropismFunction(0,0)

  virtual Vector2d getHeading(const Vector3d& pos, Matrix3d old, double dx, const Organ* leaf);
  ///<  Dices n times and takes the best shot (according to the objective function)

  /**
  * The objective function of the random optimization of getHeading(). Overwrite this function to implement a tropism.
  *
  * @param pos      current root tip position
  * @param old      rotation matrix, old(:,1) is the root tip heading
  * @param a        rotation angle alpha (angular change)
  * @param b        rotation angle beta (radial change)
  * @param dx       small ditance to look ahead
  * @param root     points to the root that called getHeading, just in case something else is needed (i.e. iheading for exotropism)
  *
  * \return         the value minimized by getHeading(), it should be in [0,1], in this way combination of various tropisms will be easier
  */
  virtual double leaftropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* leaf) { std::cout << "default\n"; return 0; }
  ///< The objective function of the random optimization of getHeading().

  static Vector3d getPosition(const Vector3d& pos, Matrix3d old, double a, double b, double dx);
  //< Auxiliary function: Applies angles a and b and goes dx [cm] into the new direction

  // random numbers
  void setSeed(double seed) { gen.seed(seed); } ///< Sets the seed of the random number generator
  double rand() { return UD(gen); } ///< Uniformly distributed random number (0,1)
  double randn() { return ND(gen); } ///< Normally distributed random number (0,1)

protected:
  double n; ///< Number of trials
  double sigma; ///< Standard deviation

private:
  std::mt19937 gen = std::mt19937(std::chrono::system_clock::now().time_since_epoch().count());  // random stuff
  std::normal_distribution<double> ND = std::normal_distribution<double>(0,1);
  std::uniform_real_distribution<double> UD = std::uniform_real_distribution<double>(0,1);;
};



/**
* Confines root growth to a geometry
*/
class ConfinedLeafTropism : public LeafTropismFunction
{
public:

  /**
  * ConfinedTropism confines a tropism to a geometry given by a signed distance function
  *
  * @param baseTropism   the underlying tropism
  * @param geometry      geometry, e.g. of a plant container or an obstacle
  */
  ConfinedLeafTropism(LeafTropismFunction* baseTropism, SignedDistanceFunction* geometry) : LeafTropismFunction(0,0), leaftropism(baseTropism), geometry(geometry) { }
  //< ConfinedTropism confines the baseTropism to geometry

  virtual Vector2d getHeading(const Vector3d& pos, Matrix3d old,  double dx, const Organ* leaf) override;
  ///< changes baseTropism->getHeading() in case geometric boundaries are hit

  virtual double leaftropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* leaf) override {
    throw std::invalid_argument( "ConfinedTropism::tropismObjective() should not be called directly" );
  } ///< this mehtod is not called, instead baseTropism->tropismObjective() is used by getHeading

private:
  LeafTropismFunction* leaftropism;
  SignedDistanceFunction* geometry;
  const int alphaN = 20;
  const int betaN = 5;
};



/**
* Gravitropism: the tendency to grow downwards
*/
class LeafGravitropism : public LeafTropismFunction
{

public:

  LeafGravitropism(double n, double sigma) : LeafTropismFunction(n,sigma) { } ///< @see TropismFunction

  virtual double leaftropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* leaf) override {
    old.times(Matrix3d::rotX(b));
    old.times(Matrix3d::rotZ(a));
    return 0.5*(old.column(0).z+1.); // negative values point downwards, tranformed to 0..1
  }
  ///< TropismFunction::getHeading minimizes this function, @see TropismFunction::getHeading and @see TropismFunction::tropismObjective

};



/**
* Plagiotropism: the tendency to stay in a horicontal layer
*/
class LeafPlagiotropism : public LeafTropismFunction
{

public:

  LeafPlagiotropism(double n, double sigma) : LeafTropismFunction(n,sigma) { } ///< @see TropismFunction

  virtual double leaftropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* leaf) override {
    old.times(Matrix3d::rotX(b));
    old.times(Matrix3d::rotZ(a));
    return std::abs(old.column(0).z); // 0..1
  }
  ///< getHeading() minimizes this function, @see TropismFunction

};



/**
* Exotropism: the tendency to keep the initial heading
*/
class LeafExotropism : public LeafTropismFunction
{

public:

  LeafExotropism(double n, double sigma) : LeafTropismFunction(n,sigma) { } ///< @see TropismFunction

  virtual double leaftropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* leaf) override;
  ///< getHeading() minimizes this function, @see TropismFunction
};



/**
 * Hydrotropism (or Chemotropism, ...): the tendency to grow towards a higher saturation (or concentration, ...)
 */
class LeafPhototropism : public LeafTropismFunction
{

public:

  LeafPhototropism(double n, double sigma, SoilProperty* soil) : LeafTropismFunction(n,sigma), soil(soil) { } ///< @see TropismFunction

  virtual double leaftropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* leaf) override;
  ///< getHeading() minimizes this function, @see TropismFunction

private:
  SoilProperty* soil;
};



/**
 * Combined tropisms, creates a linear combination of the respective objective functions
 */
class CombinedLeafTropism : public LeafTropismFunction
{

public:

  CombinedLeafTropism(double n, double sigma, std::vector<LeafTropismFunction*> tropisms_, std::vector<double> weights_): LeafTropismFunction(n,sigma), tropisms(tropisms_), weights(weights_) {
    assert(tropisms.size()>0);
    assert(weights.size()>0);
    assert(tropisms.size()==weights.size());
  } ///< linearly comibines the objective functions of multiple tropisms
  CombinedLeafTropism(double n, double sigma, LeafTropismFunction* t1, double w1, LeafTropismFunction* t2, double w2);
  ///< linearly comibines the objective functions of two tropism funcitons

  virtual double leaftropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* leaf) override;
  ///< getHeading() minimizes this function, @see TropismFunction

private:
  std::vector<LeafTropismFunction*> tropisms;
  std::vector<double> weights;
};


#include "Organ.h"


#endif
