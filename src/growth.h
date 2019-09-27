#ifndef GROWTH_H
#define GROWTH_H

#include <memory>

namespace CPlantBox {

class Organ;

/**
 * Abstract base class to all growth functions: currently LinearGrowth and ExponentialGrowth
 */
class GrowthFunction
{
public:
    virtual ~GrowthFunction() {};

    /**
     * Returns root length at root age t
     *
     * @param t     organ age [day]
     * @param r     initial growth rate [cm/day]
     * @param k     maximal organ length [cm]
     * @param root  points to the organ in case more information is needed
     *
     * \return      organ length [cm] at specific age
     */
    virtual double getLength(double t, double r, double k, std::shared_ptr<Organ> o) const
    { throw std::runtime_error( "getLength() not implemented" ); return 0; } ///< Returns root length at root age t

    /**
     * Returns the age of a root of length l
     *
     * @param l     organ length [cm]
     * @param r     initial growth rate [cm/day]
     * @param k     maximal root length [cm]
     * @param root  points to the organ in case more information is needed
     *
     * \return      organ age [day] at specific length
     */
    virtual double getAge(double l, double r, double k, std::shared_ptr<Organ> o) const
    { throw std::runtime_error( "getAge() not implemented" ); return 0; } ///< Returns the age of a root of length l


    virtual std::shared_ptr<GrowthFunction> copy() const { return std::make_shared<GrowthFunction>(*this); } ///< Copy the object
};



/**
 * LinearGrowth elongates at constant rate until the maximal length k is reached
 */
class LinearGrowth : public GrowthFunction
{
public:

    double getLength(double t, double r, double k, std::shared_ptr<Organ> o) const override { return std::min(k,r*t); } ///< @copydoc GrowthFunction::getLegngth

    double getAge(double l, double r, double k, std::shared_ptr<Organ> o)  const override { return l/r; } ///< @copydoc GrowthFunction::getAge

    std::shared_ptr<GrowthFunction> copy() const override { return std::make_shared<LinearGrowth>(*this); } ///< @copydoc GrowthFunction::copy

};


/**
 * ExponentialGrowth elongates initially at constant rate r and slows down towards the maximum length k
 */
class ExponentialGrowth : public GrowthFunction
{
public:

    double getLength(double t, double r, double k, std::shared_ptr<Organ> o) const override { return k*(1-exp(-(r/k)*t)); } ///< @copydoc GrowthFunction::getLegngth

    double getAge(double l, double r, double k, std::shared_ptr<Organ> o) const override { ///< @copydoc GrowthFunction::getAge

        double age = - k/r*log(1-l/k);
        if (std::isfinite(age)) { // the age can not be computed when root length approaches max length
            return age;
        } else {
            return 1.e9; // very old
        }
    } ///< @see GrowthFunction

    std::shared_ptr<GrowthFunction> copy() const override { return std::make_shared<ExponentialGrowth>(*this); }

};

} // end namespace CPlantBox

#endif
