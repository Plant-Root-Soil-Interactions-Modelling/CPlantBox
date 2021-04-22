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
    virtual double getLength(double t, double r, double k, std::shared_ptr<Organ> o, double CWGr =-1., double CWLength =-1.) const
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
    virtual double getAge(double l, double r, double k, std::shared_ptr<Organ> o, double CWdt=-1, double CWage=-1) const
    { throw std::runtime_error( "getAge() not implemented" ); return 0; } ///< Returns the age of a root of length l


    virtual std::shared_ptr<GrowthFunction> copy() const { return std::make_shared<GrowthFunction>(*this); } ///< Copy the object
};



/**
 * LinearGrowth elongates at constant rate until the maximal length k is reached
 */
class LinearGrowth : public GrowthFunction
{
public:

    double getLength(double t, double r, double k, std::shared_ptr<Organ> o, double CWGr =-1., double CWLength =-1.) const override { return std::min(k,r*t); } ///< @copydoc GrowthFunction::getLegngth

    double getAge(double l, double r, double k, std::shared_ptr<Organ> o, double CWdt=-1, double CWage=-1)  const override { return l/r; } ///< @copydoc GrowthFunction::getAge

    std::shared_ptr<GrowthFunction> copy() const override { return std::make_shared<LinearGrowth>(*this); } ///< @copydoc GrowthFunction::copy

};


/**
 * ExponentialGrowth elongates initially at constant rate r and slows down towards the maximum length k
 */
class ExponentialGrowth : public GrowthFunction
{
public:

    double getLength(double t, double r, double k, std::shared_ptr<Organ> o, double CWGr =-1., double CWLength =-1.) const override { return k*(1-exp(-(r/k)*t)); } ///< @copydoc GrowthFunction::getLegngth

    double getAge(double l, double r, double k, std::shared_ptr<Organ> o, double CWdt=-1, double CWage=-1) const override { ///< @copydoc GrowthFunction::getAge

        double age = - k/r*log(1-l/k);
        if (std::isfinite(age)) { // the age can not be computed when root length approaches max length
            return age;
        } else {
            return 1.e9; // very old
        }
    } ///< @see GrowthFunction

    std::shared_ptr<GrowthFunction> copy() const override { return std::make_shared<ExponentialGrowth>(*this); }

};


/**
 * CWLimitedGrowth uses growth given by phloem module
 */
class CWLimitedGrowth : public LinearGrowth
{
public:
    double getLength(double t, double r, double k, std::shared_ptr<Organ> o, double CWGr =-1., double CWLength =-1.) const override { 
		double CW_Gr= CWGr;
		double CW_length = CWLength;
		std::cout<<"\n grwoth::getlength "<<CW_Gr<<" "<<CW_length<< " "<<CW_Gr + CW_length;
		//double CW_length= o->getParameter("length");
		if (CW_Gr == -1.){
			std::cout<<"\n grwoth::fauil ";
			double length = LinearGrowth::getLength(t, r, k, o);
			return length;
		} else {return CW_Gr + CW_length; }
	}		///< @copydoc GrowthFunction::getLegngth

    double getAge(double l, double r, double k, std::shared_ptr<Organ> o, double CWdt=-1, double CWage=-1) const override {
		double CW_dt= CWdt;
		double CW_age= CWage;
		if ( CW_dt == -1){
			double age = LinearGrowth::getAge(l, r, k, o);
			return age;
		} else {return CW_age + CW_dt;}
	}

    std::shared_ptr<GrowthFunction> copy() const override { return std::make_shared<CWLimitedGrowth>(*this); }

};

} // end namespace CPlantBox

#endif
