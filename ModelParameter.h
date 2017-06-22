#ifndef MODELPARAMETER_H_
#define MODELPARAMETER_H_

/*
 * The model parameters consist of
 * RootTypeParameter            parameters for a root type
 * RootParameters     			parameter values for a specific root
 * PlantParameter               parameters for single plant
 */

#include <string>
#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include <assert.h>

#include "mymath.h"
#include "soil.h"

class RootParameter;



/**
 * RootTypeParameter: contains a parameter set describing a root type
 */
class RootTypeParameter
{

public:

	RootTypeParameter(); ///< default constructor
	RootTypeParameter(const RootTypeParameter& rp); ///< copy constructor
	virtual ~RootTypeParameter() { };

	void set(int type, double lb, double lbs, double la, double las, double ln, double lns, double nob, double nobs,
			double r, double rs, double a, double as,  double colorR, double colorG, double colorB, double tropismT, double tropismN, double tropsimS,
			double dx, const std::vector<int>& successor, const std::vector<double>& successorP, double theta, double thetas, double rlt, double rlts,
			int gf, const std::string& name); ///< sets all parameters

	RootParameter realize(); ///< Creates a specific root from the root parameter set
	int getLateralType(const Vector3d& pos); ///< Choose (dice) lateral type based on root parameter set
	double getK() const { return std::max(nob-1,double(0))*ln+la+lb; }  ///< returns the mean maximal root length [cm]

	// IO
	void read(std::istream & cin); ///< reads a single root parameter set
	void write(std::ostream & cout) const; ///< writes a single root parameter set
	std::string toString() const { std::stringstream ss; write(ss); return ss.str(); } ///< writes parameter to a string

	// random numbers
	void setSeed(double seed) { gen.seed(seed); } ///< Sets the seed of the random number generator
	double rand() { return UD(gen); } ///< Uniformly distributed random number (0,1)
	double randn() { return ND(gen); } ///< Normally distributed random number (0,1)

	/*
	 * Rootbox parameters per root type
	 */
	int type; 		///< Number of root type [1], this is the index within the vector +1
	double lb; 	 	///< Basal zone [cm]
	double lbs;  	///< Standard deviation basal zone [cm]
	double la;		///< Apical zone [cm];
	double las;		///< Standard deviation apical zone [cm];
	double ln; 		///< Inter-lateral distance [cm]
	double lns;  	///< Standard deviation inter-lateral distance [cm]
	double nob; 	///< Number of branches [1]
	double nobs; 	///< Standard deviation number of branches [1]
	double r;		///< Initial growth rate [cm day-1]
	double rs;		///< Standard deviation initial growth rate [cm day-1]
	double a; 		///< Root radius [cm]
	double as; 		///< Standard deviation root radius [cm]
	double colorR;	///< Root color (red)
	double colorG;	///< Root color (green)
	double colorB;	///< Root color (blue)
	int tropismT;	///< Root tropism parameter (Type)
	double tropismN;///< Root tropism parameter (number of trials)
	double tropismS;///< Root tropism parameter (mean value of expected changeg) [1/cm]
	double dx; 		///< Maximal segment size [cm]
	double theta; 	///< Angle between root and parent root (rad)
	double thetas; 	///< Standard deviation angle between root and parent root (rad)
	double rlt;		///< Root life time (days)
	double rlts;	///< Standard deviation root life time (days)
	int gf;			///< Growth function (1=negative exponential, 2=linear)
	std::string name;	///< Root type name
	std::vector<int> successor;			///< Lateral types [1]
	std::vector<double> successorP; 	///< Probabiltities of lateral type to emerge (sum of values == 1) [1]

	SoilProperty* se = new SoilProperty(); ///< scale elongation function
	SoilProperty* sa = new SoilProperty(); ///< scale angle function
	SoilProperty* sbp = new SoilProperty(); ///< scale branching probability function

private:
	std::mt19937 gen = std::mt19937(std::chrono::system_clock::now().time_since_epoch().count());  // random stuff
	std::uniform_real_distribution<double> UD = std::uniform_real_distribution<double>(0,1);
	std::normal_distribution<double> ND = std::normal_distribution<double>(0,1);

};



/**
 * IndividualRootParameter:
 * contains parameters of a single root, that are created by RootParameter:realize()
 */
class RootParameter
{

public:

	/* Constructor */
	RootParameter(); ///< Default constructor
	RootParameter(int type, double lb, double la, const std::vector<double>& ln, int nob, double r, double a, double theta, double rlt):
		type(type), lb(lb), la(la), nob(nob), r(r), a(a), theta(theta), rlt(rlt), ln(ln) { }
	///< Constructor setting all parameters

	/* Methods */
	void set(int type, double lb, double la, const std::vector<double>& ln, double nob, double r, double a, double theta, double rlt);
	///< Sets all the parameters

	double getK() const; ///< Returns the exact maximal root length of this realization [cm]

	void write(std::ostream & cout) const; ///< Writes parameters for debugging
	std::string toString() const { std::stringstream ss; write(ss); return ss.str(); } ///< String representation using RootParameter::write

	/* Rootbox parameters per root */
	int type; 			///< Index of root type [1]
	double lb; 			///< Basal zone [cm]
	double la;			///< Apical zone [cm];
	int nob; 			///< Number of branches [1]
	double r;			///< Initial growth rate [cm day-1]
	double a; 			///< Root radius [cm]
	double theta; 		///< Angle between root and parent root [rad]
	double rlt;			///< Root life time [day]
	std::vector<double> ln;    ///< Inter-lateral distances [cm]

};



/**
 * RootSystemParameter: contains all plant specific parameters like planting depth and describing the emergence times of basal and shoot borne roots
 *
 * This model is very limited in the moment, and we have to replace it, if we come up with something better
 */
class RootSystemParameter
{

public:

	RootSystemParameter(); ///< Default constructor
	virtual ~RootSystemParameter();

	virtual void set(double pd, double fB, double dB, int mB, int nC, double fSB, double dSB, double dRC, double nz, double simtime); ///< Sets all the parameters

	virtual void read(std::istream & cin); ///< Read plant parameters
	virtual void write(std::ostream & cout) const; ///< Write plant parameters
	std::string toString() const { std::stringstream ss; write(ss); return ss.str(); } ///< writes parameter to a string

	/* Plant parameters */
	Vector3d seedPos;   ///< Position of the seed [cm]

	//Basal roots (nodal roots)
	double firstB; 	///< Emergence of first basal root [day]
	double delayB; 	///< Time delay between the basal roots [day]
	int maxB; 	    ///< Maximal number of basal roots [1]

	//Shoot borne roots (crown roots)
	int nC; 		    ///< Maximal number of roots per root crown [1]
	double firstSB; 	///< First emergence of a shoot borne root [day]
	double delaySB; 	///< Time delay between the shoot borne roots [day]
	double delayRC; 	///< Delay between the root crowns [day]
	double nz; 		    ///< Distance between the root crowns along the shoot [cm]

	//Simulation parameters
	double simtime;    ///< recommended final simulation time
};


#endif
