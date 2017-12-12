// *** ADDED BY HEADER FIXUP ***
#include <cassert>
#include <istream>
// *** END ***
#ifndef MODELPARAMETER_H_
#define MODELPARAMETER_H_

#include <string>
#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include <assert.h>

#include "mymath.h"

/*
 * The model parameters consist of
 *
 * OrganTypeParameter			an organ type
 * OrganParameter				parameter values for a specific organ (base class)
 *
 * RootTypeParameter            a root type
 * RootParameter     			parameter values for a specific root
 *
 * SeedTypeParamter             a single seed
 * SeedParameter				parameter values for a specific seed
 *
 */

class GrowthFunction;
class TropismFunction;
class SoilProperty;
class SignedDistanceFunction;

/**
 *	Parameter base class for specific organs
 */
class OrganParameter
{
public:

	unsigned int subType = -1;

};

/**
 * Parameter base class for all organ types
 */
class OrganTypeParameter
{
public:

	OrganTypeParameter();
	virtual ~OrganTypeParameter() { }

	std::string name = "Unnamed organ";
	unsigned int organType;
	unsigned int subType;

	virtual OrganParameter* realize() const { return new OrganParameter(); }

	virtual void readXML(FILE* fp) { };
	virtual std::string writeXML(FILE* fp) const { return ""; };
	virtual std::string toString() const { return "OrganTypeParameter base class\n";}

	/* random numbers */
	void setSeed(double seed) const { gen.seed(seed); } ///< Sets the seed of the random number generator
	double rand() const { return UD(gen); } ///< Uniformly distributed random number (0,1)
	double randn() const { return ND(gen); } ///< Normally distributed random number (0,1)

	mutable std::mt19937 gen = std::mt19937(std::chrono::system_clock::now().time_since_epoch().count());  // random stuff
	mutable std::uniform_real_distribution<double> UD = std::uniform_real_distribution<double>(0,1);
	mutable std::normal_distribution<double> ND = std::normal_distribution<double>(0,1);

};


/**
 * Parameters of a single root (created by RootParameter:realize)
 */
class RootParameter : public OrganParameter
{
public:

	RootParameter() { subType = -1; };
	RootParameter(int type, double lb, double la, const std::vector<double>& ln, double r, double a, double theta, double rlt):
		lb(lb), la(la), r(r), a(a), theta(theta), rlt(rlt), ln(ln) { subType = type;  }
	///< Constructor setting all parameters

	double getK() const; ///< Returns the exact maximal root length of this realization [cm]

	void write(std::ostream & cout) const; ///< Writes parameters for debugging
	std::string toString() const { std::stringstream ss; write(ss); return ss.str(); } ///< String representation using RootParameter::write

	/* Rootbox parameters per root */
	double lb = 0.; 		///< Basal zone [cm]
	double la = 0.;			///< Apical zone [cm];
	double r = 0.;			///< Initial growth rate [cm day-1]
	double a = 0.; 			///< Root radius [cm]
	double theta = 0.; 		///< Angle between root and parent root [rad]
	double rlt = 0.;		///< Root life time [day]
	std::vector<double> ln = std::vector<double>();    ///< Inter-lateral distances [cm]

};

/**
 * A parameter set describing a root type
 */
class RootTypeParameter : public OrganTypeParameter
{
public:

	enum TropismTypes { tt_plagio = 0, tt_gravi = 1, tt_exo = 2, tt_hydro = 3 };  ///< root tropism
	enum GrowthFunctionTypes { gft_negexp = 1, gft_linear = 2 }; // root growth function

	RootTypeParameter(); ///< default constructor
	virtual ~RootTypeParameter();

	void set(int type, double lb, double lbs, double la, double las, double ln, double lns, double nob, double nobs,
			double r, double rs, double a, double as,  double colorR, double colorG, double colorB, int tropismT, double tropismN, double tropsimS,
			double dx, const std::vector<int>& successor, const std::vector<double>& successorP, double theta, double thetas, double rlt, double rlts,
			int gf, const std::string& name); ///< sets all parameters

	virtual OrganParameter* realize() const override; ///< Creates a specific root from the root parameter set
	int getLateralType(const Vector3d& pos); ///< Choose (dice) lateral type based on root parameter set
	double getK() const { return std::max(nob-1,double(0))*ln+la+lb; }  ///< returns the mean maximal root length [cm]

	/* IO */
	virtual void read(std::istream & in); ///< reads a single root parameter set
	virtual void write(std::ostream & out) const; ///< writes a single root parameter set
	virtual void readXML(FILE* fp) override; ///< reads a single root parameter set
	virtual std::string writeXML(FILE* fp) const override; ///< writes a single root parameter set
	virtual std::string toString() const override { return writeXML(0); } ///< writes parameter to a string

	/* Rootbox parameters per root type */
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
	std::vector<int> successor;			///< Lateral types [1]
	std::vector<double> successorP; 	///< Probabiltities of lateral type to emerge (sum of values == 1) [1]

	void createTropism(SignedDistanceFunction* geom = nullptr, SoilProperty* soil = nullptr);
	void createGrowth();

	/* call back functions */
	GrowthFunction* growth;
	TropismFunction* tropism;
	SoilProperty* se; ///< scale elongation function
	SoilProperty* sa; ///< scale angle function
	SoilProperty* sbp; ///< scale branching probability function

};



/**
 * Parameters of a specific seed
 */
class SeedParameter : public OrganParameter
{
public:

	// seed
	Vector3d seedPos;   ///< Location of the seed [cm]

	// basal roots
	double firstB; 	///< Emergence of first basal root [day]
	double delayB; 	///< Time delay between the basal roots [day]
	int maxB; ///< Maximal number of basal roots [1]

};

/**
 * Seed specific parameters like planting depth, and emergence times of basal roots
 */
class SeedTypeParameter: public OrganTypeParameter
{
public:

	SeedTypeParameter();
	virtual ~SeedTypeParameter() { };

	virtual OrganParameter* realize() const;

	virtual void read(std::istream & cin); ///< Read plant parameters
	virtual void write(std::ostream & cout) const; ///< Write plant parameters
	virtual void readXML(FILE* fp); ///< Read plant parameters
	virtual std::string writeXML(FILE* fp) const; ///< Write plant parameters
	virtual std::string toString() const { std::stringstream ss; write(ss); return ss.str(); } ///< writes parameter to a string

	/* Plant parameters */
	Vector3d seedPos = Vector3d(0,0,-3);   ///< Position of the seed [cm]
	Vector3d seedPoss = Vector3d();

	/* Basal roots */
	double firstB = 1.e9; 	///< Emergence of first basal root [day]
	double firstBs = 0.;
	double delayB = 1.e9; 	///< Time delay between the basal roots [day]
	double delayBs =0.;
	int maxB = 0; 	    ///< Maximal number of basal roots [1]
	double maxBs = 0;

};



/// *****************************************************************************************************************************************************************a replication of above parts

class StemGrowthFunction;
class StemTropismFunction;


/**
 *	Parameter base class for specific organs
 */





/**
 * Parameters of a single root (created by RootParameter:realize)
 */
class StemParameter : public OrganParameter
{
public:

	StemParameter() { subType = -1; };
	StemParameter(int type, double lb, double la, const std::vector<double>& ln, double r, double a, double theta, double rlt):
		lb(lb), la(la), r(r), a(a), theta(theta), rlt(rlt), ln(ln) { subType = type;  }
	///< Constructor setting all parameters

	double getK() const; ///< Returns the exact maximal root length of this realization [cm]

	void write(std::ostream & cout) const; ///< Writes parameters for debugging
	std::string toString() const { std::stringstream ss; write(ss); return ss.str(); } ///< String representation using RootParameter::write

	/* Rootbox parameters per root */
	double lb = 0.; 		///< Basal zone [cm]
	double la = 0.;			///< Apical zone [cm];
	double r = 0.;			///< Initial growth rate [cm day-1]
	double a = 0.; 			///< Root radius [cm]
	double theta = 0.; 		///< Angle between root and parent root [rad]
	double rlt = 0.;		///< Root life time [day]
	std::vector<double> ln = std::vector<double>();    ///< Inter-lateral distances [cm]

};

/**
 * A parameter set describing a root type
 */
class StemTypeParameter : public OrganTypeParameter
{
public:

	enum StemTropismTypes { tt_plagio = 0, tt_gravi = 1, tt_exo = 2, tt_hydro = 3 };  ///< stem tropism
	enum StemGrowthFunctionTypes { gft_negexp = 1, gft_linear = 2 }; // root growth function

	StemTypeParameter(); ///< default constructor
	virtual ~StemTypeParameter();

	void set(int type, double lb, double lbs, double la, double las, double ln, double lns, double nob, double nobs,
			double r, double rs, double a, double as,  double colorR, double colorG, double colorB, int tropismT, double tropismN, double tropsimS,
			double dx, const std::vector<int>& successor, const std::vector<double>& successorP, double theta, double thetas, double rlt, double rlts,
			int gf, const std::string& name); ///< sets all parameters

	virtual OrganParameter* realize() const override; ///< Creates a specific root from the root parameter set
	int getLateralType(const Vector3d& pos); ///< Choose (dice) lateral type based on root parameter set
	double getK() const { return std::max(nob-1,double(0))*ln+la+lb; }  ///< returns the mean maximal root length [cm]

	/* IO */
	virtual void read(std::istream & in); ///< reads a single root parameter set
	virtual void write(std::ostream & out) const; ///< writes a single root parameter set
	virtual void readXML(FILE* fp) override; ///< reads a single root parameter set
	virtual std::string writeXML(FILE* fp) const override; ///< writes a single root parameter set
	virtual std::string toString() const override { return writeXML(0); } ///< writes parameter to a string

	/* Rootbox parameters per root type */
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
	std::vector<int> successor;			///< Lateral types [1]
	std::vector<double> successorP; 	///< Probabiltities of lateral type to emerge (sum of values == 1) [1]

	void createTropism(SignedDistanceFunction* geom = nullptr, SoilProperty* soil = nullptr);
	void createGrowth();

	/* call back functions */
	StemGrowthFunction* growth;
	StemTropismFunction* tropism;
	SoilProperty* se; ///< scale elongation function
	SoilProperty* sa; ///< scale angle function
	SoilProperty* sbp; ///< scale branching probability function

};



/**
 * Parameters of a specific seed
 */


/**
 * Seed specific parameters like planting depth, and emergence times of basal roots
 */
///********************************************************Another COPY PASTE. Just to make it functioning***********************************************************************************************
class LeafGrowthFunction;
class LeafTropismFunction;


/**
 *	Parameter base class for specific organs
 */





/**
 * Parameters of a single root (created by RootParameter:realize)
 */
class LeafParameter : public OrganParameter
{
public:

	LeafParameter() { subType = -1; };
	LeafParameter(int type, double lb, double la, const std::vector<double>& ln, double r, double a, double theta, double rlt):
		lb(lb), la(la), r(r), a(a), theta(theta), rlt(rlt), ln(ln) { subType = type;  }
	///< Constructor setting all parameters

	double getK() const; ///< Returns the exact maximal root length of this realization [cm]

	void write(std::ostream & cout) const; ///< Writes parameters for debugging
	std::string toString() const { std::stringstream ss; write(ss); return ss.str(); } ///< String representation using RootParameter::write

	/* Rootbox parameters per root */
	double lb = 0.; 		///< Basal zone [cm]
	double la = 0.;			///< Apical zone [cm];
	double r = 0.;			///< Initial growth rate [cm day-1]
	double a = 0.; 			///< Root radius [cm]
	double theta = 0.; 		///< Angle between root and parent root [rad]
	double rlt = 0.;		///< Root life time [day]
	std::vector<double> ln = std::vector<double>();    ///< Inter-lateral distances [cm]

};

/**
 * A parameter set describing a root type
 */
class LeafTypeParameter : public OrganTypeParameter
{
public:

	enum LeafTropismTypes { tt_plagio = 0, tt_gravi = 1, tt_exo = 2, tt_hydro = 3 };  ///< stem tropism
	enum LeafGrowthFunctionTypes { gft_negexp = 1, gft_linear = 2 }; // root growth function

	LeafTypeParameter(); ///< default constructor
	virtual ~LeafTypeParameter();

	void set(int type, double lb, double lbs, double la, double las, double ln, double lns, double nob, double nobs,
			double r, double rs, double a, double as,  double colorR, double colorG, double colorB, int tropismT, double tropismN, double tropsimS,
			double dx, const std::vector<int>& successor, const std::vector<double>& successorP, double theta, double thetas, double rlt, double rlts,
			int gf, const std::string& name); ///< sets all parameters

	virtual OrganParameter* realize() const override; ///< Creates a specific root from the root parameter set
	int getLateralType(const Vector3d& pos); ///< Choose (dice) lateral type based on root parameter set
	double getK() const { return std::max(nob-1,double(0))*ln+la+lb; }  ///< returns the mean maximal root length [cm]

	/* IO */
	virtual void read(std::istream & in); ///< reads a single root parameter set
	virtual void write(std::ostream & out) const; ///< writes a single root parameter set
	virtual void readXML(FILE* fp) override; ///< reads a single root parameter set
	virtual std::string writeXML(FILE* fp) const override; ///< writes a single root parameter set
	virtual std::string toString() const override { return writeXML(0); } ///< writes parameter to a string

	/* Rootbox parameters per root type */
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
	std::vector<int> successor;			///< Lateral types [1]
	std::vector<double> successorP; 	///< Probabiltities of lateral type to emerge (sum of values == 1) [1]

	void createTropism(SignedDistanceFunction* geom = nullptr, SoilProperty* soil = nullptr);
	void createGrowth();

	/* call back functions */
	LeafGrowthFunction* growth;
	LeafTropismFunction* tropism;
	SoilProperty* se; ///< scale elongation function
	SoilProperty* sa; ///< scale angle function
	SoilProperty* sbp; ///< scale branching probability function

};



/**
 * Parameters of a specific seed
 */


/**
 * Seed specific parameters like planting depth, and emergence times of basal roots
 */
#endif
