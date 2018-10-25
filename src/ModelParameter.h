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
#include "tinyxml2.h"
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
class SoilLookUp;
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

	virtual void readXML(const tinyxml2::XMLElement* ele) { };
	virtual std::string writeXML(FILE* fp) const { return ""; };
	virtual std::string toString() const { return "OrganTypeParameter base class\n";}
     void getAttribute(const tinyxml2::XMLElement* ele2,const char* attr_name, const char* para_name, double &attr, double &deviation ); //alot of overloading, try template later..
     void getAttribute(const tinyxml2::XMLElement* ele2,const char* attr_name, const char* para_name, int &attr, double &deviation );
     void getAttribute(const tinyxml2::XMLElement* ele2,const char* attr_name, const char* para_name, double &attr, double &deviation, int &functiontype );
     void getAttribute(const tinyxml2::XMLElement* ele2,const char* attr_name, const char* para_name, std::vector<int> &successor, std::vector<double> &successorP );
     void getAttribute(const tinyxml2::XMLElement* ele2, const char* attr_name, const char* para_name, double &attr) ;
     void getAttribute(const tinyxml2::XMLElement* ele2, const char* attr_name, const char* para_name, int &attr) ;
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
	virtual void readXML(const tinyxml2::XMLElement* ele) override; ///< reads a single root parameter set
	virtual std::string writeXML(FILE* fp) const override; ///< writes a single root parameter set
	virtual std::string toString() const override { return writeXML(0); } ///< writes parameter to a string

	/* Rootbox parameters per root type */
	double lb; 	 	///< Basal zone [cm]
	double lbs;  	///< Standard deviation basal zone [cm]
	double la;		///< Apical zone [cm];
	double las;		///< Standard deviation apical zone [cm];
	double ln; 		///< Inter-lateral distance [cm]
	double lns;  	///< Standard deviation inter-lateral distance [cm]
	double k;
	double ks;
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

	void createTropism(SignedDistanceFunction* geom = nullptr, SoilLookUp* soil = nullptr);
	void createGrowth();

	/* call back functions */
	GrowthFunction* growth;
	TropismFunction* tropism ;
	SoilLookUp* se; ///< scale elongation function
	SoilLookUp* sa; ///< scale angle function
	SoilLookUp* sbp; ///< scale branching probability function
	
};



/**
 * Parameters of a specific seed
 */
class SeedParameter : public OrganParameter
{
public:

	// seed
	Vector3d seedPos;   ///< Location of the seed [cm]
	double firstB ; 	///< Emergence of first basal root [day]
	double firstBs;
	double delayB ; 	///< Time delay between the basal roots [day]
	double delayBs;
	int maxB ; 	    ///< Maximal number of basal roots [1]
	double maxBs ;
	double nC ;
	double nz ;
	double simtime ;
	double firstSB ;
	double firstSBs ;
	double delaySB ;
	double delaySBs ;
	double delayRC ;
	double delayRCs ;

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

	virtual void readXML(const tinyxml2::XMLElement* ele) override;
	virtual void read(std::istream & cin); ///< Read plant parameters
	virtual void write(std::ostream & cout) const; ///< Write plant parameters
	virtual std::string writeXML(FILE* fp) const; ///< Write plant parameters
	virtual std::string toString() const { std::stringstream ss; write(ss); return ss.str(); } ///< writes parameter to a string

	/* Plant parameters */
	Vector3d seedPos = Vector3d(0,0,-3);   ///< Position of the seed [cm]
	Vector3d seedPoss = Vector3d();
    double plantingdepth = -seedPos.z;
	/* Basal roots */
	double firstB = 1.e9; 	///< Emergence of first basal root [day]
	double firstBs = 0.;
	double delayB = 1.e9; 	///< Time delay between the basal roots [day]
	double delayBs =0.;
	int maxB = 0; 	    ///< Maximal number of basal roots [1]
	double maxBs = 0;
	double nC = 0;
	double nz = 0;
	double simtime = 60;
	double firstSB = 0;
	double firstSBs = 0;
	double delaySB = 0;
	double delaySBs = 0;
	double delayRC = 0;
	double delayRCs = 0;

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
	StemParameter(int type, double lb, double la, const std::vector<double>& ln, double r, double a, double theta, double rlt, int lnf):
		lb(lb), la(la), r(r), a(a), theta(theta), rlt(rlt), lnf(lnf), ln(ln) { subType = type;  }
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
	int lnf = 0;
	std::vector<double> ln = std::vector<double>();    ///< Inter-lateral distances [cm]

};

/**
 * A parameter set describing a root type
 */
class StemTypeParameter : public OrganTypeParameter
{
public:

	enum StemTropismTypes { tt_plagio = 0, tt_gravi = 1, tt_exo = 2, tt_hydro = 3, tt_antigravi = 4 };  ///< stem tropism
	enum StemGrowthFunctionTypes { gft_negexp = 1, gft_linear = 2 }; // root growth function

	StemTypeParameter(); ///< default constructor
	virtual ~StemTypeParameter();

	void set(int type, double lb, double lbs, double la, double las, double ln, double lns, int lnf, double nob, double nobs,
			double r, double rs, double a, double as,  double colorR, double colorG, double colorB, int tropismT, double tropismN, double tropsimS,
			double dx, const std::vector<int>& successor, const std::vector<double>& successorP, double theta, double thetas, double rlt, double rlts,
			int gf, const std::string& name); ///< sets all parameters

	virtual OrganParameter* realize() const override; ///< Creates a specific root from the root parameter set
	int getLateralType(const Vector3d& pos); ///< Choose (dice) lateral type based on root parameter set
	double getK() const { return std::max(nob-1,double(0))*ln+la+lb; }  ///< returns the mean maximal root length [cm]

	/* IO */
	virtual void read(std::istream & in); ///< reads a single root parameter set
	virtual void write(std::ostream & out) const; ///< writes a single root parameter set
	virtual void readXML(const tinyxml2::XMLElement* ele) override; ///< reads a single root parameter set
	virtual std::string writeXML(FILE* fp) const override; ///< writes a single root parameter set
	virtual std::string toString() const override { return writeXML(0); } ///< writes parameter to a string

	/* Rootbox parameters per root type */
	double lb; 	 	///< Basal zone [cm]
	double lbs;  	///< Standard deviation basal zone [cm]
	double la;		///< Apical zone [cm];
	double las;		///< Standard deviation apical zone [cm];
	double ln; 		///< Inter-lateral distance [cm]
	double lns;  	///< Standard deviation inter-lateral distance [cm]
	int lnf;
			double k;
	double ks;
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

	void createTropism(SignedDistanceFunction* geom = nullptr, SoilLookUp* soil = nullptr);
	void createGrowth();

	/* call back functions */
	StemGrowthFunction* growth;
	StemTropismFunction* tropism;
	SoilLookUp* se; ///< scale elongation function
	SoilLookUp* sa; ///< scale angle function
	SoilLookUp* sbp; ///< scale branching probability function

 /**
 * Parameters from Graal model
 */
//    ///General parameters
//    double deltat = 80 ; ///Simulation duration (d)
//    double a_lPAR = -0.10 ; ///Parameters for daily average PAR decrease during the vegetative period
//    double b_lPAR = 0.005 ; ///Parameters for daily average PAR decrease during the vegetative period
//    double a_hPAR = -0.14 ; ///Parameters for daily average PAR decrease during the vegetative period
//    double b_hPAR = 0.003 ; ///Parameters for daily average PAR decrease during the vegetative period
//    double g_organ_growth_effi = 0.73 ; ///Organ growth efficiency (dl)
//    double m_organ_mainte_effi = 0.0032 ; ///Organ maintenance respiration rate (g CO2 g -1 DM)
//    double M_si = 0.29 ; /// Seed initial dry mass (g DM)
//    double Q_10 = 2 ; /// Temperature coefficient (dl)
//    double R_r = 0.05 ; /// Reserve supply rate (d^-1)
//	double R_s = 0.15 ; /// Seed supply rate (d^-1)
//	double T_ref = 20 ; /// Reference temperature for maintenance respiration
//
//	///Shoot development
//	double a_B_i = 2.3 ; ///Parameters for final internode diameter (cm)
//	double b_B_i = 3.3 ; ///Parameters for final internode diameter (cm)
//	double c_B_i = -0.17 ; ///Parameters for final internode diameter (cm)
//	double n_B_i = 7 ; ///Parameters for final internode diameter (dl)
//
//	double a_Die = -5.16 ; ///Parameters for the time between primordia initiation and beginning of leaf elongation
//	double b_Die = 1.94 ; ///Parameters for the time between primordia initiation and beginning of leaf elongation
//	double n_Die = 3.65 ; ///Parameters for the time between primordia initiation and beginning of leaf elongation
//
//	double a_L_i = 0.202 ; ///Parameters for final internode length (cm)
//	double b_L_i = 1.233 ; ///Parameters for final internode length (cm)
//	double c_L_i = -0.040 ; ///Parameters for final internode length (cm)
//
//    double a_L_S = 3.077 ; ///Parameters for final sheath length (cm)
//    double b_L_S = 2.048 ; ///Parameters for final sheath length (cm)
//    double c_L_S = -0.569 ; ///Parameters for final sheath length (cm)
//    double n_L_S = 6 ; ///Parameters for final sheath length (dl)
//
//    double a_ls = 0.8 ; ///Parameters for lamina shape (dl)
//    double b_ls = 1.3 ; ///Parameters for lamina shape (dl)
//    double c_ls = -2.1 ; ///Parameters for lamina shape (dl)
//
//    double a_n_M = 5.93 ; ///parameters for final lamina length (dl)
//    double b_n_M = 0.33 ; ///parameters for final lamina length (dl)
//    double a_ala = -10.61 ; ///parameters for final lamina length (dl)
//    double b_ala = 0.25 ; ///parameters for final lamina length (dl)
//    double a_bla = -5.99 ; ///parameters for final lamina length (dl)
//    double b_bla = 0.27 ; ///parameters for final lamina length (dl)
//
//    double a_pp = -1.06 ; ///Parameters for the phyllochron number-plastochron number relationship
//    double b_pp = 0.54 ; ///Parameters for the phyllochron number-plastochron number relationship
//
//    double a_R_p = -0.00065 ; ///Parameters for phytomer initiation (phytomer \B0C d^-1)
//    double b_R_p = -0.0138 ; ///Parameters for phytomer initiation (phytomer \B0C d^-1)
//    double c_R_p = 0.00372 ; ///Parameters for phytomer initiation (phytomer \B0C d^-1)
//    double d_R_p = -0.000072 ; ///Parameters for phytomer initiation (phytomer \B0C d^-1)
//
//    double alpha = 0.04 ; ///Photosynthetic efficiency (miu mol CO_2 miu mol^-1)
//    double b_LSI = 0.0003; ///Lamina and sheath dry mass per unit area increment
//    double IMV_max = 0.140; ///Maximal internode dry mass per unit volume (g DM cm^-2)
//    double IMV_min = 0.005; ///Minimal internode dry mass per unit volume (g DM cm^-2)
//    double k_A = 0.75; ///Allometry coefficient for final lamina area (dl)
//    double k_W = 0.106; ///Allometry coefficient for final lamina area (dl)
//    double L_lm = 90 ; /// Minimal length of the longest internode (cm)
//    double L_im = 23 ; ///Final length of the longest internode (cm)
//    double LMA_max = 0.0075 ; ///Minimal lamina dry mass per unit area for lamina N_p (g DM cm^-2)
//    double LMA_min = 0.0010 ; ///Minimal lamina and sheath dry mass perunit area (g DM cm^-2)
//    double N_p = 16 ;///Total number of phytomers
//    double n_0 = 5 ; /// last internode that does not elongate
//    double P_max = 30 ; ///Photosynthetic rate at saturating light (miu mol CO_2 m^-2 s^-1)
//    double SMA_max = 0.0075 ; ///Maximal sheath dry mass per unit area for sheath N_p (g DM cm^-2)
//    double SMA_min = 0.0010 ; ///Minimal lamina and sheath dry mass per unit area (g DM cm^-2)
//    double T_base = 8 ; /// Base temperature (\B0C)
//    double T_max = 50 ; /// Maximal temperature (\B0C)
//    double T_opt = 31 ; /// Optimal temperature (\B0C)
//    double vstar_c = 0.41 ; /// Multiplicative factor for standard internode elongation rate (dl)
//    double vstar_1 = 0.564 ; ///standard leaf (lamina and sheath) elongation rate (cm\B0C d^-1)


};

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

	enum LeafTropismTypes { tt_plagio = 0, tt_gravi = 1, tt_exo = 2, tt_hydro = 3 , tt_antigravi = 4 };  ///< stem tropism
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
	virtual void readXML(const tinyxml2::XMLElement* ele) override; ///< reads a single root parameter set
	virtual std::string writeXML(FILE* fp) const override; ///< writes a single root parameter set
	virtual std::string toString() const override { return writeXML(0); } ///< writes parameter to a string

	/* Rootbox parameters per root type */
	double lb; 	 	///< Basal zone [cm]
	double lbs;  	///< Standard deviation basal zone [cm]
	double la;		///< Apical zone [cm];
	double las;		///< Standard deviation apical zone [cm];
	double ln; 		///< Inter-lateral distance [cm]
	double lns;  	///< Standard deviation inter-lateral distance [cm]
			double k;
	double ks;
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

	void createTropism(SignedDistanceFunction* geom = nullptr, SoilLookUp* soil = nullptr);
	void createGrowth();

	/* call back functions */
	LeafGrowthFunction* growth;
	LeafTropismFunction* tropism;
	SoilLookUp* se; ///< scale elongation function
	SoilLookUp* sa; ///< scale angle function
	SoilLookUp* sbp; ///< scale branching probability function

};

/**
 * override the PrintSpace function
 */

//void tinyxml2::XMLPrinter::PrintSpace( const char * indent, int depth ) //print space TODO
//{
//    if (indent == "") indent = " ";
//    if (depth == 0) depth = 4;  // 0 means default
//    for ( int i=0; i<depth; ++i ) {
//        Print( indent );  // depth here goes in single spaces
//    }
//}

/**
 * Parameters of a specific seed
 */



#endif
