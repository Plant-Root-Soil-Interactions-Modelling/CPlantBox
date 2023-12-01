// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ORGANPARAMETER_H_
#define ORGANPARAMETER_H_

#include <string>
#include <map>
#include <memory>
#include <iostream>
#include <vector>
#include <string>

#include "tinyxml2.h"
#include "tropism.h"

/**
 * This file describes the classes OrganSpecificParameter and OrganRandomParameter.
 * OrganSpecificParameter are drawn from the OrganRandomParameter class
 */

namespace CPlantBox {

class Organism; // forward declaration
class GrowthFunction;
class ExponentialGrowth;

/**
 * Parameters for a specific organ
 */
class OrganSpecificParameter {
public:

    OrganSpecificParameter(int t, double a): subType(t), a(a)  { }

    virtual ~OrganSpecificParameter() { }
	int created_linking_node = 0;
    int subType = -1; ///< sub type of the organ
    double a = 0.; ///< radius of the organ [cm]
    virtual std::string toString() const; ///< quick info for debugging

};

/**
 * Contains a parameter set describing as single sub type of an organ,
 * specific parameters are then created with realize().
 *
 * Organizes parameters in hash-maps for scalar double and scalar int values.
 * For this reason derived classes getParameter(), toString(), readXML(), and writeXML() should work out of the box.
 * For other parameter types the methods must be overwritten, see e.g. RootRandomParameter.
 *
 * The factory function copy() has to be overwritten for each specialization.
 */
class OrganRandomParameter
{
public:

    OrganRandomParameter(std::shared_ptr<Organism> plant); ///< default constructor
    virtual ~OrganRandomParameter() { };

    virtual std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant); ///< copies the root type parameter into a new plant

    virtual std::shared_ptr<OrganSpecificParameter> realize(); ///< creates a specific organ from the root parameter set

    virtual double getParameter(std::string name) const; // get a scalar parameter

    virtual std::string toString(bool verbose = true) const; ///< info for debugging

    virtual void readXML(tinyxml2::XMLElement* element, bool verbose); ///< reads a single sub type organ parameter set
	void readSuccessor(tinyxml2::XMLElement* p, bool verbose);
	void readXML(std::string name, bool verbose); ///< reads a single sub type organ parameter set
    virtual tinyxml2::XMLElement* writeXML(tinyxml2::XMLDocument& doc, bool comments = true) const; ///< writes a organ root parameter set
    void writeXML(std::string name) const; ///< writes a organ root parameter set

	int getLateralType(const Vector3d& pos, int ruleId); ///< Choose (dice) lateral type based on stem parameter set

	virtual void bindParameters(); ///<sets up class introspection
    void bindParameter(std::string name, int* i, std::string descr = "", double* dev = nullptr); ///< binds integer to parameter name
    void bindParameter(std::string name, double* d, std::string descr = "", double* dev = nullptr); ///< binds double to parameter name

    std::string name = "organ";
    int organType = 0;
    int subType = 0;
    double a = 0.1; 		///< Root radius [cm]
    double as = 0.; 		///< Standard deviation root radius [cm]
    double dx = 0.25; 		///< Maximal segment size [cm]
	double dxMin = 1e-6; 	///< threshold value, smaller segments will be skipped (otherwise stem tip direction can become NaN)
	double ldelay = -1.;     ///< Lateral emergence delay [day], used by RootDelay, @see RootDelay, RootSystem::initializeDB or if Organism->delayDefinition != Organism::dd_distance
    double ldelays = 0.;    ///< Standard deviation of lateral emergence delay [day]

    std::weak_ptr<Organism> plant;
	std::shared_ptr<Tropism> f_tf;  ///< tropism function (defined in constructor as new Tropism(plant))
    std::shared_ptr<GrowthFunction> f_gf;
	std::vector<std::vector<int> > successorST = std::vector<std::vector<int>>(0, std::vector<int> (0, 0));			///< Lateral types [1]
	std::vector<std::vector<double>> successorP = std::vector<std::vector<double>>(0, std::vector<double> (0, 0));  	///< Probabilities of lateral type to emerge (sum of values == 1) [1]
    std::vector<int>  successorNo = std::vector<int>(0);			///< Lateral types [1]
	//need to use double to distiguish between -0 and 0
	//default: vector empty == rule implemented everywhere
    std::vector<std::vector<double> > successorWhere = std::vector<std::vector<double>>(0, std::vector<double> (0, 0));  	///< Where should rule be implemented [1] or not [-1]
    std::vector<std::vector<int> > successorOT = std::vector<std::vector<int>>(0, std::vector<int> (0, 0));			///< Lateral types [1]


protected:

    /* class introspection */
    std::map<std::string, double*> dparam; ///< Parameters with type double that can be read and written
    std::map<std::string, int*> iparam; ///< Parameters with type double that can be read and written
    std::map<std::string, double*> param_sd; ///< Deviations of parameters
    std::map<std::string, std::string> description; ///< Parameter descriptions

    std::string vector2string(std::vector<int> vec) const;
    std::string vector2string(std::vector<double> vec) const;
    std::vector<int> string2vector(const char* xmlInput, int defaultVal);
    std::vector<double> string2vector(const char* xmlInput, double defaultVal);

    template <class IntOrDouble>
    void cpb_queryStringAttribute(std::vector<std::string> keyNames,IntOrDouble defaultVal,int sizeVector,
                                    bool replaceByDefault,
                                    std::vector<IntOrDouble> & vToFill, tinyxml2::XMLElement* key);


};

} // namespace

#endif
