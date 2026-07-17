// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ORGANPARAMETER_H_
#define ORGANPARAMETER_H_

#include "mymath.h"
#include "tinyxml2.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

/**
 * @file organparameter.h
 * @brief Defines OrganSpecificParameter and OrganRandomParameter, the base parameter classes for all organ types.
 *
 * OrganSpecificParameter holds the realised (per-organ-instance) parameters drawn by OrganRandomParameter::realize().
 * OrganRandomParameter holds a stochastic parameter set for a single organ sub-type and provides XML I/O,
 * class introspection via named hashmaps, successor rule management, and factory methods.
 */

namespace CPlantBox {

class Organism; // forward declaration
class GrowthFunction;
class ExponentialGrowth;
class Tropism;

/**
 * @brief Realised (instance-level) parameters for a single organ.
 *
 * One OrganSpecificParameter object is created per organ instance by OrganRandomParameter::realize().
 * Derived classes (e.g. RootSpecificParameter) extend this with additional organ-specific fields.
 */
class OrganSpecificParameter {
  public:
  
    OrganSpecificParameter(int t, double a) : subType(t), a(a) {} ///< @param t sub-type index, @param a radius [cm]

    virtual ~OrganSpecificParameter() {}
    int subType = -1;                     ///< Sub-type index of the organ; matches OrganRandomParameter::subType
    double a = 0.;                        ///< Radius of the organ [cm]
    virtual std::string toString() const; ///< Returns a short human-readable summary for debugging
};

/**
 * @brief Stochastic parameter set for a single organ sub-type.
 *
 * Contains all parameters needed to stochastically generate organ instances.
 * OrganSpecificParameter objects are created by realize(), which samples random
 * variables from their distributions.
 *
 * Parameters are registered in named hashmaps (dparam / iparam) via bindParameters(),
 * so that getParameter(), toString(), readXML(), and writeXML() work generically for
 * scalar int and double fields in derived classes.  Non-scalar or non-standard parameter
 * types require method overrides (see e.g. RootRandomParameter).
 *
 * The factory functions copy() and realize() must be overridden for each specialisation.
 *
 * @see OrganSpecificParameter, RootRandomParameter, StemRandomParameter, LeafRandomParameter
 */
class OrganRandomParameter {
  public:
    OrganRandomParameter(std::shared_ptr<Organism> plant); ///< Constructor; binds all scalar parameters to the introspection maps
    virtual ~OrganRandomParameter() {};

    virtual std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant); ///< Deep-copies this parameter set into @p plant; ownership of the result is transferred to caller

    virtual std::shared_ptr<OrganSpecificParameter> realize(); ///< Samples random variates and returns a realised OrganSpecificParameter; ownership transferred to caller

    virtual double getParameter(std::string name) const; ///< Returns a named scalar parameter; append "_dev" for its deviation, "_mean" for its mean; returns NaN if unknown

    virtual std::string toString(bool verbose = true) const; ///< Returns a human-readable summary; verbose=true prints a full parameter table, false prints a one-liner

    virtual void readXML(tinyxml2::XMLElement *element, bool verbose); ///< Reads parameters from @c <parameter> child tags; missing tags leave defaults unchanged; called by Organism::readParameters
    void readSuccessor(tinyxml2::XMLElement *p, bool verbose); ///< Parses a @c <parameter name="successor"> element and updates the successor rule vectors
    void readXML(std::string name, bool verbose); ///< Reads parameters from an XML file (convenience overload for testing/debugging)
    virtual tinyxml2::XMLElement *writeXML(tinyxml2::XMLDocument &doc, bool comments = true) const; ///< Serialises all registered parameters and successor rules to an XML element; not exposed to Python
    void writeXML(std::string name) const; ///< Writes this parameter set to an XML file (convenience overload for testing/debugging)

    int getLateralType(const Vector3d &pos, int ruleId); ///< Stochastically selects a lateral sub-type index using successorP[ruleId]; returns -1 if no lateral emerges

    virtual void bindParameters(); ///< Registers all scalar parameters in the introspection maps; called by constructor; derived classes must call base first
    void bindParameter(std::string name, int *i, std::string descr = "", double *dev = nullptr); ///< Registers an int member under @p name in the introspection maps
    void bindParameter(std::string name, double *d, std::string descr = "", double *dev = nullptr); ///< Registers a double member under @p name in the introspection maps

    std::string name = "organ"; ///< Name of the organ sub-type (e.g. "taproot", "small lateral")
    int organType = 0;          ///< Organ type identifier (0 = unspecified, 1 = seed, 2 = root, 3 = stem, 4 = leaf)
    int subType = 0;            ///< Unique sub-type identifier within the organ type
    double a = 0.1;             ///< Mean radius of the organ [cm]
    double as = 0.;             ///< Standard deviation of the radius [cm]
    double dx = 0.25;           ///< Maximal axial segment length (spatial resolution) [cm]
    double dxMin = 1e-6;        ///< Minimal axial segment length; segments shorter than this are skipped to avoid NaN tip directions [cm]
    double ldelay = -1.;        ///< Mean lateral emergence delay [day]; used by RootDelay and Organism::delayDefinition != Organism::dd_distance
                                ///< @see RootDelay, RootSystem::initializeDB
    double ldelays = 0.;        ///< Standard deviation of the lateral emergence delay [day]
	int multDelay = 1; ///< by how much multiply the delay between each consecutive lateral, see @Stem::getLatGrowthDelay

    /// @brief Position-based activation flag for each successor rule.
    ///
    /// successorWhere[ruleId] is a vector of linking-node positions at which the rule applies.
    /// A value of +1 means "apply", -1 means "do not apply".  An empty vector (default) means
    /// the rule applies at every linking node.  double is used to distinguish -0 from 0.
    std::vector<std::vector<double>> successorWhere = std::vector<std::vector<double>>(0, std::vector<double>(0, 0));

    std::vector<std::vector<int>> successorOT = std::vector<std::vector<int>>(0, std::vector<int>(0, 0)); ///< Organ types of candidate laterals for each successor rule [1]
    std::vector<std::vector<int>> successorST = std::vector<std::vector<int>>(0, std::vector<int>(0, 0)); ///< Sub-types of candidate laterals for each successor rule [1]
    std::vector<std::vector<double>> successorP = std::vector<std::vector<double>>(0, std::vector<double>(0, 0)); ///< Emergence probabilities for each candidate lateral (values must sum to ≤ 1) [1]
    std::vector<int> successorNo = std::vector<int>(0);                 ///< Number of laterals to create per successor rule [1]

    std::weak_ptr<Organism> plant;        ///< Weak reference to the owning organism
    std::shared_ptr<Tropism> f_tf;        ///< Tropism function (initialised to a default Tropism in the constructor)
    std::shared_ptr<GrowthFunction> f_gf; ///< Growth function (initialised to ExponentialGrowth in the constructor)

  protected:
    /* class introspection */
    std::map<std::string, double *> dparam;         ///< Named double parameters registered via bindParameter()
    std::map<std::string, int *> iparam;            ///< Named integer parameters registered via bindParameter()
    std::map<std::string, double *> param_sd;       ///< Deviations (e.g. standard deviations) for registered parameters
    std::map<std::string, std::string> description; ///< Human-readable descriptions for registered parameters

    std::vector<int>    string2vector(const char *xmlInput, int defaultVal);    ///< Parses a comma-separated XML attribute string into a vector of int
    std::vector<double> string2vector(const char *xmlInput, double defaultVal); ///< Parses a comma-separated XML attribute string into a vector of double
    std::string vector2string(std::vector<int> vec) const;                      ///< Converts an int vector to a comma-separated string for XML output
    std::string vector2string(std::vector<double> vec) const;                   ///< Converts a double vector to a comma-separated string for XML output

    /// Queries a string XML attribute under one of several candidate @p keyNames and appends parsed values to @p vToFill; falls back to @p defaultVal if not found and @p replaceByDefault is true
    template <class IntOrDouble>
    void cpb_queryStringAttribute(std::vector<std::string> keyNames, IntOrDouble defaultVal, int sizeVector, bool replaceByDefault,
                                  std::vector<IntOrDouble> &vToFill, tinyxml2::XMLElement *key);
};

} // namespace CPlantBox

#endif
