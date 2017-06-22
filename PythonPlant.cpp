#ifndef PY_ROOTBOX_H_
#define PY_ROOTBOX_H_

// copy paste for daniel
// 1. export LD_LIBRARY_PATH=~/boost_1_62_0/stage/lib
// 2.  g++ -std=c++11 -O3 -fpic -shared -o py_rootbox.so -Wl,-soname,"py_rootbox.so" PythonRootSystem.cpp -I/usr/include/python3.5 -lboost_python-py35 Debug/ModelParameter.o Debug/Root.o Debug/RootSystem.o Debug/analysis.o Debug/sdf.o Debug/tropism.o Debug/examples/Exudation/gauss_legendre.o
// 2b  g++ -std=c++11 -O3 -fpic -shared -o py_rootbox.so -Wl,-soname,"py_rootbox.so" PythonRootSystem.cpp -I/usr/include/python3.5 -lboost_python-py35 Debug/examples/Exudation/gauss_legendre.o

/**
 *  A Python module for CRootbox based on boost.python
 *
 *  build a shared library from this file
 *  put comment to line 16 to ignore this file
 */
// #define PYTHON_WRAPPER // UNCOMMENT TO BUILD SHARED LIBRARY

#ifdef PYTHON_WRAPPER

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/call.hpp>

#include "mymath.h"
#include "sdf.h"
#include "Plant.h"
#include "analysis.h"
#include "examples/Exudation/example_exudation.h"

using namespace boost::python;

/*
 * Functions overloading (by hand, there are also macros available)
 *
 * Tutorial example:
 * bool    (X::*fx1)(int)              = &X::f;
 * bool    (X::*fx2)(int, double)      = &X::f;
 * bool    (X::*fx3)(int, double, char)= &X::f;
 * int     (X::*fx4)(int, int, int)    = &X::f;
 *
 */
Vector3d (Vector3d::*times1)(const double) const = &Vector3d::times;
double (Vector3d::*times2)(const Vector3d&) const = &Vector3d::times;

void (Matrix3d::*times3)(const Matrix3d&) = &Matrix3d::times;
Vector3d (Matrix3d::*times4)(const Vector3d&) const = &Matrix3d::times;

std::string (SignedDistanceFunction::*writePVPScript)() const = &SignedDistanceFunction::writePVPScript; // because of default value

void (RootSystem::*simulate1)(double dt, bool silence) = &RootSystem::simulate;
void (RootSystem::*simulate2)() = &RootSystem::simulate;

void (SegmentAnalyser::*addSegments1)(const RootSystem& rs) = &SegmentAnalyser::addSegments;
void (SegmentAnalyser::*addSegments2)(const SegmentAnalyser& a) = &SegmentAnalyser::addSegments;
void (SegmentAnalyser::*filter1)(int st, double min, double max) = &SegmentAnalyser::filter;
void (SegmentAnalyser::*filter2)(int st, double value) = &SegmentAnalyser::filter;
double (SegmentAnalyser::*getSummed1)(int st) const = &SegmentAnalyser::getSummed;
double (SegmentAnalyser::*getSummed2)(int st, SignedDistanceFunction* geometry) const = &SegmentAnalyser::getSummed;
std::vector<double> (SegmentAnalyser::*distribution_1)(int st, double top, double bot, int n, bool exact) const = &SegmentAnalyser::distribution;
std::vector<SegmentAnalyser> (SegmentAnalyser::*distribution_2)(double top, double bot, int n) const = &SegmentAnalyser::distribution;
std::vector<std::vector<double>> (SegmentAnalyser::*distribution2_1)(int st, double top, double bot, double left, double right, int n, int m, bool exact) const = &SegmentAnalyser::distribution2;
std::vector<std::vector<SegmentAnalyser>> (SegmentAnalyser::*distribution2_2)(double top, double bot, double left, double right, int n, int m) const = &SegmentAnalyser::distribution2;
SegmentAnalyser (SegmentAnalyser::*cut1)(const SDF_HalfPlane& plane) const = &SegmentAnalyser::cut;
// static Vector3d (SegmentAnalyser::*cut2)(Vector3d in, Vector3d out, SignedDistanceFunction* geometry) = &SegmentAnalyser::cut; // not working, dont know why, problem with static?


/**
 * Default arguments: no idea how to do it by hand,  magic everywhere...
 */
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(initialize_overloads,initialize,0,2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(openFile_overloads,openFile,1,2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(simulate1_overloads,simulate,1,2);

/**
 * Virtual functions (not sure if needed, or only if we derive classes from it in python?), not working...
 *
 * it seems a bit tricky to make polymorphism work, we have to wrap the classes
 *
 * Tutorial example:
 * struct BaseWrap : Base, wrapper<Base>
 * {
 *    int f()
 *    {
 *        if (override f = this->get_override("f"))
 *            return f(); // *note*
 *        return Base::f();
 *    }
 *
 *    int default_f() { return this->Base::f(); }
 * };
 */
//struct SignedDistanceFunction_Wrap : SignedDistanceFunction, wrapper<SignedDistanceFunction>
//{
//	double getDist(const Vector3d& v) const {
//		if (override getDist = this->get_override("getDist"))
//			return getDist(v);
//		return SignedDistanceFunction::getDist(v);
//	}
//	double default_getDist(const Vector3d& v) { return this->SignedDistanceFunction::getDist(v); }
//};
//	class_<SignedDistanceFunction_Wrap, boost::noncopyable>("SignedDistanceFunction")
//	    .def("getDist", &SignedDistanceFunction_Wrap::getDist, &SignedDistanceFunction_Wrap::default_getDist)
//	; // TODO how does polymorphism work... (everything works fine, dont ask why)
// tricky booom boom (?)


struct SoilProperty_Wrap : SoilProperty, wrapper<SoilProperty> {

    double getValue(const Vector3d& pos, const Root* root = nullptr) const {

    	return this->get_override("getValue")(pos, root);
    }

    std::string toString() const {
    	return this->get_override("toString")();
    }

};

/**
 * Expose classes to Python module
 */
BOOST_PYTHON_MODULE(py_rootbox)
{
    /*
     * general
     */
    class_<std::vector<double>>("std_vector_double_")
        .def(vector_indexing_suite<std::vector<double>>() )
	;
    class_<std::vector<int>>("std_vector_int_")
        .def(vector_indexing_suite<std::vector<int>>() )
	;
	/*
	 * mymath.h
	 */
  	class_<Vector2i>("Vector2i", init<>())
			.def(init<int,int>())
			.def(init<Vector2i&>())
			.def_readwrite("x",&Vector2i::x)
			.def_readwrite("y",&Vector2i::y)
			.def("__str__",&Vector2i::toString)
	;
    class_<std::vector<Vector2i>>("std_vector_Vector2i_")
        .def(vector_indexing_suite<std::vector<Vector2i>>() )
	;
	class_<Vector3d>("Vector3d", init<>())
			.def(init<double,double,double>())
			.def(init<Vector3d&>())
			.def_readwrite("x",&Vector3d::x)
			.def_readwrite("y",&Vector3d::y)
			.def_readwrite("z",&Vector3d::z)
			.def("normalize",&Vector3d::normalize)
			.def("times",times1)
			.def("times",times2)
			.def("length",&Vector3d::length)
			.def("plus",&Vector3d::plus)
			.def("minus",&Vector3d::minus)
			.def("cross",&Vector3d::cross)
			.def("__str__",&Vector3d::toString)
			.def("__rep__",&Vector3d::toString)
	;
    class_<std::vector<Vector3d>>("std_vector_Vector3d_")
        .def(vector_indexing_suite<std::vector<Vector3d>>() )
	;
	class_<Matrix3d>("Matrix3d", init<>())
			.def(init<double,double,double,double,double,double,double,double,double>())
			.def(init<Vector3d&, Vector3d&, Vector3d&>())
			.def(init<Matrix3d&>())
			.def_readwrite("r0",&Matrix3d::r0)
			.def_readwrite("r1",&Matrix3d::r1)
			.def_readwrite("r2",&Matrix3d::r2)
			.def("rotX",&Matrix3d::rotX)
			.def("rotY",&Matrix3d::rotY)
			.def("rotZ",&Matrix3d::rotZ)
			.def("ons",&Matrix3d::ons)
			.def("det",&Matrix3d::det)
			.def("inverse",&Matrix3d::inverse)
			.def("column",&Matrix3d::column)
			.def("row",&Matrix3d::row)
			.def("times",times3)
			.def("times",times4)
			.def("__str__",&Matrix3d::toString)
			.def("__rep__",&Matrix3d::toString)
	;
	/*
	 * sdf.h
	 */
	class_<SignedDistanceFunction>("SignedDistanceFunction")
			.def("getDist",&SignedDistanceFunction::getDist)
			.def("writePVPScript", writePVPScript)
			.def("__str__",&SignedDistanceFunction::toString)
	;
    class_<std::vector<SignedDistanceFunction*>>("std_vector_SDF_")
        .def(vector_indexing_suite<std::vector<SignedDistanceFunction*>>() )
	;
	class_<SDF_PlantBox, bases<SignedDistanceFunction>>("SDF_PlantBox",init<double,double,double>())
			.def("getDist",&SDF_PlantBox::getDist)
                        .def("__str__",&SDF_PlantBox::toString)
	;
	class_<SDF_PlantContainer, bases<SignedDistanceFunction>>("SDF_PlantContainer",init<>())
			.def(init<double,double,double,double>())
			.def("getDist",&SDF_PlantContainer::getDist)
                        .def("__str__",&SDF_PlantContainer::toString)
	;
	class_<SDF_RotateTranslate, bases<SignedDistanceFunction>>("SDF_RotateTranslate",init<SignedDistanceFunction*,double,int,Vector3d&>())
			.def(init<SignedDistanceFunction*,Vector3d&>())
			.def("getDist",&SDF_RotateTranslate::getDist)
                        .def("__str__",&SDF_RotateTranslate::toString)
	;
	enum_<SDF_RotateTranslate::SDF_Axes>("SDF_Axis")
	    .value("xaxis", SDF_RotateTranslate::SDF_Axes::xaxis)
	    .value("yaxis", SDF_RotateTranslate::SDF_Axes::yaxis)
	    .value("zaxis", SDF_RotateTranslate::SDF_Axes::zaxis)
	;
	class_<SDF_Intersection, bases<SignedDistanceFunction>>("SDF_Intersection",init<std::vector<SignedDistanceFunction*>>())
			.def(init<SignedDistanceFunction*,SignedDistanceFunction*>())
			.def("getDist",&SDF_Intersection::getDist)
			.def("__str__",&SDF_Intersection::toString)
	;
	class_<SDF_Union, bases<SDF_Intersection>>("SDF_Union",init<std::vector<SignedDistanceFunction*>>())
			.def(init<SignedDistanceFunction*,SignedDistanceFunction*>())
			.def("getDist",&SDF_Union::getDist)
			.def("__str__",&SDF_Union::toString)
	;
	class_<SDF_Difference, bases<SDF_Intersection>>("SDF_Difference",init<std::vector<SignedDistanceFunction*>>())
			.def(init<SignedDistanceFunction*,SignedDistanceFunction*>())
			.def("getDist",&SDF_Difference::getDist)
			.def("__str__",&SDF_Difference::toString)
	;
	class_<SDF_Complement, bases<SignedDistanceFunction>>("SDF_Complement",init<SignedDistanceFunction*>())
			.def("getDist",&SDF_Complement::getDist)
			.def("__str__",&SDF_Complement::toString)
	;
	class_<SDF_HalfPlane, bases<SignedDistanceFunction>>("SDF_HalfPlane",init<Vector3d&,Vector3d&>())
			.def(init<Vector3d&,Vector3d&,Vector3d&>())
			.def("getDist",&SDF_HalfPlane::getDist)
			.def_readwrite("o", &SDF_HalfPlane::o)
			.def_readwrite("n", &SDF_HalfPlane::n)
			.def_readwrite("p1", &SDF_HalfPlane::p1)
			.def_readwrite("p2", &SDF_HalfPlane::p2)
			.def("__str__",&SDF_HalfPlane::toString)
	;
	/*
	 * soil.h
	 */
	class_<SoilProperty_Wrap,boost::noncopyable>("SoilProperty",init<>())
			.def("getValue",pure_virtual(&SoilProperty::getValue))
			.def("__str__",pure_virtual(&SoilProperty::toString))
	;
	class_<SoilPropertySDF, bases<SoilProperty>>("SoilPropertySDF",init<>())
			.def(init<SignedDistanceFunction*, double, double, double>())
			.def_readwrite("sdf", &SoilPropertySDF::sdf)
			.def_readwrite("fmax", &SoilPropertySDF::fmax)
			.def_readwrite("fmin", &SoilPropertySDF::fmin)
			.def_readwrite("slope", &SoilPropertySDF::slope)
			.def("__str__",&SoilPropertySDF::toString)
	;
	/*
	 * ModelParameter.h
	 */
	class_<RootTypeParameter>("RootTypeParameter", init<>())
			.def(init<RootTypeParameter&>())
			.def("realize",&RootTypeParameter::realize)
			.def("getLateralType",&RootTypeParameter::getLateralType)
			.def("getK",&RootTypeParameter::getK)
			.def_readwrite("type", &RootTypeParameter::type)
			.def_readwrite("lb", &RootTypeParameter::lb)
			.def_readwrite("lbs", &RootTypeParameter::lbs)
			.def_readwrite("la", &RootTypeParameter::la)
			.def_readwrite("las", &RootTypeParameter::las)
			.def_readwrite("ln", &RootTypeParameter::ln)
			.def_readwrite("lns", &RootTypeParameter::lns)
			.def_readwrite("nob", &RootTypeParameter::nob)
			.def_readwrite("nobs", &RootTypeParameter::nobs)
			.def_readwrite("r", &RootTypeParameter::r)
			.def_readwrite("rs", &RootTypeParameter::rs)
			.def_readwrite("a", &RootTypeParameter::a)
			.def_readwrite("a_s", &RootTypeParameter::as) // as is a keyword in python
			.def_readwrite("colorR", &RootTypeParameter::colorR)
			.def_readwrite("colorG", &RootTypeParameter::colorG)
			.def_readwrite("colorB", &RootTypeParameter::colorB)
			.def_readwrite("tropismT", &RootTypeParameter::tropismT)
			.def_readwrite("tropismN", &RootTypeParameter::tropismN)
			.def_readwrite("tropismS", &RootTypeParameter::tropismS)
			.def_readwrite("dx", &RootTypeParameter::dx)
			.def_readwrite("theta", &RootTypeParameter::theta)
			.def_readwrite("thetas", &RootTypeParameter::thetas)
			.def_readwrite("rlt", &RootTypeParameter::rlt)
			.def_readwrite("rlts", &RootTypeParameter::rlts)
			.def_readwrite("gf", &RootTypeParameter::gf)
			.def_readwrite("name", &RootTypeParameter::name)
			.def_readwrite("successor", &RootTypeParameter::successor)
			.def_readwrite("successorP", &RootTypeParameter::successorP)
			.def_readwrite("se", &RootTypeParameter::se)
			.def_readwrite("sa", &RootTypeParameter::sa)
			.def_readwrite("sbp", &RootTypeParameter::sbp)
			.def("__str__",&RootTypeParameter::toString)
	;
	class_<RootParameter>("RootParameter", init<>())
			.def(init<int , double, double, const std::vector<double>&, int, double, double, double, double>())
			.def("set",&RootParameter::set)
			.def_readwrite("type", &RootParameter::type)
			.def_readwrite("lb", &RootParameter::lb)
			.def_readwrite("la", &RootParameter::la)
			.def_readwrite("ln", &RootParameter::ln)
			.def_readwrite("nob", &RootParameter::nob)
			.def_readwrite("r", &RootParameter::r)
			.def_readwrite("a", &RootParameter::a)
			.def_readwrite("theta", &RootParameter::theta)
			.def_readwrite("rlt", &RootParameter::rlt)
			.def("getK",&RootParameter::toString)
			.def("__str__",&RootParameter::toString)
	;
	class_<RootSystemParameter>("RootSystemParameter", init<>())
			.def("set",&RootSystemParameter::set)
			.def_readwrite("seedPos", &RootSystemParameter::seedPos)
			.def_readwrite("firstB", &RootSystemParameter::firstB)
			.def_readwrite("delayB", &RootSystemParameter::delayB)
			.def_readwrite("maxB", &RootSystemParameter::maxB)
			.def_readwrite("nC", &RootSystemParameter::nC)
			.def_readwrite("firstSB", &RootSystemParameter::firstSB)
			.def_readwrite("delaySB", &RootSystemParameter::delaySB)
			.def_readwrite("delayRC", &RootSystemParameter::delayRC)
			.def_readwrite("nz", &RootSystemParameter::nz)
			.def("__str__",&RootSystemParameter::toString)
	;
	/**
	 * Root.h (no members, just data)
	 */
    class_<Root, Root*>("Root", init<RootSystem*, int, Vector3d, double, Root*, double, int>())
		.def("__str__",&Root::toString)
	    .def_readwrite("rootsystem", &Root::rootsystem)
	    .def_readwrite("param", &Root::param)
	    .def_readwrite("id", &Root::id)
	    .def_readwrite("parent_base_length", &Root::parent_base_length)
	    .def_readwrite("parent_ni", &Root::parent_ni)
	    .def_readwrite("alive", &Root::alive)
	    .def_readwrite("active", &Root::active)
	    .def_readwrite("age", &Root::age)
	    .def_readwrite("length", &Root::length)
	    .def_readwrite("parent", &Root::parent)
	    .def_readwrite("laterals", &Root::laterals)
    ;
    class_<std::vector<Root*>>("std_vector_Root_")
        .def(vector_indexing_suite<std::vector<Root*>>() )
	;
    class_<std::vector<std::vector<Vector3d>>>("std_vector_vector_Vector3d_")
		.def(vector_indexing_suite<std::vector<std::vector<Vector3d>>>() )
	;
    class_<std::vector<std::vector<double>>>("std_vector_vector_double_")
		.def(vector_indexing_suite<std::vector<std::vector<double>>>() )
	;
	/*
	 * RootSystem.h
	 */
    class_<RootSystem>("RootSystem")
		.def("setRootTypeParameter", &RootSystem::setRootTypeParameter)
		.def("getRootTypeParameter", &RootSystem::getRootTypeParameter, return_value_policy<reference_existing_object>())
		.def("setRootSystemParameter", &RootSystem::setRootSystemParameter)
		.def("getRootSystemParameter", &RootSystem::getRootSystemParameter, return_value_policy<reference_existing_object>()) // tutorial: "naive (dangerous) approach"
		.def("openFile", &RootSystem::openFile, openFile_overloads())
		.def("setGeometry", &RootSystem::setGeometry)
		.def("setSoil", &RootSystem::setSoil)
		.def("reset", &RootSystem::reset)
		.def("initialize", &RootSystem::initialize, initialize_overloads())
		.def("simulate",simulate1, simulate1_overloads())
		.def("simulate",simulate2)
		.def("getSimTime", &RootSystem::getSimTime)
		.def("getNumberOfNodes", &RootSystem::getNumberOfNodes)
		.def("getNumberOfSegments", &RootSystem::getNumberOfSegments)
		.def("getRoots", &RootSystem::getRoots)
		.def("getBaseRoots", &RootSystem::getBaseRoots)
		.def("getNodes", &RootSystem::getNodes)
		.def("getPolylines", &RootSystem::getPolylines)
		.def("getSegments", &RootSystem::getSegments)
		.def("getSegmentsOrigin", &RootSystem::getSegmentsOrigin)
		.def("getNETimes", &RootSystem::getNETimes)
		.def("getScalar", &RootSystem::getScalar)
		.def("getRootTips", &RootSystem::getRootTips)
		.def("getRootBases", &RootSystem::getRootBases)
		.def("write", &RootSystem::write)
		.def("setSeed",&RootSystem::setSeed)
	;
    enum_<RootSystem::TropismTypes>("TropismType")
    	.value("plagio", RootSystem::TropismTypes::tt_plagio)
		.value("gravi", RootSystem::TropismTypes::tt_gravi)
    	.value("exo", RootSystem::TropismTypes::tt_exo)
    	.value("hydro", RootSystem::TropismTypes::tt_hydro)
    ;
    enum_<RootSystem::GrowthFunctionTypes>("GrowthFunctionType")
    	.value("negexp", RootSystem::GrowthFunctionTypes::gft_negexp)
		.value("linear", RootSystem::GrowthFunctionTypes::gft_linear)
    ;
	enum_<RootSystem::ScalarTypes>("ScalarType")
	    .value("type", RootSystem::ScalarTypes::st_type)
	    .value("radius", RootSystem::ScalarTypes::st_radius)
	    .value("order", RootSystem::ScalarTypes::st_order)
	    .value("time", RootSystem::ScalarTypes::st_time)
	    .value("length", RootSystem::ScalarTypes::st_length)
	    .value("surface", RootSystem::ScalarTypes::st_surface)
		.value("one", RootSystem::ScalarTypes::st_one)
		.value("userdata1", RootSystem::ScalarTypes::st_userdata1)
		.value("userdata2", RootSystem::ScalarTypes::st_userdata2)
		.value("userdata3", RootSystem::ScalarTypes::st_userdata3)
		.value("parenttype", RootSystem::ScalarTypes::st_parenttype)
	;
    /*
     * analysis.h
     */
    class_<SegmentAnalyser>("SegmentAnalyser")
    .def(init<RootSystem&>())
    .def(init<SegmentAnalyser&>())
	.def("addSegments",addSegments1)
	.def("addSegments",addSegments2)
	.def("crop", &SegmentAnalyser::crop)
	.def("filter", filter1)
	.def("filter", filter2)
	.def("pack", &SegmentAnalyser::pack)
	.def("getScalar", &SegmentAnalyser::getScalar)
	.def("getSegmentLength", &SegmentAnalyser::getSegmentLength)
	.def("getSummed", getSummed1)
	.def("getSummed", getSummed2)
	.def("distribution", distribution_1)
	.def("distribution", distribution_2)
	.def("distribution2", distribution2_1)
	.def("distribution2", distribution2_2)
    .def("getRoots", &SegmentAnalyser::getRoots)
	.def("getNumberOfRoots", &SegmentAnalyser::getNumberOfRoots)
	.def("cut", cut1)
	.def("addUserData", &SegmentAnalyser::addUserData)
	.def("clearUserData", &SegmentAnalyser::clearUserData)
	.def("write", &SegmentAnalyser::write)
	// .def("cut", cut2) // not working, see top definition of cut2
    ;
    class_<std::vector<SegmentAnalyser>>("std_vector_SegmentAnalyser_")
        .def(vector_indexing_suite<std::vector<SegmentAnalyser>>() )
	;
    /*
     * example_exudation.h (rather specific for Cheng)
     */
    class_<ExudationParameters>("ExudationParameters")
		.def_readwrite("M", &ExudationParameters::M)
		.def_readwrite("Dt", &ExudationParameters::Dt)
		.def_readwrite("Dl", &ExudationParameters::Dl)
		.def_readwrite("theta", &ExudationParameters::theta)
		.def_readwrite("R", &ExudationParameters::R)
		.def_readwrite("lambda_", &ExudationParameters::lambda_)
		.def_readwrite("age_r", &ExudationParameters::age_r)
		.def_readwrite("tip", &ExudationParameters::tip)
		.def_readwrite("v", &ExudationParameters::v)
		.def_readwrite("pos", &ExudationParameters::pos)
		;
     def("getExudateConcentration", getExudateConcentration);
}

/*
 *  currently not exposed.. (because not needed)
 *
	class_<Vector2d>("Vector2d", init<>())
			.def(init<double,double>())
			.def(init<Vector2d&>())
			.def_readwrite("x",&Vector2d::x)
			.def_readwrite("y",&Vector2d::y)
			.def("__str__",&Vector2d::toString)
	;
 */

#endif /* PYTHON_WRAPPER */

#endif /* PY_ROOTBOX_H_ */

