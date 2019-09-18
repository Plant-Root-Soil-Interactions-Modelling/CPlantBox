// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef PY_ROOTBOX_H_
#define PY_ROOTBOX_H_

/**
 *  A Python module for CRootbox based on boost.python
 */

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/call.hpp>

#include "mymath.h"
#include "sdf.h"
#include "leafparameter.h"
#include "stemparameter.h"
#include "RootSystem.h"
#include "sdf_rs.h"
#include "analysis.h"

//#include "../examples/example_exudation.h"

namespace CRootBox {

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

OrganRandomParameter* (Organism::*getOrganRandomParameter1)(int otype, int subType) const = &Organism::getOrganRandomParameter;
std::vector<OrganRandomParameter*> (Organism::*getOrganRandomParameter2)(int otype) const = &Organism::getOrganRandomParameter;
void (OrganRandomParameter::*readXML1)(std::string name) = &OrganRandomParameter::readXML;
void (OrganRandomParameter::*writeXML1)(std::string name) const = &OrganRandomParameter::writeXML;
void (OrganRandomParameter::*bindDoubleParameter)(std::string name, double* d, std::string descr, double* dev) = &OrganRandomParameter::bindParameter;
void (OrganRandomParameter::*bindIntParameter)(std::string name, int* i, std::string descr, double* dev) = &OrganRandomParameter::bindParameter;

void (Organ::*addNode1)(Vector3d n, double t) = &Organ::addNode;
void (Organ::*addNode2)(Vector3d n, int id, double t)= &Organ::addNode;
std::vector<Organ*> (Organ::*getOrgans1)(int otype) = &Organ::getOrgans;
void (Organ::*getOrgans2)(int otype, std::vector<Organ*>& v) = &Organ::getOrgans;

void (RootSystem::*simulate1)(double dt, bool silence) = &RootSystem::simulate;
void (RootSystem::*simulate2)() = &RootSystem::simulate;
void (RootSystem::*simulate3)(double dt, double maxinc, ProportionalElongation* f_se, bool silence) = &RootSystem::simulate;
void (RootSystem::*initialize1)() = &RootSystem::initialize;
void (RootSystem::*initialize2)(int basal, int shootborne) = &RootSystem::initialize;

RootRandomParameter* (RootSystem::*getRootTypeParameter1)(int subType) const = &RootSystem::getRootTypeParameter;
std::vector<RootRandomParameter*> (RootSystem::*getRootTypeParameter2)() const = &RootSystem::getRootTypeParameter;

void (SegmentAnalyser::*addSegments1)(const Organism& plant) = &SegmentAnalyser::addSegments;
void (SegmentAnalyser::*addSegments2)(const SegmentAnalyser& a) = &SegmentAnalyser::addSegments;
void (SegmentAnalyser::*filter1)(std::string name, double min, double max) = &SegmentAnalyser::filter;
void (SegmentAnalyser::*filter2)(std::string name, double value) = &SegmentAnalyser::filter;
double (SegmentAnalyser::*getSummed1)(std::string name) const = &SegmentAnalyser::getSummed;
double (SegmentAnalyser::*getSummed2)(std::string name, SignedDistanceFunction* geometry) const = &SegmentAnalyser::getSummed;
std::vector<double> (SegmentAnalyser::*distribution_1)(std::string name, double top, double bot, int n, bool exact) const = &SegmentAnalyser::distribution;
std::vector<SegmentAnalyser> (SegmentAnalyser::*distribution_2)(double top, double bot, int n) const = &SegmentAnalyser::distribution;
std::vector<std::vector<double>> (SegmentAnalyser::*distribution2_1)(std::string name, double top, double bot, double left, double right, int n, int m, bool exact) const = &SegmentAnalyser::distribution2;
std::vector<std::vector<SegmentAnalyser>> (SegmentAnalyser::*distribution2_2)(double top, double bot, double left, double right, int n, int m) const = &SegmentAnalyser::distribution2;
SegmentAnalyser (SegmentAnalyser::*cut1)(const SDF_HalfPlane& plane) const = &SegmentAnalyser::cut;

/**
 * Default arguments: no idea how to do it by hand, magic everywhere...
 */
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getOrgans_overloads, getOrgans, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getParameter_overloads, getParameter, 1, 3);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getSummed_overloads, getSummed, 1, 2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getNumberOfSegments_overloads, getNumberOfSegments, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getPolylines_overloads, getPolylines, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getPolylineCTs_overloads, getPolylineCTs, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getSegments_overloads, getSegments, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getSegmentCTs_overloads, getSegmentCTs, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getSegmentOrigins_overloads, getSegmentOrigins, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getNewSegments_overloads, getNewSegments, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getNewSegmentOrigins_overloads, getNewSegmentOrigins, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(openFile_overloads,openFile,1,2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(simulate1_overloads,simulate,1,2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(simulate3_overloads,simulate,3,4);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getValue_overloads,getValue,1,2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(tropismObjective_overloads,tropismObjective,5,6);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getNumberOfRoots_overloads,getNumberOfRoots,0,1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(toString_overloads, toString, 0, 1);
// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(bindParameter_overloads, bindParameter, 2, 4);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(readParameters_overloads, readParameters, 1, 2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(writeParameters_overloads, writeParameters, 1, 3);


/**
 * Virtual functions
 */
class SoilLookUp_Wrap : public SoilLookUp, public wrapper<SoilLookUp> {
public:

    virtual double getValue(const Vector3d& pos, const Organ* o = nullptr) const override {
        return this->get_override("getValue")(pos, o);
    }

    virtual std::string toString() const override {
        return this->get_override("toString")();
    }

};

//class Tropism_Wrap : public Tropism, public wrapper<Tropism> {
//public:
//
//    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Organ* o = nullptr) override {
//        return this->get_override("tropismObjective")(pos, old, a, b, dx, o);
//    }
//
//    virtual Tropism* copy() override {
//        return this->get_override("copy")();
//    }
//
//}; // dont know how that works...



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
    class_<SignedDistanceFunction, SignedDistanceFunction*>("SignedDistanceFunction")
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
    class_<SoilLookUp_Wrap, SoilLookUp_Wrap*, boost::noncopyable>("SoilLookUp",init<>())
        .def("getValue",&SoilLookUp_Wrap::getValue)
        .def("__str__",&SoilLookUp_Wrap::toString)
        ;
    class_<SoilLookUpSDF, SoilLookUpSDF*, bases<SoilLookUp>>("SoilLookUpSDF",init<>())
        .def(init<SignedDistanceFunction*, double, double, double>())
        .def_readwrite("sdf", &SoilLookUpSDF::sdf)
        .def_readwrite("fmax", &SoilLookUpSDF::fmax)
        .def_readwrite("fmin", &SoilLookUpSDF::fmin)
        .def_readwrite("slope", &SoilLookUpSDF::slope)
        .def("__str__",&SoilLookUpSDF::toString)
        ;
    class_<MultiplySoilLookUps, MultiplySoilLookUps*, bases<SoilLookUp>>("MultiplySoilLookUps",init<SoilLookUp*, SoilLookUp*>())
        .def(init<std::vector<SoilLookUp*>>())
        .def("getValue", &MultiplySoilLookUps::getValue, getValue_overloads())
        .def("__str__",&MultiplySoilLookUps::toString)
        ;
    class_<ProportionalElongation, ProportionalElongation*, bases<SoilLookUp>>("ProportionalElongation",init<>())
        .def("getValue", &ProportionalElongation::getValue, getValue_overloads())
        .def("setScale", &ProportionalElongation::setScale)
        .def("setBaseLookUp", &ProportionalElongation::setBaseLookUp)
        .def("__str__",&ProportionalElongation::toString)
        ;
    class_<Grid1D, Grid1D*, bases<SoilLookUp>>("Grid1D",init<>())
        .def(init<size_t, std::vector<double>, std::vector<double>>())
        .def("getValue", &Grid1D::getValue, getValue_overloads())
        .def("__str__",&Grid1D::toString)
        .def("map",&Grid1D::map)
        .def_readwrite("n", &Grid1D::n)
        .def_readwrite("grid", &Grid1D::grid)
        .def_readwrite("data", &Grid1D::data)
        ;
    class_<EquidistantGrid1D, EquidistantGrid1D*, bases<Grid1D>>("EquidistantGrid1D",init<double, double, size_t>())
        .def(init<double, double, std::vector<double>>())
        .def("getValue", &EquidistantGrid1D::getValue, getValue_overloads())
        .def("map",&EquidistantGrid1D::map)
        .def_readwrite("n", &EquidistantGrid1D::n)
        .def_readwrite("grid", &EquidistantGrid1D::grid)
        .def_readwrite("data", &EquidistantGrid1D::data)
        .def("__str__",&EquidistantGrid1D::toString)
        ;
    /*
     * organparameter.h
     */
    class_<OrganSpecificParameter, OrganSpecificParameter*>("OrganSpecificParameter")
        .def_readwrite("subType",&OrganSpecificParameter::subType)
        .def("__str__",&OrganSpecificParameter::toString)
        ;
    class_<OrganRandomParameter, OrganRandomParameter*>("OrganRandomParameter", init<Organism*>())
        .def("copy",&OrganRandomParameter::copy, return_value_policy<reference_existing_object>())
        .def("realize",&OrganRandomParameter::realize, return_value_policy<reference_existing_object>())
        .def("getParameter",&OrganRandomParameter::getParameter)
        .def("writeXML",writeXML1)
        .def("readXML",readXML1)
//        .def("bindIntParameter",bindIntParameter, bindParameter_overloads()) // not working, can't pass int*
//        .def("bindDoubleParameter",bindDoubleParameter, bindParameter_overloads()) // not working, can't pass double*
        .def_readwrite("name",&OrganRandomParameter::name)
        .def_readwrite("organType",&OrganRandomParameter::organType)
        .def_readwrite("subType",&OrganRandomParameter::subType)
        .def("__str__",&OrganRandomParameter::toString, toString_overloads())
        ;
    class_<std::vector<OrganRandomParameter*>>("std_vector_OrganRandomParameter_")
        .def(vector_indexing_suite<std::vector<OrganRandomParameter*>>() )
        ;
    /**
     * Organ.h
     */
    class_<Organ, Organ*>("Organ", init<Organism*, Organ*, int, int, double>())
        .def(init<int, OrganSpecificParameter*, bool, bool, double, double, bool, int>())
        .def("copy",&Organ::copy, return_value_policy<reference_existing_object>()) // manage_new_object
        .def("organType",&Organ::organType)
        .def("simulate",&Organ::simulate, simulate1_overloads())
        .def("setParent",&Organ::setParent)
        .def("getParent",&Organ::getParent, return_value_policy<reference_existing_object>())
        .def("setOrganism",&Organ::setOrganism)
        .def("addChild",&Organ::addChild)
        .def("getId",&Organ::getId)
        .def("getParam",&Organ::getParam, return_value_policy<reference_existing_object>())
        .def("getOrganRandomParameter",&Organ::getOrganRandomParameter, return_value_policy<reference_existing_object>())
        .def("isAlive",&Organ::isAlive)
        .def("isActive",&Organ::isActive)
        .def("getAge",&Organ::getAge)
        .def("getLength",&Organ::getLength)
        .def("getNumberOfNodes",&Organ::getNumberOfNodes)
        .def("getNode",&Organ::getNode)
        .def("getNodeId",&Organ::getNodeId)
        .def("getNodeCT",&Organ::getNodeCT)
        .def("addNode",addNode1)
        .def("addNode",addNode2)
        .def("getSegments",&Organ::getSegments)
        .def("hasMoved",&Organ::hasMoved)
        .def("getOldNumberOfNodes",&Organ::getOldNumberOfNodes)
        .def("getOrgans", getOrgans1, getOrgans_overloads())
        .def("getOrgans", getOrgans2)
        .def("getParameter",&Organ::getParameter)
        .def("__str__",&Organ::toString)
        ;
    class_<std::vector<Organ*>>("std_vector_Organ_")
        .def(vector_indexing_suite<std::vector<Organ*>>() )
        ;
    /*
     * Organism.h
     */
    class_<Organism, Organism*>("Organism", init<>())
        .def(init<Organism&>())
        .def("organTypeNumber", &Organism::organTypeNumber)
        .def("organTypeName", &Organism::organTypeName)
        .def("getOrganRandomParameter", getOrganRandomParameter1, return_value_policy<reference_existing_object>())
        .def("getOrganRandomParameter", getOrganRandomParameter2)
        .def("setOrganRandomParameter", &Organism::setOrganRandomParameter)

        .def("addOrgan", &Organism::addOrgan)
        .def("initialize", &Organism::initialize)
        .def("simulate", &Organism::simulate, simulate1_overloads())
        .def("getSimTime", &Organism::getSimTime)

        .def("getOrgans", &Organism::getOrgans, getOrgans_overloads())
        .def("getParameter", &Organism::getParameter, getParameter_overloads())
        .def("getSummed", &Organism::getSummed, getSummed_overloads())

        .def("getNumberOfOrgans", &Organism::getNumberOfOrgans)
        .def("getNumberOfNodes", &Organism::getNumberOfNodes)
        .def("getNumberOfSegments", &Organism::getNumberOfSegments, getNumberOfSegments_overloads())
        .def("getPolylines", &Organism::getPolylines, getPolylines_overloads())
        .def("getPolylineCTs", &Organism::getPolylineCTs, getPolylineCTs_overloads())
        .def("getNodes", &Organism::getNodes)
        .def("getNodeCTs", &Organism::getNodeCTs)
        .def("getSegments", &Organism::getSegments, getSegments_overloads())
        .def("getSegmentCTs", &Organism::getSegmentCTs, getSegmentCTs_overloads())
        .def("getSegmentOrigins", &Organism::getSegmentOrigins,  getSegmentOrigins_overloads())

        .def("getNumberOfNewNodes", &Organism::getNumberOfNewNodes)
        .def("getNumberOfNewOrgans", &Organism::getNumberOfNewOrgans)
        .def("getUpdatedNodeIndices", &Organism::getUpdatedNodeIndices)
        .def("getUpdatedNodes", &Organism::getUpdatedNodes)
        .def("getUpdatedNodeCTs", &Organism::getUpdatedNodeCTs)

        .def("getNewNodes", &Organism::getNewNodes)
        .def("getNewNodeCTs", &Organism::getNewNodeCTs)
        .def("getNewSegments", &Organism::getNewSegments,  getNewSegments_overloads())
        .def("getNewSegmentOrigins", &Organism::getNewSegmentOrigins,  getNewSegmentOrigins_overloads())

        .def("readParameters", &Organism::readParameters, readParameters_overloads())
        .def("writeParameters", &Organism::writeParameters, writeParameters_overloads())
        .def("writeRSML", &Organism::writeRSML)
        .def("getRSMLSkip", &Organism::getRSMLSkip)
        .def("setRSMLSkip", &Organism::setRSMLSkip)
        .def("getRSMLProperties", &Organism::getRSMLProperties, return_value_policy<copy_non_const_reference>())

        .def("getOrganIndex", &Organism::getOrganIndex)
        .def("getNodeIndex", &Organism::getNodeIndex)

        .def("setSeed", &Organism::setSeed)
        .def("rand", &Organism::rand)
        .def("randn", &Organism::randn)
        .def("__str__",&Organism::toString)
        ;
    enum_<Organism::OrganTypes>("OrganTypes")
        .value("organ", Organism::OrganTypes::ot_organ)
        .value("seed", Organism::OrganTypes::ot_seed)
        .value("root", Organism::OrganTypes::ot_root)
        .value("stem", Organism::OrganTypes::ot_stem)
        .value("leaf", Organism::OrganTypes::ot_leaf)
        ;
    /*
     * sdf_rs.h
     */
    class_<SDF_RootSystem,SDF_RootSystem*, bases<SignedDistanceFunction>>("SDF_RootSystem", init<std::vector<Vector3d>, std::vector<Vector2i>, std::vector<double>, double>())
        .def(init<Root&, double>())
        .def(init<RootSystem&, double>())
        .def("getDist",&SDF_RootSystem::getDist)
        .def("__str__",&SDF_RootSystem::toString)
    ;
    /**
     * tropism.h
     */
//    class_<Tropism_Wrap, Tropism_Wrap*, boost::noncopyable>("Tropism",init<>())
//                .def("getHeading",&Tropism_Wrap::getHeading)
//                .def("tropismObjective",&Tropism_Wrap::tropismObjective, tropismObjective_overloads())
//                .def("copy",&Tropism_Wrap::copy, return_value_policy<reference_existing_object>())
//                .def("setTropismParameter",&Tropism_Wrap::setTropismParameter)
//                .def("setGeometry",&Tropism_Wrap::setGeometry) // todo dont know how that works
                ;
    class_<Tropism, Tropism*>("TropismBase",init<Organism*>()) // Base class for the following tropisms
                .def(init<Organism*, double, double>())
                .def("getHeading",&Tropism::getHeading)
                .def("tropismObjective",&Tropism::tropismObjective, tropismObjective_overloads())
                .def("copy",&Tropism::copy, return_value_policy<reference_existing_object>())
                .def("setTropismParameter",&Tropism::setTropismParameter)
                .def("setGeometry",&Tropism::setGeometry)

                ;
    class_<Gravitropism, Gravitropism*, bases<Tropism>>("Gravitropism",init<Organism*, double, double>())
        ;
    class_<Plagiotropism, Plagiotropism*, bases<Tropism>>("Plagiotropism",init<Organism*,double, double>())
        ;
    class_<Exotropism, Exotropism*, bases<Tropism>>("Exotropism",init<Organism*,double, double>())
        ;
    class_<Hydrotropism, Hydrotropism*, bases<Tropism>>("Hydrotropism",init<Organism*,double, double, SoilLookUp*>())
        ;
    //	class_<CombinedTropism, CombinedTropism*, bases<Tropism>>("CombinedTropism",init<>()) // Todo needs some extra work
    //	;
    /*
     * analysis.h
     */
    class_<SegmentAnalyser, SegmentAnalyser*>("SegmentAnalyser")
        .def(init<RootSystem&>())
        .def(init<SegmentAnalyser&>())
        .def("addSegments",addSegments1)
        .def("addSegments",addSegments2)
        .def("crop", &SegmentAnalyser::crop)
        .def("filter", filter1)
        .def("filter", filter2)
        .def("pack", &SegmentAnalyser::pack)
        .def("getParameter", &SegmentAnalyser::getParameter)
        .def("getSegmentLength", &SegmentAnalyser::getSegmentLength)
        .def("getSummed", getSummed1)
        .def("getSummed", getSummed2)
        .def("distribution", distribution_1)
        .def("distribution", distribution_2)
        .def("distribution2", distribution2_1)
        .def("distribution2", distribution2_2)
        .def("getOrgans", &SegmentAnalyser::getOrgans)
        .def("getNumberOfOrgans", &SegmentAnalyser::getNumberOfOrgans)
        .def("cut", cut1)
        .def("addUserData", &SegmentAnalyser::addUserData)
        .def("clearUserData", &SegmentAnalyser::clearUserData)
        .def("write", &SegmentAnalyser::write)
        .def_readwrite("nodes", &SegmentAnalyser::nodes)
        .def_readwrite("segments", &SegmentAnalyser::segments)
        .def_readwrite("segCTs", &SegmentAnalyser::segCTs)
        // .def("cut", cut2) // not working, see top definition of cut2
        ;
    class_<std::vector<SegmentAnalyser>>("std_vector_SegmentAnalyser_")
            .def(vector_indexing_suite<std::vector<SegmentAnalyser>>() )
            ;
    /*
     * rootparameter.h
     */
    class_<RootRandomParameter, RootRandomParameter*, bases<OrganRandomParameter>>("RootRandomParameter", init<Organism*>())
                .def("getLateralType",&RootRandomParameter::getLateralType)
                .def("getK",&RootRandomParameter::getK)
                .def_readwrite("type", &RootRandomParameter::subType)
                .def_readwrite("lb", &RootRandomParameter::lb)
                .def_readwrite("lbs", &RootRandomParameter::lbs)
                .def_readwrite("la", &RootRandomParameter::la)
                .def_readwrite("las", &RootRandomParameter::las)
                .def_readwrite("ln", &RootRandomParameter::ln)
                .def_readwrite("lns", &RootRandomParameter::lns)
                .def_readwrite("nob", &RootRandomParameter::nob)
                .def_readwrite("nobs", &RootRandomParameter::nobs)
                .def_readwrite("r", &RootRandomParameter::r)
                .def_readwrite("rs", &RootRandomParameter::rs)
                .def_readwrite("a", &RootRandomParameter::a)
                .def_readwrite("a_s", &RootRandomParameter::as) // as is a keyword in python
                .def_readwrite("colorR", &RootRandomParameter::colorR)
                .def_readwrite("colorG", &RootRandomParameter::colorG)
                .def_readwrite("colorB", &RootRandomParameter::colorB)
                .def_readwrite("tropismT", &RootRandomParameter::tropismT)
                .def_readwrite("tropismN", &RootRandomParameter::tropismN)
                .def_readwrite("tropismS", &RootRandomParameter::tropismS)
                .def_readwrite("dx", &RootRandomParameter::dx)
                .def_readwrite("theta", &RootRandomParameter::theta)
                .def_readwrite("thetas", &RootRandomParameter::thetas)
                .def_readwrite("rlt", &RootRandomParameter::rlt)
                .def_readwrite("rlts", &RootRandomParameter::rlts)
                .def_readwrite("gf", &RootRandomParameter::gf)
                .def_readwrite("name", &RootRandomParameter::name)
                .def_readwrite("successor", &RootRandomParameter::successor)
                .def_readwrite("successorP", &RootRandomParameter::successorP)
                .def_readwrite("f_tf", &RootRandomParameter::f_tf)
                .def_readwrite("f_gf", &RootRandomParameter::f_gf)
                .def_readwrite("f_se", &RootRandomParameter::f_se)
                .def_readwrite("f_sa", &RootRandomParameter::f_sa)
                .def_readwrite("f_sbp", &RootRandomParameter::f_sbp)
                .def("__str__",&RootRandomParameter::toString, toString_overloads())
                ;
    class_<std::vector<RootRandomParameter*>>("std_vector_RootRandomParameter_")
        .def(vector_indexing_suite<std::vector<RootRandomParameter*>>() )
        ;
    class_<RootSpecificParameter, RootSpecificParameter*, bases<OrganSpecificParameter> >("RootSpecificParameter", init<>())
                .def(init<int , double, double, const std::vector<double>&, int, double, double, double, double>())
                .def_readwrite("subType", &RootSpecificParameter::subType)
                .def_readwrite("lb", &RootSpecificParameter::lb)
                .def_readwrite("la", &RootSpecificParameter::la)
                .def_readwrite("ln", &RootSpecificParameter::ln)
                .def_readwrite("nob", &RootSpecificParameter::nob)
                .def_readwrite("r", &RootSpecificParameter::r)
                .def_readwrite("a", &RootSpecificParameter::a)
                .def_readwrite("theta", &RootSpecificParameter::theta)
                .def_readwrite("rlt", &RootSpecificParameter::rlt)
                .def("getK",&RootSpecificParameter::getK)
                .def("__str__",&RootSpecificParameter::toString)
                ;
    /*
     * seedparameter.h
     */
    class_<SeedRandomParameter, SeedRandomParameter*, bases<OrganRandomParameter>>("SeedRandomParameter", init<Organism*>())
            .def_readwrite("type", &SeedRandomParameter::subType)
            .def_readwrite("seedPos", &SeedRandomParameter::seedPos)
            .def_readwrite("seedPoss", &SeedRandomParameter::seedPoss)
            .def_readwrite("firstB", &SeedRandomParameter::firstB)
            .def_readwrite("firstBs", &SeedRandomParameter::firstBs)
            .def_readwrite("delayB", &SeedRandomParameter::delayB)
            .def_readwrite("delayBs", &SeedRandomParameter::delayBs)
            .def_readwrite("maxB", &SeedRandomParameter::maxB)
            .def_readwrite("maxBs", &SeedRandomParameter::maxBs)
            .def_readwrite("nC", &SeedRandomParameter::nC)
            .def_readwrite("nCs", &SeedRandomParameter::nCs)
            .def_readwrite("firstSB", &SeedRandomParameter::firstSB)
            .def_readwrite("firstSBs", &SeedRandomParameter::firstSBs)
            .def_readwrite("delaySB", &SeedRandomParameter::delaySB)
            .def_readwrite("delaySBs", &SeedRandomParameter::delaySBs)
            .def_readwrite("delayRC", &SeedRandomParameter::delayRC)
            .def_readwrite("delayRCs", &SeedRandomParameter::delayRCs)
            .def_readwrite("nz", &SeedRandomParameter::nz)
            .def_readwrite("nzs", &SeedRandomParameter::nzs)
            .def_readwrite("simtime", &SeedRandomParameter::simtime)
            .def_readwrite("simtime", &SeedRandomParameter::simtimes)
            .def("__str__",&SeedRandomParameter::toString, toString_overloads())
            ;
    class_<SeedSpecificParameter, SeedSpecificParameter*, bases<OrganSpecificParameter>>("SeedSpecificParameter", init<>())
          .def(init<int, Vector3d , double, int, int, int, double, double, double, double, double>())
          .def_readwrite("seedPos", &SeedSpecificParameter::seedPos)
          .def_readwrite("firstB", &SeedSpecificParameter::firstB)
          .def_readwrite("delayB", &SeedSpecificParameter::delayB)
          .def_readwrite("maxB", &SeedSpecificParameter::maxB)
          .def_readwrite("nC", &SeedSpecificParameter::nC)
          .def_readwrite("firstSB", &SeedSpecificParameter::firstSB)
          .def_readwrite("delaySB", &SeedSpecificParameter::delaySB)
          .def_readwrite("delayRC", &SeedSpecificParameter::delayRC)
          .def_readwrite("nz", &SeedSpecificParameter::nz)
          .def("__str__",&SeedSpecificParameter::toString)
          ;
    /*
     * leafparameter.h
     */
    class_<LeafRandomParameter, LeafRandomParameter*, bases<OrganRandomParameter>>("LeafRandomParameter", init<Organism*>())
                        .def("getLateralType",&LeafRandomParameter::getLateralType)
                        .def("getK",&LeafRandomParameter::getK)
                        .def_readwrite("type", &LeafRandomParameter::subType)
                        .def_readwrite("lb", &LeafRandomParameter::lb)
                        .def_readwrite("lbs", &LeafRandomParameter::lbs)
                        .def_readwrite("la", &LeafRandomParameter::la)
                        .def_readwrite("las", &LeafRandomParameter::las)
                        .def_readwrite("ln", &LeafRandomParameter::ln)
                        .def_readwrite("lns", &LeafRandomParameter::lns)
                        .def_readwrite("lnf", &LeafRandomParameter::lnf)
                        .def_readwrite("nob", &LeafRandomParameter::nob)
                        .def_readwrite("nobs", &LeafRandomParameter::nobs)
                        .def_readwrite("k", &LeafRandomParameter::k)
                        .def_readwrite("ks", &LeafRandomParameter::ks)
                        .def_readwrite("r", &LeafRandomParameter::r)
                        .def_readwrite("rs", &LeafRandomParameter::rs)
                        .def_readwrite("a", &LeafRandomParameter::a)
                        .def_readwrite("a_s", &LeafRandomParameter::as) // as is a keyword in python
						.def_readwrite("RotBeta", &LeafRandomParameter::RotBeta)
						.def_readwrite("BetaDev", &LeafRandomParameter::BetaDev)
						.def_readwrite("InitBeta", &LeafRandomParameter::InitBeta)
                        .def_readwrite("tropismT", &LeafRandomParameter::tropismT)
                        .def_readwrite("tropismN", &LeafRandomParameter::tropismN)
                        .def_readwrite("tropismS", &LeafRandomParameter::tropismS)
                        .def_readwrite("dx", &LeafRandomParameter::dx)
                        .def_readwrite("theta", &LeafRandomParameter::theta)
                        .def_readwrite("thetas", &LeafRandomParameter::thetas)
                        .def_readwrite("rlt", &LeafRandomParameter::rlt)
                        .def_readwrite("rlts", &LeafRandomParameter::rlts)
                        .def_readwrite("gf", &LeafRandomParameter::gf)
                        .def_readwrite("name", &LeafRandomParameter::name)
                        .def_readwrite("successor", &LeafRandomParameter::successor)
                        .def_readwrite("successorP", &LeafRandomParameter::successorP)
                        .def_readwrite("f_tf", &LeafRandomParameter::f_tf)
                        .def_readwrite("f_gf", &LeafRandomParameter::f_gf)
                        .def_readwrite("f_se", &LeafRandomParameter::f_se)
                        .def_readwrite("f_sa", &LeafRandomParameter::f_sa)
                        .def_readwrite("f_sbp", &LeafRandomParameter::f_sbp)
            .def("__str__",&LeafRandomParameter::toString, toString_overloads())
            ;
    class_<LeafSpecificParameter, LeafSpecificParameter*, bases<OrganSpecificParameter>>("LeafSpecificParameter", init<>())
          .def(init<int , double, double, const std::vector<double>&, double, double, double, double, double>())
          .def_readwrite("subType", &LeafSpecificParameter::subType)
          .def_readwrite("lb", &LeafSpecificParameter::lb)
          .def_readwrite("la", &LeafSpecificParameter::la)
          .def_readwrite("ln", &LeafSpecificParameter::ln)
          .def_readwrite("r", &LeafSpecificParameter::r)
          .def_readwrite("a", &LeafSpecificParameter::a)
          .def_readwrite("theta", &LeafSpecificParameter::theta)
          .def_readwrite("rlt", &LeafSpecificParameter::rlt)
          .def_readwrite("lnf", &LeafSpecificParameter::lnf)
          .def("getK",&LeafSpecificParameter::getK)
          .def("__str__",&SeedSpecificParameter::toString)
          ;
    /*
     * stemparameter.h
     */
    class_<StemRandomParameter, StemRandomParameter*, bases<OrganRandomParameter>>("StemRandomParameter", init<Organism*>())
                        .def("getLateralType",&StemRandomParameter::getLateralType)
                        .def("getK",&StemRandomParameter::getK)
                        .def_readwrite("type", &StemRandomParameter::subType)
                        .def_readwrite("lb", &StemRandomParameter::lb)
                        .def_readwrite("lbs", &StemRandomParameter::lbs)
                        .def_readwrite("la", &StemRandomParameter::la)
                        .def_readwrite("las", &StemRandomParameter::las)
                        .def_readwrite("ln", &StemRandomParameter::ln)
                        .def_readwrite("lns", &StemRandomParameter::lns)
                        .def_readwrite("lnf", &StemRandomParameter::lnf)
                        .def_readwrite("nob", &StemRandomParameter::nob)
                        .def_readwrite("nobs", &StemRandomParameter::nobs)
                        .def_readwrite("r", &StemRandomParameter::r)
                        .def_readwrite("rs", &StemRandomParameter::rs)
                        .def_readwrite("a", &StemRandomParameter::a)
                        .def_readwrite("a_s", &StemRandomParameter::as) // as is a keyword in python
						.def_readwrite("RotBeta", &StemRandomParameter::RotBeta)
						.def_readwrite("BetaDev", &StemRandomParameter::BetaDev)
						.def_readwrite("InitBeta", &StemRandomParameter::InitBeta)
                        .def_readwrite("tropismT", &StemRandomParameter::tropismT)
                        .def_readwrite("tropismN", &StemRandomParameter::tropismN)
                        .def_readwrite("tropismS", &StemRandomParameter::tropismS)
                        .def_readwrite("dx", &StemRandomParameter::dx)
                        .def_readwrite("theta", &StemRandomParameter::theta)
                        .def_readwrite("thetas", &StemRandomParameter::thetas)
                        .def_readwrite("rlt", &StemRandomParameter::rlt)
                        .def_readwrite("rlts", &StemRandomParameter::rlts)
                        .def_readwrite("gf", &StemRandomParameter::gf)
                        .def_readwrite("name", &StemRandomParameter::name)
                        .def_readwrite("successor", &StemRandomParameter::successor)
                        .def_readwrite("successorP", &StemRandomParameter::successorP)
                        .def_readwrite("f_tf", &StemRandomParameter::f_tf)
                        .def_readwrite("f_gf", &StemRandomParameter::f_gf)
                        .def_readwrite("f_se", &StemRandomParameter::f_se)
                        .def_readwrite("f_sa", &StemRandomParameter::f_sa)
                        .def_readwrite("f_sbp", &StemRandomParameter::f_sbp)
            .def("__str__",&StemRandomParameter::toString, toString_overloads())
            ;
    class_<StemSpecificParameter, StemSpecificParameter*, bases<OrganSpecificParameter>>("StemSpecificParameter", init<>())
          .def(init<int , double, double, const std::vector<double>&, double, double, double, double, double>())
          .def_readwrite("subType", &StemSpecificParameter::subType)
          .def_readwrite("lb", &StemSpecificParameter::lb)
          .def_readwrite("la", &StemSpecificParameter::la)
          .def_readwrite("ln", &StemSpecificParameter::ln)
          .def_readwrite("r", &StemSpecificParameter::r)
          .def_readwrite("a", &StemSpecificParameter::a)
          .def_readwrite("theta", &StemSpecificParameter::theta)
          .def_readwrite("rlt", &StemSpecificParameter::rlt)
          .def("getK",&StemSpecificParameter::getK)
          .def("__str__",&SeedSpecificParameter::toString)
          ;

    /**
     * Root.h
     */
    class_<Root, Root*, bases<Organ>>("Root", init<Organism*, int, Vector3d, double, Root*, double, int>())
            .def(init<int, OrganSpecificParameter*, bool, bool, double, double, Vector3d, double, int, bool, int >()) // with_custodian_and_ward<M,N> ?? dont know how to use it
            .def("calcCreationTime", &Root::calcCreationTime)
            .def("calcLength", &Root::calcLength)
            .def("calcAge", &Root::calcAge)
            .def("getRootTypeParameter", &Root::getRootTypeParameter, return_value_policy<reference_existing_object>())
            .def("param", &Root::param, return_value_policy<reference_existing_object>())
            .def("dx", &Root::dx)
            .def_readwrite("parent_base_length", &Root::parentBaseLength)
            .def_readwrite("parent_ni", &Root::parentNI)
            .def("__str__",&Root::toString)
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
    class_<RootSystem, RootSystem*, bases<Organism>>("RootSystem", init<>()) // bases<PlantBase>
             .def(init<RootSystem&>())
             .def("getRootTypeParameter", getRootTypeParameter1, return_value_policy<reference_existing_object>())
             .def("getRootTypeParameter", getRootTypeParameter2)
             .def("setRootSystemParameter", &RootSystem::setRootSystemParameter)
             .def("getRootSystemParameter", &RootSystem::getRootSystemParameter, return_value_policy<reference_existing_object>()) // tutorial: "naive (dangerous) approach"
             .def("openFile", &RootSystem::openFile, openFile_overloads())
             .def("setGeometry", &RootSystem::setGeometry)
             .def("setSoil", &RootSystem::setSoil)
             .def("reset", &RootSystem::reset)
             .def("initialize", initialize1)
             .def("initialize", initialize2)
             .def("setTropism", &RootSystem::setTropism)
             .def("simulate",simulate1, simulate1_overloads())
             .def("simulate",simulate2)
             .def("simulate",simulate3, simulate3_overloads())
             .def("getSimTime", &RootSystem::getSimTime)
             .def("getNumberOfNodes", &RootSystem::getNumberOfNodes)
             .def("getRoots", &RootSystem::getRoots)
             .def("getNumberOfSegments", &RootSystem::getNumberOfSegments, getNumberOfSegments_overloads())
             .def("getNumberOfRoots", &RootSystem::getNumberOfRoots, getNumberOfRoots_overloads())
             .def("getBaseRoots", &RootSystem::getBaseRoots)
             .def("getShootSegments", &RootSystem::getShootSegments)
             .def("getSegmentOrigins", &RootSystem::getSegmentOrigins)
             .def("getRootTips", &RootSystem::getRootTips)
             .def("getRootBases", &RootSystem::getRootBases)
             .def("getNumberOfNewNodes",&RootSystem::getNumberOfNewNodes)
             .def("push",&RootSystem::push)
             .def("pop",&RootSystem::pop)
             .def("write", &RootSystem::write)
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
    /*
     * exudation.h
     */
//    class_<ExudationModel, ExudationModel*>("ExudationModel", init<double, double, int, RootSystem&>())
//            .def(init<double, double, double, int, int, int, RootSystem&>())
//		    .def("calculate", &ExudationModel::calculate)
//		    .def_readwrite("Q", &ExudationModel::Q)
//		    .def_readwrite("Dl", &ExudationModel::Dl)
//		    .def_readwrite("theta", &ExudationModel::theta)
//		    .def_readwrite("R", &ExudationModel::R)
//		    .def_readwrite("k", &ExudationModel::k)
//		    .def_readwrite("l", &ExudationModel::l)
//		    .def_readwrite("type", &ExudationModel::type)
//		    .def_readwrite("n0", &ExudationModel::n0)
//		    .def_readwrite("thresh13", &ExudationModel::thresh13)
//            .def_readwrite("calc13", &ExudationModel::calc13)
//            .def_readwrite("observationRadius", &ExudationModel::observationRadius)
//            ;
//    enum_<ExudationModel::IntegrationType>("IntegrationType")
//            .value("mps_straight", ExudationModel::IntegrationType::mps_straight)
//            .value("mps", ExudationModel::IntegrationType::mps)
//            .value("mls", ExudationModel::IntegrationType::mls)
//            ;

}


} // end namespace CRootBox

#endif /* PY_ROOTBOX_H_ */
