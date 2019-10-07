// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "external/pybind11/include/pybind11/pybind11.h"
#include "external/pybind11/include/pybind11/stl.h"
namespace py = pybind11;

/**
 * A Python binding based on pybind11
 */
#include "mymath.h"
#include "sdf.h"
#include "organparameter.h"
#include "Organ.h"
#include "Organism.h"
#include "soil.h"
#include "tropism.h"

#include "rootparameter.h"
#include "seedparameter.h"
#include "leafparameter.h"
#include "stemparameter.h"
#include "Root.h"
#include "Seed.h"
#include "Leaf.h"
#include "Stem.h"

#include "RootSystem.h"
#include "Plant.h"

#include "sdf_rs.h" // todo to revise ...

namespace CPlantBox {

/**
 * Trampoline classes
 * are required for all base classes that can be derived in Python (TODO)
 */
class PyTropism : public Tropism {
public:

    using Tropism::Tropism; /* Inherit the constructors */

    std::shared_ptr<Tropism> copy(std::shared_ptr<Organism> plant) override
        { PYBIND11_OVERLOAD( std::shared_ptr<Tropism>, Tropism, copy, plant); }

    double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const std::shared_ptr<Organ> organ) override
        { PYBIND11_OVERLOAD( double, Tropism, tropismObjective, pos, old, a, b, dx, organ ); }

    // TODO do all virtuals

};





/**
 * plantbox
 */
PYBIND11_MODULE(plantbox, m) {
    /*
     * mymath
     */
    py::class_<Vector2i>(m, "Vector2i", py::buffer_protocol())
        .def(py::init<>())
        .def(py::init<int, int>())
        .def(py::init<const Vector2i&>())
        .def("__str__", &Vector2i::toString)
        .def_readwrite("x", &Vector2i::x)
        .def_readwrite("y", &Vector2i::y)
        .def_buffer([](Vector2i &v) -> py::buffer_info { /* enables numpy conversion with np.array(vector2i_instance, copy = False) */
            return py::buffer_info(
                    &v.x,                                    /* Pointer to buffer */
                    sizeof(int),                            /* Size of one scalar */
                    py::format_descriptor<int>::format(),   /* Python struct-style format descriptor */
                    1,                                       /* Number of dimensions */
                    { 2 },                                   /* Buffer dimensions */
                    { sizeof(int) }                         /* Strides (in bytes) for each index */
                );
        });
    py::class_<Vector2d>(m, "Vector2d", py::buffer_protocol())
        .def(py::init<>())
        .def(py::init<double, double>())
        .def(py::init<const Vector2d&>())
        .def("__str__", &Vector2d::toString)
        .def_readwrite("x", &Vector2d::x)
        .def_readwrite("y", &Vector2d::y)
        .def_buffer([](Vector2d &v) -> py::buffer_info { /* enables numpy conversion with np.array(vector2d_instance, copy = False) */
            return py::buffer_info(
                &v.x,                                    /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                1,                                       /* Number of dimensions */
                { 2 },                                   /* Buffer dimensions */
                { sizeof(double) }                       /* Strides (in bytes) for each index */
            );
        });
    py::class_<Vector3d>(m, "Vector3d", py::buffer_protocol())
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def(py::init<const Vector3d&>())
        .def("rotAB", &Vector3d::rotAB)
        .def("normalize", &Vector3d::normalize)
        .def("times", (double (Vector3d::*)(const Vector3d&) const) &Vector3d::times) // overloads
        .def("length", &Vector3d::length)
        .def("times", (Vector3d (Vector3d::*)(const double) const) &Vector3d::times) // overloads
        .def("plus", &Vector3d::plus)
        .def("minus", &Vector3d::minus)
        .def("cross", &Vector3d::cross)
        .def("__str__", &Vector3d::toString)
        .def_readwrite("x", &Vector3d::x)
        .def_readwrite("y", &Vector3d::y)
        .def_readwrite("z", &Vector3d::z)
        .def_buffer([](Vector3d &v) -> py::buffer_info { /* enables numpy conversion with np.array(vector3d_instance, copy = False) */
            return py::buffer_info(
                &v.x,                                    /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                1,                                       /* Number of dimensions */
                { 3 },                                   /* Buffer dimensions */
                { sizeof(double) }                       /* Strides (in bytes) for each index */
            );
        });
    py::class_<Matrix3d>(m, "Matrix3d", py::buffer_protocol())
        .def(py::init<>())
        .def(py::init<double, double, double, double, double, double, double, double, double>())
        .def(py::init<const Vector3d&, const Vector3d&, const Vector3d&>())
        .def(py::init<const Matrix3d&>())
        .def("rotX", &Matrix3d::rotX)
        .def("rotY", &Matrix3d::rotY)
        .def("rotZ", &Matrix3d::rotZ)
        .def("ons", &Matrix3d::ons)
        .def("det", &Matrix3d::det)
        .def("inverse", &Matrix3d::inverse)
        .def("column", &Matrix3d::column)
        .def("row", &Matrix3d::row)
        .def("times",(void (Matrix3d::*)(const Matrix3d&)) &Matrix3d::times) // overloads
        .def("times",(Vector3d (Matrix3d::*)(const Vector3d&) const) &Matrix3d::times) // overloads
        .def("__str__", &Matrix3d::toString)
        .def_readwrite("r0", &Matrix3d::r0)
        .def_readwrite("r1", &Matrix3d::r1)
        .def_readwrite("r2", &Matrix3d::r2)
        .def_buffer([](Matrix3d &m_) -> py::buffer_info { /* enables numpy conversion with np.array(matrix3d_instance, copy = False) */
            return py::buffer_info(
                &m_.r0.x,                               /* Pointer to buffer */
                sizeof(float),                          /* Size of one scalar */
                py::format_descriptor<double>::format(),/* Python struct-style format descriptor */
                2,                                      /* Number of dimensions */
                { 3, 3 },                               /* Buffer dimensions */
                { sizeof(double) * 3,  sizeof(double) } /* Strides (in bytes) for each index */
                );
            });
    /*
     * sdf
     */
    py::class_<SignedDistanceFunction, std::shared_ptr<SignedDistanceFunction>>(m,"SignedDistanceFunction")
        .def(py::init<>())
        .def("getDist",&SignedDistanceFunction::getDist)
        .def("writePVPScript", (std::string (SignedDistanceFunction::*)() const) &SignedDistanceFunction::writePVPScript) // overloads
        .def("getGradient",  &SignedDistanceFunction::getGradient, py::arg("p"), py::arg("eps") = 5.e-4) // defaults
        .def("__str__",&SignedDistanceFunction::toString);
    py::class_<SDF_PlantBox, SignedDistanceFunction, std::shared_ptr<SDF_PlantBox>>(m, "SDF_PlantBox")
        .def(py::init<double,double,double>());
    py::class_<SDF_PlantContainer, SignedDistanceFunction, std::shared_ptr<SDF_PlantContainer>>(m,"SDF_PlantContainer")
        .def(py::init<>())
        .def(py::init<double,double,double,double>());
    py::class_<SDF_RotateTranslate, SignedDistanceFunction, std::shared_ptr<SDF_RotateTranslate>>(m, "SDF_RotateTranslate")
        .def(py::init<SignedDistanceFunction*,double,int,Vector3d&>())
        .def(py::init<SignedDistanceFunction*,Vector3d&>());
    py::enum_<SDF_RotateTranslate::SDF_Axes>(m, "SDF_Axis")
        .value("xaxis", SDF_RotateTranslate::SDF_Axes::xaxis)
        .value("yaxis", SDF_RotateTranslate::SDF_Axes::yaxis)
        .value("zaxis", SDF_RotateTranslate::SDF_Axes::zaxis)
        .export_values();
    py::class_<SDF_Intersection, SignedDistanceFunction, std::shared_ptr<SDF_Intersection>>(m,"SDF_Intersection")
        .def(py::init<std::vector<SignedDistanceFunction*>>())
        .def(py::init<SignedDistanceFunction*,SignedDistanceFunction*>());
    py::class_<SDF_Union, SDF_Intersection, std::shared_ptr<SDF_Union>>(m, "SDF_Union")
        .def(py::init<std::vector<SignedDistanceFunction*>>())
        .def(py::init<SignedDistanceFunction*,SignedDistanceFunction*>());
    py::class_<SDF_Difference, SDF_Intersection, std::shared_ptr<SDF_Difference>>(m, "SDF_Difference")
        .def(py::init<std::vector<SignedDistanceFunction*>>())
        .def(py::init<SignedDistanceFunction*,SignedDistanceFunction*>());
    py::class_<SDF_Complement, SignedDistanceFunction, std::shared_ptr<SDF_Complement>>(m, "SDF_Complement")
        .def(py::init<SignedDistanceFunction*>());
    py::class_<SDF_HalfPlane, SignedDistanceFunction, std::shared_ptr<SDF_HalfPlane>>(m, "SDF_HalfPlane")
        .def(py::init<Vector3d&,Vector3d&>())
        .def(py::init<Vector3d&,Vector3d&,Vector3d&>())
        .def_readwrite("o", &SDF_HalfPlane::o)
        .def_readwrite("n", &SDF_HalfPlane::n)
        .def_readwrite("p1", &SDF_HalfPlane::p1)
        .def_readwrite("p2", &SDF_HalfPlane::p2);
    /*
     * organparameter.h
     */
    py::class_<OrganSpecificParameter, std::shared_ptr<OrganSpecificParameter>>(m,"OrganSpecificParameter")
        .def(py::init<>())
        .def_readwrite("subType",&OrganSpecificParameter::subType)
        .def("__str__",&OrganSpecificParameter::toString);
    py::class_<OrganRandomParameter, std::shared_ptr<OrganRandomParameter>>(m,"OrganRandomParameter")
        .def(py::init<std::shared_ptr<Organism>>())
        .def("copy",&OrganRandomParameter::copy)
        .def("realize",&OrganRandomParameter::realize)
        .def("getParameter",&OrganRandomParameter::getParameter)
        .def("__str__",&OrganRandomParameter::toString, py::arg("verbose") = true) // default
        .def("writeXML",(void (OrganRandomParameter::*)(std::string name) const) &OrganRandomParameter::writeXML) // overloads
        .def("readXML", (void (OrganRandomParameter::*)(std::string name)) &OrganRandomParameter::readXML) // overloads
        .def("bindIntParameter", (void (OrganRandomParameter::*)(std::string, int*, std::string, double*))  &OrganRandomParameter::bindParameter,
            py::arg("name"), py::arg("i"), py::arg("descr") = "", py::arg("dev") = (double*) nullptr) // overloads, defaults
        .def("bindDoubleParameter", (void (OrganRandomParameter::*)(std::string, double*, std::string, double*))  &OrganRandomParameter::bindParameter,
            py::arg("name"), py::arg("i"), py::arg("descr") = "", py::arg("dev") = (double*) nullptr) // overloads, defaults
        .def_readwrite("name",&OrganRandomParameter::name)
        .def_readwrite("organType",&OrganRandomParameter::organType)
        .def_readwrite("subType",&OrganRandomParameter::subType)
        .def_readwrite("plant",&OrganRandomParameter::plant);
    /**
     * Organ.h
     */
    py::class_<Organ, std::shared_ptr<Organ>>(m, "Organ")
        .def(py::init<std::shared_ptr<Organism>, std::shared_ptr<Organ>, int, int, double>())
        .def(py::init<int, std::shared_ptr<const OrganSpecificParameter>, bool, bool, double, double, bool, int>())
        .def("copy",&Organ::copy)
        .def("organType",&Organ::organType)
        .def("simulate",&Organ::simulate,py::arg("dt"), py::arg("verbose") = bool(false) ) // default

        .def("setParent",&Organ::setParent)
        .def("getParent",&Organ::getParent)
        .def("setOrganism",&Organ::setOrganism)
        .def("addChild",&Organ::addChild)

        .def("getId",&Organ::getId)
        .def("getParam",&Organ::getParam)
        .def("getOrganRandomParameter",&Organ::getOrganRandomParameter)
        .def("isAlive",&Organ::isAlive)
        .def("isActive",&Organ::isActive)
        .def("getAge",&Organ::getAge)
        .def("getLength",&Organ::getLength)

        .def("getNumberOfNodes",&Organ::getNumberOfNodes)
        .def("getNode",&Organ::getNode)
        .def("getNodeId",&Organ::getNodeId)
        .def("getNodeCT",&Organ::getNodeCT)
        .def("addNode",(void (Organ::*)(Vector3d n, double t)) &Organ::addNode) // overloads
        .def("addNode",(void (Organ::*)(Vector3d n, int id, double t)) &Organ::addNode) // overloads
        .def("getSegments",&Organ::getSegments)

        .def("hasMoved",&Organ::hasMoved)
        .def("getOldNumberOfNodes",&Organ::getOldNumberOfNodes)

        .def("getOrgans", (std::vector<std::shared_ptr<Organ>> (Organ::*)(int otype)) &Organ::getOrgans, py::arg("ot")=-1) //overloads, default
        .def("getOrgans", (void (Organ::*)(int otype, std::vector<std::shared_ptr<Organ>>& v)) &Organ::getOrgans)
        .def("getParameter",&Organ::getParameter)
        .def("__str__",&Organ::toString);
    /*
     * Organism.h
     */
    py::class_<Organism, std::shared_ptr<Organism>>(m, "Organism")
        .def(py::init<>())
        .def("copy", &Organism::copy)
        .def("organTypeNumber", &Organism::organTypeNumber)
        .def("organTypeName", &Organism::organTypeName)
        .def("getOrganRandomParameter", (std::shared_ptr<OrganRandomParameter> (Organism::*)(int, int) const)  &Organism::getOrganRandomParameter) //overloads
        .def("getOrganRandomParameter", (std::vector<std::shared_ptr<OrganRandomParameter>> (Organism::*)(int) const) &Organism::getOrganRandomParameter) //overloads
        .def("setOrganRandomParameter", &Organism::setOrganRandomParameter)

        .def("addOrgan", &Organism::addOrgan)
        .def("initialize", &Organism::initialize)
        .def("simulate", &Organism::simulate, py::arg("dt"), py::arg("verbose") = false) //default
        .def("getSimTime", &Organism::getSimTime)

        .def("getOrgans", &Organism::getOrgans, py::arg("ot") = -1) // default
        .def("getParameter", &Organism::getParameter, py::arg("name"), py::arg("ot") = -1, py::arg("organs") = std::vector<std::shared_ptr<Organ>>(0)) // default
        .def("getSummed", &Organism::getSummed, py::arg("name"), py::arg("ot") = -1) // default

        .def("getNumberOfOrgans", &Organism::getNumberOfOrgans)
        .def("getNumberOfNodes", &Organism::getNumberOfNodes)
        .def("getNumberOfSegments", &Organism::getNumberOfSegments, py::arg("ot") = -1) // default
        .def("getPolylines", &Organism::getPolylines, py::arg("ot") = -1) // default
        .def("getPolylineCTs", &Organism::getPolylineCTs, py::arg("ot") = -1) // default
        .def("getNodes", &Organism::getNodes)
        .def("getNodeCTs", &Organism::getNodeCTs)
        .def("getSegments", &Organism::getSegments, py::arg("ot") = -1) // default
        .def("getSegmentCTs", &Organism::getSegmentCTs, py::arg("ot") = -1) // default
        .def("getSegmentOrigins", &Organism::getSegmentOrigins, py::arg("ot") = -1) // default

        .def("getNumberOfNewNodes", &Organism::getNumberOfNewNodes)
        .def("getNumberOfNewOrgans", &Organism::getNumberOfNewOrgans)
        .def("getUpdatedNodeIndices", &Organism::getUpdatedNodeIndices)
        .def("getUpdatedNodes", &Organism::getUpdatedNodes)
        .def("getUpdatedNodeCTs", &Organism::getUpdatedNodeCTs)

        .def("getNewNodes", &Organism::getNewNodes)
        .def("getNewNodeCTs", &Organism::getNewNodeCTs)
        .def("getNewSegments", &Organism::getNewSegments, py::arg("ot") = -1)  // default
        .def("getNewSegmentOrigins", &Organism::getNewSegmentOrigins, py::arg("ot") = -1)  // default

        .def("initializeReader", &Organism::initializeReader)
        .def("readParameters", &Organism::readParameters, py::arg("name"), py::arg("basetag") = "plant")  // default
        .def("writeParameters", &Organism::writeParameters, py::arg("name"), py::arg("basetag") = "plant", py::arg("comments") = true)  // default
        .def("writeRSML", &Organism::writeRSML)
        .def("getRSMLSkip", &Organism::getRSMLSkip)
        .def("setRSMLSkip", &Organism::setRSMLSkip)
        .def("getRSMLProperties", &Organism::getRSMLProperties) //todo policy

        .def("getOrganIndex", &Organism::getOrganIndex)
        .def("getNodeIndex", &Organism::getNodeIndex)

        .def("setSeed", &Organism::setSeed)
        .def("rand", &Organism::rand)
        .def("randn", &Organism::randn)
        .def("__str__",&Organism::toString);

    py::enum_<Organism::OrganTypes>(m, "OrganTypes")
        .value("organ", Organism::OrganTypes::ot_organ)
        .value("seed", Organism::OrganTypes::ot_seed)
        .value("root", Organism::OrganTypes::ot_root)
        .value("stem", Organism::OrganTypes::ot_stem)
        .value("leaf", Organism::OrganTypes::ot_leaf)
        .export_values();
    /*
     * soil
     */
    py::class_<SoilLookUp, std::shared_ptr<SoilLookUp>>(m, "SoilLookUp")   /// <--- TODO WRAP (?)
        .def(py::init<>())
        .def("getValue",&SoilLookUp::getValue) /// <--- TODO defaults
        .def("__str__",&SoilLookUp::toString);
    py::class_<SoilLookUpSDF, SoilLookUp, std::shared_ptr<SoilLookUpSDF>>(m,"SoilLookUpSDF")
        .def(py::init<>())
        .def(py::init<SignedDistanceFunction*, double, double, double>())
        .def_readwrite("sdf", &SoilLookUpSDF::sdf)
        .def_readwrite("fmax", &SoilLookUpSDF::fmax)
        .def_readwrite("fmin", &SoilLookUpSDF::fmin)
        .def_readwrite("slope", &SoilLookUpSDF::slope);
    py::class_<MultiplySoilLookUps, SoilLookUp, std::shared_ptr<MultiplySoilLookUps>>(m, "MultiplySoilLookUps")
        .def(py::init<SoilLookUp*, SoilLookUp*>())
        .def(py::init<std::vector<SoilLookUp*>>());
    py::class_<ProportionalElongation, SoilLookUp, std::shared_ptr<ProportionalElongation>>(m, "ProportionalElongation")
        .def(py::init<>())
        .def("setScale", &ProportionalElongation::setScale)
        .def("setBaseLookUp", &ProportionalElongation::setBaseLookUp)
        .def("__str__",&ProportionalElongation::toString);
    py::class_<Grid1D, SoilLookUp, std::shared_ptr<Grid1D>>(m, "Grid1D")
        .def(py::init<>())
        .def(py::init<size_t, std::vector<double>, std::vector<double>>())
        .def("map",&Grid1D::map)
        .def_readwrite("n", &Grid1D::n)
        .def_readwrite("grid", &Grid1D::grid)
        .def_readwrite("data", &Grid1D::data);
    py::class_<EquidistantGrid1D, Grid1D, std::shared_ptr<EquidistantGrid1D>>(m, "EquidistantGrid1D")
        .def(py::init<double, double, size_t>())
        .def(py::init<double, double, std::vector<double>>())
        .def("map",&EquidistantGrid1D::map)
        .def_readwrite("n", &EquidistantGrid1D::n)
        .def_readwrite("grid", &EquidistantGrid1D::grid)
        .def_readwrite("data", &EquidistantGrid1D::data);
    /**
     * tropism.h
     */
    py::class_<Tropism, PyTropism, std::shared_ptr<Tropism>>(m, "Tropism")
        .def(py::init<std::shared_ptr<Organism>>())
        .def(py::init<std::shared_ptr<Organism>, double, double>())
        .def("copy",&Tropism::copy) // todo policy
        .def("setGeometry",&Tropism::setGeometry)
        .def("setTropismParameter",&Tropism::setTropismParameter)
        .def("getHeading",&Tropism::getHeading)
        .def("getUCHeading",&Tropism::getUCHeading)
        .def("tropismObjective",&Tropism::tropismObjective)
        .def("getPosition",&Tropism::getPosition);
    py::class_<Gravitropism, Tropism, std::shared_ptr<Gravitropism>>(m, "Gravitropism")
        .def(py::init<std::shared_ptr<Organism>, double, double>());
    py::class_<Plagiotropism, Tropism, std::shared_ptr<Plagiotropism>>(m, "Plagiotropism")
        .def(py::init<std::shared_ptr<Organism>,double, double>());
    py::class_<Exotropism, Tropism, std::shared_ptr<Exotropism>>(m, "Exotropism")
        .def(py::init<std::shared_ptr<Organism>,double, double>());
    py::class_<Hydrotropism, Tropism, std::shared_ptr<Hydrotropism>>(m, "Hydrotropism")
        .def(py::init<std::shared_ptr<Organism>,double, double, SoilLookUp*>());
//    py::class_<CombinedTropism, Tropism>(m, "CombinedTropism") // Todo constructors needs some extra work (?)
//        .def(py::init<>());
    // todo antigravi, twist ...
   /*
    * analysis.h
    */
    py::class_<SegmentAnalyser, std::shared_ptr<SegmentAnalyser>>(m, "SegmentAnalyser")
       .def(py::init<>())
       .def(py::init<Organism&>())
       .def(py::init<SegmentAnalyser&>())
       .def("addSegments",(void (SegmentAnalyser::*)(const Organism&)) &SegmentAnalyser::addSegments) //overloads
       .def("addSegments",(void (SegmentAnalyser::*)(const SegmentAnalyser&)) &SegmentAnalyser::addSegments) //overloads
       .def("crop", &SegmentAnalyser::crop)
       .def("filter", (void (SegmentAnalyser::*)(std::string, double, double)) &SegmentAnalyser::filter) //overloads
       .def("filter", (void (SegmentAnalyser::*)(std::string, double)) &SegmentAnalyser::filter) //overloads
       .def("pack", &SegmentAnalyser::pack)
       .def("getParameter", &SegmentAnalyser::getParameter)
       .def("getSegmentLength", &SegmentAnalyser::getSegmentLength)
       .def("getSummed", (double (SegmentAnalyser::*)(std::string) const) &SegmentAnalyser::getSummed) //overloads
       .def("getSummed", (double (SegmentAnalyser::*)(std::string, SignedDistanceFunction*) const) &SegmentAnalyser::getSummed) //overloads
       .def("distribution", (std::vector<double> (SegmentAnalyser::*)(std::string, double, double, int, bool) const) &SegmentAnalyser::distribution) //overloads
       .def("distribution", (std::vector<SegmentAnalyser> (SegmentAnalyser::*)(double, double, int) const) &SegmentAnalyser::distribution) //overloads
       .def("distribution2", (std::vector<std::vector<double>> (SegmentAnalyser::*)(std::string, double, double, double, double, int, int, bool) const) &SegmentAnalyser::distribution2) //overloads
       .def("distribution2", (std::vector<std::vector<SegmentAnalyser>> (SegmentAnalyser::*)(double, double, double, double, int, int) const) &SegmentAnalyser::distribution2) //overloads
       .def("getOrgans", &SegmentAnalyser::getOrgans)
       .def("getNumberOfOrgans", &SegmentAnalyser::getNumberOfOrgans)
       .def("cut", (SegmentAnalyser (SegmentAnalyser::*)(const SDF_HalfPlane&) const) &SegmentAnalyser::cut)
       .def("addUserData", &SegmentAnalyser::addUserData)
       .def("clearUserData", &SegmentAnalyser::clearUserData)
       .def("write", &SegmentAnalyser::write)
       .def_readwrite("nodes", &SegmentAnalyser::nodes)
       .def_readwrite("segments", &SegmentAnalyser::segments)
       .def_readwrite("segCTs", &SegmentAnalyser::segCTs)
       .def_readwrite("segO", &SegmentAnalyser::segO);
    /*
     * rootparameter.h
     */
    py::class_<RootRandomParameter, OrganRandomParameter, std::shared_ptr<RootRandomParameter>>(m, "RootRandomParameter")
        .def(py::init<std::shared_ptr<Organism>>())
        .def("getLateralType",&RootRandomParameter::getLateralType)
        .def("getK",&RootRandomParameter::getK)
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
        .def_readwrite("successor", &RootRandomParameter::successor)
        .def_readwrite("successorP", &RootRandomParameter::successorP)
        .def_readwrite("f_tf", &RootRandomParameter::f_tf)
        .def_readwrite("f_gf", &RootRandomParameter::f_gf)
        .def_readwrite("f_se", &RootRandomParameter::f_se)
        .def_readwrite("f_sa", &RootRandomParameter::f_sa)
        .def_readwrite("f_sbp", &RootRandomParameter::f_sbp);
    py::class_<RootSpecificParameter, OrganSpecificParameter, std::shared_ptr<RootSpecificParameter>>(m, "RootSpecificParameter")
        .def(py::init<>())
        .def(py::init<int , double, double, const std::vector<double>&, int, double, double, double, double>())
        .def_readwrite("lb", &RootSpecificParameter::lb)
        .def_readwrite("la", &RootSpecificParameter::la)
        .def_readwrite("ln", &RootSpecificParameter::ln)
        .def_readwrite("nob", &RootSpecificParameter::nob)
        .def_readwrite("r", &RootSpecificParameter::r)
        .def_readwrite("a", &RootSpecificParameter::a)
        .def_readwrite("theta", &RootSpecificParameter::theta)
        .def_readwrite("rlt", &RootSpecificParameter::rlt)
        .def("getK",&RootSpecificParameter::getK);
    /*
     * seedparameter.h
     */
    py::class_<SeedRandomParameter, OrganRandomParameter, std::shared_ptr<SeedRandomParameter>>(m, "SeedRandomParameter")
        .def(py::init<std::shared_ptr<Organism>>())
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
        .def_readwrite("simtime", &SeedRandomParameter::simtime);
    py::class_<SeedSpecificParameter, OrganSpecificParameter, std::shared_ptr<SeedSpecificParameter>>(m, "SeedSpecificParameter")
        .def(py::init<>())
        .def(py::init<int, Vector3d , double, int, int, int, double, double, double, double, int, double>())
        .def_readwrite("seedPos", &SeedSpecificParameter::seedPos)
        .def_readwrite("firstB", &SeedSpecificParameter::firstB)
        .def_readwrite("delayB", &SeedSpecificParameter::delayB)
        .def_readwrite("maxB", &SeedSpecificParameter::maxB)
        .def_readwrite("nC", &SeedSpecificParameter::nC)
        .def_readwrite("firstSB", &SeedSpecificParameter::firstSB)
        .def_readwrite("delaySB", &SeedSpecificParameter::delaySB)
        .def_readwrite("delayRC", &SeedSpecificParameter::delayRC)
        .def_readwrite("nz", &SeedSpecificParameter::nz)
        .def_readwrite("maxTil", &SeedSpecificParameter::maxTil);
    /*
     * leafparameter.h
     */
    py::class_<LeafRandomParameter, OrganRandomParameter, std::shared_ptr<LeafRandomParameter>>(m, "LeafRandomParameter")
        .def(py::init<std::shared_ptr<Organism>>())
        .def("getLateralType",&LeafRandomParameter::getLateralType)
        .def("getK",&LeafRandomParameter::getK)
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
        .def_readwrite("successor", &LeafRandomParameter::successor)
        .def_readwrite("successorP", &LeafRandomParameter::successorP)
        .def_readwrite("f_tf", &LeafRandomParameter::f_tf)
        .def_readwrite("f_gf", &LeafRandomParameter::f_gf)
        .def_readwrite("f_se", &LeafRandomParameter::f_se)
        .def_readwrite("f_sa", &LeafRandomParameter::f_sa)
        .def_readwrite("f_sbp", &LeafRandomParameter::f_sbp);
    py::class_<LeafSpecificParameter, OrganSpecificParameter, std::shared_ptr<LeafSpecificParameter>>(m, "LeafSpecificParameter")
        .def(py::init<>())
        .def(py::init<int , double, double, const std::vector<double>&, double, double, double, double>())
        .def_readwrite("lb", &LeafSpecificParameter::lb)
        .def_readwrite("la", &LeafSpecificParameter::la)
        .def_readwrite("ln", &LeafSpecificParameter::ln)
        .def_readwrite("r", &LeafSpecificParameter::r)
        .def_readwrite("a", &LeafSpecificParameter::a)
        .def_readwrite("theta", &LeafSpecificParameter::theta)
        .def_readwrite("rlt", &LeafSpecificParameter::rlt)
        .def("getK",&LeafSpecificParameter::getK);
    /*
     * stemparameter.h
     */
    py::class_<StemRandomParameter, OrganRandomParameter, std::shared_ptr<StemRandomParameter>>(m, "StemRandomParameter")
        .def(py::init<std::shared_ptr<Organism>>())
        .def("getLateralType",&StemRandomParameter::getLateralType)
        .def("getK",&StemRandomParameter::getK)
        .def_readwrite("lb", &StemRandomParameter::lb)
        .def_readwrite("lbs", &StemRandomParameter::lbs)
        .def_readwrite("la", &StemRandomParameter::la)
        .def_readwrite("las", &StemRandomParameter::las)
        .def_readwrite("k", &StemRandomParameter::k)
        .def_readwrite("ks", &StemRandomParameter::ks)
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
        .def_readwrite("successor", &StemRandomParameter::successor)
        .def_readwrite("successorP", &StemRandomParameter::successorP)
        .def_readwrite("f_tf", &StemRandomParameter::f_tf)
        .def_readwrite("f_gf", &StemRandomParameter::f_gf)
        .def_readwrite("f_se", &StemRandomParameter::f_se)
        .def_readwrite("f_sa", &StemRandomParameter::f_sa)
        .def_readwrite("f_sbp", &StemRandomParameter::f_sbp);
    py::class_<StemSpecificParameter, OrganSpecificParameter, std::shared_ptr<StemSpecificParameter>>(m, "StemSpecificParameter")
        .def(py::init<>())
        .def(py::init<int , double, double, const std::vector<double>&, double, double, double, double, double>())
        .def_readwrite("lb", &StemSpecificParameter::lb)
        .def_readwrite("la", &StemSpecificParameter::la)
        .def_readwrite("ln", &StemSpecificParameter::ln)
        .def_readwrite("r", &StemSpecificParameter::r)
        .def_readwrite("a", &StemSpecificParameter::a)
        .def_readwrite("theta", &StemSpecificParameter::theta)
        .def_readwrite("rlt", &StemSpecificParameter::rlt)
        .def("getK",&StemSpecificParameter::getK);
    /**
      * Root.h
      */
    py::class_<Root, Organ, std::shared_ptr<Root>>(m, "Root")
        .def(py::init<std::shared_ptr<Organism>, int, Vector3d, double, std::shared_ptr<Organ>, double, int>())
        .def(py::init<int, std::shared_ptr<OrganSpecificParameter>, bool, bool, double, double, Vector3d, double, int, bool, int>())
        .def("calcCreationTime", &Root::calcCreationTime)
        .def("calcLength", &Root::calcLength)
        .def("calcAge", &Root::calcAge)
        .def("getRootRandomParameter", &Root::getRootRandomParameter)
        .def("param", &Root::param)
        .def("dx", &Root::dx)
        .def_readwrite("iHeading", &Root::iHeading)
        .def_readwrite("parentBaseLength", &Root::parentBaseLength)
        .def_readwrite("parentNI", &Root::parentNI);
     /**
      * Seed.h
      */
    py::class_<Seed, Organ, std::shared_ptr<Seed>>(m, "Seed")
        .def(py::init<std::shared_ptr<Organism>>())
        .def("initialize", &Seed::initialize)
        .def("param", &Seed::param)
        .def("getNumberOfRootCrowns", &Seed::getNumberOfRootCrowns)
        .def("baseOrgans", &Seed::baseOrgans)
        .def("copyBaseOrgans", &Seed::copyBaseOrgans)
        .def("createRoot", &Seed::createRoot)
        .def("createStem", &Seed::createStem)
        .def_readwrite("tapRootType", &Seed::tapRootType)
        .def_readwrite("basalType", &Seed::basalType)
        .def_readwrite("shootborneType", &Seed::shootborneType)
        .def_readwrite("mainStemType", &Seed::mainStemType)
        .def_readwrite("tillerType", &Seed::tillerType);
    /**
     * Leaf.h
     */
    py::class_<Leaf, Organ, std::shared_ptr<Leaf>>(m, "Leaf")
        .def(py::init<std::shared_ptr<Organism>, int, Vector3d, double, std::shared_ptr<Organ>, double, int>())
        .def(py::init<int, std::shared_ptr<OrganSpecificParameter>, bool, bool, double, double, Vector3d, double, int, bool, int>())
        .def("calcCreationTime", &Leaf::calcCreationTime)
        .def("calcLength", &Leaf::calcLength)
        .def("calcAge", &Leaf::calcAge)
        .def("getLeafRandomParameter", &Leaf::getLeafRandomParameter)
        .def("param", &Leaf::param)
        .def("dx", &Leaf::dx)
        .def_readwrite("iHeading", &Leaf::iHeading)
        .def_readwrite("parentBaseLength", &Leaf::parentBaseLength)
        .def_readwrite("parentNI", &Leaf::parentNI);
   /**
    * Stem.h
    */
   py::class_<Stem, Organ, std::shared_ptr<Stem>>(m, "Stem")
       .def(py::init<std::shared_ptr<Organism>, int, Vector3d, double, std::shared_ptr<Organ>, double, int>())
       .def(py::init<int, std::shared_ptr<OrganSpecificParameter>, bool, bool, double, double, Vector3d, double, int, bool, int>())
       .def("calcCreationTime", &Stem::calcCreationTime)
       .def("calcLength", &Stem::calcLength)
       .def("calcAge", &Stem::calcAge)
       .def("getStemRandomParameter", &Stem::getStemRandomParameter)
       .def("param", &Stem::param)
       .def("dx", &Stem::dx)
       .def_readwrite("iHeading", &Stem::iHeading)
       .def_readwrite("parentBaseLength", &Stem::parentBaseLength)
       .def_readwrite("parentNI", &Stem::parentNI);
    /*
     * RootSystem.h
     */
    py::class_<RootSystem, Organism, std::shared_ptr<RootSystem>>(m, "RootSystem")
        .def(py::init<>())
        .def("getRootRandomParameter", (std::shared_ptr<RootRandomParameter> (RootSystem::*)(int) const) &RootSystem::getRootRandomParameter)
        .def("getRootRandomParameter", (std::vector<std::shared_ptr<RootRandomParameter>> (RootSystem::*)() const) &RootSystem::getRootRandomParameter)
        .def("setRootSystemParameter", &RootSystem::setRootSystemParameter)
        .def("getRootSystemParameter", &RootSystem::getRootSystemParameter)
        .def("openFile", &RootSystem::openFile, py::arg("filename"), py::arg("subdir") = "modelparameter/")
        .def("setGeometry", &RootSystem::setGeometry)
        .def("setSoil", &RootSystem::setSoil)
        .def("reset", &RootSystem::reset)
        .def("initialize", (void (RootSystem::*)()) &RootSystem::initialize)
        .def("initialize", (void (RootSystem::*)(int,int)) &RootSystem::initialize)
        .def("setTropism", &RootSystem::setTropism)
        .def("simulate",(void (RootSystem::*)(double,bool)) &RootSystem::simulate, py::arg("dt"), py::arg("verbose") = false)
        .def("simulate",(void (RootSystem::*)()) &RootSystem::simulate)
        .def("simulate",(void (RootSystem::*)(double, double, ProportionalElongation*, bool)) &RootSystem::simulate)
        .def("getRoots", &RootSystem::getRoots)
        .def("initCallbacks", &RootSystem::initCallbacks)
        .def("createTropismFunction", &RootSystem::createTropismFunction)
        .def("createGrowthFunction", &RootSystem::createGrowthFunction)
        .def("getNumberOfRoots", &RootSystem::getNumberOfRoots, py::arg("all") = false)
        .def("getBaseRoots", &RootSystem::getBaseRoots)
        .def("getShootSegments", &RootSystem::getShootSegments)
        .def("getRootTips", &RootSystem::getRootTips)
        .def("getRootBases", &RootSystem::getRootBases)
        .def("push",&RootSystem::push)
        .def("pop",&RootSystem::pop)
        .def("write", &RootSystem::write);
     /*
      * Plant.h
      */
     py::class_<Plant, Organism, std::shared_ptr<Plant>>(m, "Plant")
        .def(py::init<>())
        .def("setGeometry", &Plant::setGeometry)
        .def("reset", &Plant::reset)
        .def("setTropism", &Plant::setTropism)
        .def("simulate",(void (Plant::*)(double,bool)) &Plant::simulate, py::arg("dt"), py::arg("verbose") = false) // in Organism::simulate
        .def("simulate",(void (Plant::*)()) &Plant::simulate)
        .def("initCallbacks", &Plant::initCallbacks)
        .def("createTropismFunction", &Plant::createTropismFunction)
        .def("createGrowthFunction", &Plant::createGrowthFunction);
     py::enum_<Plant::TropismTypes>(m, "TropismType")
        .value("plagio", Plant::TropismTypes::tt_plagio)
        .value("gravi", Plant::TropismTypes::tt_gravi)
        .value("exo", Plant::TropismTypes::tt_exo)
        .value("hydro", Plant::TropismTypes::tt_hydro)
        .value("antigravi", Plant::TropismTypes::tt_antigravi)
        .value("twist", Plant::TropismTypes::tt_twist)
        .export_values();
     py::enum_<Plant::GrowthFunctionTypes>(m, "GrowthFunctionType")
         .value("negexp", Plant::GrowthFunctionTypes::gft_negexp)
         .value("linear", Plant::GrowthFunctionTypes::gft_linear)
         .export_values();

    //   /*
    //    * sdf_rs.h todo
    //    */
    //   py::class_<SDF_RootSystem, SignedDistanceFunction>(m, "SDF_RootSystem")
    //       .def(py::init<std::vector<Vector3d>, std::vector<Vector2i>, std::vector<double>, double>())
    //       .def(py::init<Root&, double>())
    //       .def(py::init<Organism&, double>());
}

}
