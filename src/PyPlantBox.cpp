// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "external/pybind11/include/pybind11/pybind11.h"
#include "external/pybind11/include/pybind11/stl.h"
#include "external/pybind11/include/pybind11/functional.h"
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
#include "MappedOrganism.h"

// functional
#include "Perirhizal.h"
#include "XylemFlux.h"
#include "ExudationModel.h"
#include "Photosynthesis.h"
#include "PiafMunch/runPM.h"

#include "PlantHydraulicParameters.h"
#include "PlantHydraulicModel.h"

// visualisation
#include "Quaternion.h"
#include "CatmullRomSpline.h"
#include "PlantVisualiser.h"

#include "sdf_rs.h" // todo to revise ...

namespace CPlantBox {

/**
 * Trampoline classes
 *
 * Are required for all base classes that can be derived in Python, currently
 * Tropism
 * SoilLookUp
 */
class PyTropism : public Tropism {
public:

    using Tropism::Tropism; /* Inherit the constructors */

    std::shared_ptr<Tropism> copy(std::shared_ptr<Organism> plant) override
    { PYBIND11_OVERLOAD( std::shared_ptr<Tropism>, Tropism, copy, plant); }

    Vector2d getHeading(const Vector3d& pos, const Matrix3d& old,  double dx, const std::shared_ptr<Organ> organ = nullptr, int nodeIdx = -1) override
    { PYBIND11_OVERLOAD( Vector2d, Tropism, getHeading, pos, old, dx, organ, nodeIdx ); }

    Vector2d getUCHeading(const Vector3d& pos, const Matrix3d& old, double dx, const std::shared_ptr<Organ> organ, int nodeIdx) override
    { PYBIND11_OVERLOAD( Vector2d, Tropism, getUCHeading, pos, old, dx, organ, nodeIdx ); }

    double tropismObjective(const Vector3d& pos, const Matrix3d& old, double a, double b, double dx, const std::shared_ptr<Organ> organ = nullptr) override
    { PYBIND11_OVERLOAD( double, Tropism, tropismObjective, pos, old, a, b, dx, organ ); }

};

class PySoilLookUp : public SoilLookUp {
public:

    using SoilLookUp::SoilLookUp; /* Inherit the constructors */

    std::shared_ptr<SoilLookUp> copy() override
    { PYBIND11_OVERLOAD( std::shared_ptr<SoilLookUp>, SoilLookUp, copy); }

    double getValue(const Vector3d& pos, const std::shared_ptr<Organ> organ = nullptr) const override
    {  PYBIND11_OVERLOAD( double, SoilLookUp, getValue, pos, organ ); }

    std::string toString() const override
    { PYBIND11_OVERLOAD( std::string, SoilLookUp, toString); }

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
            .def(py::init<const std::vector<int>&>())
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
            .def(py::init<const std::vector<double>&>())
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
            .def(py::init<const std::vector<double>&>())
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

        py::class_<Quaternion>(m, "Quaternion", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<double, double, double, double>())
      .def(py::init<const Quaternion&>())
      .def(py::init<double, const Vector3d&>())
      .def("__add__", [](const Quaternion& q1, const Quaternion& q2) { return q1 + q2; }, py::is_operator())
      .def("__sub__", [](const Quaternion& q1, const Quaternion& q2) { return q1 - q2; }, py::is_operator())
      .def("__mul__", [](const Quaternion& q1, const Quaternion& q2) { return q1 * q2; }, py::is_operator())
      .def("__mul__", [](const Quaternion& q1, double d) { return q1 * d; }, py::is_operator())
      .def("__mul__", [](double d, const Quaternion& q1) { return q1 * d; }, py::is_operator())
      .def("__truediv__", [](const Quaternion& q1, double d) { return q1 / d; }, py::is_operator())
      .def("norm", &Quaternion::norm)
      .def("normalize", &Quaternion::normalize)
      .def("inverse", &Quaternion::inverse)
      .def("Forward", &Quaternion::Forward)
      .def("Up", &Quaternion::Up)
      .def("Right", &Quaternion::Right)
      .def("Rotate", (Vector3d (Quaternion::*)(const Vector3d&) const) &Quaternion::Rotate)
      .def_static("LookAt", &Quaternion::LookAt)
      .def_static("SphericalInterpolation", &Quaternion::SphericalInterpolation)
      .def_static("geodesicRotation", &Quaternion::geodesicRotation)
      .def_buffer([](Quaternion &q) -> py::buffer_info { /* enables numpy conversion with np.array(quaternion_instance, copy = False) */
        // this only really works if the variables are stored in a contiguous block of memory
        return py::buffer_info(
            &q.w,
            sizeof(double),
            py::format_descriptor<double>::format(),
            1,
            { 4 },
            { sizeof(double) }
        );
      })
    ;
    /*
     * sdf
     */
    py::class_<SignedDistanceFunction, std::shared_ptr<SignedDistanceFunction>>(m,"SignedDistanceFunction")
            .def(py::init<>())
            .def("getDist",&SignedDistanceFunction::getDist)
            .def("dist",&SignedDistanceFunction::dist)
            .def("writePVPScript", (std::string (SignedDistanceFunction::*)() const) &SignedDistanceFunction::writePVPScript) // overloads
            .def("getGradient",  &SignedDistanceFunction::getGradient, py::arg("p"), py::arg("eps") = 5.e-4) // defaults
            .def("__str__",&SignedDistanceFunction::toString);
    py::class_<SDF_PlantBox, SignedDistanceFunction, std::shared_ptr<SDF_PlantBox>>(m, "SDF_PlantBox")
            .def(py::init<double,double,double>());
    py::class_<SDF_Cuboid, SignedDistanceFunction, std::shared_ptr<SDF_Cuboid>>(m, "SDF_Cuboid")
            .def(py::init<Vector3d,Vector3d>());
    py::class_<SDF_PlantContainer, SignedDistanceFunction, std::shared_ptr<SDF_PlantContainer>>(m,"SDF_PlantContainer")
            .def(py::init<>())
            .def(py::init<double,double,double,double>());
    py::class_<SDF_RotateTranslate, SignedDistanceFunction, std::shared_ptr<SDF_RotateTranslate>>(m, "SDF_RotateTranslate")
            .def(py::init<std::shared_ptr<SignedDistanceFunction>,double,int,Vector3d&>())
            .def(py::init<std::shared_ptr<SignedDistanceFunction>,Vector3d&>());
    py::enum_<SDF_RotateTranslate::SDF_Axes>(m, "SDF_Axis")
            .value("xaxis", SDF_RotateTranslate::SDF_Axes::xaxis)
            .value("yaxis", SDF_RotateTranslate::SDF_Axes::yaxis)
            .value("zaxis", SDF_RotateTranslate::SDF_Axes::zaxis)
            .export_values();
    py::class_<SDF_Intersection, SignedDistanceFunction, std::shared_ptr<SDF_Intersection>>(m,"SDF_Intersection")
            .def(py::init<std::vector<std::shared_ptr<SignedDistanceFunction>>>())
            .def(py::init<std::shared_ptr<SignedDistanceFunction>,std::shared_ptr<SignedDistanceFunction>>());
    py::class_<SDF_Union, SDF_Intersection, std::shared_ptr<SDF_Union>>(m, "SDF_Union")
            .def(py::init<std::vector<std::shared_ptr<SignedDistanceFunction>>>())
            .def(py::init<std::shared_ptr<SignedDistanceFunction>,std::shared_ptr<SignedDistanceFunction>>());
    py::class_<SDF_Difference, SDF_Intersection, std::shared_ptr<SDF_Difference>>(m, "SDF_Difference")
            .def(py::init<std::vector<std::shared_ptr<SignedDistanceFunction>>>())
            .def(py::init<std::shared_ptr<SignedDistanceFunction>,std::shared_ptr<SignedDistanceFunction>>());
    py::class_<SDF_Complement, SignedDistanceFunction, std::shared_ptr<SDF_Complement>>(m, "SDF_Complement")
            .def(py::init<std::shared_ptr<SignedDistanceFunction>>());
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
            .def(py::init<int, double>())
            .def_readwrite("subType",&OrganSpecificParameter::subType)
			.def_readwrite("a",&OrganSpecificParameter::a)
            .def("__str__",&OrganSpecificParameter::toString);
    py::class_<OrganRandomParameter, std::shared_ptr<OrganRandomParameter>>(m,"OrganRandomParameter")
            .def(py::init<std::shared_ptr<Organism>>())
            .def("copy",&OrganRandomParameter::copy)
            .def("realize",&OrganRandomParameter::realize)
            .def("getParameter",&OrganRandomParameter::getParameter)
            .def("__str__",&OrganRandomParameter::toString, py::arg("verbose") = true) // default
            .def("writeXML",(void (OrganRandomParameter::*)(std::string name) const) &OrganRandomParameter::writeXML) // overloads
            .def("readXML", (void (OrganRandomParameter::*)(std::string name, bool verbose)) &OrganRandomParameter::readXML, py::arg("name"), py::arg("verbose") = false) // overloads
            .def("bindParameters",&OrganRandomParameter::bindParameters)
            .def("bindIntParameter", (void (OrganRandomParameter::*)(std::string, int*, std::string, double*)) &OrganRandomParameter::bindParameter, py::arg("name"), py::arg("i"), py::arg("descr") = "", py::arg("dev") = (double*) nullptr) // overloads, defaults
            .def("bindDoubleParameter", (void (OrganRandomParameter::*)(std::string, double*, std::string, double*))  &OrganRandomParameter::bindParameter, py::arg("name"), py::arg("i"), py::arg("descr") = "", py::arg("dev") = (double*) nullptr) // overloads, defaults
            .def_readwrite("name",&OrganRandomParameter::name)
			.def_readwrite("organType",&OrganRandomParameter::organType)
			.def_readwrite("subType",&OrganRandomParameter::subType)
            .def_readwrite("a", &OrganRandomParameter::a)
            .def_readwrite("a_s", &OrganRandomParameter::as) // as is a keyword in python
            .def_readwrite("dx", &OrganRandomParameter::dx)
            .def_readwrite("dxMin", &OrganRandomParameter::dxMin)
            .def_readwrite("successor", &OrganRandomParameter::successorST)//for backward compatibility
            .def_readwrite("successorOT", &OrganRandomParameter::successorOT)
            .def_readwrite("successorST", &OrganRandomParameter::successorST)
            .def_readwrite("successorWhere", &OrganRandomParameter::successorWhere)
            .def_readwrite("successorNo", &OrganRandomParameter::successorNo)
            .def_readwrite("successorP", &OrganRandomParameter::successorP)
            .def_readwrite("ldelay", &OrganRandomParameter::ldelay)
            .def_readwrite("ldelays", &OrganRandomParameter::ldelays)
            .def_readwrite("plant",&OrganRandomParameter::plant)
            .def_readwrite("f_gf", &OrganRandomParameter::f_gf)
            .def_readwrite("f_tf", &OrganRandomParameter::f_tf);
    /**
     * Organ.h
     */
    py::class_<Organ, std::shared_ptr<Organ>>(m, "Organ")
            .def(py::init<std::shared_ptr<Organism>, std::shared_ptr<Organ>, int, int, double, int>())
            .def(py::init<int, std::shared_ptr<const OrganSpecificParameter>, bool, bool, double, double, Vector3d, int, bool, int>())
            .def("copy",&Organ::copy)
            .def("organType",&Organ::organType)
            .def("simulate",&Organ::simulate, py::arg("dt"), py::arg("verbose") = bool(false) ) // default
			.def("getNumberOfLaterals", &Organ::getNumberOfLaterals)
			.def("setParent",&Organ::setParent)
            .def("getParent",&Organ::getParent)
            .def("setOrganism",&Organ::setOrganism)
            .def("getOrganism",&Organ::getOrganism)
            .def("addChild",&Organ::addChild)
            .def("getNumberOfChildren",&Organ::getNumberOfChildren)
            .def("getChild",&Organ::getChild)
			.def("calcCreationTime", &Organ::calcCreationTime)
            .def("getId",&Organ::getId)
            .def("getParam",&Organ::param)//backward compatibility
            .def("param",&Organ::param)
            .def("getOrganRandomParameter",&Organ::getOrganRandomParameter)
            .def("isAlive",&Organ::isAlive)
            .def("isActive",&Organ::isActive)
            .def("getAge",&Organ::getAge)
            .def("getLength", (double (Organ::*)(bool realized) const) &Organ::getLength, py::arg("realized") = true)
            .def("getLength", (double (Organ::*)(int i) const) &Organ::getLength)
			.def("getEpsilon",&Organ::getEpsilon)
			.def("getNumberOfNodes",&Organ::getNumberOfNodes)
			.def("getNumberOfSegments",&Organ::getNumberOfSegments)
            .def("getOrigin",&Organ::getOrigin)
			.def("getNode",&Organ::getNode)
            .def("getNodeId",&Organ::getNodeId)
            .def("getNodeIds",&Organ::getNodeIds)
            .def("getNodeCT",&Organ::getNodeCT)
            .def("addNode",(void (Organ::*)(Vector3d n, double t)) &Organ::addNode,  py::arg("n"), py::arg("t")) // overloads
            .def("addNode",(void (Organ::*)(Vector3d n, int id, double t)) &Organ::addNode,  py::arg("n"),  py::arg("id"),py::arg("t")) // overloads
            .def("getSegments",&Organ::getSegments)
            .def("dx",&Organ::dx)
            .def("dxMin",&Organ::dxMin)
            .def("hasMoved",&Organ::hasMoved)
            .def("getOldNumberOfNodes",&Organ::getOldNumberOfNodes)
            .def("getNodes",&Organ::getNodes)
            .def("getOrgans", (std::vector<std::shared_ptr<Organ>> (Organ::*)(int otype, bool all)) &Organ::getOrgans, py::arg("ot")=-1, py::arg("all")=false) //overloads, default
            .def("getOrgans", (void (Organ::*)(int otype, std::vector<std::shared_ptr<Organ>>& v, bool all)) &Organ::getOrgans)
            .def("getParameter",&Organ::getParameter)
            .def("__str__",&Organ::toString)
            .def("orgVolume",&Organ::orgVolume, py::arg("length_")=-1, py::arg("realized")=false)
			.def("orgVolume2Length",&Organ::orgVolume2Length)
            .def("getiHeading0", &Organ::getiHeading0)
            .def_readwrite("parentNI", &Organ::parentNI);

    /*
     * Organism.h
     */
    py::class_<Organism, std::shared_ptr<Organism>>(m, "Organism")
            .def(py::init<unsigned int>(),  py::arg("seednum") = 0)
            .def("copy", &Organism::copy)
            .def("organTypeNumber", &Organism::organTypeNumber)
            .def("organTypeName", &Organism::organTypeName)
            .def("getOrganRandomParameter", (std::shared_ptr<OrganRandomParameter> (Organism::*)(int, int) const)  &Organism::getOrganRandomParameter) //overloads
            .def("getOrganRandomParameter", (std::vector<std::shared_ptr<OrganRandomParameter>> (Organism::*)(int) const) &Organism::getOrganRandomParameter) //overloads
            .def("setOrganRandomParameter", &Organism::setOrganRandomParameter)
            .def("getSeed", &Organism::getSeed)

            .def("addOrgan", &Organism::addOrgan)
            .def("initialize", &Organism::initialize, py::arg("verbose") = true)
            .def("simulate", &Organism::simulate, py::arg("dt"), py::arg("verbose") = false) //default
            .def("getSimTime", &Organism::getSimTime)

            .def("getOrgans", &Organism::getOrgans, py::arg("ot") = -1, py::arg("allOrgs")=false) // default
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
            .def("readParameters", &Organism::readParameters, py::arg("name"), py::arg("basetag") = "plant", py::arg("fromFile") = true, py::arg("verbose") = false)  // default
            .def("writeParameters", &Organism::writeParameters, py::arg("name"), py::arg("basetag") = "plant", py::arg("comments") = true)  // default
            .def("writeRSML", &Organism::writeRSML)
            .def("getRSMLSkip", &Organism::getRSMLSkip)
            .def("setRSMLSkip", &Organism::setRSMLSkip)
            .def("getRSMLProperties", &Organism::getRSMLProperties) //todo policy

            .def("getOrganIndex", &Organism::getOrganIndex)
            .def("getNodeIndex", &Organism::getNodeIndex)

            .def("setMinDx", &Organism::setMinDx)
            .def("setSeed", &Organism::setSeed)
            .def("getSeedVal", &Organism::getSeedVal)
            .def_readwrite("plantId", &Organism::plantId)

            .def("rand", &Organism::rand)
            .def("randn", &Organism::randn)
            //        .def_readwrite("seed_nC_", &Organism::seed_nC_)
            //        .def_readwrite("seed_nZ_", &Organism::seed_nZ_)
            .def("__str__",&Organism::toString);


    py::enum_<Organism::OrganTypes>(m, "OrganTypes")
            .value("organ", Organism::OrganTypes::ot_organ)
            .value("seed", Organism::OrganTypes::ot_seed)
            .value("root", Organism::OrganTypes::ot_root)
            .value("stem", Organism::OrganTypes::ot_stem)
            .value("leaf", Organism::OrganTypes::ot_leaf)
            .export_values();
    /*
     * soil.h
     */
    py::class_<SoilLookUp, PySoilLookUp, std::shared_ptr<SoilLookUp>>(m, "SoilLookUp")
            .def(py::init<>())
            .def("getValue",&SoilLookUp::getValue, py::arg("pos"), py::arg("organ") = (std::shared_ptr<Organ>) nullptr )
            .def("__str__",&SoilLookUp::toString);
    py::class_<SoilLookUpSDF, SoilLookUp, std::shared_ptr<SoilLookUpSDF>>(m,"SoilLookUpSDF")
            .def(py::init<>())
            .def(py::init<std::shared_ptr<SignedDistanceFunction>, double, double, double>())
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
            .def("getScale", &ProportionalElongation::getScale)
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
    py::class_<RectilinearGrid3D, SoilLookUp, std::shared_ptr<RectilinearGrid3D>>(m, "RectilinearGrid3D")
            .def(py::init<Grid1D*,Grid1D*,Grid1D*>())
            .def("map",&RectilinearGrid3D::map)
			.def("getData",&RectilinearGrid3D::getData)
			.def("setData",&RectilinearGrid3D::setData)
			.def("getGridPoint",&RectilinearGrid3D::getGridPoint)
            .def_readwrite("xgrid", &RectilinearGrid3D::xgrid)
            .def_readwrite("ygrid", &RectilinearGrid3D::ygrid)
            .def_readwrite("zgrid", &RectilinearGrid3D::zgrid)
            .def_readwrite("nx", &RectilinearGrid3D::nx)
            .def_readwrite("ny", &RectilinearGrid3D::ny)
            .def_readwrite("nz", &RectilinearGrid3D::nz)
            .def_readwrite("data", &RectilinearGrid3D::data);
    py::class_<EquidistantGrid3D, RectilinearGrid3D, std::shared_ptr<EquidistantGrid3D>>(m, "EquidistantGrid3D")
		.def(py::init<>())
		.def(py::init<double, double, double, int, int, int>())
		.def(py::init<double, double, int, double, double, int, double, double, int>());
    /**
     * tropism.h
     */
    py::class_<Tropism, PyTropism, std::shared_ptr<Tropism>>(m, "Tropism")
            .def(py::init<std::shared_ptr<Organism>>())
            .def(py::init<std::shared_ptr<Organism>, double, double>())
            .def("copy",&Tropism::copy) // todo policy
            .def("setGeometry",&Tropism::setGeometry)
            .def("setTropismParameter",&Tropism::setTropismParameter)
            .def("setSigma",&Tropism::setSigma)
            .def("getHeading",&Tropism::getHeading)
            .def("getUCHeading",&Tropism::getUCHeading)
            .def("tropismObjective",&Tropism::tropismObjective)
            .def("getPosition",&Tropism::getPosition)
            .def_readwrite("alphaN", &Tropism::alphaN)
            .def_readwrite("betaN", &Tropism::betaN)
            .def("isExpired",&Tropism::isExpired)
            .def("getPlant",&Tropism::getPlant);

    py::class_<Gravitropism, Tropism, std::shared_ptr<Gravitropism>>(m, "Gravitropism")
            .def(py::init<std::shared_ptr<Organism>, double, double>());
    py::class_<Plagiotropism, Tropism, std::shared_ptr<Plagiotropism>>(m, "Plagiotropism")
            .def(py::init<std::shared_ptr<Organism>,double, double>());
    py::class_<Exotropism, Tropism, std::shared_ptr<Exotropism>>(m, "Exotropism")
            .def(py::init<std::shared_ptr<Organism>,double, double>());
    py::class_<Hydrotropism, Tropism, std::shared_ptr<Hydrotropism>>(m, "Hydrotropism")
            .def(py::init<std::shared_ptr<Organism>,double, double, std::shared_ptr<SoilLookUp>>());
    //    py::class_<CombinedTropism, Tropism>(m, "CombinedTropism") // Todo constructors needs some extra work (?)
    //        .def(py::init<>());
    // todo antigravi, twist ...
    /*
     * analysis.h
     */
    py::class_<SegmentAnalyser, std::shared_ptr<SegmentAnalyser>>(m, "SegmentAnalyser")
           .def(py::init<>())
           .def(py::init<std::vector<Vector3d>, std::vector<Vector2i>, std::vector<double>, std::vector<double>>())
           .def(py::init<Organism&>())
           .def(py::init<SegmentAnalyser&>())
           .def(py::init<MappedSegments&>())
           .def("addSegments",(void (SegmentAnalyser::*)(const Organism&)) &SegmentAnalyser::addSegments) //overloads
           .def("addSegments",(void (SegmentAnalyser::*)(const SegmentAnalyser&)) &SegmentAnalyser::addSegments) //overloads
           .def("addSegment", &SegmentAnalyser::addSegment, py::arg("seg"), py::arg("ct"), py::arg("radius"), py::arg("insert") = false)
           .def("addAge", &SegmentAnalyser::addAge)
           .def("addConductivities", &SegmentAnalyser::addConductivities, py::arg("rs"), py::arg("simTime"), py::arg("kr_max") = 1.e6, py::arg("kx_max") = 1.e6)
           .def("addHydraulicConductivities", &SegmentAnalyser::addHydraulicConductivities, py::arg("rs"), py::arg("simTime"), py::arg("kr_max") = 1.e6, py::arg("kx_max") = 1.e6)
           .def("addFluxes", &SegmentAnalyser::addFluxes)
           .def("addCellIds", &SegmentAnalyser::addCellIds)
           .def("crop", &SegmentAnalyser::crop)
           .def("cropDomain", &SegmentAnalyser::cropDomain)
           .def("filter", (void (SegmentAnalyser::*)(std::string, double, double)) &SegmentAnalyser::filter) //overloads
           .def("filter", (void (SegmentAnalyser::*)(std::string, double)) &SegmentAnalyser::filter) //overloads
           .def("pack", &SegmentAnalyser::pack)
           .def("getMinBounds", &SegmentAnalyser::getMinBounds)
           .def("getMaxBounds", &SegmentAnalyser::getMaxBounds)
           .def("getParameter", &SegmentAnalyser::getParameter, py::arg("name"), py::arg("def") = std::numeric_limits<double>::quiet_NaN())
           .def("getSegmentLength", &SegmentAnalyser::getSegmentLength)
           .def("getSummed", (double (SegmentAnalyser::*)(std::string) const) &SegmentAnalyser::getSummed) //overloads
           .def("getSummed", (double (SegmentAnalyser::*)(std::string, std::shared_ptr<SignedDistanceFunction>) const) &SegmentAnalyser::getSummed) //overloads
           .def("distribution", (std::vector<double> (SegmentAnalyser::*)(std::string, double, double, int, bool) const) &SegmentAnalyser::distribution) //overloads
           .def("distribution", (std::vector<SegmentAnalyser> (SegmentAnalyser::*)(double, double, int) const) &SegmentAnalyser::distribution) //overloads
           .def("distribution2", (std::vector<std::vector<double>> (SegmentAnalyser::*)(std::string, double, double, double, double, int, int, bool) const) &SegmentAnalyser::distribution2) //overloads
           .def("distribution2", (std::vector<std::vector<SegmentAnalyser>> (SegmentAnalyser::*)(double, double, double, double, int, int) const) &SegmentAnalyser::distribution2) //overloads
           .def("mapPeriodic", &SegmentAnalyser::mapPeriodic)
           .def("map2D", &SegmentAnalyser::map2D)
           .def("getOrgans", &SegmentAnalyser::getOrgans, py::arg("ot") = -1)
           .def("getNumberOfOrgans", &SegmentAnalyser::getNumberOfOrgans)
           .def("cut", (SegmentAnalyser (SegmentAnalyser::*)(const SignedDistanceFunction&) const) &SegmentAnalyser::cut)
           .def("addData", &SegmentAnalyser::addData)
           .def("write", &SegmentAnalyser::write, py::arg("name"), py::arg("types") = std::vector<std::string>({"radius", "subType", "creationTime", "organType"}))
           .def_readwrite("nodes", &SegmentAnalyser::nodes)
           .def_readwrite("segments", &SegmentAnalyser::segments)
           .def_readwrite("segO", &SegmentAnalyser::segO)
           .def_readwrite("data", &SegmentAnalyser::data);
    /*
     * rootparameter.h
     */
    py::class_<RootRandomParameter, OrganRandomParameter, std::shared_ptr<RootRandomParameter>>(m, "RootRandomParameter")
            .def(py::init<std::shared_ptr<Organism>>())
            .def("getLateralType",&RootRandomParameter::getLateralType)
            .def("nob",&RootRandomParameter::nob)
            .def("nobs",&RootRandomParameter::nobs)
            .def_readwrite("lb", &RootRandomParameter::lb)
            .def_readwrite("lbs", &RootRandomParameter::lbs)
            .def_readwrite("la", &RootRandomParameter::la)
            .def_readwrite("las", &RootRandomParameter::las)
            .def_readwrite("ln", &RootRandomParameter::ln)
            .def_readwrite("lns", &RootRandomParameter::lns)
            .def_readwrite("lmax", &RootRandomParameter::lmax)
            .def_readwrite("lmaxs", &RootRandomParameter::lmaxs)
            .def_readwrite("r", &RootRandomParameter::r)
            .def_readwrite("rs", &RootRandomParameter::rs)
            .def_readwrite("tropismT", &RootRandomParameter::tropismT)
            .def_readwrite("tropismN", &RootRandomParameter::tropismN)
            .def_readwrite("tropismS", &RootRandomParameter::tropismS)
            .def_readwrite("theta", &RootRandomParameter::theta)
            .def_readwrite("thetas", &RootRandomParameter::thetas)
            .def_readwrite("rlt", &RootRandomParameter::rlt)
            .def_readwrite("rlts", &RootRandomParameter::rlts)
            .def_readwrite("gf", &RootRandomParameter::gf)
            .def_readwrite("lnk", &RootRandomParameter::lnk)
            .def_readwrite("f_se", &RootRandomParameter::f_se)
            .def_readwrite("f_sa", &RootRandomParameter::f_sa)
            .def_readwrite("f_sbp", &RootRandomParameter::f_sbp)
			.def_readwrite("hairsElongation", &RootRandomParameter::hairsElongation)
			.def_readwrite("hairsZone", &RootRandomParameter::hairsZone)
			.def_readwrite("hairsLength", &RootRandomParameter::hairsLength);

    py::class_<RootSpecificParameter, OrganSpecificParameter, std::shared_ptr<RootSpecificParameter>>(m, "RootSpecificParameter")
            .def(py::init<>())
            .def(py::init<int , double, double, const std::vector<double>&, double, double, double, double, bool>()) // <---------------------------------------------------
            .def_readwrite("lb", &RootSpecificParameter::lb)
            .def_readwrite("la", &RootSpecificParameter::la)
            .def_readwrite("ln", &RootSpecificParameter::ln)
            .def_readwrite("r", &RootSpecificParameter::r)
            .def_readwrite("a", &RootSpecificParameter::a)
            .def_readwrite("theta", &RootSpecificParameter::theta)
            .def_readwrite("rlt", &RootSpecificParameter::rlt)
            .def_readwrite("laterals", &RootSpecificParameter::laterals)
            .def("getK",&RootSpecificParameter::getK)
            .def("nob", &RootSpecificParameter::nob);
    /*
     * seedparameter.h
     */
    py::class_<SeedRandomParameter, OrganRandomParameter, std::shared_ptr<SeedRandomParameter>>(m, "SeedRandomParameter")
            .def(py::init<std::shared_ptr<Organism>>())
            .def_readwrite("seedPos", &SeedRandomParameter::seedPos)
            .def_readwrite("seedPoss", &SeedRandomParameter::seedPoss)
            .def_readwrite("delayDefinition", &SeedRandomParameter::delayDefinition)
            .def_readwrite("delayDefinitionShoot", &SeedRandomParameter::delayDefinitionShoot)
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
            .def_readwrite("delayTil", &SeedRandomParameter::delayTil)
            .def_readwrite("firstTil", &SeedRandomParameter::firstTil)
            .def_readwrite("maxTil", &SeedRandomParameter::maxTil)
            .def_readwrite("maxTils", &SeedRandomParameter::maxTils)
            .def_readwrite("simtime", &SeedRandomParameter::simtime)
            .def_readwrite("simtimes", &SeedRandomParameter::simtimes);
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
            .def_readwrite("maxTil", &SeedSpecificParameter::maxTil)
            .def_readwrite("simtime", &SeedSpecificParameter::simtime);
    /*
     * leafparameter.h
     */
    py::class_<LeafRandomParameter, OrganRandomParameter, std::shared_ptr<LeafRandomParameter>>(m, "LeafRandomParameter")
      .def(py::init<std::shared_ptr<Organism>>())
      .def("createGeometry", &LeafRandomParameter::createGeometry)

      .def("createLeafGeometry",&LeafRandomParameter::createLeafGeometry) // don't use directly: set parametrisationType, leafGeometryPhi, leafGeometryX, geometryN and use createGeometry()
      .def("createLeafRadialGeometry",&LeafRandomParameter::createLeafRadialGeometry) // don't use directly: set parametrisationType, leafGeometryX, geometryN and use createGeometry()

      .def("getLateralType",&LeafRandomParameter::getLateralType)
      .def("nob",&LeafRandomParameter::nob)
      .def("nobs",&LeafRandomParameter::nobs)
      .def("leafLength",&LeafRandomParameter::leafLength)
      .def("leafMid",&LeafRandomParameter::leafMid)
      .def_readwrite("lb", &LeafRandomParameter::lb)
      .def_readwrite("lbs", &LeafRandomParameter::lbs)
      .def_readwrite("la", &LeafRandomParameter::la)
      .def_readwrite("las", &LeafRandomParameter::las)
      .def_readwrite("ln", &LeafRandomParameter::ln)
      .def_readwrite("lns", &LeafRandomParameter::lns)
      .def_readwrite("lnf", &LeafRandomParameter::lnf)
      .def_readwrite("lmax", &LeafRandomParameter::lmax)
      .def_readwrite("lmaxs", &LeafRandomParameter::lmaxs)
      .def_readwrite("areaMax", &LeafRandomParameter::areaMax)
      .def_readwrite("areaMaxs", &LeafRandomParameter::areaMaxs)
      .def_readwrite("r", &LeafRandomParameter::r)
      .def_readwrite("rs", &LeafRandomParameter::rs)
      .def_readwrite("RotBeta", &LeafRandomParameter::rotBeta)
      .def_readwrite("BetaDev", &LeafRandomParameter::betaDev)
      .def_readwrite("InitBeta", &LeafRandomParameter::initBeta)
      .def_readwrite("rotBeta", &LeafRandomParameter::rotBeta)
      .def_readwrite("betaDev", &LeafRandomParameter::betaDev)
      .def_readwrite("initBeta", &LeafRandomParameter::initBeta)
      .def_readwrite("tropismT", &LeafRandomParameter::tropismT)
      .def_readwrite("tropismN", &LeafRandomParameter::tropismN)
      .def_readwrite("tropismS", &LeafRandomParameter::tropismS)
      .def_readwrite("theta", &LeafRandomParameter::theta)
      .def_readwrite("thetas", &LeafRandomParameter::thetas)
      .def_readwrite("rlt", &LeafRandomParameter::rlt)
      .def_readwrite("rlts", &LeafRandomParameter::rlts)
      .def_readwrite("gf", &LeafRandomParameter::gf)
      .def_readwrite("f_se", &LeafRandomParameter::f_se)
      .def_readwrite("f_sa", &LeafRandomParameter::f_sa)
      .def_readwrite("f_sbp", &LeafRandomParameter::f_sbp)
      .def_readwrite("tropismAge", &LeafRandomParameter::tropismAge)
      .def_readwrite("tropismAges", &LeafRandomParameter::tropismAges)
      .def_readwrite("Width_blade", &LeafRandomParameter::Width_blade)
      .def_readwrite("geometryN", &LeafRandomParameter::geometryN)
      .def_readwrite("parametrisationType", &LeafRandomParameter::parametrisationType)
      .def_readwrite("leafGeometry", &LeafRandomParameter::leafGeometry)
      .def_readwrite("leafGeometryPhi", &LeafRandomParameter::leafGeometryPhi)
      .def_readwrite("leafGeometryX", &LeafRandomParameter::leafGeometryX)
      .def_readwrite("shapeType", &LeafRandomParameter::shapeType);

    py::class_<LeafSpecificParameter, OrganSpecificParameter, std::shared_ptr<LeafSpecificParameter>>(m, "LeafSpecificParameter")
            .def(py::init<>())
            .def(py::init<int , double, double, const std::vector<double>&, double, double, double, double, double, bool, double, double>())
            .def_readwrite("lb", &LeafSpecificParameter::lb)
            .def_readwrite("la", &LeafSpecificParameter::la)
            .def_readwrite("ln", &LeafSpecificParameter::ln)
            .def_readwrite("r", &LeafSpecificParameter::r)
            .def_readwrite("a", &LeafSpecificParameter::a)
            .def_readwrite("theta", &LeafSpecificParameter::theta)
            .def_readwrite("rlt", &LeafSpecificParameter::rlt)
            .def_readwrite("leafArea", &LeafSpecificParameter::areaMax)
            .def_readwrite("laterals", &LeafSpecificParameter::laterals)
            .def_readwrite("laterals", &LeafSpecificParameter::width_blade)
            .def_readwrite("laterals", &LeafSpecificParameter::width_petiole)
			.def("getK",&LeafSpecificParameter::getK)
            .def("nob",&LeafSpecificParameter::nob);
    /*
     * stemparameter.h
     */
    py::class_<StemRandomParameter, OrganRandomParameter, std::shared_ptr<StemRandomParameter>>(m, "StemRandomParameter")
            .def(py::init<std::shared_ptr<Organism>>())
            .def("getLateralType",&StemRandomParameter::getLateralType)
            .def("nob",&StemRandomParameter::nob)
            .def("nobs",&StemRandomParameter::nobs)
            .def_readwrite("lb", &StemRandomParameter::lb)
            .def_readwrite("lbs", &StemRandomParameter::lbs)
            .def_readwrite("la", &StemRandomParameter::la)
            .def_readwrite("las", &StemRandomParameter::las)
            .def_readwrite("ln", &StemRandomParameter::ln)
            .def_readwrite("lns", &StemRandomParameter::lns)
            .def_readwrite("lnf", &StemRandomParameter::lnf)
            .def_readwrite("lmax", &StemRandomParameter::lmax)
            .def_readwrite("lmaxs", &StemRandomParameter::lmaxs)
            .def_readwrite("r", &StemRandomParameter::r)
            .def_readwrite("rs", &StemRandomParameter::rs)
            .def_readwrite("RotBeta", &StemRandomParameter::rotBeta)
            .def_readwrite("BetaDev", &StemRandomParameter::betaDev)
            .def_readwrite("InitBeta", &StemRandomParameter::initBeta)
            .def_readwrite("rotBeta", &StemRandomParameter::rotBeta)
            .def_readwrite("betaDev", &StemRandomParameter::betaDev)
            .def_readwrite("initBeta", &StemRandomParameter::initBeta)
            .def_readwrite("tropismT", &StemRandomParameter::tropismT)
            .def_readwrite("tropismN", &StemRandomParameter::tropismN)
            .def_readwrite("tropismS", &StemRandomParameter::tropismS)
            .def_readwrite("theta", &StemRandomParameter::theta)
            .def_readwrite("thetas", &StemRandomParameter::thetas)
            .def_readwrite("rlt", &StemRandomParameter::rlt)
            .def_readwrite("rlts", &StemRandomParameter::rlts)
            .def_readwrite("gf", &StemRandomParameter::gf)
            .def_readwrite("f_se", &StemRandomParameter::f_se)
            .def_readwrite("f_sa", &StemRandomParameter::f_sa)
            .def_readwrite("f_sbp", &StemRandomParameter::f_sbp)
            .def_readwrite("nodalGrowth", &StemRandomParameter::nodalGrowth)
            .def_readwrite("delayNGStart", &StemRandomParameter::delayNGStart)
            .def_readwrite("delayNGStarts", &StemRandomParameter::delayNGStarts)
            .def_readwrite("delayNGEnd", &StemRandomParameter::delayNGEnd)
            .def_readwrite("delayNGEnds", &StemRandomParameter::delayNGEnds)
            .def_readwrite("tropismAge", &StemRandomParameter::tropismAge)
            .def_readwrite("tropismAges", &StemRandomParameter::tropismAges);
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
            .def("getK",&StemSpecificParameter::getK)
            .def("nob", &StemSpecificParameter::nob)
            .def_readwrite("delayNGStart", &StemSpecificParameter::delayNGStart)
            .def_readwrite("delayNGEnd", &StemSpecificParameter::delayNGEnd);
    /**
     * Root.h
     */
    py::class_<Root, Organ, std::shared_ptr<Root>>(m, "Root")
            .def(py::init<std::shared_ptr<Organism>, int, double, std::shared_ptr<Organ>, int>())
            .def(py::init<int, std::shared_ptr<OrganSpecificParameter>, bool, bool, double, double, Vector3d, int, bool, int>())
			.def("calcLength", &Root::calcLength)
            .def("calcAge", &Root::calcAge)
            .def("getRootRandomParameter", &Root::getRootRandomParameter)
            .def("param", &Root::param);
    py::class_<StaticRoot, Root, std::shared_ptr<StaticRoot>>(m, "StaticRoot")
            .def(py::init<int, std::shared_ptr<OrganSpecificParameter>, double, int>())
            .def("initializeLaterals", &StaticRoot::initializeLaterals)
            .def("addLateral", &StaticRoot::addLateral);

    /**
     * Seed.h
     */
    py::class_<Seed, Organ, std::shared_ptr<Seed>>(m, "Seed")
            .def(py::init<int, std::shared_ptr<const OrganSpecificParameter>, bool, bool, double, double, bool, int>())
            .def(py::init<std::shared_ptr<Organism>>())
            .def("initialize", &Seed::initialize, py::arg("verbose") = true)
            .def("param", &Seed::param)
            .def("getNumberOfRootCrowns", &Seed::getNumberOfRootCrowns)
            .def("baseOrgans", &Seed::baseOrgans)
            .def("copyBaseOrgans", &Seed::copyBaseOrgans)
            .def("createRoot", &Seed::createRoot)
            .def("createStem", &Seed::createStem)
            .def_readwrite("tapType", &Seed::tapType)
            .def_readwrite("basalType", &Seed::basalType)
            .def_readwrite("shootborneType", &Seed::shootborneType)
            .def_readwrite("mainStemType", &Seed::mainStemType)
            .def_readwrite("tillerType", &Seed::tillerType);
    /**
     * Leaf.h
     */
    py::class_<Leaf, Organ, std::shared_ptr<Leaf>>(m, "Leaf")
            .def(py::init<std::shared_ptr<Organism>, int, double, std::shared_ptr<Organ>, int>())
            .def(py::init<int, std::shared_ptr<OrganSpecificParameter>, bool, bool, double, double, Vector3d, int, bool, int>())
			.def("getLeafVis", &Leaf::getLeafVis)
			.def("getLeafVisX", &Leaf::getLeafVisX)
            .def("calcLength", &Leaf::calcLength)
            .def("calcAge", &Leaf::calcAge)
            .def("getLeafRandomParameter", &Leaf::getLeafRandomParameter)
            .def("orgVolume",&Leaf::orgVolume, py::arg("length_")=-1, py::arg("realized")=false)
			.def("orgVolume2Length",&Leaf::orgVolume2Length)
            .def("param", &Leaf::param)
            .def("leafArea", &Leaf::leafArea, py::arg("realized")=false, py::arg("withPetiole")=false)
            .def("leafCenter", &Leaf::leafCenter, py::arg("realized")=false)
            .def("leafLength", &Leaf::leafLength, py::arg("realized")=false);
    /**
     * Stem.h
     */
    py::class_<Stem, Organ, std::shared_ptr<Stem>>(m, "Stem")
           .def(py::init<std::shared_ptr<Organism>, int, double, std::shared_ptr<Organ>, int>())
           .def(py::init<int, std::shared_ptr<OrganSpecificParameter>, bool, bool, double, double, Vector3d, int, bool, int>())
           .def("calcLength", &Stem::calcLength)
           .def("calcAge", &Stem::calcAge)
           .def("getStemRandomParameter", &Stem::getStemRandomParameter)
           .def("param", &Stem::param);
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
            .def("initialize", (void (RootSystem::*)(bool)) &RootSystem::initialize, py::arg("verbose") = true)
            .def("initializeLB", (void (RootSystem::*)(int, int, bool)) &RootSystem::initializeLB, py::arg("basal"), py::arg("shootborne"), py::arg("verbose") = true)
            .def("initializeDB", (void (RootSystem::*)(int, int, bool)) &RootSystem::initializeDB, py::arg("basal"), py::arg("shootborne"), py::arg("verbose") = true)
			.def("setTropism", &RootSystem::setTropism)
            .def("simulate",(void (RootSystem::*)(double,bool)) &RootSystem::simulate, py::arg("dt"), py::arg("verbose") = false)
            .def("simulate",(void (RootSystem::*)()) &RootSystem::simulate)
            .def("simulate",(void (RootSystem::*)(double, double, std::shared_ptr<ProportionalElongation>, bool)) &RootSystem::simulate)
            .def("getRoots", &RootSystem::getRoots)
            .def("initCallbacks", &RootSystem::initCallbacks)
            .def("createTropismFunction", &RootSystem::createTropismFunction)
            .def("createGrowthFunction", &RootSystem::createGrowthFunction)
            .def("getNumberOfRoots", &RootSystem::getNumberOfRoots, py::arg("all") = false)
            .def("getBaseRoots", &RootSystem::getBaseRoots)
            .def("getShootSegments", &RootSystem::getShootSegments)
            .def("getRootTips", &RootSystem::getRootTips)
            .def("getRootBases", &RootSystem::getRootBases)
            .def("write", &RootSystem::write);
    /*
     * MappedOrganism.h
     */
    py::class_<MappedSegments, std::shared_ptr<MappedSegments>>(m, "MappedSegments")
        .def(py::init<>())
        .def(py::init<std::vector<Vector3d>, std::vector<double>, std::vector<Vector2i>, std::vector<double>, std::vector<int>,  std::vector<int>>())
        .def(py::init<std::vector<Vector3d>, std::vector<double>, std::vector<Vector2i>, std::vector<double>, std::vector<int>>())
        .def(py::init<std::vector<Vector3d>, std::vector<Vector2i>, std::vector<double>>())
        .def("setRadius", &MappedSegments::setRadius)
        .def("setTypes", &MappedSegments::setSubTypes) //kept for backward compatibility
        .def("setSubTypes", &MappedSegments::setSubTypes)
        .def("setSoilGrid", (void (MappedSegments::*)(const std::function<int(double,double,double)>&, bool)) &MappedSegments::setSoilGrid, py::arg("s"), py::arg("noChanges") = false)
        .def("setSoilGrid", (void (MappedSegments::*)(const std::function<int(double,double,double)>&, Vector3d, Vector3d, Vector3d, bool, bool)) &MappedSegments::setSoilGrid,
        		py::arg("s"), py::arg("min"), py::arg("max"), py::arg("res"), py::arg("cut") = true, py::arg("noChanges") = false)
        .def("setRectangularGrid", &MappedSegments::setRectangularGrid, py::arg("min"), py::arg("max"), py::arg("res"),
				py::arg("cut") = true, py::arg("noChanges") = false)
        .def("mapSegments",  &MappedSegments::mapSegments)
        .def("cutSegments", &MappedSegments::cutSegments)
        .def_readwrite("soil_index", &MappedSegments::soil_index)
        .def("sort",&MappedSegments::sort)
		.def("segLength",&MappedSegments::segLength)
		.def("getHs",&MappedSegments::getHs)
        .def("getSegmentZ",&MappedSegments::getSegmentZ)
        .def("matric2total",&MappedSegments::matric2total)
        .def("total2matric",&MappedSegments::total2matric)
		.def("getNumberOfMappedSegments",&MappedSegments::getNumberOfMappedSegments)
		.def("getNumberOfMappedNodes",&MappedSegments::getNumberOfMappedNodes)
        .def("getSegmentMapper",&MappedSegments::getSegmentMapper)
        .def("getEffectiveRadius",&MappedSegments::getEffectiveRadius)
        .def("getEffectiveRadii",&MappedSegments::getEffectiveRadii)
        .def("getMinBounds",&MappedSegments::getMinBounds)
		.def("calcExchangeZoneCoefs",&MappedSegments::calcExchangeZoneCoefs)
        .def_readwrite("exchangeZoneCoefs", &MappedPlant::exchangeZoneCoefs)
        .def_readwrite("distanceTip", &MappedPlant::distanceTip)
        .def_readwrite("nodes", &MappedSegments::nodes)
        .def_readwrite("nodeCTs", &MappedSegments::nodeCTs)
        .def_readwrite("segments", &MappedSegments::segments)
        .def_readwrite("radii", &MappedSegments::radii)
        .def_readwrite("organTypes", &MappedSegments::organTypes)
        .def_readwrite("Types", &MappedSegments::subTypes) //kept for backward compatibility
        .def_readwrite("subTypes", &MappedSegments::subTypes)
        .def_readwrite("seg2cell", &MappedSegments::seg2cell)
        .def_readwrite("cell2seg", &MappedSegments::cell2seg)
        .def_readwrite("minBound", &MappedSegments::minBound)
        .def_readwrite("maxBound", &MappedSegments::maxBound)
        .def_readwrite("resolution", &MappedSegments::resolution)
//		.def_readwrite("organParam", &MappedSegments::plantParam)
	    .def_readwrite("segO", &MappedSegments::segO)
		.def("sumSegFluxes",&MappedSegments::sumSegFluxes)
		.def("splitSoilFluxes",&MappedSegments::splitSoilFluxes, py::arg("soilFluxes"), py::arg("type") = 0);

    /*
     * Plant.h
     */
    py::class_<Plant, Organism, std::shared_ptr<Plant>>(m, "Plant")
            .def(py::init<unsigned int>(),  py::arg("seednum")=0)
            .def("initialize", &Plant::initialize, py::arg("verbose") = true)
			.def("initializeLB", &Plant::initialize, py::arg("verbose") = true)
            .def("initializeDB", &Plant::initialize, py::arg("verbose") = true)
			.def("setGeometry", &Plant::setGeometry)
            .def("setSoil", &Plant::setSoil)
            .def("reset", &Plant::reset)
            .def("openXML", &Plant::openXML)
            .def("setTropism", &Plant::setTropism)
            .def("simulate",(void (Plant::*)(double,bool)) &Plant::simulate, py::arg("dt"), py::arg("verbose") = false)
            .def("simulate",(void (Plant::*)()) &Plant::simulate)
            .def("simulate",(void (Plant::*)(double, double, std::shared_ptr<ProportionalElongation>, bool)) &Plant::simulate)
            .def("simulateLimited", &Plant::simulateLimited)
            .def("initCallbacks", &Plant::initCallbacks)
            .def("createTropismFunction", &Plant::createTropismFunction)
            .def("createGrowthFunction", &Plant::createGrowthFunction)
            .def("write", &Plant::write)
            .def("abs2rel", &Plant::abs2rel)
            .def("rel2abs", &Plant::rel2abs);


	py::class_<MappedPlant, Plant, MappedSegments,  std::shared_ptr<MappedPlant>>(m, "MappedPlant")
			.def(py::init<unsigned int>(),  py::arg("seednum")=0)
			.def("mappedSegments", &MappedPlant::mappedSegments)
			.def("printNodes",  &MappedPlant::printNodes)
			.def("plant", &MappedPlant::plant)
			.def("getSegmentIds",&MappedPlant::getSegmentIds)
			.def("disableExtraNode",&MappedPlant::disableExtraNode)
            .def("enableExtraNode",&MappedPlant::enableExtraNode)
			.def_readwrite("leafBladeSurface",  &MappedPlant::leafBladeSurface)
			.def_readwrite("bladeLength",  &MappedPlant::bladeLength)
			.def("getNodeIds",&MappedPlant::getNodeIds);

	/**
	 * Perirhizal.h
	 */
    py::class_<Perirhizal, std::shared_ptr<Perirhizal>> (m, "Perirhizal")
            .def(py::init<>())
            .def(py::init<std::shared_ptr<MappedSegments>>())
            .def("segOuterRadii",&Perirhizal::segOuterRadii,  py::arg("type"), py::arg("vols") = std::vector<double>(0))
            .def("adapt_values",&Perirhizal::adapt_values,  py::arg("val_new_"), py::arg("minVal_"), py::arg("maxVal_")=-1., py::arg("volumes_"), py::arg("divideEqually_"),  py::arg("verbose_"))
            .def("distributeValSolute_",&Perirhizal::distributeValSolute_,  py::arg("seg_values_content"), py::arg("volumes"), py::arg("source"),
                   py::arg("dt"))
            .def("distributeValWater_",&Perirhizal::distributeValWater_,  py::arg("seg_values_perVol"), py::arg("volumes"), py::arg("source"),
                   py::arg("dt"), py::arg("theta_S"), py::arg("theta_wilting_point"))
            .def("splitSoilVals_",&Perirhizal::splitSoilVals,  py::arg("soilVals"), py::arg("cellIds"), py::arg("isWater"),
                   py::arg("seg_values"), py::arg("seg_volume"), py::arg("dt"), py::arg("theta_S"), py::arg("theta_wilting_point"))
			.def("sumSegFluxes",&Perirhizal::sumSegFluxes)
			.def("splitSoilFluxes",&Perirhizal::splitSoilFluxes, py::arg("soilFluxes"), py::arg("type") = 0)
            .def_readwrite("ms",  &Perirhizal::ms);

    /*
     * XylemFlux.h
     */
    py::class_<XylemFlux, std::shared_ptr<XylemFlux>>(m, "XylemFlux")
            .def(py::init<std::shared_ptr<CPlantBox::MappedSegments>>())
            .def(py::init<std::shared_ptr<CPlantBox::MappedPlant>>())
            .def("setKr",py::overload_cast<std::vector<double>, std::vector<double>, bool> (&XylemFlux::setKr),
					py::arg("values"), py::arg("age") = std::vector<double>(0), py::arg("verbose")=false)
            .def("setKx",py::overload_cast<std::vector<double>, std::vector<double>, bool> (&XylemFlux::setKx),
                py::arg("values"), py::arg("age") = std::vector<double>(0), py::arg("verbose")=false)
            .def("setKr",py::overload_cast<std::vector<std::vector<double>>,std::vector<std::vector<double>>, double, bool> (&XylemFlux::setKr),
                py::arg("values"), py::arg("age") = std::vector<std::vector<double>>(0),py::arg("kr_length_") = -1.0, py::arg("verbose")=false)
            .def("setKx",py::overload_cast<std::vector<std::vector<double>>,std::vector<std::vector<double>>, bool> (&XylemFlux::setKx),
                py::arg("values"), py::arg("age") = std::vector<std::vector<double>>(0), py::arg("verbose")=false)
            .def("setKrTables",py::overload_cast<std::vector<std::vector<double>>, std::vector<std::vector<double>>, bool, bool> (&XylemFlux::setKrTables),
                py::arg("values"), py::arg("age"), py::arg("verbose")=false, py::arg("ageBased")=true)
            .def("setKxTables",py::overload_cast<std::vector<std::vector<double>>, std::vector<std::vector<double>>, bool> (&XylemFlux::setKxTables),
                py::arg("values"), py::arg("age"), py::arg("verbose")=false)
            .def("setKrTables",py::overload_cast<std::vector<std::vector< std::vector<double> >>, std::vector<std::vector<std::vector<double>>>, bool, bool> (&XylemFlux::setKrTables),
					py::arg("values"), py::arg("age"), py::arg("verbose")=false, py::arg("ageBased")=true)
            .def("setKxTables",py::overload_cast< std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<std::vector<double>>>, bool> (&XylemFlux::setKxTables), py::arg("values"), py::arg("age"), py::arg("verbose")=false)
            .def("setKrValues", &XylemFlux::setKrValues)
            .def("setKxValues", &XylemFlux::setKxValues)
            .def("getEffKr", &XylemFlux::getEffKr)
            .def("getKr", &XylemFlux::getKr)
            .def("getKx", &XylemFlux::getKx)
            .def("linearSystem",&XylemFlux::linearSystem, py::arg("simTime") , py::arg("sx") , py::arg("cells") = true,
                    py::arg("soil_k") = std::vector<double>(), py::arg("verbose")=false)
            .def("soilFluxes",&XylemFlux::soilFluxes, py::arg("simTime"), py::arg("rx"), py::arg("sx"), py::arg("approx") = false,
                    py::arg("soil_k") = std::vector<double>())
            .def("segFluxes",&XylemFlux::segFluxes, py::arg("simTime"), py::arg("rx"), py::arg("sx"), py::arg("approx") = false,
                    py::arg("cells") = false, py::arg("soil_k") = std::vector<double>(), py::arg("verbose")=false)
            .def("sumSegFluxes",&XylemFlux::sumSegFluxes)
            .def("splitSoilFluxes",&XylemFlux::splitSoilFluxes, py::arg("soilFluxes"), py::arg("type") = 0)
            .def("kr_f_cpp", &XylemFlux::kr_f_wrapped, py::arg("seg_ind"), py::arg("age"), py::arg("st"), py::arg("ot"), py::arg("cells")=false)
            .def_readonly("kx_f_cpp", &XylemFlux::kx_f)
            .def_readwrite("aI", &XylemFlux::aI)
            .def_readwrite("aJ", &XylemFlux::aJ)
            .def_readwrite("aV", &XylemFlux::aV)
            .def_readwrite("aB", &XylemFlux::aB)
            .def_readwrite("kr", &XylemFlux::kr)
            .def_readwrite("kx", &XylemFlux::kx)
            .def_readwrite("rs", &XylemFlux::rs)
            .def_readwrite("psi_air", &XylemFlux::psi_air);
    /*
     * PlantHydraulicParameters.h
     */
    py::class_<PlantHydraulicParameters, std::shared_ptr<PlantHydraulicParameters>>(m, "PlantHydraulicParameters")
            .def(py::init<>())
            .def(py::init<std::shared_ptr<CPlantBox::MappedSegments>>())
            .def("setMode", &PlantHydraulicParameters::setMode)
            .def("setKrConst", &PlantHydraulicParameters::setKrConst, py::arg("v"), py::arg("subType"), py::arg("organType") = Organism::ot_root, py::arg("kr_length") = -1.)
            .def("setKxConst", &PlantHydraulicParameters::setKxConst, py::arg("v"), py::arg("subType"), py::arg("organType") = Organism::ot_root)
            .def("setKrAgeDependent", &PlantHydraulicParameters::setKrAgeDependent, py::arg("age"), py::arg("values"), py::arg("subType"), py::arg("organType")= Organism::ot_root)
            .def("setKxAgeDependent", &PlantHydraulicParameters::setKxAgeDependent, py::arg("age"), py::arg("values"), py::arg("subType"), py::arg("organType")= Organism::ot_root)
            .def("setKrDistanceDependent", &PlantHydraulicParameters::setKrDistanceDependent, py::arg("distance"), py::arg("values"), py::arg("subType"), py::arg("organType")= Organism::ot_root)
            .def("setKxDistanceDependent", &PlantHydraulicParameters::setKxDistanceDependent, py::arg("distance"), py::arg("values"), py::arg("subType"), py::arg("organType")= Organism::ot_root)
            .def("setKrValues", &PlantHydraulicParameters::setKrValues)
            .def("setKxValues", &PlantHydraulicParameters::setKxValues)
            .def("getKr", &PlantHydraulicParameters::getKr)
            .def("getEffKr", &PlantHydraulicParameters::getEffKr)
            .def("getKx", &PlantHydraulicParameters::getKx)
            .def_readonly("kr_f", &PlantHydraulicParameters::kr_f)
            .def_readonly("kx_f", &PlantHydraulicParameters::kx_f)
            //.def_readonly("kr_f_wrapped", &PlantHydraulicParameters::kr_f_wrapped)
            .def_readonly("krMode", &PlantHydraulicParameters::krMode)
            .def_readonly("kxMode", &PlantHydraulicParameters::kxMode)
            .def_readonly("maxSubTypes", &PlantHydraulicParameters::maxSubTypes)
            .def_readwrite("krValues", &PlantHydraulicParameters::krValues)
            .def_readwrite("kxValues", &PlantHydraulicParameters::kxValues)
            .def_readwrite("kr_ages", &PlantHydraulicParameters::kr_ages)
            .def_readwrite("kr_values", &PlantHydraulicParameters::kr_values)
            .def_readwrite("kx_ages", &PlantHydraulicParameters::kx_ages)
            .def_readwrite("kx_values", &PlantHydraulicParameters::kx_values)
            .def_readwrite("ms", &PlantHydraulicParameters::ms)
            .def_readwrite("psi_air", &PlantHydraulicParameters::psi_air);
        /*
         * PlantHydraulicModel.h
         */
        py::class_<PlantHydraulicModel, std::shared_ptr<PlantHydraulicModel>>(m, "PlantHydraulicModel")
            .def(py::init<std::shared_ptr<MappedSegments>, std::shared_ptr<PlantHydraulicParameters>>())
            .def("linearSystemMeunier",&PlantHydraulicModel::linearSystemMeunier, py::arg("simTime") , py::arg("sx") , py::arg("cells") = true, py::arg("soil_k") = std::vector<double>())
            .def("getRadialFluxes", &PlantHydraulicModel::getRadialFluxes)
            .def("sumSegFluxes", &PlantHydraulicModel::sumSegFluxes)
            .def_readwrite("ms", &PlantHydraulicModel::ms)
            .def_readwrite("params", &PlantHydraulicModel::params)
            .def_readwrite("aI", &PlantHydraulicModel::aI)
            .def_readwrite("aJ", &PlantHydraulicModel::aJ)
            .def_readwrite("aV", &PlantHydraulicModel::aV)
            .def_readwrite("aB", &PlantHydraulicModel::aB);

	/*
     * Photosynthesis.h
     */
   py::class_<Photosynthesis, PlantHydraulicModel, std::shared_ptr<Photosynthesis>>(m, "Photosynthesis")
            .def(py::init<std::shared_ptr<CPlantBox::MappedPlant>, std::shared_ptr<PlantHydraulicParameters>,double, double>(),  py::arg("plant"),py::arg("params"),  py::arg("psiXylInit")=-500.0 ,  py::arg("ciInit")= 350e-6)
			.def("solve_photosynthesis",&Photosynthesis::solve_photosynthesis,py::arg("sim_time"),
					py::arg("sxx")  , py::arg("ea"),py::arg("es"),py::arg("TleafK")  ,
					py::arg("cells") = true,py::arg("soil_k") = std::vector<double>(),
					py::arg("doLog")=false, py::arg("verbose")=true,    py::arg("outputDir")="")
            .def_readwrite("PhotoType", &Photosynthesis::PhotoType)
            .def_readwrite("psiXyl_old", &Photosynthesis::psiXyl_old)
            .def_readwrite("psiXyl4Phloem", &Photosynthesis::psiXyl4Phloem)
            .def_readwrite("maxErrLim", &Photosynthesis::maxErrLim)
            .def_readwrite("maxErrAbsLim", &Photosynthesis::maxErrAbsLim)
            .def_readwrite("maxErr", &Photosynthesis::maxErr)
            .def_readwrite("maxErrAbs", &Photosynthesis::maxErrAbs)
			.def_readwrite("psiXyl", &Photosynthesis::psiXyl)
            .def_readwrite("An", &Photosynthesis::An)
            .def_readwrite("kp25", &Photosynthesis::kp25)
            .def_readwrite("kp", &Photosynthesis::kp)
            .def_readwrite("Vp", &Photosynthesis::Vp)
            .def_readwrite("Vc", &Photosynthesis::Vc)
            .def_readwrite("Vcrefmax", &Photosynthesis::Vcrefmax)
            .def_readwrite("Jrefmax", &Photosynthesis::Jrefmax)
            .def_readwrite("Jmax", &Photosynthesis::Jmax)
            .def_readwrite("J", &Photosynthesis::J)
            .def_readwrite("Vj", &Photosynthesis::Vj)
            .def_readwrite("fw", &Photosynthesis::fw)
            .def_readwrite("fwr", &Photosynthesis::fwr)
            .def_readwrite("fw_cutoff", &Photosynthesis::fw_cutoff)
            .def_readwrite("sh", &Photosynthesis::sh)
            .def_readwrite("p_lcrit", &Photosynthesis::p_lcrit)
            .def_readwrite("ci", &Photosynthesis::ci)
            .def_readwrite("deltagco2", &Photosynthesis::deltagco2)
            .def_readwrite("delta", &Photosynthesis::delta)
            .def_readwrite("oi", &Photosynthesis::oi)
            .def_readwrite("Rd", &Photosynthesis::Rd)
            .def_readwrite("gco2", &Photosynthesis::gco2)
            .def_readwrite("es", &Photosynthesis::es)
            .def_readwrite("ea", &Photosynthesis::ea)
            .def_readwrite("gm",&Photosynthesis::gm)
            .def_readwrite("PVD",&Photosynthesis::PVD)
            .def_readwrite("EAL",&Photosynthesis::EAL)
            .def_readwrite("hrelL",&Photosynthesis::hrelL)
            .def_readwrite("pg",&Photosynthesis::pg)
            .def_readwrite("Qlight", &Photosynthesis::Qlight)
            .def_readwrite("Jw", &Photosynthesis::Jw)
            .def_readwrite("Ev", &Photosynthesis::Ev)
            .def_readwrite("plant", &Photosynthesis::plant)
            .def_readwrite("Ag4Phloem", &Photosynthesis::Ag4Phloem)
            .def_readwrite("minLoop", &Photosynthesis::minLoop)
            .def_readwrite("maxLoop", &Photosynthesis::maxLoop)
            .def_readwrite("loop", &Photosynthesis::loop)
            .def_readwrite("Patm", &Photosynthesis::Patm)
            .def_readwrite("cs", &Photosynthesis::cs)
            .def_readwrite("TleafK", &Photosynthesis::TleafK)
            .def_readwrite("Chl", &Photosynthesis::Chl)
            .def_readwrite("g0", &Photosynthesis::g0)
            .def_readwrite("g_bl", &Photosynthesis::g_bl)
            .def_readwrite("g_canopy", &Photosynthesis::g_canopy)
            .def_readwrite("g_air", &Photosynthesis::g_air)
            .def_readwrite("theta", &Photosynthesis::theta)
            .def_readwrite("alpha", &Photosynthesis::alpha)
            .def_readwrite("a1", &Photosynthesis::a1)
            .def_readwrite("a2_stomata", &Photosynthesis::a2_stomata)
            .def_readwrite("a2_bl", &Photosynthesis::a2_bl)
            .def_readwrite("a2_canopy", &Photosynthesis::a2_canopy)
            .def_readwrite("a2_air", &Photosynthesis::a2_air)
            .def_readwrite("a3", &Photosynthesis::a3)
            .def_readwrite("Rd_ref", &Photosynthesis::Rd_ref)
            .def_readwrite("Kc_ref", &Photosynthesis::Kc_ref)
            .def_readwrite("Q10_photo", &Photosynthesis::Q10_photo)
            .def_readwrite("VcmaxrefChl1", &Photosynthesis::VcmaxrefChl1)
            .def_readwrite("VcmaxrefChl2", &Photosynthesis::VcmaxrefChl2)
            .def_readwrite("outputFlux", &Photosynthesis::outputFlux)
            .def_readwrite("outputFlux_old", &Photosynthesis::outputFlux_old)
            .def_readwrite("doLog", &Photosynthesis::doLog)
            .def_readonly("rho_h2o", &Photosynthesis::rho_h2o)
            .def_readonly("R_ph", &Photosynthesis::R_ph)
            .def_readonly("Mh2o", &Photosynthesis::Mh2o);

    py::enum_<Photosynthesis::PhotoTypes>(m, "PhotoTypes")
            .value("C3", Photosynthesis::PhotoTypes::C3)
            .value("C4", Photosynthesis::PhotoTypes::C4)
            .export_values();

	/*
     * runPM.h
     */
    py::class_<PhloemFlux, Photosynthesis, std::shared_ptr<PhloemFlux>>(m, "PhloemFlux")
            .def(py::init<std::shared_ptr<CPlantBox::MappedPlant>, std::shared_ptr<PlantHydraulicParameters>, double, double>(),  py::arg("plant"),py::arg("params"),
			py::arg("psiXylInit"),  py::arg("ciInit") )
            .def("waterLimitedGrowth",&PhloemFlux::waterLimitedGrowth)
            .def("setKr_st",&PhloemFlux::setKr_st, py::arg("values"), py::arg("kr_length") = -1.0, py::arg("verbose") = false)
            .def("setKx_st",&PhloemFlux::setKx_st, py::arg("values"), py::arg("verbose") = false)
            .def("setRmax_st",&PhloemFlux::setRmax_st, py::arg("values"), py::arg("verbose") = false)
            .def("setAcross_st",&PhloemFlux::setAcross_st, py::arg("values"), py::arg("verbose") = false)
            .def("setRhoSucrose",&PhloemFlux::setRhoSucrose, py::arg("values"), py::arg("verbose") = false)
            .def("setKrm1",&PhloemFlux::setKrm1, py::arg("values"), py::arg("verbose") = false)
            .def("setKrm2",&PhloemFlux::setKrm2, py::arg("values"), py::arg("verbose") = false)
			.def("startPM",&PhloemFlux::startPM)
            .def_readonly("rhoSucrose_f",&PhloemFlux::rhoSucrose_f)
            .def_readwrite("psiMin", &PhloemFlux::psiMin)
            .def_readwrite("Q_out",&PhloemFlux::Q_outv)
            .def_readwrite("Q_init",&PhloemFlux::Q_init)
            .def_readwrite("Q_out_dot",&PhloemFlux::Q_out_dotv)
            .def_readwrite("a_ST",&PhloemFlux::a_STv)
            .def_readwrite("vol_ST",&PhloemFlux::vol_STv)
            .def_readwrite("vol_Meso",&PhloemFlux::vol_Mesov)
            .def_readwrite("CSTimin",&PhloemFlux::CSTimin)
            .def_readwrite("C_ST",&PhloemFlux::C_STv)
            .def_readwrite("r_ST_ref",&PhloemFlux::r_ST_refv)
            .def_readwrite("r_ST",&PhloemFlux::r_STv)
            .def_readwrite("update_viscosity",&PhloemFlux::update_viscosity_)
            .def_readwrite("AgPhl",&PhloemFlux::Agv)
            .def_readwrite("Q_Grmax",&PhloemFlux::Q_Grmaxv)
            .def_readwrite("Q_Exudmax",&PhloemFlux::Q_Exudmaxv)
            .def_readwrite("Q_Rmmax",&PhloemFlux::Q_Rmmaxv)
            .def_readwrite("Fl",&PhloemFlux::Flv)
            //.def_readwrite("KMgr",&PhloemFlux::KMgr)
            .def_readwrite("KMfu",&PhloemFlux::KMfu)
            .def_readwrite("CsoilDefault",&PhloemFlux::CsoilDefault)
            .def_readwrite("Csoil_seg",&PhloemFlux::Csoil_seg)
            .def_readwrite("Csoil_node",&PhloemFlux::Csoil_node)
            .def_readwrite("deltaSucOrgNode",&PhloemFlux::deltaSucOrgNode_)
            .def_readwrite("usePsiXyl",&PhloemFlux::usePsiXyl)
            //.def_readwrite("expression",&PhloemFlux::expression)
            .def_readwrite("JW_ST",&PhloemFlux::JW_STv)
            .def_readwrite("Gr_Y",&PhloemFlux::Gr_Y)
            .def_readwrite("solver",&PhloemFlux::solver)
            .def_readwrite("atol",&PhloemFlux::atol_double)
            .def_readwrite("rtol",&PhloemFlux::rtol_double)
            .def_readwrite("surfMeso",&PhloemFlux::surfMeso)
			.def_readwrite("sameVolume_meso_st",&PhloemFlux::sameVolume_meso_st)
			.def_readwrite("sameVolume_meso_seg",&PhloemFlux::sameVolume_meso_seg)
			//.def_readwrite("Cobj_ST",&PhloemFlux::Cobj_ST)
			.def_readwrite("Vmaxloading",&PhloemFlux::Vmaxloading)
			.def_readwrite("useCWGr",&PhloemFlux::useCWGr)
			.def_readwrite("beta_loading",&PhloemFlux::beta_loading)
			.def_readwrite("Mloading",&PhloemFlux::Mloading)
			.def_readwrite("withInitVal",&PhloemFlux::withInitVal)
			.def_readwrite("initValST",&PhloemFlux::initValST)
			.def_readwrite("initValMeso",&PhloemFlux::initValMeso)
			.def_readwrite("doTroubleshooting",&PhloemFlux::doTroubleshooting)
			//.def_readwrite("Q_GrUnbornv_i",&PhloemFlux::Q_GrUnbornv_i)
			//.def_readwrite("Q_GrmaxUnbornv_i",&PhloemFlux::Q_GrmaxUnbornv_i)
			.def_readwrite("Fpsi",&PhloemFlux::Fpsi)
			.def_readwrite("Flen",&PhloemFlux::Flen)
			.def_readwrite("Q10",&PhloemFlux::Q10)
			.def_readwrite("TrefQ10",&PhloemFlux::TrefQ10)
			.def_readwrite("leafGrowthZone",&PhloemFlux::leafGrowthZone)
			.def_readwrite("StemGrowthPerPhytomer",&PhloemFlux::StemGrowthPerPhytomer)
			.def_readwrite("GrowthZone",&PhloemFlux::GrowthZone)
			.def_readwrite("GrowthZoneLat",&PhloemFlux::GrowthZoneLat)
			.def_readwrite("psi_osmo_proto",&PhloemFlux::psi_osmo_proto)
			.def_readwrite("psi_p_symplasm",&PhloemFlux::psi_p_symplasm)
            .def_readwrite("C_targ",&PhloemFlux::C_targ)
            .def_readwrite("C_targMesophyll",&PhloemFlux::C_targMesophyll)
            .def_readwrite("Vmax_S_ST",&PhloemFlux::Vmax_S_ST)
            .def_readwrite("kM_S_ST",&PhloemFlux::kM_S_ST)
            .def_readwrite("kHyd_S_ST",&PhloemFlux::kHyd_S_ST)
            .def_readwrite("k_mucil",&PhloemFlux::k_mucil)
            .def_readwrite("k_S_ST",&PhloemFlux::k_S_ST)
            .def_readwrite("Vmax_S_Mesophyll",&PhloemFlux::Vmax_S_Mesophyll)
            .def_readwrite("kM_S_Mesophyll",&PhloemFlux::kM_S_Mesophyll)
            .def_readwrite("kHyd_S_Mesophyll",&PhloemFlux::kHyd_S_Mesophyll)
            .def_readwrite("k_S_Mesophyll",&PhloemFlux::k_S_Mesophyll)
            .def_readwrite("k_mucil_",&PhloemFlux::k_mucil_)
            .def_readwrite("kr_st",&PhloemFlux::kr_st)
            .def_readwrite("kx_st",&PhloemFlux::kx_st)
            .def_readwrite("Across_st",&PhloemFlux::Across_st)
            .def_readwrite("Perimeter_st",&PhloemFlux::Perimeter_st)
            .def_readwrite("Rmax_st",&PhloemFlux::Rmax_st)
            .def_readwrite("rhoSucrose",&PhloemFlux::rhoSucrose)
            .def_readwrite("krm1v",&PhloemFlux::krm1v)
            .def_readwrite("krm2v",&PhloemFlux::krm2v);

    py::class_<PlantVisualiser, std::shared_ptr<PlantVisualiser>>(m, "PlantVisualiser")
        .def(py::init<>())
        .def(py::init<std::shared_ptr<MappedPlant>>())
        .def("ComputeGeometryForOrgan",&PlantVisualiser::ComputeGeometryForOrgan, py::arg("organId"))
        .def("ComputeGeometryForOrganType",&PlantVisualiser::ComputeGeometryForOrganType, py::arg("organType"), py::arg("clearFirst") = true)
        .def("ComputeGeometry",&PlantVisualiser::ComputeGeometry)
        .def("GetGeometry",&PlantVisualiser::GetGeometry)
        .def("GetGeometryColors",&PlantVisualiser::GetGeometryColors)
        .def("GetGeometryNormals",&PlantVisualiser::GetGeometryNormals)
        .def("GetGeometryIndices",&PlantVisualiser::GetGeometryIndices)
        .def("GetGeometryTextureCoordinates",&PlantVisualiser::GetGeometryTextureCoordinates)
        .def("GetGeometryNodeIds", &PlantVisualiser::GetGeometryNodeIds)
        .def("SetGeometryResolution",&PlantVisualiser::SetGeometryResolution, py::arg("resolution"))
        .def("SetLeafResolution",&PlantVisualiser::SetLeafResolution, py::arg("resolution"))
        .def("SetComputeMidlineInLeaf", &PlantVisualiser::SetComputeMidlineInLeaf, py::arg("inCompute"))
        .def("HasGeometry", &PlantVisualiser::HasGeometry)
        .def("ResetGeometry", &PlantVisualiser::ResetGeometry)
    ;

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
             .value("CWLim", Plant::GrowthFunctionTypes::gft_CWLim)
             .export_values();

    py::class_<ExudationModel, std::shared_ptr<ExudationModel>>(m, "ExudationModel")
            .def(py::init<double, double, int, std::shared_ptr<RootSystem>>())
            .def(py::init<double, double, double, int, int, int, std::shared_ptr<RootSystem>>())
            .def_readwrite("Q", &ExudationModel::Q)
            .def_readwrite("Dl", &ExudationModel::Dl)
            .def_readwrite("theta", &ExudationModel::theta)
            .def_readwrite("R", &ExudationModel::R)
            .def_readwrite("k", &ExudationModel::k)
            .def_readwrite("l", &ExudationModel::l)
            .def_readwrite("type", &ExudationModel::type)
            .def_readwrite("n0", &ExudationModel::n0)
            .def_readwrite("thresh13", &ExudationModel::thresh13)
            .def_readwrite("calc13", &ExudationModel::calc13)
            .def_readwrite("observationRadius", &ExudationModel::observationRadius)
            .def("calculate",  &ExudationModel::calculate, py::arg("tend"), py::arg("i0") = 0, py::arg("iend")=-1);
    py::enum_<ExudationModel::IntegrationType>(m, "IntegrationType")
            .value("mps_straight", ExudationModel::IntegrationType::mps_straight )
            .value("mps", ExudationModel::IntegrationType::mps )
            .value("mls", ExudationModel::IntegrationType::mls )
            .export_values();



}

}
