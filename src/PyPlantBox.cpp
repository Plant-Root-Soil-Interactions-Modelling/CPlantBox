// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "pybind11/include/pybind11/pybind11.h"
namespace py = pybind11;

/**
 * A Python binding based on pybind11
 */

#include "mymath.h"
//#include "sdf.h"
//#include "soil.h"
#include "Organism.h"
#include "organparameter.h"

namespace CRootBox {

/**
 * bybind11 Python bindings
 */
PYBIND11_MODULE(py_plantbox, m) {
    /* mymath */
    py::class_<Vector2i>(m, "Vector2i")
        .def(py::init<>())
        .def(py::init<int, int>())
        .def(py::init<const Vector2i&>())
        .def("__str__", &Vector2i::toString)
        .def_readwrite("x", &Vector2i::x)
        .def_readwrite("y", &Vector2i::y);
    py::class_<Vector2d>(m, "Vector2d")
        .def(py::init<>())
        .def(py::init<double, double>())
        .def(py::init<const Vector2d&>())
        .def("__str__", &Vector2d::toString)
        .def_readwrite("x", &Vector2d::x)
        .def_readwrite("y", &Vector2d::y);
    py::class_<Vector3d>(m, "Vector3d")
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
        .def_readwrite("z", &Vector3d::z);
    py::class_<Matrix3d>(m, "Matrix3d")
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
        .def_readwrite("r2", &Matrix3d::r2);
//    /* organparameter.h  */
//    py::class_<OrganSpecificParameter>(m,"OrganSpecificParameter")
//        .def(py::init<>())
//        .def_readwrite("subType",&OrganSpecificParameter::subType)
//        .def("__str__",&OrganSpecificParameter::toString);
//    py::class_<OrganRandomParameter>(m,"OrganRandomParameter")
//        .def(py::init<Organism*>())
////        .def("copy",&OrganRandomParameter::copy, return_value_policy<reference_existing_object>())
////        .def("realize",&OrganRandomParameter::realize, return_value_policy<reference_existing_object>())
//        .def("getParameter",&OrganRandomParameter::getParameter)
////        .def("writeXML",writeXML1)
////        .def("readXML",readXML1)
////        .def("bindIntParameter",bindIntParameter, bindParameter_overloads()) // not working, can't pass int*
////        .def("bindDoubleParameter",bindDoubleParameter, bindParameter_overloads()) // not working, can't pass double*
//        .def_readwrite("name",&OrganRandomParameter::name)
//        .def_readwrite("organType",&OrganRandomParameter::organType)
//        .def_readwrite("subType",&OrganRandomParameter::subType)
////        .def("__str__",&OrganRandomParameter::toString, toString_overloads())
//        ;


    //    /* sdf */
//    py::class_<SignedDistanceFunction>(m,"SignedDistanceFunction")
//    .def(py::init<>());
//        //.def("getDist",&SignedDistanceFunction::getDist)
//        .def("writePVPScript", (std::string (SignedDistanceFunction::*)() const) &SignedDistanceFunction::writePVPScript) // overloads
// .def("getGradient",  &SignedDistanceFunction::getGradient, py::arg("p"), py::arg("eps") = 5.e-4) // defaults
      //  .def("__str__",&SignedDistanceFunction::toString);
//    py::class_<SDF_PlantBox, SignedDistanceFunction>(m, "SDF_PlantBox")
//        .def(py::init<double,double,double>());
////        .def("getDist",&SDF_PlantBox::getDist)
////        .def("__str__",&SDF_PlantBox::toString) // should be ok?
//    py::class_<SDF_PlantContainer, SignedDistanceFunction>(m,"SDF_PlantContainer")
////        .def(py::init<>())
//        .def(py::init<double,double,double,double>());
////        .def("getDist",&SDF_PlantContainer::getDist)
////        .def("__str__",&SDF_PlantContainer::toString);
//    py::class_<SDF_RotateTranslate, SignedDistanceFunction>(m, "SDF_RotateTranslate")
//        .def(py::init<SignedDistanceFunction*,double,int,Vector3d&>())
//        .def(py::init<SignedDistanceFunction*,Vector3d&>());
////        .def("getDist",&SDF_RotateTranslate::getDist)
////        .def("__str__",&SDF_RotateTranslate::toString);
//    py::enum_<SDF_RotateTranslate::SDF_Axes>(m, "SDF_Axis")
//        .value("xaxis", SDF_RotateTranslate::SDF_Axes::xaxis)
//        .value("yaxis", SDF_RotateTranslate::SDF_Axes::yaxis)
//        .value("zaxis", SDF_RotateTranslate::SDF_Axes::zaxis)
//        .export_values();;
//    py::class_<SDF_Intersection, SignedDistanceFunction>(m,"SDF_Intersection")
//            .def(py::init<std::vector<SignedDistanceFunction*>>())
//            .def(py::init<SignedDistanceFunction*,SignedDistanceFunction*>());
////        .def("getDist",&SDF_Intersection::getDist)
////        .def("__str__",&SDF_Intersection::toString)
////        ;
//    py::class_<SDF_Union, SDF_Intersection>(m, "SDF_Union")
//        .def(py::init<std::vector<SignedDistanceFunction*>>())
//        .def(py::init<SignedDistanceFunction*,SignedDistanceFunction*>());
////        .def("getDist",&SDF_Union::getDist)
////        .def("__str__",&SDF_Union::toString)
////        ;
//    py::class_<SDF_Difference, SDF_Intersection>(m, "SDF_Difference")
//        .def(py::init<std::vector<SignedDistanceFunction*>>())
//        .def(py::init<SignedDistanceFunction*,SignedDistanceFunction*>());
////        .def("getDist",&SDF_Difference::getDist)
////        .def("__str__",&SDF_Difference::toString)
////        ;
//    py::class_<SDF_Complement, SignedDistanceFunction>(m, "SDF_Complement")
//        .def(py::init<SignedDistanceFunction*>());
////        .def("getDist",&SDF_Complement::getDist)
////        .def("__str__",&SDF_Complement::toString)
////        ;
//    py::class_<SDF_HalfPlane, SignedDistanceFunction>(m, "SDF_HalfPlane")
//        .def(py::init<Vector3d&,Vector3d&>())
//        .def(py::init<Vector3d&,Vector3d&,Vector3d&>())
////        .def("getDist",&SDF_HalfPlane::getDist)
//        .def_readwrite("o", &SDF_HalfPlane::o)
//        .def_readwrite("n", &SDF_HalfPlane::n)
//        .def_readwrite("p1", &SDF_HalfPlane::p1)
//        .def_readwrite("p2", &SDF_HalfPlane::p2);
////        .def("__str__",&SDF_HalfPlane::toString)
//        ;
}

}
