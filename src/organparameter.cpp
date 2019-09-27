// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "organparameter.h"

#include "Organism.h"

#include <limits>
#include <iostream>
#include <exception>

namespace CPlantBox {

/**
 * @return Quick info about the object for debugging
 */
std::string OrganSpecificParameter::toString() const
{
    std::stringstream str;
    str << "subType\t" << subType << std::endl;
    return str.str();
}

/**
 * Default constructor using default values for the parameters.
 *
 * @param plant     the organism managing this organ type parameter
 */
OrganRandomParameter::OrganRandomParameter(std::weak_ptr<Organism> p): plant(p)
{
    iparam["organType"] = &organType; // set up class introspection
    iparam["subType"] = &subType;
    description["name"]  = "Name of the sub type of the organ, e.g. small lateral";
    description["organType"]  = "Organ type (unspecified organ = 0, seed = 1, root = 2, stem = 3, leaf = 4)";
    description["subType"]  = "Unique identifier of this sub type";
}

/**
 * Copies the OrganTypeParameter and into another plant.
 *
 * @param plant     target organism managing this organ type parameter
 *
 * @return A new organ type parameter object, with same values as this, ownership is passed to the caller
 */
std::shared_ptr<OrganRandomParameter> OrganRandomParameter::copy(std::weak_ptr<Organism> p)
{
    auto o = std::make_shared<OrganRandomParameter>(p); // copy constructor would break class introspection
    o->name == name;
    o->organType = organType;
    o->subType = subType;
    return o;
}

/**
 * An organ specific parameter set (OrganParameter) is created from the general parameter set (OrganTypeParameter).
 * In the specializations of this class, random variables are drawn from their distributions.
 *
 * Normally, both OrganParameter, and OrganTypeParameter are constant over the organ life time.
 *
 * @return The organ specific parameter set (each organ has exactly one set), ownership is passed to caller
 */
std::shared_ptr<OrganSpecificParameter> OrganRandomParameter::realize()
{
    auto op = std::make_shared<OrganSpecificParameter>();
    op->subType = subType;
    return op;
}

/**
 * Returns a single scalar parameter.
 * Works for int and double in derived classes, if the parameters were added to the respective hashmaps
 *
 * @param name      name of the parameter, add "_dev" to obtain the parameter's deviation (usually standard deviation),
 *                  optionally, add "_mean" to obtain the mean value (to avoid naming conflicts with the specific parameters).
 *
 * @return A single parameter value, if unknown NaN
 */
double OrganRandomParameter::getParameter(std::string name) const
{
    if ((name.length()>4) && (name.substr(name.length()-4)=="_dev")) {// looking for standard deviation?
        std::string n = name.substr(0,name.length()-4);
        if (param_sd.count(n)) {
            return *param_sd.at(n);
        }
    }
    if ((name.length()>5) && (name.substr(name.length()-5)=="_mean")) {// looking the mean value?
        name = name.substr(0,name.length()-5);
    }
    if (iparam.count(name)) { // looking for an int parameter
        return (double)*iparam.at(name);
    }
    if (dparam.count(name)) { // looking for a double parameter
        return *dparam.at(name);
    }
    return std::numeric_limits<double>::quiet_NaN(); // default if name is unknown
}

/**
 * Quick info about the object for debugging
 *
 * @param verbose     if true, all parameters names, values, and descriptions are written
 */
std::string OrganRandomParameter::toString(bool verbose) const
{
    std::stringstream str;
    if (verbose) {
        str << "[Parameters of "<< name << "]" << std::endl;
        str << "Variable\t\t" << "Value\t\t" << "Deviation\t\t" << "Description" << std::endl;
        str << "===" << std::endl;
        for (auto& ip : iparam) { // int valued parameters
            std::string key = ip.first;
            std::string keystr = key;
            while (keystr.length()<18) {
                keystr.append(" ");
            }
            str << keystr << "\t" << *ip.second<< "\t\t";
            if (param_sd.count(key)) { // deviations
                str << *param_sd.at(key)<< "\t\t\t";
            } else {
                str << "-"<< "\t\t\t";
            }
            if (description.count(key)) { // descriptions
                str << description.at(key);
            } else {
                str << "-";
            }
            str << std::endl;
        }
        for (auto& dp : dparam) { // double valued parameters
            std::string key = dp.first;
            std::string keystr = key;
            while (keystr.length()<18) {
                keystr.append(" ");
            }
            str << keystr << "\t" << *dp.second<< "\t\t";
            if (param_sd.count(key)) { // deviations
                str << *param_sd.at(key)<< "\t\t\t";
            } else {
                str << "-"<< "\t\t\t";
            }
            if (description.count(key)) { // descriptions
                str << description.at(key);
            } else {
                str << "-";
            }
            str << std::endl;
        }
        str << "===" << std::endl;
        return str.str();
    } else {
        str << "Name: " << name << ", " << "organType: "<< organType << ", " << "subType" << ", " << subType;
        return str.str();
    }
}

/**
 * Reads a XML Tag with parameter values for this organ type parameter.
 * This works in derived classes for int and double parameters including deviation.
 *
 * The reader does not check if all parameters are in the XML Tag,
 * missing parameters are not altered.
 *
 * Not exposed to Python
 *
 * @param element  The XML element containing the parameter tags
 */
void OrganRandomParameter::readXML(tinyxml2::XMLElement* element)
{
    std::string type = element->Name();
    this->organType =  Organism::organTypeNumber(type);
    subType = element->IntAttribute("subType");
    name = element->Attribute("name");
    auto p = element->FirstChildElement("parameter");
    while(p) {
        std::string key = p->Attribute("name");
        if (iparam.count(key)>0) {
            *iparam[key] = p->IntAttribute("value");
        } else if (dparam.count(key)>0) {
            *dparam[key] = p->DoubleAttribute("value");
        }
        if (param_sd.count(key)>0) {
            *param_sd[key] = p->DoubleAttribute("dev", 0); // optional value
        }
        p = p->NextSiblingElement("parameter");
    }
}

/**
 * Reads a single XML Tag with parameter values for this organ type parameter from a file.
 * Mainly for testing and debugging.
 *
 * @param name      file name
 */
void OrganRandomParameter::readXML(std::string name)
{
    tinyxml2::XMLDocument doc;
    doc.LoadFile(name.c_str());
    if(doc.ErrorID() == 0) {
        tinyxml2::XMLElement* element = doc.FirstChildElement();
        readXML(element);
    } else {
        std::cout << "readXML(): could not open file\n" << std::flush;
    }
}

/**
 * Creates a XML tag representing this organ type parameter
 * This works in derived classes for int and double parameters including deviation.*
 *
 * Not exposed to Python
 *
 * @param doc         the XML document you are writing
 * @param comments    if true, add the parameter descriptions as comments
 *
 * @return the xml tag representing this organ type parameter
 */
tinyxml2::XMLElement* OrganRandomParameter::writeXML(tinyxml2::XMLDocument& doc, bool comments) const
{
    std::string typeName = Organism::organTypeName(organType);
    tinyxml2::XMLElement* organ = doc.NewElement(typeName.c_str());
    organ->SetAttribute("name", name.c_str());
    organ->SetAttribute("subType", subType);
    for (auto& ip : iparam) { // int valued parameters
        std::string key = ip.first;
        if (!(key.compare("subType")==0 || key.compare("organType")==0)) { // already written in organ attributes
            tinyxml2::XMLElement* p = doc.NewElement("parameter");
            p->SetAttribute("name", key.c_str());
            p->SetAttribute("value", *ip.second);
            if (param_sd.count(key)) { // deviations
                double d = *(param_sd.at(key));
                p->SetAttribute("dev", d);
            }
            organ->InsertEndChild(p);
            if (comments && description.count(key)) { // descriptions
                std::string str = description.at(key);
                tinyxml2::XMLComment* c = doc.NewComment(str.c_str());
                organ->InsertEndChild(c);
            }
        }
    }
    for (auto& dp : dparam) { // double valued parameters
        std::string key = dp.first;
        if (!(key.compare("subType")==0 || key.compare("organType")==0)) { // already written in organ attributes
            tinyxml2::XMLElement* p = doc.NewElement("parameter");
            p->SetAttribute("name", key.c_str());
            p->SetAttribute ("value", float(*dp.second)); // float output is much nicer
            if (param_sd.count(key)) { // deviations
                p->SetAttribute("dev", *(param_sd.at(key)));
            }
            organ->InsertEndChild(p);
            if (comments && description.count(key)) { // descriptions
                std::string str = description.at(key);
                tinyxml2::XMLComment* c = doc.NewComment(str.c_str());
                organ->InsertEndChild(c);
            }
        }
    }
    return organ;
}

/**
 * Writes a single XML tag representing this organ type parameter into a file.
 * Mainly for testing and debugging.
 *
 * @param name      file name
 */
void OrganRandomParameter::writeXML(std::string name) const
{
    tinyxml2::XMLDocument doc;
    tinyxml2::XMLElement* organ = this->writeXML(doc, true);
    doc.InsertEndChild(organ);
    doc.SaveFile(name.c_str());
}

/**
 *  Binds a parameter name to its int member variable,
 *  the values can then be accessed with OrganTypeParameter::getParameter,
 *  and are used in OrganTypeParameter::toString, and OrganTypeParameter::writeXML.
 *
 *  @param name         the parameter name
 *  @param i            a pointer to the integer member variable
 *  @param descr        a parameter description
 *  @param dev          optionally, a deviation of the parameter like standard deviation.
 */
void OrganRandomParameter::bindParameter(std::string name, int* i, std::string descr, double* dev)
{
    iparam[name] = i;
    if (!descr.empty()) {
        description[name] = descr;
    }
    if (dev!=nullptr) {
        param_sd[name] = dev;
    }
}

/**
 *  Binds a parameter name to its double member variable,
 *  the values can then be accessed with OrganTypeParameter::getParameter,
 *  and are used in OrganTypeParameter::toString, and  OrganTypeParameter::writeXML.
 *
 *  @param name         the parameter name
 *  @param i            a pointer to the integer member variable
 *  @param descr        a parameter description
 *  @param dev          optionally, a deviation of the parameter like standard deviation.
 */
void OrganRandomParameter::bindParameter(std::string name, double* d, std::string descr, double* dev)
{
    dparam[name] = d;
    if (!descr.empty()) {
        description[name] = descr;
    }
    if (dev!=nullptr) {
        param_sd[name] = dev;
    }
}

} // namespace
