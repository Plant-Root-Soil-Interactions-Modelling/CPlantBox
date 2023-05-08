// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "organparameter.h"

#include "Organism.h"
#include "growth.h"

#include <limits>
#include <iostream>
#include <exception>
#include <numeric>

namespace CPlantBox {

/**
 * @return Quick info about the object for debugging
 */
std::string OrganSpecificParameter::toString() const
{
    std::stringstream str;
    str << "subType " << subType << ", radius " << a << " cm " << std::endl;
    return str.str();
}

/**
 * Default constructor using default values for the parameters.
 *
 * @param plant     the organism managing this organ type parameter
 */
OrganRandomParameter::OrganRandomParameter(std::shared_ptr<Organism> p): plant(p)
{
    bindParameters();
    f_gf = std::make_shared<ExponentialGrowth>(); ///< growth function
}

/**
 * Copies the OrganTypeParameter and into another plant.
 *
 * @param plant     target organism managing this organ type parameter
 *
 * @return A new organ type parameter object, with same values as this, ownership is passed to the caller
 */
std::shared_ptr<OrganRandomParameter> OrganRandomParameter::copy(std::shared_ptr<Organism> p)
{
    // std::cout << "OrganRandomParameter::copy\n"<< std::flush;
    auto o = std::make_shared<OrganRandomParameter>(p); // copy constructor would break class introspection
    o->name = name;
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
    auto p = plant.lock();
    double a_ = std::max(a + p->randn()*as, 0.); // radius
    return std::make_shared<OrganSpecificParameter>(subType, a_);
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
	throw std::runtime_error("OrganRandomParameter::getParameter: parameter "+name+" not found. Send NaN as output.");
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
        str << "successor organ types\t";
        for (int i=0; i<successorOT.size(); i++) {
            for (int j=0; j<successorOT.at(i).size(); j++) {
				str << successorOT.at(i).at(j) << " ";
			}
			str << "; ";
        }
        str << "\t" << description.at("successorOT") << std::endl;
        str << "successor sub types\t";
        for (int i=0; i<successorST.size(); i++) {
            for (int j=0; j<successorST.at(i).size(); j++) {
				str << successorST.at(i).at(j) << " ";
			}
			str << "; ";
        }
        str << "\t" << description.at("successorST") << std::endl;
        str << "successor lateral number\t";
        for (int i=0; i<successorNo.size(); i++) {
            str << successorNo.at(i) << "; ";
        }
        str << "\t" << description.at("successorNo") << std::endl;
        str << "successorP\t";
        for (int i=0; i<successorP.size(); i++) {
            for (int j=0; j<successorP.at(i).size(); j++) {
				str << successorP.at(i).at(j) << " ";
			}
			str << "; ";
        }
        str << "\t" << description.at("successorP") << std::endl;
        str << "successorWhere\t";
        for (int i=0; i<successorWhere.size(); i++) {
            for (int j=0; j<successorWhere.at(i).size(); j++) {
				str << successorWhere.at(i).at(j) << " ";
			}
			str << "; ";
        }
        str << "\t" << description.at("successorWhere") << std::endl;

        return str.str();
    } else {
        str << "name: " << name << ", " << "organType: "<< organType << ", " << "subType: " << subType << ".";
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
 * Not directly exposed to Python (called by Organism::readParameters)
 *
 * @param element  The XML element containing the parameter tags
 */
void OrganRandomParameter::readXML(tinyxml2::XMLElement* element, bool verbose)
{
    subType = element->IntAttribute("subType", 0);
    const char* cname = element->Attribute("name");
    if (cname!=nullptr) {
        name = std::string(cname);
    } else {
        name = "undefined*";
    }
	successorOT.resize(0, std::vector<int>(0));//2D, int
    successorST.resize(0, std::vector<int>(0));//2D, int
    successorP.resize(0, std::vector<double>(0));//2D, double
    successorNo.resize(0);//1D, int
    successorWhere.resize(0, std::vector<double>(0));//2D, double
    auto p = element->FirstChildElement("parameter");
    while(p!=nullptr) {
        const char* str = p->Attribute("name");
        if (str!=nullptr) {
            std::string key = std::string(str);
            int i = 0;
            if (iparam.count(key)>0) {
                *iparam[key] = p->IntAttribute("value", *iparam[key]);  // if not found, leave default
                i++;
            } else if (dparam.count(key)>0) {
                *dparam[key] = p->DoubleAttribute("value", *dparam[key]);  // if not found, leave default
                i++;
            }
            if (param_sd.count(key)>0) {
                *param_sd[key] = p->DoubleAttribute("dev", *param_sd[key]); // if not found, leave default
                i++;
            }
			if (key.compare("successor")==0) {
				readSuccessor(p, verbose);
				i++;
			}
            if (i == 0) {
                if ((key.compare("leafGeometry")!=0)||(organType!=Organism::ot_leaf)) {
					if(verbose){
						std::cout << "OrganRandomParameter::readXML: warning! parameter " << key <<
                        " is defined in the xml, but not available in organ " << Organism::organTypeName(organType) << "\n" << std::flush;
					}
                }
            }
            p = p->NextSiblingElement("parameter");
        } else {
            std::cout << "OrganRandomParameter::readXML: warning! tag has no attribute 'name' \n" << std::flush;
            p = p->NextSiblingElement("parameter");
        }
    }
}

/**
 * Reads a single XML Tag with parameter values for this organ type parameter from a file.
 * Mainly for testing and debugging.
 *
 * @param name      file name
 */
void OrganRandomParameter::readXML(std::string name, bool verbose)
{
    tinyxml2::XMLDocument doc;
    doc.LoadFile(name.c_str());
    if(doc.ErrorID() == 0) {
        tinyxml2::XMLElement* element = doc.FirstChildElement();
        readXML(element, verbose);
    } else {
        std::cout << "readXML(): could not open file\n" << std::flush;
    }
}

void OrganRandomParameter::readSuccessor(tinyxml2::XMLElement* p, bool verbose)
{
	double p_ , defaultVald;
	int defaultVal, defaultSize;
	std::vector<std::string> lookfor ;
	bool replaceByDefaultValue;
	std::string key = p->Attribute("name");
	if (key.compare("successor")==0)  {

		int ruleId = std::max(0,int(successorOT.size()-1));
		int	success = p->QueryIntAttribute("ruleId",&ruleId);
		if(success != tinyxml2::XML_SUCCESS){
			success = p->QueryIntAttribute("number",&ruleId);
			if(success != tinyxml2::XML_SUCCESS){
				if(verbose)
				{
					std::cout<<"OrganRandomParameter::readSuccessor: for parameter of organ "+std::to_string(organType)
					+", subType "+std::to_string(subType)+", 'ruleId' (and 'number') not found in successor definition. ";
					std::cout<<"Use defeault ruleId instead: "+ std::to_string(ruleId)+"\n";
				}
			}
		}

		if(ruleId>= successorNo.size() )//create new rules
		{
			int toAdd = ruleId + 1;
			successorOT.resize(toAdd,std::vector<int>());
			successorST.resize(toAdd,std::vector<int>());
			successorWhere.resize(toAdd,std::vector<double>());
			successorP.resize(toAdd,std::vector<double>());
			successorNo.resize(toAdd,1);//default == make one lateral
		}

		int numLat;
		success = p->QueryIntAttribute("numLat",&numLat);
		if(success == tinyxml2::XML_SUCCESS){successorNo.at(ruleId) = numLat;}


		//default == empty vector == apply rule to all the linking nodes
		replaceByDefaultValue = true;//replace by default value if not found (if false: throw an error)
		lookfor = std::vector<std::string>{"where"};//parameter name
		defaultVald = -0.0;
		defaultSize = 0;//how many time do we need to repeat the value. Here: leave vectore empty == apply everywhere

		cpb_queryStringAttribute(lookfor,
					defaultVald,defaultSize, replaceByDefaultValue,
					successorWhere.at(ruleId), p);//name, default value, vector to fill, accept not found

		replaceByDefaultValue = false;lookfor = std::vector<std::string>{"subType","subtype","type"};//subtype (or type for backarwad compatibility)
		defaultVal = 1.0;defaultSize = 0;
		cpb_queryStringAttribute(lookfor,
					defaultVal,defaultSize, replaceByDefaultValue,
					successorST.at(ruleId), p);



		replaceByDefaultValue = true;
		lookfor = std::vector<std::string>{"probability","percentage"};//probability (or percentage for backward compatibility)
		defaultVald = successorNo.at(ruleId)/successorST.at(ruleId).size();//default = number of laterals / option for laterals
		defaultSize = (successorST.at(ruleId).size() - successorP.at(ruleId).size());
		if(!std::isfinite(defaultVald))
		{
			for(int k = 0;k < successorST.at(ruleId).size(); k++){
				std::cout<<successorST.at(ruleId).at(k)<<" ";}
			std::cout<<std::endl;
			throw std::runtime_error("OrganRandomParameter::readSuccessor !std::isfinite(default growth probability) "+
			std::to_string(defaultVald)+" , successorST.at(ruleId).size() == "+ std::to_string(successorST.at(ruleId).size()));
		}
		cpb_queryStringAttribute(lookfor,
					defaultVald,defaultSize, replaceByDefaultValue,
					successorP.at(ruleId), p);

		replaceByDefaultValue = true;lookfor = std::vector<std::string>{"organType","organtype"};
		if((successorST.at(ruleId).at(0) == 2)&&(this->organType == Organism::ot_stem)){
			if(verbose)
			{
				std::cout<<"OrganRandomParameter::readSuccessor: gave a stem a successor of subtype 2 and did not specify type.";
				std::cout<<" For backward compatibility, this will be considered as a leaf successor"<<std::endl;
			}
			//for backward compatibility =>
			//if (no organtype given) + (parent is stem) + (subtype == 2) == we want a leaf
			defaultVal = Organism::ot_leaf;
		}else{defaultVal = this->organType;}
		defaultSize = (successorST.at(ruleId).size() - successorOT.at(ruleId).size());
		cpb_queryStringAttribute(lookfor,
					defaultVal,defaultSize, replaceByDefaultValue,
					successorOT.at(ruleId), p);

		//sum(p_) <= 1.
		p_ = std::accumulate(successorP.at(ruleId).begin(), successorP.at(ruleId).end(), 0.);
		if(p_ > (1. + 1e-6)) //can be < 1 but not > 1
		{
			throw std::runtime_error("OrganRandomParameter::readSuccessor: probability of emergence for rule "+ std::to_string(ruleId)
			+" is > 1: "+ std::to_string(p_));
		}
	}
}




template <class IntOrDouble>
void OrganRandomParameter::cpb_queryStringAttribute(std::vector<std::string> keyNames,IntOrDouble defaultVal,int sizeVector,
	bool replaceByDefault,
	std::vector<IntOrDouble> & vToFill, tinyxml2::XMLElement* key)
{
	int success = -1;
	std::vector<IntOrDouble> dummy;
	for(int i = 0; (i < keyNames.size())&&(tinyxml2::XML_SUCCESS != success); i++){
		const char* cckey;
		success = key->QueryStringAttribute(keyNames.at(i).c_str(),&cckey);
		if (tinyxml2::XML_SUCCESS == success){
			//use default val to defin type of element in vector
			dummy = string2vector(cckey, defaultVal);
		}

	}
	if(tinyxml2::XML_SUCCESS != success){
		assert(replaceByDefault &&
		"mymath::queryStringAttribute: key not found in xml file without default value");
		dummy = std::vector<IntOrDouble>(sizeVector, defaultVal);
	};
	if(vToFill.size() == 0){
	    vToFill = dummy;
	} else {
		vToFill.insert( vToFill.end(), dummy.begin(), dummy.end() );
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
                if (d!=0) {
                    p->SetAttribute("dev", float(d));
                }
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
                double d = *(param_sd.at(key));
                if (d!=0) {
                    p->SetAttribute("dev", float(d));
                }
            }
            organ->InsertEndChild(p);
            if (comments && description.count(key)) { // descriptions
                std::string str = description.at(key);
                tinyxml2::XMLComment* c = doc.NewComment(str.c_str());
                organ->InsertEndChild(c);
            }
        }
    }
	double p_ = 0.;
	for (int j = 0; j<successorST.size(); j++) {
			tinyxml2::XMLElement* p = doc.NewElement("parameter");
			p->SetAttribute("name", "successor");
			p->SetAttribute("ruleId",j);
			if(successorNo.size()>j){p->SetAttribute("numLat", successorNo.at(j));}
			if(successorWhere.size()>j){p->SetAttribute("Where", vector2string(successorWhere.at(j)).c_str());}
			p->SetAttribute("subType", vector2string(successorST.at(j)).c_str());
			if(successorOT.size()>j){p->SetAttribute("organType", vector2string(successorOT.at(j)).c_str());}
			p->SetAttribute("percentage", vector2string(successorP.at(j)).c_str());
			organ->InsertEndChild(p);
			if (comments) {
				std::string str = description.at("successorST");
				tinyxml2::XMLComment* c = doc.NewComment(str.c_str());
				organ->InsertEndChild(c);
			}
			p_ += std::accumulate(successorP.at(j).begin(), successorP.at(j).end(), 0.);
			if (p_>1) {
				std::cout << "OrganRandomParameter::writeXML: Warning! percentages "<< p_ <<"  > 1.\n";
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
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void OrganRandomParameter::bindParameters()
{
    bindParameter("organType", &organType, "Organ type (unspecified organ = 0, seed = 1, root = 2, stem = 3, leaf = 4)");
    bindParameter("subType", &subType, "Unique identifier of this sub type");
    bindParameter("a", &a, "radius [cm]", &as);
    bindParameter("dx", &dx, "Axial resolution [cm] (maximal segment size)");
    bindParameter("dxMin", &dxMin, "Axial resolution [cm] (minimal segment size)");
	bindParameter("ldelay", &ldelay, "Lateral emergence delay [day]", &ldelays);
    // other parameters (descriptions only)
    description["name"]  = "Name of the sub type of the organ, e.g. small lateral";
    // other parameters (descriptions only)
    description["successorNo"] = "number of lateral ";
    description["successorWhere"] = "linking node id of the lateral (no value = apply everywhere)";
    description["successorOT"] = "organ type of lateral ";
    description["successorST"] = "Sub type of lateral ";
    description["successorP"] = "Probability of each sub type to occur";
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


/**
 *  converts string to vector of double see @LeafRandomParameter::readXML
 *	used for leaf shape
 *  @param xmlInput     input array in xml
 *
 * @return xmlInput converted to vector double
 */
std::vector<int> OrganRandomParameter::string2vector(const char* xmlInput, int defaultVal)
{
	std::string buf;                 // Have a buffer string
	std::stringstream ss(xmlInput);       // Insert the string into a stream
	std::vector<int> tokens; // Create vector to hold our words
	while (std::getline(ss, buf, ',')) {
	    tokens.push_back(std::stod(buf));
	}
	return tokens;
}
std::vector<double> OrganRandomParameter::string2vector(const char* xmlInput, double defaultVal)
{
    std::string buf;                 // Have a buffer string
    std::stringstream ss(xmlInput);       // Insert the string into a stream
    std::vector<double> tokens; // Create vector to hold our words
    while (std::getline(ss, buf, ',')) {
        tokens.push_back(std::stod(buf));
    }
    return tokens;
}


std::string OrganRandomParameter::vector2string(std::vector<int> vec) const
{
    std::stringstream ss;
    for (auto it = vec.begin(); it != vec.end(); it++)    {
        if (it != vec.begin()) {
            ss << ", ";
        }
        ss << *it;
    }
	return ss.str() ;
}
std::string OrganRandomParameter::vector2string(std::vector<double> vec) const
{
    std::stringstream ss;
    for (auto it = vec.begin(); it != vec.end(); it++)    {
        if (it != vec.begin()) {
            ss << ", ";
        }
        ss << *it;
    }
    return ss.str() ;
}


/**
 * Choose (dice) lateral type based on stem parameters successor and successorP
 *
 * @param pos       spatial position (for coupling to a soil model)
 * @return          stem sub type of the lateral stem
 */
int OrganRandomParameter::getLateralType(const Vector3d& pos, int ruleId)//
{
	 assert(successorST.at(ruleId).size()==successorP.at(ruleId).size()
        && "StemTypeParameter::getLateralType: Successor sub type and probability vector does not have the same size");
    if (successorP.at(ruleId).size()>0) { // at least 1 successor type
        double d = plant.lock()->rand(); // in [0,1]
        int i=0;
        double p=successorP.at(ruleId).at(i);
        i++;
        while ((p<d) && (i<successorP.at(ruleId).size())) {
            p+=successorP.at(ruleId).at(i);
            i++;
        }
        if (p>=d) { // success
            // std::cout << "lateral type " << successor.at(i-1) << "\n" << std::flush;
            return i-1;//successor.at(i-1);
        } else { // no successors
            // std::cout << "no lateral type " << std::flush;
            return -1;
        }
    } else {
        return -1; // no successors
    }
}

} // namespace
