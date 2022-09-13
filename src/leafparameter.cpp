// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "leafparameter.h"

#include "Organism.h"

#include <cmath>
#include <iostream>
#include <chrono>
#include <assert.h>
#include <algorithm>

namespace CPlantBox {

/**
 * @return Mean maximal leaf length of this leaf type
 */
double LeafSpecificParameter::getK() const {
	double l = std::accumulate(ln.begin(), ln.end(), 0.); // branchning zone
	return l+la+lb;
}

/**
 * @copydoc OrganParameter::toString()
 */
std::string LeafSpecificParameter::toString() const
{
	std::stringstream str;
	str << "Leaf" << std::endl;
	str << "subType\t" << subType << std::endl;
	str << "lb\t" << lb << std::endl << "la\t" << la << std::endl;
	str << "r\t" << r << std::endl << "a\t" << a << std::endl;
	str << "theta\t" << theta << std::endl << "rlt\t" << rlt << std::endl;
	str << "areaMax\t" << areaMax << std::endl << "laterals\t" << laterals << std::endl;
	str << "ln\t";
	for (int i=0; i<ln.size(); i++) {
		str << ln[i] << " ";
	}
	str << std::endl;
	return str.str();
}



/**
 * Default constructor sets up hashmaps for class introspection
 */
LeafRandomParameter::LeafRandomParameter(std::shared_ptr<Organism> plant) :OrganRandomParameter(plant)
{
	// base class default values
	name = "undefined";
	organType = Organism::ot_leaf;
	subType = -1;
	f_tf = std::make_shared<Tropism>(plant);
	bindParameters();
}

/**
 * @copydoc OrganTypeParameter::copy()
 */
std::shared_ptr<OrganRandomParameter> LeafRandomParameter::copy(std::shared_ptr<Organism> plant)
{
	// std::cout << "LeafRandomParameter::copy\n"<< std::flush;
	auto r = std::make_shared<LeafRandomParameter>(*this); // copy constructor breaks class introspection
	r->plant = plant;
	r->bindParameters(); // fix class introspection
	r->f_tf = f_tf->copy(plant); // copy call back classes
	r->f_gf = f_gf->copy();
	r->f_se = f_se->copy();
	r->f_sa = f_sa->copy();
	r->f_sbp = f_sbp->copy();
	return r;
}

/**
 * @copydoc OrganTypeParameter::realize()
 *
 * Creates a specific leaf from the leaf type parameters.
 * @return Specific leaf parameters derived from the leaf type parameters
 */
std::shared_ptr<OrganSpecificParameter> LeafRandomParameter::realize()
{
	bool verbose = true;
	bool hasLaterals = (successorST.size()>0);
	auto p = plant.lock();
	//define the parameters outside fo the if functions:
	double lb_;
	double la_;
	std::vector<double> ln_; // stores the inter-distances
	double res;
	int nob_real = 0; //number of branching nodes
	if((shapeType == shape_2D)&&(leafGeometryPhi.size()==0))//because default value of shapeType is shape_2D
	{
		std::cout<<"shapeType set to shape_2D but no 2D shape data given, set shapeType to cylinder"<<std::endl;
		shapeType = shape_cylinder;
	}
	if (dx <= dxMin){
		std::cout<<"dx <= dxMin, dxMin set to dx/2"<<std::endl;
		this->dxMin = dx/2;
	}
	if (!hasLaterals) { // no laterals

    	//length of sheath or petiole => basal zone has a meaning even without lateral
		lb_ = std::max(lb + p->randn()*lbs,double(0)); 
		double lmax_ = std::max(lmax + p->randn()*lmaxs, 0.) ; // adapt total length
		res = lmax_-floor(lmax_ / dx)*dx;
		if(res < dxMin && res != 0){
			if(res <= dxMin/2){ lmax_ -= res;
			}else{lmax_ =  floor(lmax_ / dx)*dx + dxMin;}   
		}			//make la_ compatible with dx() and dxMin()
		la_ = lmax_ - lb_; 
    } else {
	lb_ = std::max(lb + p->randn()*lbs,double(0)); // length of basal zone
	la_ = std::max(la + p->randn()*las,double(0)); // length of apical zone
	nob_real = std::max((lmax-la_-lb_)/ln+1, 1.);// real maximal number of branches 
	res = lb_ - floor(lb_/dx)* dx;	
	if (res < dxMin) {
		if (res <= dxMin/2){
			lb_ -= res;
		} else {
			lb_ =  floor(lb_ / dx)*dx + dxMin;
		}
	}	

	res = la_-floor(la_ / dx)*dx;	
	if(res < dxMin && res != 0){
		if(res <= dxMin/2){ la_ -= res;
		}else{la_ =  floor(la_ / dx)*dx + dxMin;}
	}
	double ln_mean = ln;
	if (ln < dxMin*0.99 && ln !=0){
		std::cout<<"\nLeafRandomParameter::realize inter-lateral distance (ln) "<<ln<<" below minimum resolution (dxMin) "<<dxMin<<". ln set to dxMin"<<std::endl;
		ln_mean = dxMin;
	}

	int nob1 = std::max((lmax-la_-lb_)/ln_mean+1, 1.);//use new la_, lb_ and ln_mean
    int nob_ = std::min(std::max(round(nob1 + p->randn()*nobs()), 1.),double(nob_real)); // maximal number of branches +1
	int latMissing = nob_real - nob_;
	int latExtra1 = floor(latMissing/nob_);//mean number of extra laterals per branching point to keep correct number
	int latExtra2 = latMissing - latExtra1*(nob_);
	int latExtra2_ = latExtra2;
	//at end of basal zone
	for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
	if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}
	if(verbose)
	{
	std::cout<<"leaf=> lnf is "<<lnf<<" "<<subType<<std::endl;
	std::cout<<"leaf=> in create lat branch "<<nob1<< " " <<nob_<<" "<<latMissing<<" "<<latExtra1<<" "<<latExtra2<<std::endl;
	}
	switch(lnf) {
	case 0: // homogeneously distributed leaf nodes
		for (int i = 0; i<nob_-1; i++) { // create inter-leaf distances
			double d = std::max(ln_mean + p->randn()*lns,dxMin); //Normal function of equal internode distance
			res = d -floor(d / dx)*dx;
			if(res < dxMin && res != 0){
				if(res <= dxMin/2){d -= res;
				}else{d = floor(d / dx)*dx + dxMin;}

			} //make ln compatible with dx() and dxMin().

			ln_.push_back(d);
			for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
			if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}
		}
		break;
	case 1: //nodes distance increase linearly TODO
		for (int i = 0; i<nob_*2-1; i++) { // create inter-leaf distances
			double d =  std::max(ln_mean*(1+i) + p->randn()*lns,dxMin); //std::max(  );//ln + randn()*lns,1e-9);
			res = d -floor(d / dx)*dx;
			if(res < dxMin && res != 0){
				if(res <= dxMin/2){d -= res;
				}else{d = floor(d / dx)*dx + dxMin;}

			} //make ln compatible with dx() and dxMin().

			ln_.push_back(d);
			for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
			if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}	
		}
		break;
	case 2: //nodes distance decrease linearly TODO
		for (int i = 0; i<nob_-1; i++) { // create inter-leaf distances
			double d =  std::max(ln_mean*(1+i) + p->randn()*lns,dxMin); //std::max(  );//ln + randn()*lns,1e-9);
			res = d -floor(d / dx)*dx;
			if(res < dxMin && res != 0){
				if(res <= dxMin/2){d -= res;
				}else{d = floor(d / dx)*dx + dxMin;}

			} //make ln compatible with dx() and dxMin().

			ln_.push_back(d);
			for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
			if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}	
		}
		break;
	case 3: //nodes distance increase exponential TODO
		for (int i = 0; i<nob_-1; i++) { // create inter-leaf distances
			double d =  std::max(ln_mean + p->randn()*lns,dxMin); //std::max(  );//ln + randn()*lns,1e-9);
			res = d -floor(d / dx)*dx;
			if(res < dxMin && res != 0){
				if(res <= dxMin/2){d -= res;
				}else{d = floor(d / dx)*dx + dxMin;}

			} //make ln compatible with dx() and dxMin().

			ln_.push_back(d);
			for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
			if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}	
			// latExtra2_-=0.5;}			   
			// ln_.push_back(0);
			// for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
			// if (latExtra2_> 0) {ln_.push_back(0);latExtra2_-=0.5;}
		}
		break;
	case 4://nodes distance decrease exponential TODO
		for (int i = 0; i<nob_*2-1; i++) { // create inter-leaf distances
			double d =  std::max(ln_mean/(1+i) + p->randn()*lns,dxMin); //std::max(  );//ln + randn()*lns,1e-9);
			res = d -floor(d / dx)*dx;
			if(res < dxMin && res != 0){
				if(res <= dxMin/2){d -= res;
				}else{d = floor(d / dx)*dx + dxMin;}

			} //make ln compatible with dx() and dxMin().

			ln_.push_back(d);
			for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
			if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}		
		}
		break;
	case 5://nodes distance decrease exponential
		for (int i = 0; i<nob_*2-1; i++) { // create inter-leaf distances
			double d =  std::max(ln_mean/(1+i) + p->randn()*lns,dxMin); //std::max(  );//ln + randn()*lns,1e-9);
			res = d -floor(d / dx)*dx;
			if(res < dxMin && res != 0){
				if(res <= dxMin/2){d -= res;
				}else{d = floor(d / dx)*dx + dxMin;}

			} //make ln compatible with dx() and dxMin().

			ln_.push_back(d);
					for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
					if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;} 
																											 
   
		}
		break;
	default:
		throw 1; // TODO make a nice one
	}}
	double r_ = std::max(r + p->randn()*rs, 0.); // initial elongation
	double  Width_blade_ = std::max( Width_blade + p->randn()* Width_blades, 0.); // initial elongation
	double  Width_petiole_ = std::max( Width_petiole + p->randn()* Width_petioles, 0.); // initial elongation
	double a_ = std::max(a + p->randn()*as, 0.); // radius
	double theta_ = std::max(theta + p->randn()*thetas, 0.); // initial elongation
	double rlt_ = std::max(rlt + p->randn()*rlts, 0.); // leaf life time
	double leafArea_ = std::max(areaMax + p->randn()*areaMaxs, 0.); // radius
	double delayLat_ = std::max(delayLat + p->randn()*delayLats, 0.);
	if(verbose)
	{
	std::cout<<"to std::make_shared<LeafSpecificParameter> ";
	std::cout<<subType<<" "<<lb_<<" "<<la_<<" "<<r_<<" "<<hasLaterals<<std::endl;
	std::cout<<"ln: "<<std::endl;
	for (int j = 0; j<ln_.size(); j++) {std::cout<< ln_.at(j)<<" ";}
	std::cout<<std::endl;
	}
	return std::make_shared<LeafSpecificParameter>(subType,lb_,la_,ln_,r_,a_,theta_,rlt_,leafArea_, 
			hasLaterals, Width_blade_, Width_petiole_, delayLat_);
}

/**
 * Choose (dice) lateral type based on leaf parameters successor and successorP
 *
 * @param pos       spatial position (for coupling to a soil model)
 * @return          leaf sub type of the lateral leaf
 */
int LeafRandomParameter::getLateralType(const Vector3d& pos, int ruleId)
{
	assert(successorST.at(ruleId).size()==successorP.at(ruleId).size()
        && "LeafTypeParameter::getLateralType: Successor sub type and probability vector does not have the same size");
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
			return i-1;
		} else { // no successors
			// std::cout << "no lateral type " << std::flush;
			return -1;
		}
	} else {
		return -1; // no successors
	}
}

/**
 * todo docme
 *
 * todo I have no idea why this holds...
 */
double LeafRandomParameter::nobs() const
{
	double nobs = (lmaxs/lmax - lns/ln)*lmax/ln; // error propagation
	if (la>0) {
		nobs -= (las/la - lns/ln)*la/ln;
	}
	if (lb>0) {
		nobs -= (lbs/lb - lns/ln)*lb/ln;
	}
	return std::max(nobs,0.);
}

/**
 * @copydoc OrganTypeParameter::toString()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 */
std::string LeafRandomParameter::toString(bool verbose) const {

	if (verbose) {
		std::string s = OrganRandomParameter::toString(true);
		std::stringstream str;
		str << "successor\t";
		for (int i=0; i<successorST.size(); i++) {
            for (int j=0; i<successorST.at(i).size(); j++) {
				str << successorST[i][j] << " ";
			}
			str << "; ";
        }
		str << "\t" << description.at("successor") << std::endl;
		str << "successorP\t";
		for (int i=0; i<successorP.size(); i++) {
            for (int j=0; i<successorP.at(i).size(); j++) {
				str << successorP[i][j] << " ";
			}
			str << "; ";
        }
		str << "\t" << description.at("successorP") << std::endl;
		///
		
		str << "leafGeometryPhi\t";
		for (int i=0; i<leafGeometryPhi.size(); i++) {
			str << leafGeometryPhi[i] << " ";
		}
		str << "\t" << description.at("leafGeometryPhi") << std::endl;
		str << "leafGeometryX\t";
		for (int i=0; i<leafGeometryX.size(); i++) {
			str << leafGeometryX[i] << " ";
		}
		str << "\t" << description.at("leafGeometryX") << std::endl;
		
		return s.insert(s.length()-4, str.str());
	} else {
		return OrganRandomParameter::toString(false);
	}

}

/**
 * @copydoc OrganTypeParameter::readXML()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 *
 * If the parameter successor or successorP are not in the element, they are set to zero size.
 */
void LeafRandomParameter::readXML(tinyxml2::XMLElement* element)
{
	OrganRandomParameter::readXML(element);
	tinyxml2::XMLElement* p = element->FirstChildElement("parameter");
	successorOT.resize(0, std::vector<int>(0));//2D, int
    successorST.resize(0, std::vector<int>(0));//2D, int
    successorP.resize(0, std::vector<double>(0));//2D, double
    successorNo.resize(0);//1D, int
    successorWhere.resize(0, std::vector<int>(0));//2D, int
	double p_ =0., defaultVald;
	int defaultVal, defaultSize, numLat, success;
	std::vector<std::string> lookfor ;
	bool replaceByDefaultValue;
	int ruleId = 0;
	bool verbose = true;
	while(p) {
		std::string key = p->Attribute("name");
		if(verbose){
			std::cout<<"LeafRandomParameter::readXML ";
			std::cout<<key<<" "<<key.compare("successor")<<std::endl;
		}
        if (key.compare("successor")==0)  {
			if(successorNo.size()< (ruleId+1))
			{
				successorOT.push_back(std::vector<int>());
				successorST.push_back(std::vector<int>());
				successorWhere.push_back(std::vector<int>());
				successorP.push_back(std::vector<double>());
				successorNo.push_back(1);
			}
			//default == make one lateral
			success = p->QueryIntAttribute("numLat",&numLat);
			if(success == tinyxml2::XML_SUCCESS){successorNo.at(ruleId) = numLat;}
			
			
			//default == empty vector == apply rule to all the linking nodes
			replaceByDefaultValue = true;lookfor = std::vector<std::string>{"where"};
			defaultVal = -1;defaultSize = int(successorWhere.at(ruleId).size() ==0);
			cpb_queryStringAttribute(lookfor,
						defaultVal,defaultSize, replaceByDefaultValue,
						successorWhere.at(ruleId), p);//name, default value, vector to fill, accept not found
			
			replaceByDefaultValue = false;lookfor = std::vector<std::string>{"subType","type"};
			defaultVal = 1.0;defaultSize = 0;
			cpb_queryStringAttribute(lookfor,
						defaultVal,defaultSize, replaceByDefaultValue,
						successorST.at(ruleId), p);
			
			
			
			replaceByDefaultValue = true;lookfor = std::vector<std::string>{"percentage"};
			defaultVald = 1./successorST.at(ruleId).size();
			defaultSize = (successorST.at(ruleId).size() - successorP.at(ruleId).size());
			
			if(!std::isfinite(defaultVal))
			{
				std::cout<<"!std::isfinite(defaultVal) "<<defaultVald<<" ";
				std::cout<<" did sucST "<< successorST.at(ruleId).size() <<std::endl;
				for(int k = 0;k < successorST.at(ruleId).size(); k++){
					std::cout<<successorST.at(ruleId).at(k)<<" ";}
				std::cout<<std::endl;
				assert(false);
			}
			cpb_queryStringAttribute(lookfor,
						defaultVald,defaultSize, replaceByDefaultValue,
						successorP.at(ruleId), p);
			
			replaceByDefaultValue = true;lookfor = std::vector<std::string>{"organType"};
			if(successorST.at(ruleId).at(0) == 2){
				//for backward compatibility => if not organtype and subtype == 2, then we want a leaf
				defaultVal = Organism::ot_leaf; 
			}else{defaultVal = this->organType;}
			defaultSize = (successorST.at(ruleId).size() - successorOT.at(ruleId).size());
			cpb_queryStringAttribute(lookfor,
						defaultVal,defaultSize, replaceByDefaultValue,
						successorOT.at(ruleId), p);
			
			//sum(p_) <= 1. 
			p_ = std::accumulate(successorP.at(ruleId).begin(), successorP.at(ruleId).end(), 0.);
			if(verbose){
			std::cout<<"sum p_ is "<<p_<<" "<< (p_ == 1.) <<" "<< ruleId <<std::endl;
				std::cout<<" did sucST "<< successorST.size() <<std::endl;
				std::cout<<" did sucST "<< successorST.at(ruleId).size() <<std::endl;
				std::cout<<" did sucP "<< successorP.size() <<std::endl;
				std::cout<<" did sucP "<< successorP.at(ruleId).size() <<std::endl;
				std::cout<<" did sucQT "<< successorOT.size() <<std::endl;
				std::cout<<" did sucQT "<< successorOT.at(ruleId).size() <<std::endl;
				std::cout<<" did sucNo "<< successorNo.size() <<std::endl;
				std::cout<<" did sucWhere "<< successorWhere.size() <<std::endl;
				std::cout<<" did sucWhere "<< successorWhere.at(ruleId).size() <<std::endl;
			}
			if(p_ == 1.) //we gathered on group of successor
			{
				ruleId ++;
				if(verbose){
				std::cout<<"new ruleId "<<ruleId<<std::endl;
				}
				// successorP.push_back(sucP); 	sucP.resize(0);
				// successorOT.push_back(sucOT); 	sucOT.resize(0);
				// successorST.push_back(sucST); 	sucST.resize(0);
				// successorWhere.push_back(sucW); sucW.resize(0);
			}
        }
		p = p->NextSiblingElement("parameter");
	}
	//double p_ = std::accumulate(successorP.begin(), successorP.end(), 0.);
	if  ((p_<1) && (p_!=0))  {
		std::cout << "LeafRandomParameter::readXML: Warning! percentages to not add up to 1. \n";
	}
	//go back to start of the list?
	p = element->FirstChildElement("parameter");						 
	leafGeometryPhi.resize(0);
	leafGeometryX.resize(0);
	while(p) {
		std::string key = p->Attribute("name");
		if (key.compare("leafGeometry")==0)  {
			double dummy = 0.0;
			//leafGeometryPhi.push_back(p->DoubleAttribute("phi"));
			leafGeometryPhi = string2vector(p->Attribute("phi"), dummy);
			//leafGeometryX.push_back(p->DoubleAttribute("x"));
			leafGeometryX = string2vector(p->Attribute("x"), dummy);
		}
		p = p->NextSiblingElement("parameter");
	}
	// create geometry
	if (parametrisationType==0) {
		createLeafRadialGeometry(leafGeometryPhi, leafGeometryX, geometryN);
	} else if (parametrisationType==1) {
		createLeafGeometry(leafGeometryPhi, leafGeometryX, geometryN);
	} else {
		std::cout << "LeafRandomParameter::readXML: Warning! unknown parametrisation type, could not create leaf geometry \n";
	}
}

/**
 * @copydoc OrganTypeParameter::writeXML()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 */
tinyxml2::XMLElement* LeafRandomParameter::writeXML(tinyxml2::XMLDocument& doc, bool comments) const
{
	assert(successorST.size()==successorP.size() &&
			"LeafTypeParameter::writeXML: Successor sub type and probability vector does not have the same size" );
	tinyxml2::XMLElement* element = OrganRandomParameter::writeXML(doc, comments);
	double p_ = 0.;
	
    for (int j = 0; j<successorST.size(); j++) {	
			tinyxml2::XMLElement* p = doc.NewElement("parameter");
			p->SetAttribute("name", "successor");
			if(successorNo.size()>j){p->SetAttribute("numLat", successorNo.at(j));}
			if(successorWhere.size()>j){p->SetAttribute("Where", vector2string(successorWhere.at(j)).c_str());}
			p->SetAttribute("subType", vector2string(successorST.at(j)).c_str());
			if(successorOT.size()>j){p->SetAttribute("organType", vector2string(successorOT.at(j)).c_str());}
			p->SetAttribute("percentage", vector2string(successorP.at(j)).c_str());
			element->InsertEndChild(p);
			if (comments) {
				std::string str = description.at("successor");
				tinyxml2::XMLComment* c = doc.NewComment(str.c_str());
				element->InsertEndChild(c);
			}
			p_ += std::accumulate(successorP.at(j).begin(), successorP.at(j).end(), 0.);
			if ((p_<1) && (p_!=0)) {
				std::cout << "LeafRandomParameter::writeXML: Warning! percentages do not add up to 1. = " << p_ << "\n";
			}

		
	}
	// double p_ = std::accumulate(successorP.begin(), successorP.end(), 0.);
	// if ((p_<1) && (p_!=0)) {
		// std::cout << "LeafRandomParameter::writeXML: Warning! percentages do not add up to 1. = " << p_ << "\n";
	// }
	for (int i = 0; i<leafGeometryPhi.size(); i++) {
		tinyxml2::XMLElement* p = doc.NewElement("parameter");
		p->SetAttribute("name", "leafGeometry");
		p->SetAttribute("phi", leafGeometryPhi[i]);
		p->SetAttribute("x", leafGeometryX[i]);
		element->InsertEndChild(p);
		if (comments) {
			std::string str = description.at("leafGeometry");
			tinyxml2::XMLComment* c = doc.NewComment(str.c_str());
			element->InsertEndChild(c);
		}
	}
	return element;
}

/**
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void LeafRandomParameter::bindParameters()
{
	OrganRandomParameter::bindParameters();
	bindParameter("shapeType", &shapeType, "leaf shape type for definition of area and volume (cylinder = 0, cuboid = 1, user-defined 2D shape= 2)");
	bindParameter("Width_petiole",  &Width_petiole, "width petiole [cm]", &Width_petioles);//when leaf defined as a cuboid
	bindParameter("Width_blade", &Width_blade, "width blade [cm]", &Width_blades);//when leaf defined as a cuboid
	bindParameter("lb", &lb, "Basal zone [cm]", &lbs);
	bindParameter("la", &la, "Apical zone [cm]", &las);
	bindParameter("ln", &ln, "Inter-lateral distance [cm]", &lns);
	bindParameter("lmax", &lmax, "Maximal leaf length [cm]", &lmaxs);
	bindParameter("areaMax", &areaMax, "Maximal leaf area [cm2]", &areaMaxs);
	bindParameter("r", &r, "Initial growth rate [cm day-1]", &rs);
	bindParameter("RotBeta", &rotBeta, "RevRotation of the leaf"); /// todo improve description, start lower letter
	bindParameter("BetaDev", &betaDev, "RevRotation deviation"); /// todo improve description, start lower letter
	bindParameter("InitBeta", &initBeta, "Initial RevRotation"); /// todo improve description, start lower letter
	bindParameter("tropismT", &tropismT, "Type of leaf tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)");
	bindParameter("tropismN", &tropismN, "Number of trials of leaf tropism");
	bindParameter("tropismS", &tropismS, "Mean value of expected change of leaf tropism [1/cm]");
	bindParameter("tropismAge", &tropismAge, "Age at which organ switch tropism", &tropismAges);
	bindParameter("theta", &theta, "Angle between leaf and parent leaf [rad]", &thetas);
	bindParameter("rlt", &rlt, "Leaf life time [day]", &rlts);
	bindParameter("gf", &gf, "Growth function number [1]");
	bindParameter("lnf", &lnf, "Type of inter-branching distance (0 homogeneous, 1 linear inc, 2 linear dec, 3 exp inc, 4 exp dec)");
	bindParameter("parametrisationType", &parametrisationType, "Leaf geometry parametrisation type");
	bindParameter("geometryN", &geometryN, "leaf geometry resolution");
	bindParameter("delayLat", &delayLat, "delay between latteral creation and start of nodal growth", &delayLats);
    // other parameters (descriptions only)
	description["successorST"] = "Sub type of lateral leaf veins";
	description["successorP"] = "Probability of each sub type to occur";
	description["leafGeometryPhi"] = "Leaf geometry parametrisation parameter";
	description["leafGeometryX"] = "Leaf geometry parametrisation";
}

/**
 * resamples incoming leaf geometry (radial star domain) to N y-coordinates axis (2d xy)
 * and normalizes (see normalize())
 */
void LeafRandomParameter::createLeafRadialGeometry(std::vector<double> phi, std::vector<double> l, int N) {
	if (phi.size()>0 && (phi.size()==l.size())) {
		leafGeometry.resize(N);
		auto y_ = Function::linspace(0., leafLength(), N);
		for (int i = 0; i<N; i++) {
			std::vector<double> x = intersections(y_[i], phi, l);
			std::sort(x.begin(), x.end());
			leafGeometry.at(i) = x;
		}
		normalizeLeafNodes();
	} else {
		std::cout << "LeafRandomParameter::createLeafRadialGeometry: Warning! parametrisation vectors y and l are empty or differ in size\n";
	}
}

/**
 * resamples incoming leaf geometry
 */
void LeafRandomParameter::createLeafGeometry(std::vector<double> y, std::vector<double> l, int N) {
	if (y.size()>0 && (y.size()==l.size())) {
		double midy = leafMid();
		leafGeometry.resize(N);
		auto y_ = Function::linspace(-midy, leafLength()-midy, N);
		for (int i = 0; i<N; i++) {
			leafGeometry.at(i).push_back(Function::interp1(y_[i], y, l));
		}
		normalizeLeafNodes();
	} else {
		std::cout << "LeafRandomParameter::createLeafGeometry: Warning! parametrisation vectors y and l are empty or differ in size\n";
	}
}

/**
 *  normalizes points to obtain a leaf area of 1, assuming a length of 1
 */
void LeafRandomParameter::normalizeLeafNodes() {
	double s = 0.;
	for (int i = 0; i< leafGeometry.size(); i++) {
		auto g = leafGeometry.at(i);
		for (auto g_ : g) {
			s += g_;
		}
	}
	s = 2*s/leafGeometry.size(); // leaf area (leafGeometry describes only half a leaf)
	for (int i = 0; i< leafGeometry.size(); i++) {
		auto& g = leafGeometry.at(i);
		for (auto& g_ : g) {
			g_ = g_/s;
		}
	}
}

/**
 * returns the intersection of a horizontal line at y-coordinate with the leaf geometry
 */
std::vector<double> LeafRandomParameter::intersections(double y, std::vector<double> phi, std::vector<double> l) {
	int N = 1000; // high resolution, since we just use nearest neighbor to find the zero
	double midy = leafMid();
	std::vector<double> ix = std::vector<double>();
	bool below = true;
	for (int i = 0; i<N; i++) {
		double p = (double)i/(N-1)*M_PI - M_PI/2.; // [-pi/2, pi/2]
		double l_ = Function::interp1(p, phi, l);
		if (below) {
			if (l_*sin(p)+midy>=y) {
				below = false;
				ix.push_back(l_*cos(p));
			}
		} else {
			if (l_*sin(p)+midy<=y) {
				below = true;
				ix.push_back(l_*cos(p));
			}
		}
	}
	return ix;
}

} // end namespace CPlantBox
