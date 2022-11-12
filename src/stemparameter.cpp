// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "stemparameter.h"

#include "Organism.h"

#include <cmath>
#include <iostream>
#include <chrono>
#include <assert.h>

namespace CPlantBox {

/**
 * @return Mean maximal stem length of this stem type
 */
double StemSpecificParameter::getK() const {
    double l = std::accumulate(ln.begin(), ln.end(), 0.);
    return l+la+lb;
}

/**
 * @copydoc OrganParameter::toString()
 */
std::string StemSpecificParameter::toString() const
{
    std::stringstream str;
    str << "subType\t" << subType << std::endl;
    str << "lb\t" << lb << std::endl << "la\t" << la << std::endl;
    str << "nob\t" << nob() << std::endl << "r\t" << r << std::endl << "a\t" << a << std::endl;
    str << "theta\t" << theta << std::endl << "rlt\t" << rlt << std::endl;
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
StemRandomParameter::StemRandomParameter(std::shared_ptr<Organism> plant) :OrganRandomParameter(plant)
{
    // base class default values
    name = "undefined";
    organType = Organism::ot_stem;
    subType = -1;
    f_tf = std::make_shared<Tropism>(plant);
    bindParmeters();
}

/**
 * @copydoc OrganTypeParameter::copy()
 */
std::shared_ptr<OrganRandomParameter> StemRandomParameter::copy(std::shared_ptr<Organism> plant)
{
    auto r = std::make_shared<StemRandomParameter>(*this); // copy constructor breaks class introspection
    r->plant = plant;
    r->bindParmeters(); // fix class introspection
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
 * Creates a specific stem from the stem type parameters.
 * @return Specific stem parameters derived from the stem type parameters
 */
std::shared_ptr<OrganSpecificParameter> StemRandomParameter::realize()
{
	bool verbose = false;
	if(verbose){
		std::cout<<"StemRandomParameter::realize() "<<successorST.size()<<" "<<lb<<" "<<la<<" "<<ln<<std::endl;
    }
	auto p = plant.lock();
	double lb_;
    double la_;
    std::vector<double> ln_; // stores the inter-distances
	double res;
	int nob_real = 0;
	bool hasLaterals = (successorST.size()>0);
	if (dx <= dxMin){
		std::cout<<"dx <= dxMin, dxMin set to dx/2"<<std::endl;
		this->dxMin = dx/2;
	}
	if (!hasLaterals) { // no laterals

    	lb_ = 0;
        la_ = std::max(lmax + p->randn()*lmaxs, 0.); // la, and lb is ignored
		res = la_-floor(la_ / dx)*dx;
		if(res < dxMin && res != 0){
			if(res <= dxMin/2){ la_ -= res;
			}else{la_ =  floor(la_ / dx)*dx + dxMin;}
		}			//make la_ compatible with dx() and dxMin()

    } else {
    lb_ = std::max(lb + p->randn()*lbs, 0.); // length of basal zone
	la_ = std::max(la + p->randn()*las, 0.); // length of apical zone
	nob_real = std::max((lmax-la_-lb_)/ln+1, 1.);//std::max(round(nob() + p->randn()*nobs()), 1.); // real maximal number of branches 			  
	res = lb_ - floor(lb_/dx)* dx;	
	if((res < dxMin) && (res != 0)){
		if(res <= dxMin/2){ lb_ -= res;
		}else{lb_ =  floor(lb_ / dx)*dx + dxMin;}
	}	
    res = la_-floor(la_ / dx)*dx;	
	if(res < dxMin && res != 0){
		if(res <= dxMin/2){ la_ -= res;
		}else{la_ =  floor(la_ / dx)*dx + dxMin;}
	}	
	double ln_mean = ln;
	if(ln < dxMin*0.99 && ln !=0){
		std::cout<<"\nStemRandomParameter::realize inter-lateral distance (ln) "<<ln<<" below minimum resolution (dxMin) "<<dxMin<<". ln set to dxMin"<<std::endl;
		ln_mean = dxMin;
	}
	
	//adapt number of laterals by branching point to keep same total number of lats
	//in spite of dxMin
	int nob1 = std::max((lmax-la_-lb_)/ln_mean+1, 1.);//use new la_, lb_ and ln_mean
    int nob_ = std::min(std::max(round(nob1 + p->randn()*nobs()), 1.),double(nob_real));// maximal number of branches
	int latMissing = nob_real - nob_;
	int latExtra1 = floor(latMissing/nob_);//mean number of extra laterals per branching point to keep correct number
	int latExtra2 = latMissing - latExtra1*(nob_);
	int latExtra2_ = latExtra2;
	if(verbose){
		std::cout<<"in create lat branch "<<nob_real<<" ";
		std::cout<<nob1<< " " <<nob_<<" "<<latMissing<<" "<<latExtra1<<" "<<latExtra2<<std::endl;
	}
		//at end of basal zone
		for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
		if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}
		
		switch(lnf) {
		case 0: // homogeneously distributed stem nodes
		for (int i = 0; i<nob_-1; i++) { // create inter-stem distances
			double d = std::max(ln_mean +p->randn()*lns,1.e-5); //Normal function of equal internode distance
			res = d -floor(d / dx)*dx;
			if(res < dxMin && res != 0){
				if(res <= dxMin/2){d -= res;
				}else{d = floor(d / dx)*dx + dxMin;}
				
				} //make ln compatible with dx() and dxMin().
			
			ln_.push_back(d);
			for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
			if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}  


		};break;
		case 1: //  nodes distance increase linearly
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d =  std::max(ln_mean*(1+i) +p->randn()*lns,1.e-5); //std::max(  );//ln +p->randn()*lns,1e-9);
			res = d -floor(d / dx)*dx;
			if(res < dxMin && res != 0){
				if(res <= dxMin/2){d -= res;
				}else{d = floor(d / dx)*dx + dxMin;}
				
				} //make ln compatible with dx() and dxMin().
			
			ln_.push_back(d);
			for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
			if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}

		};break;
		case 2: //nodes distance decrease linearly
		for (int i = 0; i<nob_-1; i++) { // create inter-stem distances
			double d =  std::max(ln_mean*(1+i) +p->randn()*lns,1.e-5); //std::max(  );//ln +p->randn()*lns,1e-9);
			res = d -floor(d / dx)*dx;
			if(res < dxMin && res != 0){
				if(res <= dxMin/2){d -= res;
				}else{d = floor(d / dx)*dx + dxMin;}
				
				} //make ln compatible with dx() and dxMin().
			
			ln_.push_back(d);
			for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
			if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}

		};break;
		case 3: //nodes distance cst
		for (int i = 0; i<nob_-1; i++) { // create inter-stem distances
			double d =  std::max(ln_mean +p->randn()*lns,1.e-5); //std::max(  );//ln +p->randn()*lns,1e-9);
			res = d -floor(d / dx)*dx;
			if(res < dxMin && res != 0){
				if(res <= dxMin/2){d -= res;
				}else{d = floor(d / dx)*dx + dxMin;}
				
				} //make ln compatible with dx() and dxMin().
			
			ln_.push_back(d);
			for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
			if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}

		};break;

		case 4://nodes distance decrease exponential
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d =  std::max(ln_mean/(1+i) +p->randn()*lns,1.e-5); //std::max(  );//ln +p->randn()*lns,1e-9);
			res = d -floor(d / dx)*dx;
			if(res < dxMin && res != 0){
				if(res <= dxMin/2){d -= res;
				}else{d = floor(d / dx)*dx + dxMin;}
				
				} //make ln compatible with dx() and dxMin().
			
			ln_.push_back(d);
			for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
			if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}		
		}; break;
		case 5://nodes distance decrease exponential
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d =  std::max(ln_mean/(1+i) +p->randn()*lns,1.e-5); //std::max(  );//ln +p->randn()*lns,1e-9);
			res = d -floor(d / dx)*dx;
			if(res < dxMin && res != 0){
				if(res <= dxMin/2){d -= res;
				}else{d = floor(d / dx)*dx + dxMin;}
				
				} //make ln compatible with dx() and dxMin().
				
			ln_.push_back(d);
			for (int j = 0; j<latExtra1; j++) { ln_.push_back(0);}
			if (latExtra2_> 0) {ln_.push_back(0);latExtra2_--;}	
		};break;
default:
		throw std::runtime_error("StemRandomParameter::realize type of inter-branching distance not recognized"); 
}}
    double r_ = std::max(r + p->randn()*rs, 0.); // initial elongation
    double a_ = std::max(a + p->randn()*as, 0.); // radius
    double theta_ = std::max(theta + p->randn()*thetas, 0.); // initial elongation
    double rlt_ = std::max(rlt + p->randn()*rlts, 0.); // stem life time
	double delayNGStart_ = std::max(delayNGStart + p->randn()*delayNGStarts, 0.);
	double delayNGEnd_ = std::max(delayNGEnd + p->randn()*delayNGEnds, 0.);
	if(delayNGEnd_ < delayNGStart_){
		std::cout<<"StemRandomParameter::realize() : delayNGEnd_ < delayNGStart_ \n";
		std::cout<<"set delayNGEnd_ = delayNGStart_ = "<<delayNGStart_<<std::endl;
		delayNGEnd_ = delayNGStart_;
	}
	double delayLat_ = std::max(delayLat + p->randn()*delayLats, 0.);
	if(verbose){
		std::cout<<"to std::make_shared<StemSpecificParameter> ";
		std::cout<<subType<<" "<<lb_<<" "<<la_<<" "<<r_<<" "<<hasLaterals<<" "<<this->nodalGrowth;
		std::cout<<" "<<delayNGStart_<<" "<< delayNGEnd_<<" "<< delayLat_<<std::endl;
		std::cout<<"ln: "<<std::endl;
		for (int j = 0; j<ln_.size(); j++) {std::cout<< ln_.at(j)<<" ";}
		std::cout<<std::endl;
	}
    return std::make_shared<StemSpecificParameter>(subType,lb_,la_,ln_,r_,a_,theta_,rlt_,hasLaterals, 
			this->nodalGrowth, delayNGStart_, delayNGEnd_, delayLat_);
}

/**
 * Choose (dice) lateral type based on stem parameters successor and successorP
 *
 * @param pos       spatial position (for coupling to a soil model)
 * @return          stem sub type of the lateral stem
 */
int StemRandomParameter::getLateralType(const Vector3d& pos, int ruleId)//
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

/**
 * todo docme
 *
 * todo I have no idea why this holds...
 */
double StemRandomParameter::nobs() const
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
std::string StemRandomParameter::toString(bool verbose) const {

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
void StemRandomParameter::readXML(tinyxml2::XMLElement* element)
{
    OrganRandomParameter::readXML(element);
    tinyxml2::XMLElement* p = element->FirstChildElement("parameter");
    successorOT.resize(0, std::vector<int>(0));//2D, int
    successorST.resize(0, std::vector<int>(0));//2D, int
    successorP.resize(0, std::vector<double>(0));//2D, double
    successorNo.resize(0);//1D, int
    successorWhere.resize(0, std::vector<int>(0));//2D, int
	double p_ , defaultVald;
	int defaultVal, defaultSize, numLat, success;
	std::vector<std::string> lookfor ;
	bool replaceByDefaultValue;
	int ruleId = 0;
	bool verbose = false;
    while(p) {
        std::string key = p->Attribute("name");
		if(verbose){
			std::cout<<"StemRandomParameter::readXML ";
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
	//p_ == 1 or p_ == successorP.sze()
    //double p_ = std::accumulate(successorP.begin(), successorP.end(), 0.);//need to sum to 1 over the group or over all the successor
    if  ((p_<1)&&(p_>0))  {
        std::cout<<p_ << "StemRandomParameter::readXML: Warning! percentages do not add up to 1. \n";
    }
    assert(successorST.size()==successorP.size() &&
        "StemTypeParameter::readXML: Successor sub type and probability vector does not have the same size" );
	if(verbose){
		std::cout << "successorP "<<successorP.size()<<std::endl;
	for(int k = 0;k < successorP.size(); k++){//<< " "<<successorP.at(k).size()
		for(int o = 0;o < successorP.at(k).size(); o++)
		{
			
			std::cout<<successorP.at(k).at(o)<<" ";
		}std::cout<<std::endl;
	 }std::cout<<std::endl;std::cout<<std::endl;
		std::cout << "successorOT "<<successorOT.size()<<std::endl;
	for(int k = 0;k < successorOT.size(); k++){
		for(int o = 0;o < successorOT.at(k).size(); o++)
		{
			std::cout<<"at k o "<<k<<" "<<o<<" weget ";//<<std::endl;
			std::cout<<successorOT.at(k).at(o)<<" ";
		}std::cout<<std::endl;
	}std::cout<<std::endl;std::cout<<std::endl;
		std::cout << "successorST "<<successorST.size()<<std::endl;
	for(int k = 0;k < successorST.size(); k++){
		for(int o = 0;o < successorST.at(k).size(); o++)
		{
			std::cout<<successorST.at(k).at(o)<<" ";
		}std::cout<<std::endl;
	}std::cout<<std::endl;std::cout<<std::endl;
		std::cout << "successorWhere "<<successorWhere.size()<<std::endl;
	for(int k = 0;k < successorWhere.size(); k++){
		for(int o = 0;o < successorWhere.at(k).size(); o++)
		{
			std::cout<<successorWhere.at(k).at(o)<<" ";
		}std::cout<<std::endl;
	}std::cout<<std::endl;std::cout<<std::endl;
	}
}

/**
 * @copydoc OrganTypeParameter::writeXML()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 */
tinyxml2::XMLElement* StemRandomParameter::writeXML(tinyxml2::XMLDocument& doc, bool comments) const
{
    assert(successorST.size()==successorP.size() &&
        "StemTypeParameter::writeXML: Successor sub type and probability vector does not have the same size" );
    tinyxml2::XMLElement* element = OrganRandomParameter::writeXML(doc, comments);
	double p_ = 0.;
    for (int j = 0; j<successorST.size(); j++) {
			std::cout<<"StemRandomParameter::writeXML "<<j<<" "<<successorNo.size()<<" "
			<<successorWhere.size()<<" "<<successorST.size()<<" "<<successorP.size()
			<< " "<<successorOT.size()<<std::endl;
			tinyxml2::XMLElement* p = doc.NewElement("parameter");
			p->SetAttribute("name", "successor");
			std::cout<<"StemRandomParameter::writeXMLA"<<std::endl;
			if(successorNo.size()>j){p->SetAttribute("numLat", successorNo.at(j));}
			std::cout<<"StemRandomParameter::writeXMLB"<<std::endl;
			if(successorWhere.size()>j){p->SetAttribute("Where", vector2string(successorWhere.at(j)).c_str());}
			std::cout<<"StemRandomParameter::writeXMLC"<<std::endl;
			p->SetAttribute("subType", vector2string(successorST.at(j)).c_str());
			std::cout<<"StemRandomParameter::writeXMLD"<<std::endl;
			if(successorOT.size()>j){p->SetAttribute("organType", vector2string(successorOT.at(j)).c_str());}
			std::cout<<"StemRandomParameter::writeXMLE"<<std::endl;
			p->SetAttribute("percentage", vector2string(successorP.at(j)).c_str());
			std::cout<<"StemRandomParameter::writeXMLF"<<std::endl;
			element->InsertEndChild(p);
			std::cout<<"StemRandomParameter::writeXMLG"<<std::endl;
			if (comments) {
				std::cout<<"StemRandomParameter::writeXMLH"<<std::endl;
				std::string str = description.at("successorST");
			std::cout<<"StemRandomParameter::writeXMLI"<<std::endl;
				tinyxml2::XMLComment* c = doc.NewComment(str.c_str());
			std::cout<<"StemRandomParameter::writeXMLJ"<<std::endl;
				element->InsertEndChild(c);
			std::cout<<"StemRandomParameter::writeXMLK"<<std::endl;
			}
			std::cout<<"StemRandomParameter::writeXMLL"<<std::endl;
			p_ += std::accumulate(successorP.at(j).begin(), successorP.at(j).end(), 0.);
			std::cout<<"StemRandomParameter::writeXMLM"<<std::endl;
			if ((p_<1) && (p_!=0)) {
				std::cout << "LeafRandomParameter::writeXML: Warning! percentages do not add up to 1. = " << p_ << "\n";
			std::cout<<"StemRandomParameter::writeXMLN"<<std::endl;
			}
			std::cout<<"StemRandomParameter::writeXMLO"<<std::endl;

		
	}
			std::cout<<"StemRandomParameter::writeXMLP"<<std::endl;
    return element;
}

/**
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void StemRandomParameter::bindParmeters()
{
    OrganRandomParameter::bindParameters();
    bindParameter("lb", &lb, "Basal zone [cm]", &lbs);
    bindParameter("la", &la, "Apical zone [cm]", &las);
    bindParameter("ln", &ln, "Inter-lateral distance [cm]", &lns);
    bindParameter("lnf", &lnf, "Type of inter-branching distance (0 homogeneous, 1 linear inc, 2 linear dec, 3 exp inc, 4 exp dec)");
    bindParameter("lmax", &lmax, "Maximal stem length [cm]", &lmaxs);
    bindParameter("r", &r, "Initial growth rate [cm day-1]", &rs);
    bindParameter("a", &a, "Stem radius [cm]", &as);
    bindParameter("RotBeta", &rotBeta, "RevRotation of the stem");  /// todo improve description, start lower letter
    bindParameter("BetaDev", &betaDev, "RevRotation deviation");  /// todo improve description, start lower letter
    bindParameter("InitBeta", &initBeta, "Initial RevRotation");  /// todo improve description, start lower letter
    bindParameter("tropismT", &tropismT, "Type of stem tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)");
    bindParameter("tropismN", &tropismN, "Number of trials of stem tropism");
    bindParameter("tropismS", &tropismS, "Mean value of expected change of stem tropism [1/cm]");
	bindParameter("tropismAge", &tropismAge, "Age at which organ switch tropism", &tropismAges);
    bindParameter("theta", &theta, "Angle between stem and parent stem [rad]", &thetas);
    bindParameter("rlt", &rlt, "Stem life time [day]", &rlts);
    bindParameter("gf", &gf, "Growth function number [1]", &rlts);
	bindParameter("nodalGrowth", &nodalGrowth, "nodal growth function (sequential = 0, equal = 0)");
    bindParameter("delayNGStart", &delayNGStart, "delay between stem creation and start of nodal growth", &delayNGStarts);
    bindParameter("delayNGEnd", &delayNGEnd, "delay between stem creation and start of nodal growth", &delayNGEnds);
    bindParameter("delayLat", &delayLat, "delay between latteral creation and start of nodal growth", &delayLats);
     // other parameters (descriptions only)
    description["successorST"] = "Sub type of lateral stems";
    description["successorP"] = "Probability of each sub type to occur";
}

} // end namespace CPlantBox
