#include "ModelParameter.h"
#include "external/tinyxml2.h"
#include "Plant.h"

#include <cmath>

namespace CPlantBox {




// auxiliary function for xml parsing
void readXMLvs(tinyxml2::XMLElement* el_, std::string name, double* v, double* s)
{
	tinyxml2::XMLElement* el = el_->FirstChildElement(name.c_str());
	if (el!=0) {
		el->QueryDoubleText(v);
		tinyxml2::XMLElement* sd = el->FirstChildElement("sd");
		if (sd!=0) {
			sd->QueryDoubleText(s);
		}
	}
}



OrganRandomOrganParameter::OrganRandomOrganParameter()
{
	organType = Organ::ot_organ;
	subType = -1; // means undefined
}

void OrganRandomOrganParameter::getAttribute(const tinyxml2::XMLElement* param, const char* attr_name, const char* para_name, double &attr, double &deviation )  { //parse through different subtype

        for (; param != 0; param = param->NextSiblingElement(para_name) )
        {   if (param->Attribute("name" , attr_name)) {
                    param->QueryDoubleAttribute("value" , &attr);
                    param->QueryDoubleAttribute("dev" , &deviation);
//                       std::cout << attr_name<<" is " << attr <<"  " << attr_name<<" dev is " << deviation << "\n";
        }
             }
}

void OrganRandomOrganParameter::getAttribute(const tinyxml2::XMLElement* param, const char* attr_name, const char* para_name, double &attr, double &deviation, int &functiontype )  { //parse through different subtype

        for (; param != 0; param = param->NextSiblingElement(para_name) )
        {   if (param->Attribute("name" , attr_name)) {
                    param->QueryDoubleAttribute("value" , &attr);
                    param->QueryDoubleAttribute("dev" , &deviation);
		    param->QueryIntAttribute("functiontype" , &functiontype);
//                       std::cout << attr_name<<" is " << attr <<"  " << attr_name<<" dev is " << deviation << attr_name<<" functiontype is " << functiontype << "\n";
        }
             }
}




void OrganRandomOrganParameter::getAttribute(const tinyxml2::XMLElement* param, const char* attr_name, const char* para_name, int &attr, double &deviation )  { //parse through different subtype

        for (; param != 0; param = param->NextSiblingElement(para_name) )
        {   if (param->Attribute("name" , attr_name)) {
                    param->QueryIntAttribute("value" , &attr);
                    param->QueryDoubleAttribute("dev" , &deviation);
//                       std::cout << attr_name<<" is " << attr <<"  " << attr_name<<" dev is " << deviation << "\n";
        }
             }
}


void OrganRandomOrganParameter::getAttribute(const tinyxml2::XMLElement* param, const char* attr_name, const char* para_name, double &attr)  { //parse through different subtype

        for (; param != 0; param = param->NextSiblingElement(para_name) )
        {   if (param->Attribute("name" , attr_name)) {
                    param->QueryDoubleAttribute("value" , &attr);
//                       std::cout << attr_name<<" is " << attr <<"  " << attr_name<<"\n";
        }
             }
}

void OrganRandomOrganParameter::getAttribute(const tinyxml2::XMLElement* param, const char* attr_name, const char* para_name, int &attr)  { //parse through different subtype

        for (; param != 0; param = param->NextSiblingElement(para_name) )
        {   if (param->Attribute("name" , attr_name)) {
                    param->QueryIntAttribute("value" , &attr);
//                       std::cout << attr_name<<" is " << attr <<"  " << attr_name<<"\n";
        }
             }
}

void OrganRandomOrganParameter::getAttribute(const tinyxml2::XMLElement* param,const char* attr_name, const char* para_name, std::vector<int> &successor, std::vector<double> &successorP )  { //parse through different subtype
        int n = 0;
        for (; param != 0; param = param->NextSiblingElement(para_name) )
        {   if (param->Attribute("name" , attr_name)) {
                    param->QueryIntAttribute("number", &n);
                    successor.push_back(0);
                    successorP.push_back(0);
                    param->QueryIntAttribute("type" , &successor.at(n));
                    param->QueryDoubleAttribute("percentage" , &successorP.at(n));
//                       std::cout << attr_name<<" number " << successor.at(n) <<" percentage " << successorP.at(n)<<"\n";
        }
             }
}



double RootParameter::getK() const
{
	double l = 0;
	for (auto const& dl : ln) {
		l += dl;
	}
	return l+la+lb;
}

void RootParameter::write(std::ostream & cout) const
{
	cout << "# Root Parameters \n";
	cout << "type\t" << subType << "\n"  << "lb\t"<< lb <<"\n" << "la\t"<< la <<"\n" << "ln\t";
	for (size_t i=0; i<ln.size(); i++) {
		cout << ln[i] << "\t";
	}
	cout << "\n" << "nob\t"<< ln.size() <<"\n" << "r\t"<< r <<"\n" << "a\t" << a << "\n" << "theta\t" << theta << "\n" << "rlt\t" << rlt << "\n";
}



RootRandomOrganParameter::RootRandomOrganParameter()
{
	organType = Organ::ot_root;
	subType = -1; // means undefined
	tropism = new TropismFunction(0,0);
	growth = new GrowthFunction();
	se = new SoilLookUp();
	sa = new SoilLookUp();
	sbp = new SoilLookUp();
	set(-1, 0., 0., 10., 0., 1., 0., 0., 0., 1., 0, 0.1, 0., 150./255.,150./255.,50./255., 1, 1. ,0.2, 0.1,
			successor, successorP, 1.22, 0., 1.e9, 0., 1, "undefined");
}

RootRandomOrganParameter::~RootRandomOrganParameter()
{
	delete tropism;
	delete growth;
	delete f_se;
	delete sa;
	delete sbp;
}

/**
 * todo comment
 */
void RootRandomOrganParameter::set(int type, double lb, double lbs, double la, double las, double ln, double lns, double nob, double nobs,
		double r, double rs, double a, double as,  double RotBeta, double BetaDev, double InitBeta, int tropismT, double tropismN, double tropismS,
		double dx, const std::vector<int>& successor, const std::vector<double>& successorP, double theta, double thetas, double rlt, double rlts,
		int gf, const std::string& name)
{
	this->subType = type;
	this->lb = lb;
	this->lbs = lbs;
	this->la = la;
	this->las = las;
	this->ln = ln;
	this->lns = lns;
	this->nob = nob;
	this->nobs = nobs;
	this->r = r;
	this->rs = rs;
	this->a = a;
	this->as = as;
	this->RotBeta = RotBeta;
	this->BetaDev = BetaDev;
	this->InitBeta = InitBeta;
	this->tropismT = tropismT;
	this->tropismN = tropismN;
	this->tropismS = tropismS;
	this->dx=dx;
	this->successor = successor; // vector
	this->successorP = successorP; // vector
	assert(successor.size()==successorP.size());
	this->theta = theta;
	this->thetas = thetas;
	this->rlt = rlt;
	this->rlts = rlts;
	this->gf = gf;
	this->name = name; // string

	createTropism();
	createGrowth();
}

void RootRandomOrganParameter::createTropism(SignedDistanceFunction* geom, SoilLookUp* soil)
{
	//delete tropism;
	TropismFunction* t;
	switch (tropismT) {
	case tt_plagio: {
		t = new Plagiotropism(tropismN,tropismS);
		break;
	}
	case tt_gravi: {
		t = new Gravitropism(tropismN,tropismS);
		break;
	}
	case tt_exo: {
		t = new Exotropism(tropismN,tropismS);
		break;
	}
	case tt_hydro: {
		TropismFunction* gt =  new Gravitropism(tropismN,tropismS);
		TropismFunction* ht = new Hydrotropism(tropismN,tropismS,soil);
		t = new CombinedTropism(tropismN,tropismS,ht,10.,gt,1.); // does only use the objective functions from gravitropism and hydrotropism
		break;
	}
	default: throw std::invalid_argument( "RootSystem::createTropismFunction() tropism type not implemented" );
	}

	if (geom!=nullptr) {
		tropism = new ConfinedTropism(t,geom);
	} else {
		tropism = t;
	}

}

void RootRandomOrganParameter::createGrowth()
{
	switch (gf) {
	case gft_negexp: {
		growth = new ExponentialGrowth();
		break;
	}
	case gft_linear: {
		growth =  new LinearGrowth();
		break;
	}
	default: throw std::invalid_argument( "RootSystem::createGrowthFunction() growth function type not implemented" );
	}
}

/**
 * Creates a specific root from the root type parameters.
 * The unique root id is not set, but must be set from outside.
 * (called by Root::Root())
 *
 * minimal ln distance is 1.e-9
 *
 * \return Specific root parameters derived from the root type parameters
 */
SpecificOrganParamter* RootRandomOrganParameter::realize() const
{
	// type does not change
	double lb_ = std::max(lb + randn()*lbs,double(0)); // length of basal zone
	double la_ = std::max(la + randn()*las,double(0)); // length of apical zone
	std::vector<double> ln_; // stores the inter-distances
	int nob_ = std::max(round(nob + randn()*nobs),double(0)); // maximal number of branches
	for (int i = 0; i<nob_-1; i++) { // create inter-root distances
		double d = std::max(ln + randn()*lns,1e-9);
		ln_.push_back(d);
	}
	double r_ = std::max(r + randn()*rs,double(0)); // initial elongation
	double a_ = std::max(a + randn()*as,double(0)); // radius
	double theta_ = std::max(theta + randn()*thetas,double(0)); // initial elongation
	double rlt_ = std::max(rlt + randn()*rlts,double(0)); // root life time
	RootParameter* p =  new RootParameter(subType,lb_,la_,ln_,r_,a_,theta_,rlt_);
	return p;
}

/**
 * Choose (dice) lateral type based on root parameter set,
 * (since there can be more than one lateral type)
 *
 * @param pos       spatial position (for coupling to a soil model)
 */
int RootRandomOrganParameter::getLateralType(const Vector3d& pos)
{
	assert(successor.size()==successorP.size());
	double scale = sbp->getValue(pos);  //the current model makes not a lot of sense, we may come up with something more clever
	if (successorP.size()>0) { // at least 1 successor typese
		if (successorP.size()>1) { // if there are more than one lateral we have to dice
			double d = rand();
			int i=0;
			double p=successorP.at(i)*scale;
			while (p<=d) {
				i++;
				p+=successorP.at(i)*scale;
			}
			return successor.at(i);
		} else {
			return successor.front();
		}
	} else { // no successors
		return -1;
	}
}

void RootRandomOrganParameter::write(std::ostream & os) const
{
	os << "# Root type parameter for " << name << "\n";
	os << "type\t" << subType << "\n" << "name\t" << name << "\n" << "lb\t"<< lb <<"\t"<< lbs << "\n" << "la\t"<< la <<"\t"<< las << "\n"
			<< "ln\t" << ln << "\t" << lns << "\n" << "nob\t"<< nob <<"\t"<< nobs << "\n" << "r\t"<< r <<"\t"<< rs << "\n" <<
			"a\t" << a << "\t" << as << "\n" << "color\t"<< RotBeta <<"\t"<< BetaDev << "\t" << InitBeta << "\n"
			<< "tropism\t"<< tropismT <<"\t"<< tropismN << "\t" << tropismS << "\n" << "dx\t" << dx << "\n" << "successor\t" << successor.size() << "\t";
	for (size_t i=0; i<successor.size(); i++) {
		os << successor[i] << "\t";
	}
	os << "\n" << "successorP\t" << successorP.size() <<"\t";
	for (size_t i=0; i<successorP.size(); i++) {
		os << successorP[i] << "\t";
	}
	os << "\n" << "theta\t" << theta << "\t" << thetas << "\n" << "rlt\t" << rlt << "\t" << rlts << "\n" << "gf\t" << gf << "\n";
	std::cout << "RootRandomOrganParameter::write is deprecated, use RootRandomOrganParameter::writeXML instead\n";
}



void RootRandomOrganParameter::read(std::istream & is)
{
	char ch[256]; // dummy
	is.getline(ch,256);
	std::string s; // dummy

	is >> s >> subType >> s >> name >> s >> lb >> lbs >> s >> la >> las >> s >> ln >> lns >> s >> k >> ks;
	is >> s >> r >> rs >> s >> a >> as >> s >> RotBeta >> BetaDev >> InitBeta >> s >> tropismT >> tropismN >> tropismS >> s >> dx;
	if (ln > 0) {
		nob=  (k-la-lb)/ln+1;   //conversion, because the input file delivers the lmax value and not the nob value
		nob = std::max(nob,0.);
		nobs = (ks/k - lns/ln)*k/ln; // error propagation
		if (la>0) {
			nobs -= (las/la - lns/ln)*la/ln;
		}
		if (lb>0) {
			nobs -= (lbs/lb - lns/ln)*lb/ln;
		}
		nobs = std::max(nobs,0.);
		if (std::isnan(nobs)) {
			std::cout << "RootRandomOrganParameter::read() nobs is nan \n";
			nobs =0;
		}
	} else {
		nob=0;
		nobs = 0;
	}
	int n;
	is  >> s >> n;
	successor.clear();
	int is_;
	for (int i=0; i<n; i++) {
		is >> is_;
		successor.push_back(is_);
	}
	is >> s >> n;
	successorP.clear();
	double ds;
	for (int i=0; i<n; i++) {
		is >> ds;
		successorP.push_back(ds);
	}
	is >> s >> theta >> thetas >> s >> rlt >> rlts >> s >> gf >> s;
	createTropism();
	createGrowth();
	std::cout << "RootRandomOrganParameter::read is deprecated, use RootRandomOrganParameter::readXML instead\n";
}

void RootRandomOrganParameter::readXML(const tinyxml2::XMLElement* ele) //read subtype parameter from different organ type, used by Plant::openXML
{
   const tinyxml2::XMLElement* ele_param = ele->FirstChildElement("parameter"); //XML elements for parameters
   const char* name; //name of the suborgan,
   ele->QueryUnsignedAttribute("subType", &subType);
   ele->QueryStringAttribute("name", &name);
   getAttribute(ele_param, "lb", "parameter", lb, lbs);
   getAttribute(ele_param, "la", "parameter", la, las);
   getAttribute(ele_param, "ln", "parameter", ln, lns);
   getAttribute(ele_param, "lmax", "parameter", k, ks);
   getAttribute(ele_param, "nob", "parameter", nob, nobs);
   	if (ln > 0) {
		nob=  (k-la-lb)/ln+1;   //conversion, because the input file delivers the lmax value and not the nob value
		nob = std::max(nob,0.);
		nobs = (ks/k - lns/ln)*k/ln; // error propagation
		if (la>0) {
			nobs -= (las/la - lns/ln)*la/ln;
		}
		if (lb>0) {
			nobs -= (lbs/lb - lns/ln)*lb/ln;
		}
		nobs = std::max(nobs,0.);
		if (std::isnan(nobs)) {
			std::cout << "RootRandomOrganParameter::read() nobs is nan \n";
			nobs =0;
		}
	} else {
		nob=0;
		nobs = 0;
	}
   getAttribute(ele_param, "r", "parameter", r, rs);
   getAttribute(ele_param, "a", "parameter", a, as);
   getAttribute(ele_param, "RotBeta", "parameter", RotBeta);
   getAttribute(ele_param, "BetaDev", "parameter", BetaDev);
   getAttribute(ele_param, "InitBeta", "parameter", InitBeta);
   getAttribute(ele_param, "tropismN", "parameter", tropismN);
   getAttribute(ele_param, "tropismT", "parameter", tropismT);
   getAttribute(ele_param, "tropismS", "parameter", tropismS);
      getAttribute(ele_param, "dx", "parameter", dx);
   getAttribute(ele_param, "theta", "parameter", theta, thetas);
   getAttribute(ele_param, "rlt", "parameter", rlt, rlts);
   getAttribute(ele_param, "gf", "parameter", gf);
   getAttribute(ele_param, "successor", "parameter", successor, successorP);
   	createTropism();
	createGrowth();
   std::cout<<"subType "<<subType<<"\n";


}


std::string RootRandomOrganParameter::writeXML(FILE* fp) const
{
tinyxml2::XMLPrinter printer( fp, false, 0 ); //compact mode false, and 0 indent
        printer.OpenElement("organ", false); //compact mode false
        printer.PushAttribute("type","root");

	    switch (subType) {
	case 1 :  printer.PushAttribute("name","taproot"); printer.PushAttribute("subType",subType);// See

	break;
    case 2 :  printer.PushAttribute("name","lateral1"); printer.PushAttribute("subType",subType); // See
	break;
    case 3 :  printer.PushAttribute("name","lateral2"); printer.PushAttribute("subType",subType);// See
	break;
	case 4 :  printer.PushAttribute("name","nodal_root"); printer.PushAttribute("subType",subType); // See
	break;
    case 5 :  printer.PushAttribute("name","shoot_borne_root"); printer.PushAttribute("subType",subType); // See
	break;
	case 6:  printer.PushAttribute("name", "lateral3"); printer.PushAttribute("subType", subType); // See
		break;
    }

    printer.OpenElement("parameter");

        printer.PushAttribute("name","lb"); printer.PushAttribute("value",lb); printer.PushAttribute("dev",lbs);  printer.CloseElement(); printer.PushComment("Basal zone [cm]");	//		 	///< Basal zone [cm]
	printer.OpenElement("parameter");

        printer.PushAttribute("name","la"); printer.PushAttribute("value",la); printer.PushAttribute("dev",las);  printer.CloseElement(); printer.PushComment("Apical zone [cm];");	///< Apical zone [cm];
	printer.OpenElement("parameter");

        printer.PushAttribute("name","ln"); printer.PushAttribute("value",ln); printer.PushAttribute("dev",lns); printer.CloseElement();printer.PushComment("Inter-lateral distance [cm];");		///< Inter-lateral distance [cm];
	printer.OpenElement("parameter");//	printer.PushComment("Number of branches [1];");

	printer.PushAttribute("name","lmax"); printer.PushAttribute("value",k); printer.PushAttribute("dev",ks);  printer.CloseElement();printer.PushComment("Inter-lateral distance [cm];");		///< Inter-lateral distance [cm];
	printer.OpenElement("parameter");//	printer.PushComment("Number of branches [1];");

        printer.PushAttribute("name","nob"); printer.PushAttribute("value",nob); printer.PushAttribute("dev",nobs);  printer.CloseElement();		///< Standard deviation apical zone [cm];
	printer.OpenElement("parameter");//	printer.PushComment("Initial growth rate [cm day-1]");
        printer.PushAttribute("name","r"); printer.PushAttribute("value",r); printer.PushAttribute("dev",rs); printer.CloseElement();		///< Initial growth rate [cm day-1]
	printer.OpenElement("parameter");//	printer.PushComment("Root radius [cm]");
        printer.PushAttribute("name","a"); printer.PushAttribute("value",a); printer.PushAttribute("dev",as);   printer.CloseElement();		///< Root radius [cm]

	printer.OpenElement("parameter");//
        printer.PushAttribute("name","RotBeta");  printer.PushAttribute("value",RotBeta); printer.CloseElement();
    printer.OpenElement("parameter");
        printer.PushAttribute("name","BetaDev");  printer.PushAttribute("value",BetaDev); printer.CloseElement();
    printer.OpenElement("parameter");
        printer.PushAttribute("name","InitBeta");  printer.PushAttribute("value",InitBeta); printer.CloseElement();	///< Root color (red)

	printer.OpenElement("parameter");//	printer.PushComment("Root tropism parameter (Type, number of trials, mean vale of expected change)");
        printer.PushAttribute("name","tropismT");  printer.PushAttribute("value",tropismT);  printer.CloseElement();
    printer.OpenElement("parameter");
        printer.PushAttribute("name","tropismN");  printer.PushAttribute("value",tropismN);  printer.CloseElement();
    printer.OpenElement("parameter");
        printer.PushAttribute("name","tropismS");  printer.PushAttribute("value",tropismS);  printer.CloseElement();	///< Root tropism parameter (Type)

	printer.OpenElement("parameter");//	printer.PushComment("Maximal segment size");
        printer.PushAttribute("name","dx"); printer.PushAttribute("value",dx);  printer.CloseElement();		///< Maximal segment size [cm]
	printer.OpenElement("parameter");//	printer.PushComment("Angle between root and parent root");
        printer.PushAttribute("name","theta"); printer.PushAttribute("value",theta); printer.PushAttribute("dev",thetas); printer.CloseElement();	///< Angle between root and parent root (rad)
	printer.OpenElement("parameter");//	printer.PushComment("Root life time");
        printer.PushAttribute("name","rlt"); printer.PushAttribute("value",rlt); printer.PushAttribute("dev",rlts);  printer.CloseElement();	///< Root life time (days)
	printer.OpenElement("parameter");//	printer.PushComment("Growth function");
        printer.PushAttribute("name","gf"); printer.PushAttribute("value",gf); printer.CloseElement();		///< Growth function (1=negative exponential, 2=linear)


	for (int successorCount = 0; successorCount < successor.size(); successorCount ++ ){
            std::string s = std::to_string(successorCount);
            char const *schar = s.c_str();
        printer.OpenElement("parameter"); printer.PushAttribute("name","successor"); printer.PushAttribute("number",schar); printer.PushAttribute("type",successor[successorCount]); printer.PushAttribute("percentage",successorP[successorCount]);    printer.CloseElement();
//        printer.OpenElement("parameter"); printer.PushAttribute("name","successorP"); printer.PushAttribute("number",schar); printer.PushAttribute("percentage",successorP[successorCount]);   printer.CloseElement();
//
    }
//    printer.PrintSpace();
    printer.CloseElement(false); // close element compact mode false

	return std::string(printer.CStr());


}



SeedRandomOrganParameter::SeedRandomOrganParameter() {
	organType = Organ::ot_seed;
	subType = 0; // there is currently only one type
}

void SeedRandomOrganParameter::read(std::istream & is) {

	std::string s; // dummy
	is  >>  s >> plantingdepth;
	is >> s >> firstB >> s >> delayB >> s >> maxB >> s >> nC >> s >> firstSB >> s >> delaySB >> s >> delayRC >> s >> nz >> s >> simtime;
	seedPos = Vector3d(0,0,-plantingdepth);
	std::cout << "SeedTypeParamter::read is deprecated, use SeedTypeParamter::readXML instead\n";
}


void SeedRandomOrganParameter::write(std::ostream & cout) const {
	double pd = -seedPos.z;
	cout <<  "plantingdepth\t" << pd << "\n" <<  "firstB\t" << firstB << "\n" <<  "delayB\t" << delayB << "\n"
			<<  "maxB\t" << maxB << "\n" <<  "nC\t" << 0 << "\n" <<  "firstSB\t" << 0 << "\n"
			<<  "delaySB\t" << 0 << "\n" <<  "delayRC\t" << 0 << "\n" <<  "nz\t" << 0 << "\n" << "simulationTime\t" << 0 << "\n";
	std::cout << "SeedTypeParamter::write is deprecated, use SeedTypeParamter::writeXML instead\n";
}

void SeedRandomOrganParameter::readXML(const tinyxml2::XMLElement* ele) //read subtype parameter from different organ type, used by Plant::openXML
{
   const tinyxml2::XMLElement* ele_param = ele->FirstChildElement("parameter"); //XML elements for parameters
   const char* name;
   ele->QueryUnsignedAttribute("subType", &subType);
   ele->QueryStringAttribute("name", &name);

   getAttribute(ele_param, "seedPos.x", "parameter", seedPos.x);
   getAttribute(ele_param, "seedPos.y", "parameter", seedPos.y);
   getAttribute(ele_param, "seedPos.z", "parameter", seedPos.z);
   getAttribute(ele_param, "plantingdepth", "parameter", plantingdepth);
   getAttribute(ele_param, "firstB", "parameter", firstB, firstBs);
   getAttribute(ele_param, "firstSB", "parameter", firstSB, firstSBs);
   getAttribute(ele_param, "delayB", "parameter", delayB, delayBs);
   getAttribute(ele_param, "delaySB", "parameter", delaySB, delaySBs);
   getAttribute(ele_param, "delayRC", "parameter", delayRC, delayRCs);
   getAttribute(ele_param, "nC", "parameter", nC);
   getAttribute(ele_param, "nz", "parameter", nz);
   getAttribute(ele_param, "maxB", "parameter", maxB, maxBs);
    getAttribute(ele_param, "maxTi", "parameter", maxTi, maxTis);
   getAttribute(ele_param, "simulationTime", "parameter", simtime);
   std::cout<<"seed subType "<<subType<<"\n";
  //

}


std::string SeedRandomOrganParameter::writeXML(FILE* fp) const {

tinyxml2::XMLPrinter printer( fp, false, 0 ); //compact mode false, and 0 indent
    printer.OpenElement("organ", false); //compact mode false
    printer.PushAttribute("type","seed"); printer.PushAttribute("name", &name);printer.PushAttribute("subType",subType);
        printer.OpenElement("parameter");    printer.PushAttribute("name", "seedPos.x");    printer.PushAttribute("value", seedPos.x);  printer.CloseElement(false);
        printer.OpenElement("parameter");    printer.PushAttribute("name", "seedPos.y");    printer.PushAttribute("value", seedPos.y);  printer.CloseElement(false);
        printer.OpenElement("parameter");    printer.PushAttribute("name", "seedPos.z");    printer.PushAttribute("value", seedPos.z);  printer.CloseElement(false);

        printer.OpenElement("parameter");    printer.PushAttribute("name", "plantingdepth");    printer.PushAttribute("value", plantingdepth);  printer.CloseElement(false);
        printer.OpenElement("parameter");    printer.PushAttribute("name", "firstB");    printer.PushAttribute("value", firstB);  printer.PushAttribute("dev", firstBs); printer.CloseElement(false);
        printer.OpenElement("parameter");    printer.PushAttribute("name", "delayB");    printer.PushAttribute("value", delayB);  printer.PushAttribute("dev", delayBs); printer.CloseElement(false);
		printer.OpenElement("parameter");    printer.PushAttribute("name", "maxB");    printer.PushAttribute("value", maxB);  printer.PushAttribute("dev", maxBs); printer.CloseElement(false);
		printer.OpenElement("parameter");    printer.PushAttribute("name", "maxTi");    printer.PushAttribute("value", maxTi);  printer.PushAttribute("dev", maxTis); printer.CloseElement(false);
		printer.OpenElement("parameter");    printer.PushAttribute("name", "nC");    printer.PushAttribute("value", nC);  printer.CloseElement(false);
        printer.OpenElement("parameter");    printer.PushAttribute("name", "nz");    printer.PushAttribute("value", nz);  printer.CloseElement(false);
		printer.OpenElement("parameter");    printer.PushAttribute("name", "firstSB");    printer.PushAttribute("value", firstSB);  printer.PushAttribute("dev", firstSBs); printer.CloseElement(false);
		printer.OpenElement("parameter");    printer.PushAttribute("name", "delaySB");    printer.PushAttribute("value", delaySB);  printer.PushAttribute("dev", delaySBs); printer.CloseElement(false);
		printer.OpenElement("parameter");    printer.PushAttribute("name", "delayRC");    printer.PushAttribute("value", delayRC);  printer.PushAttribute("dev", delayRCs); printer.CloseElement(false);
        printer.OpenElement("parameter");    printer.PushAttribute("name", "simulationTime");    printer.PushAttribute("value", simtime);  printer.PushAttribute("dev", maxBs); printer.CloseElement(false);
    printer.CloseElement(false); // close element compact mode false

	return std::string(printer.CStr());

}

/**
 *
 */
SpecificOrganParamter* SeedRandomOrganParameter::realize() const
{
    //std::cout<<"seed name "<< <<"\n";


	SeedParameter* sp = new SeedParameter();
	sp->firstB = firstB + randn()*firstBs;
	sp->delayB = delayB + randn()*firstBs;
	sp->maxB = std::round(maxB + randn()*maxBs);
	sp->firstB = std::max(0.,firstB);
	sp->delayB = std::max(0.,delayB);
	sp->maxTi = std::round(maxTi + randn()*maxTis);
	sp->seedPos = seedPos;
	sp->nC = nC;
	sp->nz = nz;

	return sp;
}







/*************************************** Stem *****************************/


double StemParameter::getK() const
{
	double l = 0;
	for (auto const& dl : ln) {
		l += dl;
	}
	return l+la+lb;
}

void StemParameter::write(std::ostream & cout) const
{
	cout << "# Root Parameters \n";
	cout << "type\t" << subType << "\n"  << "lb\t"<< lb <<"\n" << "la\t"<< la <<"\n" << "ln\t";
	for (size_t i=0; i<ln.size(); i++) {
		cout << ln[i] << "\t";
	}
	cout << "\n" << "nob\t"<< ln.size() <<"\n" << "r\t"<< r <<"\n" << "a\t" << a << "\n" << "theta\t" << theta << "\n" << "rlt\t" << rlt << "\n";
}



StemRandomOrganParameter::StemRandomOrganParameter()
{
	organType = Organism::ot_stem;
	subType = -1; // means undefined
	tropism = new StemTropismFunction(0,0);
	growth = new StemGrowthFunction();
	se = new SoilLookUp();
	sa = new SoilLookUp();
	sbp = new SoilLookUp();
	set(-1, 0., 0., 10., 0., 1., 0., 0, 0., 0., 1., 0, 0.1, 0., 150./255.,150./255.,50./255., 1, 1. ,0.2, 0.1,
			successor, successorP, 1.22, 0., 1.e9, 0., 1, "undefined");
}

StemRandomOrganParameter::~StemRandomOrganParameter()
{
	delete tropism;
	delete growth;
	delete f_se;
	delete sa;
	delete sbp;
}

/**
 * todo comment
 */
void StemRandomOrganParameter::set(int type, double lb, double lbs, double la, double las, double ln, double lns, int inf, double nob, double nobs,
		double r, double rs, double a, double as,  double RotBeta, double BetaDev, double InitBeta, int tropismT, double tropismN, double tropismS,
		double dx, const std::vector<int>& successor, const std::vector<double>& successorP, double theta, double thetas, double rlt, double rlts,
		int gf, const std::string& name)
{
	this->subType = type;
	this->lb = lb;
	this->lbs = lbs;
	this->la = la;
	this->las = las;
	this->lns = lns;
	this->ln = ln;
	this->lnf = lnf;
	this->nob = nob;
	this->nobs = nobs;
	this->r = r;
	this->rs = rs;
	this->a = a;
	this->as = as;
	this->RotBeta = RotBeta;
	this->BetaDev = BetaDev;
	this->InitBeta = InitBeta;
	this->tropismT = tropismT;
	this->tropismN = tropismN;
	this->tropismS = tropismS;
	this->dx=dx;
	this->successor = successor; // vector
	this->successorP = successorP; // vector
	assert(successor.size()==successorP.size());
	this->theta = theta;
	this->thetas = thetas;
	this->rlt = rlt;
	this->rlts = rlts;
	this->gf = gf;
	this->name = name; // string

	createTropism();
	createGrowth();
}

void StemRandomOrganParameter::createTropism(SignedDistanceFunction* geom, SoilLookUp* soil)
{
	delete tropism;
	StemTropismFunction* t;
	switch (tropismT) {
	case tt_plagio: {
		t = new StemPlagiotropism(tropismN,tropismS);
		break;
	}
	case tt_gravi: {
		t = new StemGravitropism(tropismN,tropismS);
		break;
	}
	case tt_exo: {
		t = new StemExotropism(tropismN,tropismS);
		break;
	}
	case tt_hydro: {
		StemTropismFunction* gt =  new StemAntiGravitropism(tropismN,tropismS);
		StemTropismFunction* ht = new StemPlagiotropism(tropismN,tropismS);
		t = new CombinedStemTropism(tropismN,tropismS,ht,10.,gt,1.); // does only use the objective functions from gravitropism and hydrotropism
		break;
	}
    case tt_antigravi: {
		t = new StemAntiGravitropism(tropismN,tropismS);
		break;
    }
    case tt_twist: {
		t = new StemTwist(tropismN,tropismS);
		break;
    }
	default: throw std::invalid_argument( "StemSystem::createTropismFunction() tropism type not implemented" );
	}

	if (geom!=nullptr) {
		tropism = new ConfinedStemTropism(t,geom);
	} else {
		tropism = t;
	}

}

void StemRandomOrganParameter::createGrowth()
{
	switch (gf) {
	case gft_negexp: {
		growth = new StemExponentialGrowth();
		break;
	}
	case gft_linear: {
		growth =  new StemLinearGrowth();
		break;
	}
	default: throw std::invalid_argument( "RootSystem::createGrowthFunction() growth function type not implemented" );
	}
}

/**
 * Creates a specific root from the root type parameters.
 * The unique root id is not set, but must be set from outside.
 * (called by Stem::Stem())
 *
 * minimal ln distance is 1.e-9
 *
 * \return Specific stem parameters derived from the Stem type parameters
 */
SpecificOrganParamter* StemRandomOrganParameter::realize() const
{
	// type does not change
	double lb_ = std::max(lb + randn()*lbs,double(0)); // length of basal zone
	double la_ = std::max(la + randn()*las,double(0)); // length of apical zone
	std::vector<double> ln_; // stores the inter-distances
	int nob_ = std::max(round(nob + randn()*nobs),double(0)); // maximal number of branches
	switch(lnf) {
		case 0: // homogeneously distributed stem nodes
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d = std::max(ln + randn()*lns,1e-9); //Normal function of equal internode distance
			ln_.push_back(d);
			ln_.push_back(0);

		};break;
		case 1: //  nodes distance increase linearly
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d =  std::max(ln*(1+i) + randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);
			ln_.push_back(0);

		};break;
		case 2: //nodes distance decrease linearly
		for (int i = 0; i<nob_-1; i++) { // create inter-stem distances
			double d =  std::max(ln*(1+i) + randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);

		};break;
		case 3: //nodes distance decrease exponential
		for (int i = 0; i<nob_-1; i++) { // create inter-stem distances
			double d =  std::max(ln + randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);

		};break;

		case 4://nodes distance decrease exponential
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d =  std::max(ln/(1+i) + randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);
			ln_.push_back(0);
		}; break;
		case 5://nodes distance decrease exponential
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d =  std::max(ln/(1+i) + randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);
		};break;

	}
	double r_ = std::max(r + randn()*rs,double(0)); // initial elongation
	double a_ = std::max(a + randn()*as,double(0)); // radius
	double theta_ = std::max(theta + randn()*thetas,double(0)); // initial elongation
	double rlt_ = std::max(rlt + randn()*rlts,double(0)); // root life time
	int lnf_ = lnf;
	StemParameter* stem_p =  new StemParameter(subType,lb_,la_,ln_,r_,a_,theta_,rlt_,lnf_);
	return stem_p;
	}



/**
 * Choose (dice) lateral type based on root parameter set,
 * (since there can be more than one lateral type)
 *
 * @param pos       spatial position (for coupling to a soil model)
 */
int StemRandomOrganParameter::getLateralType(const Vector3d& pos)
{
	assert(successor.size()==successorP.size());
	double scale = sbp->getValue(pos);  //the current model makes not a lot of sense, we may come up with something more clever
	if (successorP.size()>0) { // at least 1 successor type
		if (successorP.size()>1) { // if there are more than one lateral we have to dice
			double d = rand();
			int i=0;
			double stem_p=successorP.at(i)*scale;
			while (stem_p<=d) {
				i++;
				stem_p+=successorP.at(i)*scale;
			}
			return successor.at(i);
		} else {
			return successor.front();
		}
	} else { // no successors
		return -1;
	}
}

void StemRandomOrganParameter::write(std::ostream & os) const
{
	os << "# Root type parameter for " << name << "\n";
	os << "type\t" << subType << "\n" << "name\t" << name << "\n" << "lb\t"<< lb <<"\t"<< lbs << "\n" << "la\t"<< la <<"\t"<< las << "\n"
			<< "ln\t" << ln << "\t" << lns << "\n" << "nob\t"<< nob <<"\t"<< nobs << "\n" << "r\t"<< r <<"\t"<< rs << "\n" <<
			"a\t" << a << "\t" << as << "\n" << "color\t"<< RotBeta <<"\t"<< BetaDev << "\t" << InitBeta << "\n"
			<< "tropism\t"<< tropismT <<"\t"<< tropismN << "\t" << tropismS << "\n" << "dx\t" << dx << "\n" << "successor\t" << successor.size() << "\t";
	for (size_t i=0; i<successor.size(); i++) {
		os << successor[i] << "\t";
	}
	os << "\n" << "successorP\t" << successorP.size() <<"\t";
	for (size_t i=0; i<successorP.size(); i++) {
		os << successorP[i] << "\t";
	}
	os << "\n" << "theta\t" << theta << "\t" << thetas << "\n" << "rlt\t" << rlt << "\t" << rlts << "\n" << "gf\t" << gf << "\n";
	std::cout << "RootRandomOrganParameter::write is deprecated, use RootRandomOrganParameter::writeXML instead\n";
}

void StemRandomOrganParameter::read(std::istream & is)
{
	char ch[256]; // dummy
	is.getline(ch,256);
	std::string s; // dummy
	is >> s >> subType >> s >> name >> s >> lb >> lbs >> s >> la >> las >> s >> ln >> lns >> s >> k >> ks;
	is >> s >> r >> rs >> s >> a >> as >> s >> RotBeta >> BetaDev >> InitBeta >> s >> tropismT >> tropismN >> tropismS >> s >> dx;
	if (ln > 0) {
		nob=  (k-la-lb)/ln+1;   //conversion, because the input file delivers the lmax value and not the nob value
		nob = std::max(nob,0.);
		nobs = (ks/k - lns/ln)*k/ln; // error propagation
		if (la>0) {
			nobs -= (las/la - lns/ln)*la/ln;
		}
		if (lb>0) {
			nobs -= (lbs/lb - lns/ln)*lb/ln;
		}
		nobs = std::max(nobs,0.);
		if (std::isnan(nobs)) {
			std::cout << "StemRandomOrganParameter::read() nobs is nan \n";
			nobs =0;
		}
	} else {
		nob=0;
		nobs = 0;
	}
	int n;
	is  >> s >> n;
	successor.clear();
	int is_;
	for (int i=0; i<n; i++) {
		is >> is_;
		successor.push_back(is_);
	}
	is >> s >> n;
	successorP.clear();
	double ds;
	for (int i=0; i<n; i++) {
		is >> ds;
		successorP.push_back(ds);
	}
	is >> s >> theta >> thetas >> s >> rlt >> rlts >> s >> gf >> s;
	createTropism();
	createGrowth();
	std::cout << "StemRandomOrganParameter::read is deprecated, use StemRandomOrganParameter::readXML instead\n";
}

void StemRandomOrganParameter::readXML(const tinyxml2::XMLElement* ele) //read subtype parameter from different organ type, used by Plant::openXML
{
   const tinyxml2::XMLElement* ele_param = ele->FirstChildElement("parameter"); //XML elements for parameters
   const char* name;
   ele->QueryUnsignedAttribute("subType", &subType);
   ele->QueryStringAttribute("name", &name);
   getAttribute(ele_param, "lb", "parameter", lb, lbs);
   getAttribute(ele_param, "la", "parameter", la, las);
   getAttribute(ele_param, "ln", "parameter", ln, lns, lnf);
   getAttribute(ele_param, "lmax", "parameter", k, ks);
   getAttribute(ele_param, "nob", "parameter", nob, nobs);
   	if (ln > 0) {
		nob=  (k-la-lb)/ln+1;   //conversion, because the input file delivers the lmax value and not the nob value
		nob = std::max(nob,0.);
		nobs = (ks/k - lns/ln)*k/ln; // error propagation
		if (la>0) {
			nobs -= (las/la - lns/ln)*la/ln;
		}
		if (lb>0) {
			nobs -= (lbs/lb - lns/ln)*lb/ln;
		}
		nobs = std::max(nobs,0.);
		if (std::isnan(nobs)) {
			std::cout << "RootRandomOrganParameter::read() nobs is nan \n";
			nobs =0;
		}
	} else {
		nob=0;
		nobs = 0;
	}
      // if (lnf = 1){ nob = nob/3;}
   getAttribute(ele_param, "r", "parameter", r, rs);
   getAttribute(ele_param, "a", "parameter", a, as);
   getAttribute(ele_param, "RotBeta", "parameter", RotBeta);
   getAttribute(ele_param, "BetaDev", "parameter", BetaDev);
   getAttribute(ele_param, "InitBeta", "parameter", InitBeta);
   getAttribute(ele_param, "tropismN", "parameter", tropismN);
   getAttribute(ele_param, "tropismT", "parameter", tropismT);
   getAttribute(ele_param, "tropismS", "parameter", tropismS);
      getAttribute(ele_param, "dx", "parameter", dx);
   getAttribute(ele_param, "theta", "parameter", theta, thetas);
   getAttribute(ele_param, "rlt", "parameter", rlt, rlts);
   getAttribute(ele_param, "gf", "parameter", gf);
   getAttribute(ele_param, "successor", "parameter", successor, successorP);
   std::cout<<"subType "<<subType<<"\n";
      std::cout<<"theta "<<theta<<"\n";
   	createTropism();
	createGrowth();
}


std::string StemRandomOrganParameter::writeXML(FILE* fp) const
{
tinyxml2::XMLPrinter printer( fp, false, 0 ); //compact mode false, and 0 indent
        printer.OpenElement("organ", false); //compact mode false
        printer.PushAttribute("type","stem");

	    switch (subType) {
	case 1 :  printer.PushAttribute("name","mainstem"); printer.PushAttribute("subType",subType);// See

	break;
    case 2 :  printer.PushAttribute("name","1stbranch"); printer.PushAttribute("subType",subType); // See
	break;
    case 3 :  printer.PushAttribute("name","2ndbranch"); printer.PushAttribute("subType",subType);// See
	break;
	case 4 :  printer.PushAttribute("name","otherbranch1"); printer.PushAttribute("subType",subType); // See
	break;
    case 5 :  printer.PushAttribute("name","otherbranch2"); printer.PushAttribute("subType",subType); // See
	break;
    }

    printer.OpenElement("parameter");

        printer.PushAttribute("name","lb"); printer.PushAttribute("value",lb); printer.PushAttribute("dev",lbs);  printer.CloseElement(); printer.PushComment("Basal zone [cm]");	//		 	///< Basal zone [cm]
	printer.OpenElement("parameter");

        printer.PushAttribute("name","la"); printer.PushAttribute("value",la); printer.PushAttribute("dev",las);  printer.CloseElement(); printer.PushComment("Apical zone [cm];");	///< Apical zone [cm];
	printer.OpenElement("parameter");

        printer.PushAttribute("name","ln"); printer.PushAttribute("value",ln); printer.PushAttribute("dev",lns); printer.PushAttribute("functiontype",lnf); printer.CloseElement();printer.PushComment("Inter-lateral distance [cm];");		///< Inter-lateral distance [cm];
	printer.OpenElement("parameter");//	printer.PushComment("Number of branches [1];");

	printer.PushAttribute("name","lmax"); printer.PushAttribute("value",k); printer.PushAttribute("dev",ks);  printer.CloseElement();printer.PushComment("Inter-lateral distance [cm];");		///< Inter-lateral distance [cm];
	printer.OpenElement("parameter");//	printer.PushComment("Number of branches [1];");

        printer.PushAttribute("name","nob"); printer.PushAttribute("value",nob); printer.PushAttribute("dev",nobs);  printer.CloseElement();		///< Standard deviation apical zone [cm];
	printer.OpenElement("parameter");//	printer.PushComment("Initial growth rate [cm day-1]");
        printer.PushAttribute("name","r"); printer.PushAttribute("value",r); printer.PushAttribute("dev",rs); printer.CloseElement();		///< Initial growth rate [cm day-1]
	printer.OpenElement("parameter");//	printer.PushComment("Root radius [cm]");
        printer.PushAttribute("name","a"); printer.PushAttribute("value",a); printer.PushAttribute("dev",as);   printer.CloseElement();		///< Root radius [cm]

	printer.OpenElement("parameter");//
        printer.PushAttribute("name","RotBeta");  printer.PushAttribute("value",RotBeta); printer.CloseElement();
    printer.OpenElement("parameter");
        printer.PushAttribute("name","BetaDev");  printer.PushAttribute("value",BetaDev); printer.CloseElement();
    printer.OpenElement("parameter");
        printer.PushAttribute("name","InitBeta");  printer.PushAttribute("value",InitBeta); printer.CloseElement();	///< Root color (red)

	printer.OpenElement("parameter");//	printer.PushComment("Root tropism parameter (Type, number of trials, mean vale of expected change)");
        printer.PushAttribute("name","tropismT");  printer.PushAttribute("value",tropismT);  printer.CloseElement();
    printer.OpenElement("parameter");
        printer.PushAttribute("name","tropismN");  printer.PushAttribute("value",tropismN);  printer.CloseElement();
    printer.OpenElement("parameter");
        printer.PushAttribute("name","tropismS");  printer.PushAttribute("value",tropismS);  printer.CloseElement();	///< Root tropism parameter (Type)

	printer.OpenElement("parameter");//	printer.PushComment("Maximal segment size");
        printer.PushAttribute("name","dx"); printer.PushAttribute("value",dx);  printer.CloseElement();		///< Maximal segment size [cm]
	printer.OpenElement("parameter");//	printer.PushComment("Angle between root and parent root");
        printer.PushAttribute("name","theta"); printer.PushAttribute("value",theta); printer.PushAttribute("dev",thetas); printer.CloseElement();	///< Angle between root and parent root (rad)
	printer.OpenElement("parameter");//	printer.PushComment("Root life time");
        printer.PushAttribute("name","rlt"); printer.PushAttribute("value",rlt); printer.PushAttribute("dev",rlts);  printer.CloseElement();	///< Root life time (days)
	printer.OpenElement("parameter");//	printer.PushComment("Growth function");
        printer.PushAttribute("name","gf"); printer.PushAttribute("value",gf); printer.CloseElement();		///< Growth function (1=negative exponential, 2=linear)


	for (int successorCount = 0; successorCount < successor.size(); successorCount ++ ){
            std::string s = std::to_string(successorCount);
            char const *schar = s.c_str();
        printer.OpenElement("parameter"); printer.PushAttribute("name","successor"); printer.PushAttribute("number",schar); printer.PushAttribute("type",successor[successorCount]); printer.PushAttribute("percentage",successorP[successorCount]);    printer.CloseElement();
//        printer.OpenElement("parameter"); printer.PushAttribute("name","successorP"); printer.PushAttribute("number",schar); printer.PushAttribute("percentage",successorP[successorCount]);   printer.CloseElement();
//
    }
//    printer.PrintSpace();
    printer.CloseElement(false); // close element compact mode false

	return std::string(printer.CStr());


}



///**********************************************************Another COPY PASTE to create leaf******************************************************************************************

double LeafParameter::getK() const
{
	double l = 0;
	for (auto const& dl : ln) {
		l += dl;
	}
	return l+la+lb;
}

void LeafParameter::write(std::ostream & cout) const
{
	cout << "# Root Parameters \n";
	cout << "type\t" << subType << "\n"  << "lb\t"<< lb <<"\n" << "la\t"<< la <<"\n" << "ln\t";
	for (size_t i=0; i<ln.size(); i++) {
		cout << ln[i] << "\t";
	}
	cout << "\n" << "nob\t"<< ln.size() <<"\n" << "r\t"<< r <<"\n" << "a\t" << a << "\n" << "theta\t" << theta << "\n" << "rlt\t" << rlt << "\n";
}



LeafRandomOrganParameter::LeafRandomOrganParameter()
{
	organType = Organ::ot_leafe;
	subType = -1; // means undefined
	tropism = new LeafTropismFunction(0,0);
	growth = new LeafGrowthFunction();
	se = new SoilLookUp();
	sa = new SoilLookUp();
	sbp = new SoilLookUp();
	set(-1, 0., 0., 10., 0., 1., 0., 0, 0., 0., 1., 0, 0.1, 0., 150./255.,150./255.,50./255., 1, 1. ,0.2, 0.1,
			successor, successorP, 1.22, 0., 1.e9, 0., 1, "undefined", "sth");
}

LeafRandomOrganParameter::~LeafRandomOrganParameter()
{
	delete tropism;
	delete growth;
	delete f_se;
	delete sa;
	delete sbp;
}

/**
 * todo comment
 */
void LeafRandomOrganParameter::set(int type, double lb, double lbs, double la, double las, double ln, double lns, int inf, double nob, double nobs,
		double r, double rs, double a, double as,  double RotBeta, double BetaDev, double InitBeta, int tropismT, double tropismN, double tropismS,
		double dx, const std::vector<int>& successor, const std::vector<double>& successorP, double theta, double thetas, double rlt, double rlts,
		int gf, const std::string& name,const char* organName)
{
	this->subType = type;
	this->lb = lb;
	this->lbs = lbs;
	this->la = la;
	this->las = las;
	this->ln = ln;
	this->lnf = lnf;
	this->lns = lns;
	this->nob = nob;
	this->nobs = nobs;
	this->r = r;
	this->rs = rs;
	this->a = a;
	this->as = as;
	this->RotBeta = RotBeta;
	this->BetaDev = BetaDev;
	this->InitBeta = InitBeta;
	this->tropismT = tropismT;
	this->tropismN = tropismN;
	this->tropismS = tropismS;
	this->dx=dx;
	this->successor = successor; // vector
	this->successorP = successorP; // vector
	assert(successor.size()==successorP.size());
	this->theta = theta;
	this->thetas = thetas;
	this->rlt = rlt;
	this->rlts = rlts;
	this->gf = gf;
	this->name = name; // string
	std::cout<<this->name<<"\n";
    this->organName =organName;
	createTropism();
	createGrowth();
}

void LeafRandomOrganParameter::createTropism(SignedDistanceFunction* geom, SoilLookUp* soil)
{
	delete tropism;
	LeafTropismFunction* t;
	switch (tropismT) {
	case tt_plagio: {
		t = new LeafPlagiotropism(tropismN,tropismS);
		break;
	}
	case tt_gravi: {
		t = new LeafGravitropism(tropismN,tropismS);
		break;
	}
	case tt_exo: {
		t = new LeafExotropism(tropismN,tropismS);
		break;
	}
	case tt_hydro: {
		LeafTropismFunction* gt =  new LeafGravitropism(tropismN,tropismS);
		LeafTropismFunction* ht = new LeafPhototropism(tropismN,tropismS,soil);
		t = new CombinedLeafTropism(tropismN,tropismS,ht,10.,gt,1.); // does only use the objective functions from gravitropism and hydrotropism
		break;
	}
	case tt_antigravi: {
		t = new LeafAntiGravitropism(tropismN,tropismS);
		break;
    }
	default: throw std::invalid_argument( "LeafSystem::createTropismFunction() tropism type not implemented" );
	}

	if (geom!=nullptr) {
		tropism = new ConfinedLeafTropism(t,geom);
	} else {
		tropism = t;
	}

}

void LeafRandomOrganParameter::createGrowth()
{
	switch (gf) {
	case gft_negexp: {
		growth = new LeafExponentialGrowth();
		break;
	}
	case gft_linear: {
		growth =  new LeafLinearGrowth();
		break;
	}
	default: throw std::invalid_argument( "Plant::createGrowthFunction() growth function type not implemented" );
	}
}

/**
 * Creates a specific root from the root type parameters.
 * The unique root id is not set, but must be set from outside.
 * (called by Root::Root())
 *
 * minimal ln distance is 1.e-9
 *
 * \return Specific root parameters derived from the root type parameters
 */
SpecificOrganParamter* LeafRandomOrganParameter::realize() const
{
	// type does not change
	double lb_ = std::max(lb + randn()*lbs,double(0)); // length of basal zone
	double la_ = std::max(la + randn()*las,double(0)); // length of apical zone
	std::vector<double> ln_; // stores the inter-distances

	// stores the inter-distances
	int nob_ = std::max(round(nob + randn()*nobs),double(0)); // maximal number of branches
	switch(lnf) {
		case 0: // homogeneously distributed stem nodes
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d = std::max(ln + randn()*lns,1e-9); //Normal function of equal internode distance
			ln_.push_back(d);
			ln_.push_back(0);

		};
		case 1: //  nodes distance increase linearly
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d =  std::max(ln*(1+i) + randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);
			ln_.push_back(0);

		};
		case 2: //nodes distance decrease linearly
		for (int i = 0; i<nob_-1; i++) { // create inter-stem distances
			double d =  std::max(ln*(1+i) + randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);

		};
		case 3: //nodes distance increase exponential
		for (int i = 0; i<nob_-1; i++) { // create inter-stem distances
			double d =  std::max(ln + randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);

		};

		case 4://nodes distance decrease exponential
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d =  std::max(ln/(1+i) + randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);
			ln_.push_back(0);
		};
		case 5://nodes distance decrease exponential
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d =  std::max(ln/(1+i) + randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);
		};

	}


	double r_ = std::max(r + randn()*rs,double(0)); // initial elongation
	double a_ = std::max(a + randn()*as,double(0)); // radius
	double theta_ = std::max(theta + randn()*thetas,double(0)); // initial elongation
	double rlt_ = std::max(rlt + randn()*rlts,double(0)); // root life time
	int lnf_ = lnf;
	LeafParameter* leaf_p =  new LeafParameter(subType,lb_,la_,ln_,r_,a_,theta_,rlt_,lnf_);
	return leaf_p;

//for (int i = 0; i<nob_-1; i++) { // create inter-root distances
//		double d = 0.01*2*i; //std::max(  );//ln + randn()*lns,1e-9);
//		ln_.push_back(d);
//	}

}

/**
 * Choose (dice) lateral type based on root parameter set,
 * (since there can be more than one lateral type)
 *
 * @param pos       spatial position (for coupling to a soil model)
 */
int LeafRandomOrganParameter::getLateralType(const Vector3d& pos)
{
	assert(successor.size()==successorP.size());
	double scale = sbp->getValue(pos);  //the current model makes not a lot of sense, we may come up with something more clever
	if (successorP.size()>0) { // at least 1 successor type
		if (successorP.size()>1) { // if there are more than one lateral we have to dice
			double d = rand();
			int i=0;
			double leaf_p=successorP.at(i)*scale;
			while (leaf_p<=d) {
				i++;
				leaf_p+=successorP.at(i)*scale;
			}
			return successor.at(i);
		} else {
			return successor.front();
		}
	} else { // no successors
		return -1;
	}
}

void LeafRandomOrganParameter::write(std::ostream & os) const
{
	os << "# Root type parameter for " << name << "\n";
	os << "type\t" << subType << "\n" << "name\t" << name << "\n" << "lb\t"<< lb <<"\t"<< lbs << "\n" << "la\t"<< la <<"\t"<< las << "\n"
			<< "ln\t" << ln << "\t" << lns << "\n" << "nob\t"<< nob <<"\t"<< nobs << "\n" << "r\t"<< r <<"\t"<< rs << "\n" <<
			"a\t" << a << "\t" << as << "\n" << "color\t"<< RotBeta <<"\t"<< BetaDev << "\t" << InitBeta << "\n"
			<< "tropism\t"<< tropismT <<"\t"<< tropismN << "\t" << tropismS << "\n" << "dx\t" << dx << "\n" << "successor\t" << successor.size() << "\t";
	for (size_t i=0; i<successor.size(); i++) {
		os << successor[i] << "\t";
	}
	os << "\n" << "successorP\t" << successorP.size() <<"\t";
	for (size_t i=0; i<successorP.size(); i++) {
		os << successorP[i] << "\t";
	}
	os << "\n" << "theta\t" << theta << "\t" << thetas << "\n" << "rlt\t" << rlt << "\t" << rlts << "\n" << "gf\t" << gf << "\n";
	std::cout << "RootRandomOrganParameter::write is deprecated, use RootRandomOrganParameter::writeXML instead\n";
}

void LeafRandomOrganParameter::read(std::istream & is)
{
	char ch[256]; // dummy
	is.getline(ch,256);
	std::string s; // dummy
	is >> s >> subType >> s >> name >> s >> lb >> lbs >> s >> la >> las >> s >> ln >> lns >> s >> k >> ks;
	is >> s >> r >> rs >> s >> a >> as >> s >> RotBeta >> BetaDev >> InitBeta >> s >> tropismT >> tropismN >> tropismS >> s >> dx;
	if (ln > 0) {
		nob=  (k-la-lb)/ln+1;   //conversion, because the input file delivers the lmax value and not the nob value
		nob = std::max(nob,0.);
		nobs = (ks/k - lns/ln)*k/ln; // error propagation
		if (la>0) {
			nobs -= (las/la - lns/ln)*la/ln;
		}
		if (lb>0) {
			nobs -= (lbs/lb - lns/ln)*lb/ln;
		}
		nobs = std::max(nobs,0.);
		if (std::isnan(nobs)) {
			std::cout << "LeafRandomOrganParameter::read() nobs is nan \n";
			nobs =0;
		}
	} else {
		nob=0;
		nobs = 0;
	}
	int n;
	is  >> s >> n;
	successor.clear();
	int is_;
	for (int i=0; i<n; i++) {
		is >> is_;
		successor.push_back(is_);
	}
	is >> s >> n;
	successorP.clear();
	double ds;
	for (int i=0; i<n; i++) {
		is >> ds;
		successorP.push_back(ds);
	}
	is >> s >> theta >> thetas >> s >> rlt >> rlts >> s >> gf >> s;
	createTropism();
	createGrowth();
	std::cout << "LeafRandomOrganParameter::read is deprecated, use LeafRandomOrganParameter::readXML instead\n";
}

void LeafRandomOrganParameter::readXML(const tinyxml2::XMLElement* ele) //read subtype parameter from different organ type, used by Plant::openXML
{
   const tinyxml2::XMLElement* ele_param = ele->FirstChildElement("parameter"); //XML elements for parameters

   ele->QueryUnsignedAttribute("subType", &subType);
   ele->QueryStringAttribute("name", &organName);
   std::cout<<"ltp xml read is"<< organName<<"\n";
   name = std::string(organName);
   getAttribute(ele_param, "lb", "parameter", lb, lbs);
   getAttribute(ele_param, "la", "parameter", la, las);
   getAttribute(ele_param, "ln", "parameter", ln, lns, lnf);
   getAttribute(ele_param, "lmax", "parameter", k, ks);
   getAttribute(ele_param, "nob", "parameter", nob, nobs);
   	if (ln > 0) {
		nob=  (k-la-lb)/ln+1;   //conversion, because the input file delivers the lmax value and not the nob value
		nob = std::max(nob,0.);
		nobs = (ks/k - lns/ln)*k/ln; // error propagation
		if (la>0) {
			nobs -= (las/la - lns/ln)*la/ln;
		}
		if (lb>0) {
			nobs -= (lbs/lb - lns/ln)*lb/ln;
		}
		nobs = std::max(nobs,0.);
		if (std::isnan(nobs)) {
			std::cout << "RootRandomOrganParameter::read() nobs is nan \n";
			nobs =0;
		}
	} else {
		nob=0;
		nobs = 0;
	}
   getAttribute(ele_param, "r", "parameter", r, rs);
   getAttribute(ele_param, "a", "parameter", a, as);
   getAttribute(ele_param, "RotBeta", "parameter", RotBeta);
   getAttribute(ele_param, "BetaDev", "parameter", BetaDev);
   getAttribute(ele_param, "InitBeta", "parameter", InitBeta);
   getAttribute(ele_param, "tropismN", "parameter", tropismN);
   getAttribute(ele_param, "tropismT", "parameter", tropismT);
   getAttribute(ele_param, "tropismS", "parameter", tropismS);
      getAttribute(ele_param, "dx", "parameter", dx);
   getAttribute(ele_param, "theta", "parameter", theta, thetas);
   getAttribute(ele_param, "rlt", "parameter", rlt, rlts);
   getAttribute(ele_param, "gf", "parameter", gf);
   getAttribute(ele_param, "successor", "parameter", successor, successorP);
   std::cout<<"subType "<<subType<<"\n";
   	createTropism();
	createGrowth();
}


std::string LeafRandomOrganParameter::writeXML(FILE* fp) const
{
tinyxml2::XMLPrinter printer( fp, false, 0 ); //compact mode false, and 0 indent
        printer.OpenElement("organ", false); //compact mode false
        printer.PushAttribute("type","leaf");

	    switch (subType) {
	case 1 :  printer.PushAttribute("name",organName); printer.PushAttribute("subType",subType);// See

	break;
    case 2 :  printer.PushAttribute("name",organName); printer.PushAttribute("subType",subType); // See
	break;
    case 3 :  printer.PushAttribute("name",organName); printer.PushAttribute("subType",subType);// See
	break;
	case 4 :  printer.PushAttribute("name",organName); printer.PushAttribute("subType",subType); // See
	break;
    case 5 :  printer.PushAttribute("name",organName); printer.PushAttribute("subType",subType); // See
	break;
    }

    printer.OpenElement("parameter");

        printer.PushAttribute("name","lb"); printer.PushAttribute("value",lb); printer.PushAttribute("dev",lbs);  printer.CloseElement(); printer.PushComment("Basal zone [cm]");	//		 	///< Basal zone [cm]
	printer.OpenElement("parameter");

        printer.PushAttribute("name","la"); printer.PushAttribute("value",la); printer.PushAttribute("dev",las);  printer.CloseElement(); printer.PushComment("Apical zone [cm];");	///< Apical zone [cm];
	printer.OpenElement("parameter");

        printer.PushAttribute("name","ln"); printer.PushAttribute("value",ln); printer.PushAttribute("dev",lns); printer.PushAttribute("functiontype",lnf); printer.CloseElement();printer.PushComment("Inter-lateral distance [cm];");		///< Inter-lateral distance [cm];
	printer.OpenElement("parameter");//	printer.PushComment("Number of branches [1];");

	printer.PushAttribute("name","lmax"); printer.PushAttribute("value",k); printer.PushAttribute("dev",ks);  printer.CloseElement();printer.PushComment("Inter-lateral distance [cm];");		///< Inter-lateral distance [cm];
	printer.OpenElement("parameter");//	printer.PushComment("Number of branches [1];");

        printer.PushAttribute("name","nob"); printer.PushAttribute("value",nob); printer.PushAttribute("dev",nobs);  printer.CloseElement();		///< Standard deviation apical zone [cm];
	printer.OpenElement("parameter");//	printer.PushComment("Initial growth rate [cm day-1]");
        printer.PushAttribute("name","r"); printer.PushAttribute("value",r); printer.PushAttribute("dev",rs); printer.CloseElement();		///< Initial growth rate [cm day-1]
	printer.OpenElement("parameter");//	printer.PushComment("Root radius [cm]");
        printer.PushAttribute("name","a"); printer.PushAttribute("value",a); printer.PushAttribute("dev",as);   printer.CloseElement();		///< Root radius [cm]

	printer.OpenElement("parameter");//
        printer.PushAttribute("name","RotBeta");  printer.PushAttribute("value",RotBeta); printer.CloseElement();
    printer.OpenElement("parameter");
        printer.PushAttribute("name","BetaDev");  printer.PushAttribute("value",BetaDev); printer.CloseElement();
    printer.OpenElement("parameter");
        printer.PushAttribute("name","InitBeta");  printer.PushAttribute("value",InitBeta); printer.CloseElement();	///< Root color (red)

	printer.OpenElement("parameter");//	printer.PushComment("Root tropism parameter (Type, number of trials, mean vale of expected change)");
        printer.PushAttribute("name","tropismT");  printer.PushAttribute("value",tropismT);  printer.CloseElement();
    printer.OpenElement("parameter");
        printer.PushAttribute("name","tropismN");  printer.PushAttribute("value",tropismN);  printer.CloseElement();
    printer.OpenElement("parameter");
        printer.PushAttribute("name","tropismS");  printer.PushAttribute("value",tropismS);  printer.CloseElement();	///< Root tropism parameter (Type)

	printer.OpenElement("parameter");//	printer.PushComment("Maximal segment size");
        printer.PushAttribute("name","dx"); printer.PushAttribute("value",dx);  printer.CloseElement();		///< Maximal segment size [cm]
	printer.OpenElement("parameter");//	printer.PushComment("Angle between root and parent root");
        printer.PushAttribute("name","theta"); printer.PushAttribute("value",theta); printer.PushAttribute("dev",thetas); printer.CloseElement();	///< Angle between root and parent root (rad)
	printer.OpenElement("parameter");//	printer.PushComment("Root life time");
        printer.PushAttribute("name","rlt"); printer.PushAttribute("value",rlt); printer.PushAttribute("dev",rlts);  printer.CloseElement();	///< Root life time (days)
	printer.OpenElement("parameter");//	printer.PushComment("Growth function");
        printer.PushAttribute("name","gf"); printer.PushAttribute("value",gf); printer.CloseElement();		///< Growth function (1=negative exponential, 2=linear)


	for (int successorCount = 0; successorCount < successor.size(); successorCount ++ ){
            std::string s = std::to_string(successorCount);
            char const *schar = s.c_str();
        printer.OpenElement("parameter"); printer.PushAttribute("name","successor"); printer.PushAttribute("number",schar); printer.PushAttribute("type",successor[successorCount]); printer.PushAttribute("percentage",successorP[successorCount]);    printer.CloseElement();
//        printer.OpenElement("parameter"); printer.PushAttribute("name","successorP"); printer.PushAttribute("number",schar); printer.PushAttribute("percentage",successorP[successorCount]);   printer.CloseElement();
//
    }
//    printer.PrintSpace();
    printer.CloseElement(false); // close element compact mode false

	return std::string(printer.CStr());


}



} // namespace CPlantBox
