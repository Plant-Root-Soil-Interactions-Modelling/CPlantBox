#include "ModelParameter.h"
#include "tinyxml2.h"
#include "Plant.h"

#include <cmath>



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



OrganTypeParameter::OrganTypeParameter()
{
	organType = Organ::ot_organ;
	subType = -1; // means undefined
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



RootTypeParameter::RootTypeParameter()
{
	organType = Organ::ot_root;
	subType = -1; // means undefined
	tropism = new TropismFunction(0,0);
	growth = new GrowthFunction();
	se = new SoilProperty();
	sa = new SoilProperty();
	sbp = new SoilProperty();
	set(-1, 0., 0., 10., 0., 1., 0., 0., 0., 1., 0, 0.1, 0., 150./255.,150./255.,50./255., 1, 1. ,0.2, 0.1,
			successor, successorP, 1.22, 0., 1.e9, 0., 1, "undefined");
}

RootTypeParameter::~RootTypeParameter()
{
	delete tropism;
	delete growth;
	delete se;
	delete sa;
	delete sbp;
}

/**
 * todo comment
 */
void RootTypeParameter::set(int type, double lb, double lbs, double la, double las, double ln, double lns, double nob, double nobs,
		double r, double rs, double a, double as,  double colorR, double colorG, double colorB, int tropismT, double tropismN, double tropismS,
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
	this->colorR = colorR;
	this->colorG = colorG;
	this->colorB = colorB;
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

void RootTypeParameter::createTropism(SignedDistanceFunction* geom, SoilProperty* soil)
{
	delete tropism;
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

void RootTypeParameter::createGrowth()
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
OrganParameter* RootTypeParameter::realize() const
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
int RootTypeParameter::getLateralType(const Vector3d& pos)
{
	assert(successor.size()==successorP.size());
	double scale = sbp->getValue(pos);  //the current model makes not a lot of sense, we may come up with something more clever
	if (successorP.size()>0) { // at least 1 successor type
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

void RootTypeParameter::write(std::ostream & os) const
{
	os << "# Root type parameter for " << name << "\n";
	os << "type\t" << subType << "\n" << "name\t" << name << "\n" << "lb\t"<< lb <<"\t"<< lbs << "\n" << "la\t"<< la <<"\t"<< las << "\n"
			<< "ln\t" << ln << "\t" << lns << "\n" << "nob\t"<< nob <<"\t"<< nobs << "\n" << "r\t"<< r <<"\t"<< rs << "\n" <<
			"a\t" << a << "\t" << as << "\n" << "color\t"<< colorR <<"\t"<< colorG << "\t" << colorB << "\n"
			<< "tropism\t"<< tropismT <<"\t"<< tropismN << "\t" << tropismS << "\n" << "dx\t" << dx << "\n" << "successor\t" << successor.size() << "\t";
	for (size_t i=0; i<successor.size(); i++) {
		os << successor[i] << "\t";
	}
	os << "\n" << "successorP\t" << successorP.size() <<"\t";
	for (size_t i=0; i<successorP.size(); i++) {
		os << successorP[i] << "\t";
	}
	os << "\n" << "theta\t" << theta << "\t" << thetas << "\n" << "rlt\t" << rlt << "\t" << rlts << "\n" << "gf\t" << gf << "\n";
	std::cout << "RootTypeParameter::write is deprecated, use RootTypeParameter::writeXML instead\n";
}

void RootTypeParameter::read(std::istream & is)
{
	char ch[256]; // dummy
	is.getline(ch,256);
	std::string s; // dummy
	double k;
	double ks;
	is >> s >> subType >> s >> name >> s >> lb >> lbs >> s >> la >> las >> s >> ln >> lns >> s >> k >> ks;
	is >> s >> r >> rs >> s >> a >> as >> s >> colorR >> colorG >> colorB >> s >> tropismT >> tropismN >> tropismS >> s >> dx;
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
			std::cout << "RootTypeParameter::read() nobs is nan \n";
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
	std::cout << "RootTypeParameter::read is deprecated, use RootTypeParameter::readXML instead\n";
}

void RootTypeParameter::readXML(FILE* fp)
{

}


std::string RootTypeParameter::writeXML(FILE* fp) const
{
tinyxml2::XMLPrinter printer( fp, false, 0 ); //compact mode false, and 0 indent
        printer.OpenElement("organ", false); //compact mode false
        printer.PushAttribute("type","root");

	    switch (subType) {
	case 1 :  printer.PushAttribute("name","taproot"); printer.PushAttribute("id","id to discuss");// See

	break;
    case 2 :  printer.PushAttribute("name","lateral1"); printer.PushAttribute("id","id to discuss"); // See
	break;
    case 3 :  printer.PushAttribute("name","lateral2"); printer.PushAttribute("id","id to discuss");// See
	break;
	case 4 :  printer.PushAttribute("name","nodal_root"); printer.PushAttribute("id","id to discuss"); // See
	break;
    case 5 :  printer.PushAttribute("name","shoot_borne_root"); printer.PushAttribute("id","id to discuss"); // See
	break;
    }

    printer.OpenElement("parameter");
        printer.PushComment("Basal zone [cm]");
        printer.PushAttribute("name","lb"); printer.PushAttribute("value",lb); printer.PushAttribute("dev",lbs);  printer.CloseElement();	 	///< Basal zone [cm]
	printer.OpenElement("parameter");
		printer.PushComment("Apical zone [cm];");
        printer.PushAttribute("name","la"); printer.PushAttribute("value",la); printer.PushAttribute("dev",las);  printer.CloseElement();	///< Apical zone [cm];
	printer.OpenElement("parameter");
		printer.PushComment("Inter-lateral distance [cm];");
        printer.PushAttribute("name","ln"); printer.PushAttribute("value",ln); printer.PushAttribute("dev",lns);  printer.CloseElement();		///< Inter-lateral distance [cm];
	printer.OpenElement("parameter");//	printer.PushComment("Number of branches [1];");
        printer.PushAttribute("name","nob"); printer.PushAttribute("value",nob); printer.PushAttribute("dev",nobs);  printer.CloseElement();		///< Standard deviation apical zone [cm];
	printer.OpenElement("parameter");//	printer.PushComment("Initial growth rate [cm day-1]");
        printer.PushAttribute("name","r"); printer.PushAttribute("value",r); printer.PushAttribute("dev",rs); printer.CloseElement();		///< Initial growth rate [cm day-1]
	printer.OpenElement("parameter");//	printer.PushComment("Root radius [cm]");
        printer.PushAttribute("name","a"); printer.PushAttribute("value",a); printer.PushAttribute("dev",as);   printer.CloseElement();		///< Root radius [cm]
	printer.OpenElement("parameter");//	printer.PushComment("Root color (red, green, blue)");
        printer.PushAttribute("colorR",colorR); printer.PushAttribute("colorG",colorG); printer.PushAttribute("colorB",colorB); printer.CloseElement();	///< Root color (red)
	printer.OpenElement("parameter");//	printer.PushComment("Root tropism parameter (Type, number of trials, mean vale of expected change)");
        printer.PushAttribute("tropismT",tropismT); printer.PushAttribute("tropismN",tropismN); printer.PushAttribute("tropismS",tropismS);  printer.CloseElement();	///< Root tropism parameter (Type)
	printer.OpenElement("parameter");//	printer.PushComment("Maximal segment size");
        printer.PushAttribute("name","dx"); printer.PushAttribute("value",dx);  printer.CloseElement();		///< Maximal segment size [cm]
	printer.OpenElement("parameter");//	printer.PushComment("Angle between root and parent root");
        printer.PushAttribute("name","theta"); printer.PushAttribute("value",theta); printer.PushAttribute("dev",thetas); printer.CloseElement();	///< Angle between root and parent root (rad)
	printer.OpenElement("parameter");//	printer.PushComment("Root life time");
        printer.PushAttribute("name","rlt"); printer.PushAttribute("value",rlt); printer.PushAttribute("dev",rlts);  printer.CloseElement();	///< Root life time (days)
	printer.OpenElement("parameter");//	printer.PushComment("Growth function");
        printer.PushAttribute("name","dx"); printer.PushAttribute("value",dx); printer.CloseElement();		///< Growth function (1=negative exponential, 2=linear)


	for (int successorCount = 0; successorCount < successor.size(); successorCount ++ ){
            std::string s = std::to_string(successorCount);
            char const *schar = s.c_str();
        printer.OpenElement("successor"); printer.PushAttribute("name",schar); printer.PushAttribute("location",successor[successorCount]);    printer.CloseElement();
        printer.OpenElement("successorP"); printer.PushAttribute("name",schar); printer.PushAttribute("percentage",successorP[successorCount]);   printer.CloseElement();
//
    }
//    printer.PrintSpace();
    printer.CloseElement(false); // close element compact mode false

	return std::string(printer.CStr());


}



SeedTypeParameter::SeedTypeParameter() {
	organType = Organ::ot_seed;
	subType = 0; // there is currently only one type
}

void SeedTypeParameter::read(std::istream & is) {
	double plantingdepth;
	std::string s; // dummy
	is  >>  s >> plantingdepth;
	is >> s >> firstB >> s >> delayB >> s >> maxB >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s >> s;
	seedPos = Vector3d(0,0,-plantingdepth);
	std::cout << "SeedTypeParamter::read is deprecated, use SeedTypeParamter::readXML instead\n";
}


void SeedTypeParameter::write(std::ostream & cout) const {
	double pd = -seedPos.z;
	cout <<  "plantingdepth\t" << pd << "\n" <<  "firstB\t" << firstB << "\n" <<  "delayB\t" << delayB << "\n"
			<<  "maxB\t" << maxB << "\n" <<  "nC\t" << 0 << "\n" <<  "firstSB\t" << 0 << "\n"
			<<  "delaySB\t" << 0 << "\n" <<  "delayRC\t" << 0 << "\n" <<  "nz\t" << 0 << "\n" << "simulationTime\t" << 0 << "\n";
	std::cout << "SeedTypeParamter::write is deprecated, use SeedTypeParamter::writeXML instead\n";
}

void SeedTypeParameter::readXML(FILE* fp) {
//	tinyxml2::XMLDocument doc( fp );
//	tinyxml2::XMLElement* seed = doc.FirstChildElement( "Seed" );
//	seedPos = Vector3d(0,0,-3);
//	seedPoss = Vector3d();
//	tinyxml2::XMLElement* pos = seed->FirstChildElement("Location");
//	if (pos!=0) {
//		readXMLvs(pos,"x",&seedPos.x,&seedPoss.x);
//		readXMLvs(pos,"y",&seedPos.y,&seedPoss.y);
//		readXMLvs(pos,"z",&seedPos.z,&seedPoss.z);
//	}
//
//	firstB = 1.e9;
//	firstBs = 0.;
//	delayB = 1.e9;
//	delayBs =0.;
//	maxB = 0;
//	maxBs = 0;
//	tinyxml2::XMLElement* basal = seed->FirstChildElement("Basal roots");
//	if (basal!=0) {
//		readXMLvs(basal,"First",&firstB,&firstBs);
//		readXMLvs(basal,"Delay",&delayB,&delayBs);
//		double dmB=0;
//		readXMLvs(basal,"Maximum",&dmB,&maxBs);
//		maxB = int(dmB);
//	}
}

std::string SeedTypeParameter::writeXML(FILE* fp) const {

tinyxml2::XMLPrinter printer( fp, false, 0 ); //compact mode false, and 0 indent
    printer.OpenElement("organ", false); //compact mode false
    printer.PushAttribute("type","seed");
        printer.OpenElement("parameter");
        printer.PushAttribute("location_x", "0");
        printer.PushAttribute("location_y", "0");
        printer.PushAttribute("location_z", "0");
    printer.CloseElement(false); // close element compact mode false
    printer.CloseElement(false); // close element compact mode false

	return std::string(printer.CStr());
//    tinyxml2::XMLPrinter printer( fp );
//	printer.OpenElement("Seed");
//
//	printer.OpenElement("Location");
//	printer.OpenElement("x");
//	printer.PushText(seedPos.x);
//	if (seedPoss.x!=0) {
//		printer.OpenElement("sd"); printer.PushText(seedPoss.x); printer.CloseElement();
//	}
//	printer.CloseElement();
//	printer.OpenElement("y");
//	printer.PushText(seedPos.y);
//	if (seedPoss.y!=0) {
//		printer.OpenElement("sd"); printer.PushText(seedPoss.y); printer.CloseElement();
//	}
//	printer.CloseElement();
//	printer.OpenElement("z");
//	printer.PushText(seedPos.z);
//	if (seedPoss.y!=0) {
//		printer.OpenElement("sd"); printer.PushText(seedPoss.z); printer.CloseElement();
//	}
//	printer.CloseElement();
//	printer.CloseElement(); // Location
//
//	printer.OpenElement("Basal");
//	printer.OpenElement("First");
//	printer.PushText(firstB);
//	if (firstBs!=0) {
//		printer.OpenElement("sd"); printer.PushText(firstBs); printer.CloseElement();
//	}
//	printer.CloseElement();
//	printer.OpenElement("Delay");
//	printer.PushText(delayB);
//	if (delayBs!=0) {
//		printer.OpenElement("sd"); printer.PushText(delayBs); printer.CloseElement();
//	}
//	printer.CloseElement();
//	printer.OpenElement("Maximum");
//	printer.PushText(maxB);
//	if (maxBs!=0) {
//		printer.OpenElement("sd"); printer.PushText(maxBs); printer.CloseElement();
//	}
//	printer.CloseElement();
//	printer.CloseElement(); // Basal Roots
//
//	printer.CloseElement(); // Seed
//
//	return std::string(printer.CStr());
}

/**
 *
 */
OrganParameter* SeedTypeParameter::realize() const
{
	SeedParameter* sp = new SeedParameter();
	sp->firstB = firstB + randn()*firstBs;
	sp->delayB = delayB + randn()*firstBs;
	sp->maxB = std::round(maxB + randn()*maxBs);
	sp->firstB = std::max(0.,firstB);
	sp->delayB = std::max(0.,delayB);
	sp->maxB = std::max(0,maxB);
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



StemTypeParameter::StemTypeParameter()
{
	organType = Organ::ot_stem;
	subType = -1; // means undefined
	tropism = new StemTropismFunction(0,0);
	growth = new StemGrowthFunction();
	se = new SoilProperty();
	sa = new SoilProperty();
	sbp = new SoilProperty();
	set(-1, 0., 0., 10., 0., 1., 0., 0., 0., 1., 0, 0.1, 0., 150./255.,150./255.,50./255., 1, 1. ,0.2, 0.1,
			successor, successorP, 1.22, 0., 1.e9, 0., 1, "undefined");
}

StemTypeParameter::~StemTypeParameter()
{
	delete tropism;
	delete growth;
	delete se;
	delete sa;
	delete sbp;
}

/**
 * todo comment
 */
void StemTypeParameter::set(int type, double lb, double lbs, double la, double las, double ln, double lns, double nob, double nobs,
		double r, double rs, double a, double as,  double colorR, double colorG, double colorB, int tropismT, double tropismN, double tropismS,
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
	this->colorR = colorR;
	this->colorG = colorG;
	this->colorB = colorB;
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

void StemTypeParameter::createTropism(SignedDistanceFunction* geom, SoilProperty* soil)
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
		StemTropismFunction* gt =  new StemGravitropism(tropismN,tropismS);
		StemTropismFunction* ht = new StemPhototropism(tropismN,tropismS,soil);
		t = new CombinedStemTropism(tropismN,tropismS,ht,10.,gt,1.); // does only use the objective functions from gravitropism and hydrotropism
		break;
	}
    case tt_antigravi: {
		t = new StemAntiGravitropism(tropismN,tropismS);
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

void StemTypeParameter::createGrowth()
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
 * (called by Root::Root())
 *
 * minimal ln distance is 1.e-9
 *
 * \return Specific root parameters derived from the root type parameters
 */
OrganParameter* StemTypeParameter::realize() const
{
	// type does not change
	double lb_ = std::max(lb + randn()*lbs,double(0)); // length of basal zone
	double la_ = std::max(la + randn()*las,double(0)); // length of apical zone
	std::vector<double> ln_; // stores the inter-distances
	int nob_ = std::max(round(nob + randn()*nobs),double(0)); // maximal number of branches
	for (int i = 0; i<nob_-1; i++) { // create inter-root distances
		double d =  ln*(1+i) + randn()*lns; //std::max(  );//ln + randn()*lns,1e-9);
		ln_.push_back(d);
	}
	double r_ = std::max(r + randn()*rs,double(0)); // initial elongation
	double a_ = std::max(a + randn()*as,double(0)); // radius
	double theta_ = std::max(theta + randn()*thetas,double(0)); // initial elongation
	double rlt_ = std::max(rlt + randn()*rlts,double(0)); // root life time
	StemParameter* stem_p =  new StemParameter(subType,lb_,la_,ln_,r_,a_,theta_,rlt_);
	return stem_p;


    for (int i = 0; i<nob_-1; i++) { // create inter-root distances
		double d = 1 +2*i; //std::max(  );//ln + randn()*lns,1e-9);
		ln_.push_back(d);
	}

}

/**
 * Choose (dice) lateral type based on root parameter set,
 * (since there can be more than one lateral type)
 *
 * @param pos       spatial position (for coupling to a soil model)
 */
int StemTypeParameter::getLateralType(const Vector3d& pos)
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

void StemTypeParameter::write(std::ostream & os) const
{
	os << "# Root type parameter for " << name << "\n";
	os << "type\t" << subType << "\n" << "name\t" << name << "\n" << "lb\t"<< lb <<"\t"<< lbs << "\n" << "la\t"<< la <<"\t"<< las << "\n"
			<< "ln\t" << ln << "\t" << lns << "\n" << "nob\t"<< nob <<"\t"<< nobs << "\n" << "r\t"<< r <<"\t"<< rs << "\n" <<
			"a\t" << a << "\t" << as << "\n" << "color\t"<< colorR <<"\t"<< colorG << "\t" << colorB << "\n"
			<< "tropism\t"<< tropismT <<"\t"<< tropismN << "\t" << tropismS << "\n" << "dx\t" << dx << "\n" << "successor\t" << successor.size() << "\t";
	for (size_t i=0; i<successor.size(); i++) {
		os << successor[i] << "\t";
	}
	os << "\n" << "successorP\t" << successorP.size() <<"\t";
	for (size_t i=0; i<successorP.size(); i++) {
		os << successorP[i] << "\t";
	}
	os << "\n" << "theta\t" << theta << "\t" << thetas << "\n" << "rlt\t" << rlt << "\t" << rlts << "\n" << "gf\t" << gf << "\n";
	std::cout << "RootTypeParameter::write is deprecated, use RootTypeParameter::writeXML instead\n";
}

void StemTypeParameter::read(std::istream & is)
{
	char ch[256]; // dummy
	is.getline(ch,256);
	std::string s; // dummy
	double k;
	double ks;
	is >> s >> subType >> s >> name >> s >> lb >> lbs >> s >> la >> las >> s >> ln >> lns >> s >> k >> ks;
	is >> s >> r >> rs >> s >> a >> as >> s >> colorR >> colorG >> colorB >> s >> tropismT >> tropismN >> tropismS >> s >> dx;
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
			std::cout << "StemTypeParameter::read() nobs is nan \n";
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
	std::cout << "StemTypeParameter::read is deprecated, use StemTypeParameter::readXML instead\n";
}

void StemTypeParameter::readXML(FILE* fp)
{

}


std::string StemTypeParameter::writeXML(FILE* fp) const
{
tinyxml2::XMLPrinter printer( fp, false, 0 ); //compact mode false, and 0 indent
        printer.OpenElement("organ", false); //compact mode false
        printer.PushAttribute("type","stem");

	    switch (subType) {
	case 1 :  printer.PushAttribute("name","taproot"); printer.PushAttribute("id","id to discuss");// See

	break;
    case 2 :  printer.PushAttribute("name","lateral1"); printer.PushAttribute("id","id to discuss"); // See
	break;
    case 3 :  printer.PushAttribute("name","lateral2"); printer.PushAttribute("id","id to discuss");// See
	break;
	case 4 :  printer.PushAttribute("name","nodal_root"); printer.PushAttribute("id","id to discuss"); // See
	break;
    case 5 :  printer.PushAttribute("name","shoot_borne_root"); printer.PushAttribute("id","id to discuss"); // See
	break;
    }

    printer.OpenElement("parameter");
//    printer.PushComment("Basal zone [cm]");
    printer.PushAttribute("name","lb"); printer.PushAttribute("value",lb); printer.PushAttribute("dev",lbs);  printer.CloseElement();	 	///< Basal zone [cm]
	printer.OpenElement("parameter");
//	printer.PushComment("Apical zone [cm];");
	printer.PushAttribute("name","la"); printer.PushAttribute("value",la); printer.PushAttribute("dev",las);  printer.CloseElement();	///< Apical zone [cm];
	printer.OpenElement("parameter");
//	printer.PushComment("Inter-lateral distance [cm];");
	printer.PushAttribute("name","ln"); printer.PushAttribute("value",ln); printer.PushAttribute("dev",lns);  printer.CloseElement();		///< Inter-lateral distance [cm];
	printer.OpenElement("parameter");
//	printer.PushComment("Number of branches [1];");
	printer.PushAttribute("name","nob"); printer.PushAttribute("value",nob); printer.PushAttribute("dev",nobs);  printer.CloseElement();		///< Standard deviation apical zone [cm];
	printer.OpenElement("parameter");
//	printer.PushComment("Initial growth rate [cm day-1]");
	printer.PushAttribute("name","r"); printer.PushAttribute("value",r); printer.PushAttribute("dev",rs); printer.CloseElement();		///< Initial growth rate [cm day-1]
	printer.OpenElement("parameter");
//	printer.PushComment("Root radius [cm]");
	printer.PushAttribute("name","a"); printer.PushAttribute("value",a); printer.PushAttribute("dev",as);   printer.CloseElement();		///< Root radius [cm]
	printer.OpenElement("parameter");
//	printer.PushComment("Root color (red, green, blue)");
	printer.PushAttribute("colorR",colorR); printer.PushAttribute("colorG",colorG); printer.PushAttribute("colorB",colorB); printer.CloseElement();	///< Root color (red)
	printer.OpenElement("parameter");
//	printer.PushComment("Root tropism parameter (Type, number of trials, mean vale of expected change)");
	printer.PushAttribute("tropismT",tropismT); printer.PushAttribute("tropismN",tropismN); printer.PushAttribute("tropismS",tropismS);  printer.CloseElement();	///< Root tropism parameter (Type)
	printer.OpenElement("parameter");
//	printer.PushComment("Maximal segment size");
	printer.PushAttribute("name","dx"); printer.PushAttribute("value",dx);  printer.CloseElement();		///< Maximal segment size [cm]
	printer.OpenElement("parameter");
//	printer.PushComment("Angle between root and parent root");
	printer.PushAttribute("name","theta"); printer.PushAttribute("value",theta); printer.PushAttribute("dev",thetas); printer.CloseElement();	///< Angle between root and parent root (rad)
	printer.OpenElement("parameter");
//	printer.PushComment("Root life time");
	printer.PushAttribute("name","rlt"); printer.PushAttribute("value",rlt); printer.PushAttribute("dev",rlts);  printer.CloseElement();	///< Root life time (days)
	printer.OpenElement("parameter");
//	printer.PushComment("Growth function");
	printer.PushAttribute("name","dx"); printer.PushAttribute("value",dx); printer.CloseElement();		///< Growth function (1=negative exponential, 2=linear)


	for (int successorCount = 0; successorCount < successor.size(); successorCount ++ ){
            std::string s = std::to_string(successorCount);
            char const *schar = s.c_str();
        printer.OpenElement("successor"); printer.PushAttribute("name",schar); printer.PushAttribute("location",successor[successorCount]);    printer.CloseElement();
        printer.OpenElement("successorP"); printer.PushAttribute("name",schar); printer.PushAttribute("percentage",successorP[successorCount]);   printer.CloseElement();
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



LeafTypeParameter::LeafTypeParameter()
{
	organType = Organ::ot_leafe;
	subType = -1; // means undefined
	tropism = new LeafTropismFunction(0,0);
	growth = new LeafGrowthFunction();
	se = new SoilProperty();
	sa = new SoilProperty();
	sbp = new SoilProperty();
	set(-1, 0., 0., 10., 0., 1., 0., 0., 0., 1., 0, 0.1, 0., 150./255.,150./255.,50./255., 1, 1. ,0.2, 0.1,
			successor, successorP, 1.22, 0., 1.e9, 0., 1, "undefined");
}

LeafTypeParameter::~LeafTypeParameter()
{
	delete tropism;
	delete growth;
	delete se;
	delete sa;
	delete sbp;
}

/**
 * todo comment
 */
void LeafTypeParameter::set(int type, double lb, double lbs, double la, double las, double ln, double lns, double nob, double nobs,
		double r, double rs, double a, double as,  double colorR, double colorG, double colorB, int tropismT, double tropismN, double tropismS,
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
	this->colorR = colorR;
	this->colorG = colorG;
	this->colorB = colorB;
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

void LeafTypeParameter::createTropism(SignedDistanceFunction* geom, SoilProperty* soil)
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

void LeafTypeParameter::createGrowth()
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
OrganParameter* LeafTypeParameter::realize() const
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
	LeafParameter* leaf_p =  new LeafParameter(subType,lb_,la_,ln_,r_,a_,theta_,rlt_);
	return leaf_p;
}

/**
 * Choose (dice) lateral type based on root parameter set,
 * (since there can be more than one lateral type)
 *
 * @param pos       spatial position (for coupling to a soil model)
 */
int LeafTypeParameter::getLateralType(const Vector3d& pos)
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

void LeafTypeParameter::write(std::ostream & os) const
{
	os << "# Root type parameter for " << name << "\n";
	os << "type\t" << subType << "\n" << "name\t" << name << "\n" << "lb\t"<< lb <<"\t"<< lbs << "\n" << "la\t"<< la <<"\t"<< las << "\n"
			<< "ln\t" << ln << "\t" << lns << "\n" << "nob\t"<< nob <<"\t"<< nobs << "\n" << "r\t"<< r <<"\t"<< rs << "\n" <<
			"a\t" << a << "\t" << as << "\n" << "color\t"<< colorR <<"\t"<< colorG << "\t" << colorB << "\n"
			<< "tropism\t"<< tropismT <<"\t"<< tropismN << "\t" << tropismS << "\n" << "dx\t" << dx << "\n" << "successor\t" << successor.size() << "\t";
	for (size_t i=0; i<successor.size(); i++) {
		os << successor[i] << "\t";
	}
	os << "\n" << "successorP\t" << successorP.size() <<"\t";
	for (size_t i=0; i<successorP.size(); i++) {
		os << successorP[i] << "\t";
	}
	os << "\n" << "theta\t" << theta << "\t" << thetas << "\n" << "rlt\t" << rlt << "\t" << rlts << "\n" << "gf\t" << gf << "\n";
	std::cout << "RootTypeParameter::write is deprecated, use RootTypeParameter::writeXML instead\n";
}

void LeafTypeParameter::read(std::istream & is)
{
	char ch[256]; // dummy
	is.getline(ch,256);
	std::string s; // dummy
	double k;
	double ks;
	is >> s >> subType >> s >> name >> s >> lb >> lbs >> s >> la >> las >> s >> ln >> lns >> s >> k >> ks;
	is >> s >> r >> rs >> s >> a >> as >> s >> colorR >> colorG >> colorB >> s >> tropismT >> tropismN >> tropismS >> s >> dx;
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
			std::cout << "LeafTypeParameter::read() nobs is nan \n";
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
	std::cout << "LeafTypeParameter::read is deprecated, use LeafTypeParameter::readXML instead\n";
}

void LeafTypeParameter::readXML(FILE* fp)
{

}


std::string LeafTypeParameter::writeXML(FILE* fp) const
{
tinyxml2::XMLPrinter printer( fp, false, 0 ); //compact mode false, and 0 indent
        printer.OpenElement("organ", false); //compact mode false
        printer.PushAttribute("type","leaf");

	    switch (subType) {
	case 1 :  printer.PushAttribute("name","taproot"); printer.PushAttribute("id","id to discuss");// See

	break;
    case 2 :  printer.PushAttribute("name","lateral1"); printer.PushAttribute("id","id to discuss"); // See
	break;
    case 3 :  printer.PushAttribute("name","lateral2"); printer.PushAttribute("id","id to discuss");// See
	break;
	case 4 :  printer.PushAttribute("name","nodal_root"); printer.PushAttribute("id","id to discuss"); // See
	break;
    case 5 :  printer.PushAttribute("name","shoot_borne_root"); printer.PushAttribute("id","id to discuss"); // See
	break;
    }

    printer.OpenElement("parameter");
//    printer.PushComment("Basal zone [cm]");
    printer.PushAttribute("name","lb"); printer.PushAttribute("value",lb); printer.PushAttribute("dev",lbs);  printer.CloseElement();	 	///< Basal zone [cm]
	printer.OpenElement("parameter");
//	printer.PushComment("Apical zone [cm];");
	printer.PushAttribute("name","la"); printer.PushAttribute("value",la); printer.PushAttribute("dev",las);  printer.CloseElement();	///< Apical zone [cm];
	printer.OpenElement("parameter");
//	printer.PushComment("Inter-lateral distance [cm];");
	printer.PushAttribute("name","ln"); printer.PushAttribute("value",ln); printer.PushAttribute("dev",lns);  printer.CloseElement();		///< Inter-lateral distance [cm];
	printer.OpenElement("parameter");
//	printer.PushComment("Number of branches [1];");
	printer.PushAttribute("name","nob"); printer.PushAttribute("value",nob); printer.PushAttribute("dev",nobs);  printer.CloseElement();		///< Standard deviation apical zone [cm];
	printer.OpenElement("parameter");
//	printer.PushComment("Initial growth rate [cm day-1]");
	printer.PushAttribute("name","r"); printer.PushAttribute("value",r); printer.PushAttribute("dev",rs); printer.CloseElement();		///< Initial growth rate [cm day-1]
	printer.OpenElement("parameter");
//	printer.PushComment("Root radius [cm]");
	printer.PushAttribute("name","a"); printer.PushAttribute("value",a); printer.PushAttribute("dev",as);   printer.CloseElement();		///< Root radius [cm]
	printer.OpenElement("parameter");
//	printer.PushComment("Root color (red, green, blue)");
	printer.PushAttribute("colorR",colorR); printer.PushAttribute("colorG",colorG); printer.PushAttribute("colorB",colorB); printer.CloseElement();	///< Root color (red)
	printer.OpenElement("parameter");
//	printer.PushComment("Root tropism parameter (Type, number of trials, mean vale of expected change)");
	printer.PushAttribute("tropismT",tropismT); printer.PushAttribute("tropismN",tropismN); printer.PushAttribute("tropismS",tropismS);  printer.CloseElement();	///< Root tropism parameter (Type)
	printer.OpenElement("parameter");
//	printer.PushComment("Maximal segment size");
	printer.PushAttribute("name","dx"); printer.PushAttribute("value",dx);  printer.CloseElement();		///< Maximal segment size [cm]
	printer.OpenElement("parameter");
//	printer.PushComment("Angle between root and parent root");
	printer.PushAttribute("name","theta"); printer.PushAttribute("value",theta); printer.PushAttribute("dev",thetas); printer.CloseElement();	///< Angle between root and parent root (rad)
	printer.OpenElement("parameter");
//	printer.PushComment("Root life time");
	printer.PushAttribute("name","rlt"); printer.PushAttribute("value",rlt); printer.PushAttribute("dev",rlts);  printer.CloseElement();	///< Root life time (days)
	printer.OpenElement("parameter");
//	printer.PushComment("Growth function");
	printer.PushAttribute("name","dx"); printer.PushAttribute("value",dx); printer.CloseElement();		///< Growth function (1=negative exponential, 2=linear)


	for (int successorCount = 0; successorCount < successor.size(); successorCount ++ ){
            std::string s = std::to_string(successorCount);
            char const *schar = s.c_str();
        printer.OpenElement("successor"); printer.PushAttribute("name",schar); printer.PushAttribute("location",successor[successorCount]);    printer.CloseElement();
        printer.OpenElement("successorP"); printer.PushAttribute("name",schar); printer.PushAttribute("percentage",successorP[successorCount]);   printer.CloseElement();
//
    }
//    printer.PrintSpace();
    printer.CloseElement(false); // close element compact mode false

	return std::string(printer.CStr());


}

