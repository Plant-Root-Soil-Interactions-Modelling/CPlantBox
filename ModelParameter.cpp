#include "ModelParameter.h"
#include <cmath>

/*
 * RootParameter: parameters for a root type
 */

/**
 * Default constructor
 */
RootTypeParameter::RootTypeParameter() {
	set(-1, 0., 0., 10., 0., 1., 0., 0., 0., 1., 0, 0.1, 0., 150./255.,150./255.,50./255., 1, 1. ,0.2, 0.1,
			successor, successorP, 1.22, 0., 1.e9, 0., 1, "undefined");
}

/**
 * Copy constructor
 */
RootTypeParameter::RootTypeParameter(const RootTypeParameter& rp) :type(rp.type), lb(rp.lb), lbs(rp.lbs), la(rp.la), las(rp.las),
		ln(rp.ln), lns(rp.lns), nob(rp.nob), nobs(rp.nobs), r(rp.r), rs(rp.rs), a(rp.a), as(rp.as), colorR(rp.colorR), colorG(rp.colorG),
		colorB(rp.colorB), tropismT(rp.tropismT), tropismN(rp.tropismN), tropismS(rp.tropismS), dx(rp.dx), theta(rp.theta), thetas(rp.thetas),
		rlt(rp.rlt), rlts(rp.rlts), gf(rp.gf), name(rp.name), successor(rp.successor),
		successorP(rp.successorP) { }

void RootTypeParameter::set(int type, double lb, double lbs, double la, double las, double ln, double lns, double nob, double nobs,
		double r, double rs, double a, double as,  double colorR, double colorG, double colorB, double tropismT, double tropismN, double tropismS,
		double dx, const std::vector<int>& successor, const std::vector<double>& successorP, double theta, double thetas, double rlt, double rlts,
		int gf, const std::string& name) {
	this->type = type;
	this->lb = lb;		this->lbs = lbs;
	this->la = la;		this->las = las;
	this->ln = ln;		this->lns = lns;
	this->nob = nob;	this->nobs = nobs;
	this->r = r;		this->rs = rs;
	this->a = a;		this->as = as;
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
	this->theta = theta;	this->thetas = thetas;
	this->rlt = rlt;		this->rlts = rlts;
	this->gf = gf;
	this->name = name; // string
}

/**
 * Creates a specific root from the root type parameters.
 * The unique root id is not set, but must be set from outside.
 * (called by Root::Root())
 *
 * minimal ln distance is 1.e-9 per cut off
 *
 * \return Specific root parameters derived from the root type parameters
 */
RootParameter RootTypeParameter::realize() {
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
	RootParameter p(type,lb_,la_,ln_,nob_,r_,a_,theta_,rlt_);
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

void RootTypeParameter::write(std::ostream & cout) const {
	cout << "# Root type parameter for " << name << "\n";
	cout << "type\t" << type << "\n" << "name\t" << name << "\n" << "lb\t"<< lb <<"\t"<< lbs << "\n" << "la\t"<< la <<"\t"<< las << "\n"
			<< "ln\t" << ln << "\t" << lns << "\n" << "nob\t"<< nob <<"\t"<< nobs << "\n" << "r\t"<< r <<"\t"<< rs << "\n" <<
			"a\t" << a << "\t" << as << "\n" << "color\t"<< colorR <<"\t"<< colorG << "\t" << colorB << "\n"
			<< "tropism\t"<< tropismT <<"\t"<< tropismN << "\t" << tropismS << "\n" << "dx\t" << dx << "\n" << "successor\t" << successor.size() << "\t";
	for (size_t i=0; i<successor.size(); i++) {
		cout << successor[i] << "\t";
	}
	cout << "\n" << "successorP\t" << successorP.size() <<"\t";
	for (size_t i=0; i<successorP.size(); i++) {
		cout << successorP[i] << "\t";
	}
	cout << "\n" << "theta\t" << theta << "\t" << thetas << "\n" << "rlt\t" << rlt << "\t" << rlts << "\n" << "gf\t" << gf << "\n";
}

void RootTypeParameter::read(std::istream & cin) {
	char ch[256]; // dummy
	cin.getline(ch,256);
	std::string s; // dummy
	double k;
	double ks;
	cin >> s >> type >> s >> name >> s >> lb >> lbs >> s >> la >> las >> s >> ln >> lns >> s >> k >> ks;
	cin >> s >> r >> rs >> s >> a >> as >> s >> colorR >> colorG >> colorB >> s >> tropismT >> tropismN >> tropismS >> s >> dx;
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
	cin  >> s >> n;
	successor.clear();
	int is;
	for (int i=0; i<n; i++) {
		cin >> is;
		successor.push_back(is);
	}
	cin >> s >> n;
	successorP.clear();
	double ds;
	for (int i=0; i<n; i++) {
		cin >> ds;
		successorP.push_back(ds);
	}
	cin >> s >> theta >> thetas >> s >> rlt >> rlts >> s >> gf >> s;
	//successorP.push_back(1); // to be on the safe side, in case P does not sum up to 1
	//successor.push_back(-1); // -1 means no lateral
}

/*
 * class RootParameter
 */

/**
 * Default constructor
 */
RootParameter:: RootParameter() {
	std::vector<double> ln;
	set(-1,0.,0.,ln, 0,0.,0.,0.,0.);
}

double RootParameter::getK() const {
	double l = 0;
	for (auto const& dl : ln) {
		l += dl;
	}
	return l+la+lb;
}

void RootParameter::set(int type, double lb, double la, const std::vector<double>& ln, double nob, double r, double a, double theta, double rlt) {
	this->type=type;
	this->lb=lb;
	this->la=la;
	this->ln=ln;
	this->nob=nob;
	this->r=r;
	this->a=a;
	this->theta=theta;
	this->rlt=rlt;
}

void RootParameter::write(std::ostream & cout) const {
	cout << "# Root Parameters \n";
	cout << "type\t" << type << "\n"  << "lb\t"<< lb <<"\n" << "la\t"<< la <<"\n" << "ln\t";
	for (size_t i=0; i<ln.size(); i++) {
		cout << ln[i] << "\t";
	}
	cout << "\n" << "nob\t"<< nob <<"\n" << "r\t"<< r <<"\n" << "a\t" << a << "\n" << "theta\t" << theta << "\n" << "rlt\t" << rlt << "\n";
}



/*
 * class RootSystemParameter
 */

/**
 * Default constructor: No basal roots, not shoot borne roots, planting depth 3 [cm]
 */
RootSystemParameter::RootSystemParameter() {
	set(3.,1.e9,0.,0, //pd, fB, dB, mB,
			0,1.e9,1.e9,0.,0.,30.);  // nC, fSB, dSB, dRC, nz
}


RootSystemParameter::~RootSystemParameter() {
}


void RootSystemParameter::read(std::istream & cin) {
	double plantingdepth;
	std::string s; // dummy
	cin  >>  s >> plantingdepth;
	cin >> s >> firstB >> s >> delayB >> s >> maxB >> s >> nC >> s >> firstSB >> s >> delaySB >> s >> delayRC >> s >> nz >> s >> simtime;
	seedPos = Vector3d(0,0,-plantingdepth);
}


void RootSystemParameter::write(std::ostream & cout) const {
	double pd = -seedPos.z;
	cout <<  "plantingdepth\t" << pd << "\n" <<  "firstB\t" << firstB << "\n" <<  "delayB\t" << delayB << "\n"
			<<  "maxB\t" << maxB << "\n" <<  "nC\t" << nC << "\n" <<  "firstSB\t" << firstSB << "\n"
			<<  "delaySB\t" << delaySB << "\n" <<  "delayRC\t" << delayRC << "\n" <<  "nz\t" << nz << "\n" << "simulationTime\t" << simtime << "\n";
}

/**
 * Set all plant parameters
 *
 * @param pd          Planting depth [cm]
 * @param fB          Emergence of first basal root [day]
 * @param dB          Time delay between the basal roots [day]
 * @param mB          Maximal number of basal roots [1]
 * @param nC          Maximal number of roots per root crown [1]
 * @param fSB         First emergence of a shoot borne root [day]
 * @param dSB         Time delay between the shoot borne roots [day]
 * @param dRC         Delay between the root crowns [day]
 * @param nz          Distance between the root crowns along the shoot [cm]
 * @param simtime 	  Recommended final simulation time (e.g. used in the web interface)
 */
void RootSystemParameter::set(double pd, double fB, double dB, int mB, int nC, double fSB, double dSB, double dRC, double nz, double simtime) {
	seedPos=Vector3d(0,0,-pd);
	firstB=fB;
	delayB=dB;
	maxB=mB;
	this->nC=nC;
	firstSB=fSB;
	delaySB=dSB;
	delayRC=dRC;
	this->nz=nz;
	this->simtime=simtime;
}
