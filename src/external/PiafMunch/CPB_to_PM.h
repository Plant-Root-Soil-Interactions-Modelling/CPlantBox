#ifndef CPB2PM_H_
#define CPB2PM_H_


#include "MappedOrganism.h"
//#include "runPM.h"

class PhloemFlux;

class CPB_to_PM: public std::enable_shared_from_this<CPB_to_PM>
{

	public:
	
	CPB_to_PM(std::shared_ptr<CPlantBox::MappedPlant> plant, 
				std::shared_ptr<PhloemFlux> piafmunch): plant(plant), piafmunch(piafmunch){};
    std::shared_ptr<CPB_to_PM> CPB_to_PM_ptr() { return shared_from_this(); }; // up-cast for Python binding
	virtual ~CPB_to_PM() { }
	
	void waterLimitedGrowth(std::shared_ptr<CPlantBox::Organ> org, 
							std::vector<int> localGrowingNodesId, double t, double deltaVol);
	double adaptDt(std::shared_ptr<CPlantBox::Organ> org, double dt);
	void assertUsedCReserves(std::shared_ptr<CPlantBox::Organ> org);
	double getMaxVolumicGrowth(std::shared_ptr<CPlantBox::Organ> org, double t, double rmax, double Lmax);
	std::vector<int> getGrowingNodes(std::shared_ptr<CPlantBox::Organ> org);
	void computeFpsi();
	void computeFlen(std::shared_ptr<CPlantBox::Organ> org, std::vector<int> growingNodesId);
	void updateBudStage(double EndTime);
	void organToNodeData(double dt);
	void getRmaxLmax_st_f(std::shared_ptr<CPlantBox::Organ> org, int st, int ot, double & rmax, double & Lmax);
	void setAuxinSource(std::shared_ptr<CPlantBox::Organ> org);
	void computeAuxinValue(std::shared_ptr<CPlantBox::Organ> org);
	
	std::shared_ptr<CPlantBox::MappedPlant> plant;	
	std::shared_ptr<PhloemFlux> piafmunch;	
	
	std::map<int, double> BackUpMaxGrowth;//to check at runtime if growth is correct. is a map because max(orgID) > len(orgs)
	
    double computeBerth(double ss, double aa){return (aa+Berthkaa)/(ss+Berthkss); };
		
	//		from plant shape	
	void setKr_st(std::vector<std::vector<double>> values, double kr_length_, bool verbose = false); ///< sets a callback for kr_suc:=kr_suc(ot,type), 
	void setKx_st(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type), 
	void setRmax_st(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type), 
	void setAcross_st(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type), 
	void setPerimeter_st(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type),  
	void setRhoSucrose(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type),  
	void setKrm1(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type), 
	void setKrm2(std::vector<std::vector<double>> values, bool verbose = false); ///< sets a callback for kx_suc:=kx_suc(ot,type), 
	
    std::function<double(int, int, int)> kr_st_f = []( int type, int orgtype, int si) {
		throw std::runtime_error("kr_st_f not implemented"); 
		return 0.; };
    std::function<double(int,int)> kx_st_f = [](int type, int orgtype) {
		throw std::runtime_error("kx_st_f not implemented"); 
		return 0.; };
    std::function<double(int, int)> Across_st_f = [](int type, int orgtype) {//cross-sectional area of all the sieve tubes in segment
		throw std::runtime_error("get_Across_st not implemented"); 
		return 0.; };
    std::function<double(int, int)> Perimeter_st_f = [](int type, int orgtype) {//cross-sectional area of all the sieve tubes in segment
		throw std::runtime_error("get_Across_st not implemented"); 
		return 0.; };
    std::function<double(int,int)> Rmax_st_f = []( int type, int orgtype){//maximum initial growth rate use by phloem module
		throw std::runtime_error("get_Rmax not implemented"); 
		return 0.; };
    std::function<double(int,int)> rhoSucrose_f = []( int type, int orgtype){//maximum initial growth rate use by phloem module
		throw std::runtime_error("get_rhoSucrose not implemented"); 
		return 0.; };
    std::function<double(int,int)> krm1_f = []( int type, int orgtype){//maximum initial growth rate use by phloem module
		throw std::runtime_error("get_rhoSucrose not implemented"); 
		return 0.; };
    std::function<double(int,int)> krm2_f = []( int type, int orgtype){//maximum initial growth rate use by phloem module
		throw std::runtime_error("get_rhoSucrose not implemented"); 
		return 0.; };

	
	// per type
    std::vector<std::vector<double>> kr_st;//  [mmol hPa-1 day-1]
	std::vector<std::vector<double>> kx_st; //  [cm3 hPa-1 day-1]
    std::vector<std::vector<double>> Across_st; // [cm2]
    std::vector<std::vector<double>> Perimeter_st; // [cm]
	std::vector<std::vector<double>> Rmax_st; // [cm day-1]
	std::vector<std::vector<double>> rhoSucrose;
	std::vector<std::vector<double>> krm1v; 
	std::vector<std::vector<double>> krm2v;
	
	
    //for growth
	double psi_osmo_proto = -4000*1.0197;
    double psiMin = 2000*1.0197;//limit wat. pot. in xylem for water-limited growth, [cm]
    double leafGrowthZone = 2;//cm
    bool StemGrowthPerPhytomer = true;
	std::vector<double> Fpsi;//water scarcity factor for growth, (-)
    std::vector<double> Flen;
    std::vector<int> GrowthZone;
    std::vector<int> GrowthZoneLat;
	std::vector<std::map<int,double>> deltaSucOrgNode;//maximal sucrose need for growth per node, (mmol Suc d-1)
    std::vector<double> psi_p_symplasm;
	bool do_gf_warning = true; // get a warning in case we have an unexpected growth function (gives warning only once)
		
	// auxin
	double CSTthreshold = 0.3;
	std::vector<double> AuxinSource;
    double L_dead_threshold = 2.;
    double limLenActive;
    double BerthLim = -1;
	double Berthkss = 0.2;
	double Berthkaa = 1.;
	
	//plant data
	int Nt; //number of nodes
	//std::vector<int> budTypes = {2,3};
	std::vector<double> maxLBudDormant;
	std::vector<double> maxLBud;
	double budGR;
	
	protected:
	
	
	
	//retrieve tissue-specific parameters
	double kr_st_const(  int type, int organType, int si) { return kr_st.at(0).at(0); }  //constant
    double kr_st_perOrgType( int type, int organType, int si) { return kr_st.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double kr_st_perType( int type, int organType, int si) {return kr_st.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	double kr_st_RootExchangeZonePerType(int type, int organType, int si)
	{ 
		if (organType == CPlantBox::Organism::ot_root){
			double coef = plant->exchangeZoneCoefs.at(si);//% of segment length in the root exchange zone, see MappedPlant::simulate
			return coef * kr_st.at(organType - 2).at(type); 
		}
		return kr_st.at(organType - 2).at(type);  
	} //subtype, type and depend on distance to tip for roots
	
	double kx_st_const( int type, int organType) { return kx_st.at(0).at(0); }  //constant
    double kx_st_perOrgType( int type, int organType) { return kx_st.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double kx_st_perType( int type, int organType) {return kx_st.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	
	double Across_st_const( int type, int organType) { return Across_st.at(0).at(0); }  //constant
    double Across_st_perOrgType( int type, int organType) { return Across_st.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double Across_st_perType( int type, int organType) {return Across_st.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
    
	double Perimeter_st_const( int type, int organType) { return Perimeter_st.at(0).at(0); }  //constant
    double Perimeter_st_perOrgType( int type, int organType) { return Perimeter_st.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double Perimeter_st_perType( int type, int organType) {return Perimeter_st.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	
	double Rmax_st_const( int type, int organType) { return Rmax_st.at(0).at(0); }  //constant
    double Rmax_st_perOrgType( int type, int organType) { return Rmax_st.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double Rmax_st_perType( int type, int organType) {return Rmax_st.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	
	double rhoSucrose_const( int type, int organType) { return rhoSucrose.at(0).at(0); }  //constant
    double rhoSucrose_perOrgType( int type, int organType) { return rhoSucrose.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double rhoSucrose_perType( int type, int organType) {return rhoSucrose.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	
	double krm1_const( int type, int organType) { return krm1v.at(0).at(0); }  //constant
    double krm1_perOrgType( int type, int organType) { return krm1v.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double krm1_perType( int type, int organType) {return krm1v.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	
	double krm2_const( int type, int organType) { return krm2v.at(0).at(0); }  //constant
    double krm2_perOrgType( int type, int organType) { return krm2v.at(organType - 2).at(0); } //per organ type (goes from 2 (root) to 4 (leaf))
    double krm2_perType( int type, int organType) {return krm2v.at(organType - 2).at(type); }//per subtype and organ type (goes from 2 (root) to 4 (leaf))
	
	
};
#endif