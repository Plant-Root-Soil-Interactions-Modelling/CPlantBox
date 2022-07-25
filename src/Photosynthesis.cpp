// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

#include "Photosynthesis.h"
//#include <armadillo>
//#include <algorithm>
#include <set>
#include <external/Eigen/Dense>
#include <external/Eigen/Sparse> 
#include <iostream>
#include <fstream>
																							
																		

namespace CPlantBox {

/**
 * A "photosyntehsis" object, as needed for water flux computations + sucrose assimilation + stomatale opening
 *
 * @param plant_     	MappedPlant object
 * @param psiXylInit   	Initial guess for the value of xylem wat. pot [cm]
 * @param ciInit      	Initial guess for intracellular CO2 partial pressure [mol mol-1]
 */
Photosynthesis::Photosynthesis(std::shared_ptr<CPlantBox::MappedPlant> plant_, double psiXylInit, double ciInit): 
	XylemFlux(std::shared_ptr<CPlantBox::MappedSegments>(plant_)), plant(plant_), psiXylInit(psiXylInit), ciInit(ciInit)
{//check when plant and planphotosyn are diff
	//std::cout<<"alive creation "<<std::endl;
	//this->seg_leaves_idx = plant->getNodeIds(4);//ids of leaf segments
	this->seg_leaves_idx = plant->getSegmentIds(4);//ids of leaf segments
	psiXyl.resize(plant->nodes.size(), psiXylInit);//-500
	An.resize(seg_leaves_idx.size(), 0.);
	gco2.resize(seg_leaves_idx.size(), 0.);
	ci.resize(seg_leaves_idx.size(), ciInit);//350e-6
	pg.resize(seg_leaves_idx.size(), psiXylInit);
											  
	k_stomatas.resize(seg_leaves_idx.size(), 0.);
										  
	//std::cout<<"segleafnum "<<seg_leaves_idx.size()<<std::endl;
}


	/* solves the coupled water flux + carbon assimilation and stomatal opening, 
		@param sim_time_ [day]           simulation time
		@param sxx_ [cm]                 soil matric potentials given per segment or per soil cell
		@param cells_                    indicates if the matric potentials are given per cell (True) or by segments (False)
		@param soil_k [day-1]           optionally, soil conductivities can be prescribed per segment, 
										conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)
		@param doLog_                    indicates if computed values should be printed in a text file (True) or not (False) 
		@param verbose_                  print at runtime nothing (0), sparsly (1), many outputs (2) 
		@param RH_ [%]                	relative ari humidity
		@param TairC_ [°C]               leaf temperature (mean)     
	*/
	
void Photosynthesis::solve_photosynthesis(double sim_time_,std::vector<double> sxx_, bool cells_ , 
	 std::vector<double> soil_k_, bool doLog_ , int verbose_ , double RH_, double TairC_)
{
	//		save Environmental and other input variables
	doLog = doLog_; verbose_photosynthesis = verbose_;	
					   
	loop = 0;        
	this->stop = false;
	this->seg_leaves_idx = plant->getSegmentIds(4);//ids of leaf segments
	this->TairC = TairC_;
	this->TleafK = TairC_ + 273.15;
	this->es = 0.61078 * std::exp(17.27 * this->TairC / (this->TairC + 237.3)) /10;//hPa
	this->ea = this->es * RH_ ;//hPa
	//cm = log(-) * (mg cm-3) * (hPa cm3 K−1 mmol−1) * K * (1/[mg mmol-1]) * (cm/hPa)
	this->psi_air = std::log(RH_) * rho_h2o * R_ph * (this->TairC + 237.3)/Mh2o * (1/0.9806806)  ; //in cm
	assert(((plant->kr_length < 0)||(plant->exchangeZoneCoefs.size()==plant->segments.size()))&&"(plant->exchangeZoneCoefs.size()==plant->segments.size()) while kr_length >0");
	//		creat first guesses arrays + "old" values
	psiXyl= std::vector<double>( plant->nodes.size(), psiXylInit);//-500
	An= std::vector<double>( seg_leaves_idx.size(), 0.);
	gco2= std::vector<double>( seg_leaves_idx.size(), 0.);
	ci= std::vector<double>( seg_leaves_idx.size(), ciInit);
	pg= std::vector<double>( seg_leaves_idx.size(), 0.);
															 
	k_stomatas= std::vector<double>( seg_leaves_idx.size(), 0.);
														 
	orgsVec = plant->getOrgans(-1);
 
	An_old = An; gco2_old = gco2; ci_old = ci;
	outputFlux_old.resize(psiXyl.size(), 0.);
	outputFlux.resize(plant->segments.size(), 0.);
	psiXyl_old = psiXyl; pg_old = pg; k_stomatas_old = k_stomatas;
	
	k_stomatas.clear();//to not take k_stomatas into account in @Photosynthesis::initCalcs
	assert(k_stomatas.empty() &&"Photosynthesis::initStruct: k_stomatas not empty");
	
	//		compute parameters which do not change between the loops
	initCalcs(sim_time_);
	k_stomatas = k_stomatas_old ;											
	while(!this->stop){  
		std::fill(maxErr.begin(), maxErr.end(), 0.);//re-initialize error vector
		loopCalcs(sim_time_) ;//compute photosynthesis outputs
		if((verbose_photosynthesis > 1)){std::cout<<"to linearSystem"<<std::endl;}
		linearSystemSolve(sim_time_, sxx_, cells_, soil_k_); //compute psiXyl
		if((verbose_photosynthesis > 1)){std::cout<<"to outputFlux"<<std::endl;}
		//usefull only to know whether we reached convergence
		outputFlux = segFluxes(sim_time_, this->psiXyl, sxx_, false, cells_, std::vector<double>());//approx = false
		if((verbose_photosynthesis > 1)){std::cout<<"to getError"<<std::endl;}
		getError(sim_time_);
		
		loop++ ;
		this->stop = ((loop > maxLoop) || ((maxMaxErr < limMaxErr)&&(loop>minLoop)));//reached convergence or max limit of loops?
		
		if(!this->stop){
			outputFlux_old = outputFlux;
			k_stomatas_old = k_stomatas; ci_old = this->ci; pg_old = this->pg;
			psiXyl_old = this->psiXyl; An_old = this->An; gco2_old = this->gco2;
		}
		
		if((verbose_photosynthesis > 0)){
			std::cout<<"leuning computation module at "<<(loop-1)<<" trials. "
			"Sum of max relative error calculated at the "
			"last two trials: "<<maxMaxErr<<std::endl;
			std::cout<<"each val "<<maxErr[0]<<" "<<maxErr[1]<<" "<<maxErr[2]<<" "<<maxErr[3]<<" "<<maxErr[4];
			std::cout<<" "<<maxErr[5]<<" "<<maxErr[6]<<" "<<maxErr[7]<<" "<<maxErr[8]<<std::endl;
		}
		
	}
	
		
	outputFlux = segFluxes(sim_time_, this->psiXyl, sxx_, false, cells_, std::vector<double>());//approx = false
	loop++ ;
	
	// for phloem flow
	getAg4Phloem();
	doAddGravity(); 
	
	
	if((verbose_photosynthesis  > 0))
	{
		std::cout<<"leuning computation module stopped after "<<(loop-1)<<" trials. "
		"Sum of max relative error calculated at the "
		"last two trials: "<<maxMaxErr<<std::endl;
		std::cout<<"each val "<<maxErr[0]<<" "<<maxErr[1]<<" "<<maxErr[2]<<" "<<maxErr[3]<<" "<<maxErr[4];
		std::cout<<" "<<maxErr[5]<<" "<<maxErr[6]<<" "<<maxErr[7]<<" "<<maxErr[8]<<std::endl;
	}
}


/**
 * Solves the linear system filled by @see XylemFlux::linearSystem
 *
 * @param simTime[day]  	current simulation time, needed for age dependent conductivities,
 *                  		to calculate the age from the creation times (age = sim_time - segment creation time).
 * @param sx [cm]			soil matric potential in the cells or around the segments, given per cell or per segment
 * @param cells 			sx per cell (true), or segments (false)
 * @param soil_k [day-1]    optionally, soil conductivities can be prescribed per segment,
 *                          conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)
 */
void Photosynthesis::linearSystemSolve(double simTime_, const std::vector<double>& sxx_, bool cells_, const std::vector<double> soil_k_)
{
	//get "tripletList" and "b"
	bool withEigen_ = true;
	linearSystem(simTime_, sxx_, cells_, soil_k_, withEigen_); //see XylemFlux::linearSystem
    int N = rs->nodes.size(); // number of nodes
	Eigen::SparseMatrix<double> mat(N,N);
	mat.reserve(Eigen::VectorXi::Constant(N,2));
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	mat.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
	lu.compute(mat);
	
	if(lu.info() != Eigen::Success){
		std::cout << "XylemFlux::linearSystem  matrix Compute with Eigen failed: " << lu.info() << std::endl;
		assert(false);
	}
	
	Eigen::VectorXd v2;
	try{ 
		v2= lu.solve(b);
	}catch(...){
		assert(false&&"XylemFlux::linearSystem error when solving wat. pot. xylem with Eigen ");
	}
	std::vector<double> v3(&v2[0], v2.data()+v2.cols()*v2.rows());
	psiXyl = v3;
}

/* 
		Computes error % for each segment for each of the variables of interestes.
		@param sim_time_ [day]           simulation time, needed when doLog = true
*/

void Photosynthesis::getError(double simTime)
{
	std::ofstream myfile1;
	if(doLog){
		std::string name1 = "errphoto";
		std::string name2 = ".txt";
		std::string namefile = name1 +"_"+ std::to_string(simTime)+ "_"+std::to_string(loop) + name2;
		myfile1.open (namefile);
		myfile1 <<"psi err "<<std::endl;
		myfile1 <<std::endl;}
	
	for (int i = 0; i < this->psiXyl.size(); i++) 
	{ 
		// maxErr[0] =std::max(std::abs((this->psiXyl[i]-psiXyl_old[i])/this->psiXyl[i]),maxErr[0]);
		// if(doLog){
			// myfile1 << "maxErr[1] "<<i <<", "<<plant->organTypes[i]<<", "<<plant->subTypes[i]<<", "<<this->psiXyl[i]<<" "<<psiXyl_old[i];
		// }
		assert(this->psiXyl[i]==this->psiXyl[i] && "Photosynthesis xylnew is nan");
		assert(!std::isnan(this->psiXyl[i]) && "Photosynthesis xylcurrent is nan" );
		assert(!std::isnan(psiXyl_old[i]) && "Photosynthesis xylold is nan" );
		//std::cout<<this->psiXylnew[i]<<" "<<(this->psiXylnew[i]<=0)<<" "<<std::flush;
		if(this->psiXyl[i]>0){
			std::cout<< "Photosynthesis xyl >0, "<<i<<" "<<this->psiXyl[i]<<std::endl;
			assert(false);}
			
		double tempVal = 1;
		
		if(this->psiXyl[i] !=0){tempVal= std::min(std::abs(this->psiXyl[i]), std::abs(this->psiXyl_old[i]));}else{tempVal=1.;}
		maxErr[0] =std::max(std::abs((this->psiXyl[i]-psiXyl_old[i])/tempVal),maxErr[0]);
		if( i < (this->psiXyl.size() - 1)){
			if(this->outputFlux[i] !=0){tempVal=  std::min(std::abs(this->outputFlux[i]), std::abs(this->outputFlux_old[i]));
				}else{tempVal=1.;}
			maxErr[5] =std::max(std::abs((this->outputFlux[i]-outputFlux_old[i])/tempVal),maxErr[5]);
			maxErr[7] += std::abs((this->outputFlux[i]-outputFlux_old[i])/tempVal);//becasue sum of fluxes important 
		}
		
		if(doLog){
			myfile1 <<i<<" "<<", maxErr[5] "<<plant->organTypes[std::max(i-1,0)]<<" "<<maxErr[0];
			myfile1<<" "<<psiXyl[i]<<" "<<psiXyl_old[i]<<" "<<i<<" "<< (this->psiXyl.size() - 1) ;
			if( i < (this->psiXyl.size() - 1)){
				myfile1<<" "<<std::abs((this->outputFlux[i]-outputFlux_old[i])/this->outputFlux[i]);
				myfile1 <<" "<<outputFlux[i] <<" "<<outputFlux_old[i]<<" "<<maxErr[5];
			}
			myfile1 <<std::endl;
		}
	}
	myfile1 <<std::endl<<std::endl<<std::endl;
	for (int i = 0; i < this->An.size(); i++) 
	{ //f = f < 0 ? -f : f;
		
		double tempVal = 1;
		if(this->An[i] !=0){tempVal= std::min(std::abs(this->An[i]), std::abs(An_old[i]));}else{tempVal=1.;}
		maxErr[1] =std::max(std::abs((this->An[i]-An_old[i])/tempVal),maxErr[1]);
		maxErr[8] += std::abs((this->An[i]-An_old[i])/tempVal);
		
		if(this->gco2[i] !=0){tempVal=std::min(std::abs(this->gco2[i]), std::abs(this->gco2_old[i]));}else{tempVal=1.;}
		maxErr[2] =std::max(std::abs((this->gco2[i]-gco2_old[i])/tempVal),maxErr[2]);
		
		if(this->ci[i] !=0){tempVal=std::min(std::abs(this->ci[i]), std::abs(ci_old[i]));}else{tempVal=1.;}
		maxErr[3] =std::max(std::abs((this->ci[i]-ci_old[i])/tempVal),maxErr[3]);
		
		if(this->pg[i] !=0){tempVal=std::min(std::abs(this->pg[i]), std::abs(this->pg[i]));}else{tempVal=1.;}
		maxErr[4] =std::max(std::abs((this->pg[i]-pg_old[i])/tempVal),maxErr[4]);
		
		//	if(doLog){myfile1 <<"in if(plant->organTypes[i] == 4){"<<std::endl;} }
		
		if(this->k_stomatas[i] !=0){tempVal=std::min(std::abs(this->k_stomatas[i]),std::abs(this->k_stomatas_old[i])) ;}else{tempVal=1.;}
		maxErr[6] =std::max(std::abs((this->k_stomatas[i]-k_stomatas_old[i])/tempVal),maxErr[6]);
		if(doLog){
			myfile1 <<i<<" "<<maxErr[1] <<" "<<maxErr[2]<<" "<<maxErr[3]<<" "<<maxErr[4] <<" "<<maxErr[5] <<" "<<maxErr[6]<<" an: ";
			myfile1  <<this->An[i]<<" an_old: "<<An_old[i]<<" gco2: "<<this->gco2[i]<<" gco2_old> "<<gco2_old[i];
			myfile1 <<" ci: "<<this->ci[i]<<" ci_old:"<<ci_old[i]<<" "<<pg[i]<<" "<<pg_old[i]<<" ks: ";
			myfile1 << k_stomatas[i]<<" ksold: "<< k_stomatas_old[i] <<std::endl;
		}
		assert(!std::isnan(pg[i]) && "Photosynthesis psi old guard cell is nan");
		assert(!std::isnan(pg_old[i]) && "Photosynthesis psi old guard cell is nan");
	}
	maxMaxErr = *std::max_element(maxErr.begin(), maxErr.end());
}	



	/* 
		Converts An (assimilation rate per unit of surface) [mol CO2 m-2 s-1] to Ag4Phloem (assimilation rate per segment) [mmol Suc d-1]
	*/
	
void Photosynthesis::getAg4Phloem()
{	
	this->Ag4Phloem = std::vector<double>(this->plant->nodes.size(), 0.);
	for(int idl_seg = 0; idl_seg < seg_leaves_idx.size(); idl_seg ++){
		int idl_yNode = seg_leaves_idx.at(idl_seg) + 1;
		int si = seg_leaves_idx.at(idl_seg) ;
		//from mol CO2 m-2 s-1
		//to mmol Suc cm-2 d-1
		double temp = (this->An.at(idl_seg) )/12 * 24*60*60 *1000/10000;//+ this->Rd
		double sideSurface = plant->leafBladeSurface.at(si)*2;//transpiration on both side; 2 * M_PI * a * l ;
		//to mmol Suc d-1
		this->Ag4Phloem.at(idl_yNode) =  sideSurface*temp;//to mmol Suc d-1
	}
}

	/* give all the values that only need to be computed once and not a each loop
		@param sim_time_ [day]           simulation time, needed time-dependent conductivity
	*/
void Photosynthesis::initCalcs(double sim_time_){
	
	initStruct(sim_time_);
	initVcVjRd();		
}
	
	/*  Computes variables constants at each loop linked needed for computing  the transpiration (Ev)
		@param sim_time_ [day]           simulation time, needed time-dependent conductivity
	*/
void Photosynthesis::initStruct(double sim_time_){
	lengths =  this->plant->segLength();
	assert((plant->leafBladeSurface.size() == plant->segments.size())&&"plant->leafBladeSurface.size() != plant->segments.size()");
	//sideSurface_leaf = this->plant->leafSurfacePerSeg();
	double kr, kx,  l;//a
	 int ot = Organism::ot_leaf; 
	if(fv.size() != seg_leaves_idx.size()){
		fv.resize(seg_leaves_idx.size(), 0.);
		tauv.resize(seg_leaves_idx.size(), 0.);
		dv.resize(seg_leaves_idx.size(), 0.);
		//sideSurface.resize(seg_leaves_idx.size(), 0.);
	}
	for(int li_ = 0; li_ < this->seg_leaves_idx.size(); li_++){
		int li = this->seg_leaves_idx.at(li_);
		l = lengths.at(li);		//a = plant->radii[li]; xylem data. issue: need other radius
		int st = this->plant->subTypes.at(li);
		kx = 0.;
		kr = 0.;
        try {
            kx = kx_f(li, sim_time_, st, ot);
            kr = kr_f(li, sim_time_, st, ot, li_);
        } catch(...) {
            std::cout << "\n Photosynthesis::initStruct: conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<ot<< " subtype " << st <<std::flush;
        }
		double perimeter = plant->leafBladeSurface.at(li)/l*2; 
		fv.at(li_) = -perimeter*kr;
		tauv.at(li_) = std::sqrt(perimeter*kr/kx);
		dv.at(li_) = std::exp(-tauv.at(li_)*l)-std::exp(tauv.at(li_)*l);
		if((kr<0)||(kx<=0))
		{
			std::cout<<"pt, st "<<ot<<" "<<st<<" "<<li<<" "<<li_<<std::endl;
			std::cout<<"kr etc "<<kx<<" "<<kr<<" "<<perimeter<<" "<<fv.at(li_)<<" "<<tauv.at(li_)<<" "<<dv.at(li_)<<std::endl;
			
			assert((kr>0)&&"kr<=0"); assert((kx>0)&&"kx<=0");
			assert((perimeter>0)&&"perimeter<=0");
			assert((fv[li_]>0)&&"fv[li_]<=0");
			assert((tauv[li_]>0)&&"tauv[li_]<=0");
			assert((dv[li_]>0)&&"dv[li_]<=0");
		}
	}
	
}

	/* Computes variables constants at each loop needed for computing carboxylation (Vc) and photon flux rates (J)
	*/
	
void Photosynthesis::initVcVjRd(){
	if(Vcmax.size() != seg_leaves_idx.size()){
		Vcmax.resize(seg_leaves_idx.size(), 0.);
		deltagco2.resize(seg_leaves_idx.size(), 0.);
		Vc.resize(seg_leaves_idx.size(), 0.);
		Vj.resize(seg_leaves_idx.size(), 0.);
		//An.resize(seg_leaves_idx.size(), 0.);
		fw.resize(seg_leaves_idx.size(), 0.);
		//gco2.resize(seg_leaves_idx.size(), 0.);
		Jw.resize(seg_leaves_idx.size(), 0.);
		Ev.resize(seg_leaves_idx.size(), 0.);
		//outputFluxL.resize(seg_leaves_idx.size(), 0.);
	}
	double Chl_;
	for(int li_ = 0; li_ < this->seg_leaves_idx.size(); li_++){
		//carboxylation rate
		//Vc25max
		
		if(Chl.size() != seg_leaves_idx.size()){
			Chl_ = Chl.at(0);
		}else{Chl_ = Chl.at(li_);}
		//prewprint from qian replaced with actuall article
		double Vcrefmax = (VcmaxrefChl1* Chl_ + VcmaxrefChl2)*1e-6 ;//double mol m-2 s-1
		//std::cout<<Vcrefmax <<" "<< VcmaxrefChl1<<" "<< Chl_ <<" "<< VcmaxrefChl2<<" "<<std::endl;
		//Vcmax
		double expo1 = std::exp(Eav /(R_ph*0.1*Tref)*(1. - Tref/TleafK));
		double expo2 = std::exp((S* TleafK - Edv)/(R_ph*0.1* TleafK));
		Vcmax.at(li_) = Vcrefmax * expo1 / (1. + expo2); //Eq 11
		//for Vc
		//mmol mmol-1 * exp(mJ mmol-1/(hPa cm3K−1mmol−1 *(mJ/(hPa/cm3))*K)*(-))=mmol mmol-1 * exp(-)
		Ko = Ko_ref * std::exp(Eao/(R_ph*0.1*Tref)*(1.-Tref/TleafK)); //Eq 9
		Kc = Kc_ref * std::exp(Eac/(R_ph*0.1*Tref)*(1.-Tref/TleafK)) ;//Eq 9
		delta = gamma0* (1.+ gamma1*(TleafK - Tref) + gamma2*std::pow((TleafK - Tref),2.) ) ;//Eq 10
		
		//electron transport rate
		//Jrefmax
		double Jrefmax = Vcrefmax * a3 ;//Eq 25
		//Jmax
		expo1 = std::exp(Eaj /(R_ph*0.1*Tref)*(1. - Tref/TleafK));
		expo2 = std::exp((S * TleafK - Edj)/(R_ph *0.1* TleafK));
		double Jmax = std::min(Jrefmax * expo1 / (1. + expo2), Jrefmax); //Eq 24
		//J
		double coefa = theta;
		double coefb = -(alpha * Qlight + Jmax);
		double coefc = alpha * Qlight * Jmax;
		double dis = std::pow(coefb,2.) - (4.*coefa*coefc);
		J =  ((-coefb- std::sqrt(dis))/(2.*coefa));
		
		
		Rd = Rd_ref * std::exp(Eard/(R_ph*Tref*0.1)*(1.-Tref/TleafK));
		deltagco2.at(li_) = (delta + Kc*Rd*(1. + oi/Ko)/Vcmax.at(li_))/(1-Rd/Vcmax.at(li_));
	}
}


	/*
		add gravitational wat. pot to total wat. pot. (used in phloem module)
	*/
void Photosynthesis::doAddGravity()
{
	Vector3d minB = plant->getMinBounds();
	psiXyl4Phloem.resize(psiXyl.size(), 0.);
	for(int i = 0; i<  psiXyl.size(); i++)
	{
		psiXyl4Phloem.at(i) =psiXyl.at(i) +(  plant->nodes.at(i).z - minB.z);//in cm 
	}
}


	/*
		Computes the output variables => ci, go2, An, Ev
		@param simtime		
	*/
void Photosynthesis::loopCalcs(double simTime){
	std::ofstream myfile4;
	if(doLog)
	{
		std::string name1 = "loopphoto";
		std::string name2 = ".txt";
		std::string namefile = name1 +"_"+std::to_string(simTime)+"_"+std::to_string(loop) + name2;
		myfile4.open (namefile);
	}
	for(int i = 0; i<seg_leaves_idx.size();i++)
	{
		int idl= seg_leaves_idx.at(i);
		if((verbose_photosynthesis ==2)){std::cout<<"in loopcalcs "<<i<<" "<<idl<<std::endl;}
		double l = lengths.at(idl);
		//double a = plant->radii[idl];
		double sideArea = plant->leafBladeSurface.at(idl)*2;//
		//double perimeter = sideArea/l*2;
		if((verbose_photosynthesis ==2)){std::cout<<"geom "<<l<<" "<<sideArea<<std::endl;}
		double eps = 0.;
		double rxi = this->psiXyl.at(plant->segments.at(idl).x);
		if(rxi == 0) //node just got created
		{
			int newxid = plant->segments.at(seg_leaves_idx.at(i-1)).x;
			rxi = this->psiXyl.at(newxid);
			if(doLog ){myfile4<<"new rxi "<<plant->segments.at(idl).x<<" "<< newxid <<" "<<rxi<<std::endl;}
		}
		double rxj = this->psiXyl[plant->segments.at(idl).y];
		if(rxj == 0)//node just got created
		{
			auto n1 = plant->nodes[plant->segments.at(idl).x];
			auto n2 = plant->nodes[plant->segments.at(idl).y];
			auto v = n2.minus(n1);
			double vz = v.z / l; // normed direction
			rxj = this->psiXyl.at(plant->segments.at(idl).x)-vz;
					
			if(doLog ){myfile4<<"new rxj "<<plant->segments.at(idl).y<<" "<< vz <<" "<<rxj<<std::endl;}
		}
		double p_lhPa =(rxi + rxj)*0.5*0.9806806;// cm => hPa
		//(mg mmol-1)* hPa /((hPa cm3K−1mmol−1) mg cm-3 K) =(-)
		double HRleaf = std::exp(Mh2o*this->pg.at(i)*0.9806806 /(rho_h2o*R_ph*TleafK)) ;//fractional relative humidity in the intercellular spaces
		//double ea = es - VPD;
		double ea_leaf = es * HRleaf;//hPa
		if((verbose_photosynthesis ==2)){std::cout<<"git to leaf "<<ea_leaf<<" "<<HRleaf<<std::endl;}
		if(std::abs(fv.at(i)) > 1e-16)//i.e., perimeter * kr > 1e-16 like for @see Xylem::solveLinear
		{
			if((verbose_photosynthesis ==2)){std::cout<<"in compute gco2 "<<sideArea<<" "<<(sideArea > 1e-16)<<std::endl;}
			//sideArea = 2. * M_PI * a ;//
			//carboxylation and electron transport  rate
			Vc.at(i) = std::min(std::max(Vcmax.at(i) * (ci.at(i) - delta) / (ci.at(i) + Kc*(1. + oi/Ko)),0.),Vcmax.at(i)); //Eq 8
			
			if(ci.at(i) == 2. * delta){eps = 0.001*delta ;}
			Vj.at(i) = std::max(J/4. * (ci.at(i) - delta)/ (ci.at(i) - 2. * delta+eps), 0.) ;//Eq 22
			
			//An mol m-2 s-1
			An.at(i) = std::min(Vc.at(i), Vj.at(i)) - Rd;//Eq 6
			//fw (-)
			fw.at(i) = fwr + (1.- fwr)*std::exp(-std::exp(-sh*(p_lhPa*0.0001 - p_lcrit)*10228.)) ;//Eq 5
			// mol CO2 m-2 s-1
			gco2.at(i) = g0 + fw.at(i) * a1 *( An.at(i) + Rd)/(ci.at(i) - deltagco2.at(i));//tuzet2003
			// mol H2O m-2 s-1 MPa-1
			//double k_stomate_1 = (gco2.at(i) * a2) / Patm;
			//(mol m-2 s-1)*(mmol/mol)*(hPa/hPa) * (mg mmol-1) /(mg cm-3) *(h/d)*(s/h)*(m2 m-2) =  ( cm3)/d*(cm-2)
			Jw.at(i) = (gco2.at(i) * a2) *1000* (ea_leaf - ea)/Patm * Mh2o/rho_h2o * 24.*3600*1e-4 ;//in cm3 cm-2 d-1
			Ev.at(i) = Jw.at(i)* sideArea; //in cm3 d-1
            //double f = -2*a*M_PI*kr; // flux is proportional to f // *rho*g
			
			//ci and pg
			//gruard cell wat. pot. to havee water flux from xylem to gard cell. kr = permeability of xylem membrane only.
			this->pg.at(i) = (-1/2.)*((Ev.at(i))/(-fv.at(i)*(1./(tauv.at(i)*dv.at(i)))
				*(2.-std::exp(-tauv.at(i)*l)-std::exp(tauv.at(i)*l))) - (rxi + rxj)) ;//cm
				
			k_stomatas.at(i) = Jw.at(i)/(this->pg.at(i) - psi_air);
			if((verbose_photosynthesis ==2)){
																				
				std::cout<<"sizes "<<An.size()<<" "<< gco2.size()<<" "<<ci.size()<<" "<<ci_old.size() <<std::endl;

			}
			ci.at(i) = (cs*a1*fw.at(i) +deltagco2.at(i))/(1+a1* fw.at(i)) ;//Eq 26	
			if((!std::isfinite(this->pg.at(i)))||(!std::isfinite(k_stomatas.at(i)))) {
			std::cout<<"shape leaf "<<idl<<" "<<sideArea<<" "<<ci_old.at(i)<<" "<<ci.at(i)<<std::endl;
			std::cout<<"an calc "<<An.at(i)<<" "<<Vc.at(i)<<" "<< Vj.at(i)<<" "<<J<<" "<<Vcmax.at(i)<<" "<<Kc<<" "<<Ko<<" ";
			std::cout<<" "<<delta<<" "<<oi<<" "<<eps<<std::endl;
			std::cout<<"forgco2 "<<gco2.at(i) <<" "<< g0<<" "<<  fw.at(i) <<" "<<  a1 <<" "<< An.at(i)<<" "<< Rd<<" "<< deltagco2.at(i)<<std::endl;
			std::cout<<"forJW, Jw "<<Jw.at(i)<<" drout_in "<<(this->pg.at(i) - (rxi + rxj)/2)<<" "<<ea_leaf <<" "<< ea<<" "<<Patm<<" "<<Mh2o<<" "<<rho_h2o <<std::endl;
			std::cout<<"forpg "<<sideArea <<" "<< Jw.at(i)<<" "<<fv.at(i)<<" "<<tauv.at(i)<<" "<<dv.at(i)<<" "<<l<<" "<<rxi <<" "<< rxj<<" "<<this->pg.at(i)<<" numleaf: "<<i <<std::endl;
			std::cout<<"diff Ev and lat fluw: "<<Ev.at(i)<<std::endl;//<<" "<<outputFluxL.at(i)
				throw std::runtime_error("Phtotosynthesis: nan or Inf k_stomatas.at(i) of pg.at(i)");
			}
			
		}else{ci.at(i) = 0.0;}
		if(doLog ){
			if((verbose_photosynthesis ==2)){std::cout<<"dolog loop calcs "<<std::endl;}
			myfile4<<"shape leaf "<<idl<<" "<<sideArea<<" "<<ci_old.at(i)<<" "<<ci.at(i)<<std::endl;
			myfile4<<"an calc "<<An.at(i)<<" "<<Vc.at(i)<<" "<< Vj.at(i)<<" "<<J<<" "<<Vcmax.at(i)<<" "<<Kc<<" "<<Ko<<" ";
			myfile4<<" "<<delta<<" "<<oi<<" "<<eps<<std::endl;
			myfile4<<"forgco2 "<<gco2.at(i) <<" "<< g0<<" "<<  fw.at(i) <<" "<<  a1 <<" "<< An.at(i)<<" "<< Rd<<" "<< deltagco2.at(i)<<std::endl;
			myfile4<<"forJW, Jw "<<Jw.at(i)<<" drout_in "<<(this->pg.at(i) - (rxi + rxj)/2)<<" "<<ea_leaf <<" "<< ea<<" "<<Patm<<" "<<Mh2o<<" "<<rho_h2o <<std::endl;
			myfile4<<"forpg "<<sideArea <<" "<< Jw.at(i)<<" "<<fv.at(i)<<" "<<tauv.at(i)<<" "<<dv.at(i)<<" "<<l<<" "<<rxi <<" "<< rxj<<" "<<this->pg.at(i)<<" numleaf: "<<i <<std::endl;
			myfile4<<"diff Ev and lat fluw: "<<Ev.at(i)<<std::endl;//<<" "<<outputFluxL.at(i)
		}
	}
	if(doLog){myfile4.close();}

}
			   

	/* Computes water-limited growth*/

}//end namespace
