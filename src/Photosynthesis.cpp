// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

#ifdef USE_PHOTOSYNTHESIS
#include "Photosynthesis.h"
//#include <armadillo>
//#include <algorithm>
#include <set>
#include <Eigen/Dense>
#include <Eigen/Sparse> 
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
	if((plant->kr_length > 0)&&(plant->exchangeZoneCoefs.size()!=plant->segments.size())){
		plant->calcExchangeZoneCoefs();
											 
											 
	 
											 
											
	 
			
						
	}
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
	std::cout<<"set back stomatas"<<std::endl;
	k_stomatas = k_stomatas_old ;
	std::cout<<"did back stomatas"<<std::endl;
	bool withEigen_ = true; //solve water flux with Eigen library, @see XylemFlux::linearSystem
					
 

	while(!this->stop){  
		//std::cout<<"in loop "<<loop<<std::endl;
		std::fill(maxErr.begin(), maxErr.end(), 0.);//re-initialize error vector
  
						 
	
					  
			
	
		loopCalcs(sim_time_) ;//compute photosynthesis outputs
		if((verbose_photosynthesis > 1)){std::cout<<"to linearSystem"<<std::endl;}
					  
							 
													 
																		  
   
	 
							   
							 
	 
		linearSystem(sim_time_, sxx_, cells_, soil_k_, withEigen_); //compute psiXyl
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
		int idl_yNode = seg_leaves_idx[idl_seg] + 1;
		int si = seg_leaves_idx[idl_seg] ;
		//from mol CO2 m-2 s-1
		//to mmol Suc cm-2 d-1
		double temp = (this->An[idl_seg] )/12 * 24*60*60 *1000/10000;//+ this->Rd
		double sideSurface = plant->leafBladeSurface[si]*2;//transpiration on both side; 2 * M_PI * a * l ;
		//to mmol Suc d-1
		this->Ag4Phloem[idl_yNode] =  sideSurface*temp;//to mmol Suc d-1
	}
}
	
														 
																	
																									
 
				  
				
				   
						   
		  
				 
														  
														  
																
 
 

	/* give all the values that only need to be computed once and not a each loop
		@param sim_time_ [day]           simulation time, needed time-dependent conductivity
	*/
void Photosynthesis::initCalcs(double sim_time_){
	
	std::cout<<"initStructandN "<<std::endl;
	initStruct(sim_time_);
	std::cout<<"initVcVjRd "<<std::endl;
	initVcVjRd();		
	std::cout<<"end initVcVjRd "<<std::endl;
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
		std::cout<<	"initStructandN "<<	li_<<" "<<seg_leaves_idx.size()<<" "<<lengths.size()<<" "<<ot;
		int li = this->seg_leaves_idx[li_];
		std::cout<<" "<<li<<std::endl;
		l = lengths[li];		//a = plant->radii[li]; xylem data. issue: need other radius
		int st = this->plant->subTypes[li];
		kx = 0.;
		kr = 0.;
        try {
            kx = kx_f(li, sim_time_, st, ot);
            kr = kr_f(li, sim_time_, st, ot, li_);
        } catch(...) {
            std::cout << "\n Photosynthesis::initStruct: conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<ot<< " subtype " << st <<std::flush;
        }
		std::cout<<"get leaf blade "<<plant->leafBladeSurface.size()<<" "<<li<<std::endl;
		double perimeter = plant->leafBladeSurface[li]/l*2; 
		std::cout<<"got leaf blade "<< perimeter<<std::endl;
		fv[li_] = -perimeter*kr;
		tauv[li_] = std::sqrt(perimeter*kr/kx);
		dv[li_] = std::exp(-tauv[li_]*l)-std::exp(tauv[li_]*l);	
		//sideSurface_leaf[li_] = 2.*a*M_PI*l;
		std::cout<<	"	initStructandN_end "<<	li_<<" "<<l<<" "<<kr<<" "<<kx<<" "<<st<<" "<<ot<<std::endl;
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
			Chl_ = Chl[0];
		}else{Chl_ = Chl[li_];}
		//prewprint from qian replaced with actuall article
		double Vcrefmax = (VcmaxrefChl1* Chl_ + VcmaxrefChl2)*1e-6 ;//double mol m-2 s-1
		//std::cout<<Vcrefmax <<" "<< VcmaxrefChl1<<" "<< Chl_ <<" "<< VcmaxrefChl2<<" "<<std::endl;
		//Vcmax
		double expo1 = std::exp(Eav /(R_ph*0.1*Tref)*(1. - Tref/TleafK));
		double expo2 = std::exp((S* TleafK - Edv)/(R_ph*0.1* TleafK));
		Vcmax[li_] = Vcrefmax * expo1 / (1. + expo2); //Eq 11
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
		deltagco2[li_] = (delta + Kc*Rd*(1. + oi/Ko)/Vcmax[li_])/(1-Rd/Vcmax[li_]);
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
		psiXyl4Phloem[i] =psiXyl[i] +(  plant->nodes[i].z - minB.z);//in cm 
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
		int idl= seg_leaves_idx[i];
		if((verbose_photosynthesis ==2)){std::cout<<"in loopcalcs "<<i<<" "<<idl<<std::endl;}
		double l = lengths[idl];
		//double a = plant->radii[idl];
		double sideArea = plant->leafBladeSurface[idl]*2;//
		//double perimeter = sideArea/l*2;
		if((verbose_photosynthesis ==2)){std::cout<<"geom "<<l<<" "<<sideArea<<std::endl;}
		double eps = 0.;
		double rxi = this->psiXyl[plant->segments[idl].x];
		if(rxi == 0) //node just got created
		{
			int newxid = plant->segments[seg_leaves_idx[i-1]].x;
			rxi = this->psiXyl[newxid];
			if(doLog ){myfile4<<"new rxi "<<plant->segments[idl].x<<" "<< newxid <<" "<<rxi<<std::endl;}
		}
		double rxj = this->psiXyl[plant->segments[idl].y];
		if(rxj == 0)//node just got created
		{
			auto n1 = plant->nodes[plant->segments[idl].x];
			auto n2 = plant->nodes[plant->segments[idl].y];
			auto v = n2.minus(n1);
			double vz = v.z / l; // normed direction
			rxj = this->psiXyl[plant->segments[idl].x]-vz;
					
			if(doLog ){myfile4<<"new rxj "<<plant->segments[idl].y<<" "<< vz <<" "<<rxj<<std::endl;}
		}
		double p_lhPa =(rxi + rxj)*0.5*0.9806806;// cm => hPa
		//(mg mmol-1)* hPa /((hPa cm3K−1mmol−1) mg cm-3 K) =(-)
		double HRleaf = std::exp(Mh2o*this->pg[i]*0.9806806 /(rho_h2o*R_ph*TleafK)) ;//fractional relative humidity in the intercellular spaces
		//double ea = es - VPD;
		double ea_leaf = es * HRleaf;//hPa
		if((verbose_photosynthesis ==2)){std::cout<<"git to leaf "<<ea_leaf<<" "<<HRleaf<<std::endl;}
		if(sideArea > 1e-16)
		{
			if((verbose_photosynthesis ==2)){std::cout<<"in compute gco2 "<<sideArea<<" "<<(sideArea > 1e-16)<<std::endl;}
			//sideArea = 2. * M_PI * a ;//
			//carboxylation and electron transport  rate
			Vc[i] = std::min(std::max(Vcmax[i] * (ci[i] - delta) / (ci[i] + Kc*(1. + oi/Ko)),0.),Vcmax[i]); //Eq 8
			
			if(ci[i] == 2. * delta){eps = 0.001*delta ;}
			Vj[i] = std::max(J/4. * (ci[i] - delta)/ (ci[i] - 2. * delta+eps), 0.) ;//Eq 22
			
			//An mol m-2 s-1
			An[i] = std::min(Vc[i], Vj[i]) - Rd;//Eq 6
			//fw (-)
			fw[i] = fwr + (1.- fwr)*std::exp(-std::exp(-sh*(p_lhPa*0.0001 - p_lcrit)*10228.)) ;//Eq 5
			// mol CO2 m-2 s-1
			gco2[i] = g0 + fw[i] * a1 *( An[i] + Rd)/(ci[i] - deltagco2[i]);//tuzet2003
			// mol H2O m-2 s-1 MPa-1
			//double k_stomate_1 = (gco2[i] * a2) / Patm;
			//(mol m-2 s-1)*(mmol/mol)*(hPa/hPa) * (mg mmol-1) /(mg cm-3) *(h/d)*(s/h)*(m2 m-2) =  ( cm3)/d*(cm-2)
			Jw[i] = (gco2[i] * a2) *1000* (ea_leaf - ea)/Patm * Mh2o/rho_h2o * 24.*3600*1e-4 ;//in cm3 cm-2 d-1
			Ev[i] = Jw[i]* sideArea; //in cm3 d-1
            //double f = -2*a*M_PI*kr; // flux is proportional to f // *rho*g
			
			//ci and pg
			//gruard cell wat. pot. to havee water flux from xylem to gard cell. kr = permeability of xylem membrane only.
			this->pg[i] = (-1/2.)*((Ev[i])/(-fv[i]*(1./(tauv[i]*dv[i]))
				*(2.-std::exp(-tauv[i]*l)-std::exp(tauv[i]*l))) - (rxi + rxj)) ;//cm
				
			k_stomatas[i] = Jw[i]/(this->pg[i] - psi_air);
			if((verbose_photosynthesis ==2)){
																				
				std::cout<<"sizes "<<An.size()<<" "<< gco2.size()<<" "<<ci.size()<<" "<<ci_old.size() <<std::endl;

			}
			ci[i] = (cs*a1*fw[i] +deltagco2[i])/(1+a1* fw[i]) ;//Eq 26	
			
		}else{ci[i] = 0.0;}
							
	  
								   
													  
																	  
						 
	  
		if(doLog ){
			if((verbose_photosynthesis ==2)){std::cout<<"dolog loop calcs "<<std::endl;}
			myfile4<<"shape leaf "<<idl<<" "<<sideArea<<" "<<ci_old[i]<<" "<<ci[i]<<std::endl;
			myfile4<<"an calc "<<An[i]<<" "<<Vc[i]<<" "<< Vj[i]<<" "<<J<<" "<<Vcmax[i]<<" "<<Kc<<" "<<Ko<<" ";
			myfile4<<" "<<delta<<" "<<oi<<" "<<eps<<std::endl;
			myfile4<<"forgco2 "<<gco2[i] <<" "<< g0<<" "<<  fw[i] <<" "<<  a1 <<" "<< An[i]<<" "<< Rd<<" "<< deltagco2[i]<<std::endl;
			myfile4<<"forJW, Jw "<<Jw[i]<<" drout_in "<<(this->pg[i] - (rxi + rxj)/2)<<" "<<ea_leaf <<" "<< ea<<" "<<Patm<<" "<<Mh2o<<" "<<rho_h2o <<std::endl;
			myfile4<<"forpg "<<sideArea <<" "<< Jw[i]<<" "<<fv[i]<<" "<<tauv[i]<<" "<<dv[i]<<" "<<l<<" "<<rxi <<" "<< rxj<<" "<<this->pg[i]<<" numleaf: "<<i <<std::endl;
			myfile4<<"diff Ev and lat fluw: "<<Ev[i]<<std::endl;//<<" "<<outputFluxL[i]
		}
	}
	if(doLog){myfile4.close();}

}

   
																				  
														   
  
																						  
																											
																											
														  
																						  
																														 
   
																																   
 
														  
					  
					  
					  
												   
				 
										
										  
										 
										 
			   
					   
 
									
							  
						
														 
						   
						   
											
												
					  
 
 
					  
								
							
																							 
								   
																								 
								   

	/* Computes water-limited growth
									  
																	 
					 
											  
				
															
											 
										 
																			   
																								 
				   
																					
																					
					
	 
							   
								 
											 
						
									 
				 
					
								
																			   
			 
				
													  
		 
  
   
  
																				  
																													
												 
										  
						
						 
																									
								  
								  
							  
							   

			 
													
													
					  
																								
																							
				 
   
												 
																																				
														
																		
	
  
																   
							   
							   
									   
																				   
			  
															  
	 
																																	
																						   
	
			  
														
	 
	
   
						
													 
													 
									 
																					  
											
				  
						
		 

							  
											
		 

					   
																									
															   
		 
												 

							

													 
					
																			  
							   
   
																	 
																				   
										 
																								   
													  
									   
											   
						 
						  
						  
					  
					
														   
																																	   
																
	
					   
									 
																																		 
																 
	   
			 
																																										   
	
   

								
																														 
																															   
		   

													
												   
										   
																							
									  
		 
											 
																							
									  
			   

				   
								  
																						
													 
													  
  
										   
																							
									  
			   
											 
																							
									  
			   
	 
																
																		 
												   
														
				  
															 
					  
	
	*/
												  
				 
   
											   
													
					
 
 
								 
																  
											   
				
  
 
					   
	  
				  
			 
			
																		 
				  
											   
  
		   
														   
				  
																  
						   
 




std::vector<std::map<int,double>> Photosynthesis::waterLimitedGrowth(double t)
																			   
{
																		
					   
		   
									
							 
																				   
						   
																		  
					   
																						 
															
													

								   
								   
										   
																								   
					 
				 
															
										  
   
										 
																			   
				   
																					
																				 
					
	 
	
							   
								 
											 
						
									 
				 
					
								
			 
				
													  
		 
																		 


																										  
											  
									   

					   
					   
			 
												   
												   
					  
																						
																							
		 
																			   

  
																   
							   
							   
									   
																																   
	
   
						
													 
													 
									 
															 
												 
   
																					  
										  
				 
						
		 
  
							  
										  
		 
																		
  
  
									 
						
					
																									  
																							 
									   
									   
											   
																	 
						  
			 
		  
																		 
													 
   
																			 
							   
   
																			
																							  
																	 

																	  
															   
																																		   

														  
																		   
							  
									
																																										 
				  
	
																				   
   
			
																																		  
																																				   
  
   
				 
																												   
																	 
												
	
														  
																										 
					  
				
							
						  
   
			
																																					 
  
   
  
 
																																			  
																																					
				   
	
					 
																				
														 
														 
							
																										 
										   
	
   
																																								  
																			   
								  
								  
																						  
																				   
									   
									   
																								
																					
												   
																												 
																																									
																																									
																						   
																						  
														   
	 
															  
	
																			   
	
			 
																														 
																																			  
																							  
																	 
	
																			 
			  
																						
								
								
									 
									 
												  
																							
																						
																					
																						  
																						 
														  
																  
												   
	
			 
																		  
																					  
																				  
																						
																					   
														
																
	
															  
	
   
   
   
	 
							
				  
 

									
																								
 
	std::ofstream out("calcDeltaVolOrgNode.txt");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	if(doLog){
		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
	}
	
	if((verbose_photosynthesis ==2)||doLog){
		std::cout<<"calcDeltaVolOrgNode"<<std::endl;
	}
	int Nr = plant->nodes.size();
	std::map<int,double> initMap = { { -1,0. } };//-1 : total need per node
	std::vector<std::map<int,double>> deltaVolOrgNode(Nr, initMap);//, 0.);
	bool allOrgs = true;
	auto orgs = plant->getOrgans(-1,allOrgs);//also org with length - Epsilon == 0
	//std::vector<double> Lseg = plant->segLength();//seg
	double Flen, Fpsi;
	int ot;
	for(auto org: orgs)
	{
		//organ alive and active at this time step
		if(((org->getAge()+t)>0)&&(org->isAlive())&&(org->isActive()))
		{
			
			//what happens with nodal growth? (L_realized > L_init). still need grmax > 0
			int f_gf_ind = org->getParameter("gf");
			if(f_gf_ind == 3){f_gf_ind = 1;}
			auto f_gf =  plant->createGrowthFunction(f_gf_ind);
			double Linit = org->getLength(false);//theoretical length 
			double age_ = f_gf->getAge(Linit, org->getParameter("r"), org->getParameter("k"), org->shared_from_this());//age if max growth
			
			//(length,getLeafRandomParameter()->r,param()->getK(),shared_from_this());
			double dt__= t; // time step
			double dt_; // time step
			ot = org->organType();
			
			
			double delayNG_ = 0.; //delay between organ creation and expansion
			
			double ageAtLb =f_gf->getAge(org->getParameter("lb"), org->getParameter("r"), org->getParameter("k"), org->shared_from_this());
			bool belowage = (age_+t< ageAtLb);
			//faire en sorte que la 1ere fois, on eneleve age a waitedNG
			
			bool finishedWaiting = (org->getAge() >=org->getParameter("delayNG"));
			if( (ot == 3)&&(!belowage)&&(!finishedWaiting) ){
				
				if((verbose_photosynthesis ==2)||doLog){
					std::cout<<"(ot==3)orgID_ is "<<org->getId()<<" "<<ot<<std::endl;//<<" "<<targetlength<<" "<<org->orgVolume(targetlength, false) ;
					std::cout<<"(ot==3)age_ and so on: "<< age_<<" "<< t <<" "<<ageAtLb<<" "<<org->getParameter("delayNG")<<std::endl;
					std::cout<<org->getAge()<<" "<<t<<std::endl;			
				}
				dt__=0;
				if((age_< ageAtLb)){dt__= ageAtLb - age_;}
				//dt__ = std::min(0.,org->getAge()+t-ageAtLb -org->getParameter("delayNG") );
				//if(((org->getAge()+t) >= ageAtLb)){
					//std::cout<<"(org->getAge()+t) > ageAtLb) "<<age_<<" "<<t<<" "<< ageAtLb<<std::endl;
					//std::cout<<org->getParameter("delayNG") <<" "<<org->waitedNG<<std::endl;
					double dt___ = std::max(0.,std::min(t-dt__,t-dt__- (org->getParameter("delayNG") -org->getAge())));
					//org->waitedNG += t - dt___;
					dt__ += dt___;
					//std::cout<<"changed dt__: "<<dt__<<" "<<org->getAge() <<std::endl;
				//}//
			}
			if ((org->getAge()+ dt__)<dt__) { // the root emerged in this time step, adjust time step
				//std::cout<<"((org->getAge()+t-delayNG_)<t) "<<org->getAge()<<" "<<dt__<<std::endl;
				
				if ((org->getAge()+dt__)>org->getParameter("rlt")) { // root life time
					//std::cout<<"((org->getAge()+t-delayNG_)>org->getParameter(rlt)) "<<std::endl;
					dt_=std::max(0.,org->getParameter("rlt")-delayNG_);//dt = life span
				}else{dt_= std::max(0.,org->getAge()+dt__);}
			} else {
				
				if ((org->getAge()+dt__)>org->getParameter("rlt")) { // root life time
				dt_= org->getParameter("rlt") -org->getAge(); // remaining life span
				}else{dt_= std::max(0.,dt__);}
			}
			if((verbose_photosynthesis ==2)||doLog){
				std::cout<<"(ot==?)orgID_ is "<<org->getId()<<" "<<ot<<" "<< (ot == 3)<<std::endl;			
			}
			
			if((not((org->getOrganRandomParameter()->f_gf->CW_Gr.empty()) || 
					(org->getOrganRandomParameter()->f_gf->CW_Gr.count(org->getId()) ==0) ||
					(org->getOrganRandomParameter()->f_gf->CW_Gr.find(org->getId())->second<0.)))&&
					org->isActive())
			{
				std::cout<<"error cw_gr2"<<std::endl;
				std::cout<<org->getId()<<" "<<org->getOrganRandomParameter()->f_gf->CW_Gr.find(org->getId())->second<<std::endl;
				std::cout<<org->calcLength(1)<<" "<< ot <<" "<<org->getAge()<<std::endl;
				std::cout<<"error cw_gr2"<<std::endl;
				assert(false);
			}
			if((verbose_photosynthesis ==2)||doLog){
				std::cout<<"orgID_ is "<<org->getId()<<" "<<ot<<" "<<f_gf_ind<<" "<<org->getLength(true)<<" "<<Linit<<std::endl;//<<" "<<targetlength<<" "<<org->orgVolume(targetlength, false) ;
				std::cout<<"age_ and so on: "<< age_<<" "<< dt_ <<" "<<delayNG_<<" "<<(age_+dt_ -delayNG_)<<" "<<org->getParameter("r")<<" "<<org->getParameter("lb")<<std::endl;
				std::cout<<org->getAge()<<" "<<t<<std::endl;			
			}
			//double targetlength = org->calcLength(age_+dt_ *increase);
			double targetlength = f_gf->getLength(age_+dt_ , org->getParameter("r"), org->getParameter("k"), org->shared_from_this());
			//double targetlength = targetlength_;//+ org->getEpsilon();

			double e = targetlength-Linit; // unimpeded elongation in time step dt
			//double Rmax_org = e/t;
			std::vector<int> nodeIds_;// = org->getNodeIds();
			if((e + Linit)> org->getParameter("k")){
				std::cout<<"Photosynthesis::rmaxSeg: target length too high "<<e<<" "<<t<<" "<<Linit;
				std::cout<<" "<<org->getParameter("k")<<" "<<org->getId()<<std::endl;
				assert(false);
			}
			//delta_length to delta_vol
			double deltavol = org->orgVolume(targetlength, false) - org->orgVolume(Linit, false);//volume from theoretical length
			
			int nNodes = org->getNumberOfNodes();
			if ((nNodes==1)||(ot == 2)) {//organ not represented because below dx limit or is root
				//rather do age >0 and getlength <0 no?
				nodeIds_.push_back(-1); //because count start at 1
				nodeIds_.push_back(org->getNodeId(nNodes-1));	//globalID of parent node
				//std::cout<<"	data small org "<<o->getParameter("radius")<<" "<<o->getParameter("r")<<" "<<o->parentNI<<std::endl;
				//Linit = 1; //to avoid division by zero when computing => not needed, as Flen set to 1 
								
			}else{nodeIds_ = org->getNodeIds();}
			
			double Linit_realized = org->getLength(true);
			if((verbose_photosynthesis ==2)||doLog){
				std::cout<<"orgID is "<<org->getId()<<" "<<ot<<" "<<f_gf_ind<<" "<<Linit_realized<<" "<<Linit<<" "<<targetlength<<" "<<org->orgVolume(targetlength, false) ;
				std::cout<<" "<< org->orgVolume(Linit, false)<<" "<<deltavol<<" "<<nNodes<<" "<<org->getEpsilon()<<std::endl;
				std::cout<<"nodeIds_.size() "<<nodeIds_.size()<<std::endl;
			
				for(int k=0;k< nodeIds_.size();k++){std::cout<<nodeIds_[k]<<" ";}
				std::cout<<std::endl;
			}
			double Flen_tot = 0.;
			double deltaVol_tot = 0.;
			//auto orgNodes = org->getNodes();
			for(int k=1;k< nodeIds_.size();k++)//don't take  first node
			{
				//int segId = nodeIds_[k] - 1;//segID = seg.y -1
				int nodeId = nodeIds_[k];
				int nodeId_h = -1;
				double Lseg = 0.;
				
				bool isRootTip = ((ot==2)&&(k==(nodeIds_.size()-1)));
				if((nNodes==1)||isRootTip){Flen = 1.;
				}else{
					
					nodeId_h = nodeIds_[k-1];
					Lseg = org->getLength(k) - org->getLength(k-1);//getLength uses local node ID
					Flen = (Lseg/Linit_realized) * double(ot != 2) ;
					if((verbose_photosynthesis ==2)||doLog){
						std::cout<<"		long stem or leaf "<<nodeId<<" "<<nodeId_h<<" "<<org->getLength(k) <<" "<< org->getLength(k-1)<<" "<<Lseg<<std::endl;
					}
					auto nodei = org->getNode(k-1);
					auto nodej = org->getNode(k);
					double length2 = nodej.minus(nodei).length();
					if((verbose_photosynthesis ==2)||doLog){
						std::cout<<"		"<<nodei.toString()<<" "<<nodej.toString()<<" "<<nodej.minus(nodei).toString() <<" "<<length2<<std::endl;
					}
				}//att! for this division we need the realized length, not the theoretical one
				
				if(psiXyl.size()>0){
					Fpsi = (std::min(psiXyl[nodeId], psiMax) - psiMin)/(psiMax - psiMin);
				}else{Fpsi = 1.;}
				double deltavolSeg = deltavol * Flen * Fpsi;
				if((deltavolSeg<0.)||(deltavolSeg != deltavolSeg)){
					//could be error of pressision (if l = Lmax)
					// or that, because of nodal growth and dxMin, org->getEpsilon() <0
					if(((deltavolSeg>-1e-5)||(targetlength>0.))
							&&(not(deltavolSeg != deltavolSeg))){
						deltavolSeg=0.;	//within margin of error
					}else{
						std::cout<<org->getId()<<" t:"<<t<<" ot:"<<ot<<" Li:"<<Linit<<" Le:"<<targetlength<<std::endl;
						std::cout<<"		k "<<k<<" "<<" id:"<<nodeId<<" Flen:"<<Flen <<" Fpsi:"<< Fpsi;
						std::cout<<" Rtip:"<<isRootTip<<" Lseg:"<<Lseg<<" rorg"<<e<<" "<<deltavolSeg;
						std::cout<<" "<<e<<" "<<(targetlength-Linit)<<std::endl;
						std::cout<<"coucou"<<std::endl;
						throw std::runtime_error("(deltavolSeg<0.)||(deltavolSeg != deltavolSeg)");
						//assert(false);
					}
					
					
				}
				deltaVolOrgNode[nodeId][org->getId()] = deltavolSeg;
				deltaVolOrgNode[nodeId][-1] += deltavolSeg;
				Flen_tot += Flen;
				deltaVol_tot += deltavolSeg;
				if((verbose_photosynthesis ==2)||doLog){
					std::cout<<"		k "<<k<<" "<<" id:"<<nodeId<<" idh:"<<nodeId_h<<" Flen:"<<Flen <<" Fpsi:"<< Fpsi;
					std::cout<<" Rtip:"<<isRootTip<<" Lseg:"<<Lseg<<" "<<deltavolSeg<<std::endl;
					
						std::cout<<"		"<<org->getId()<<" t:"<<t<<" ot:"<<ot<<" Li:"<<Linit<<" Le:"<<targetlength;
						std::cout<<" rorg "<<e<<" "<<deltavol<<" "<<(targetlength-Linit)<<std::endl;
						std::cout<<"		"<<org->getLength(true)<<" "<<org->getEpsilon()<<std::endl;
				}
			}
			if((verbose_photosynthesis ==2)||doLog){
				std::cout<<"Flen_tot "<<Flen_tot<<" "<<(Flen_tot == 1.)<<std::endl;
				std::cout<<"Flen_tot "<<( 1. - Flen_tot )<<std::endl;
				std::cout<<"Flen_tot "<<(Flen_tot < 1.)<<" "<<(Flen_tot > 1.)<<std::endl;
				std::cout<<"Flen_tot "<<(Flen_tot <= 1.)<<" "<<(Flen_tot >= 1.)<<std::endl;
			}
			assert((std::abs(Flen_tot - 1.)<1e-10)&&"wrong tot Flen");
			if((verbose_photosynthesis ==2)||doLog){
				std::cout<<"vol_tot "<<deltaVol_tot<<" "<<deltavol<<" "<<(deltaVol_tot -deltavol)<<std::endl;
			}
			assert(((deltaVol_tot -deltavol)<1e-10)&&"deltavol_tot too high");//deltaVol_tot <=deltavol
		}
	
	}
	
	if(doLog){std::cout<<"\n\n\n"<<std::endl;std::cout.rdbuf(coutbuf);} //reset to standard output again
		
	return deltaVolOrgNode;
}


   
											
  
																							  
																										  
													
															
																	  
  
																			
   
																															  
												  
 
																		 
																			
   
																															   
																				  
 

// void Photosynthesis::r_forPhloem(double lightTimeRatio, int ot)
// {
	
	// bool allOrgs = true;
	// auto orgs = plant->getOrgans(ot,allOrgs);//also org with length - Epsilon == 0
	// for(auto org: orgs)
	// {
		// double oldR = org->getParameter("r");
		// org->getParam()->setParam("r",lightTimeRatio);
		// double newR = org->getParameter("r");
		// if(std::abs(newR/oldR - lightTimeRatio) > 1e-14){
			// std::cout<<"lightratio "<<oldR<<" "<<newR<<" "<<lightTimeRatio<<" "<<std::abs(newR/oldR - lightTimeRatio)<<std::endl;
			// assert(false);
		// }
		// std::cout<<"lightratio "<<org->organType()<<" "<<oldR<<" "<<newR<<" "<<lightTimeRatio<<" "<<std::abs(newR/oldR - lightTimeRatio)<<std::endl;
	// }
									  
																											 
								   
												  
							
																							   
		 
													  
																		
		 
	 
				 
 


	  
												
									   
	  
					  
																			
					
																											   
	

					  
																			
					   
																											   
// }

															   
 
 
					 
																			   
					
  
									   
												
									   
												   
																														  
				 
   
																																				
  
 
}//end namespace
#endif