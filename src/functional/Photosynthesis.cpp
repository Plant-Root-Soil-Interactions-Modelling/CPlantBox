// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

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
Photosynthesis::Photosynthesis(std::shared_ptr<CPlantBox::MappedPlant> plant_, std::shared_ptr<CPlantBox::PlantHydraulicParameters> params, double psiXylInit, double ciInit):
	PlantHydraulicModel(std::shared_ptr<CPlantBox::MappedSegments>(plant_), params), plant(plant_), psiXylInit(psiXylInit), ciInit(ciInit)
{//check when plant and planphotosyn are diff
	//std::cout<<"alive creation "<<std::endl;
	//this->seg_leaves_idx = plant->getNodeIds(4);//ids of leaf segments
	this->seg_leaves_idx = plant->getSegmentIds(4);//ids of leaf segments
	psiXyl.resize(plant->nodes.size(), psiXylInit);//-500
	An.resize(seg_leaves_idx.size(), 0.);
	gco2.resize(seg_leaves_idx.size(), 0.);
	ci.resize(seg_leaves_idx.size(), ciInit);//350e-6
	pg.resize(seg_leaves_idx.size(), psiXylInit);
    gtotOx.resize(seg_leaves_idx.size(), 0.);
	PVD.resize(seg_leaves_idx.size(), 0.);
	hrelL.resize(seg_leaves_idx.size(), 0.);
	EAL.resize(seg_leaves_idx.size(), 0.);
}


	/* solves the coupled water flux + carbon assimilation and stomatal opening,
		@param sim_time_ [day]           simulation time
		@param sxx_ [cm]                 soil matric potentials given per segment or per soil cell
		@param es [hPa]                 atmospheric humidity at saturation
		@param ea [hPa]                 atmospheric humidity
		@param cells_                    indicates if the matric potentials are given per cell (True) or by segments (False)
		@param soil_k [day-1]           optionally, soil conductivities can be prescribed per segment,
										conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)
		@param doLog_                    indicates if computed values should be printed in a text file (True) or not (False)
		@param verbose_                  print at runtime nothing (0), sparsly (1), many outputs (2)
		@param TleafK_ [°K]               leaf temperature (mean)
	*/

void Photosynthesis::solve_photosynthesis( double sim_time_ , std::vector<double> sxx_, 
				double ea_,double es_, std::vector<double> TleafK_, 
				bool cells_ ,std::vector<double> soil_k_, bool doLog_ , int verbose_ , 
				std::string outputDir_)
{
	//		save Environmental and other input variables
	doLog = doLog_; verbose_photosynthesis = verbose_;

	loop = 0;
	this->stop = false;
	this->seg_leaves_idx = plant->getSegmentIds(4);//ids of leaf segments
	this->TleafK = TleafK_;
	this->es = es_;
	this->ea = ea_;
	//cm = log(-) * (mg cm-3) * (hPa cm3 K−1 mmol−1) * K * (1/[mg mmol-1]) * (cm/hPa)
    double RH_ = this->ea /this->es;
    if((RH_>=1)||(RH_<=0))
    {
        throw std::runtime_error("Photosynthesis::solve_photosynthesis : (RH_>=1)||(RH_<=0)");
    }
	double TleafK_mean = std::accumulate(this->TleafK.begin(), this->TleafK.end(), 0.0) / this->TleafK.size();
	this->params->psi_air = std::log(RH_) * rho_h2o * R_ph * TleafK_mean/Mh2o * (1/0.9806806)  ; //in cm
	assert(((plant->kr_length < 0)||(plant->exchangeZoneCoefs.size()==plant->segments.size()))&&"(plant->exchangeZoneCoefs.size()==plant->segments.size()) while kr_length >0");
	//		creat first guesses arrays + "old" values
	psiXyl= std::vector<double>( plant->nodes.size(), psiXylInit);//-500
	An= std::vector<double>( seg_leaves_idx.size(), 0.);
	gco2= std::vector<double>( seg_leaves_idx.size(), 0.);
	ci= std::vector<double>( seg_leaves_idx.size(), ciInit);
	pg= std::vector<double>( seg_leaves_idx.size(), 0.);
	gtotOx= std::vector<double>( seg_leaves_idx.size(), 0.);
	PVD= std::vector<double>( seg_leaves_idx.size(), 0.);
	hrelL= std::vector<double>( seg_leaves_idx.size(), 0.);
	EAL= std::vector<double>( seg_leaves_idx.size(), 0.);


	orgsVec = plant->getOrgans(-1);

	An_old = An; gco2_old = gco2; ci_old = ci;
	outputFlux_old.resize(psiXyl.size(), 0.);
	outputFlux.resize(plant->segments.size(), 0.);
	psiXyl_old = psiXyl; pg_old = pg; //k_stomatas_old = k_stomatas;

	//		compute parameters which do not change between the loops
	initCalcs(sim_time_);
	//k_stomatas = k_stomatas_old ;
	while(!this->stop){
		std::fill(maxErrAbs.begin(), maxErrAbs.end(), 0.);//re-initialize absolute error vector
		std::fill(maxErr.begin(), maxErr.end(), 0.);//re-initialize relative error vector
		loopCalcs(sim_time_, sxx_, cells_) ;//compute photosynthesis outputs
		if((verbose_photosynthesis > 1)){std::cout<<"to linearSystem"<<std::endl;}
		linearSystemSolve(sim_time_, sxx_, cells_, soil_k_); //compute psiXyl
		if((verbose_photosynthesis > 1)){std::cout<<"to outputFlux"<<std::endl;}
		//usefull only to know whether we reached convergence
		outputFlux = getRadialFluxes(sim_time_, this->psiXyl, sxx_, false, cells_, soil_k_);//approx = false
		if((verbose_photosynthesis > 1)){std::cout<<"to getError"<<std::endl;}
		getError(sim_time_);
		this->stop = canStop();

		loop++ ;
		
		if(!this->stop){
			outputFlux_old = outputFlux;
			ci_old = this->ci; pg_old = this->pg;
			psiXyl_old = this->psiXyl; An_old = this->An; gco2_old = this->gco2;
		}

		if((verbose_photosynthesis > 0)){
			std::cout<<"leuning computation module at "<<(loop-1)<<" trials. "
			"Sum of max relative error calculated at the "
			"last two trials: stop? "<<(this->stop)<<std::endl;
			std::cout<<"each val "<<maxErr[0]<<" "<<maxErr[1]<<" "<<maxErr[2]<<" "<<maxErr[3]<<" "<<maxErr[4];
			std::cout<<" "<<maxErr[5]<<" "<<maxErr[6]<<" "<<maxErr[7]<<" "<<maxErr[8]<<std::endl;
			std::cout<<"each val abs "<<maxErrAbs[0]<<" "<<maxErrAbs[1]<<" "<<maxErrAbs[2]<<" "<<maxErrAbs[3]<<" "<<maxErrAbs[4];
			std::cout<<" "<<maxErrAbs[5]<<" "<<maxErrAbs[6]<<" "<<maxErrAbs[7]<<" "<<maxErrAbs[8]<<std::endl;
		}

	}


	 
	//loopCalcs(sim_time_, sxx_, cells_) ;//compute photosynthesis outputs. when fw ~ 0, need to do it one last time to be sure that seg_flux_leaf = Ev # TODO: add as error matrix because then the sum of fluxes becomes wrong
	// outputFlux = getRadialFluxes(sim_time_, this->psiXyl, sxx_, false, cells_, soil_k_);//approx = false
	loop++ ;

	// for phloem flow
	getAg4Phloem();
	doAddGravity();

	if(loop>maxLoop)
	{
		throw std::runtime_error("photosynthesis::solve: did not reach convergence");
	}

	if((verbose_photosynthesis  > 0))
	{
		std::cout<<"leuning computation module stopped after "<<(loop-1)<<" trials. "
		"Sum of max relative error calculated at the "
		"last two trials: "<<(this->stop)<<std::endl;
		std::cout<<"each val "<<maxErr[0]<<" "<<maxErr[1]<<" "<<maxErr[2]<<" "<<maxErr[3]<<" "<<maxErr[4];
		std::cout<<" "<<maxErr[5]<<" "<<maxErr[6]<<" "<<maxErr[7]<<" "<<maxErr[8]<<std::endl;
			std::cout<<"each val abs "<<maxErrAbs[0]<<" "<<maxErrAbs[1]<<" "<<maxErrAbs[2]<<" "<<maxErrAbs[3]<<" "<<maxErrAbs[4];
			std::cout<<" "<<maxErrAbs[5]<<" "<<maxErrAbs[6]<<" "<<maxErrAbs[7]<<" "<<maxErrAbs[8]<<std::endl;
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
	
	int Ns = ms->segments.size(); // number of segments
    int N = ms->nodes.size(); // number of nodes
	tripletList.clear();
	tripletList.reserve(Ns*4);
	b = Eigen::VectorXd(N);
	//get "tripletList" and "b"
	linearSystemMeunier(simTime_, sxx_, cells_, soil_k_); //see XylemFlux::linearSystem
	Eigen::SparseMatrix<double> mat(N,N);
	mat.reserve(Eigen::VectorXi::Constant(N,2));
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	mat.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
	lu.compute(mat);

	if(lu.info() != Eigen::Success){
		std::cout << "XylemFlux::linearSystem  matrix Compute with Eigen failed: " << lu.info() << std::endl;
		throw std::runtime_error("XylemFlux::linearSystem  matrix Compute with Eigen failed" );
	}

	Eigen::VectorXd v2;
	try{
		v2= lu.solve(b);
	}catch(...){
		 throw std::runtime_error("XylemFlux::linearSystem error when solving wat. pot. xylem with Eigen ");
	}
	std::vector<double> v3(&v2[0], v2.data()+v2.cols()*v2.rows());
	psiXyl = v3;
}

/**
 * fill the matrices to be solved. overloads @see XylemFlux::fillVectors
 * @param k				index for the row- and column-index vectors
 * @param i, j			indexes of the non-zero elements of the sparse matrix
 * @param psi_s 		outer water potential [cm]
 * @param bi			value of variable b at row i [cm3/d]
 * @param cii			value of variable c at row i col i [cm2/d]
 * @param cij			value of variable c at row i col j [cm2/d]
 * @return k			next index for the row- and column-index vectors
 */
size_t Photosynthesis::fillVectors(size_t k, int i, int j, double bi, double cii, double cij, double psi_s) 
{
	typedef Eigen::Triplet<double> Tri;
	aB[i] += ( bi + cii * psi_s +cij * psi_s) ;

	b(i) = aB[i];
	tripletList.push_back(Tri(i,i,cii));
	k += 1;
	tripletList.push_back(Tri(i,j,cij));
	k += 1;

	int ii = i;
	i = j;  j = ii; // edge ji
	aB[i] += ( -bi + cii * psi_s +cij * psi_s) ; // (-bi) Eqn (14) with changed sign

	b(i) = aB[i];
	tripletList.push_back(Tri(i,i,cii));
	k += 1;

	tripletList.push_back(Tri(i,j,cij));
	k += 1;
	return k;
}


/**
 *  give outer water potential [cm] overloads by @see Xylem::getPsiOut
 * @param cells 		sx per cell (true), or segments (false)
 * @param si 			segment index
 * @param sx        [cm] soil matric potential for each cell
 */
 
double Photosynthesis::getPsiOut(bool cells, int si, const std::vector<double>& sx_) const
{
	int organType = plant->organTypes.at(si);
    double psi_s;
	if (cells) { // soil matric potential given per cell
		int cellIndex = plant->seg2cell[si];
		if (cellIndex>=0) {
			if(organType ==Organism::ot_leaf){ 
				if(verbose_photosynthesis>0)
				{
					std::cout<<"Photosynthesis::linearSystem: Leaf segment n#"<<si<<" below ground. OrganType: ";
					std::cout<<organType<<" cell Index: "<<cellIndex<<std::endl;
				}
				psi_s = 0;
				//throw std::runtime_error("Photosynthesis::linearSystem: Leaf segment is belowground.");
			}
			if(sx_.size()>1) {
				psi_s = sx_.at(cellIndex);
			} else {
				psi_s = sx_.at(0);
			}
		} else 
		{
			switch(organType) {
				case Organism::ot_root: 
						if(verbose_photosynthesis>0)
					{
						std::cout<<"Photosynthesis::linearSystem: Root segment n#"<<si<<" aboveground. OrganType: ";
						std::cout<<organType<<" cell Index: "<<cellIndex<<std::endl;
					}
					psi_s = 0;
						//throw std::runtime_error("Photosynthesis::linearSystem: root segment is aboveground.");
					break;
				case  Organism::ot_stem: 
					psi_s = params->psi_air;
					break;
				case Organism::ot_leaf: 
					psi_s = pg.at(plant->getSegment2leafId(si));
					break;
				default:
					throw std::runtime_error("Photosynthesis::getPsiOut: organType not recognized.");
			}
		}
	}else{ //TODO: use better the psi_s as given by the air_ModelsPlant == boundary layer
		switch(organType) {
			case Organism::ot_root: 
				psi_s = sx_.at(si); // j-1 = segIdx = s.y-1
				//throw std::runtime_error("Photosynthesis::linearSystem: root segment is aboveground.");
				break;
			case  Organism::ot_stem: 
				psi_s = sx_.at(si); // j-1 = segIdx = s.y-1. should be == psi_air
				break;
			case Organism::ot_leaf: 
				psi_s = pg.at(plant->getSegment2leafId(si));//sx_.at(si) is used in loop function
				break;
			default:
				throw std::runtime_error("Photosynthesis::getPsiOut: organType not recognized.");
		}
	}		
	return psi_s;
}

/*
		Check if convergence reached
*/

bool Photosynthesis::canStop()
{
    bool canStop = false;
    if(loop > maxLoop)
    {
        canStop = true;
    }else{
        if(loop > minLoop)
        {
            canStop = true;
            for(int tt = 0; tt<maxErr.size();tt++)
            {
                if((maxErr.at(tt)>maxErrLim.at(tt))&&(maxErrAbs.at(tt)>maxErrAbsLim.at(tt)))
                {
                    canStop = false;
                }
            }
        }
    }
    return canStop;
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
		myfile1.open (outputDir+namefile);
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
			std::cout<< "Photosynthesis xyl >0, "<<i<<" "<<this->psiXyl[i]<<std::flush;
			throw std::runtime_error("Photosynthesis::getError : xyl >0");
		}

		double tempVal = 1;

		if(this->psiXyl[i] !=0){tempVal= std::min(std::abs(this->psiXyl[i]), std::abs(this->psiXyl_old[i]));}else{tempVal=1.;}
		maxErrAbs[0] =std::max(std::abs((this->psiXyl[i]-psiXyl_old[i])),maxErrAbs[0]);
		maxErr[0] =std::max(std::abs((this->psiXyl[i]-psiXyl_old[i])/tempVal),maxErr[0]);
		if( i < (this->psiXyl.size() - 1)){
			if(this->outputFlux[i] !=0){tempVal=  std::min(std::abs(this->outputFlux[i]), std::abs(this->outputFlux_old[i]));
				}else{tempVal=1.;}
			maxErrAbs[5] =std::max(std::abs((this->outputFlux[i]-outputFlux_old[i])),maxErrAbs[5]);
			maxErrAbs[7] +=(this->outputFlux[i]);//becasue sum of fluxes important
			maxErr[5] =std::max(std::abs((this->outputFlux[i]-outputFlux_old[i])/tempVal),maxErr[5]);
			maxErr[7] += this->outputFlux[i]/tempVal;//becasue sum of fluxes important
		}

		if(doLog){
			myfile1 <<i<<" "<<", maxErr[5] "<<plant->organTypes[std::max(i-1,0)]<<" "<<maxErr[0];
			myfile1<<" "<<psiXyl[i]<<" "<<psiXyl_old[i]<<" "<<i<<" "<< (this->psiXyl.size() - 1) ;
			if( i < (this->psiXyl.size() - 1)){
				myfile1<<" "<<std::abs((this->outputFlux[i]-outputFlux_old[i])/this->outputFlux[i]);
				myfile1 <<" "<<outputFlux[i] <<" "<<outputFlux_old[i]<<" "<<maxErr[5]<<" "<<maxErr[7];
			}
			myfile1 <<std::endl;
		}
	}
	maxErrAbs[7] = std::abs(maxErrAbs[7]);
	maxErr[7] = std::abs(maxErr[7]);
	myfile1 <<std::endl<<std::endl<<std::endl;
	for (int i = 0; i < this->An.size(); i++)
	{ //f = f < 0 ? -f : f;

		double tempVal = 1;
		if(this->An[i] !=0){tempVal= std::min(std::abs(this->An[i]), std::abs(An_old[i]));}else{tempVal=1.;}
		maxErrAbs[1] =std::max(std::abs((this->An[i]-An_old[i])),maxErrAbs[1]);
		maxErrAbs[8] += std::abs((this->An[i]-An_old[i]));

		maxErr[1] =std::max(std::abs((this->An[i]-An_old[i])/tempVal),maxErr[1]);
		maxErr[8] += std::abs((this->An[i]-An_old[i])/tempVal);

		if(this->gco2[i] !=0){tempVal=std::min(std::abs(this->gco2[i]), std::abs(this->gco2_old[i]));}else{tempVal=1.;}
		maxErrAbs[2] =std::max(std::abs((this->gco2[i]-gco2_old[i])),maxErrAbs[2]);
		maxErr[2] =std::max(std::abs((this->gco2[i]-gco2_old[i])/tempVal),maxErr[2]);

		if(this->ci[i] !=0){tempVal=std::min(std::abs(this->ci[i]), std::abs(ci_old[i]));}else{tempVal=1.;}
		maxErrAbs[3] =std::max(std::abs((this->ci[i]-ci_old[i])),maxErrAbs[3]);
		maxErr[3] =std::max(std::abs((this->ci[i]-ci_old[i])/tempVal),maxErr[3]);

		if(this->pg[i] !=0){tempVal=std::min(std::abs(this->pg[i]), std::abs(this->pg[i]));}else{tempVal=1.;}
		maxErrAbs[4] =std::max(std::abs((this->pg[i]-pg_old[i])),maxErrAbs[4]);
		maxErr[4] =std::max(std::abs((this->pg[i]-pg_old[i])/tempVal),maxErr[4]);

		int idl= seg_leaves_idx.at(i); // global segment indx
		maxErrAbs[6] = std::abs((this->outputFlux[idl]-this->Ev.at(i)));
		maxErr[6] = std::abs((this->outputFlux[idl]-this->Ev.at(i)))/std::min(std::abs(this->outputFlux[i]),std::abs(this->Ev.at(i)));
		
		if(doLog){
			myfile1 <<i<<" abs: "<<maxErrAbs[1] <<" "<<maxErrAbs[1] <<" "<<maxErrAbs[2]<<" "<<maxErrAbs[3]<<" "<<maxErrAbs[4] <<" "<<maxErrAbs[5] <<" "<<maxErrAbs[6]<<" "<<maxErrAbs[7]<<" "<<maxErrAbs[8]<<" an: ";
			myfile1  <<this->An[i]<<" an_old: "<<An_old[i]<<" gco2: "<<this->gco2[i]<<" gco2_old> "<<gco2_old[i];
			myfile1 <<" ci: "<<this->ci[i]<<" ci_old:"<<ci_old[i]<<" "<<pg[i]<<" "<<pg_old[i]<<std::endl;
            
			myfile1 <<i<<" "<<maxErr[1] <<" "<<maxErr[2]<<" "<<maxErr[3]<<" "<<maxErr[4] <<" "<<maxErr[5] <<" "<<maxErr[6]<<" an: ";
			myfile1  <<this->An[i]<<" an_old: "<<An_old[i]<<" gco2: "<<this->gco2[i]<<" gco2_old> "<<gco2_old[i];
			myfile1 <<" ci: "<<this->ci[i]<<" ci_old:"<<ci_old[i]<<" "<<pg[i]<<" "<<pg_old[i]<<std::endl;
		}
		assert(!std::isnan(pg[i]) && "Photosynthesis psi old guard cell is nan");
		assert(!std::isnan(pg_old[i]) && "Photosynthesis psi old guard cell is nan");
	}
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
	for(int i = 0; i < this->seg_leaves_idx.size(); i++){
		int li = this->seg_leaves_idx.at(i);
		l = lengths.at(li);		//a = plant->radii[li]; xylem data. issue: need other radius
		int st = this->plant->subTypes.at(li);
		kx = 0.;
		kr = 0.;
        try {
            kx = this->params->kx_f(li, sim_time_, st, ot);
            kr = kr_f(li, sim_time_, st, ot);
        } catch(...) {
            std::cout << "\n Photosynthesis::initStruct: conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<ot<< " subtype " << st <<std::flush;
        }
		double perimeter = plant->leafBladeSurface.at(li)/l*2;
		fv.at(i) = -perimeter*kr;
		tauv.at(i) = std::sqrt(perimeter*kr/kx);
		dv.at(i) = std::exp(-tauv.at(i)*l)-std::exp(tauv.at(i)*l);
		if((kr<0)||(kx<=0))
		{
			std::cout<<"pt, st "<<ot<<" "<<st<<" "<<li<<" "<<i<<std::endl;
			std::cout<<"kr etc "<<kx<<" "<<kr<<" "<<perimeter<<" "<<fv.at(i)<<" "<<tauv.at(i)<<" "<<dv.at(i)<<std::endl;

			assert((kr>0)&&"kr<=0"); assert((kx>0)&&"kx<=0");
			assert((perimeter>0)&&"perimeter<=0");
			assert((fv[i]>0)&&"fv[i]<=0");
			assert((tauv[i]>0)&&"tauv[i]<=0");
			assert((dv[i]>0)&&"dv[i]<=0");
		}
	}

}

	/* Computes variables constants at each loop needed for computing carboxylation (Vc) and photon flux rates (J)
	*/

void Photosynthesis::initVcVjRd(){
	if(Vcmax.size() != seg_leaves_idx.size()){
		Vcmax.resize(seg_leaves_idx.size(), 0.);
		deltagco2.resize(seg_leaves_idx.size(), 0.);
		delta.resize(seg_leaves_idx.size(), 0.);
		Vc.resize(seg_leaves_idx.size(), 0.);
		Vj.resize(seg_leaves_idx.size(), 0.);
		fw.resize(seg_leaves_idx.size(), 0.);
		Jw.resize(seg_leaves_idx.size(), 0.);
		Ev.resize(seg_leaves_idx.size(), 0.);
		Rd_ref.resize(seg_leaves_idx.size(), 0.);
		Rd.resize(seg_leaves_idx.size(), 0.);
        Vcrefmax = std::vector<double>(seg_leaves_idx.size(), 0.);
        Jrefmax = std::vector<double>(seg_leaves_idx.size(), 0.); 
		Ko.resize(seg_leaves_idx.size(), 0.);
		Kc.resize(seg_leaves_idx.size(), 0.);
		
		
		//C3
		Jmax.resize(seg_leaves_idx.size(), 0.);
		J.resize(seg_leaves_idx.size(), 0.);
		//for C4		
		Vp.resize(seg_leaves_idx.size(), 0.);
		kp.resize(seg_leaves_idx.size(), 0.);
		kp25.resize(seg_leaves_idx.size(), 0.);
	}
	for(int i = 0; i < this->seg_leaves_idx.size(); i++){
		//carboxylation rate
		//Vc25max
		double Chl_ = getMeanOrSegData(Chl, i);		
			//prewprint from qian replaced with actuall article
		Vcrefmax.at(i) = (VcmaxrefChl1* Chl_ + VcmaxrefChl2)*1e-6 ;//double mol m-2 s-1
		
		//mmol mmol-1 * exp(mJ mmol-1/(hPa cm3K−1mmol−1 *(mJ/(hPa/cm3))*K)*(-))=mmol mmol-1 * exp(-)
		Ko.at(i) = Ko_ref * Arrhenius(i, Eao); //Eq 9
		Kc.at(i) = Kc_ref * Arrhenius(i, Eac);//Eq 9
		delta.at(i) = delta_ref * Arrhenius(i, Ead);//gamma0* (1.+ gamma1*(TleafK.at(i)- Tref) + gamma2*std::pow((TleafK.at(i) - Tref),2.) ) ;//Eq 10
		
		//compute Rd before deltagco2
		//std::cout<<"Photosynthesis::initVcVjRd "<<PhotoType<<" "<<C3<<" "<<C4<<" "<<(PhotoType == C3)<<" "<<(PhotoType == C4)<<std::endl;
		switch(PhotoType)
		{
			case C3: {photoC3_init(i);break;}		
			case C4: {photoC4_init(i);break;}	
			default:{throw std::runtime_error("Photosynthesis::initVcVjRd: PhotoType not recognised  ");}
		}


		deltagco2.at(i) = std::max(0.,(delta.at(i) + Kc.at(i)*Rd.at(i)*(1. + oi/Ko.at(i))/Vcmax.at(i))/(1-Rd.at(i)/Vcmax.at(i)));
	}
}


	/*
		add gravitational wat. pot to total wat. pot. (used in phloem module)
	*/
void Photosynthesis::doAddGravity()
{
	psiXyl4Phloem.resize(psiXyl.size(), 0.);
	for(int i = 0; i<  psiXyl.size(); i++)
	{
		psiXyl4Phloem.at(i) =psiXyl.at(i) +  plant->nodes.at(i).z ;//in cm
	}

}

double Photosynthesis::thermalBreakdown(int index, double Ed)
{
	 //only evaluate denominator as the nominator is ~ 1
	double TleafK_ =  getMeanOrSegData(TleafK, index);	
	return 1/(1 + std::exp((S * TleafK_ - Ed)/(R_ph *0.1* TleafK_)));
}

double Photosynthesis::Arrhenius(int index, double Ea)
{
	return std::exp(Ea /(R_ph*0.1*Tref)*(1. - Tref/getMeanOrSegData(TleafK, index)));
}




void Photosynthesis::photoC4_init(int i)
{
	Rd_ref.at(i) = 0.025 * Vcrefmax.at(i); //Bonan2019Chap11

	//Vcmax
	double TleafK_ =  getMeanOrSegData(TleafK, i);	
	Vcmax.at(i) = Vcrefmax.at(i) * Q10f(i) / (1 + std::exp(s1 * ( TleafK_ - s2) ) ) / (1 + std::exp(s3 * ( s4 - TleafK_) ) ); //Bonan2019Chap11
	Rd.at(i) = Rd_ref.at(i) * Q10f(i) / (1 + std::exp(s5 * ( TleafK_ - s6) ) );
	kp25.at(i) = 0.02 * Vcrefmax.at(i); //from Bonan2019Chap11
	kp.at(i) = kp25.at(i) * Q10f(i);
		
}

void Photosynthesis::photoC3_init(int i)
{		
	Rd_ref.at(i) = 0.015 * Vcrefmax.at(i); //Bonan2019Chap11
	
	//Vcmax
	Vcmax.at(i) = Vcrefmax.at(i) * Arrhenius(i, Eav) * thermalBreakdown(i, Edv); //Eq 11
	Rd.at(i) = Rd_ref.at(i) * Arrhenius(i, Eard);// * thermalBreakdown(i, Edrd);
	
	//electron transport rate
	//Jrefmax
	Jrefmax.at(i) = Vcrefmax.at(i) * a3 ;//Eq 25
	//Jmax
	Jmax.at(i) = std::min(Jrefmax.at(i) * Arrhenius(i, Eaj) * thermalBreakdown(i, Edv), Jrefmax.at(i)); //Eq 24
	//J
	double Qlight_ = getMeanOrSegData(Qlight, i);
	double coefa = theta;
	double coefb = -(alpha * Qlight_ + Jmax.at(i));
	double coefc = alpha * Qlight_ * Jmax.at(i);
	double dis = std::pow(coefb,2.) - (4.*coefa*coefc);
	if (dis < 0) {
			throw std::runtime_error("Photosynthesis::initVcVjRd : root for J not found");
	  }
	J.at(i) =  ((-coefb- std::sqrt(dis))/(2.*coefa));//rostamza2020, Bonan2019Chap11
	if (J.at(i) < 0) {
			throw std::runtime_error("Photosynthesis::loopCalcs : J < 0");
	  }
	
}


void Photosynthesis::photoC3_loop( int i)
{
	//carboxylation and electron transport  rate
	Vc.at(i) = std::max(std::min(std::max(Vcmax.at(i) * (ci.at(i) - delta.at(i)) / (ci.at(i) + Kc.at(i)*(1. + oi/Ko.at(i))),
						0.),
						Vcmax.at(i)),0.); //Eq 8

	Vj.at(i) = std::max(J.at(i)/4. * (ci.at(i) - delta.at(i))/ (ci.at(i) + 2. * delta.at(i)), 0.) ;//Eq 22

	//An mol m-2 s-1
	An.at(i) = std::min(Vc.at(i), Vj.at(i)) - Rd.at(i);//Eq 6
	
}


void Photosynthesis::photoC4_loop(int i)
{
	Vc.at(i) = Vcmax.at(i);
	double Qlight_ = getMeanOrSegData(Qlight, i);
	
	Vj.at(i) = alpha * Qlight_;
	//std::cout<<"Photosynthesis::photoC4_loop "<<Vj.at(i)<<" "<< alpha <<" "<< Qlight_<<std::endl;
	
	Vp.at(i) = kp.at(i) * (std::max(ci.at(i) - deltagco2.at(i),0.)* 1e6);//[mol CO2 m-2 s-1] * [mumol mol-1]
	//An mol m-2 s-1
	An.at(i) = std::min(std::min(Vc.at(i), Vj.at(i)), Vp.at(i)) - Rd.at(i);//Eq 6
}



	/*
		Computes the output variables => ci, go2, An, Ev
		@param simtime
	*/
void Photosynthesis::loopCalcs(double simTime, std::vector<double> sxx_, bool cells_){
	std::ofstream myfile4;
	if(doLog)
	{
		std::string name1 = "loopphoto";
		std::string name2 = ".txt";
		std::string namefile = name1 +"_"+std::to_string(simTime)+"_"+std::to_string(loop) + name2;
		myfile4.open (outputDir+namefile);
	}
	for(int i = 0; i<seg_leaves_idx.size();i++)
	{
		int idl= seg_leaves_idx.at(i); // global segment indx
		
		double TleafK_ =  getMeanOrSegData(TleafK, i);	
		double cs_ = getMeanOrSegData(cs, i);	
		double g_bl_ = getMeanOrSegData(g_bl, i);	
		double g_canopy_ = getMeanOrSegData(g_canopy, i);	
		double g_air_ = getMeanOrSegData(g_air, i);	
		
		double ea_;
		
		if(cells_){ea_ = ea; //mean air ea
		}else{ // TODO: change how that is handled
			double psi_air_ = sxx_.at(idl);
			ea_ = std::exp(psi_air_/(rho_h2o * R_ph * TleafK_/Mh2o * (1/0.9806806) ))*this->es  ; //in cm
			
		}
		
		
		if((verbose_photosynthesis ==2)){std::cout<<"in loopcalcs "<<i<<" "<<idl<<std::endl;}
		double l = lengths.at(idl);
		
		double sideArea = plant->leafBladeSurface.at(idl)*2;//
		
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
        auto n1 = plant->nodes[plant->segments.at(idl).x].z;
        auto n2 = plant->nodes[plant->segments.at(idl).y].z;
		//(mg mmol-1)* hPa /((hPa cm3K−1mmol−1) mg cm-3 K) =(-)
		double HRleaf = std::exp(Mh2o*(this->pg.at(i) + (n1+n2)/2)*0.9806806 /(rho_h2o*R_ph*TleafK_)) ;//fractional relative humidity in the intercellular spaces
		//double ea = es - VPD;
		double ea_leaf = es * HRleaf;//hPa
		if((verbose_photosynthesis ==2)){std::cout<<"git to leaf "<<ea_leaf<<" "<<HRleaf<<std::endl;}
		if(std::abs(fv.at(i)) > 1e-16)//i.e., perimeter * kr > 1e-16 like for @see Xylem::solveLinear
		{
            double p_lhPa = this->pg.at(i)*0.9806806;// cm => hPa +(n1+n2)/2
            
			fw.at(i) = fwr +std::max(0., (1.- fwr)*(1+std::exp(sh*p_lcrit))/(1+std::exp(sh*(p_lcrit-p_lhPa))) - fw_cutoff);
            
			ci.at(i) = std::max((cs_*a1*fw.at(i) +deltagco2.at(i))/(1+a1* fw.at(i)),deltagco2.at(i)) ;
			if((verbose_photosynthesis ==2)){std::cout<<"in compute gco2 "<<sideArea<<" "<<(sideArea > 1e-16)<<std::endl;}
			
			switch(PhotoType)
			{
				case C3: {photoC3_loop(i);break;}		
				case C4: {photoC4_loop(i);break;}	
				default:{throw std::runtime_error("Photosynthesis::loopCalcs: PhotoType not recognised");}
			}
			
			// mol CO2 m-2 s-1
			gco2.at(i) = std::max(g0 + fw.at(i) * a1 *( An.at(i) + Rd.at(i))/(ci.at(i) - deltagco2.at(i)),0.);//tuzet2003
			// mol H2O m-2 s-1 MPa-1
			
			//(mol m-2 s-1)*(mmol/mol)*(hPa/hPa) * (mg mmol-1) /(mg cm-3) *(h/d)*(s/h)*(m2 m-2) =  ( cm3)/d*(cm-2)
				if((gco2.at(i)>0.)&&(gm*fw.at(i)>0.)&&(g_bl_>0.)&&(g_canopy_>0.)&&(g_air_>0.))//
            {//
                gtotOx.at(i) = 1/(1/(gco2.at(i) * a2_stomata) +1/(gm*fw.at(i) )+ 1/(g_bl_ * a2_bl) + 1/(g_canopy_ * a2_canopy) + 1/(g_air_ * a2_air) );
            }else
            {
                gtotOx.at(i) = 0.;
            }																							  
			Jw.at(i) = gtotOx.at(i) *1000* (ea_leaf - ea_)/Patm * Mh2o/rho_h2o * 24.*3600*1e-4 ;//in cm3 cm-2 d-1
			Ev.at(i) = Jw.at(i)* sideArea; //in cm3 d-1
            PVD.at(i) =  ea_leaf - ea_ ;
            EAL.at(i) =  ea_leaf;
            hrelL.at(i) =  HRleaf;						 
            //double f = -2*a*M_PI*kr; // flux is proportional to f // *rho*g

			//ci and pg
			//gruard cell wat. pot. to havee water flux from xylem to gard cell. kr = permeability of xylem membrane only.
			this->pg.at(i) = (-1/2.)*((Ev.at(i))/(-fv.at(i)*(1./(tauv.at(i)*dv.at(i)))
				*(2.-std::exp(-tauv.at(i)*l)-std::exp(tauv.at(i)*l))) - (rxi + rxj)) ;//cm

			if((verbose_photosynthesis ==2)){

				std::cout<<"sizes "<<An.size()<<" "<< gco2.size()<<" "<<ci.size()<<" "<<ci_old.size() <<std::endl;

			}
			bool erroHappened = (!std::isfinite(this->pg.at(i)))||(!std::isfinite(ci.at(i)))||(ci.at(i)<0)||(fw.at(i)-1>1e-10)||(fw.at(i)<0);
			if(erroHappened) {
			std::cout<<"shape leaf "<<idl<<" "<<sideArea<<" "<<ci_old.at(i)<<" "<<ci.at(i)<<std::endl<<std::flush;
			std::cout<<"an calc "<<An.at(i)<<" "<<Vc.at(i)<<" "<< Vj.at(i)<<" "<<J.at(i)<<" "<<Vcmax.at(i)<<" "<<Kc.at(i)<<" "<<Ko.at(i)<<" "<<std::flush;
			std::cout<<" "<<delta.at(i)<<" "<<oi<<" "<<eps<<std::endl<<std::flush;
			std::cout<<"forgco2 "<<gco2.at(i) <<" "<< g0<<" "<< ( fw.at(i) -1)<<" "<<  fw.at(i) <<" "<<  a1 <<" "<< An.at(i)<<" "<< Rd.at(i)<<" "<< deltagco2.at(i)<<std::endl<<std::flush;
			std::cout<<"forJW, Jw "<<Jw.at(i)<<" drout_in "<<(this->pg.at(i) - (rxi + rxj)/2)<<" "<<ea_leaf <<" "<< ea<<" "<<Patm<<" "<<Mh2o<<" "<<rho_h2o <<std::endl;
			std::cout<<"forpg "<<sideArea<<" "<<fv.at(i)<<" "<<tauv.at(i)<<" "<<dv.at(i)<<" "<<l<<" "<<rxi <<" "<< rxj<<" "<<this->pg.at(i)<<" numleaf: "<<i <<std::endl<<std::flush;
			std::cout<<"diff Ev and lat fluw: "<<Ev.at(i)<<std::endl<<std::flush;//<<" "<<outputFluxL.at(i)
			std::cout<<"cause of the error: "<<(!std::isfinite(this->pg.at(i)))<<" "<< (!std::isfinite(ci.at(i)))<<" "<<(ci.at(i)<0)<<" "<< (fw.at(i)>1)<<" "<< (fw.at(i)<0)<<std::endl<<std::flush;
                
				throw std::runtime_error("Phtotosynthesis: nan or Inf  pg.at(i)");
			}

		}else{ci.at(i) = 0.0;}
		if(doLog ){
			if((verbose_photosynthesis ==2)){std::cout<<"dolog loop calcs "<<std::endl;}
			myfile4<<"shape leaf "<<idl<<" "<<sideArea<<" "<<ci_old.at(i)<<" "<<ci.at(i)<<std::endl;
			myfile4<<"an calc "<<An.at(i)<<" "<<Vc.at(i)<<" "<< Vj.at(i)<<" "<<J.at(i)<<" "<<Vcmax.at(i)<<" "<<Kc.at(i)<<" "<<Ko.at(i)<<" ";
			myfile4<<" "<<delta.at(i)<<" "<<oi<<" "<<eps<<std::endl;
			myfile4<<"forgco2 "<<gco2.at(i) <<" "<< g0<<" "<<  fw.at(i) <<" "<<  a1 <<" "<< An.at(i)<<" "<< Rd.at(i)<<" "<< deltagco2.at(i)<<std::endl;
			myfile4<<"forJW, Jw "<<Jw.at(i)<<" drout_in "<<(this->pg.at(i) - (rxi + rxj)/2)<<" "<<ea_leaf <<" "<< ea<<" "<<Patm<<" "<<Mh2o<<" "<<rho_h2o <<std::endl;
			myfile4<<"forpg "<<sideArea <<" "<<fv.at(i)<<" "<<tauv.at(i)<<" "<<dv.at(i)<<" "<<l<<" "<<rxi <<" "<< rxj<<" "<<this->pg.at(i)<<" numleaf: "<<i <<std::endl;
			myfile4<<"diff Ev and lat fluw: "<<Ev.at(i)<<std::endl;//<<" "<<outputFluxL.at(i)
		}
	}
	if(doLog){myfile4.close();}

}



double Photosynthesis::kr_f(int si, double age, int subType, int organType)
{
	//try{
		//std::cout<<"Photosynthesis::kr_f "<<si<<" "<<subType<<" "<< organType<<" "<<plant->seg2cell.size()<<std::endl;
		int cellIndex = plant->seg2cell.at(si);
		if (((cellIndex>=0)&&(organType ==Organism::ot_leaf))||((cellIndex < 0)&&(organType ==Organism::ot_root)))
		{ 
				return 0.;
		}else
		{ 
			//std::cout<<"to XylemFlux::kr_f"<<std::endl;
			return params->kr_f(si, age, subType, organType);
		}
	//} catch(...) { 
	//	return XylemFlux::kr_f(si, age, subType, organType);
	//}
	
}

}//end namespace
