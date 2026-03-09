#include <PiafMunch/CPB_to_PM.h>
#include "runPM.h"

	/*water limited deltaSuc per node*/
void CPB_to_PM::organToNodeData(double dt) //waterLimitedGrowth
{
	// initialise data
	Nt = plant->nodes.size();
	std::map<int,double> initMap = { { -1,0. } };//-1 : total need per node
	deltaSucOrgNode = std::vector<std::map<int,double>>(Nt, initMap);//, 0.);
	bool allOrgs = true;
	auto orgs = plant->getOrgans(-1,allOrgs);//also org with length - Epsilon == 0	
	BackUpMaxGrowth = {};//std::vector<double>(orgs.size(), 0.); //for checking in @see CPB_to_PM::computeOrgGrowth
	Flen = std::vector<double>(Nt,0.);
	AuxinSource = std::vector<double>(Nt , 0.)  ;//reset


	// start computation
	computeFpsi();
	for(auto org: orgs)
	{
		int ot = org->organType();
		int st = org->getParameter("subType");
		assert(org != nullptr);	
		std::vector<int> localGrowingNodesId = getGrowingNodes(org); // localIds	
		double rmax, Lmax;
		getRmaxLmax_st_f(org,st,ot, rmax, Lmax);	//TODO: change, does not work to have a diff rmax	
		double deltaVol = getMaxVolumicGrowth(org, dt, rmax, Lmax); // cm3				
		computeFlen(org, localGrowingNodesId);
		waterLimitedGrowth(org, localGrowingNodesId, dt, deltaVol);
		setAuxinSource(org);		
	}	
}


void CPB_to_PM::setAuxinSource(std::shared_ptr<CPlantBox::Organ> org)
{
	if(org->isActive()) // isActive False in Lmax or decapitated
	{
		int youngExpanding = (org->getLength(false) < org->getParameter("k")*limLenActive);//TODO: check
		if((org->organType()==4)&&(youngExpanding)) // base of young expading leaves
		{
			int nodeSourceId = org->getNodeId(0);//= global pni
			AuxinSource.at(nodeSourceId) += 1; //can carry several expanding leaves == several sources

		}
		if(org->organType()==3)// tip of active stems
		{
			AuxinSource.at(org->getNodeId(org->getNumberOfNodes()-1)) += 1;
		}
	}                
}


void CPB_to_PM::computeAuxinValue(std::shared_ptr<CPlantBox::Organ> org)
{
	double externalAuxin;
	if(org->getParent()->getNumberOfNodes() > (org->parentNI+1) )
	{
		externalAuxin = piafmunch->C_Auxinv.at(org->getParent()->getNodeId(org->parentNI+1) );//concentration above
	}else{
		if(org->getNumberOfNodes() >1)//concentration base - concentration next to it
		{
			externalAuxin = piafmunch->C_Auxinv.at(org->getNodeId(0) ) -  piafmunch->C_Auxinv.at(org->getNodeId(1) ) ;
		}else{
			externalAuxin = piafmunch->C_Auxinv.at(org->getNodeId(0) )  ;//no solutions
		}
	}
	if(org->auxTested > externalAuxin) // theory of canal created and maintained
	{
		org->setAuxTested(externalAuxin);
	}
}

void CPB_to_PM::updateBudStage(double EndTime)// TODO: check at each time step
{
	bool allOrgs = true;
	auto orgs = plant->getOrgans(3,allOrgs);//also org with length - Epsilon == 0
    for(auto org: orgs)
	{
		if(org->getBudStage() != 2)//org->getOrganRandomParameter()->isBudType)
		{
			int bu_bs = org->getBudStage();
			//use the org_mean instead?
			double suc =  piafmunch->C_STv.at(org->getNodeId(org->getNumberOfNodes()-1)) ;
			switch(org->getBudStage())
			{
				case 1://active bud
				{
					computeAuxinValue(org);                        

					//concentration in node above (not affected by own outflux)
					double lenFact = org->getLength(false); // TODO: update to change directly kx of bud
					double suc_ =  suc * lenFact;// go from concentration to somewhat content ///plant->maxLBud);

						
					org->setBerthFact(computeBerth(suc_, org->auxTested));
					
					if(org->BerthFact <= BerthLim){org->setBudStage(2);}//congrats! you are a branch
					if(org->BerthFact > L_dead_threshold){org->setBudStage(-1);}//sorry! you are dead
					org->setSucTested(suc);
					break;
				}
				case 0://dormant
				{
					
						computeAuxinValue(org); 
						
					double lenFact = org->getLength(false);
					double suc_ =  suc * lenFact;// * (org->getLength(false)/plant->maxLBudDormant);//content and not concentration is important
						org->setSucTested(suc);
						
					if(suc_ >= CSTthreshold)
					{
						org->setBudStage(1);
						org->setAge(0.); //reset age
						
					}//congrats! you are released
					break;
				}

			}
			if(bu_bs != org->getBudStage())
			{ org->setBudStageChange(org->getBudStage(),EndTime);}
		
		}
    }
}

double CPB_to_PM::adaptDt(std::shared_ptr<CPlantBox::Organ> org, double dt)
{
	double age = org->getAge();
	int ot = org->organType(); 
	double orgLT = org->getParameter("rlt");
	
	if (age+dt>orgLT) { // root life time
		dt=orgLT-age; // remaining life span
	}
	
	//no probabilistic branching models and no other scaling via getRootRandomParameter()->f_se->getValue(nodes.back(), shared_from_this());
	if(ot == CPlantBox::Organism::ot_stem){
		
		double delayNGStart = org->getParameter("delayNGStart");
		double delayNGEnd = org->getParameter("delayNGEnd");
		
		if((age+dt) > delayNGStart){//simulation ends after start of growth pause
			if((age+dt)  < delayNGEnd){dt = 0;//during growth pause
			}else{
				double beforePause = 0.; double afterPause =0.;
				if(age < delayNGStart){
					beforePause = std::min(dt,delayNGStart - age); //reaches start of growth pause and then stopes
				}
				if((age+dt)  > delayNGEnd){
					afterPause = std::min(dt, std::max(age +dt  - delayNGEnd, 0.)); //part of the simulation after end of pause
				}
				dt = beforePause + afterPause;//part of growth pause during simulation				
			}
		}
	}
	if(dt < 0){
		throw std::runtime_error("CPB_to_PM::adaptDt: dt <0");
	}
	return dt;
	
}


void CPB_to_PM::assertUsedCReserves(std::shared_ptr<CPlantBox::Organ> org)
{
	if((!((org->getOrganRandomParameter()->f_gf->CW_Gr.empty()) || 
					(org->getOrganRandomParameter()->f_gf->CW_Gr.count(org->getId()) ==0) ||
					(org->getOrganRandomParameter()->f_gf->CW_Gr.find(org->getId())->second<0.)))&&
					org->isActive())
	{
		std::cout<<org->getId()<<" "<<org->getOrganRandomParameter()->f_gf->CW_Gr.find(org->getId())->second<<std::endl;
		std::cout<<org->calcLength(1)<<" "<< org->organType() <<" "<<org->getAge()<<std::endl;
		throw std::runtime_error("CPB_to_PM::assertUsedCReserves: sucrose for growth has not been used at last time step");
	}
}

void CPB_to_PM::getRmaxLmax_st_f(std::shared_ptr<CPlantBox::Organ> org, int st, int ot, double & rmax, double & Lmax)
{
    if(org->getBudStage() != 2)//org->getParameter("isBudType"))
	{		
		double maxLBudDormant_ = maxLBudDormant.at(st);// todo check
		double maxLBud_ = maxLBud.at(st);// todo check
		switch(org->getBudStage()) 
		{
				case -1:{rmax =  Rmax_st_f(st,ot); Lmax = org->getLength(false);break;}
				case 0:{rmax = budGR;//1 mm/d
						Lmax = maxLBudDormant_; 
						break;}
				case 1 :{rmax = budGR;Lmax = maxLBud_;break;}//1 mm/d
				//case 2 :{rmax =  Rmax_st_f(st,ot);//org->getParameter("r");
				//		 Lmax = org->getParameter("k");break;}//1 mm/d
				default:{std::cout<<"stem::simulate: budStage not recognised "<< org->getBudStage() <<std::flush;
						assert(false);}
		}
	}else{
		rmax =  Rmax_st_f(st,ot);//org->getParameter("r");
		Lmax = org->getParameter("k");
	}
}


double CPB_to_PM::getMaxVolumicGrowth(std::shared_ptr<CPlantBox::Organ> org, double t, double rmax, double Lmax)
{
	double Linit = org->getLength(false);//theoretical length 
			
	if((org->getParameter("gf") != 3)&&(org->getParameter("gf") != 4)&&do_gf_warning)
	{
		std::cout<<"CPB_to_PM::getMaxVolumicGrowth: organ(s) do(es) not use carbon-limited growth"<<std::endl;
		do_gf_warning = false;
	}
	auto f_gf = plant->createGrowthFunction(org->getParameter("gf"));
	if(org->getParameter("gf") == 4)
	{
		f_gf =  plant->createGrowthFunction(2);
	}
	if(org->getParameter("gf") == 3)
	{
		f_gf =  plant->createGrowthFunction(1);
	}
		
	double age_ = f_gf->getAge(Linit, rmax, Lmax, org->shared_from_this());			
	assertUsedCReserves(org);			
	double dt = adaptDt(org, t);				
	//	params to compute growth
	//double LinitTemp = f_gf->getLength(age_ , rmax, Lmax, org->shared_from_this());
	double targetlength = f_gf->getLength(age_ + dt , rmax, org->getParameter("k"), org->shared_from_this());
	double e = std::max(0.,targetlength-Linit);// unimpeded elongation in time step dt
	
	if((e + org->getLength(false))> org->getParameter("k") + 1e-10){
		std::cout<<"CPB_to_PM::getMaxLengthGrowth: target length too high "<<e<<" "<<dt<<" "<<Linit;
		std::cout<<" "<<org->getParameter("k")<<" "<<org->getId()<<std::endl;
		assert(false);
	}
	BackUpMaxGrowth[org->getId()] = Linit + e; // to compare with final growth length in @see CPB_to_PM::computeOrgGrowth
	double deltavol = std::max(0.,org->orgVolume(Linit + e, false) - org->orgVolume(Linit, false));//volume from theoretical length
	return deltavol;
}

/*
	//roots grow at their tip
	//leaves grow in their groth zones
	//stems grow at each phytomere nodes
*/

std::vector<int> CPB_to_PM::getGrowingNodes(std::shared_ptr<CPlantBox::Organ> org)
{
	int ot = org->organType(); 
	int nNodes = org->getNumberOfNodes();
	std::vector<int> nodeLocalIds;

	if ((nNodes==1)||(ot == 2)) { //organ not represented because below dx limit or is root
		nodeLocalIds.push_back(nNodes - 1);
	}else{
		if(ot == 4)//only for nodes in the leaf growth zone
		{
			for(int k = 0; k < nNodes; k++)
			{
				nodeLocalIds.push_back(k);
				if(org->getLength(k) > this->leafGrowthZone)
				{
					break;
				}
			}
		}
		if(ot == 3)
		{
			int nn1 = 0;
			if((StemGrowthPerPhytomer)&&(org->getNumberOfChildren() > 0))
			{
				// only for nodes in the growing phytomers. Phytomer ends each time we have a lateral.
				auto stemParams = std::static_pointer_cast<CPlantBox::Stem>(org)->param();
				bool foundPhytoIdx = false;
				int PhytoIdx = 0; // to remember PhytoIdx outside of loop
				for(; ((PhytoIdx < stemParams->ln.size())&(!foundPhytoIdx));PhytoIdx++)
				{
					double maxPhytoLen = stemParams->ln.at(PhytoIdx);
					double currentPhytoLen = org->getLength(org->getChild(PhytoIdx+1)->parentNI) - org->getLength(org->getChild(std::max(0,PhytoIdx))->parentNI);
					foundPhytoIdx = (maxPhytoLen - currentPhytoLen > 1e-10);//first still growing phytomere
					if(maxPhytoLen - currentPhytoLen < -1e-10){throw std::runtime_error("maxPhytoLen - currentPhytoLen < -1e10;");}
					
				}
				assert((nn1 < nNodes)&&"nn1 > nNodes");
			}
			for (int i = nn1; i < nNodes; ++i) {
				nodeLocalIds.push_back(i);
			}
		}
	}
	return nodeLocalIds;
}

void CPB_to_PM::computeFpsi()
{
	Fpsi = std::vector<double>(Nt,1.); 
	psi_p_symplasm = std::vector<double>(Nt,1.); 
	for(std::size_t nodeId = 0; nodeId < piafmunch->psiXyl.size(); nodeId++)
	{
		psi_p_symplasm.at(nodeId) =  piafmunch->psiXyl.at(nodeId) - psi_osmo_proto;
		if(-psi_osmo_proto - psiMin>1e-5){
			Fpsi.at(nodeId) = std::max((std::max(psi_p_symplasm.at(nodeId), psiMin) - psiMin)/(-psi_osmo_proto - psiMin),0.);
		}else{Fpsi.at(nodeId) = 0.;}
		
		if((Fpsi.at(nodeId) > (1+1e-6))||(Fpsi.at(nodeId) < -1e-6))
		{
			throw std::runtime_error("Fpsi < 0 or Fpsi > 1");
		}
	}
}

void CPB_to_PM::computeFlen(std::shared_ptr<CPlantBox::Organ> org, std::vector<int> growingNodesId)
{
	double Flen_tot = 0.;
	if (growingNodesId.size() == 1)
	{
		Flen.at(org->getNodeId(growingNodesId.at(0))) = 1.;
	}else{
		double L_n0 = org->getLength(growingNodesId.at(0)); 
		double L_growth_length = org->getLength(growingNodesId.at(growingNodesId.size() - 1)) - L_n0;
		for(std::size_t k = 1; k < growingNodesId.size(); k++)//int nodeId : growingNodesId)
		{
			int nodeGlobalId = org->getNodeId(growingNodesId.at(k));
			auto nodei = plant->nodes.at(org->getNodeId(growingNodesId.at(k - 1)));
			auto nodej = plant->nodes.at(nodeGlobalId);
			double Lseg = nodej.minus(nodei).length();
			Flen.at(nodeGlobalId) = (Lseg/L_growth_length);
			Flen_tot += Flen.at(nodeGlobalId);
		}
		if((std::abs(Flen_tot - 1.)>1e-10)||(std::abs(Flen_tot - 1.)<-1e-10))
		{
			throw std::runtime_error("CPB_to_PM::computeFlen: (Flen_tot != 1.)&&wrong tot Flen");
		}
	}
}

void CPB_to_PM::waterLimitedGrowth(std::shared_ptr<CPlantBox::Organ> org, 
							std::vector<int> localGrowingNodesId, double t, double deltaVol)
{
	double dt = t; //copy
	if(((org->getAge()+dt)>0)&&(org->isAlive())&&(org->isActive()))
	{ //organ alive and active at this time step
		int ot = org->organType(); 
		int st = org->getParameter("subType");
		double deltaSucGrowth_tot = 0.;
		double rhoSucrose_double = rhoSucrose_f(st,ot);
		for(int localNodeId : localGrowingNodesId) // fill out the sink power at each node
		{
			int nodeId = org->getNodeId(localNodeId);
			double deltaSucGrowth_per_node = deltaVol * Flen.at(nodeId) * Fpsi.at(nodeId) * rhoSucrose_double;
			
			if((deltaSucGrowth_per_node<0.)||(deltaSucGrowth_per_node != deltaSucGrowth_per_node)){
				// could be error of pressision (if l = Lmax)
				// or that, because of nodal growth and dxMin, org->getEpsilon() <0
				if(((deltaSucGrowth_per_node>-1e-5)||(deltaVol>0.))
						&&(deltaSucGrowth_per_node == deltaSucGrowth_per_node)){
					deltaSucGrowth_per_node=0.;	//within margin of error
				}else{
					throw std::runtime_error("(deltaSucGrowth_per_node<0.)||(deltaSucGrowth_per_node != deltaSucGrowth_per_node)");
				}					
			}				
			assert(nodeId >= 0 && nodeId < (int)deltaSucOrgNode.size());
			deltaSucOrgNode.at(nodeId)[org->getId()] = deltaSucGrowth_per_node;
			deltaSucOrgNode.at(nodeId).at(-1) += deltaSucGrowth_per_node;								
			deltaSucGrowth_tot += deltaSucGrowth_per_node;
		}
		assert(((deltaSucGrowth_tot * rhoSucrose_double - deltaVol)<1e-10)&&"deltaSucGrowth_tot too high");//deltaSucGrowth_tot * rhoSucrose_double <= deltaVol
	}else{
		BackUpMaxGrowth[org->getId()] = org->getLength(false);
	}
}


/**
 *  Sets the radial conductivity in [1 day-1]
 * in case of organ_type specific kr
 * @param values 		kr per organ pr/and organ type) or/and per age [cm-1]
 * @param age 			ages if kr per age
 * @param kr_length_ 	exchange zone in root, where kr > 0 [cm from root tip], default = -1.0, i.e., no kr_length
 */
//type/subtype dependent
void CPB_to_PM::setKr_st(std::vector<std::vector<double>> values, double kr_length_, bool verbose) {
    kr_st = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			kr_st_f = std::bind(&CPB_to_PM::kr_st_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
            if (verbose)
            {
                std::cout << "Kr_st is constant " << values[0][0] << " 1 day-1 \n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			kr_st_f = std::bind(&CPB_to_PM::kr_st_perOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
			if (verbose)
            {
                std::cout << "Kr_st is constant per organ type, organ type 2 (root) = " << values[0][0] << " 1 day-1 \n";
            }
		} else {
			if (verbose)
            {
                std::cout << "Exchange zone in roots: kr_st > 0 until "<< kr_length_<<"cm from root tip "<<(kr_length_ > 0)<<" "<<(kr_length_ > 0.)<<std::endl;
            }
			if(kr_length_ > 0.){
				if (verbose)
            {
                std::cout << "Exchange zone in roots: kr > 0 until "<< kr_length_<<"cm from root tip"<<std::endl;
                }
				plant->kr_length = kr_length_; //in MappedPlant. define distance to root tipe where kr > 0 as cannot compute distance from age in case of carbon-limited growth
				plant->calcExchangeZoneCoefs();	
				kr_st_f  = std::bind(&CPB_to_PM::kr_st_RootExchangeZonePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
			}else{
				kr_st_f  = std::bind(&CPB_to_PM::kr_st_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
			}
			if (verbose)
            {
                std::cout << "Kr_st is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << " 1 day-1 \n";
            }
		}
	}
    
}	
// type/subtype dependent
void CPB_to_PM::setKx_st(std::vector<std::vector<double>> values, bool verbose) {
    kx_st = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			kx_st_f = std::bind(&CPB_to_PM::kx_st_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Kx_st is constant " << values[0][0] << " cm3 day-1 \n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			kx_st_f = std::bind(&CPB_to_PM::kx_st_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Kx_st is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm3 day-1 \n";
            }
		} else {
			kx_st_f  = std::bind(&CPB_to_PM::kx_st_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Kx_st is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << " cm3 day-1 \n";
            }
		}
	}
}	


//type/subtype dependent
void CPB_to_PM::setAcross_st(std::vector<std::vector<double>> values, bool verbose) {
    Across_st = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			Across_st_f = std::bind(&CPB_to_PM::Across_st_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Across_st is constant " << values[0][0] << " cm2\n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			Across_st_f = std::bind(&CPB_to_PM::Across_st_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Across_st is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm2 \n";
            }
		} else {
			Across_st_f  = std::bind(&CPB_to_PM::Across_st_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Across_st is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << " cm2 \n";
            }
		}
	}
}	


//type/subtype dependent ==> to compute Rhat Fhat
void CPB_to_PM::setPerimeter_st(std::vector<std::vector<double>> values, bool verbose) {
    Perimeter_st = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			Perimeter_st_f = std::bind(&CPB_to_PM::Perimeter_st_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Perimeter_st is constant " << values[0][0] << " cm\n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			Perimeter_st_f = std::bind(&CPB_to_PM::Perimeter_st_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Perimeter_st is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm\n";
            }
		} else {
			Perimeter_st_f  = std::bind(&CPB_to_PM::Perimeter_st_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Perimeter_st is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << " cm\n";
            }
		}
	}
}	

// type/subtype dependent
void CPB_to_PM::setRmax_st(std::vector<std::vector<double>> values, bool verbose) {
    Rmax_st = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			Rmax_st_f = std::bind(&CPB_to_PM::Rmax_st_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Rmax_st is constant " << values[0][0] << " cm day-1 \n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			Rmax_st_f = std::bind(&CPB_to_PM::Rmax_st_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Rmax_st is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm day-1 \n";
            }
		} else {
			Rmax_st_f  = std::bind(&CPB_to_PM::Rmax_st_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "Rmax_st is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << " cm day-1 \n";
            }
		}
	}
}	


//either age or type/subtype dependent
void CPB_to_PM::setRhoSucrose(std::vector<std::vector<double>> values, bool verbose) {
    rhoSucrose = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			rhoSucrose_f = std::bind(&CPB_to_PM::rhoSucrose_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "rhoSucrose is constant " << values[0][0] << " mmol cm-3\n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			rhoSucrose_f = std::bind(&CPB_to_PM::rhoSucrose_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "rhoSucrose is constant per organ type, organ type 2 (root) = " << values[0][0] << " mmol cm-3\n";}
		} else {
			rhoSucrose_f  = std::bind(&CPB_to_PM::rhoSucrose_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "rhoSucrose is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << " mmol cm-3\n";}
		}
	}
}	


//either age or type/subtype dependent
void CPB_to_PM::setKrm1(std::vector<std::vector<double>> values, bool verbose) {
    krm1v = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			krm1_f = std::bind(&CPB_to_PM::krm1_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "krm1 is constant " << values[0][0] << " -\n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			krm1_f = std::bind(&CPB_to_PM::krm1_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "krm1 is constant per organ type, organ type 2 (root) = " << values[0][0] << " -\n";
            }
		} else {
			krm1_f  = std::bind(&CPB_to_PM::krm1_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "krm1 is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << "-\n";
            }
		}
	}
}	

//either age or type/subtype dependent
void CPB_to_PM::setKrm2(std::vector<std::vector<double>> values, bool verbose) {
    krm2v = values;
	if (values.size()==1) {
		if (values[0].size()==1) {
			krm2_f = std::bind(&CPB_to_PM::krm2_const, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "krm2 is constant " << values[0][0] << " -\n";
            }
		} 
	} else {
		if (values[0].size()==1) {
			krm2_f = std::bind(&CPB_to_PM::krm2_perOrgType, this, std::placeholders::_1, std::placeholders::_2);
			std::cout << "krm2 is constant per organ type, organ type 2 (root) = " << values[0][0] << " -\n";
		} else {
			krm2_f  = std::bind(&CPB_to_PM::krm2_perType, this, std::placeholders::_1, std::placeholders::_2);
			if (verbose)
            {
                std::cout << "krm2 is constant per subtype of organ type, for root, subtype 1 = " << values[0].at(1) << "-\n";
            }
		}
	}
}	