#include "Stem.h"

#include "Leaf.h"
#include "Root.h"
#include "Plant.h"

namespace CPlantBox {

std::vector<int> Stem::phytomerId = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/**
 * Constructs a root from given data.
 * The organ tree must be created, @see Organ::setPlant, Organ::setParent, Organ::addChild
 * Organ geometry must be created, @see Organ::addNode, ensure that this->getNodeId(0) == parent->getNodeId(pni)
 *
 * @param id        		the organ's unique id (@see Organ::getId)
 * @param param     		the organs parameters set, ownership transfers to the organ
 * @param alive     		indicates if the organ is alive (@see Organ::isAlive)
 * @param active    		indicates if the organ is active (@see Organ::isActive)
 * @param age       		the current age of the organ (@see Organ::getAge)
 * @param length    		the current length of the organ (@see Organ::getLength)
 * @param partialIHeading 	the initial partial heading of this root
 * @param pbl       		base length of the parent root, where this root emerges
 * @param pni       		local node index, where this root emerges
 * @deprecated moved		indicates if nodes were moved in the previous time step (default = false)
 * @param oldNON    		the number of nodes of the previous time step (default = 0)
 */
Stem::Stem(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
		Vector3d partialIHeading_, int pni, bool moved, int oldNON)
:Organ(id, param, alive, active, age, length, 
Matrix3d(Vector3d(0., 0., 1.), Vector3d(0., 1., 0.), Vector3d(1., 0., 0.)),  
pni, moved,  oldNON), partialIHeading(partialIHeading_)
{}

/**
 * Constructor
 * This is a Copy Paste of the Root.cpp but it works independently, it has its own parameter file (in .stparam file) tropism, growth function, txt and vtp writing system.
 * All of those can be modified to fit the real growth of the Plant.
 *
 * Typically called by the Plant::Plant(), or Stem::createNewStem().
 * For stem the initial node and node emergence time (netime) must be set from outside
 *
 * @param plant 		points to the plant
 * @param parent 		points to the parent organ
 * @param subtype		sub type of the stem
 * @param delay 		delay after which the organ starts to develop (days)
 * @param rheading		relative heading (within parent organ)
 * @param pni			parent node index
 * @param pbl			parent base length
 */
Stem::Stem(std::shared_ptr<Organism> plant, int type, Matrix3d iHeading, double delay,  std::shared_ptr<Organ> parent, int pni)
:Organ(plant, parent, Organism::ot_stem, type, delay, iHeading, pni)
{
	assert(parent!=nullptr && "Stem::Stem parent must be set");
	auto p = this->param();
	addPhytomerId(p->subType);
	double beta = getphytomerId(p->subType)*M_PI*getStemRandomParameter()->rotBeta +
			M_PI*plant->rand()*getStemRandomParameter()->betaDev;
	beta = beta + getStemRandomParameter()->initBeta*M_PI;
	if (getStemRandomParameter()->initBeta >0 && getphytomerId(p->subType)==0 ){
		beta = beta + getStemRandomParameter()->initBeta*M_PI;
	}
	double theta = p->theta;//M_PI*p->theta;
	if (parent->organType()!=Organism::ot_seed) { // scale if not a base organ, to delete?
		double scale = getStemRandomParameter()->f_sa->getValue(parent->getNode(pni), parent);
		theta *= scale;
	}
	//used when computing actual heading, @see Stem::getIHeading
	this->partialIHeading = Vector3d::rotAB(theta,beta);
	if (parent->organType()!=Organism::ot_seed) { // initial node
		//if lateral of stem, initial creation time: 
		//time when stem reached end of basal zone (==CT of parent node of first lateral) + delay
		// @see stem::createLateral
		double creationTime;
		if (parent->getNumberOfChildren() == 0){creationTime = parent->getNodeCT(pni)+delay;
		}else{creationTime = parent->getChild(0)->getParameter("creationTime") + delay;}
		
		addNode(Vector3d(0.,0.,0.), parent->getNodeId(pni), creationTime);
	//std::cout<<"stem created "<<getId()<<" "<<getNodeId(0)<<" "<<pni<<" "<<parent->getNodeId(pni)<<std::endl;
	}
}

/**
 * Deep copies the organ into the new plant @param rs.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
std::shared_ptr<Organ> Stem::copy(std::shared_ptr<Organism> p)
{
	auto s = std::make_shared<Stem>(*this); // shallow copy
	s->parent = std::weak_ptr<Organ>();
	s->plant = p;
	s->param_ = std::make_shared<StemSpecificParameter>(*param()); // copy parameters
	for (size_t i=0; i< children.size(); i++) {
		s->children[i] = children[i]->copy(p); // copy laterals
		s->children[i]->setParent(s);
	}
	return s;
}

/**
 * Simulates growth of this stem for a time span dt
 *
 * @param dt       time step [day]
 * @param verbose  indicates if status messages are written to the console (cout) (default = false)
 */
void Stem::simulate(double dt, bool verbose)
{
	if(verbose||(active &&(getParameter("subType")==2)&&plant.lock()->verboseh))
	{
		std::cout<<"          stem::simulate "<<getId()<<" "<<this->param_->subType
				<<" "<<dt<<" "<<age<<" ";
				std::cout<<"length "<<length<<" "<<getLength(int(nodes.size()-1))<<" "
				<<getLength(true)<<" "<<getLength(false)<<std::endl;
				//std::cout<<" kid/node "<< children.size()<<" "<< nodes.size()<<std::endl;
				
				
// 		std::cout<<"get kid pni";
// 		for(int u = 0; u < children.size();u++)
// 		{
// 			std::cout<<u<<") "<<children.at(u)->parentNI<<" "<<getLength(children.at(u)->parentNI)<<", ";
// 		}std::cout<<std::endl;
// 		std::cout<<"get node"<<std::endl;
// 		for(int u = 0; u < nodes.size();u++)
// 		{
// 			std::cout<<u<<") "<<getLength(u)<<", ";
// 		}std::cout<<std::endl;
				
				
	}
	if(!getOrganism()->hasRelCoord()){
		throw std::runtime_error("organism no set in rel coord");
	}
	const StemSpecificParameter& p = *param(); // rename
	
	if(children.size()>0)
	{
		double l_ = getLength(children.at(0)->parentNI);
		if( p.lb - l_ > 1e-11){std::cout<<"START first ked below lb "<<getId()<<" "<<l_<<" "<<p.lb<<std::endl;assert(false);}
	}
	
	firstCall = true;
	oldNumberOfNodes = nodes.size();
	auto p_all = plant.lock();
	auto p_stem = p_all->getOrganRandomParameter(Organism::ot_stem);



	int additional_childern;
	additional_childern = 0;
	

	if (alive) { // dead roots wont grow

		// increase age
		if (age+dt>p.rlt) { // root life time
			dt=p.rlt-age; // remaining life span
			alive = false; // this root is dead
		}
		age+=dt;

		// probabilistic branching model (todo test)
		if ((age>0) && (age-dt<=0)) { // the root emerges in this time step
			//use relative coordinates for this function. Delete as it s not a root?
			double P = getStemRandomParameter()->f_sbp->getValue(nodes.back(),shared_from_this());
			if (P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
				double p = 1.-std::pow((1.-P), dt); //probability of emergence in this time step
				if (plant.lock()->rand()>p) { // not rand()<p
					age -= dt; // the root does not emerge in this time step
				}
			}
		}

		if (age>0) { // unborn  roots have no children

			// children first (lateral roots grow even if base root is inactive)
			if(children.size()>0)
			{
				double l_ = getLength(children.at(0)->parentNI);
				if( p.lb - l_ > 1e-11){std::cout<<"before kid first ked below lb "<<l_<<" "<<p.lb<<std::endl;assert(false);}
			}
			for (auto l:children) {
				l->simulate(dt,verbose);
			}
			if(children.size()>0)
			{
				double l_ = getLength(children.at(0)->parentNI);
				if( p.lb - l_ > 1e-11){std::cout<<"after kid first ked below lb "<<getId()<<" "<<l_<<" "<<p.lb<<std::endl;assert(false);}
			}

			if(verbose){std::cout<<"stem::simulateAFTERCHILDREN "<<getId()<<" "<<this->param_->subType
				<<" "<<dt<<" "<<age<<" ";
				std::cout<<"length "<<length<<" "<<getLength(int(nodes.size()-1))<<" "
				<<getLength(true)<<" "<<getLength(false);
				std::cout<<" kid/node "<< children.size()<<" "<< nodes.size()<<std::endl;
				
				
				std::cout<<"get kid pni"<<std::endl;
				for(int u = 0; u < children.size();u++)
				{
					std::cout<<u<<") "<<children.at(u)->parentNI<<" "<<getLength(children.at(u)->parentNI)<<", ";
				}std::cout<<std::endl;
				std::cout<<"get node"<<std::endl;
				for(int u = 0; u < nodes.size();u++)
				{
					std::cout<<u<<") "<<getLength(u)<<", ";
				}std::cout<<std::endl;
						
				
				}
				
				
			if(children.size()>0)
			{
				double l_ = getLength(children.at(0)->parentNI);
				if( p.lb - l_ > 1e-11){std::cout<<"afterkid first ked below lb "<<getId()<<" "<<l_<<" "<<p.lb<<std::endl;assert(false);}
			}
			if (active) {
                
                double dt_; // time step
                if (age<dt) { // the root emerged in this time step, adjust time step
                    dt_= age;
                } else {
                    dt_=dt;
                }

				// length increment
				double age__ = age;
				if(age > p.delayNGStart){//simulation ends after start of growth pause
					if(age < p.delayNGEnd){age__ =p.delayNGStart;//during growth pause
					}else{
						age__ = age - (p.delayNGEnd - p.delayNGStart);//simulation ends after end of growth pause
					}
				}//delay to apply 
				/*as we currently do not implement impeded growth for stem and leaves
				*we can use directly the organ's age to cumpute the target length
				*/
                double e;
                
                if(verbose || ((getParameter("subType")==2)&&plant.lock()->verboseh))
                {
                    std::cout<<"stem::simulate, check growth rate "<<getOrganism()->useCWGr<<" "<<this->getOrganRandomParameter()->f_gf->CW_Gr.empty()<<" "<<age__<<std::endl;
                }
                double rmax, Lmax;
                switch(budStage) 
                {
                    case -1:{Lmax = length; break;}
                    case 0:{rmax = plant.lock()->budGR;//1 mm/d
                            Lmax = plant.lock()->maxLBudDormant; 
                            if(parentLinkingNode == 1)//2nd bud
                            {
                                Lmax = plant.lock()->maxLBudDormant_1; 
                            }
                            break;}
                    case 1 :{rmax = plant.lock()->budGR;Lmax = plant.lock()->maxLBud;break;}//1 mm/d
                    case 2 :{rmax = getParameter("r");
                             Lmax = getParameter("k");break;}//1 mm/d
                    default:{std::cout<<"stem::simulate: budStage not recognised "<< budStage<<std::flush;
                            assert(false);}
                }
                
                if(this->getOrganRandomParameter()->f_gf->CW_Gr.empty())
                {
                    
                    if(Lmax - length < -1e-10)
                    {
                        std::cout<<"Lmax - length < -1e-10 "<<getId()<<" "<<Lmax<<" "<<
                            length<<" "<<this->epsilonDx<<std::flush;
                        assert(false);
                    }  
                    
                    double age_ = getStemRandomParameter()->f_gf->getAge(length, rmax, Lmax, shared_from_this());
                    double LinitTemp = getStemRandomParameter()->f_gf->getLength(age_  , rmax, Lmax, shared_from_this());
                    double targetlength = getStemRandomParameter()->f_gf->getLength(age_ +dt_ , rmax, Lmax, shared_from_this()) + this->epsilonDx;
                    e = std::max(0.,targetlength-LinitTemp); // unimpeded elongation in time step dt
                    
                    if(verbose || ((getParameter("subType")==2)&&plant.lock()->verboseh))
                    {
                        std::cout<<"yes empty "<<LinitTemp<<" "<<targetlength<<" "<<e  <<" "<<length<<" "<<this->epsilonDx<<std::endl;
                    }
                    
                    if((e + getLength(true)) - Lmax> 1e-10){
                        std::cout<<"Stem::simulate: target length too high "<<e<<" "<<dt<<" "<<getLength(false);
                        std::cout<<" "<<Lmax<<" "<<getId()<<" "<<budStage<<std::endl;
                        assert(false);
                    }
                }else{
                    double targetlength = calcLength(age__)+ this->epsilonDx;//length from piafmunch module
                    e = targetlength-length; // store value of elongation to add
                    if(verbose || ((getParameter("subType")==2)&&plant.lock()->verboseh))
                    {
                        std::cout<<"NO has was "<<targetlength<<" "<<e  <<" "<<length<<" "<<this->epsilonDx<<std::endl;
                    }
                }
                assert(((budStage !=1)||(length - plant.lock()->maxLBud < 1e-10))&&"!(((budStage !=1)||(length - 1 < 1e-10)))");
                assert(((budStage !=0)||(parentLinkingNode == 1)||(length - plant.lock()->maxLBudDormant < 1e-10))&&"!(((budStage !=0)||(length < 1e-10)))");
                assert(((budStage !=0)||(parentLinkingNode != 1)||(length - plant.lock()->maxLBudDormant_1 < 1e-10))&&"!(((budStage !=0)||(length < 1e-10)))");
				double dl = e;//length increment = calculated length + increment from last time step too small to be added
				length = getLength(true);
				this->epsilonDx = 0.; // now it is "spent" on length (no need for -this->epsilonDx in the following)
				// create geometry
				if (p.laterals||bool(additional_childern)) { // stem has laterals
					//std::cout<<"sim seed nC is"<< nC<<"\n";
					//std::cout<<"sim seed nZ is"<< nZ<<"\n";
					/*
                    shoot born root
// 					 */
// 					if ((dl>0)&&(length< nZ)) {
// 						if (length+dl <= nZ) {
// 							createSegments(dl,verbose);
// 							length+=dl ;
// 							dl=0;
// 						} else {
// 							double ddx = nZ - length;
// 							createSegments(ddx,verbose);//should it not be ddx here?

// 							dl-=ddx;
// 							shootBorneRootGrow(verbose);
// 							length = nZ;
// 						}


// 					}

					/* basal zone */
					if ((dl>0)&&(length<p.lb)) { // length is the current length of the root
						if (length+dl<=p.lb) {
							createSegments(dl,verbose);
							length+=dl;
							dl=0;
						} else {
							double ddx = p.lb-length;
							createSegments(ddx,verbose);
							dl-=ddx; // ddx already has been created
							length=p.lb;
							
							//if(this->epsilonDx != 0){//this sould not happen as p.lb was redefined in rootparameter::realize to avoid this
							//	throw std::runtime_error("Stem::simulate: p.lb - length < dxMin");
							//}
						}
					}
					/* branching zone */
					//go into branching zone if organ has laterals and has reached 
					//the end of the basal zone
                    
						if(verbose){std::cout<<"to branching zone "<<std::endl;}
						double objectifLength = std::accumulate(p.ln.begin(),p.ln.begin()+std::max(0,created_linking_node-1),p.lb); 
						bool longEnough = ((length- objectifLength) > -1e-11);
						bool needLat = (created_linking_node < (p.ln.size() + 1));
						bool canGrow = ((dl>0)||(getStemRandomParameter()->createFirst && needLat) );
					//internodal elongation, if the basal zone of the stem is created and still has to grow
					if(verbose){
						std::cout<<"create nodes? "<<getId()<<" "<<created_linking_node <<" "<< (p.ln.size())<<" length "
						<<length<<" "<<p.lb<<" "<<(length>=p.lb)<<" non "<< nodes.size()<<" objLen ";
						std::cout<< objectifLength <<" "<< longEnough <<" needl "<< needLat<<" canGrow "<< canGrow <<std::endl;
					}
						while((length>=p.lb) && canGrow && (!longEnough || needLat))
					{
						
						if(verbose){
							std::cout<<"create lat?"<<getId()<<" "<<" "<<this->param_->subType<<" "<< longEnough <<" "<< needLat<<" "
							<<created_linking_node <<" "<<  (p.ln.size() + 1) <<" "<<
							 length<<" "<< objectifLength<<" "<< length- objectifLength<<std::endl;
						}
						//current objective length before creating next linkging nodes
						//while (needGrow && canGrow) 
						//{
							//
							// create seg
							if(longEnough && needLat){
								//std::cout<<"create segment "<<created_linking_node<<" "<<(created_linking_node>0)<<std::endl;
								
								//std::cout<<"get kid pni"<<std::endl;
								//for(int u = 0; u < children.size();u++)
								//{
								//	std::cout<<u<<") "<<children.at(u)->parentNI<<" "<<getLength(children.at(u)->parentNI)<<", ";
								//}std::cout<<std::endl;
								//std::cout<<"get node"<<std::endl;
								//for(int u = 0; u < nodes.size();u++)
								//{
								//	std::cout<<u<<") "<<getLength(u)<<", ";
								//}std::cout<<std::endl;
								if(created_linking_node>0){
									firstCall = false;
									createSegments(this->dxMin(),verbose, -1);
									dl-=this->dxMin();
									length+=this->dxMin();
								//std::cout<<"create segment ! "<<created_linking_node<<" "<<(created_linking_node>0)<<std::endl;
								
								// std::cout<<"get kid pni"<<std::endl;
								// for(int u = 0; u < children.size();u++)
								// {
								// 	std::cout<<u<<") "<<children.at(u)->parentNI<<" "<<getLength(children.at(u)->parentNI)<<", ";
								// }std::cout<<std::endl;
								// std::cout<<"get node"<<std::endl;
								// for(int u = 0; u < nodes.size();u++)
								// {
								// 	std::cout<<u<<") "<<getLength(u)<<", ";
								// }std::cout<<std::endl;
								// 
                                }
								if(verbose){
									std::cout<<"to create lateral "<<getId()<<" "<<this->param_->subType<<" "
									<<children.size()<<" "<<created_linking_node 
									<<" "<< nodes.size()<<std::endl;
								}
								
								//create lat
								createLateral(verbose);
								
								//s += p.ln.at(created_linking_node - 2);//new obj length
								objectifLength = std::accumulate(p.ln.begin(),p.ln.begin()+(created_linking_node-1),p.lb); 
								longEnough = ((length- objectifLength) > -1e-11);
								needLat = (created_linking_node < (p.ln.size() + 1));
								if(verbose){
									std::cout<<"after lat "<<getId()<<" "<<this->param_->subType<<" "
									<< longEnough <<" "<< needLat<<" "
									<<created_linking_node <<" "<<  (p.ln.size() + 1) <<" "<<
									 length<<" "<< objectifLength<<" "<< length- objectifLength
									<<" "<<this->param_->subType <<std::endl;
								}
							}
							if(verbose){
								std::cout<<"to create segment after lateral "<<children.size()<<" "<<this->param_->subType
								<<" "<<created_linking_node 
								<<" "<<p.ln.size()<<" "<< nodes.size()<<std::endl;
							}
							//}
							//growth 
							//double maxInternodeDistance = s - p.lb;
							if(verbose){
								std::cout<<"nodalgrowth? "<<getId()<<" "<<this->param_->subType<<" "<<objectifLength <<" "<<p.getK()
								<<" "<<dl<<" "<< nodes.size()<<" "<<children.size()<<" "<<
								created_linking_node<<" "<<p.ln.size()<<" "<<
								length<<" "<< length- objectifLength<<" "<<this->param_->subType<<std::endl;
								std::cout<<"get kid pni"<<std::endl;
								for(int u = 0; u < children.size();u++)
								{
									std::cout<<u<<") "<<children.at(u)->parentNI<<" "<<getLength(children.at(u)->parentNI)<<", ";
								}std::cout<<std::endl;
								std::cout<<"get node"<<std::endl;
								for(int u = 0; u < nodes.size();u++)
								{
									std::cout<<u<<") "<<getLength(u)<<", ";
								}std::cout<<std::endl;
							}
							 
							if((dl>0)&&(!longEnough))
							{
								int nn = children.at(children.size()-1)->parentNI; //node carrying the last lateral == end of branching zone
								double currentLength = getLength(nn);// - p.lb; //actual length of branching zone
								//double maxInternodeDistance = std::accumulate(p.ln.begin(),p.ln.begin()+created_linking_node -1,p.lb);
								double ddx = std::min(objectifLength-currentLength, dl);//length to add to branching zone 
								
								if(verbose){
									std::cout<<"elongation? "<<" "<<this->param_->subType<<" "
									<<objectifLength<<" "<<p.getK()<<" "<<dl<< " ";
									std::cout<<nn<<" "<< nodes.size()<<" "<<currentLength<<" "<<ddx<<std::endl;
								}
								
								if(ddx > 0){
									internodalGrowth(ddx, verbose);
									dl -= ddx;
								length += ddx;
									
								}
								
								
							}
									
							canGrow = ((dl>0)||(getStemRandomParameter()->createFirst && needLat));
							longEnough = ((length- objectifLength) > -1e-11);
							
							//check for next loop
							//needGrow = (created_linking_node < (p.ln.size()));
							
						//}
						//not needed anymore?
						//createLateral(verbose);
						// if (getStemRandomParameter()->getLateralType(getNode(nodes.size()-1))==2){
										// leafGrow(verbose);
						// }
					}
					if((length>=p.lb)&&(((p.ln.size()+1)!=created_linking_node))&&(getStemRandomParameter()->createFirst)){
						std::stringstream errMsg;
						errMsg <<"Stem::simulate(): different number of realized laterals ("<<created_linking_node<<
						") and max laterals ("<<p.ln.size()+1<<")";
						throw std::runtime_error(errMsg.str().c_str());
					}
					
					
					/* apical zone */
					//only grows once the basal and branching nodes are developped
					//double maxInternodeDistance = std::accumulate(p.ln.begin(),p.ln.end(),0.);
					if(verbose){
								
						std::cout<<"to apical growth? "<<" "<<this->param_->subType<<" "<<dl<<" "<<length
						<<" "<<objectifLength<<" "<<p.lb<<" "<<
						(dl>0)<<" "<<(length-objectifLength)<<std::endl;
					}
					if ((length>=p.lb) && canGrow && (longEnough || !needLat))
					{
					//((dl>0)&&(length-(maxInternodeDistance + p.lb) > -1e-11)) {
						createSegments(dl,verbose);
						length+=dl;
						dl = 0;
					} 
				} else { // no laterals
					if (dl>0) {
						createSegments(dl,verbose);
						length+=dl;
						dl = 0;
						
						}
				} // if lateralgetLengths
			if(dl <0){ //to keep in memory that realised length is too long, as created nodes to carry children
									
				this->epsilonDx += dl;//targetlength + e - length;
				length += this->epsilonDx;//go back to having length = theoratical length
			}
                active = getLength(false)<=(p.getK()*(1 - 1e-11)); // become inactive, if final length is nearly reached
			} // if active
			//set limit below 1e-10, as the test files see if correct length 
			//once rounded at the 10th decimal
			//@see test/test_stem_ng.py
		}
	} // if alive
	if(verbose){std::cout<<"stem id "<<getId()<<" "<<this->param_->subType<<" end simulate "<<length<<" "<< epsilonDx <<" ";
				std::cout<<p.ln.size()<<" "<< p.laterals <<" "<< children.size()<<std::endl;
				std::cout<<nodes.size()<<" "<<getLength(int(nodes.size() - 1))<<std::endl;
			std::cout<<"get kid pni"<<std::endl;
			for(int u = 0; u < children.size();u++)
			{
				std::cout<<u<<") "<<children.at(u)->parentNI<<" "<<getLength(children.at(u)->parentNI)<<", ";
			}std::cout<<std::endl;
			std::cout<<"get node"<<std::endl;
			for(int u = 0; u < nodes.size();u++)
			{
				std::cout<<u<<") "<<getLength(u)<<", ";
			}std::cout<<std::endl;
				}
	if(children.size()>0)
	{
		double l_ = getLength(children.at(0)->parentNI);
		if( p.lb - l_ > 1e-11){std::cout<<"END first ked below lb "<<getId()<<" "<<l_<<" "<<p.lb<<std::endl;assert(false);}
	}
}



/**
 * Simulates internodal growth of dl for this stem
 * divid total stem growth between the phytomeres
 * currently two option:
 * growth devided equally between the phytomeres or 
 * the phytomere grow sequentially
 *
 * @param 	dl			total length of the segments that are created [cm]
 * @param	verbose		print information
 */
void Stem::internodalGrowth(double dl, bool verbose)
{
	if(verbose){
		std::cout<<"Stem::internodalGrowth "<<std::endl;
	}
	const StemSpecificParameter& p = *param(); // rename
	std::vector<double> toGrow(created_linking_node - 1);//p.ln.size());
	double dl_;
	//wont be needed anymore no? only if stochasticity makes ssome ln at 0
	//number of laterals wich grow on smae branching point as the one before
	const int ln_0 = std::count(p.ln.cbegin(), p.ln.cbegin() + created_linking_node, 0);
	if(p.nodalGrowth==0){//sequentiall growth
		toGrow[0] = dl;
		std::fill(toGrow.begin()+1,toGrow.end(),0) ;
	}
	if(p.nodalGrowth ==1)
	{//equal growth
		std::fill(toGrow.begin(),toGrow.end(),dl/(toGrow.size()-ln_0)) ; 
	}
	
	int numPhyto = 0; int numPhyto_ =0;
	int firstLat = 0; int secLat;
	while( (dl > 1e-12)&&(firstLat < (children.size()-1))) {//(numPhyto_ <(p.ln.size()*2)) 
		//if the phytomere can do a growth superior to the mean phytomere growth, we add the value of "missing" 
		//(i.e., length left to grow to get the predefined total growth of the branching zone)
		int nn1 = children.at(firstLat)->parentNI; //node at the beginning of phytomere
		secLat = firstLat + 1;
		for (secLat = firstLat + 1;(secLat < children.size())&&(children.at(secLat)->parentNI == nn1 );secLat ++){
			if(verbose){
				std::cout<<"add seclat "<<nn1<<" "<<children.at(secLat)->parentNI<<" "<<firstLat<<" "<<secLat<<" "<<children.size()<<std::endl;
			}
			//secLat++;
			}//look for id of next linking node
		int nn2 = children.at(secLat)->parentNI; //node at end of phytomere
		if(nn2 != nn1)
		{
			double availableForGrowth = p.ln.at(numPhyto) -( getLength(nn2) - getLength(nn1) ) ;//difference between maximum and current length of the phytomere
			dl_ = std::min(std::min(toGrow[numPhyto],availableForGrowth), dl);
			if((numPhyto+1)< toGrow.size()){//p.ln.size()){
				toGrow[numPhyto + 1] +=  toGrow[numPhyto] - dl_ ;
			}
			if(verbose){
				std::cout<<"in nodal growth elongation_a "<< numPhyto <<" "<<nn1<<" "<<nn2<<" "<<firstLat
				<<" "<<secLat<<" "<<p.ln.at(numPhyto) <<" "<< getLength(nn2) <<" "<< getLength(nn1) <<" "
				<<availableForGrowth<<std::endl;
			}
			if(dl_ > 0){
				if(verbose){
									std::cout<<"in nodal growth elongation "
									<<dl_ <<" "<<(numPhyto+1)<<std::endl;
								}
				createSegments(dl_,verbose, nn2); dl -= dl_;}
				if(verbose){
					std::cout<<"back to internodalGrowth "<<numPhyto<<" "<< numPhyto_<<" "<<firstLat <<" "<<secLat;
					std::cout<<" "<<p.ln.size()<<" "<<children.size()<<" "<<	dl<<std::endl;
				}
			//back to internodalGrowth 1 1 1 2 2 2.22045e-16
			//do the loop at most twice
			numPhyto ++;numPhyto_ ++;
			if(numPhyto>=toGrow.size()){numPhyto=0;}
		}
			firstLat = secLat;
	}
	if(std::abs(dl)> 1e-6){//this sould not happen as computed dl to be <= sum(availableForGrowth)
		std::stringstream errMsg;
		errMsg <<"Stem::internodalGrowth length left to grow: "<<dl;
		throw std::runtime_error(errMsg.str().c_str());
	}
}
/**
 * Returns a parameter per organ
 *
 * @param name 		parameter name (returns nan if not available)
 *
 */
double Stem::getParameter(std::string name) const
{
	if (name=="lb") { return param()->lb; } // basal zone [cm]
	if (name=="delayNGStart") { return param()->delayNGStart; } // delay for nodal growth [day]
	if (name=="delayNGEnd") { return param()->delayNGEnd; } // delay for nodal growth [day]
	if (name=="la") { return param()->la; } // apical zone [cm]
	if (name=="nob") { return param()->nob(); } // number of branches
	if (name=="r"){ return param()->r; }  // initial growth rate [cm day-1]
	if (name=="radius") { return param()->a; } // root radius [cm]
	if (name=="a") { return param()->a; } // root radius [cm]
	if (name=="theta") { return param()->theta; } // angle between root and parent root [rad]
	if (name=="rlt") { return param()->rlt; } // root life time [day]
	if (name=="k") { return param()->getK(); }; // maximal root length [cm]
	if (name=="lnMean") { // mean lateral distance [cm]
		auto& v =param()->ln;
		return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
	}
	if (name=="lnDev") { // standard deviation of lateral distance [cm]
		auto& v =param()->ln;
		double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
		double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
		return std::sqrt(sq_sum / v.size() - mean * mean);
	}
	if (name=="volume") { return param()->a*param()->a*M_PI*getLength(true); } // // root volume [cm^3]
	if (name=="surface") { return 2*param()->a*M_PI*getLength(true); }
	if (name=="type") { return this->param_->subType; }  // in CPlantBox the subType is often called just type
	if (name=="parentNI") { return parentNI; } // local parent node index where the lateral emerges
	return Organ::getParameter(name);
}


/**
 * Analytical creation (=emergence) time of a node at a length along the stem
 *
 * @param length   length of the stem [cm]
 */
double Stem::calcCreationTime(double length)
{
	assert(length >= 0 && "Stem::getCreationTime() negative length");
	double stemage = calcAge(length);
	stemage = std::min(stemage, age);
	assert(stemage >= 0 && "Stem::getCreationTime() negative stem age");
	return stemage+nodeCTs[0];
}

/**
 * Analytical length of the stem at a given age
 *
 * @param age          age of the stem [day]
 */
double Stem::calcLength(double age)
{
	assert(age>=0 && "Stem::calcLength() negative root age");
	return getStemRandomParameter()->f_gf->getLength(age,getStemRandomParameter()->r,param()->getK(),shared_from_this());
}

/**
 * Analytical age of the stem at a given length
 * no scaling of organ growth , so can return age directly
 * otherwise cannot compute exact age between delayNGStart and delayNGEnd
 * @param length   length of the stem [cm]
 */
double Stem::calcAge(double length)
{
	assert(length>=0 && "Stem::calcAge() negative root age");
	double age__ = getStemRandomParameter()->f_gf->getAge(length,getStemRandomParameter()->r,param()->getK(),shared_from_this());
	if(age__ >param()->delayNGStart ){age__ += (param()->delayNGEnd - param()->delayNGStart);}
	return age__;
}

/**
 * Creates a new lateral by calling Stem::createNewstem().
 *
 * Overwrite this method to implement more specialized stem classes.
 */
void Stem::createLateral(bool verbose)
{ 
	auto rp = getStemRandomParameter(); // rename
	auto sp = param(); // rename
	Matrix3d h = Matrix3d(); //not needed anymore
	int ot;
	if(verbose){
		std::cout<<"Stem::createLateral "<<getId()<<" "<<created_linking_node<<" "<<children.size()<<std::endl;
		}
	for(int i = 0; i < rp->successorST.size(); i++){//go through each successor rule
		//found id
		bool foundId, applyAll;
		if(rp->successorWhere.size()>i){
			foundId = (std::find (rp->successorWhere.at(i).begin(), rp->successorWhere.at(i).end(), created_linking_node)
				!= rp->successorWhere.at(i).end());
			//applies to all linking nodes
			applyAll = (std::find (rp->successorWhere.at(i).begin(), rp->successorWhere.at(i).end(), -1)
				!=rp->successorWhere.at(i).end());
		}else{applyAll = true;} //default
		
		if(verbose){
			std::cout<<"at i "<<i<<" found "<<foundId<<" "<<applyAll<<std::endl;
			}
		if(foundId || applyAll  )//aply here
		{
			int numlats = 1;//how many laterals? default = 1
			if(rp->successorNo.size()>i){numlats =  rp->successorNo.at(i);}
			if(verbose){std::cout<<"numlats "<<numlats <<std::endl;}
			for(int nn = 0; nn < numlats; nn++)
			{
				// if(
				// (std::adjacent_find( myvector.begin(), myvector.end(), std::not_equal_to<>() ) 
				// == myvector.end())&&(successorP[i][0]==1.))
				// {p=nn;}else{
					const Vector3d& pos = Vector3d();
					int p_id = rp->getLateralType(pos, i);//if probabilistic branching
					//}//, soccessorOT[i], successorST[i])
				if(rp->successorOT.size()>i){ot = rp->successorOT.at(i).at(p_id);
				}else{ot = getParameter("organType");}//default
				int st = rp->successorST.at(i).at(p_id);
				double ageLN = this->calcAge(length); // sp->lb age of stem when lateral node is created
				double delay = sp->delayLat * created_linking_node;
				if(verbose){std::cout<<"from stem data 4 lat "<<getId()<<" "<<ot<<" "<<st<<" "<<p_id<<" "<< (nodes.size() - 1)<<" ageLN "<<ageLN<<" "<<age;
				
						std::cout<<" "<<getNodeId(nodes.size() - 1)<<" "<<getNodeId(0)<<std::endl;}
                        int created_linking_node_ = created_linking_node;
                        if(getParameter("subType") != 1){created_linking_node_ = this->parentLinkingNode;}
				switch(ot){
					case Organism::ot_root:{
						double  beta = i*M_PI*getStemRandomParameter()->rotBeta;
						Vector3d newHeading = iHeading.times(Vector3d::rotAB(0,beta));
						auto lateral = std::make_shared<Root>(plant.lock(), st, newHeading, delay, shared_from_this(),  nodes.size() - 1);
						children.push_back(lateral);
                        lateral->parentLinkingNode = created_linking_node_ ; 
						lateral->simulate(age-ageLN,verbose); 
						break;}
					case Organism::ot_stem:{
						auto lateral = std::make_shared<Stem>(plant.lock(), st, h, delay, shared_from_this(),  nodes.size() - 1);
						children.push_back(lateral);
                        lateral->parentLinkingNode = created_linking_node_ ; 
                        if(st == 2)//auxillary bud
                        {lateral->budStage = 0; //dormant bud
                        }
						lateral->simulate(age-ageLN,verbose); 
						break;}
					case Organism::ot_leaf:{
						auto lateral = std::make_shared<Leaf>(plant.lock(), st, h, delay, shared_from_this(),  nodes.size() - 1);
						children.push_back(lateral);
                        lateral->parentLinkingNode = created_linking_node_; 
						lateral->simulate(age-ageLN,verbose); 
						break;}
				}
				
			}
		}
		
	}
	created_linking_node ++;
	if(verbose){
		std::cout<<"Stem::createLateral_END "<<getId()<<" "<<created_linking_node<<" "<<children.size()<<std::endl;
		}
	
}


/*
 * Searches for subtype of first leaf
 * wont be needed anymore
 */
int Stem::getLeafSubType()
{
	auto orp = plant.lock()->getOrganRandomParameter(Organism::ot_leaf);
	for(int st_ = 1; st_ < orp.size();st_++)//skipe st_ ==0, never an organ st
	{
		if(orp[st_] != NULL) {
			return orp[st_]->subType;
		}
	}
	return -1;
}

/*
 *
 */
void Stem::shootBorneRootGrow(bool verbose)
{
	//const double maxT = 365.; // maximal simulation time
	auto sp = this->param(); // rename
	auto stem_p = this->param();
	auto p = plant.lock();
	auto p_seed = p->getOrganRandomParameter(Organism::ot_seed,0);
	int st = p->getParameterSubType(Organism::ot_root, "shootborne");
	if (st>0) {
		shootborneType = st;
	} // otherwise stick with default
	try {
		p->getOrganRandomParameter(Organism::ot_root, shootborneType); // if the type is not defined an exception is thrown
	} catch (...) {
		if (verbose) {
			std::cout << "Seed::initialize:Shootborne root type #" << shootborneType << " was not defined, using tap root parameters instead\n";
		}
		shootborneType =1;
	}
	//        std::cout <<"subtype ="<<stem_p->subType <<"stem getPhytomerId =" <<getphytomerId(stem_p->subType)<< "\n";
	addPhytomerId(stem_p->subType);

	int nC = getPlant()->getSeed()->param()->nC;
	double nZ = getPlant()->getSeed()->param()->nz;
	double res = nZ-floor(nZ / dx())*dx();
	if(res < dxMin() && res != 0){
		if(res <= dxMin()/2){ nZ -= res;
		}else{nZ =  floor(nZ / dx())*dx() + dxMin();}
		if(verbose){std::cout<<"\nStem::shootBorneRootGrow: nZ changed to "<<nZ<<" for compatibility with dx and dxMin"<<std::endl;}
	}
	if ( nC>0 ) { // only if there are any shootborne roots && (p_seed->firstSB+p_seed->delaySB<maxT)
		std::cout<<"seed nC is "<< nC <<"\n";
		std::cout<<"seed nZ is "<< nZ <<"\n";
		double ageLN = this->calcAge(getLength(true)); // age of stem when lateral node is created
		double ageLG = this->calcAge(getLength(true)+sp->la); // age of the stem, when the lateral starts growing (i.e when the apical zone is developed)
		double delay = ageLG-ageLN; // time the lateral has to wait
		for (int i=0; i< nC; i++) {
			double  beta = i*M_PI*getStemRandomParameter()->rotBeta;
			Vector3d newHeading = iHeading.times(Vector3d::rotAB(0,beta));
			auto shootBorneRoot = std::make_shared<Root>(plant.lock() , shootborneType, newHeading, delay,
					shared_from_this(), nodes.size() - 1);
			children.push_back(shootBorneRoot);
			shootBorneRoot->simulate(age-ageLN,verbose);
			std::cout<<"root grow number "<<i<<"\n";
		}
	}
}

/**
 * Returns the increment of the next segments
 *
 *  @param p       coordinates of previous node (in absolute coordinates)
 *  @param sdx     length of next segment [cm]
 *  @param n       index of the node at the beginning of the segment
 *  @return        the vector representing the increment
 */
Vector3d Stem::getIncrement(const Vector3d& p, double sdx, int n)
{
	Vector3d h = heading(n);
	Matrix3d ons = Matrix3d::ons(h);
	Vector2d ab = getStemRandomParameter()->f_tf->getHeading(p, ons, dx(), shared_from_this(), n+1);
	Vector3d sv = ons.times(Vector3d::rotAB(ab.x,ab.y));
	return sv.times(sdx);
}


/**
 * @return Current absolute heading of the organ at node n, based on initial heading, or direction of the segment going from node n-1 to node n
 */
Vector3d Stem::heading(int n ) const
{
	if(n<0){n=nodes.size()-1 ;}
	if ((nodes.size()>1)&&(n>0)) {
		n = std::min(int(nodes.size()),n);
		Vector3d h = getNode(n).minus(getNode(n-1));
		h.normalize();
		return h;
	} else {
		return getiHeading0();
	}
}

/**
 * @return Current absolute heading of the organ at node n, based on initial heading, or segment before
 */
Vector3d Stem::getiHeading0()  const
{	
	Matrix3d iHeading;
	if (getParent()->organType()==Organism::ot_seed) { // from seed?
		iHeading = Matrix3d(Vector3d(0, 0, 1), Vector3d(0, 1, 0), Vector3d(1, 0, 0));
	}else{
		Vector3d vIHeading = getParent()->heading(parentNI);
		iHeading = Matrix3d::ons(vIHeading);};
	auto heading = iHeading.column(0);
	Vector3d new_heading = Matrix3d::ons(heading).times(this->partialIHeading);
	return Matrix3d::ons(new_heading).column(0);
}
/**
 *  Creates nodes and node emergence times for a length l
 *
 *  Checks that each new segments length is <= dx but >= smallDx
 *
 *  @param l        total length of the segments that are created [cm]
 *  @param verbose  turns console output on or off
 *  @param nodeIdx   index of the node to shift. default = -1 (last node)
 */
void Stem::createSegments(double l, bool verbose, int nodeIdx)
{
	if (l==0) {
		std::cout << "Stem::createSegments: zero length encountered \n";
		return;
	}
	if (l<0) {
		std::cout << "Stem::createSegments: negative length encountered \n";
	}

	// shift first node to axial resolution
	double shiftl = 0; // length produced by shift
	int nn = nodes.size();
	if(verbose){
		std::cout<<"Stem::createSegments "<<l<<" "<<nodeIdx<<std::endl;
	}
  if( nodeIdx >= 0){ //if we are doing internodal growth,  nodeIdx >= 0.
		//autonodeIdx o = children.at(nodeIdx);
		nn = nodeIdx +1; //shift the last node of the phytomere n° nodeIdx instead of the last node of the organ
		if(verbose){
			std::cout<<"nodeIdx >= 0 "<<nn<<std::endl;
		}
	}else
	{
		if(verbose){
			std::cout<<"nodeIdx < 0 "<<nn<< " "<<firstCall<<" "<<(firstCall||(nodeIdx >= 0))<<std::endl;
		}
	}		
	if (firstCall||(nodeIdx >= 0)) { // first call of createSegments (in Root::simulate)
		firstCall = false;
		if(verbose){
			std::cout<<"in make move nn:"<<nn<<" "<<(nn>1)<<" ";
		}
		if ((nn>1)) { // don't move first node
			Vector3d h = nodes[nn-1];
			double olddx = h.length(); // length of last segment
			if(verbose){
				std::cout<<"len "<<h.length()<<" dx "<<dx()<<" "<<olddx<<" "<<dx()*0.99;
			}
			if (dx()-olddx > 1e-13) { // shift node instead of creating a new node
				shiftl = std::min(dx()-olddx, l);
				double sdx = olddx + shiftl; // length of new segment
				h.normalize();  
				nodes[nn-1] =  h.times(sdx);
				double et = this->calcCreationTime(getLength(true)+shiftl);
				nodeCTs[nn-1] = et; // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
				l -= shiftl;
				if (l<=0) { // ==0 should be enough
					return;
				}
				if(verbose){
					std::cout<<" shift "<<shiftl<<" "<<sdx<<std::endl;
				}
			} 
		} 
	}
	// create n+1 new nodes
	double sl = 0; // summed length of created segment
	int n = floor(l/dx());
	if(verbose){std::cout<<"seg to create "<<l<<" "<<n<<" "<<dx()<<" "<<this->epsilonDx<<"\n";}
	for (int i = 0; i < n + 1; i++) {

		double sdx; // segment length (<=dx)
		if (i<n) {  // normal case
			sdx = dx();
		} else { // last segment
			sdx = l-n*dx();
			if (sdx<dxMin()*0.99 ) { // quit if l is too small
				if (verbose&& sdx != 0) {
					std::cout << "length increment below dx threshold ("<< sdx <<" < "<< dxMin() << ") and kept in memory\n";
				}
				if( nodeIdx >= 0){
					this->epsilonDx += sdx;
				}else{this->epsilonDx += sdx;}
				return;
			}
			//this->epsilonDx = 0; //no residual
		}
		sl += sdx;
		Vector3d newnode = Vector3d(sdx, 0., 0.);
		double et = this->calcCreationTime(getLength(true)+shiftl+sl);
		// in case of impeded growth the node emergence time is not exact anymore,
		// but might break down to temporal resolution
		bool shift = (nodeIdx >= 0); //node will be insterted between 2 nodes. only happens if we have internodal growth (PhytoIdx >= 0)
		if(verbose){
			std::cout<<"stem addNode "<<getId()<<" "<<newnode.toString()<<" "<<size_t(nn+i)<<" "<<shift<<" "<<nodeIdx
			<<" "<<this->epsilonDx<<std::endl;
		}
		addNode(newnode, et, size_t(nn+i), shift);
	}
}
/**
 * @return the organs length from start node up to the node with index @param i.
 */
double Stem::getLength(int i) const 
{
	double l = 0.; // length until node i
	if(getOrganism()->hasRelCoord()){//is currently using relative coordinates?
		for (int j = 0; j<i; j++) {
			l += nodes.at(j+1).length(); // relative length equals absolute length
		}
	}else{
		for (int j = 0; j<i; j++) {
			l += nodes.at(j+1).minus(nodes.at(j)).length(); // relative length equals absolute length
		}
	}
	return l;
}

 /* @param realized	FALSE:	get theoretical organ length, INdependent from spatial resolution (dx() and dxMin()) 
 *					TRUE:	get realized organ length, dependent from spatial resolution (dx() and dxMin())
 *					DEFAULT = TRUE
 * @return 			The chosen type of organ length (realized or theoretical).
 */
double Stem::getLength(bool realized) const
{
	if (realized) {
		return length - this->epsilonDx;
	} else {
		return length;
	}
}

/**
 * convert the nodes' positions from relative to absolute coordinates
 */
void Stem::rel2abs() 
{
	
	nodes[0] = getOrigin(); //recompute postiion of the first node
	
	for(size_t i=1; i<nodes.size(); i++){
		Vector3d newdx = nodes[i];
		if((i>=oldNumberOfNodes)|| active){//if we have a new node or have nodal growth, need to update the tropism effect
			double sdx = nodes[i].length();
			newdx = getIncrement(nodes[i-1], sdx, i-1); //add tropism
		}
		nodes[i] = nodes[i-1].plus(newdx); //replace relative by absolute position
		
		
	}
	//if carry children, update their pos
	
	for(size_t i=0; i<children.size(); i++){
		(children[i])->rel2abs();
	}
	
}

/**
 *  convert the nodes' positions from absolute to relative coordinates
 */
void Stem::abs2rel()
{
	for (int j = nodes.size(); j>1; j--) {
		nodes[j-1] = nodes.at(j-1).minus(nodes.at(j-2));
	}
	nodes[0] = Vector3d(0.,0.,0.);
	for(size_t i=0; i<children.size(); i++){
		(children[i])->abs2rel();
	}//if carry children, update their pos
	
}


/**
 * @return The RootTypeParameter from the plant
 */
std::shared_ptr<StemRandomParameter> Stem::getStemRandomParameter() const
{
	return std::static_pointer_cast<StemRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_stem, param_->subType));
}

/**
 * @return Parameters of the specific root
 */
std::shared_ptr<const StemSpecificParameter> Stem::param() const
{
	return std::static_pointer_cast<const StemSpecificParameter>(param_);
}

/*
 * Quick info about the object for debugging
 * additionally, use getParam()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Stem::toString() const
{
	std::stringstream newstring;
	newstring << "; initial heading: " << iHeading.toString() << ", parent node index" << parentNI << ".";
	return Organ::toString()+newstring.str();
}

} // namespace CPlantBox
