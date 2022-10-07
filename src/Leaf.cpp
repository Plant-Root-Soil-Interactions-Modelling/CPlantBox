#include "Stem.h"
#include "Leaf.h"
#include "Root.h"
#include "Plant.h"

namespace CPlantBox {

/**
 * Constructs a leaf from given data.
 * The organ tree must be created, @see Organ::setPlant, Organ::setParent, Organ::addChild
 * Organ geometry must be created, @see Organ::addNode, ensure that this->getNodeId(0) == parent->getNodeId(pni)
 *
 * @param id        	the organ's unique id (@see Organ::getId)
 * @param param     	the organs parameters set, ownership transfers to the organ
 * @param alive     	indicates if the organ is alive (@see Organ::isAlive)
 * @param active    	indicates if the organ is active (@see Organ::isActive)
 * @param age       	the current age of the organ (@see Organ::getAge)
 * @param length    	the current length of the organ (@see Organ::getLength)
 * @param iheading  	the initial heading of this leaf
 * @param pbl       	base length of the parent leaf, where this leaf emerges
 * @param pni       	local node index, where this leaf emerges
 * @deprecated moved	as long as stem is active, nodes are assumed to have moved (@see Organism::getUpdatedNodes)
 * @param oldNON    	the number of nodes of the previous time step (default = 0)
 */
Leaf::Leaf(int id, const std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length,
		Matrix3d iHeading,int pni, bool moved, int oldNON)
		:Organ(id, param, alive, active, age, length, iHeading, pni, moved,  oldNON )
{}

/**
 * Constructor
 * Typically called by the Plant::Plant(), or Leaf::createNewLeaf().
 * For leaf the initial node (in nodes) and node emergence time (in nodeCTs) must be set from outside
 *
 * @param plant 		points to the plant
 * @param parent 		points to the parent organ
 * @param subtype		sub type of the leaf
 * @param delay 		delay after which the organ starts to develop (days)
 * @param rheading		relative heading (within parent organ)
 * @param pni			parent node index
 * @param pbl			parent base length
 */
Leaf::Leaf(std::shared_ptr<Organism> plant, int type, Matrix3d iHeading, double delay,  std::shared_ptr<Organ> parent, int pni)
:Organ(plant, parent, Organism::ot_leaf, type, delay, iHeading, pni)
{
	assert(parent!=nullptr && "Leaf::Leaf parent must be set");
	addleafphytomerID(param()->subType);
	ageDependentTropism = getLeafRandomParameter()->f_tf->ageSwitch > 0;
	// Calculate the rotation of the leaves. The code begins here needs to be rewritten, because another following project will work on the leaves. The code here is just temporally used to get some nice visualizations. When someone rewrites the code, please take "gimbal lock" into consideration.  
	//Rewritten Begin: 															 
	double beta = getleafphytomerID(param()->subType)*M_PI*getLeafRandomParameter()->rotBeta
			+ M_PI*plant->rand()*getLeafRandomParameter()->betaDev ;  //+ ; //2 * M_PI*plant->rand(); // initial rotation
	beta = beta + getLeafRandomParameter()->initBeta*M_PI;
	if (getLeafRandomParameter()->initBeta >0 && getLeafRandomParameter()->subType==2 && getLeafRandomParameter()->lnf==5 && getleafphytomerID(2)%4==2) {
		beta = beta + getLeafRandomParameter()->initBeta*M_PI;
	} else if (getLeafRandomParameter()->initBeta >0 && getLeafRandomParameter()->subType==2 && getLeafRandomParameter()->lnf==5 && getleafphytomerID(2)%4==3) {
		beta = beta + getLeafRandomParameter()->initBeta*M_PI + M_PI;
	}
	double theta = param()->theta;
	if (parent->organType()!=Organism::ot_seed) { // scale if not a base leaf
		double scale = getLeafRandomParameter()->f_sa->getValue(parent->getNode(pni), parent);
		theta *= scale;
	}
	//used when computing actual heading, @see LEaf::getIHeading
	this->partialIHeading = Vector3d::rotAB(theta,beta);
	// Rewritten ends 
	if (parent->organType()!=Organism::ot_seed) { // if not base organ
	
		double creationTime;
		if (parent->organType()==Organism::ot_stem) {
			//if lateral of stem, initial creation time: 
			//time when stem reached end of basal zone (==CT of parent node of first lateral) + delay
			// @see stem::leafGrow
			if (parent->getNumberOfChildren() == 0){creationTime = parent->getNodeCT(pni)+delay;
			}else{creationTime = parent->getChild(0)->getParameter("creationTime") + delay;}
		}else{
			creationTime = parent->getNodeCT(pni)+delay;
		}
		addNode(Vector3d(0.,0.,0.), parent->getNodeId(pni), creationTime);//create first node. relative coordinate = (0,0,0)
	}
	//std::cout<<"leaf created "<<getId()<<" "<<getNodeId(0)<<" "<<pni<<" "<<parent->getNodeId(pni)<<std::endl;

}

/**
 * Deep copies the organ into the new plant @param rs.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
std::shared_ptr<Organ> Leaf::copy(std::shared_ptr<Organism> p)
{
	auto l = std::make_shared<Leaf>(*this); // shallow copy
	l->parent = std::weak_ptr<Organ>();
	l->plant = p;
	l->param_ = std::make_shared<LeafSpecificParameter>(*param()); // copy parameters
	for (size_t i=0; i< children.size(); i++) {
		l->children[i] = children[i]->copy(p); // copy laterals
		l->children[i]->setParent(l);
	}
	return l;
}

/**
 * Simulates f_gf of this leaf for a time span dt
 *
 * @param dt       time step [day]
 * @param verbose  indicates if status messages are written to the console (cout) (default = false)
 */
void Leaf::simulate(double dt, bool verbose)
{
	
	if(verbose){std::cout<<"leaf::simulate "<<getId()<<" "<<this->param_->subType
				<<" "<<dt<<" "<<age<<" ";
				std::cout<<length<<" "<< children.size()<<" "<< nodes.size()<<std::endl;}
	firstCall = true;
	oldNumberOfNodes = nodes.size();

	const LeafSpecificParameter& p = *param(); // rename

	if (alive) { // dead leafs wont grow

		// increase age
		if (age+dt>p.rlt) { // leaf life time
			dt=p.rlt-age; // remaining life span
			alive = false; // this leaf is dead
		}
		age+=dt;

		// probabilistic branching model (todo test)
		if ((age>0) && (age-dt<=0)) { // the leaf emerges in this time step
			//currently, does not use absolute coordinates for these function. 
			double P = getLeafRandomParameter()->f_sbp->getValue(nodes.back(),shared_from_this());
			if (P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
				double p = 1.-std::pow((1.-P), dt); //probability of emergence in this time step
				if (plant.lock()->rand()>p) { // not rand()<p
					age -= dt; // the leaf does not emerge in this time step
				}
			}
		}

		if (age>0) { // unborn leafs have no children

			// children first (lateral leafs grow even if base leaf is inactive)
			for (auto l:children) {
				l->simulate(dt,verbose);
			}

			if (active) {

				// length increment
				double age_ = calcAge(length); // leaf age as if grown unimpeded (lower than real age)
				double dt_; // time step
				if (age<dt) { // the leaf emerged in this time step, adjust time step
					dt_= age;
				} else {
					dt_=dt;
				}

				double targetlength = calcLength(age_+dt_)+ this->epsilonDx;
				double e = targetlength-length; // unimpeded elongation in time step dt
				double dl = std::max(e, 0.);// length increment = calculated length + increment from last time step too small to be added
				length = getLength(true);
				this->epsilonDx = 0.; // now it is "spent" on targetlength (no need for -this->epsilonDx in the following)
				// create geometry
				if (p.laterals) { // leaf has laterals
					/* basal zone */
					if ((dl>0)&&(length<p.lb)) { // length is the current length of the leaf
						if (length+dl<=p.lb) {
							createSegments(dl,verbose);
							length+=dl;
							dl=0;
						} else {
							double ddx = p.lb-length;
							createSegments(ddx,verbose);
							dl-=ddx; // ddx already has been created
							length=p.lb;
						}
					}
					double s = p.lb; // summed length
					/* branching zone */
						if(verbose){
							std::cout<<"to branching zone? "<<p.ln.size()<<" "<<created_linking_node
							<<" "<<getLength(true)<<" "<<dl<<" "<<nodes.size()<<std::endl;
						}
					if ((dl>0)&&(length>=p.lb)) {
						for (size_t i=0; ((i<p.ln.size()) && (dl>0)); i++) {
							if(verbose){
							std::cout<<"elongate before lat "<<i<<" "<<p.ln.size()<<" "<<created_linking_node
							<<" "<<s<<" "<<p.ln.at(i)<<" "<<length<<" "<<nodes.size()<<std::endl;
						}
							s+=p.ln.at(i);
							if (length<s) {
								if(verbose){std::cout<<"if (length<s) {"<<std::endl;}
								if (i== created_linking_node ) { // new lateral
								if(verbose){std::cout<<"createLateral"<<std::endl;}
									createLateral(verbose);
								}
								if (length+dl<=s) { // finish within inter-lateral distance i
								if(verbose){std::cout<<"if (length+dl<=s) {"<<std::endl;}
									createSegments(dl,verbose);
									length+=dl;//- this->epsilonDx;
									dl=0;
								} else { // grow over inter-lateral distance i
								if(verbose){std::cout<<"if (length+dl>s) {"<<std::endl;}
									double ddx = s-length;
									createSegments(ddx,verbose);
									dl-=ddx;
									length=s;
								}
							}
						}
						if(verbose){
							std::cout<<"to last lat? "<<p.ln.size()<<" "<<created_linking_node
							<<" "<<getLength(true)<<" "<<s<<" "<<nodes.size()<<std::endl;
						}
						if ((p.ln.size()==created_linking_node)&& (getLength(true)>=s)) { // new lateral (the last one)
							if(verbose){
								std::cout<<"do last lat! ";
							}
							createLateral(verbose);
						}
					}
					/* apical zone */
					if (dl>0) {
						createSegments(dl,verbose);//y not with dt_?
						length+=dl;//- this->epsilonDx;
					}
				} else { // no laterals
					if (dl>0) {
						createSegments(dl,verbose);
						length+=dl;//- this->epsilonDx;
					}
				} // if lateralgetLengths
			} // if active
			//level of precision = 1e-10 to not create an error in the test files
			active = getLength(false)<=(p.getK()*(1 - 1e-11)); // become inactive, if final length is nearly reached
		}
	} // if alive
	if(verbose){std::cout<<"leaf id "<<getId()<<" end simulate "<<length<<" "<< epsilonDx <<" ";
				std::cout<<p.ln.size()<<" "<< p.laterals <<" "<< children.size()<<std::endl;
				std::cout<<nodes.size()<<" "<<getLength(int(nodes.size() - 1))<<std::endl;
				}

}

/**
 *
 */
double Leaf::getParameter(std::string name) const {
	if (name=="shapeType") { return getLeafRandomParameter()->shapeType; } // definition type of the leaf shape 
	if (name=="Width_petiole") { return param()->Width_petiole; } // [cm]
	if (name=="Width_blade") { return param()->Width_blade; } // [cm]
	if (name=="lb") { return param()->lb; } // basal zone [cm]
	if (name=="la") { return param()->la; } // apical zone [cm]
	//if (name=="nob") { return param()->nob; } // number of branches
	if (name=="r"){ return param()->r; }  // initial growth rate [cm day-1]
	if (name=="radius") { return param()->a; } // leaf radius or thickness [cm]
	if (name=="a") { return param()->a; } // leaf radius or thickness [cm]
	if (name=="theta") { return param()->theta; } // angle between leaf and parent root [rad]
	if (name=="rlt") { return param()->rlt; } // leaf life time [day]
	if (name=="k") { return param()->getK(); }; // maximal leaf length [cm]
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
	if (name=="volume_th") { return orgVolume(-1, false); } // // theoretical leaf volume [cm^3]
	if (name=="surface_th") { return leafArea(false); } // // theoretical leaf surface [cm^2]
	if (name=="volume_realized") { return orgVolume(-1, true); } // // realized leaf volume [cm^3]
	if (name=="surface_realized") { return leafArea(true); } // // realized leaf surface [cm^2]
	if (name=="volume") { return orgVolume(-1, true); } // // realized leaf volume [cm^3]
	if (name=="surface") { return leafArea(true); } // // realized leaf surface [cm^2]
	if (name=="type") { return this->param_->subType; }  // in CPlantBox the subType is often called just type
	if (name=="parentNI") { return parentNI; } // local parent node index where the lateral emerges
	return Organ::getParameter(name);
}


/**
 * in case there are no lateral leafs return leaf surface area [cm2]
 * upper side only. If used for photosynthesis, 
 * with C3 plants (stomata on upper + lower side) need to do * 2
 * @param realized		use realized (true) or theoretical (false) length and area (default = false)
 * @param withPetiole	take into account leaf petiole or sheath (true) or not (false). Default = false (for computation of transpiration)
 * @return 	total leaf blade Area  (withPetiole == false) or total leaf Area (withPetiole == true) [cm2] 
 */
double Leaf::leafArea(bool realized, bool withPetiole) const
{																			 
	double length_ = getLength(realized);
	double surface_ = 0;
	double surfacePetiole = 0;
	if (param()->laterals) {//no blade
		return 0.;
	} else {
		int shapeType = getLeafRandomParameter()->shapeType;
		switch(shapeType) 
		{
			case LeafRandomParameter::shape_cuboid:{ 
				double Width_blade = getParameter("Width_blade") ;
				double Width_petiole = getParameter("Width_petiole") ;
				if (length_ <= param()->lb) {
					surfacePetiole =  Width_petiole * length_ ; 
				} else {
					//surface of basal zone
					surfacePetiole = Width_petiole *param()->lb  ;
					//surface rest of leaf
											  
					length_ -= param()->lb;
											
					double surfaceBlade =  Width_blade * length_ ;
					surface_ =  surfaceBlade;				
				}
				if(withPetiole){surface_ += surfacePetiole;}
				return surface_;
			
			} break;
			case LeafRandomParameter::shape_cylinder:{
				// divide by two to get only upper side of leaf
				double perimeter =  2 * M_PI * param()->a;
				if (length_ <= param()->lb) {
					surfacePetiole =  perimeter  * length_ /2; 
				} else {
					//surface of basal zone
					surfacePetiole = perimeter  *param()->lb /2 ;
					//surface rest of leaf
											  
					length_ -= param()->lb;
											
					double surfaceBlade =  perimeter  * length_ /2;
					surface_ =  surfaceBlade;				
				}
				if(withPetiole){surface_ += surfacePetiole;}
				return surface_;
				
			} break;
			case LeafRandomParameter::shape_2D:{
				// how to take into account possible petiole area? add perimeter  *param()->lb /2 ?
				return param()->areaMax * (leafLength(realized)/param()->leafLength());
			} break;
			
			default:
				throw  std::runtime_error("Leaf::leafArea: undefined leaf shape type");
		}
	}
	return 0.;
};

/**
 * leaf BLADE Area at segment n°localSegId
 * upper side only. If used for photosynthesis, 
 * with C3 plants (stomata on upper + lower side) need to do * 2
 * see @XylemFlux::segFluxes and @XylemFlux::linearSystem 
 * @param localSegId	index for which evaluate area == nodey_localid + 1
 * @param realized		use realized (true) or theoretical (false) length and area (default = false)
 * @param withPetiole	take into account leaf petiole or sheath (true) or not (false). Default = false (for computation of transpiration)
 * @return 	leaf area at segment n°localSegId [cm2]
 */
double Leaf::leafAreaAtSeg(int localSegId, bool realized, bool withPetiole)
{
	double surface_ = 0.;
	if (param()->laterals) {//no blade
		return 0.;
	} else {
		int shapeType = getLeafRandomParameter()->shapeType;
		auto n1 = nodes.at(localSegId);
		auto n2 = nodes.at(localSegId + 1);
		auto v = n2.minus(n1);
		double length_ = v.length();
		double lengthAt_x = getLength(localSegId);
		double lengthInPetiole = std::min(length_,std::max(param()->lb - lengthAt_x,0.));//petiole or sheath
		double lengthInBlade = std::max(length_ - lengthInPetiole, 0.);
		assert(((lengthInBlade+lengthInPetiole)==length_)&&"leafAreaAtSeg: lengthInBlade+lengthInPetiole !=lengthSegment");
		switch(shapeType) 
		{
			case LeafRandomParameter::shape_cuboid:{ 
				double Width_blade = getParameter("Width_blade") ;
				
				double surfaceBlade =  Width_blade * lengthInBlade ;
				surface_ =  surfaceBlade ;
				double surfacePetiole = 0;
				if(withPetiole)
				{ 
					double Width_petiole = getParameter("Width_petiole") ;
					surfacePetiole =   Width_petiole * lengthInPetiole ;
					surface_ +=  surfacePetiole;
				}
			} break;
			case LeafRandomParameter::shape_cylinder:{
				// divide by two to get only upper side of leaf
				surface_ =  2 * M_PI * lengthInBlade * param()->a / 2; 
				if(withPetiole)
				{ 
					surface_ +=  2 * M_PI * lengthInPetiole * param()->a / 2; 
				}				
				
			} break;
			case LeafRandomParameter::shape_2D:{
				//TODO: compute it better later? not sur how to do it if the leaf is not convex
				// how to take into account possible petiole area? add perimeter  *lengthInPetiole /2 ?
				surface_ = (lengthInBlade / leafLength(realized)) * leafArea(realized);
			} break;
			
			default:
				throw  std::runtime_error("Leaf::leafAreaAtSeg: undefined leaf shape type");
		}
	}
	if(surface_ < 1e-15){ surface_ = 0;}
	return surface_;
};

/**
 * leaf BLADE Area at segment n°localSegId
 * upper side only. If used for photosynthesis, 
 * with C3 plants (stomata on upper + lower side) need to do * 2
 * see @XylemFlux::segFluxes and @XylemFlux::linearSystem 
 * @param localSegId	index for which evaluate area == nodey_localid + 1
 * @param realized		use realized (true) or theoretical (false) length and area (default = false)
 * @param withPetiole	take into account leaf petiole or sheath (true) or not (false). Default = false (for computation of transpiration)
 * @return 	leaf area at segment n°localSegId [cm2]
 */
double Leaf::leafLengthAtSeg(int localSegId, bool withPetiole)
{
	
	
	double length_out = 0.;
	if (!(param()->laterals)) {
		auto n1 = nodes.at(localSegId);
		auto n2 = nodes.at(localSegId + 1);
		auto v = n2.minus(n1);
		double length_ = v.length();
		double lengthAt_x = getLength(localSegId);
		double lengthInPetiole = std::min(length_,std::max(param()->lb - lengthAt_x,0.));//petiole or sheath
		double lengthInBlade = std::max(length_ - lengthInPetiole, 0.);
		assert(((lengthInBlade+lengthInPetiole)==length_)&&"leafAreaAtSeg: lengthInBlade+lengthInPetiole !=lengthSegment");
		length_out = lengthInBlade;
		if(withPetiole){length_out += lengthInPetiole;}
	}
	return length_out;
};



/**
 * leaf BLADE Area at segment n°localSegId
 * upper side only. If used for photosynthesis, 
 * with C3 plants (stomata on upper + lower side) need to do * 2
 * see @XylemFlux::segFluxes and @XylemFlux::linearSystem 
 * @param localSegId	index for which evaluate area == nodey_localid + 1
 * @param realized		use realized (true) or theoretical (false) length and area (default = false)
 * @param withPetiole	take into account leaf petiole or sheath (true) or not (false). Default = false (for computation of transpiration)
 * @return 	leaf area at segment n°localSegId [cm2]
 */
double Leaf::leafVolAtSeg(int localSegId,bool realized, bool withPetiole)
{
	double vol_ = 0.;
	//if (param()->laterals) {
	//	return 0.;
	//} else {
		int shapeType = getLeafRandomParameter()->shapeType;
		auto n1 = nodes.at(localSegId);
		auto n2 = nodes.at(localSegId + 1);
		auto v = n2.minus(n1);
		double length_ = v.length();
		double lengthAt_x = getLength(localSegId);
		double lengthInPetiole = std::min(length_,std::max(param()->lb - lengthAt_x,0.));//petiole or sheath
		double lengthInBlade = std::max(length_ - lengthInPetiole, 0.);
		double a = getParameter("a") ;//radius or thickness
		assert(((lengthInBlade+lengthInPetiole)==length_)&&"leafVolAtSeg: lengthInBlade+lengthInPetiole !=lengthSegment");
		switch(shapeType) 
		{
			case LeafRandomParameter::shape_cuboid:{ 
				double Width_blade = getParameter("Width_blade") ;
				
				double volBlade =  Width_blade * lengthInBlade *a;
				vol_ =  volBlade ;
				double volPetiole = 0;
				if(withPetiole)
				{ 
					double Width_petiole = getParameter("Width_petiole") ;
					volPetiole =   Width_petiole * lengthInPetiole *a;
					vol_ +=  volPetiole;
				}
			} break;
			case LeafRandomParameter::shape_cylinder:{
				// divide by two to get only upper side of leaf
				vol_ =  M_PI * lengthInBlade * param()->a * param()->a; 
				if(withPetiole)
				{ 
					vol_ +=  M_PI * lengthInPetiole * param()->a * param()->a; 
				}				
				
			} break;
			case LeafRandomParameter::shape_2D:{
				//TODO: compute it better later? not sur how to do it if the leaf is not convex
				// how to take into account possible petiole area? add perimeter  *lengthInPetiole /2 ?
				vol_ = (lengthInBlade / leafLength(realized)) * leafArea(realized) *a;
			} break;
			
			default:
				throw  std::runtime_error("Leaf::leafVolAtSeg: undefined leaf shape type");
		}
	//}
	if(!std::isfinite(vol_ )){
		std::stringstream errMsg;
		errMsg <<"Leaf::leafVolAtSeg:END computation of leaf volume failed "<<localSegId<<" "<<getParameter("Width_petiole")
			<<" "<<getLeafRandomParameter()->shapeType<<" "<<getParameter("Width_blade") <<" "<<param()->laterals<<"\n";
		throw std::runtime_error(errMsg.str().c_str());
	}
	return vol_;
};



/**
 * @param length_	total leaf length for which to evaluate volume. default = -1 (i.e., use current volume)
 *					for phloem module, need to compute volume for other lengths
 * @param realized		use realized (true) or theoretical (false) length and area (default = false)
 * @return leaf volume [cm3]
 */
double Leaf::orgVolume(double length_, bool realized) const
{
	double vol_;
	const LeafSpecificParameter& p = *param(); 
	int shapeType = getLeafRandomParameter()->shapeType;																						 
	if(length_ == -1){length_ = getLength(realized);}//theoretical
	switch(shapeType) 
	{
		case LeafRandomParameter::shape_cuboid:{ 
			double Width_blade = getParameter("Width_blade") ;
			double Width_petiole = getParameter("Width_petiole") ;
			if ((p.laterals)||(length_ <= p.lb)) {
				vol_ =  Width_petiole * length_ * p.a;
			} else {
				//volume of basal zone
				double volPetiole = Width_petiole * p.lb * p.a ;//assume p.a is thickness
				//volume rest of leaf
				length_ -= p.lb;
				double volBlade =  Width_blade * length_ * p.a; //assume p.a is thickness
				vol_ =  volBlade + volPetiole;
			}
		}break;
		case LeafRandomParameter::shape_cylinder:{
			vol_ = length_ * p.a * p.a * M_PI;
		} break;
		case LeafRandomParameter::shape_2D:{
			// how to take into account possible petiole volume? add lengthPetiole * p.a * p.a * M_PI;  ?
			vol_ = leafArea(realized) * p.a ;//assume p.a is thickness
		} break;
		default:
			throw  std::runtime_error("Leaf::orgVolume: undefined leaf shape type");
	}
	return vol_;
};

/**
 * @param volume_	total leaf length for which to evaluate volume.
 *					for phloem module, need to compute lengths for different volumes
 * @return leaf length [cm]
 */
double Leaf::orgVolume2Length(double volume_) 
{
	const LeafSpecificParameter& p = *param(); 
	double length_;
	int shapeType = getLeafRandomParameter()->shapeType;	
	switch(shapeType) 
		{
			case LeafRandomParameter::shape_cuboid:{ 
				double Width_blade = getParameter("Width_blade") ;
				double Width_petiole = getParameter("Width_petiole") ;
				double volPetiole = Width_petiole * p.lb * p.a;//assume p.a is thickness
				if(volume_ <= volPetiole){
					length_ = volume_/( Width_petiole * p.a);//assume p.a is thickness
				}else{
					double lengthBlade = (volume_ - volPetiole)/p.a/Width_blade;
					length_ = p.lb + lengthBlade;
				}
			} break;
			case LeafRandomParameter::shape_cylinder:
			{
				length_ = volume_/(p.a*p.a*M_PI);
			} break;
			case LeafRandomParameter::shape_2D:{
				double area_ = volume_ / p.a;//assume p.a is thickness
				bool realized = false;
				length_ = leafLength(realized ) * (area_ / p.areaMax); //assume area/areaMax = length / lengthmax
			} break;
			default:
				throw  std::runtime_error("Leaf::orgVolume2Length: undefined leaf shape type");
	}
	return length_;
};

/**
 * indicates if the node is in the leaf surface are and should be viusalized as polygon
 *
 * leaf base (false), branched leaf (false), or leaf surface area (true)
 */
bool Leaf::nodeLeafVis(double l)
{
	if (param()->laterals) {
		return false;
	} else {
		return l >= param()->lb; // true if not in basal zone
	}
}

/**
 * Parameterization x value, at position l along the leaf axis
 */
std::vector<double> Leaf::getLeafVisX_(double l) {
	auto& lg = getLeafRandomParameter()->leafGeometry;
	int n = lg.size();
	int ind = int( ((l - param()->lb) /leafLength())*(n-1) + 0.5); // index within precomputed normalized geometry
	auto x_ = lg.at(ind); // could be more than one point for non-convex geometries
	return x_;
}

/**
 * for Python binding
 */
std::vector<double> Leaf::getLeafVisX(int i) {
	return getLeafVisX_(getLength(i));
}

/**
 * Scales unit leaf shape to the specific leaf,
 * and returns leaf shape coordinates per node (normally 2 points, for convex domain it could be more points)
 * see used by vtk_plot.py to create a polygon representation of the leaf area
 */
std::vector<Vector3d> Leaf::getLeafVis(int i)
{
	double l = getLength(i);
	if (nodeLeafVis(l)) {
		auto& lg = getLeafRandomParameter()->leafGeometry;
		int n = lg.size();
		if (n>0) {
			std::vector<Vector3d> coords;
			auto x_ = getLeafVisX_(l);
			Vector3d x1= getiHeading0();
			x1.normalize();
			Vector3d y1 = Vector3d(0,0,-1).cross(x1); // todo angle between leaf - halfs
			y1.normalize(); 
			
			double a  = leafArea() / leafLength(); // scale radius
			for (double x :x_) {
				coords.push_back(getNode(i).plus(y1.times(x*a)));
			}
			for (double x :x_) {
				coords.push_back(getNode(i).minus(y1.times(x*a)));
			}
			return coords;
		} else {
			std::cout << "Leaf::getLeafVis: WARNING leaf geometry was not set \n";
			return std::vector<Vector3d>();
		}
	} else { // no need for polygonal visualisation
		return std::vector<Vector3d>();
	}
}

/**
 * Analytical creation (=emergence) time of a node at a length along the leaf
 *
 * @param length   length of the leaf [cm]
 */
double Leaf::calcCreationTime(double length)
{
	assert(length >= 0 && "Leaf::getCreationTime() negative length");
	double leafage = calcAge(length);
	leafage = std::min(leafage, age);
	assert(leafage >= 0 && "Leaf::getCreationTime() negative leaf age");
	return leafage+nodeCTs[0];
}

/**
 * Analytical length of the leaf at a given age
 *
 * @param age          age of the leaf [day]
 */
double Leaf::calcLength(double age)
{
	assert(age>=0  && "Leaf::calcLength() negative root age");
	return getLeafRandomParameter()->f_gf->getLength(age,getLeafRandomParameter()->r,param()->getK(),shared_from_this());
}

/**
 * Analytical age of the leaf at a given length
 *
 * @param length   length of the leaf [cm]
 */
double Leaf::calcAge(double length)
{
	assert(length>=0 && "Leaf::calcAge() negative root length");
	return getLeafRandomParameter()->f_gf->getAge(length,getLeafRandomParameter()->r,param()->getK(),shared_from_this());
}

/**
 *
 */
void Leaf::minusPhytomerId(int subtype)
{
	getPlant()->leafphytomerID[subtype]--;
}

/**
 *
 */
int Leaf::getleafphytomerID(int subtype)
{
	return getPlant()->leafphytomerID[subtype];
}

/**
 *
 */
void Leaf::addleafphytomerID(int subtype)
{
	getPlant()->leafphytomerID.at(subtype)++;
}

/**
 * Creates a new lateral by calling Leaf::createNewleaf().
 *
 * Overwrite this method to implement more specialized leaf classes.
 * 
 * This method was done in a rush, because there will be a following project specifically
 * focus on the leaves. This should be rewritten instead of kept. 
 * 
 */
void Leaf::createLateral(bool silence)
{
	bool verbose = silence;
	auto rp = getLeafRandomParameter(); // rename
	auto sp = param(); // rename
	Matrix3d h = Matrix3d(); //not needed anymore
	int ot;
	if((this->param_->subType == 1)&&verbose){
		std::cout<<"Leaf::createLateral "<<getId()<<" "<<created_linking_node<<" "<<children.size()
		<<" "<<nodes.size()<<std::endl;
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
		if(foundId || applyAll  )//aply here
		{
			int numlats = 1;//how many laterals? default = 1
			if(rp->successorNo.size()>i){numlats =  rp->successorNo.at(i);}
			if((this->param_->subType == 1)&&verbose){std::cout<<"numlats "<<numlats <<std::endl;}
			for(int nn = 0; nn < numlats; nn++)
			{
				const Vector3d& pos = Vector3d();
				int p_id = rp->getLateralType(pos, i);
				if(rp->successorOT.size()>i){ot = rp->successorOT.at(i).at(p_id);
				}else{ot = getParameter("organType");}//default
				int st = rp->successorST.at(i).at(p_id);
				double ageLN = this->calcAge(sp->lb); // age of stem when lateral node is created
				double delay = sp->delayLat * created_linking_node;
				if(verbose){std::cout<<"from leaf data 4 lat "<<getId()<<" "<<ot<<" "<<st<<" "<<p_id<<" "<< (nodes.size() - 1) ;
				
						std::cout<<" "<<getNodeId(nodes.size() - 1)<<" "<<getNodeId(0)<<std::endl;
						}
				switch(ot){
					case Organism::ot_root:{
						double  beta = i*M_PI*getLeafRandomParameter()->rotBeta;
						Vector3d newHeading = iHeading.times(Vector3d::rotAB(0,beta));
						auto lateral = std::make_shared<Root>(plant.lock(), st, newHeading, delay, shared_from_this(),  nodes.size() - 1);
						children.push_back(lateral);
						lateral->simulate(age-ageLN,verbose); 
						break;}
					case Organism::ot_stem:{
						auto lateral = std::make_shared<Stem>(plant.lock(), st, h, delay, shared_from_this(),  nodes.size() - 1);
						children.push_back(lateral);
						lateral->simulate(age-ageLN,verbose); 
						break;}
					case Organism::ot_leaf:{
						auto lateral = std::make_shared<Leaf>(plant.lock(), st, h, delay, shared_from_this(),  nodes.size() - 1);
						children.push_back(lateral);
						lateral->simulate(age-ageLN,verbose); 
						break;}
				}
				
			}
		}
		
	}
	created_linking_node ++;
}

/**
 * computes absolute coordinates from relative coordinates
 * when this function is called, the parent organ has already
 * its absolute coordinates
 * called by @see Plant::rel2abs
 */
void Leaf::rel2abs() 
{
	double ageSwitch = getLeafRandomParameter()->f_tf->ageSwitch; //rename
	bool tropismChange = (age > ageSwitch);
	nodes[0] = getOrigin();//get absolute coordinates of first node via coordinates of parent
	for(size_t i=1; i<nodes.size(); i++){
		Vector3d newdx = nodes[i];
		//if new node or has an age-dependent tropism + reached age at which tropism changes. Might need to update the conditions if do new tropism functions
		//i.e., gradual change according to age
		if((i>= oldNumberOfNodes )|| (this->ageDependentTropism&& tropismChange)){
			double sdx = nodes[i].length();
			newdx = getIncrement(nodes[i-1], sdx, i-1);
		}
		nodes[i] = nodes[i-1].plus(newdx);
		
	}
	if(this->ageDependentTropism && tropismChange){
		this->ageDependentTropism = false; //switch done
	}
	for(size_t i=0; i<children.size(); i++){
		(children[i])->rel2abs();
	}//if carries children, update their coordinates from relative to absolute
	
	
}

/**
 * computes relative coordinates from absolute coordinates
 * when this function is called, the parent organ has already
 * its relative coordinates
 * called by @see Plant::abs2rel
 */
void Leaf::abs2rel()
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
 * Returns the increment of the next segments
 *
 *  @param p       coordinates of previous node (in absolute coordinates)
 *  @param sdx     length of next segment [cm]
 *  @param n       index at which new node is to be inserted
 *  @return        the vector representing the increment
 */
Vector3d Leaf::getIncrement(const Vector3d& p, double sdx, int n)
{
	Vector3d h = heading(n);
	Matrix3d ons = Matrix3d::ons(h);
	//use dx() rather rhan sdx to compute heading
	//to make tropism independante from growth rate
	Vector2d ab = getLeafRandomParameter()->f_tf->getHeading(p, ons, dx(),shared_from_this());
	Vector3d sv = ons.times(Vector3d::rotAB(ab.x,ab.y));
	return sv.times(sdx);
}


/**
 * @return Current absolute heading of the organ at node n, based on initial heading, or segment before
 */
Vector3d Leaf::heading(int n) const
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
Vector3d Leaf::getiHeading0()  const
{
	Vector3d vIHeading = getParent()->heading(parentNI);
	Matrix3d iHeading = Matrix3d::ons(vIHeading);
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
 */
void Leaf::createSegments(double l, bool verbose)
{
	if (l==0) {
		std::cout << "Leaf::createSegments: zero length encountered \n";
		return;
	}
	if (l<0) {
		std::cout << "Leaf::createSegments: negative length encountered \n";
	}

	// shift first node to axial resolution
	double shiftl = 0; // length produced by shift
	int nn = nodes.size();
	if (firstCall) { // first call of createSegments (in Leaf::simulate)
		firstCall = false;
		if (nn>1){ // don't move first node of organ
			Vector3d h = nodes[nn-1];
			double olddx = h.length(); // length of last segment
			if (dx()- olddx > 1e-13) { // shift node instead of creating a new node
				shiftl = std::min(dx()-olddx, l);
				double sdx = olddx + shiftl; // length of new segment
				h.normalize();  
				nodes[nn-1] = h.times(sdx);
				double et = this->calcCreationTime(getLength(true)+shiftl);
				nodeCTs[nn-1] = et; // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
				l -= shiftl;
				if (l<=0) { // ==0 should be enough
					return;
				}
			}
		}
	}
	// create n+1 new nodes
	double sl = 0; // summed length of created segment
	int n = floor(l/dx());
	for (int i = 0; i < n + 1; i++) {
		double sdx; // segment length (<=dx)
		if (i<n) {  // normal case
			sdx = dx();
		} else { // last segment
			sdx = l-n*dx();
			if (sdx<dxMin()*0.99) { // quit if l is too small
				if (verbose&& sdx != 0) {
					std::cout <<  "Leaf::createSegments(): length increment below dxMin threshold ("<< sdx <<" < "<< dxMin() << ") and kept in memory\n";
				}
				this->epsilonDx = sdx;
				return;
			}
			this->epsilonDx = 0; //no residual
		}
		sl += sdx;
		Vector3d newnode = Vector3d(sdx, 0., 0.);//set relative position (@see rel2abs for addition of tropism)
		double et = this->calcCreationTime(getLength(true)+shiftl+sl);
		// in case of impeded growth the node emergence time is not exact anymore,
		// but might break down to temporal resolution
		addNode(newnode, et);
	}
}


/**
 * @return The LeafTypeParameter from the plant
 */
std::shared_ptr<LeafRandomParameter> Leaf::getLeafRandomParameter() const
{
	return std::static_pointer_cast<LeafRandomParameter>(plant.lock()->getOrganRandomParameter(Organism::ot_leaf, param_->subType));
}

/**
 * @return Parameters of the specific leaf
 */
std::shared_ptr<const LeafSpecificParameter>  Leaf::param() const
{
	return std::static_pointer_cast<const LeafSpecificParameter>(param_);
}

/**
 * Quick info about the object for debugging
 * additionally, use getParam()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Leaf::toString() const
{
	std::stringstream newstring;
	newstring << "; initial heading: " << getiHeading0().toString()  << ", parent node index " << parentNI << ".";
	return  Organ::toString()+newstring.str();
}

/**
 * @return the organs length from start node up to the node with index @param i.
 */
double Leaf::getLength(int i) const 
{
	double l = 0.; // length until node i
	if(getOrganism()->hasRelCoord()){
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

 /* 
 * @param realized	FALSE:	get theoretical organ length, INdependent from spatial resolution (dx() and dxMin()) 
 *					TRUE:	get realized organ length, dependent from spatial resolution (dx() and dxMin())
 *					DEFAULT = TRUE
 * @return 			The chosen type of organ length (realized or theoretical).
 */
double Leaf::getLength(bool realized) const
{
	if (realized) {
		return length - this->epsilonDx;
	} else {
		return length;
	}
}
} // namespace CPlantBox
