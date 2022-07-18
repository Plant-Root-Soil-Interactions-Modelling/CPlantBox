// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Organ.h"

#include "Organism.h"
#include "Plant.h"
#include <iostream>

#include "organparameter.h"

namespace CPlantBox {

/**
 * Constructs an organ from given data.
 * The organ tree must be created, @see Organ::setPlant, Organ::setParent, Organ::addChild
 * Organ geometry must be created, @see Organ::addNode, ensure that this->getNodeId(0) == parent->getNodeId(pni)
 *
 * @param id        the organ's unique id (@see Organ::getId)
 * @param param     the organs parameters set, ownership transfers to the organ
 * @param alive     indicates if the organ is alive (@see Organ::isAlive)
 * @param active    indicates if the organ is active (@see Organ::isActive)
 * @param age       the current age of the organ (@see Organ::getAge)
 * @param length    the current length of the organ (@see Organ::getLength)
 * @param iHeading TODO
 * @param pni
 * @param moved     indicates if nodes were moved in the previous time step (default = false)
 * @param oldNON    the number of nodes of the previous time step (default = 0)
 */
Organ::Organ(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active,
		double age, double length, Matrix3d iHeading, int pni, bool moved, int oldNON)
:iHeading(iHeading), parentNI(pni), plant(), parent(), id(id), param_(param), alive(alive), active(active), age(age),
 length(length), moved(moved), oldNumberOfNodes(oldNON)
{ }

/**
 * The constructor is used for simulation.
 * The organ parameters are chosen from random distributions within the the OrganTypeParameter class.
 * The next organ id is retrieved from the plant,
 * and the organ starts growing after a delay (starts with age = -delay).
 *
 * @param plant     the plant the new organ will be part of
 * @param parent    the parent organ, equals nullptr if there is no parent
 * @param ot        organ type
 * @param st        sub type of the organ type, e.g. different root types
 * @param delay     time delay in days when the organ will start to grow
 * @param iHeading TODO
 * @param pni
 */
Organ::Organ(std::shared_ptr<Organism> plant, std::shared_ptr<Organ> parent, int ot, int st, double delay,
		Matrix3d iHeading, int pni)
:iHeading(iHeading), parentNI(pni), plant(plant), parent(parent), id(plant->getOrganIndex()),
  param_(plant->getOrganRandomParameter(ot, st)->realize()), /* root parameters are diced in the getOrganRandomParameter class */
  age(-delay)
{ }

/*
 * Deep copies this organ into the new plant @param plant.
 * All children are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 * @return          the newly created copy (ownership is passed)
 */
std::shared_ptr<Organ> Organ::copy(std::shared_ptr<Organism>  p)
{
	auto o = std::make_shared<Organ>(*this); // shallow copy
	o->parent = std::weak_ptr<Organ>();
	o->plant = p;
	o->param_ = std::make_shared<OrganSpecificParameter>(*param_); // copy parameters
	for (size_t i=0; i< children.size(); i++) {
		o->children[i] = children[i]->copy(p); // copy lateral
		o->children[i]->setParent(o);
	}
	return o;
}

/**
 * @param realized	FALSE:	get theoretical organ length, INdependent from spatial resolution (dx() and dxMin())
 *					TRUE:	get realized organ length, dependent from spatial resolution (dx() and dxMin())
 *					DEFAULT = TRUE
 * @return 			The chosen type of organ length (realized or theoretical).
 */
double Organ::getLength(bool realized) const
{
	if (realized) {
		return length - this->epsilonDx;
	} else {
		return length;
	}
}

/**
 * @return the organs length from start node up to the node with index @param i.
 */
double Organ::getLength(int i) const
{
	double l = 0.; // length until node i
	for (int j = 0; j<i; j++) {
		l += nodes.at(j+1).minus(nodes.at(j)).length(); // relative length equals absolute length
	}
	return l;
}

/**
 * @return The organ type, which is a coarse classification of the organs,
 * for a string representation see see Organism::organTypeNames
 *
 * Currently there are: ot_organ (for unspecified organs) = 0, ot_seed = 1, ot_root = 2, ot_stem = 3, and ot_leaf = 4.
 * There can be different classes with the same organ type.
 */
int Organ::organType() const
{
	return Organism::ot_organ;
}

/**
 * @return The organ type parameter is retrieved from the plant organism.
 * The Organism class manages all organs type parameters.
 */
std::shared_ptr<OrganRandomParameter> Organ::getOrganRandomParameter() const
{
	return plant.lock()->getOrganRandomParameter(this->organType(), param_->subType);
}

/**
 * Simulates the development of the organ in a time span of @param dt days.
 *
 * @param dt        time step [day]
 * @param verbose   turns console output on or off
 */
void Organ::simulate(double dt, bool verbose)
{
	// store information of this time step
	oldNumberOfNodes = nodes.size();
	moved = false;

	// if the organ is alive, manage children
	if (alive) {
		age += dt;
		for (auto& c : children)  {
			c->simulate(dt, verbose);
		}
	}
}

/**
 *
 */
std::shared_ptr<Plant> Organ::getPlant() const
{
	return std::dynamic_pointer_cast<Plant>(plant.lock());
}

/*
 * Adds a subsequent organ (e.g. a lateral root)
 *
 * @param c     the organ to add (ownership is passed)
 */
void Organ::addChild(std::shared_ptr<Organ> c)
{
	c->setParent(shared_from_this());
	children.push_back(c);
}

/**
 * Adds a node to the organ.
 *
 * For simplicity nodes can not be deleted, organs can only become deactivated or die
 *
 * @param n        new node
 * @param id       global node index
 * @param t        exact creation time of the node
 * @param index	   position were new node is to be added
 * @param shift	   do we need to shift the nodes? (i.e., is the new node inserted between existing nodes because of internodal growth?)
 */
void Organ::addNode(Vector3d n, int id, double t, size_t index, bool shift)
{
	if(!shift){//node added at the end of organ										
		nodes.push_back(n); // node
		nodeIds.push_back(id); //unique id
		nodeCTs.push_back(t); // exact creation time
	}
	else{//could be quite slow  to insert, but we won t have that many (node-)tillers (?) 
		nodes.insert(nodes.begin() + index, n);//add the node at index
		//add a global index. 
		//no need for the nodes to keep the same global index and makes the update of the nodes position for MappedPlant object more simple)
		nodeIds.push_back(id);  
		nodeCTs.insert(nodeCTs.begin() + index-1, t);
		for(auto kid : children){//if carries children after the added node, update their "parent node index"
			if(kid->parentNI >= index-1){
				kid->moveOrigin(kid->parentNI + 1);
				}
			
		}
		
	}
}

/**
 * change idx of node linking to parent organ (in case of internodal growth)
 * @see Organ::addNode
 * @param idx      new idx
 */
void Organ::moveOrigin(int idx)
{
	this->parentNI = idx;
	nodeIds.at(0) = getParent()->getNodeId(idx);
	
}

/**
 * Adds the node with the next global index to the root
 *
 * For simplicity nodes can not be deleted, organs can only become deactivated or die
 *
 * @param n        the new node
 * @param t        exact creation time of the node
 * @param index	   position were new node is to be added
 * @param shift	   do we need to shift the nodes? (i.e., is the new node inserted between existing nodes because of internodal growth?)
 */
void Organ::addNode(Vector3d n, double t, size_t index, bool shift)
{
	addNode(n,plant.lock()->getNodeIndex(),t, index, shift);
}

/**
 * By default the organ is represented by a polyline,
 * i.e. the segments of the nodes {n1, n2, n3, n4}, are { [i1,i2], [i2,i3], [i3,i4] }, where i1-i4 are node indices.
 *
 * @return A vector of line segments, where each line segment is described as two global node indices.
 * If there are less than two nodes an empty vector is returned.
 */
std::vector<Vector2i> Organ::getSegments() const
{
	if (this->nodes.size()>1) {
		std::vector<Vector2i> segs = std::vector<Vector2i>(nodes.size()-1);
		for (size_t i=0; i<nodes.size()-1; i++) {
			Vector2i s(getNodeId(i),getNodeId(i+1));
			segs[i] = s;
		}
		return segs;
	} else {
		return std::vector<Vector2i>(0);
	}
}

/**
 * returns the maximal axial resolution
 */
double Organ::dx() const
{
	return getOrganRandomParameter()->dx;
}

/**
 * returns the minimal axial resolution,
 * length overhead is stored in epsilon and realized in the next simulation step (see Organ::getEpsilon)
 */
double Organ::dxMin() const
{
	return getOrganRandomParameter()->dxMin;
}

/**
 * Returns the organs as sequential list, copies only organs with more than one node.
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default).
 * @param all       get also the organs with only one node? default: false. Sometimes true for carbon-limited growth
 *
 * @return A sequential list of organs. If there is less than one node,
 * or another organ type is expected, an empty vector is returned.
 */
std::vector<std::shared_ptr<Organ>> Organ::getOrgans(int ot, bool all)
{
	auto v = std::vector<std::shared_ptr<Organ>> ();
	this->getOrgans(ot, v, all);
	return v;
}

/**
 * Returns the organs as sequential list, if all == false,copies only organs with more than one node.
 * if all == true return all born organs (age > 0 ) except for the seed and bulb (i.e., organs which do not grow)
 * @param ot        the expected organ type, where -1 denotes all organ types (default).
 * @param all       get also the organs with only one node? default: false. Sometimes true for carbon-limited growth
 * @param v         vector of organs where the subtree is added,
 *                  only expected organ types with more than one nodes are added.
 */
void Organ::getOrgans(int ot, std::vector<std::shared_ptr<Organ>>& v, bool all)
{
	bool notBulb = !((this->organType() == Organism::ot_stem)&&(this->getParameter("subType") == 2));//do not count leaf bulb
	//might have age <0 and node.size()> 1 when adding organ manuelly @see test_organ.py
	bool forCarbon_limitedGrowth = (all && (this->getAge()>0));//when ask for "all" organs which have age > 0 even if nodes.size() == 1
	bool notSeed = ( this->organType() != Organism::ot_seed);
	
	if ((this->nodes.size()>1 || forCarbon_limitedGrowth)&& notBulb &&notSeed) {
		if ((ot<0) || (ot==this->organType())) {
			v.push_back(shared_from_this());
		}
	}
	// std::cout << "Organ::getOrgans recursive: number of children " <<  this->children.size() << "\n" << std::flush;
	for (const auto& c : this->children) {
		c->getOrgans(ot,v, all);
	}
}

/**
 * @return The number of emerged laterals (i.e. number of children with age>0)
 * @see Organ::getNumberOfChildren
 * needed for the test files
 */
int Organ::getNumberOfLaterals() const {
	int nol = 0;
	for (auto& c : children)  {
		if (c->getAge()>0) { // born
			nol ++;
		}
	}
	return nol;
}

/**
 * Returns a single scalar parameter called @param name of the organ.
 * This method is for post processing, since it is flexible but slow.
 * Overwrite to add more parameters for specific organs.
 *
 * For OrganRandomParameters add '_mean', or '_dev',
 * to avoid naming conflicts with the organ specific parameters.
 *
 * @return The parameter value, if unknown NaN
 */
double Organ::getParameter(std::string name) const {
	// specific parameters
	if (name=="subType") { return this->param_->subType; }
    if (name=="a") { return param_->a; } // root radius [cm]
	if (name=="radius") { return this->param_->a; } // root radius [cm]
	if (name=="diameter") { return 2.*this->param_->a; } // root diameter [cm]
	// organ member variables
    if (name=="iHeadingX") { return getiHeading0().x; } // root initial heading x - coordinate [cm]
    if (name=="iHeadingY") { return getiHeading0().y; } // root initial heading y - coordinate [cm]
    if (name=="iHeadingZ") { return getiHeading0().z; } // root initial heading z - coordinate [cm]
    if (name=="parentNI") { return parentNI; } // local parent node index where the lateral emerges
    if (name=="parent-node") { // local parent node index for RSML (higher order roots are missing the first node)
        if (this->parent.expired()) {
            return -1;
        }
        if (this->parent.lock()->organType()==Organism::ot_seed) { // if it is base root
    		return -1;
    	}
        auto p = this->parent.lock();
        if (p->parent.expired()) { // if parent is base root
            return parentNI;
        }
        if (p->parent.lock()->organType()==Organism::ot_seed){ // if parent is base root
            return parentNI;
        } else {
            return std::max(parentNI-1,0); // higher order roots are missing the first node
            // TODO for 0 this can be negative... (belongs to other branch in rsml)
        }
    }
    // organ member functions
	if (name=="organType") { return this->organType(); }
    if (name=="numberOfChildren") { return children.size(); }
	if (name=="id") { return getId(); }
	if (name=="alive") { return isAlive(); }
	if (name=="active") { return isActive(); }
	if (name=="age") { return getAge(); }
    if (name=="length") { return getLength(true); } //realized organ length, dependent on dxMin and dx
	if (name=="lengthTh") { return getLength(false); } //theoratical organ length, dependent on dxMin and dx
    if (name=="numberOfNodes") { return getNumberOfNodes(); }
    if (name=="numberOfSegments") { return getNumberOfSegments(); }
    if (name=="hasMoved") { return hasMoved(); }
    if (name=="oldNumberOfNodes") { return getOldNumberOfNodes(); }
	if (name=="numberOfLaterals") { return getNumberOfLaterals(); }
    // further
    if (name=="creationTime") { return getNodeCT(0); }
	if (name=="order") { // count how often it is possible to move up
		int o = 0;
		auto p = shared_from_this();
		while ((!p->parent.expired()) && (p->parent.lock()->organType()!=Organism::ot_seed)) {
			o++;
			p = p->parent.lock(); // up the organ tree
		}
		return o;
	}
	if (name=="one") { return 1; } // e.g. for counting the organs
    return this->getOrganRandomParameter()->getParameter(name); // ask the random parameter
}

/**
 * Writes the organs RSML root tag, if it has more than one node.
 *
 * Called by Organism::getRSMLScene, not exposed to Python
 *
 * @param doc          the xml document (supplies factory functions)
 * @param parent       the parent xml element, where the organ's tag is added
 */
void Organ::writeRSML(tinyxml2::XMLDocument& doc, tinyxml2::XMLElement* parent) const
{
	if (this->nodes.size()>1) {
		int nn = plant.lock()->getRSMLSkip()+1;
		// organ
		// std::string name = getOrganTypeParameter()->name; // todo where to put it
		tinyxml2::XMLElement* organ = doc.NewElement("root"); // TODO use ot to fetch tag name?
		organ->SetAttribute("ID", id);
		// geometry
		tinyxml2::XMLElement* geometry = doc.NewElement("geometry");
		organ->InsertEndChild(geometry);
		tinyxml2::XMLElement* polyline = doc.NewElement("polyline");


		int o;
		if (this->parent.expired()) { // baseRoot = 0, others = 1
		    // std::cout << this->toString() << std::flush;
		    o = 0;
		} else {
		    if (this->parent.lock()->organType()==Organism::ot_seed) {
		        o = 0;
		    } else {
		        o = 1;
		    }
		}
		for (int i = o; i<getNumberOfNodes(); i+=nn) {
			auto n = getNode(i);
			tinyxml2::XMLElement* p = doc.NewElement("point");
			p->SetAttribute("x", float(n.x));
			p->SetAttribute("y", float(n.y));
			p->SetAttribute("z", float(n.z));
			polyline->InsertEndChild(p);
		}
		geometry->InsertEndChild(polyline);
		// properties
		tinyxml2::XMLElement* properties = doc.NewElement("properties");
		auto prop_names = plant.lock()->getRSMLProperties();
		for (const auto& pname : prop_names) {
			tinyxml2::XMLElement* p = doc.NewElement(pname.c_str());
			p->SetAttribute("value", float(this->getParameter(pname)));
			properties->InsertEndChild(p);
		}
		organ->InsertEndChild(properties);
		/* laterals roots */
		for (size_t i = 0; i<children.size(); i+=nn) {
			children[i]->writeRSML(doc, organ);
		}
		// functions
		tinyxml2::XMLElement* fcts = doc.NewElement("functions");
		tinyxml2::XMLElement* fun1 = doc.NewElement("function");
		fun1->SetAttribute("domain","polyline");
		fun1->SetAttribute("name","node_creation_time");
		for (int i = o; i<getNumberOfNodes(); i+=nn) {
			double ct = getNodeCT(i);
			tinyxml2::XMLElement* p = doc.NewElement("sample");
			p->SetAttribute("value", ct);
			fun1->InsertEndChild(p);
		}
		tinyxml2::XMLElement* fun2 = doc.NewElement("function");
		fun2->SetAttribute("domain","polyline");
		fun2->SetAttribute("name","node_index");
		for (int i = o; i<getNumberOfNodes(); i+=nn) {
			int nid = getNodeId(i);
			tinyxml2::XMLElement* p = doc.NewElement("sample");
			p->SetAttribute("value", nid);
			fun2->InsertEndChild(p);
		}
		fcts->InsertEndChild(fun1);
		fcts->InsertEndChild(fun2);
		organ->InsertEndChild(fcts);
		parent->InsertEndChild(organ);
	}
}

/**
 * @return Quick info about the object for debugging,
 * additionally, use getParam()->toString() and getOrganRandomParameter()->toString() to obtain all information.
 */
std::string Organ::toString() const
{
	std::stringstream str;
	str << Organism::organTypeNames.at(this->organType()) << " #"<< getId() <<": sub type "<< param_->subType
					<< ", realized length " << getLength(true)
					<< " cm , theoretic length " << getLength(false) << " cm , age " << getAge()
    				<< " days, alive " << isAlive() << ", active " << isActive() << ", number of nodes " << this->getNumberOfNodes()
					<< ", with "<< children.size() << " children";
	return str.str();
}

/**
 * @return Current absolute heading of the organ tip, based on initial heading, or last two nodes
 */
Vector3d Organ::heading() const
{
	if (nodes.size()>1) {
		int n = nodes.size();
		Vector3d h = getNode(n-1).minus(getNode(n-2));
		h.normalize();
		return h;
	} else {
		return getiHeading0();
	}
}

/**
 * Needed for carbon-limited growth: to know sucrose necessary for length increase
 * Overwritten by @Leaf::orgVolume
 * @param length   length for which volume is calculated. if length = -1, use current organ length
 * @param realized for length = -1: use current theoratical or realized length
 * @return Volume for specific or current length. Overriden for @Leaf::orgVolume
 */
double Organ::orgVolume(double length_,  bool realized) const
{
	if(length_ == -1){length_ = getLength(realized);}
	return M_PI * length_ * param_->a * param_->a;//cylinder
};
}
