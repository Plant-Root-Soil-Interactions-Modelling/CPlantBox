// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Organ.h"

#include "Organism.h"
#include "Plant.h"
#include <algorithm>
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
    double age, double length, Vector3d partialIHeading_, int pni, bool moved, int oldNON)
    : parentNI(pni),partialIHeading(partialIHeading_), plant(), parent(), id(id), param_(param), alive(alive), active(active), age(age),
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
 * @param pni       parent node index
 */
Organ::Organ(std::shared_ptr<Organism> plant, std::shared_ptr<Organ> parent, int ot, int st, double delay,
    int pni)
: parentNI(pni), plant(plant), parent(parent), id(plant->getOrganIndex()),
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
 * @return the organs length from start node up to the node with index @param i.
 */
double Organ::getLength(int i) const
{
    double l = 0.; // length until node i
    if(hasRelCoord()){//is currently using relative coordinates?
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
double Organ::getLength(bool realized) const
{
    if (realized) {
        return length - this->epsilonDx;
    } else {
        return length;
    }
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
 * Adds a node to the organ. Overriden by @see Stem::addNode
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
    nodes.push_back(n); // node
    nodeIds.push_back(id); //unique id
    nodeCTs.push_back(t); // exact creation time
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
    //deprecated: do not need bulb anymore, stems of subtype 2 are normal stems
    //bool notBulb = !((this->organType() == Organism::ot_stem)&&(this->getParameter("subType") == 2));//do not count leaf bulb
    //might have age <0 and node.size()> 1 when adding organ manuelly @see test_organ.py
    bool forCarbon_limitedGrowth = (all && (this->getAge()>0));//when ask for "all" organs which have age > 0 even if nodes.size() == 1
    bool notSeed = ( this->organType() != Organism::ot_seed);

    if ((this->nodes.size()>1 || forCarbon_limitedGrowth) && notSeed) {//&& notBulb
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
    if (name=="length") { return getLength(true); } // realized organ length, dependent on dxMin and dx
    if (name=="lengthTh") { return getLength(false); } // theoratical organ length, dependent on dxMin and dx
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
 * additionally, use param()->toString() and getOrganRandomParameter()->toString() to obtain all information.
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
 * convert the nodes' positions from relative to absolute coordinates
 */
void Organ::rel2abs()
{
    if(hasRelCoord())
    {
        nodes[0] = getOrigin(); //recompute position of the first node

        for(size_t i=1; i<nodes.size(); i++)
        {
            double sdx = nodes[i].length();
            Vector3d newdx = getIncrement(nodes[i-1], sdx, i-1); //add tropism
            nodes[i] = nodes[i-1].plus(newdx); //replace relative by absolute position
        }
        moved = true; //update position of existing nodes in MappedSegments
    }
    //if carry children, update their pos

    for(size_t i=0; i<children.size(); i++){
        (children[i])->rel2abs();//even if parent does not have relCoordinate, the laterals might
    }
}

/**
 *  convert the nodes' positions from absolute to relative coordinates
 */
void Organ::abs2rel()
{
    bool isShoot = ((organType()==Organism::ot_stem)||(organType()==Organism::ot_leaf));
    if(isShoot||(getParent()->hasRelCoord()))//convert to relative coordinate if is shoot organ or carried by shoot organs
    {
        for (int j = nodes.size(); j>1; j--) {
            double sdx = (nodes.at(j-1).minus(nodes.at(j-2))).length();
            nodes.at(j-1) = Vector3d(sdx,0.,0.);
            //nodes.at(j-1) = nodes.at(j-1).minus(nodes.at(j-2));
        }
        nodes[0] = Vector3d(0.,0.,0.);
        moved = true; //update position of existing nodes in MappedSegments
    }
    for(size_t i=0; i<children.size(); i++){
        //if((children[i])->organType()!=Organism::ot_root){
        (children[i])->abs2rel();
        //}
    }//if carry children, update their pos

}

/**
 * @return Current absolute heading of the organ at node n, based on initial heading, or segment before
 */
Vector3d Organ::getiHeading0()  const
{
    if (!getParent()) { // in case of class RootSystem base roots (tap, basal, shootborne) or Organism organs created manually have no parent
        Vector3d vparentHeading ;
        if(organType() <= Organism::ot_root){//root (2) or unrecognized organ (0)
            vparentHeading = Vector3d(0, 0, -1);
        }else{//stem (3) or leaf(4)
            vparentHeading = Vector3d(0, 0, 1);
        }
        Matrix3d parentHeading = Matrix3d::ons(vparentHeading);
        return parentHeading.times(this->partialIHeading);
    }
    Matrix3d parentHeading;
    bool isBaseOrgan = (getParent()->organType()==Organism::ot_seed);
    bool isShootBornRoot = ((getParent()->organType()==Organism::ot_stem)&&(organType()==Organism::ot_root));
    if (isBaseOrgan||isShootBornRoot) { // from seed?
        if (organType()==Organism::ot_root) {
            parentHeading = Matrix3d(Vector3d(0, 0, -1), Vector3d(0, -1, 0), Vector3d(-1, 0, 0));
        }else{
            parentHeading = Matrix3d(Vector3d(0, 0, 1), Vector3d(0, 1, 0), Vector3d(1, 0, 0));
        }
    } else {
        Vector3d vparentHeading = getParent()->heading(parentNI);
        parentHeading = Matrix3d::ons(vparentHeading);
    }
    auto heading = parentHeading.column(0);
    Vector3d new_heading = Matrix3d::ons(heading).times(this->partialIHeading);
    return Matrix3d::ons(new_heading).column(0);
}
/**
 * @return Current absolute heading of the organ at node n, based on initial heading, or direction of the segment going from node n-1 to node n
 */
Vector3d Organ::heading(int n ) const
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

/**
 * Returns the increment of the next segments
 *
 *  @param p       coordinates of previous node
 *  @param sdx     length of next segment [cm]
 *  @return        the vector representing the increment
 */
Vector3d Organ::getIncrement(const Vector3d& p, double sdx, int n)
{
    Vector3d h = heading(n);
    Matrix3d ons = Matrix3d::ons(h);
    Vector2d ab = getOrganRandomParameter()->f_tf->getHeading(p, ons, dx(), shared_from_this(), n+1); // << THIS causes segmentation fault if f_tf->plant is expired
    //for leaves: necessary?
    //Vector2d ab = getLeafRandomParameter()->f_tf->getHeading(p, ons, dx(),shared_from_this());
    Vector3d sv = ons.times(Vector3d::rotAB(ab.x,ab.y));
    return sv.times(sdx);
}

/**
 * Analytical creation (=emergence) time of a point along the already grown root
 *
 * @param length   length along the organ, where the point is located [cm]
 * @param dt 	   current time step [day]
 * @return         the analytic time when this point was reached by the growing root [day],
 * 				   if growth is impeded, the value is not exact, but approximated dependent on the temporal resolution.
 */
double Organ::calcCreationTime(double length, double dt)
{
    assert(length >= 0 && "Organ::getCreationTime() negative length");
    double age_ = calcAge(std::max(length,0.)); // organ age as if grown unimpeded (lower than real age)
    double a = std::max(age_, age-dt /*old age*/);
    a = std::min(a, age); // a in [age-dt, age]
    assert((a+nodeCTs[0] >= 0) && "Organ::getCreationTime() negative creation time");
    return a+nodeCTs[0];
}


/**
 *  Creates nodes and node emergence times for a length l
 *
 *  Checks that each new segments length is <= dx but >= parent->minDx
 *
 *  @param l        total length of the segments that are created [cm]
 *  @param dt       time step [day]
 *  @param PhytoIdx index of phytomere node to elongate (optional) [1]
 *  @param verbose  turns console output on or off
 */
void Organ::createSegments(double l, double dt, bool verbose, int PhytoIdx)
{
    if (l==0) {
        std::cout << "Organ::createSegments: zero length encountered \n";
        return;
    }
    if (l<0) {
        std::cout << "Organ::createSegments: negative length encountered \n";
    }
    // shift first node to axial resolution
    double shiftl = 0; // length produced by shift
    int nn = nodes.size();
    bool stemElongation = (PhytoIdx >= 0);//if we are doing internodal growth,  PhytoIdx >= 0.
    if( stemElongation){
        nn = PhytoIdx +1;
    }
    if (firstCall||stemElongation) { // first call of createSegments (in Organ::simulate)
        if (!stemElongation) {
            firstCall = false;
        }
        //don t move first node
        bool notFirstNode =  (nn>1);
        bool notChildBaseNode = (children.empty() || (nn-1 != std::static_pointer_cast<Organ>(children.back())->parentNI));
        // don't move a child base node unless you have phytomere expansion
        if (notFirstNode && (notChildBaseNode || stemElongation) ) {
            Vector3d h;
            Vector3d n2 = nodes.at(nn-2);

            if (hasRelCoord()) {
                h = nodes.at(nn-1);
            } else {
                Vector3d n1 = nodes.at(nn-1);
                h = n1.minus(n2);
            }
            double olddx = h.length(); // length of last segment
            if (olddx<dx()*(1 - 1e-10)) { // shift node instead of creating a new node
                shiftl = std::min(dx()-olddx, l);
                double sdx = olddx + shiftl; // length of new segment
                // Vector3d newdxv = getIncrement(n2, sdx);

                if(hasRelCoord()) {
                    nodes.at(nn-1) =  Vector3d(sdx,0.,0.);//h.times(sdx);
                }else{
                    h.normalize();
                    nodes.at(nn-1) = Vector3d(n2.plus(h.times(sdx))); // n2.plus(newdxv)
                }
                double et = this->calcCreationTime(getLength(true)+shiftl, dt);
                nodeCTs.at(nn-1) = et; // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
                moved = true;
                l -= shiftl;
                if (l<=0) { // ==0 should be enough
                    return;
                }
            } else {
                moved = false;
            }
        } else {
            moved = false;
        }
    }
    // create n+1 new nodes
    double sl = 0; // summed length of created segment
    int n = floor(l/dx());
    for (int i = 0; i < n + 1; i++) {

        double sdx; // segment length (<=dx)
        if (i<n) {  // normal case
            sdx = dx();
        } else
        { // last segment
            sdx = l-n*dx();
            if (sdx<dxMin()*(1-1e-10)) { //plant.lock()->getMinDx()) { // quit if l is too small
                //                if (verbose && (sdx != 0)) {
                //                    std::cout <<"Organ::createSegments(): "<<organType()<<" length increment below dxMin threshold ("<< sdx <<" < "<< dxMin() << ") and kept in memory\n";
                //                }
                if( PhytoIdx >= 0){
                    this->epsilonDx += sdx;
                }else{this->epsilonDx = sdx;}
                return;
            }
            this->epsilonDx = 0; //no residual
        }
        sl += sdx;
        Vector3d newnode;
        if (hasRelCoord()) {
            newnode = Vector3d(sdx, 0., 0.);
        } else {
            Vector3d newdx = getIncrement(nodes.back(), sdx);
            newnode = Vector3d(nodes.back().plus(newdx));
        }

        double et = this->calcCreationTime(getLength(true)+shiftl+sl, dt);//here length or get length? it s the same because epsilonDx was set back to 0 at beginning of simulate no?
        // in case of impeded growth the node emergence time is not exact anymore,
        // but might break down to temporal resolution
        addNode(newnode, et,size_t(nn+i),stemElongation);
    }
}


/**
 * Creates a new lateral
 *  @param dt       time step [day]
 *  @param verbose  turns console output on or off
 */
void Organ::createLateral(double dt, bool verbose)
{
    auto rp = getOrganRandomParameter(); // rename

    for(int i = 0; i < rp->successorST.size(); i++){//go through each successor rule
        //found id
        bool applyHere = getApplyHere(i);

        if(applyHere)
        {
            int numlats = 1;//how many laterals? default = 1
            if(rp->successorNo.size()>i){numlats =  rp->successorNo.at(i);}
            for(int nn = 0; nn < numlats; nn++)
            {

                const Vector3d& pos = Vector3d();
                int p_id = rp->getLateralType(pos, i);//if probabilistic branching

                if(p_id >=0)
                {
                    int ot;

                    if((rp->successorOT.size()>i)&&(rp->successorOT.at(i).size()>p_id)){
                        ot = rp->successorOT.at(i).at(p_id);
                    }else{ot = getParameter("organType");}//default

                    int st = rp->successorST.at(i).at(p_id);

                    double delay = getLatGrowthDelay(ot, st, dt);// forDelay*multiplyDelay
                    double growth_dt = getLatInitialGrowth(dt);


                    switch(ot){
                    case Organism::ot_root:{
                        auto lateral = std::make_shared<Root>(plant.lock(), st,  delay, shared_from_this(),  nodes.size() - 1);
                        children.push_back(lateral);
                        lateral->simulate(growth_dt,verbose);
                        break;}
                    case Organism::ot_stem:{
                        auto lateral = std::make_shared<Stem>(plant.lock(), st, delay, shared_from_this(),  nodes.size() - 1);
                        children.push_back(lateral);
                        lateral->simulate(growth_dt,verbose);
                        break;}
                    case Organism::ot_leaf:{
                        auto lateral = std::make_shared<Leaf>(plant.lock(), st,  delay, shared_from_this(),  nodes.size() - 1);
                        children.push_back(lateral);
                        lateral->simulate(growth_dt,verbose);//age-ageLN,verbose);
                        break;}
                    }
                }
            }
        }

    }
    created_linking_node ++;
    storeLinkingNodeLocalId(created_linking_node,verbose);//needed (currently) only for stems when doing nodal growth

}


/**
 * See if should apply successor rule at a specific linking node, @see Organ::createLateral
 *  @param i       rule id
 *  @return whether to apply the rule
 */
bool Organ::getApplyHere(int i) const
{
    bool applyHere;
    auto rp = getOrganRandomParameter(); // rename
    if((rp->successorWhere.size()>i)&&(rp->successorWhere.at(i).size()>0)){
        if(!std::signbit(rp->successorWhere.at(i).at(0)))//true if number is signed
        {//gave which linking nodes to include
            applyHere = (std::find (rp->successorWhere.at(i).begin(), rp->successorWhere.at(i).end(), created_linking_node)
                != rp->successorWhere.at(i).end());
        }else{//gave which linking nodes to ignore
            applyHere = !(std::find (rp->successorWhere.at(i).begin(), rp->successorWhere.at(i).end(), -double(created_linking_node))
                != rp->successorWhere.at(i).end());
        }

    }else{applyHere = true;} //default
    return applyHere;
}


/**
 *  @see Organ::createLateral
 *  @param dt       time step recieved by parent organ [day]
 *  @return growth period to send to lateral after creation
 */
double Organ::getLatInitialGrowth(double dt)
{
    double ageLN = this->calcAge(getLength(true)); // MINIMUM age of root when lateral node is created
    ageLN = std::max(ageLN, age-dt);
    return age-ageLN;
}


/**
 *  @see Organ::createLateral
 *  @param ot_lat       organType of lateral to create
 *  @param st_lat       subType of lateral to create
 *  @param dt       time step recieved by parent organ [day]
 *  @return emergence delay to send to lateral after creation
 */
double Organ::getLatGrowthDelay(int ot_lat, int st_lat, double dt) const //override for stems
{
    auto rp = getOrganRandomParameter(); // rename
    double growthDelay; //store necessary variables to define lateral growth delay
    int delayDefinition = getOrganism()->getDelayDefinition(ot_lat);

    assert(delayDefinition >= 0);

    switch(delayDefinition){
      case Organism::dd_distance:
      {
          double meanLn = getParameter("lnMean"); // mean inter-lateral distance
          double effectiveLa = std::max(getParameter("la")-meanLn/2, 0.); // effective apical distance, observed apical distance is in [la-ln/2, la+ln/2]
          double ageLN = this->calcAge(getLength(true)); // theoretical age of root when lateral node is created
          ageLN = std::max(ageLN, age-dt);
          double ageLG = this->calcAge(getLength(true)+effectiveLa); // age of the root, when the lateral starts growing (i.e when the apical zone is developed)
          growthDelay = ageLG-ageLN; // time the lateral has to wait
          break;
      }
      case Organism::dd_time_lat:
      {
          // time the lateral has to wait
          growthDelay = std::max(rp->ldelay + plant.lock()->randn()*rp->ldelays, 0.);
          break;
      }
      case Organism::dd_time_self:
      {

          //get delay per lateral
          auto latRp = plant.lock()->getOrganRandomParameter(ot_lat, st_lat); // random parameter of lateral to create
          growthDelay = std::max(latRp->ldelay + plant.lock()->randn()*latRp->ldelays, 0.);
          break;
      }
      default:
      {
          std::cout<<"delayDefinition "<<delayDefinition<<" "<<Organism::dd_distance<<" ";
          std::cout<< Organism::dd_time_lat<<" "<< Organism::dd_time_self<<std::endl<<std::flush;
          std::cout<<"				"<<(delayDefinition==Organism::dd_distance)<<" ";
          std::cout<<(delayDefinition== Organism::dd_time_lat)<<" "<< (delayDefinition==Organism::dd_time_self)<<std::endl<<std::flush;
          throw std::runtime_error("Delay definition type (delayDefinition) not recognised");
      }
    }
    return growthDelay;
}
/**
 * Check if the organ has relative coordinates.
 * meaning: organ is not a basal/tap root and first node is at (0,0,0) but
 */
bool Organ::hasRelCoord() const
{
    bool nullNode0 = (nodes.at(0) == Vector3d(0.,0.,0.)); // rel coordinate for node 0?
    bool isSeed = organType() == Organism::ot_seed;
    bool basalOrgan = true;
    if (getParent()) { // in case of class RootSystem base roots (tap, basal, shootborne) or Organism organs created manually have no parent
        basalOrgan = (isSeed||(getParent()->organType() == Organism::ot_seed));
    }
    bool isBasalRoot = ((organType() == Organism::ot_root)&&basalOrgan);
    return (nullNode0&&(!isBasalRoot)&&(!isSeed));
}


} // namespace CPlantBox
