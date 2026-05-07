// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// own header
#include "Organ.h"

// project headers
#include "Organism.h"
#include "Plant.h"
#include "organparameter.h"

// standard library
#include <algorithm>
#include <iostream>

namespace CPlantBox {

/**
 * @brief Constructs an organ from explicit state.
 *
 * Tree links and geometry are expected to be set consistently by the caller.
 *
 * @param id Organ identifier.
 * @param param Organ-specific parameters.
 * @param alive Initial alive state.
 * @param active Initial active state.
 * @param age Initial age [day].
 * @param length Initial length [cm].
 * @param partialIHeading_ Initial local heading.
 * @param pni Parent node index in the parent organ.
 * @param moved Whether nodes moved in the previous time step.
 * @param oldNON Number of nodes before the previous step.
 */
Organ::Organ(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length, Vector3d partialIHeading_,
             int pni, bool moved, int oldNON)
    : parentNI(pni), partialIHeading(partialIHeading_), plant(), parent(), id(id), param_(param), alive(alive), active(active), age(age), length(length),
      moved(moved), oldNumberOfNodes(oldNON) {}

/**
 * @brief Simulation constructor.
 *
 * Parameters are drawn via OrganRandomParameter::realize(); the organ starts with
 * age = -delay so it begins growing after @p delay days.
 *
 * @param plant   Organism that owns this organ.
 * @param parent  Parent organ (nullptr for base organs).
 * @param ot      Organ type.
 * @param st      Sub-type index.
 * @param delay   Growth delay [days]; organ age is initialised to -delay.
 * @param pni     Local node index in the parent where this organ attaches.
 */
Organ::Organ(std::shared_ptr<Organism> plant, std::shared_ptr<Organ> parent, int ot, int st, double delay, int pni)
    : parentNI(pni), plant(plant), parent(parent), id(plant->getOrganIndex()),
      param_(plant->getOrganRandomParameter(ot, st)->realize()), /* root parameters are diced in the getOrganRandomParameter class */
      age(-delay) {}

/**
 * @brief Deep-copies this organ and its subtree into another organism.
 * @param p Target organism for the copied organ tree.
 * @return Root of the copied subtree.
 */
std::shared_ptr<Organ> Organ::copy(std::shared_ptr<Organism> p) {
    auto o = std::make_shared<Organ>(*this); // shallow copy
    o->parent = std::weak_ptr<Organ>();
    o->plant = p;
    o->param_ = std::make_shared<OrganSpecificParameter>(*param_); // copy parameters
    for (size_t i = 0; i < children.size(); i++) {
        o->children[i] = children[i]->copy(p); // copy lateral
        o->children[i]->setParent(o);
    }
    return o;
}

/**
 * @brief Returns length from the first node up to node index @p i.
 *
 * Works correctly in both relative and absolute coordinate modes.
 *
 * @param i  Local node index (exclusive upper bound).
 * @return   Length [cm].
 */
double Organ::getLength(int i) const {
    double l = 0.;       // length until node i
    if (hasRelCoord()) { // is currently using relative coordinates?
        for (int j = 0; j < i; j++) {
            l += nodes.at(j + 1).length(); // relative length equals absolute length
        }
    } else {
        for (int j = 0; j < i; j++) {
            l += nodes.at(j + 1).minus(nodes.at(j)).length(); // relative length equals absolute length
        }
    }
    return l;
}

/**
 * @brief Returns either realized or theoretical organ length.
 * @param realized If true, subtracts unresolved residual growth epsilonDx.
 * @return Organ length [cm].
 */
double Organ::getLength(bool realized) const {
    if (realized) {
        return length - this->epsilonDx;
    } else {
        return length;
    }
}

/**
 * @brief Returns the organ type.
 *
 * Coarse classification: ot_organ=0, ot_seed=1, ot_root=2, ot_stem=3, ot_leaf=4.
 * Override in each derived class.
 * For string names see Organism::organTypeNames.
 */
int Organ::organType() const { return Organism::ot_organ; }

/**
 * @brief Returns the stochastic parameter set for this organ's sub-type.
 *
 * Delegates to the owning organism, which manages all OrganRandomParameter sets.
 */
std::shared_ptr<OrganRandomParameter> Organ::getOrganRandomParameter() const {
    return plant.lock()->getOrganRandomParameter(this->organType(), param_->subType);
}

/**
 * @brief Advances organ development by @p dt days.
 *
 * Increments age and recursively simulates all children.
 * Override in derived classes to implement organ-specific growth (segment creation, lateral initiation).
 *
 * @param dt       Time step [days].
 * @param verbose  If true, print diagnostic output.
 */
void Organ::simulate(double dt, bool verbose) {
    // store information of this time step
    oldNumberOfNodes = nodes.size();
    moved = false;

    // if the organ is alive, manage children
    if (alive) {
        age += dt;
        for (auto &c : children) {
            c->simulate(dt, verbose);
        }
    }
}

/**
 * @brief Downcast the owning organism to Plant.
 * @return Shared pointer to the Plant, or nullptr if the organism is not a Plant.
 */
std::shared_ptr<Plant> Organ::getPlant() const { return std::dynamic_pointer_cast<Plant>(plant.lock()); }

/**
 * @brief Adds a child organ and sets its parent pointer.
 * @param c Child organ.
 */
void Organ::addChild(std::shared_ptr<Organ> c) {
    c->setParent(shared_from_this());
    children.push_back(c);
}

/**
 * @brief Adds a node with an explicit global id.
 *
 * Nodes cannot be deleted; organs become deactivated or die instead.
 * Override in Stem for phytomere (internodal) handling.
 *
 * @param n      Node coordinate.
 * @param id     Global node index.
 * @param t      Creation time [days].
 * @param index  Position at which to insert the node.
 * @param shift  True if the new node is inserted between existing nodes (internodal growth).
 */
void Organ::addNode(Vector3d n, int id, double t, size_t index, bool shift) {
    nodes.push_back(n);    // node
    nodeIds.push_back(id); // unique id
    nodeCTs.push_back(t);  // exact creation time
}

/**
 * @brief Updates the parent attachment node index.
 *
 * Used during internodal (stem) growth.  Also updates the first nodeId to match
 * the new attachment node.
 *
 * @param idx  New local parent node index.
 */
void Organ::moveOrigin(int idx) {
    this->parentNI = idx;
    nodeIds.at(0) = getParent()->getNodeId(idx);
}

/**
 * @brief Adds a node using the next global index from the organism.
 *
 * Convenience wrapper around addNode(n, id, t, index, shift) that fetches
 * the next available global node index from the owning organism.
 *
 * @param n      Node coordinate.
 * @param t      Creation time [days].
 * @param index  Insertion position.
 * @param shift  True if inserting between existing nodes.
 */
void Organ::addNode(Vector3d n, double t, size_t index, bool shift) { addNode(n, plant.lock()->getNodeIndex(), t, index, shift); }

/**
 * @brief Returns all polyline segments as pairs of global node indices.
 *
 * The polyline {n0, n1, n2, n3} produces segments {[i0,i1], [i1,i2], [i2,i3]}.
 * Returns an empty vector if fewer than two nodes exist.
 */
std::vector<Vector2i> Organ::getSegments() const {
    if (this->nodes.size() > 1) {
        std::vector<Vector2i> segs = std::vector<Vector2i>(nodes.size() - 1);
        for (size_t i = 0; i < nodes.size() - 1; i++) {
            Vector2i s(getNodeId(i), getNodeId(i + 1));
            segs[i] = s;
        }
        return segs;
    } else {
        return std::vector<Vector2i>(0);
    }
}

/**
 * @brief Returns the maximal axial resolution dx.
 */
double Organ::dx() const { return getOrganRandomParameter()->dx; }

/**
 * @brief Returns the minimal axial resolution dxMin.
 *
 * Residual growth below this threshold is accumulated in epsilonDx.
 */
double Organ::dxMin() const { return getOrganRandomParameter()->dxMin; }

/**
 * @brief Returns this organ and all descendants as a flat list.
 *
 * Only organs with more than one node are included unless @p all is true.
 * If @p all is true, organs with age > 0 (but possibly only one node, e.g. carbon-limited growth)
 * are also included, except for seeds.
 *
 * @param ot   Filter by organ type; -1 returns all types.
 * @param all  If true, also include organs with a single node (age > 0).
 */
std::vector<std::shared_ptr<Organ>> Organ::getOrgans(int ot, bool all) {
    auto v = std::vector<std::shared_ptr<Organ>>();
    this->getOrgans(ot, v, all);
    return v;
}

/**
 * @brief Appends this organ and all descendants to @p v (in-place variant).
 *
 * @param ot   Filter by organ type; -1 returns all types.
 * @param v    Output vector; qualifying organs are appended.
 * @param all  If true, also include organs with a single node (age > 0).
 */
void Organ::getOrgans(int ot, std::vector<std::shared_ptr<Organ>> &v, bool all) {
    // deprecated: do not need bulb anymore, stems of subtype 2 are normal stems
    // bool notBulb = !((this->organType() == Organism::ot_stem)&&(this->getParameter("subType") == 2));//do not count leaf bulb
    // might have age <0 and node.size()> 1 when adding organ manuelly @see test_organ.py
    bool forCarbon_limitedGrowth = (all && (this->getAge() > 0)); // when ask for "all" organs which have age > 0 even if nodes.size() == 1
    bool notSeed = (this->organType() != Organism::ot_seed);

    if ((this->nodes.size() > 1 || forCarbon_limitedGrowth) && notSeed) { //&& notBulb
        if ((ot < 0) || (ot == this->organType())) {
            v.push_back(shared_from_this());
        }
    }
    // std::cout << "Organ::getOrgans recursive: number of children " <<  this->children.size() << "\n" << std::flush;
    for (const auto &c : this->children) {
        c->getOrgans(ot, v, all);
    }
}

/**
 * @brief Returns the number of children whose age > 0 (i.e. emerged laterals).
 * @see getNumberOfChildren()
 */
int Organ::getNumberOfLaterals() const {
    int nol = 0;
    for (auto &c : children) {
        if (c->getAge() > 0) { // born
            nol++;
        }
    }
    return nol;
}

/**
 * @brief Returns a named scalar parameter of the organ.
 *
 * Recognises organ-specific fields (subType, a, radius, diameter, iHeadingX/Y/Z, parentNI,
 * parent-node), member variables (organType, id, alive, active, age, length, lengthTh,
 * numberOfNodes, numberOfSegments, hasMoved, oldNumberOfNodes, numberOfLaterals, creationTime,
 * order, one), and delegates to OrganRandomParameter::getParameter() for anything else.
 *
 * Intended for post-processing; flexible but slower than direct member access.
 * Override in derived classes to expose additional parameters.
 */
double Organ::getParameter(std::string name) const {
    // specific parameters
    if (name == "subType") {
        return this->param_->subType;
    }
    if (name == "a") { // root radius [cm]
        return param_->a;
    }
    if (name == "radius") { // root radius [cm]
        return this->param_->a;
    }
    if (name == "diameter") { // root diameter [cm]
        return 2. * this->param_->a;
    }
    // organ member variables
    if (name == "iHeadingX") { // root initial heading x - coordinate [cm]
        return getiHeading0().x;
    }
    if (name == "iHeadingY") { // root initial heading y - coordinate [cm]
        return getiHeading0().y;
    }
    if (name == "iHeadingZ") { // root initial heading z - coordinate [cm]
        return getiHeading0().z;
    }
    if (name == "parentNI") { // local parent node index where the lateral emerges
        return parentNI;
    }
    if (name == "parent-node") { // local parent node index for RSML (higher order roots are missing the first node)
        if (this->parent.expired()) {
            return -1;
        }
        if (this->parent.lock()->organType() == Organism::ot_seed) { // if it is base root
            return -1;
        }
        auto p = this->parent.lock();
        if (p->parent.expired()) { // if parent is base root
            return parentNI;
        }
        if (p->parent.lock()->organType() == Organism::ot_seed) { // if parent is base root
            return parentNI;
        } else {
            return std::max(parentNI - 1, 0); // higher order roots are missing the first node
        }
    }
    // organ member functions
    if (name == "organType") {
        return this->organType();
    }
    if (name == "numberOfChildren") {
        return children.size();
    }
    if (name == "id") {
        return getId();
    }
    if (name == "alive") {
        return isAlive();
    }
    if (name == "active") {
        return isActive();
    }
    if (name == "age") {
        return getAge();
    }
    if (name == "length") { // realized organ length, dependent on dxMin and dx
        return getLength(true);
    }
    if (name == "lengthTh") {
        return getLength(false);
    } // theoratical organ length, dependent on dxMin and dx
    if (name == "numberOfNodes") {
        return getNumberOfNodes();
    }
    if (name == "numberOfSegments") {
        return getNumberOfSegments();
    }
    if (name == "hasMoved") {
        return hasMoved();
    }
    if (name == "oldNumberOfNodes") {
        return getOldNumberOfNodes();
    }
    if (name == "numberOfLaterals") {
        return getNumberOfLaterals();
    }
    // further
    if (name == "creationTime") {
        return getNodeCT(0);
    }
    if (name == "order") { // count how often it is possible to move up
        int o = 0;
        auto p = shared_from_this();
        while ((!p->parent.expired()) && (p->parent.lock()->organType() != Organism::ot_seed)) {
            o++;
            p = p->parent.lock(); // up the organ tree
        }
        return o;
    }
    if (name == "one") {
        return 1;
    } // e.g. for counting the organs
    return this->getOrganRandomParameter()->getParameter(name); // ask the random parameter
}

/**
 * @brief Writes an RSML @c <root> element for this organ.
 *
 * Skips organs with fewer than two nodes.  Writes geometry (polyline), properties,
 * lateral children, and functions (node_creation_time, node_index).
 * Called by Organism::getRSMLScene(); not exposed to Python.
 *
 * @param doc     XML document supplying element factory functions.
 * @param parent  Parent XML element to which the new @c <root> tag is appended.
 */
void Organ::writeRSML(tinyxml2::XMLDocument &doc, tinyxml2::XMLElement *parent) const {
    if (this->nodes.size() > 1) {
        int nn = plant.lock()->getRSMLSkip() + 1;
        // organ
        // std::string name = getOrganTypeParameter()->name; // todo where to put it
        tinyxml2::XMLElement *organ = doc.NewElement("root"); // TODO use ot to fetch tag name?
        organ->SetAttribute("ID", id);
        // geometry
        tinyxml2::XMLElement *geometry = doc.NewElement("geometry");
        organ->InsertEndChild(geometry);
        tinyxml2::XMLElement *polyline = doc.NewElement("polyline");
        int o;
        if (this->parent.expired()) { // baseRoot = 0, others = 1
            // std::cout << this->toString() << std::flush;
            o = 0;
        } else {
            if (this->parent.lock()->organType() == Organism::ot_seed) {
                o = 0;
            } else {
                o = 1;
            }
        }
        for (int i = o; i < getNumberOfNodes(); i += nn) {
            auto n = getNode(i);
            tinyxml2::XMLElement *p = doc.NewElement("point");
            p->SetAttribute("x", float(n.x));
            p->SetAttribute("y", float(n.y));
            p->SetAttribute("z", float(n.z));
            polyline->InsertEndChild(p);
        }
        geometry->InsertEndChild(polyline);
        // properties
        tinyxml2::XMLElement *properties = doc.NewElement("properties");
        auto prop_names = plant.lock()->getRSMLProperties();
        for (const auto &pname : prop_names) {
            tinyxml2::XMLElement *p = doc.NewElement(pname.c_str());
            p->SetAttribute("value", float(this->getParameter(pname)));
            properties->InsertEndChild(p);
        }
        organ->InsertEndChild(properties);
        /* laterals roots */
        for (size_t i = 0; i < children.size(); i += nn) {
            children[i]->writeRSML(doc, organ);
        }
        // functions
        tinyxml2::XMLElement *fcts = doc.NewElement("functions");
        tinyxml2::XMLElement *fun1 = doc.NewElement("function");
        fun1->SetAttribute("domain", "polyline");
        fun1->SetAttribute("name", "node_creation_time");
        for (int i = o; i < getNumberOfNodes(); i += nn) {
            double ct = getNodeCT(i);
            tinyxml2::XMLElement *p = doc.NewElement("sample");
            p->SetAttribute("value", ct);
            fun1->InsertEndChild(p);
        }
        tinyxml2::XMLElement *fun2 = doc.NewElement("function");
        fun2->SetAttribute("domain", "polyline");
        fun2->SetAttribute("name", "node_index");
        for (int i = o; i < getNumberOfNodes(); i += nn) {
            int nid = getNodeId(i);
            tinyxml2::XMLElement *p = doc.NewElement("sample");
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
 * @brief Returns a compact human-readable summary for debugging.
 *
 * For full details call param()->toString() and getOrganRandomParameter()->toString().
 */
std::string Organ::toString() const {
    std::stringstream str;
    str << Organism::organTypeNames.at(this->organType()) << " #" << getId() << ": sub type " << param_->subType << ", realized length " << getLength(true)
        << " cm , theoretic length " << getLength(false) << " cm , age " << getAge() << " days, alive " << isAlive() << ", active " << isActive()
        << ", number of nodes " << this->getNumberOfNodes() << ", with " << children.size() << " children";
    return str.str();
}

/**
 * @brief Converts node coordinates from relative to absolute (recursive over children).
 *
 * Recomputes absolute positions using the parent heading at each node.
 * Sets moved=true so MappedSegments can update.
 */
void Organ::rel2abs() {
    if (hasRelCoord()) {
        nodes[0] = Vector3d(getOrigin()); // recompute position of the first node
        this->has_rel_coord = false;
        for (size_t i = 1; i < nodes.size(); i++) {
            bool doAddIncrement = false;
            if (addIncrement.size() > 0) {
                auto it = find(addIncrement.begin(), addIncrement.end(), i);
                if (it != addIncrement.end()) {
                    doAddIncrement = true;
                }
            }
            if (doAddIncrement) {
                Vector3d newdx = getIncrement(nodes.at(i - 1), nodes.at(i).length(), i - 1);
                nodes[i] = Vector3d(nodes.at(i - 1).plus(newdx));
            } else {
                Vector3d h = heading(i - 1);
                Matrix3d ons = Matrix3d::ons(h);
                nodes[i] = nodes[i - 1].plus(ons.times(nodes[i]));
            }
        }
        moved = true; // update position of existing nodes in MappedSegments
    }
    addIncrement.clear();

    // if carry children, update their pos
    for (size_t i = 0; i < children.size(); i++) {
        (children[i])->rel2abs(); // even if parent does not have relCoordinate, the laterals might
    }
}

/**
 * @brief Converts node coordinates from absolute to relative (recursive over children).
 *
 * Only affects shoot organs (stems and leaves) or root organs carried by a shoot parent.
 * Sets moved=true so MappedSegments can update.
 */
void Organ::abs2rel() {
    for (size_t i = 0; i < children.size(); i++) {
        (children[i])->abs2rel();
    } // if carry children, update their pos
    bool isShoot = ((organType() == Organism::ot_stem) || (organType() == Organism::ot_leaf));
    bool hasShootParent = false;
    if (!isShoot) {
        auto ot_parent = -1;
        auto parent_ = getParent();
        while ((ot_parent != 0) && (ot_parent <= 2) &&
               (parent_ != nullptr)) // we reach either a seed (finished going up parents) or shoot organ (need to be in relative coordinates)
        {
            ot_parent = parent_->organType();
            parent_ = parent_->getParent();
        }
        hasShootParent = ot_parent > 2;
    }

    if ((isShoot || hasShootParent) && (!hasRelCoord())) // convert to relative coordinate if is shoot organ or carried by shoot organs
    {
        for (int j = nodes.size(); j > 1; j--) {
            auto nodes_j_1 = nodes.at(j - 1).minus(nodes.at(j - 2));
            Vector3d h = heading(j - 2);
            Matrix3d onsinv = Matrix3d::ons(h).inverse();
            nodes.at(j - 1) = onsinv.times(nodes_j_1);
        }
        nodes[0] = Vector3d(0., 0., 0.);
        this->has_rel_coord = true;
        moved = true; // update position of existing nodes in MappedSegments
    }
}

/**
 * @brief Returns the initial absolute heading when the organ was created.
 *
 * Reconstructed from partialIHeading and the parent's heading at parentNI.
 * Handles base organs (no parent), shoot-borne roots, and root-borne shoots as special cases.
 */
Vector3d Organ::getiHeading0() const {
    if (!getParent()) { // in case of class RootSystem base roots (tap, basal, shootborne) or Organism organs created manually have no parent
        Vector3d vparentHeading;
        if (organType() <= Organism::ot_root) { // root (2) or unrecognized organ (0)
            vparentHeading = Vector3d(0, 0, -1);
        } else { // stem (3) or leaf(4)
            vparentHeading = Vector3d(0, 0, 1);
        }
        Matrix3d parentHeading = Matrix3d::ons(vparentHeading);
        return parentHeading.times(this->partialIHeading);
    }
    Matrix3d parentHeading;
    bool isBaseOrgan = (getParent()->organType() == Organism::ot_seed);
    bool isShootBornRoot = ((getParent()->organType() == Organism::ot_stem) && (organType() == Organism::ot_root));
    bool isRootBornShoot = ((getParent()->organType() == Organism::ot_root) && ((organType() == Organism::ot_stem) || (organType() == Organism::ot_leaf)));
    if (isBaseOrgan || isShootBornRoot || isRootBornShoot) { // from seed?
        if (organType() == Organism::ot_root) {
            parentHeading = Matrix3d(Vector3d(0, 0, -1), Vector3d(0, -1, 0), Vector3d(-1, 0, 0));
        } else {
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
 * @brief Returns the absolute heading at node @p n.
 *
 * For n > 0 and at least two nodes: direction of segment [n-1, n], normalised.
 * For n == 0 or single-node organs: falls back to getiHeading0().
 * Negative @p n is treated as the last node.
 */
Vector3d Organ::heading(int n) const {

    if (n < 0) {
        n = nodes.size() - 1;
    }
    if ((nodes.size() > 1) && (n > 0)) {
        n = std::min(int(nodes.size()), n);
        Vector3d h = getNode(n).minus(getNode(n - 1));
        h.normalize();
        return h;
    } else {
        return getiHeading0();
    }
}

/**
 * @brief Returns organ volume [cm³].
 *
 * Uses a cylinder approximation: V = π·r²·L.
 * Override in Leaf (and other non-cylindrical organs), e.g. for carbon limited growth
 *
 * @param length_   Length to use; -1 means use current organ length.
 * @param realized  When length_=-1: true uses realised length, false uses theoretical length.
 */
double Organ::orgVolume(double length_, bool realized) const {
    if (length_ == -1) {
        length_ = getLength(realized);
    }
    return M_PI * length_ * param_->a * param_->a; // cylinder
};

/**
 * @brief Returns the growth increment vector for the next segment.
 *
 * Queries the tropism function (f_tf) to obtain heading angles alpha/beta,
 * then rotates the heading frame accordingly.
 * Override in Leaf.
 *
 * @param p    Coordinates of the previous node.
 * @param sdx  Length of the next segment [cm].
 * @param n    Node index used to determine the current heading (-1 = last node).
 * @return     Increment vector of length @p sdx in the growth direction.
 */
Vector3d Organ::getIncrement(const Vector3d &p, double sdx, int n) {
    Vector3d h = heading(n);
    Matrix3d ons = Matrix3d::ons(h);
    Vector2d ab = getOrganRandomParameter()->f_tf->getHeading(p, ons, dx(), shared_from_this()); // << THIS causes segmentation fault if f_tf->plant is expired
    // for leaves: necessary?
    // Vector2d ab = getLeafRandomParameter()->f_tf->getHeading(p, ons, dx(),shared_from_this());
    Vector3d sv = ons.times(Vector3d::rotAB(ab.x, ab.y));
    return sv.times(sdx);
}

/**
 * @brief Returns the analytical creation time of a node at @p length along the organ.
 *
 * If growth was impeded the result is an approximation bounded to [age-dt, age].
 *
 * @param length  Distance along the organ from its base [cm].
 * @param dt      Current simulation time step [days].
 * @return        Absolute creation time [days].
 */
double Organ::calcCreationTime(double length, double dt) {
    assert(length >= 0 && "Organ::getCreationTime() negative length");
    double age_ = calcAge(std::max(length, 0.)); // organ age as if grown unimpeded (lower than real age)
    double a = std::max(age_, age - dt /*old age*/);
    a = std::min(a, age); // a in [age-dt, age]
    assert((a + nodeCTs[0] >= 0) && "Organ::getCreationTime() negative creation time");
    return a + nodeCTs[0];
}

/**
 * @brief Creates polyline segments totalling length @p l.
 *
 * On the first call (or during internodal growth) shifts the last existing node
 * if it is shorter than dx(), then appends new nodes at intervals of dx().
 * Residual growth smaller than dxMin() is stored in epsilonDx.
 *
 * @param l         Total length of new segments to create [cm].
 * @param dt        Current time step [days] (used for creation-time calculation).
 * @param verbose   If true, print diagnostic output.
 * @param phytoIdx  Node index for internodal elongation; -1 for tip growth.
 */
void Organ::createSegments(double l, double dt, bool verbose, int phytoIdx) {
    if (l == 0) {
        std::cout << "Organ::createSegments: zero length encountered \n";
        return;
    }
    if (l < 0) {
        std::cout << "Organ::createSegments: negative length encountered \n";
        return;
    }
    // shift first node to axial resolution
    double shiftl = 0; // length produced by shift
    int nn = nodes.size();
    bool stemElongation = (phytoIdx >= 0); // if we are doing internodal growth,  phytoIdx >= 0.
    if (stemElongation) {
        nn = phytoIdx + 1;
    }
    if (firstCall || stemElongation) { // first call of createSegments (in Organ::simulate) or elongation
        if (!stemElongation) {
            firstCall = false;
        }
        // don t move first node
        bool notFirstNode = (nn > 1);
        bool notChildBaseNode = (children.empty() || (nn - 1 != std::static_pointer_cast<Organ>(children.back())->parentNI));
        // don't move a child base node unless you have phytomere expansion
        if (notFirstNode && (notChildBaseNode || stemElongation)) {
            Vector3d h;
            Vector3d n2 = nodes.at(nn - 2);

            if (hasRelCoord()) {
                h = Vector3d(nodes.at(nn - 1));
            } else {
                Vector3d n1 = nodes.at(nn - 1);
                h = n1.minus(n2);
            }
            double olddx = h.length();        // length of last segment
            if (olddx < dx() * (1 - 1e-10)) { // shift node instead of creating a new node
                shiftl = std::min(dx() - olddx, l);
                double sdx = olddx + shiftl; // length of new segment
                h.normalize();
                if (hasRelCoord()) {
                    nodes.at(nn - 1) = h.times(sdx);

                } else {

                    nodes.at(nn - 1) = Vector3d(n2.plus(h.times(sdx)));
                }
                double et = this->calcCreationTime(getLength(true) + shiftl, dt);
                nodeCTs.at(nn - 1) = et; // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
                moved = true;
                l -= shiftl;
                if (l <= 0) { // ==0 should be enough
                    return;
                }
            } else {
                moved = false;
            }
        } else {
            moved = false;
        }
    } // first call or elongation

    // create n+1 new nodes
    double sl = 0; // summed length of created segment
    int n = floor(l / dx());
    for (int i = 0; i < n + 1; i++) {

        double sdx;  // segment length (<=dx)
        if (i < n) { // normal case
            sdx = dx();
        } else { // last segment
            sdx = l - n * dx();
            if (sdx < dxMin() * (1 - 1e-10)) { // quit if l is too small
                if (phytoIdx >= 0) {
                    this->epsilonDx += sdx;
                } else {
                    this->epsilonDx = sdx;
                }
                return;
            }
            this->epsilonDx = 0; // no residual
        }
        sl += sdx;
        Vector3d newnode;
        if (hasRelCoord()) {
            newnode = Vector3d(0., 0., sdx);
            if (stemElongation) {
                for (int ii = 0; ii < addIncrement.size(); ii++) {
                    if (addIncrement.at(ii) >= nn + i) {
                        addIncrement.at(ii) += 1;
                    }
                }
            }
        } else {
            Vector3d newdx = getIncrement(nodes.back(), sdx);
            newnode = Vector3d(nodes.back().plus(newdx));
        }
        addIncrement.push_back(nn + i);
        double et = this->calcCreationTime(getLength(true) + shiftl + sl,
                                           dt); // here length or get length? it s the same because epsilonDx was set back to 0 at beginning of simulate no?
        // in case of impeded growth the node emergence time is not exact anymore,
        // but might break down to temporal resolution
        addNode(newnode, et, size_t(nn + i), stemElongation);
    }
}

/**
 * @brief Creates a new lateral organ at the current organ tip.
 *
 * Iterates over all successor rules; for each rule that applies at the current
 * linking node, draws a lateral type stochastically and constructs the
 * corresponding Root, Stem, or Leaf child.
 *
 * @param dt       Time step received by the parent organ [days].
 * @param verbose  If true, print diagnostic output.
 */
void Organ::createLateral(double dt, bool verbose) {
    auto rp = getOrganRandomParameter(); // rename

    for (int i = 0; i < rp->successorST.size(); i++) { // go through each successor rule

        bool applyHere = getApplyHere(i);

        if (applyHere) {
            int numlats = 1; // how many laterals? default = 1
            if (rp->successorNo.size() > i) {
                numlats = rp->successorNo.at(i);
            }
            for (int nn = 0; nn < numlats; nn++) {

                const Vector3d &pos = Vector3d();
                int p_id = rp->getLateralType(pos, i); // if probabilistic branching

                if (p_id >= 0) {
                    int ot;

                    if ((rp->successorOT.size() > i) && (rp->successorOT.at(i).size() > p_id)) {
                        ot = rp->successorOT.at(i).at(p_id);
                    } else {
                        ot = getOrganRandomParameter()->organType; // default
                    }

                    int st = rp->successorST.at(i).at(p_id);

                    double delay = getLatGrowthDelay(ot, st, dt); // forDelay*multiplyDelay
                    double growth_dt = getLatInitialGrowth(dt);

                    switch (ot) {
                    case Organism::ot_root: {
                        auto lateral = plant.lock()->createRoot(st, delay, shared_from_this(), nodes.size() - 1);
                        lateral->has_rel_coord = this->has_rel_coord;
                        children.push_back(lateral);
                        lateral->simulate(growth_dt, verbose);
                        break;
                    }
                    case Organism::ot_stem: {
                        auto lateral = plant.lock()->createStem(st, delay, shared_from_this(), nodes.size() - 1);
                        lateral->has_rel_coord = this->has_rel_coord;
                        children.push_back(lateral);
                        lateral->simulate(growth_dt, verbose);
                        break;
                    }
                    case Organism::ot_leaf: {
                        auto lateral = plant.lock()->createLeaf(st, delay, shared_from_this(), nodes.size() - 1);
                        lateral->has_rel_coord = this->has_rel_coord;
                        children.push_back(lateral);
                        lateral->simulate(growth_dt, verbose); // age-ageLN,verbose);
                        break;
                    }
                    }
                }
            }
        }
    }
    created_linking_node++;
    storeLinkingNodeLocalId(created_linking_node, verbose); // needed (currently) only for stems when doing nodal growth
}

/**
 * @brief Returns true if successor rule @p i applies at the current linking node.
 *
 * Checks successorWhere[i]: positive values list included nodes, negative values
 * list excluded nodes.  An empty successorWhere entry means the rule applies everywhere.
 *
 * @param i  Successor rule index.
 */
bool Organ::getApplyHere(int i) const {
    bool applyHere;
    auto rp = getOrganRandomParameter(); // rename
    if ((rp->successorWhere.size() > i) && (rp->successorWhere.at(i).size() > 0)) {
        if (!std::signbit(rp->successorWhere.at(i).at(0))) // true if number is signed
        {                                                  // gave which linking nodes to include
            applyHere = (std::find(rp->successorWhere.at(i).begin(), rp->successorWhere.at(i).end(), created_linking_node) != rp->successorWhere.at(i).end());
        } else { // gave which linking nodes to ignore
            applyHere =
                !(std::find(rp->successorWhere.at(i).begin(), rp->successorWhere.at(i).end(), -double(created_linking_node)) != rp->successorWhere.at(i).end());
        }

    } else {
        applyHere = true;
    } // default
    return applyHere;
}

/**
 * @brief Returns the initial growth period to pass to a newly created lateral.
 *
 * Computed as the time elapsed between the lateral node's creation and the
 * current simulation time, so the lateral immediately receives its "missed" growth.
 *
 * @param dt  Time step received by the parent organ [days].
 * @return    Growth period to give to the new lateral [days].
 */
double Organ::getLatInitialGrowth(double dt) {
    double ageLN = this->calcAge(getLength(true)); // MINIMUM age of root when lateral node is created
    ageLN = std::max(ageLN, age - dt);
    return age - ageLN;
}

/**
 * @brief Returns the emergence delay for a new lateral organ.
 *
 * Behaviour depends on the delay-definition mode set in the owning organism:
 * - dd_distance: lateral waits until the parent's apical zone has developed.
 * - dd_time_lat: fixed delay drawn from the parent's ldelay/ldelays parameters.
 * - dd_time_self: fixed delay drawn from the lateral's own ldelay/ldelays parameters.
 *
 * @param ot_lat  Organ type of the lateral to create.
 * @param st_lat  Sub-type of the lateral to create.
 * @param dt      Time step received by the parent organ [days].
 * @return        Emergence delay [days].
 */
double Organ::getLatGrowthDelay(int ot_lat, int st_lat, double dt) const // override for stems
{
    auto rp = getOrganRandomParameter(); // rename
    double growthDelay;                  // store necessary variables to define lateral growth delay
    int delayDefinition = getOrganism()->getDelayDefinition(ot_lat);

    assert(delayDefinition >= 0);

    switch (delayDefinition) {
    case Organism::dd_distance: {                                           // lateral has to wait for a apical zone to develop
        double meanLn = getParameter("lnMean");                             // mean inter-lateral distance
        double effectiveLa = std::max(getParameter("la") - meanLn / 2, 0.); // effective apical distance, observed apical distance is in [la-ln/2, la+ln/2]
        double ageLN = this->calcAge(getLength(true));                      // theoretical age of root when lateral node is created
        ageLN = std::max(ageLN, age - dt);
        double ageLG = this->calcAge(getLength(true) + effectiveLa); // age of the root, when the lateral starts growing (i.e when the apical zone is developed)
        growthDelay = ageLG - ageLN;                                 // time the lateral has to wait
        break;
    }
    case Organism::dd_time_lat: { // lateral has to wait for a fixed time (defined by parent organ)
        growthDelay = std::max(rp->ldelay + plant.lock()->randn() * rp->ldelays, 0.);
        break;
    }
    case Organism::dd_time_self: { // lateral has to wait for a fixed time (defined by lateral organ itself)        
        auto latRp = plant.lock()->getOrganRandomParameter(ot_lat, st_lat); // random parameter of lateral to create
        growthDelay = std::max(latRp->ldelay + plant.lock()->randn() * latRp->ldelays, 0.);
        break;
    }
    default: {
        std::cout << "delayDefinition " << delayDefinition << " " << Organism::dd_distance << " ";
        std::cout << Organism::dd_time_lat << " " << Organism::dd_time_self << std::endl << std::flush;
        std::cout << "				" << (delayDefinition == Organism::dd_distance) << " ";
        std::cout << (delayDefinition == Organism::dd_time_lat) << " " << (delayDefinition == Organism::dd_time_self) << std::endl << std::flush;
        throw std::runtime_error("Delay definition type (delayDefinition) not recognised");
    }
    }
    return growthDelay;
}
/**
 * @brief Checks whether node coordinates are stored in relative form.
 */
bool Organ::hasRelCoord() const { return has_rel_coord; }

} // namespace CPlantBox
