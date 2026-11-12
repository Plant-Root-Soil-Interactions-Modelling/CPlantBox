// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// own header
#include "Organ.h"

// project headers
#include "Organism.h"
#include "Plant.h"
#include "organparameter.h"

// standard library
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <typeinfo>

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
Organ::Organ(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length, Vector3d partialIHeading_, int pni, bool moved, int oldNON)
    : parentNI(pni), partialIHeading(partialIHeading_), plant(), parent(), id(id), param_(param), alive(alive), active(active), age(age), length(length), moved(moved), oldNumberOfNodes(oldNON) {}

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
 * @brief Returns all member variables of the organ by name.
 *
 * Links to organs are stored as ids ("parentId" is -1 without parent, "childIds" lists the children).
 * The owning organism is not stored. Derived classes extend the returned map.
 * Only valid while the plant is in absolute coordinates (between simulate() calls); relative-coordinate scratch
 * variables (has_rel_coord, addIncrement) and firstCall (reset in simulate()) are therefore not stored.
 */
OrganState Organ::getState() const {
    OrganState s;
    stateSet(s, "id", id);
    stateSet(s, "organType", organType());
    stateSet(s, "subType", param_ ? param_->subType : -1);
    if (param_) {
        stateSet(s, "param.a", param_->a);
    }
    auto p = parent.lock();
    stateSet(s, "parentId", p ? p->getId() : -1);
    std::vector<int> cids;
    for (const auto &c : children)
        cids.push_back(c->getId());
    stateSet(s, "childIds", cids);
    stateSet(s, "parentNI", parentNI);
    stateSet(s, "partialIHeading", partialIHeading);
    stateSet(s, "alive", alive);
    stateSet(s, "active", active);
    stateSet(s, "age", age);
    stateSet(s, "length", length);
    stateSet(s, "epsilonDx", epsilonDx);
    stateSet(s, "created_linking_node", static_cast<double>(created_linking_node));
    stateSet(s, "nodes", nodes);
    stateSet(s, "nodeIds", nodeIds);
    stateSet(s, "nodeCTs", nodeCTs);
    stateSet(s, "moved", moved);
    stateSet(s, "oldNumberOfNodes", oldNumberOfNodes);
    return s;
}

/**
 * @brief Restores member variables from a state created by getState().
 *
 * Keys missing in @p state are left unchanged. The organ id is const and not restored.
 * If @p organs (id -> organ) is non-empty, parent and children are looked up by id; links with unknown ids are left unchanged.
 */
void Organ::setState(const OrganState &state, const std::map<int, std::shared_ptr<Organ>> &organs) {
    if (param_ && typeid(*param_) == typeid(OrganSpecificParameter)) { // derived classes restore their own parameter type
        auto p = std::make_shared<OrganSpecificParameter>(*param_);
        stateGet(state, "subType", p->subType);
        stateGet(state, "param.a", p->a);
        param_ = p;
    }
    stateGet(state, "parentNI", parentNI);
    stateGet(state, "partialIHeading", partialIHeading);
    stateGet(state, "alive", alive);
    stateGet(state, "active", active);
    stateGet(state, "age", age);
    stateGet(state, "length", length);
    stateGet(state, "epsilonDx", epsilonDx);
    stateGet(state, "created_linking_node", created_linking_node);
    stateGet(state, "nodes", nodes);
    stateGet(state, "nodeIds", nodeIds);
    stateGet(state, "nodeCTs", nodeCTs);
    stateGet(state, "moved", moved);
    stateGet(state, "oldNumberOfNodes", oldNumberOfNodes);
    if (!organs.empty()) {
        int pid = -1;
        if (stateGet(state, "parentId", pid)) {
            auto it = organs.find(pid);
            if (it != organs.end()) {
                parent = it->second;
            } else if (pid < 0) {
                parent.reset();
            }
        }
        std::vector<int> cids;
        if (stateGet(state, "childIds", cids)) {
            std::vector<std::shared_ptr<Organ>> c;
            for (int cid : cids) {
                auto it = organs.find(cid);
                if (it == organs.end()) {
                    break;
                }
                c.push_back(it->second);
            }
            if (c.size() == cids.size()) { // only relink if all children are known
                children = c;
            }
        }
    }
}

/**
 * @brief Returns length from the first node up to node index @p i.
 *
 * @param i  Local node index (exclusive upper bound).
 * @return   Length [cm].
 */
double Organ::getLength(int i) const {
    double l = 0.; // length until node i
    for (int j = 0; j < i; j++) {
        l += nodes.at(j + 1).length(); // relative node = segment vector in the local frame, its length equals the segment length
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
 * @brief Returns the number of polyline segments.
 *
 * Segments connect consecutive nodes, so the count equals the number of
 * nodes minus one (clamped to zero for empty organs).
 *
 * @return Number of segments.
 */
int Organ::getNumberOfSegments() const {
    const int numberOfNodes = getNumberOfNodes();
    return (numberOfNodes > 0) ? (numberOfNodes - 1) : 0;
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
std::shared_ptr<OrganRandomParameter> Organ::getOrganRandomParameter() const { return plant.lock()->getOrganRandomParameter(this->organType(), param_->subType); }

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
 * @brief Appends a node with an explicit global id at the end of the organ.
 *
 * Nodes cannot be deleted; organs become deactivated or die instead.
 *
 * @param n      Node coordinate.
 * @param id     Global node index.
 * @param t      Creation time [days].
 */
void Organ::addNode(Vector3d n, int id, double t) {
    nodes.push_back(n);    // node
    nodeIds.push_back(id); // unique id
    nodeCTs.push_back(t);  // exact creation time
}

/**
 * @brief Appends a node using the next global index from the organism.
 *
 * Convenience wrapper around addNode(n, id, t) that fetches
 * the next available global node index from the owning organism.
 *
 * @param n  Node coordinate.
 * @param t  Creation time [days].
 */
void Organ::addNode(Vector3d n, double t) { addNode(n, plant.lock()->getNodeIndex(), t); }

/**
 * @brief Inserts a node between existing nodes (internodal growth).
 *
 * Coordinates and creation time are inserted at @p index; the global id is appended to the list of node ids,
 * so that existing nodes keep their global index (simplifies the update of the node positions in a MappedPlant).
 * Children attached at or after the new node get their parent node index incremented.
 * Override in Stem to also update its linking nodes.
 *
 * @param n      Node coordinate.
 * @param id     Global node index.
 * @param t      Creation time [days].
 * @param index  Position at which to insert the node.
 */
void Organ::insertNode(Vector3d n, double t, size_t index) {
    nodes.insert(nodes.begin() + index, n);
    nodeIds.push_back(plant.lock()->getNodeIndex());
    nodeCTs.insert(nodeCTs.begin() + index, t);
    for (auto kid : children) { // if carries children after the added node, update their "parent node index"
        if ((kid->getParentNI() >= index) && (kid->getParentNI() > 0)) {
            kid->incParentNI();
            kid->nodeIds.at(0) = kid->getParent()->getNodeId(kid->getParentNI());
        }
    }
}

/**
 * @brief Returns all polyline segments as pairs of global node indices.
 *
 * The polyline {n0, n1, n2, n3} produces segments {[i0,i1], [i1,i2], [i2,i3]}.
 * Returns an empty vector if fewer than two nodes exist.
 */
std::vector<Vector2i> Organ::getSegments() const {
    int n = this->getNumberOfNodes();
    if (n > 1) {
        std::vector<Vector2i> segs = std::vector<Vector2i>(n - 1);
        for (size_t i = 0; i < n - 1; i++) {
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
    bool forCarbon_limitedGrowth = (all && (this->getAge() > 0)); // when ask for "all" organs which have age > 0 even if getNumberOfNodes() == 1
    bool notSeed = (this->organType() != Organism::ot_seed);

    if ((this->getNumberOfNodes() > 1 || forCarbon_limitedGrowth) && notSeed) { //&& notBulb
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
 * @brief Returns a named scalar parameter of the organ.
 *
 * Recognises organ-specific fields (subType, a, radius, diameter, iHeadingX/Y/Z, parentNI,
 * parent-node), member variables (organType, id, alive, active, age, length, lengthTh,
 * numberOfNodes, numberOfSegments, hasMoved, oldNumberOfNodes, creationTime,
 * order, one), and delegates to OrganRandomParameter::getParameter() for anything else.
 *
 * Intended for post-processing; flexible but slower than direct member access.
 * Override in derived classes to expose additional parameters.
 */
double Organ::getParameter(std::string name, std::map<std::string, double> addParams) const {

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
        return getInitialHeading().x;
    }
    if (name == "iHeadingY") { // root initial heading y - coordinate [cm]
        return getInitialHeading().y;
    }
    if (name == "iHeadingZ") { // root initial heading z - coordinate [cm]
        return getInitialHeading().z;
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
        auto absn = getNodes();
        for (int i = o; i < getNumberOfNodes(); i += nn) {
            auto n = absn.at(i);
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
    str << Organism::organTypeNames.at(this->organType()) << " #" << getId() << ": sub type " << param_->subType << ", realized length " << getLength(true) << " cm , theoretic length "
        << getLength(false) << " cm , age " << getAge() << " days, alive " << isAlive() << ", active " << isActive() << ", number of nodes " << this->getNumberOfNodes() << ", with " << children.size()
        << " children";
    return str.str();
}

/**
 * @brief Returns the initial absolute heading
 *
 * Reconstructed from partialIHeading and the parent's heading at parentNI.
 * Handles base organs (no parent), shoot-borne roots, and root-borne shoots as special cases.
 */
Vector3d Organ::getInitialHeading() const {
    Vector3d parentHeading; // heading of the parent at the attachment point
    auto parent = getParent();
    bool isRoot = parent ? (organType() == Organism::ot_root) : (organType() <= Organism::ot_root); // without parent, unrecognized organs (0) are treated as roots
    bool isBaseOrgan = (parent == nullptr) || (parent->organType() == Organism::ot_seed); // in case of class RootSystem base roots (tap, basal, shootborne) or Organism organs created wwithout parent
    bool isShootBornRoot = (parent != nullptr) && (parent->organType() == Organism::ot_stem) && (organType() == Organism::ot_root);
    bool isRootBornShoot = (parent != nullptr) && (parent->organType() == Organism::ot_root) && ((organType() == Organism::ot_stem) || (organType() == Organism::ot_leaf));
    if (isBaseOrgan || isShootBornRoot || isRootBornShoot) {
        parentHeading = isRoot ? Vector3d(0, 0, -1) : Vector3d(0, 0, 1); // roots grow downwards, shoots upwards
    } else {
        parentHeading = parent->heading(parentNI);
    }
    Vector3d new_heading = Matrix3d::ons(parentHeading).times(this->partialIHeading);
    return Matrix3d::ons(new_heading).column(0);
}

/**
 * @brief Returns the absolute coordinate of node @p i.
 *
 * Organs with relative coordinates (Stem, Leaf, and roots carried by them) calculate it from getOrigin() and the relative nodes:
 * node j is stored as the segment vector expressed in the orthonormal system of the direction of the segment [j-2, j-1]
 * (for j = 1 of the initial heading), see segmentHeading().
 */
Vector3d Organ::getNode(int i) const {
    Vector3d p = getParent() ? getOrigin() : nodes.at(0);
    Vector3d h = getInitialHeading(); // direction of the previous segment
    for (int j = 1; j <= i; j++) {
        Vector3d s = Matrix3d::ons(h).times(nodes.at(j)); // segment [j-1, j]
        p = p.plus(s);
        h = s.normalized();
    }
    return p;
}

/**
 * @brief Returns the absolute coordinates of all nodes.
 */
std::vector<Vector3d> Organ::getNodes() const {
    std::vector<Vector3d> v;
    v.reserve(getNumberOfNodes());
    if (nodes.empty()) {
        return v;
    }
    v.push_back(getParent() ? getOrigin() : nodes.at(0));
    Vector3d h = getInitialHeading(); // direction of the previous segment
    for (size_t j = 1; j < nodes.size(); j++) {
        Vector3d s = Matrix3d::ons(h).times(nodes.at(j)); // segment [j-1, j]
        v.push_back(v.back().plus(s));
        h = s.normalized();
    }
    return v;
}

/**
 * @brief Returns the absolute heading at node @p n.
 *
 * The heading at node 0 is the initial heading, the heading at node j is the direction of the segment [j-1, j].
 * Negative @p n is treated as the last node. Relative-coordinate organs use their segment coordinates;
 * absolute-coordinate organs use the direction between their absolute node positions.
 */
Vector3d Organ::heading(int n) const {
    if (n < 0) {
        n = nodes.size() - 1;
    }
    n = std::min(n, int(nodes.size()) - 1);
    Vector3d h = getInitialHeading();
    for (int j = 1; j <= n; j++) {
        h = Matrix3d::ons(h).times(nodes[j]);
        h.normalize();
    }
    return h;
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
 * @brief Returns the organ length [cm] corresponding to a volume, assuming a cylinder.
 *
 * @param volume_ Organ volume [cm³].
 * @return Organ length [cm].
 */
double Organ::orgVolume2Length(double volume_) { return volume_ / (M_PI * param_->a * param_->a); }

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
                        children.push_back(lateral);
                        lateral->simulate(growth_dt, verbose);
                        break;
                    }
                    case Organism::ot_stem: {
                        auto lateral = plant.lock()->createStem(st, delay, shared_from_this(), nodes.size() - 1);
                        children.push_back(lateral);
                        lateral->simulate(growth_dt, verbose);
                        break;
                    }
                    case Organism::ot_leaf: {
                        auto lateral = plant.lock()->createLeaf(st, delay, shared_from_this(), nodes.size() - 1);
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
    case Organism::dd_time_self: {                                          // lateral has to wait for a fixed time (defined by lateral organ itself)
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
            applyHere = !(std::find(rp->successorWhere.at(i).begin(), rp->successorWhere.at(i).end(), -double(created_linking_node)) != rp->successorWhere.at(i).end());
        }

    } else {
        applyHere = true;
    } // default
    return applyHere;
}

/* state helpers */
void Organ::stateSet(OrganState &s, const std::string &key, const std::vector<Vector3d> &v) {
    std::vector<double> d;
    d.reserve(3 * v.size());
    for (const auto &p : v) {
        d.insert(d.end(), {p.x, p.y, p.z});
    }
    s[key] = d;
}
bool Organ::stateGet(const OrganState &s, const std::string &key, double &v) {
    auto it = s.find(key);
    if (it == s.end() || it->second.empty())
        return false;
    v = it->second[0];
    return true;
}
bool Organ::stateGet(const OrganState &s, const std::string &key, bool &v) {
    double d;
    if (!stateGet(s, key, d))
        return false;
    v = (d != 0.);
    return true;
}
bool Organ::stateGet(const OrganState &s, const std::string &key, int &v) {
    double d;
    if (!stateGet(s, key, d))
        return false;
    v = static_cast<int>(std::lround(d));
    return true;
}
bool Organ::stateGet(const OrganState &s, const std::string &key, size_t &v) {
    double d;
    if (!stateGet(s, key, d))
        return false;
    v = static_cast<size_t>(std::llround(d));
    return true;
}
bool Organ::stateGet(const OrganState &s, const std::string &key, std::vector<double> &v) {
    auto it = s.find(key);
    if (it == s.end())
        return false;
    v = it->second;
    return true;
}
bool Organ::stateGet(const OrganState &s, const std::string &key, std::vector<int> &v) {
    auto it = s.find(key);
    if (it == s.end())
        return false;
    v.clear();
    for (double d : it->second)
        v.push_back(static_cast<int>(std::lround(d)));
    return true;
}
bool Organ::stateGet(const OrganState &s, const std::string &key, Vector3d &v) {
    auto it = s.find(key);
    if (it == s.end() || it->second.size() < 3)
        return false;
    v = Vector3d(it->second[0], it->second[1], it->second[2]);
    return true;
}
bool Organ::stateGet(const OrganState &s, const std::string &key, std::vector<Vector3d> &v) {
    auto it = s.find(key);
    if (it == s.end())
        return false;
    v.clear();
    for (size_t i = 0; i + 2 < it->second.size(); i += 3)
        v.push_back(Vector3d(it->second[i], it->second[i + 1], it->second[i + 2]));
    return true;
}

} // namespace CPlantBox
