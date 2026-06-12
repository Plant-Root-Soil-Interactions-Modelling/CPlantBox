// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "GrassLeaf.h"

#include "Organism.h"
#include "Plant.h"

#include <cassert>
#include <iostream>
#include <sstream>

namespace CPlantBox {

/**
 * @brief Restoration constructor; tree links and geometry must be set by the caller.
 */
GrassLeaf::GrassLeaf(int id, std::shared_ptr<const OrganSpecificParameter> param, bool alive, bool active, double age, double length, Vector3d partialIHeading_,
                     int pni, bool moved, int oldNON)
    : Organ(id, param, alive, active, age, length, partialIHeading_, pni, moved, oldNON) { }

/**
 * @brief Simulation constructor; sets up the TurtlePolyline anchor at the parent attachment node.
 *
 * The anchor frame is aligned with the parent heading at @p pni, randomly rotated around that
 * heading, then pitched slightly upward so the leaf initially grows along the parent axis.
 */
GrassLeaf::GrassLeaf(std::shared_ptr<Organism> plant, int subtype, double delay, std::shared_ptr<Organ> parent, int pni)
    : Organ(plant, parent, Organism::ot_leaf, subtype, delay, pni) {
    assert(parent != nullptr && "GrassLeaf: parent must not be null");

    // Anchor the meristem at the attachment point with the parent's heading as
    // the initial turtle frame (+x = heading, +z = up).
    Vector3d anchorPos(0., 0., 0.); // relative coordinate; abs is resolved via rel2abs
    Matrix3d anchorFrame;           // identity = upward along x, will be overridden below

    // Anchor frame: align heading column with parent heading at pni.
    if (parent->getNumberOfNodes() > 0) {
        Vector3d h;
        if (parent->hasRelCoord()) {
            parent->rel2abs(); 
            h = parent->heading(pni); // otherwise that does not work
            parent->abs2rel();
        } else {
            h = parent->heading(pni);
        }
        anchorFrame = Matrix3d::ons(h); // heading → col 0; col 1, 2 perpendicular
    }
    auto rotX = anchorFrame.rotX(2 * M_PI * this->plant.lock()->rand());    
    anchorFrame.times(rotX); // random rotation around the heading to make the leaf orientation truly random
    turtle.setAnchor(anchorPos);
    
    // rotate the anchor frame: pitch the heading up by 90° so the leaf grows upward along the parent axis
    Turtle3D t(Vector3d(0., 0., 0.), anchorFrame);
    t.pitchUp(0.03 * M_PI / 2.0); // sheath is almost straight
    anchorFrame = t.getFrame();
    turtle.setAnchorFrame(anchorFrame);

    // Seed nodeIds/nodeCTs for the initial meristem node.
    double creationTime = parent->getNodeCT(pni) + std::max(delay, 0.0);
    nodeIds.push_back(parent->getNodeId(pni)); // anchorframe id and ct is located at index 0
    nodeCTs.push_back(creationTime);
}

/**
 * @brief Returns the absolute Cartesian position of meristem node @p i, with the anchor set to the parent attachment point.
 */
Vector3d GrassLeaf::getNode(int i) const {
    auto anchorPoint = getParent()->getNode(parentNI);
    if (i==0) {
        return anchorPoint; 
    } else {
        turtle.setAnchor(anchorPoint);
        return turtle.getNode(i-1); // -1 because the initial node at index 0 
    }
}

/**
 * @brief Deep-copies this organ into organism @p p; children are not copied (grass leaves have none).
 */
std::shared_ptr<Organ> GrassLeaf::copy(std::shared_ptr<Organism> p) {
    auto gl = std::make_shared<GrassLeaf>(*this); // shallow copy
    gl->plant = p;
    gl->parent = std::weak_ptr<Organ>();
    gl->param_ = std::make_shared<GrassLeafSpecificParameter>(*param());
    // Grass leaf has no children, so no subtree to copy.
    return gl;
}

/**
 * @brief Returns total leaf length (sheath + blade) at @p age using the growth function f_gf.
 *
 * The maximum length k = sheathLength + bladeLength and the initial growth rate
 * r = k / leafGrowthDuration are passed to f_gf.
 */
double GrassLeaf::calcLength(double age) {
    double k = param()->sheathLength + param()->bladeLength;
    double r = k / param()->leafGrowthDuration;
    return getGrassLeafRandomParameter()->f_gf->getLength(age, r, k, shared_from_this());
}

/**
 * @brief Returns the age at which the leaf reaches @p length; inverse of calcLength().
 */
double GrassLeaf::calcAge(double length) const {
    double k = param()->sheathLength + param()->bladeLength;
    double r = k / param()->leafGrowthDuration;
    return getGrassLeafRandomParameter()->f_gf->getAge(length, r, k, shared_from_this());
}


/**
 * @brief Advances the leaf simulation by @p dt days.
 *
 * Each step:
 *  1. Computes the target length from the growth function.
 *  2. Calls growLeaf() to prepend new base segments if the leaf is still growing.
 *  3. Updates sheathLength and bladeLength from the current length ratio.
 *  4. Calls fixPitch() to assign correct pitch values to sheath vs blade segments.
 *  5. Deactivates the organ once full length is reached.
 */
void GrassLeaf::simulate(double dt, bool verbose) {
    firstCall = true; // |\label{l81:incremental}|
    moved = false;
    oldNumberOfNodes = getNumberOfNodes(); // |\label{l81:incremental_end}|

    if (!alive) { // |\label{l81:alive}| 
        return; 
    }

    age += dt;
    double dt_ = (age < dt) ? age : dt; //age < dt emerged in this time step
    if (age <= 0.) {
        return; 
    } // |\label{l81:alive_end}| 
    
    const GrassLeafSpecificParameter &p = *param();
    
    if (active) {
        double age_ = calcAge(length);  // |\label{l81:age}|
        double targetLength = calcLength(age_ + dt_) + this->epsilonDx;
        double scale = getGrassLeafRandomParameter()->f_se->getValue(getNode(0), shared_from_this());
        double dl = std::max(scale * (targetLength - length), 0.); // |\label{l81:age_end}|
        if (dl > 0.) {
            growLeaf(dl, dt_); // |\label{l81:growLeaf}|
        }
        sheathLength = length * p.sheathLength / (p.sheathLength + p.bladeLength); // |\label{l81:fixPitch}|
        bladeLength = length * p.bladeLength / (p.sheathLength + p.bladeLength);
        fixPitch(); // apply blade angle to the new segments // |\label{l81:fixPitch_end}|
    }
    
    active = (length < p.sheathLength + p.bladeLength - 1.e-6) && alive;
}

/**
 * @brief Returns the named scalar parameter; handles leaf-specific names before delegating to Organ::getParameter().
 */
double GrassLeaf::getParameter(std::string name) const {
    if (name == "sheathLength")
        return param()->sheathLength;
    if (name == "leafGrowthDuration")
        return param()->leafGrowthDuration;
    if (name == "bladeAngle")
        return param()->bladeAngle;
    if (name == "bladeWidth")
        return param()->bladeWidth;
    if (name == "bladeLength")
        return param()->bladeLength;
    return Organ::getParameter(name);
}

/**
 * @brief Returns a compact human-readable summary of the leaf state.
 */
std::string GrassLeaf::toString() const {
    std::ostringstream s;
    s << "GrassLeaf id=" << getId() << " age=" << age << " sheath=" << sheathLength << "/" << param()->sheathLength << " blade=" << bladeLength << "/"
      << param()->bladeLength << " nodes=" << turtle.size() << "\n";
    return s.str();
}

/**
 * @brief Returns the stochastic parameter set for this leaf's sub-type.
 */
std::shared_ptr<GrassLeafRandomParameter> GrassLeaf::getGrassLeafRandomParameter() const {
    return std::dynamic_pointer_cast<GrassLeafRandomParameter>(getOrganism()->getOrganRandomParameter(Organism::ot_leaf, param()->subType));
}

/**
 * @brief Returns the realised (instance-level) parameters cast to GrassLeafSpecificParameter.
 */
std::shared_ptr<const GrassLeafSpecificParameter> GrassLeaf::param() const { return std::dynamic_pointer_cast<const GrassLeafSpecificParameter>(param_); }

/**
 * @brief Prepends new segments totalling @p dl [cm] at the base of the TurtlePolyline.
 *
 * Implements intercalary meristem growth: new segments are inserted before the existing
 * base so that older tissue moves toward the tip.  Each segment gets a creation time
 * interpolated from the organ age and the base creation time.  Pitch angles are left at
 * zero here and corrected afterwards by fixPitch().
 */
void GrassLeaf::growLeaf(double dl, double dt) {
    const GrassLeafSpecificParameter &p = *param(); 

    double baseCT = nodeCTs[turtle.getInitialNodeIndex()]; // organ creation time 

    int n = static_cast<int>(std::floor(dl / dx()));
    double last = dl - n * dx(); // if last piece is below dxMin, absorb it into epsilonDx and skip
    bool addLast = (last >= dxMin() * (1. - 1e-10));
    if (!addLast) {
        epsilonDx += last;
    } else {
        epsilonDx = 0.;
    }
    int nSegs = n + (addLast ? 1 : 0);

    // prepend in reverse order so that node 0 ends up as the new base
    for (int i = nSegs - 1; i >= 0; --i) {
        double sdx = (i < n) ? dx() : last;
        length += sdx;
        double age_ = calcAge(length);
        double a = std::min(std::max(age_, age - dt), age);  // in case of failing we fall back to temporal discretization   
        double ct = a + baseCT;
        int nid = plant.lock()->getNodeIndex();
        turtle.addNodeFront(sdx, 0., 0., 0.);
        nodeIds.insert(nodeIds.begin()+1, nid); // +1 because the initial node at index 0 is the anchor and does not have a nodeId
        nodeCTs.insert(nodeCTs.begin()+1, ct);
        moved = true;
    }
}

/**
 * @brief Assigns pitch angles to every turtle node based on its position along the leaf.
 *
 * Sheath nodes (cumulative arc-length < sheathLength) get pitch = 0 (straight).
 * Blade nodes get pitch = bladeBending * segment dist [rad/cm], producing a gentle outward curve.
 * The sheath/blade boundary is located with TurtlePolyline::getNodeIndexAtLength().
 *
 * The node deque is mutated in place via const_cast (safe because @c turtle is mutable).
 */
void GrassLeaf::fixPitch() {
    if (turtle.size() == 0)
        return;
    int bladeStart = turtle.getNodeIndexAtLength(sheathLength);
    auto &nodes = const_cast<std::deque<TurtlePolyline::TurtleNode> &>(turtle.getTurtleNodes());
    for (int i = 0; i < nodes.size(); i++) {
        nodes[i].pitch = (i < bladeStart) ? 0. : param()->bladeBending * nodes[i].dist;    
    }
    double p = length/(param()->sheathLength + param()->bladeLength);
    nodes[bladeStart].pitch = param()->bladeAngle * p; 
}


/**
 * @brief Returns two 3D coordinates at node @p i for VTK polygon rendering.
 *
 * @param i  node index (0-based)
 * @return   empty for sheath nodes, or {leftEdge, rightEdge} for blade nodes
 */
std::vector<Vector3d> GrassLeaf::getLeafVis(int i) {

    int turtleIdx = std::max(0, i - 1); // Cumulative arc length to node i (turtle nodes are offset by 1 vs organ nodes)
    double cumLen = turtle.getLength(turtleIdx); 
    double p = (cumLen - getSheathLength()) / getBladeLength(); 
    double a = getParent()->param()->a; // parent radius
    Vector3d node = getNode(i);
    Vector3d y1 = turtle.getNodeFrame(turtleIdx).column(1);

    if (p<=0)  {
        return { node.plus(y1.times(a)), node.minus(y1.times(a)) };
    }

    double bw; 
    double halfWidth = param()->bladeWidth / 2.;
    double halfWidth_ = std::min(length/2., halfWidth); // avoid unrealistic width at the very base of the blade  

    if (p<0.25) {
        bw = (halfWidth_-a) * (p/0.25)+a;
    } else if (p<0.6) {
        bw = halfWidth_;
    } else {
        bw = halfWidth_*(1. - (p-0.6)/0.4);
    }
        
    return { node.plus(y1.times(bw)), node.minus(y1.times(bw)) };
}

} // end namespace CPlantBox
