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
    : Organ(id, param, alive, active, age, length, partialIHeading_, pni, moved, oldNON) {}

/**
 * @brief Simulation constructor; realizes parameters from GrassLeafRandomParameter and initialises the Meristem anchor at the parent attachment point.
 */
GrassLeaf::GrassLeaf(std::shared_ptr<Organism> plant, int subtype, double delay, std::shared_ptr<Organ> parent, int pni)
    : Organ(plant, parent, Organism::ot_leaf, subtype, delay, pni) {
    assert(parent != nullptr && "GrassLeaf: parent must not be null");

    std::cout << "Creating GrassLeaf with parent " << parent->getId() << " pni " << pni << " delay " << delay << "\n" << std::flush;

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
    meristem.setAnchor(anchorPos);
    
    // rotate the anchor frame: pitch the heading up by 90° so the leaf grows upward along the parent axis
    Turtle3D t(Vector3d(0., 0., 0.), anchorFrame);
    t.pitchUp(0.1 * M_PI / 2.0);
    anchorFrame = t.getFrame();
    meristem.setAnchorFrame(anchorFrame);

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
        meristem.setAnchor(anchorPoint);
        return meristem.getNode(i-1); // -1 because the initial node at index 0 
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
 * @brief Advances the leaf by @p dt days 
 */
void GrassLeaf::simulate(double dt, bool verbose) {

    oldNumberOfNodes = getNumberOfNodes();

    if (!alive) {
        return; // leafs don't die (yet)
    }

    // std::cout << "Simulating GrassLeaf id=" << getId() << " age=" << age << " dt=" << dt << "\n" << std::flush;

    const GrassLeafSpecificParameter &p = *param();

    age += dt;
    double dt_ = (age < dt) ? age : dt; // time step; age < dt means the organ emerged in this time step

    if (age <= 0.) {
        return; // still in delay
    }

    // GrassLeaf hold no children

    if (active) {

        // double age_ = calcAge(sheathLengthGrown, p.sheathLength, p.sheathLength/p.leafGrowthDuration); 
        // double targetSheathLength = calcLength(age_ + dt_, p.sheathLength, p.sheathLength/p.leafGrowthDuration) + this->epsilonDx;        
        // double dl = targetSheathLength - sheathLengthGrown; 
        // if (dl > 0.) {
        //     this->epsilonDx = 0; // its in dl already
        //     growSheath(dl, dt);
        // }

        double age_ = calcAge(bladeLengthGrown); 
        double targetBladeLength = calcLength(age_ + dt_) + this->epsilonDx;
        double dl = targetBladeLength - bladeLengthGrown;
        if (dl > 0.) {
            growLeaf(dl, dt);
        }
        
    }

    // Sync the length field used by the base class infrastructure.
    length = sheathLengthGrown + bladeLengthGrown;
    active = (length < p.sheathLength + p.bladeLength) && alive;
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
    s << "GrassLeaf id=" << getId() << " age=" << age << " sheath=" << sheathLengthGrown << "/" << param()->sheathLength << " blade=" << bladeLengthGrown << "/"
      << param()->bladeLength << " nodes=" << meristem.size() << "\n";
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
 * @brief Grows the blade by @p dl [cm]; applies bladeAngle yaw on the first segment, then subdivides via addMeristemNodes().
 */
void GrassLeaf::growLeaf(double dl, double dt) {

    const GrassLeafSpecificParameter &p = *param(); 

    double baseCT = nodeCTs[meristem.getInitialNodeIndex()]; // organ creation time 

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
        bladeLengthGrown += sdx;
        double age_ = calcAge(bladeLengthGrown);
        double a = std::min(std::max(age_, age - dt), age);        
        double ct = a + baseCT;
        int nid = plant.lock()->getNodeIndex();
        meristem.addNodeFront(sdx, 0., 0.03, 0.);
        nodeIds.insert(nodeIds.begin()+1, nid); // +1 because the initial node at index 0 is the anchor and does not have a nodeId
        nodeCTs.insert(nodeCTs.begin()+1, ct);
    }
}

/**
 *  
 */
void GrassLeaf::elongateLeaf(double dl, double dt) {


}


/**
 * @brief Returns total leaf length (sheath + blade) at @p age_; piecewise-linear approximation.
 */
double GrassLeaf::calcLength(double age) {
    double k = param()->sheathLength + param()->bladeLength;
    double r = k / param()->leafGrowthDuration;
    return getGrassLeafRandomParameter()->f_gf->getLength(age, r, k, shared_from_this());
}

/**
 * @brief Returns the age at which the leaf reaches @p length_; inverse of calcLength().
 */
double GrassLeaf::calcAge(double length) const {
    double k = param()->sheathLength + param()->bladeLength;
    double r = k / param()->leafGrowthDuration;
    return getGrassLeafRandomParameter()->f_gf->getAge(length, r, k, shared_from_this());
}

/**
 * @brief Returns two 3D coordinates at node @p i for VTK polygon rendering.
 *
 * @param i  node index (0-based)
 * @return   empty for sheath nodes, or {leftEdge, rightEdge} for blade nodes
 */
std::vector<Vector3d> GrassLeaf::getLeafVis(int i) {
    // Cumulative arc length to node i
    double cumLen = 0.;
    for (int j = 0; j < i; j++) {
        cumLen += getNode(j + 1).minus(getNode(j)).length();
    } // TODO meristem knows length, do shortcut
    double p = cumLen/(getSheathLengthGrown() + getBladeLengthGrown());

    // Local segment direction at node i
    Vector3d dir;
    if (i > 0) {
        dir = getNode(i).minus(getNode(i - 1));
    } else if (getNumberOfNodes() > 1) {
        dir = getNode(1).minus(getNode(0));
    } else {
        return {};
    }
    if (dir.length() < 1e-12) {
        return {};
    }
    dir.normalize();

    // Lateral direction: perpendicular to segment and world-up (0,0,1)
    Vector3d up(0., 0., 1.);
    Vector3d y1 = dir.cross(up);
    if (y1.length() < 1e-6) {
        y1 = Vector3d(1., 0., 0.);  // fallback when segment is nearly vertical
    } else {
        y1.normalize();
    }

    double bw; 
    double halfWidth = param()->bladeWidth / 2.;
    double a = getParent()->param()->a; // parent radius
    if (p<0.25) {
        bw = (halfWidth-a) * (p/0.25)+a;
    } else if (p<0.6) {
        bw = halfWidth;
    } else {
        bw = halfWidth*(1. - (p-0.6)/0.4);
    }
    
    Vector3d node = getNode(i);
    return { node.plus(y1.times(bw)), node.minus(y1.times(bw)) };
}

} // end namespace CPlantBox
