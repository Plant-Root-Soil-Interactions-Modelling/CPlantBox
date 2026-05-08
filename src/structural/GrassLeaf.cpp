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
        Vector3d h = parent->heading(pni);
        anchorFrame = Matrix3d::ons(h); // heading → col 0; col 1, 2 perpendicular
    }

    meristem.setAnchor(anchorPos);
    meristem.setAnchorFrame(anchorFrame);

    // Seed nodeIds/nodeCTs for the initial meristem node.
    double creationTime = parent->getNodeCT(pni) + std::max(delay, 0.0);
    nodeIds.push_back(plant->getNodeIndex());
    nodeCTs.push_back(creationTime);
}

/**
 * @brief Returns the absolute Cartesian position of meristem node @p i, with the anchor set to the parent attachment point.
 */
Vector3d GrassLeaf::getNode(int i) const {
    meristem.setAnchor(getParent()->getNode(parentNI));
    return meristem.getNode(i);
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
 * @brief Advances the leaf by @p dt days through sheath elongation, blade delay, and blade elongation phases.
 */
void GrassLeaf::simulate(double dt, bool verbose) {

    oldNumberOfNodes = meristem.size();

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

        double age_ = calcAge(sheathLengthGrown, p.sheathLength, p.sheathLength/p.sheathDuration); 
        double targetSheathLength = calcLength(age_ + dt_, p.sheathLength, p.sheathLength/p.sheathDuration) + this->epsilonDx;
        double dl = targetSheathLength - sheathLengthGrown; 
        if (dl > 0.) {
            growSheath(dl, dt);
        }

        // // blade elongation
        // if (isSheathComplete() && isBladeEmerged()) {
        //     double bladeRate = (p.bladeDuration > 0.)
        //                     ? p.bladeLength / p.bladeDuration
        //                     : p.bladeLength;
        //     double dl = std::min(bladeRate * dt, p.bladeLength - bladeLengthGrown);
        //     if (dl > 0.) growBlade(dl, dt);
        // }
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
    if (name == "sheathDuration")
        return param()->sheathDuration;
    if (name == "bladeAngle")
        return param()->bladeAngle;
    if (name == "bladeWidth")
        return param()->bladeWidth;
    if (name == "bladeLength")
        return param()->bladeLength;
    if (name == "bladeDelay")
        return param()->bladeDelay;
    if (name == "bladeDuration")
        return param()->bladeDuration;
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
 * @brief Returns the organ age at which sheath elongation completes.
 */
double GrassLeaf::calcAgeAtSheathComplete() const {
    const GrassLeafSpecificParameter &p = *param();
    double sheathRate = (p.sheathDuration > 0.) ? p.sheathLength / p.sheathDuration : 1e9;
    return p.sheathLength / sheathRate;
}

/**
 * @brief Grows the sheath by @p dl [cm]; subdivides into dx()-sized segments via addMeristemNodes().
 */
void GrassLeaf::growSheath(double dl, double dt) {
    sheathLengthGrown += dl;
    addSheathNodes(dl, 0., dt);
}

/**
 * @brief Grows the blade by @p dl [cm]; applies bladeAngle yaw on the first segment, then subdivides via addMeristemNodes().
 */
void GrassLeaf::growBlade(double dl, double dt) {
    const GrassLeafSpecificParameter &p = *param();
    double yaw = (bladeLengthGrown == 0.) ? p.bladeAngle : 0.; // blade angle only on first segment
    bladeLengthGrown += dl;
    addSheathNodes(dl, yaw, dt);
}

/**
 * @brief Subdivides @p dl into dx()-sized segments and prepends each as a meristem node.
 *
 * Mirrors the logic of Organ::createSegments:
 *  - Full segments of length dx() are prepended first.
 *  - The remaining piece is prepended only if it is >= dxMin(); otherwise it
 *    is accumulated in epsilonDx and prepended in a future time step.
 *  - @p yaw is applied only to the innermost (first-prepended) segment; all
 *    subsequent segments are straight.
 *  - nodeId and creation time are recorded for each node, compatible with
 *    Organ::nodeIds / Organ::nodeCTs conventions.
 *
 * @param dl   Total length to add [cm].
 * @param yaw  Yaw rotation [rad] for the innermost new segment (0 = straight).
 * @param dt   Current time step [days]; used for creation-time interpolation.
 */
void GrassLeaf::addSheathNodes(double dl, double yaw, double dt) {
    if (dl <= 0.)
        return;

    double baseCT = nodeCTs.at(meristem.getInitialNodeIndex());
    double totalLen = sheathLengthGrown; // already updated by caller
    double accLen = totalLen - dl;                          // length before this call

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
        accLen += sdx;
        double age_ = calcAge(std::max(accLen, 0.));
        double a = std::min(std::max(age_, age - dt), age);
        double ct = a + baseCT;
        int nid = plant.lock()->getNodeIndex();
        double segYaw = (i == nSegs - 1) ? yaw : 0.; // yaw only on innermost segment
        meristem.addNodeFront(sdx, segYaw, 0., 0.);
        nodeIds.insert(nodeIds.begin(), nid);
        nodeCTs.insert(nodeCTs.begin(), ct);
    }
}




double GrassLeaf::calcLength(double length, double k, double r) {
    assert(age >= 0 && "GrassLeaf::calcLength() negative leaf age");
    return getGrassLeafRandomParameter()->f_gf->getLength(age, r, k, shared_from_this());
}


double GrassLeaf::calcAge(double length, double k, double r) {
    assert(length >= 0 && "GrassLeaf::calcAge() negative leaf length");
    return getGrassLeafRandomParameter()->f_gf->getAge(length, r, k, shared_from_this());
}



/**
 * @brief Returns total leaf length (sheath + blade) at @p age_; piecewise-linear approximation.
 */
double GrassLeaf::calcLength(double age_) {
    throw std::runtime_error("GrassLeaf::calcLength(age) not implemented, use calcLength(age, k, r) instead");
    const GrassLeafSpecificParameter &p = *param();
    if (age_ <= 0.)
        return 0.;

    double sheathRate = (p.sheathDuration > 0.) ? p.sheathLength / p.sheathDuration : p.sheathLength;
    double sheathLen = std::min(sheathRate * age_, p.sheathLength);
    double remaining = age_ - sheathLen / std::max(sheathRate, 1e-12);
    if (remaining <= 0.)
        return sheathLen;

    double afterDelay = remaining - p.bladeDelay;
    if (afterDelay <= 0.)
        return sheathLen; // still in blade delay

    double bladeRate = (p.bladeDuration > 0.) ? p.bladeLength / p.bladeDuration : p.bladeLength;
    double bladeLen = std::min(bladeRate * afterDelay, p.bladeLength);
    return sheathLen + bladeLen;
}

/**
 * @brief Returns the age at which the leaf reaches @p length_; inverse of calcLength().
 */
double GrassLeaf::calcAge(double length_) const {
    throw std::runtime_error("GrassLeaf::calcAge(length) not implemented, use calcAge(length, k, r) instead");
    const GrassLeafSpecificParameter &p = *param();
    if (length_ <= 0.)
        return 0.;

    double sheathRate = (p.sheathDuration > 0.) ? p.sheathLength / p.sheathDuration : 1e9;
    if (length_ <= p.sheathLength) {
        return length_ / sheathRate;
    }

    double ageAtSheath = p.sheathLength / sheathRate;
    double bladeRate = (p.bladeDuration > 0.) ? p.bladeLength / p.bladeDuration : 1e9;
    double bladeLen = length_ - p.sheathLength;
    return ageAtSheath + p.bladeDelay + bladeLen / bladeRate;
}




} // end namespace CPlantBox
