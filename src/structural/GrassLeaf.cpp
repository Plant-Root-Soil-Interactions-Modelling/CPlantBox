// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "GrassLeaf.h"

#include "Organism.h"
#include "Plant.h"

#include <cassert>
#include <iostream>
#include <sstream>

namespace CPlantBox {

// ========================================================================== //
//  Constructors
// ========================================================================== //

/**
 * Restoration constructor.
 * Tree links and geometry must be set by the caller after construction.
 */
GrassLeaf::GrassLeaf(int id,
                     std::shared_ptr<const OrganSpecificParameter> param,
                     bool alive, bool active,
                     double age, double length,
                     Vector3d partialIHeading_,
                     int pni,
                     bool moved,
                     int oldNON)
    : Organ(id, param, alive, active, age, length, partialIHeading_, pni, moved, oldNON)
{}

/**
 * Simulation constructor.
 * Called by the plant/organism during growth.  Parameters are realised from
 * the GrassLeafRandomParameter registered for @p subtype.
 *
 * The Meristem anchor is placed at the parent attachment point and given the
 * initial heading of the organ (upward along the stem, as is typical for a
 * grass sheath).
 */
GrassLeaf::GrassLeaf(std::shared_ptr<Organism> plant,
                     int subtype,
                     double delay,
                     std::shared_ptr<Organ> parent,
                     int pni)
    : Organ(plant, parent, Organism::ot_leaf, subtype, delay, pni)
{
    assert(parent != nullptr && "GrassLeaf: parent must not be null");

    std::cout << "Creating GrassLeaf with parent " << parent->getId() << " pni " << pni << " delay " << delay << "\n" << std::flush;

    // Anchor the meristem at the attachment point with the parent's heading as
    // the initial turtle frame (+x = heading, +z = up).
    Vector3d anchorPos(0., 0., 0.); // relative coordinate; abs is resolved via rel2abs
    Matrix3d anchorFrame;           // identity = upward along x, will be overridden below

    // Build a sensible initial frame from the parent heading at pni.
    if (parent->getNumberOfNodes() > 0) {
        Vector3d h = parent->heading(pni);
        anchorFrame = Matrix3d::ons(h); // heading → col 0; col 1, 2 perpendicular
    }

    meristem.setAnchor(anchorPos);
    meristem.setAnchorFrame(anchorFrame);

    // Register the initial (anchor) node in the Organ base geometry so that
    // the organism node-management infrastructure works correctly.
    double creationTime = parent->getNodeCT(pni) + std::max(delay, 0.0);
    Organ::addNode(Vector3d(0., 0., 0.), parent->getNodeId(pni), creationTime);
}

// ========================================================================== //
//  Organ interface
// ========================================================================== //

std::shared_ptr<Organ> GrassLeaf::copy(std::shared_ptr<Organism> p)
{
    auto gl = std::make_shared<GrassLeaf>(*this); // shallow copy
    gl->plant  = p;
    gl->parent = std::weak_ptr<Organ>();
    gl->param_ = std::make_shared<GrassLeafSpecificParameter>(*param());
    for (size_t i = 0; i < children.size(); i++) {
        gl->children[i] = children[i]->copy(p);
        gl->children[i]->setParent(gl);
    }
    return gl;
}

/**
 * Advances the grass leaf by @p dt days.
 *
 * Growth is divided into three sequential phases:
 *  1. Sheath elongation  – until sheathLength is reached.
 *  2. Blade delay        – waiting period bladeDelay after sheath completes.
 *  3. Blade elongation   – until bladeLength is reached.
 *
 * Each phase calls the respective helper which prepends new nodes into the
 * Meristem via addNodeFront().
 */
void GrassLeaf::simulate(double dt, bool verbose)
{
    oldNumberOfNodes = meristem.size();

    if (!alive) return;

    // std::cout << "Simulating GrassLeaf id=" << getId() << " age=" << age << " dt=" << dt << "\n" << std::flush;

    const GrassLeafSpecificParameter& p = *param();

    // Kill if life time exceeded.
    if (age + dt > p.sheathDuration + p.bladeDelay + p.bladeDuration) {
        dt    = std::max(p.sheathDuration + p.bladeDelay + p.bladeDuration - age, 0.);
        alive = false;
    }
    age += dt;

    if (age <= 0.) return; // still in delay

    // ----- Phase 1: sheath elongation ----- //
    if (sheathLengthGrown < p.sheathLength) {
        double sheathRate = (p.sheathDuration > 0.)
                            ? p.sheathLength / p.sheathDuration
                            : p.sheathLength; // instant if duration = 0
        double dl = std::min(sheathRate * dt, p.sheathLength - sheathLengthGrown);
        if (dl > 0.) growSheath(dl);
    }

    // ----- Phase 2: blade delay ----- //
    // Nothing to grow, just wait.

    // ----- Phase 3: blade elongation ----- //
    if (isSheathComplete() && isBladeEmerged()) {
        double bladeRate = (p.bladeDuration > 0.)
                           ? p.bladeLength / p.bladeDuration
                           : p.bladeLength;
        double dl = std::min(bladeRate * dt, p.bladeLength - bladeLengthGrown);
        if (dl > 0.) growBlade(dl);
    }

    // Sync the length field used by the base class infrastructure.
    length = sheathLengthGrown + bladeLengthGrown;
    active = (length < p.sheathLength + p.bladeLength) && alive;
}

double GrassLeaf::getParameter(std::string name) const
{
    if (name == "sheathLength")  return param()->sheathLength;
    if (name == "sheathDuration")return param()->sheathDuration;
    if (name == "bladeAngle")    return param()->bladeAngle;
    if (name == "bladeWidth")    return param()->bladeWidth;
    if (name == "bladeLength")   return param()->bladeLength;
    if (name == "bladeDelay")    return param()->bladeDelay;
    if (name == "bladeDuration") return param()->bladeDuration;
    return Organ::getParameter(name);
}

std::string GrassLeaf::toString() const
{
    std::ostringstream s;
    s << "GrassLeaf id=" << getId()
      << " age=" << age
      << " sheath=" << sheathLengthGrown << "/" << param()->sheathLength
      << " blade=" << bladeLengthGrown  << "/" << param()->bladeLength
      << " nodes=" << meristem.size() << "\n";
    return s.str();
}

/**
 * Returns the total leaf length (sheath + blade) at a given age.
 * Approximates growth as piecewise linear through the three phases.
 */
double GrassLeaf::calcLength(double age_)
{
    const GrassLeafSpecificParameter& p = *param();
    if (age_ <= 0.) return 0.;

    double sheathRate = (p.sheathDuration > 0.) ? p.sheathLength / p.sheathDuration : p.sheathLength;
    double sheathLen  = std::min(sheathRate * age_, p.sheathLength);
    double remaining  = age_ - sheathLen / std::max(sheathRate, 1e-12);
    if (remaining <= 0.) return sheathLen;

    double afterDelay = remaining - p.bladeDelay;
    if (afterDelay <= 0.) return sheathLen; // still in blade delay

    double bladeRate = (p.bladeDuration > 0.) ? p.bladeLength / p.bladeDuration : p.bladeLength;
    double bladeLen  = std::min(bladeRate * afterDelay, p.bladeLength);
    return sheathLen + bladeLen;
}

/**
 * Returns the age at which the leaf has reached @p length.
 * Inverse of calcLength(), piecewise linear approximation.
 */
double GrassLeaf::calcAge(double length_) const
{
    const GrassLeafSpecificParameter& p = *param();
    if (length_ <= 0.) return 0.;

    double sheathRate = (p.sheathDuration > 0.) ? p.sheathLength / p.sheathDuration : 1e9;
    if (length_ <= p.sheathLength) {
        return length_ / sheathRate;
    }

    double ageAtSheath = p.sheathLength / sheathRate;
    double bladeRate   = (p.bladeDuration > 0.) ? p.bladeLength / p.bladeDuration : 1e9;
    double bladeLen    = length_ - p.sheathLength;
    return ageAtSheath + p.bladeDelay + bladeLen / bladeRate;
}

// ========================================================================== //
//  Parameter convenience
// ========================================================================== //

std::shared_ptr<GrassLeafRandomParameter> GrassLeaf::getGrassLeafRandomParameter() const
{
    return std::dynamic_pointer_cast<GrassLeafRandomParameter>(
        getOrganism()->getOrganRandomParameter(Organism::ot_leaf, param()->subType));
}

std::shared_ptr<const GrassLeafSpecificParameter> GrassLeaf::param() const
{
    return std::dynamic_pointer_cast<const GrassLeafSpecificParameter>(param_);
}

// ========================================================================== //
//  Private helpers
// ========================================================================== //

double GrassLeaf::calcAgeAtSheathComplete() const
{
    const GrassLeafSpecificParameter& p = *param();
    double sheathRate = (p.sheathDuration > 0.) ? p.sheathLength / p.sheathDuration : 1e9;
    return p.sheathLength / sheathRate;
}

/**
 * Grows the sheath by @p dl [cm].
 *
 * Sheath segments are added via addNodeFront() so that the oldest node sits at
 * the top of the sheath (ligule) and the most recently added node sits at the
 * base (where the leaf emerges from the stem).  The sheath grows vertically
 * relative to the anchor frame (no rotation – yaw = pitch = roll = 0).
 */
void GrassLeaf::growSheath(double dl)
{
    sheathLengthGrown += dl;
    meristem.addNodeFront(dl, 0., 0., 0.); // straight upward (no rotation)
}

/**
 * Grows the blade by @p dl [cm].
 *
 * The blade segment is prepended at the base of the current polyline.  The
 * first blade node applies the blade angle (a pitch rotation relative to the
 * sheath direction) and the blade width is stored in the parameter set but not
 * encoded in the polyline geometry itself (it is used by the visualiser).
 *
 * All subsequent blade increments continue straight (no additional rotation).
 */
void GrassLeaf::growBlade(double dl)
{
    const GrassLeafSpecificParameter& p = *param();
    // On the very first blade increment apply the blade angle as a yaw.
    double yaw = (bladeLengthGrown == 0.) ? p.bladeAngle : 0.;
    bladeLengthGrown += dl;
    meristem.addNodeFront(dl, yaw, 0., 0.);
}

} // end namespace CPlantBox
