// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef SDF_RS_H
#define SDF_RS_H

#include <iostream>
#include <vector>
#include <stdexcept>

#include "mymath.h"
#include "Organism.h"
#include "Root.h"
#include "Hyphae.h"

#include "aabbcc/AABB.h"

namespace CPlantBox {

/**
 * Distance to a static root system
 *
 * segment, nodes, and radii are copied,
 * segment centers are put into a aabb tree, for fast distance lookup
 *
 * dx is the rectangular observation radius
 */
class SDF_RootSystem : public SignedDistanceFunction
{
public:

    SDF_RootSystem(const Root& r, double dx = 0.5);

    SDF_RootSystem(const Organism& plant, double dx = 0.5);

    SDF_RootSystem(const Hyphae& h, double dx = 0.5);

    SDF_RootSystem(std::vector<Vector3d> nodes, const std::vector<Vector2i> segments, std::vector<double> radii, double dx = 0.5);

    virtual double getDist(const Vector3d& p) const override;

    virtual std::string toString() const override { return "SDF_RootSystem"; }

    std::vector<Vector3d> nodes_;
    std::vector<Vector2i> segments_;
    std::vector<double> radii_;
    std::vector<int> organTypes_;
    std::vector<int> treeIds_;
    std::vector<std::weak_ptr<Organ>> segO;
    double dx_;

    int selectedOrganType = -1; // getDist will compare only to this organTypes, -1 for all
    int excludeTreeId = -1; // getDist will exclude this hyphal tree from the search (only if selectedOrganType >=0)
    mutable int distIndex = -1; // last apical node index (s.y) of getDist call (connected to segment [s.x, s.y]) , -1 if
    mutable std::weak_ptr<Organ> lastOrgan;

protected:

    void buildTree();

    mutable aabb::Tree tree = aabb::Tree();

};

} // namespace



#endif
