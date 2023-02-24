// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef SDF_RS_H
#define SDF_RS_H

#include <iostream>
#include <vector>
#include <stdexcept>

#include "aabbcc/AABB.h"

#include "sdf.h"
#include "Organism.h"
#include "mymath.h"
#include "SegmentAnalyser.h"
#include "Root.h"

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

    SDF_RootSystem(std::vector<Vector3d> nodes, const std::vector<Vector2i> segments, std::vector<double> radii, double dx = 0.5);

    virtual double getDist(const Vector3d& p) const override;

    virtual std::string toString() const override { return "SDF_RootSystem"; }

    std::vector<Vector3d> nodes_;
    std::vector<Vector2i> segments_;
    std::vector<double> radii_;
    double dx_;

protected:

    void buildTree();

    mutable aabb::Tree tree = aabb::Tree();

};



/**
 * Constructors
 */
SDF_RootSystem::SDF_RootSystem(const Root& r, double dx): dx_(dx) {
  size_t n = r.getNumberOfNodes();
  nodes_.resize(n);
  segments_.resize(n-1);
  radii_.resize(n-1);
  for (size_t i=0; i<n; i++) {
      nodes_[i] = r.getNode(i);
  }
  for (size_t i=1; i<n; i++) {
      segments_[i-1] = Vector2i(i-1,i);
      radii_[i-1] = r.param()->a;
  }
  buildTree();
}

SDF_RootSystem::SDF_RootSystem(const Organism& plant, double dx): dx_(dx) {
    auto ana = SegmentAnalyser(plant);
    nodes_ = ana.nodes;
    segments_ = ana.segments;
    radii_ = ana.getParameter("radius");
    buildTree();
}

SDF_RootSystem::SDF_RootSystem(std::vector<Vector3d> nodes, const std::vector<Vector2i> segments, const std::vector<double> radii, double dx)
    :nodes_(nodes), segments_(segments), radii_(radii), dx_(dx) {
    buildTree();
}

void SDF_RootSystem::buildTree() {
    size_t c = 0;
    for (const auto& s : segments_) { // fill the tree
        Vector3d mid = nodes_[s.x].plus(nodes_[s.y]).times(0.5);
        std::vector<double> d = { mid.x, mid.y, mid.z };
        tree.insertParticle(c, d, radii_[c]);
        c++;
    }
}


double SDF_RootSystem::getDist(const Vector3d& p) const {

    std::vector<double> a = { p.x-dx_, p.y-dx_, p.z-dx_ };
    std::vector<double> b = { p.x+dx_, p.y+dx_, p.z+dx_ };
    aabb::AABB box = aabb::AABB(a,b);
    double mdist = 1e100; // far far away
    auto indices = tree.query(box);
    // std::cout << indices.size() << " segments in range\n";
    for (int i : indices) {

        Vector3d x1 = nodes_[segments_[i].x];
        Vector3d x2 = nodes_[segments_[i].y];
        Vector3d v = x2.minus(x1);
        Vector3d w = p.minus(x1);

        double c1 = v.times(w);
        double c2 = v.times(v);

        double l;
        if (c1<=0) {
            l = w.length();
        } else if (c1>=c2) {
            l = p.minus(x2).length();
        } else {
            l = p.minus(x1.plus(v.times(c1/c2))).length();
        }
        l -= radii_[i];

        if (l < mdist) {
            mdist = l;
        }

    }
    return -mdist;
}

} // namespace


#endif
