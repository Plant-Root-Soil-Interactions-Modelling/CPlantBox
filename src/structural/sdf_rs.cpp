// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "sdf_rs.h"

#include <iostream>
#include <vector>
#include <stdexcept>

#include "sdf.h"
#include "Organism.h"
#include "mymath.h"
#include "SegmentAnalyser.h"
#include "Root.h"
#include "Hyphae.h"

namespace CPlantBox {


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

SDF_RootSystem::SDF_RootSystem(const Hyphae& h, double dx): dx_(dx) {
  size_t n = h.getNumberOfNodes();
  nodes_.resize(n);
  segments_.resize(n-1);
  radii_.resize(n-1);
  for (size_t i=0; i<n; i++) {
      nodes_[i] = h.getNode(i);
  }
  for (size_t i=1; i<n; i++) {
      segments_[i-1] = Vector2i(i-1,i);
      radii_[i-1] = h.param()->a;
  }
  buildTree();
}

SDF_RootSystem::SDF_RootSystem(const Organism& plant, double dx): dx_(dx) {
    auto ana = SegmentAnalyser(plant);
    nodes_ = ana.nodes;
    segments_ = ana.segments;
    radii_ = ana.getParameter("radius");    
    auto vd = ana.getParameter("organType");
    organTypes_.resize(vd.size());
    std::transform(vd.begin(), vd.end(), organTypes_.begin(), [](double x) { return static_cast<int>(x); });   
    vd = ana.getParameter("hyphalTreeIndex");
    treeIds_.resize(vd.size());
    std::transform(vd.begin(), vd.end(), treeIds_.begin(), [](double x) { return static_cast<int>(x); });
    segO = ana.segO; // weak pointer to the organ containing the segment
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

void SDF_RootSystem::updateTree(const Organism& plant) {
    size_t c = tree.nParticles();
    auto newNodes = plant.getNewNodes();
    nodes_.insert(nodes_.end(), newNodes.begin(),newNodes.end());
    auto newSegments = plant.getNewSegments();
    segments_.insert(segments_.end(),newSegments.begin(),newSegments.end());
    auto newsegO = plant.getNewSegmentOrigins();
    segO.insert(segO.end(),newsegO.begin(),newsegO.end());
    for (const auto& o : newsegO) {
        radii_.push_back(o->getParameter("radius"));
        organTypes_.push_back( o->organType());
        treeIds_.push_back(o->getParameter("hyphalTreeIndex"));
    }
    // For checking implementation
    assert(segments_.size() == radii_.size());
    assert(segments_.size() == organTypes_.size());
    assert(segments_.size() == treeIds_.size());
    assert(segments_.size() == segO.size());
    assert(c + newSegments.size() == segments_.size());

    for (const auto& s:newSegments) {
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
    distIndex = -1;
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
        if (i < 0 || static_cast<size_t>(i) >= treeIds_.size()) {
            std::cout << "BAD TREE INDEX: " << i << std::endl;
            throw std::runtime_error("bad tree index");
        }

        if (selectedOrganType == -1) {
			if (l < mdist) {
				mdist = l;
                distIndex = segments_[i].y;
                lastOrgan = segO[i];
			}
        } else {
        	if (excludeTreeId == -1) {
				if ((l < mdist) && (selectedOrganType == organTypes_[i]))  {
					mdist = l;
                    distIndex = segments_[i].y;
                    lastOrgan = segO[i];
				}
        	} else {
				if ((l < mdist) && (selectedOrganType == organTypes_[i]) && (excludeTreeId != treeIds_.at(i)))  {
					mdist = l;
                    distIndex = segments_[i].y;
                    lastOrgan = segO[i];
				}
        	}
        }
    }
    return -mdist;
}

} // namespace


