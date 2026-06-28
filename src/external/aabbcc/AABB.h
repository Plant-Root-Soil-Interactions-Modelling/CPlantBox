/*
  Copyright (c) 2009 Erin Catto http://www.box2d.org
  Copyright (c) 2016-2018 Lester Hedges <lester.hedges+aabbcc@gmail.com>

  This software is provided 'as-is', without any express or implied
  warranty. In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.

  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.

  3. This notice may not be removed or altered from any source distribution.

  This code was adapted from parts of the Box2D Physics Engine,
  http://www.box2d.org
*/

#ifndef _AABB_H
#define _AABB_H

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>
#include <array>
#include <unordered_map>

/// Null node flag.
const unsigned int NULL_NODE = 0xffffffff;

namespace aabb
{
    /*! \brief The axis-aligned bounding box object.

        Axis-aligned bounding boxes (AABBs) store information for the minimum
        orthorhombic bounding-box for an object. Support is provided for
        dimensions >= 2. (In 2D the bounding box is either a rectangle,
        in 3D it is a rectangular prism.)

        Class member functions provide functionality for merging AABB objects
        and testing overlap with other AABBs.
     */
    class AABB
    {
    public:
        /// Constructor.
        AABB(){};

        //! Constructor.
        /*! \param dimension
                The dimensionality of the system.
         */
        AABB(unsigned int dimension)
		{};

        //! Constructor.
        /*! \param lowerBound_
                The lower bound in each dimension.

            \param upperBound_
                The upper bound in each dimension.
         */
        AABB(const std::array<double, 3>& lowerBound_, const std::array<double, 3>& upperBound_) :
        lowerBound(lowerBound_), upperBound(upperBound_)
    {


        // Validate that the upper bounds exceed the lower bounds.
		// that is just check/trouble shooting
        // for (unsigned int i=0;i<lowerBound.size();i++)
        // {
            // // Validate the bound.
            // if (lowerBound[i] > upperBound[i])
            // {
                // throw std::invalid_argument("[ERROR]: AABB lower bound is greater than the upper bound!");
            // }
        // }

        surfaceArea = computeSurfaceArea();
        computeCentre();//centre = 
    };

        /// Compute the surface area of the box.
		inline double computeSurfaceArea() const{
			double dx = upperBound[0] - lowerBound[0];
			double dy = upperBound[1] - lowerBound[1];
			double dz = upperBound[2] - lowerBound[2];
			return 2.0 * (dx*dy + dy*dz + dx*dz);
		};

		
        /// Get the surface area of the box.
        double getSurfaceArea() const {
			return surfaceArea;
		};

        //! Merge two AABBs into this one.
        /*! \param aabb1
                A reference to the first AABB.

            \param aabb2
                A reference to the second AABB.
         */
        inline void merge(const AABB& aabb1, const AABB& aabb2)
    {
		lowerBound[0] = std::min(aabb1.lowerBound[0], aabb2.lowerBound[0]);
		lowerBound[1] = std::min(aabb1.lowerBound[1], aabb2.lowerBound[1]);
		lowerBound[2] = std::min(aabb1.lowerBound[2], aabb2.lowerBound[2]);

		upperBound[0] = std::max(aabb1.upperBound[0], aabb2.upperBound[0]);
		upperBound[1] = std::max(aabb1.upperBound[1], aabb2.upperBound[1]);
		upperBound[2] = std::max(aabb1.upperBound[2], aabb2.upperBound[2]);

        surfaceArea = computeSurfaceArea();
        computeCentre();//centre =
    };

        //! Test whether the AABB is contained within this one.
        /*! \param aabb
                A reference to the AABB.

            \return
                Whether the AABB is fully contained.
         */
        bool contains(const AABB& aabb) const
		{
			return
				aabb.lowerBound[0] >= lowerBound[0] &&
				aabb.lowerBound[1] >= lowerBound[1] &&
				aabb.lowerBound[2] >= lowerBound[2] &&
				aabb.upperBound[0] <= upperBound[0] &&
				aabb.upperBound[1] <= upperBound[1] &&
				aabb.upperBound[2] <= upperBound[2];
		};

        //! Test whether the AABB overlaps this one.
        /*! \param aabb
                A reference to the AABB.

            \param touchIsOverlap
                Does touching constitute an overlap?

            \return
                Whether the AABB overlaps.
         */
        bool overlaps(const AABB& aabb, bool touchIsOverlap) const
		{
        return
            !(aabb.upperBound[0] < lowerBound[0] ||
              aabb.lowerBound[0] > upperBound[0] ||
              aabb.upperBound[1] < lowerBound[1] ||
              aabb.lowerBound[1] > upperBound[1] ||
              aabb.upperBound[2] < lowerBound[2] ||
              aabb.lowerBound[2] > upperBound[2]);
    };

        //! Compute the centre of the AABB.
        /*! \returns
                The position vector of the AABB centre.
         */
        inline void computeCentre()
		 {
		
		//    std::vector<double> position(lowerBound.size());
		centre[0] = 0.5 * (lowerBound[0] + upperBound[0]);
		centre[1] = 0.5 * (lowerBound[1] + upperBound[1]);
		centre[2] = 0.5 * (lowerBound[2] + upperBound[2]);
        //return centre;
    };

        //! Set the dimensionality of the AABB.
        /*! \param dimension
                The dimensionality of the system.
         */
        void setDimension(unsigned int)
		{};

		std::array<double, 3> lowerBound;
		std::array<double, 3> upperBound;
		std::array<double, 3> centre;

        /// The AABB's surface area.
        double surfaceArea;
    };

    /*! \brief A node of the AABB tree.

        Each node of the tree contains an AABB object which corresponds to a
        particle, or a group of particles, in the simulation box. The AABB
        objects of individual particles are "fattened" before they are stored
        to avoid having to continually update and rebalance the tree when
        displacements are small.

        Nodes are aware of their position within in the tree. The isLeaf member
        function allows the tree to query whether the node is a leaf, i.e. to
        determine whether it holds a single particle.
     */
    struct Node
    {
        /// Constructor.
        Node(){};

        /// The fattened axis-aligned bounding box.
        AABB aabb;

        /// Index of the parent node.
        unsigned int parent;

        /// Index of the next node.
        unsigned int next;

        /// Index of the left-hand child.
        unsigned int left;

        /// Index of the right-hand child.
        unsigned int right;

        /// Height of the node. This is 0 for a leaf and -1 for a free node.
        int height;

        /// The index of the particle that the node contains (leaf nodes only).
        unsigned int particle;

        //! Test whether the node is a leaf.
        /*! \return
                Whether the node is a leaf node.
         */
        bool isLeaf() const{
			return (left == NULL_NODE);
		};
    };

    /*! \brief The dynamic AABB tree.

        The dynamic AABB tree is a hierarchical data structure that can be used
        to efficiently query overlaps between objects of arbitrary shape and
        size that lie inside of a simulation box. Support is provided for
        periodic and non-periodic boxes, as well as boxes with partial
        periodicity, e.g. periodic along specific axes.
     */
    class Tree
    {
    public:
        //! Constructor (non-periodic).
        /*! \param dimension_
                The dimensionality of the system.

            \param skinThickness_
                The skin thickness for fattened AABBs, as a fraction
                of the AABB base length.

            \param nParticles
                The number of particles (for fixed particle number systems).

            \param touchIsOverlap
                Does touching count as overlapping in query operations?
         */
        Tree(){};
			
        Tree(unsigned int nParticles, double skinThickness_ = 0.05,
            bool touchIsOverlap=true);

        //! Constructor (custom periodicity).
        /*! \param dimension_
                The dimensionality of the system.

            \param skinThickness_
                The skin thickness for fattened AABBs, as a fraction
                of the AABB base length.

            \param periodicity_
                Whether the system is periodic in each dimension.

            \param boxSize_
                The size of the simulation box in each dimension.

            \param nParticles
                The number of particles (for fixed particle number systems).

            \param touchIsOverlap
                Does touching count as overlapping in query operations?
         */
        Tree( unsigned int nParticles ,double, const std::array<bool, 3>&, const std::array<double, 3>&,
             bool touchIsOverlap=true);

        //! Set the periodicity of the simulation box.
        /*! \param periodicity_
                Whether the system is periodic in each dimension.
         */
        void setPeriodicity(const std::array<bool, 3>& periodicity_)
		{
			periodicity = periodicity_;
		};

        //! Set the size of the simulation box.
        /*! \param boxSize_
                The size of the simulation box in each dimension.
         */
        void setBoxSize(const std::array<double, 3>& boxSize_)
		{
			boxSize = boxSize_;
		};

        //! Insert a particle into the tree (point particle).
        /*! \param index
                The index of the particle.

            \param position
                The position vector of the particle.

            \param radius
                The radius of the particle.
         */
        void insertParticle(unsigned int, std::array<double, 3>&, double);

        //! Insert a particle into the tree (arbitrary shape with bounding box).
        /*! \param index
                The index of the particle.

            \param lowerBound
                The lower bound in each dimension.

            \param upperBound
                The upper bound in each dimension.
         */
        void insertParticle(unsigned int, std::array<double, 3>&, std::array<double, 3>&);

        /// Return the number of particles in the tree.
        unsigned int nParticles()
		{
			return particleMap.size();
		};

        //! Remove a particle from the tree.
        /*! \param particle
                The particle index (particleMap will be used to map the node).
         */
        // void removeParticle(unsigned int);

        /// Remove all particles from the tree.
        void removeAll();

        //! Update the tree if a particle moves outside its fattened AABB.
        /*! \param particle
                The particle index (particleMap will be used to map the node).

            \param position
                The position vector of the particle.

            \param radius
                The radius of the particle.

            \param alwaysReinsert
                Always reinsert the particle, even if it's within its old AABB (default:false)

            \return
                Whether the particle was reinserted.
         */
        // bool updateParticle(unsigned int, std::array<double, 3>&, double, bool alwaysReinsert=false);

        //! Update the tree if a particle moves outside its fattened AABB.
        /*! \param particle
                The particle index (particleMap will be used to map the node).

            \param lowerBound
                The lower bound in each dimension.

            \param upperBound
                The upper bound in each dimension.

            \param alwaysReinsert
                Always reinsert the particle, even if it's within its old AABB (default: false)
         */
        // bool updateParticle(unsigned int, std::array<double, 3>&, std::array<double, 3>&, bool alwaysReinsert=false);

        //! Query the tree to find candidate interactions for a particle.
        /*! \param particle
                The particle index.

            \return particles
                A vector of particle indices.
         */
        std::vector<unsigned int> query(unsigned int);

        //! Query the tree to find candidate interactions for an AABB.
        /*! \param particle
                The particle index.

            \param aabb
                The AABB.

            \return particles
                A vector of particle indices.
         */
        std::vector<unsigned int> query(unsigned int, const AABB&);
        std::vector<unsigned int> query(const AABB& aabb);

        //! Query the tree to find candidate interactions for an AABB.
        /*! \param aabb
                The AABB.

            \return particles
                A vector of particle indices.
         */
        std::vector<unsigned int> queryMax(const AABB& aabb)
    {
        // Make sure the tree isn't empty.
        if (particleMap.size() == 0)
        {
            return std::vector<unsigned int>();
        }

        // Test overlap of AABB against all particles.
        return query(std::numeric_limits<unsigned int>::max(), aabb);
    };

        //! Get a particle AABB.
        /*! \param particle
                The particle index.
         */
        const AABB& getAABB(unsigned int);

        //! Get the height of the tree.
        /*! \return
                The height of the binary tree.
         */
        unsigned int getHeight() const;

        //! Get the number of nodes in the tree.
        /*! \return
                The number of nodes in the tree.
         */
        unsigned int getNodeCount() const;

        //! Compute the maximum balancance of the tree.
        /*! \return
                The maximum difference between the height of two
                children of a node.
         */
        unsigned int computeMaximumBalance() const;

        //! Compute the surface area ratio of the tree.
        /*! \return
                The ratio of the sum of the node surface area to the surface
                area of the root node.
         */
        double computeSurfaceAreaRatio() const;

        /// Validate the tree.
        void validate() const;

        /// Rebuild an optimal tree.
        void rebuild();
		
		static inline double mergedSurfaceArea(const AABB& a, const AABB& b) {
			double dx = std::max(a.upperBound[0], b.upperBound[0]) - std::min(a.lowerBound[0], b.lowerBound[0]);
			double dy = std::max(a.upperBound[1], b.upperBound[1]) - std::min(a.lowerBound[1], b.lowerBound[1]);
			double dz = std::max(a.upperBound[2], b.upperBound[2]) - std::min(a.lowerBound[2], b.lowerBound[2]);
			return 2.0 * (dx*dy + dy*dz + dx*dz);
		}


    private:
        /// The index of the root node.
        unsigned int root;

        /// The dynamic tree.
        std::vector<Node> nodes;

        /// The current number of nodes in the tree.
        unsigned int nodeCount;

        /// The current node capacity.
        unsigned int nodeCapacity;

        /// The position of node at the top of the free list.
        unsigned int freeList;

        /// The dimensionality of the system.
        unsigned int dimension;

        /// Whether the system is periodic along at least one axis.
        bool isPeriodic;

        /// The skin thickness of the fattened AABBs, as a fraction of the AABB base length.
        double skinThickness;

        std::array<double, 3> size;
        /// Whether the system is periodic along each axis.
        std::array<bool, 3> periodicity;

        /// The size of the system in each dimension.
        std::array<double, 3> boxSize;

        /// The position of the negative minimum image.
        std::array<double, 3> negMinImage;

        /// The position of the positive minimum image.
        std::array<double, 3> posMinImage;

        /// A map between particle and node indices.
        std::vector<unsigned int> particleMap;

        /// Does touching count as overlapping in tree queries?
        bool touchIsOverlap;

        //! Allocate a new node.
        /*! \return
                The index of the allocated node.
         */
        unsigned int allocateNode();

        //! Free an existing node.
        /*! \param node
                The index of the node to be freed.
         */
        void freeNode(unsigned int);

        //! Insert a leaf into the tree.
        /*! \param leaf
                The index of the leaf node.
         */
        void insertLeaf(unsigned int);

        //! Remove a leaf from the tree.
        /*! \param leaf
                The index of the leaf node.
         */
        void removeLeaf(unsigned int);

        //! Balance the tree.
        /*! \param node
                The index of the node.
         */
        unsigned int balance(unsigned int);

        //! Compute the height of the tree.
        /*! \return
                The height of the entire tree.
         */
        unsigned int computeHeight() const;

        //! Compute the height of a sub-tree.
        /*! \param node
                The index of the root node.

            \return
                The height of the sub-tree.
         */
        unsigned int computeHeight(unsigned int) const;

        //! Assert that the sub-tree has a valid structure.
        /*! \param node
                The index of the root node.
         */
        void validateStructure(unsigned int) const;

        //! Assert that the sub-tree has valid metrics.
        /*! \param node
                The index of the root node.
         */
        void validateMetrics(unsigned int) const;

        //! Apply periodic boundary conditions.
        /* \param position
                The position vector.
         */
        void periodicBoundaries(std::array<double, 3>&);

        //! Compute minimum image separation.
        /*! \param separation
                The separation vector.

            \param shift
                The shift vector.

            \return
                Whether a periodic shift has been applied.
         */
        bool minimumImage(std::array<double, 3>&, std::array<double, 3>&);
    };
}

#endif /* _AABB_H */
