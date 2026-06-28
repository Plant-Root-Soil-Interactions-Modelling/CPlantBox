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

#include "AABB.h"

namespace aabb
{


    Tree::Tree(
               unsigned int nParticles,
               double skinThickness_,
               bool touchIsOverlap_) :
       isPeriodic(false), skinThickness(skinThickness_),
        touchIsOverlap(touchIsOverlap_)
    {
       
		dimension = 3;
        // Initialise the periodicity vector.
        //periodicity.resize(dimension);
        std::fill(periodicity.begin(), periodicity.end(), false);

        // Initialise the tree.
        root = NULL_NODE;
        nodeCount = 0;
        nodeCapacity = nParticles;
		particleMap.resize(nParticles, NULL_NODE);
        nodes.resize(nodeCapacity);

        // Build a linked list for the list of free nodes.
        for (unsigned int i=0;i<nodeCapacity-1;i++)
        {
            nodes[i].next = i + 1;
            nodes[i].height = -1;
        }
        nodes[nodeCapacity-1].next = NULL_NODE;
        nodes[nodeCapacity-1].height = -1;

        // Assign the index of the first free node.
        freeList = 0;
    }

    Tree::Tree(unsigned int nParticles,
               double skinThickness_,
               const std::array<bool, 3>& periodicity_,
               const std::array<double, 3>& boxSize_,               
               bool touchIsOverlap_) :
        skinThickness(skinThickness_),
        periodicity(periodicity_), boxSize(boxSize_),
        touchIsOverlap(touchIsOverlap_)
    {
        


        // Initialise the tree.
        root = NULL_NODE;
        touchIsOverlap = true;
        nodeCount = 0;
        nodeCapacity = nParticles;
		particleMap.resize(nParticles, NULL_NODE);
        nodes.resize(nodeCapacity);

        // Build a linked list for the list of free nodes.
        for (unsigned int i=0;i<nodeCapacity-1;i++)
        {
            nodes[i].next = i + 1;
            nodes[i].height = -1;
        }
        nodes[nodeCapacity-1].next = NULL_NODE;
        nodes[nodeCapacity-1].height = -1;

        // Assign the index of the first free node.
        freeList = 0;

        // Check periodicity.
        isPeriodic = false;
        // for (unsigned int i=0;i<dimension;i++)
        // {
            posMinImage[0] =  0.5*boxSize[0];
            negMinImage[0] = -0.5*boxSize[0];
            posMinImage[1] =  0.5*boxSize[1];
            negMinImage[1] = -0.5*boxSize[1];
            posMinImage[2] =  0.5*boxSize[2];
            negMinImage[2] = -0.5*boxSize[2];

			isPeriodic = (periodicity[0])||(periodicity[1])||(periodicity[2]);

        // }
    }



    unsigned int Tree::allocateNode()
    {
        // Expand the node pool as needed.
        if (freeList == NULL_NODE)
        {
            assert(nodeCount == nodeCapacity);

            // The free list is empty. Rebuild a bigger pool.
            nodeCapacity *= 2;
            nodes.resize(nodeCapacity);

            // Build a linked list for the list of free nodes.
            for (unsigned int i=nodeCount;i<nodeCapacity-1;i++)
            {
                nodes[i].next = i + 1;
                nodes[i].height = -1;
            }
            nodes[nodeCapacity-1].next = NULL_NODE;
            nodes[nodeCapacity-1].height = -1;

            // Assign the index of the first free node.
            freeList = nodeCount;
        }

        // Peel a node off the free list.
        unsigned int node = freeList;
        freeList = nodes[node].next;
        nodes[node].parent = NULL_NODE;
        nodes[node].left = NULL_NODE;
        nodes[node].right = NULL_NODE;
        nodes[node].height = 0;
        //nodes[node].aabb.setDimension(dimension);
        nodeCount++;

        return node;
    }

    void Tree::freeNode(unsigned int node)
    {
        assert(node < nodeCapacity);
        assert(0 < nodeCount);

        nodes[node].next = freeList;
        nodes[node].height = -1;
        freeList = node;
        nodeCount--;
    }

    void Tree::insertParticle(unsigned int particle,
                          std::array<double, 3>& position,
                          double radius)
{
	 if (particle >= particleMap.size())
    {
        size_t newSize = std::max(particle + 1, particle * 2);
		particleMap.resize(newSize, NULL_NODE);
    }else{
		// Make sure the particle doesn't already exist.
		if (particleMap[particle] != NULL_NODE)
		{
			throw std::invalid_argument("[ERROR]: Particle already exists in tree!");
		}
	}

    // Allocate a new node for the particle.
    unsigned int node = allocateNode();

    // AABB size in each dimension.
    std::array<double, 3> size;

    nodes[node].aabb.lowerBound[0] = position[0] - radius;
    nodes[node].aabb.upperBound[0] = position[0] + radius;
    size[0] = nodes[node].aabb.upperBound[0] - nodes[node].aabb.lowerBound[0];

    nodes[node].aabb.lowerBound[1] = position[1] - radius;
    nodes[node].aabb.upperBound[1] = position[1] + radius;
    size[1] = nodes[node].aabb.upperBound[1] - nodes[node].aabb.lowerBound[1];

    nodes[node].aabb.lowerBound[2] = position[2] - radius;
    nodes[node].aabb.upperBound[2] = position[2] + radius;
    size[2] = nodes[node].aabb.upperBound[2] - nodes[node].aabb.lowerBound[2];

    // Fatten the AABB.
    nodes[node].aabb.lowerBound[0] -= skinThickness * size[0];
    nodes[node].aabb.upperBound[0] += skinThickness * size[0];

    nodes[node].aabb.lowerBound[1] -= skinThickness * size[1];
    nodes[node].aabb.upperBound[1] += skinThickness * size[1];

    nodes[node].aabb.lowerBound[2] -= skinThickness * size[2];
    nodes[node].aabb.upperBound[2] += skinThickness * size[2];

    nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();
    nodes[node].aabb.computeCentre();

    nodes[node].height = 0;

    // Insert into BVH
    insertLeaf(node);

    // Store reverse mapping (array instead of unordered_map)
    particleMap[particle] = node;

    // Store particle in node (still useful for queries)
    nodes[node].particle = particle;
}

    void Tree::insertParticle(unsigned int particle, std::array<double, 3>& lowerBound, std::array<double, 3>& upperBound)
    {
        // Make sure the particle doesn't already exist.
        if (particle >= particleMap.size())
    {
        size_t newSize = std::max(particle + 1, particle * 2);
		particleMap.resize(newSize, NULL_NODE);
    }else{
		// Make sure the particle doesn't already exist.
		if (particleMap[particle] != NULL_NODE)
		{
			throw std::invalid_argument("[ERROR]: Particle already exists in tree!");
		}
	}

        // // Validate the dimensionality of the bounds vectors.
        // if ((lowerBound.size() != dimension) || (upperBound.size() != dimension))
        // {
            // throw std::invalid_argument("[ERROR]: Dimensionality mismatch!");
        // }

        // Allocate a new node for the particle.
        unsigned int node = allocateNode();


        // // Compute the AABB limits.
        // for (unsigned int i=0;i<dimension;i++)
        // {
            // // Validate the bound.
            // if (lowerBound[i] > upperBound[i])
            // {
                // throw std::invalid_argument("[ERROR]: AABB lower bound is greater than the upper bound!");
            // }

            nodes[node].aabb.lowerBound[0] = lowerBound[0];
            nodes[node].aabb.upperBound[0] = upperBound[0];
            size[0] = upperBound[0] - lowerBound[0];
            
			nodes[node].aabb.lowerBound[1] = lowerBound[1];
            nodes[node].aabb.upperBound[1] = upperBound[1];
            size[1] = upperBound[1] - lowerBound[1];
            
			nodes[node].aabb.lowerBound[2] = lowerBound[2];
            nodes[node].aabb.upperBound[2] = upperBound[2];
            size[2] = upperBound[2] - lowerBound[2];
        // }

        // Fatten the AABB.
        // for (unsigned int i=0;i<dimension;i++)
        // {
            nodes[node].aabb.lowerBound[0] -= skinThickness * size[0];
            nodes[node].aabb.upperBound[0] += skinThickness * size[0];
			
            nodes[node].aabb.lowerBound[1] -= skinThickness * size[1];
            nodes[node].aabb.upperBound[1] += skinThickness * size[1];
			
            nodes[node].aabb.lowerBound[2] -= skinThickness * size[2];
            nodes[node].aabb.upperBound[2] += skinThickness * size[2];
        // }
        nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();
        nodes[node].aabb.computeCentre();//nodes[node].aabb.centre = 

        // Zero the height.
        nodes[node].height = 0;

        // Insert a new leaf into the tree.
        insertLeaf(node);

        // Add the new particle to the map.
        particleMap[particle] = node;

        // Store the particle index.
        nodes[node].particle = particle;
    }

    

    // void Tree::removeParticle(unsigned int particle)
    // {
        // // Map iterator.
        // std::map<unsigned int, unsigned int>::iterator it;

        // // Find the particle.
        // it = particleMap.find(particle);

        // // The particle doesn't exist.
        // if (it == particleMap.end())
        // {
            // throw std::invalid_argument("[ERROR]: Invalid particle index!");
        // }

        // // Extract the node index.
        // unsigned int node = it->second;

        // // Erase the particle from the map.
        // particleMap.erase(it);

        // assert(node < nodeCapacity);
        // assert(nodes[node].isLeaf());

        // removeLeaf(node);
        // freeNode(node);
    // }

    void Tree::removeAll()
    {
        // Iterator pointing to the start of the particle map.
        // std::map<unsigned int, unsigned int>::iterator it = particleMap.begin();

        // Iterate over the map.
        for (unsigned int& node : particleMap)
        {

            assert(node < nodeCapacity);
            assert(nodes[node].isLeaf());

            removeLeaf(node);
            freeNode(node);

            // it++;
        }

        // Clear the particle map.
        particleMap.clear();
    }

    // bool Tree::updateParticle(unsigned int particle, std::array<double, 3>& position, double radius,
                              // bool alwaysReinsert)
    // {
        // // Validate the dimensionality of the position vector.
        // if (position.size() != dimension)
        // {
            // throw std::invalid_argument("[ERROR]: Dimensionality mismatch!");
        // }

        // // AABB bounds vectors.
        // std::array<double, 3> lowerBound;//(dimension);
        // std::array<double, 3> upperBound;//(dimension);

        // // Compute the AABB limits.
        // // for (unsigned int i=0;i<dimension;i++)
        // // {
            // lowerBound[0] = position[0] - radius;
            // upperBound[0] = position[0] + radius;
			
            // lowerBound[1] = position[1] - radius;
            // upperBound[1] = position[1] + radius;
			
            // lowerBound[2] = position[2] - radius;
            // upperBound[2] = position[2] + radius;
        // // }

        // // Update the particle.
        // return updateParticle(particle, lowerBound, upperBound, alwaysReinsert);
    // }

    // bool Tree::updateParticle(unsigned int particle, std::array<double, 3>& lowerBound,
                              // std::array<double, 3>& upperBound, bool alwaysReinsert)
    // {
        // // Validate the dimensionality of the bounds vectors.
        // if ((lowerBound.size() != dimension) && (upperBound.size() != dimension))
        // {
            // throw std::invalid_argument("[ERROR]: Dimensionality mismatch!");
        // }

        // // Map iterator.
        // std::map<unsigned int, unsigned int>::iterator it;

        // // Find the particle.
        // it = particleMap.find(particle);

        // // The particle doesn't exist.
        // if (it == particleMap.end())
        // {
            // throw std::invalid_argument("[ERROR]: Invalid particle index!");
        // }

        // // Extract the node index.
        // unsigned int node = it->second;

        // assert(node < nodeCapacity);
        // assert(nodes[node].isLeaf());

        // // AABB size in each dimension.
        
        // // GCC supports variable length arrays, but MSVC does not.
// #ifdef _MSC_VER
        // std::array<double, 3> size;//(dimension);
// #else
        // double size[dimension];
// #endif

        // // Compute the AABB limits.
        // // for (unsigned int i=0;i<dimension;i++)
        // // {
            // // // Validate the bound.
            // // if (lowerBound[i] > upperBound[i])
            // // {
                // // throw std::invalid_argument("[ERROR]: AABB lower bound is greater than the upper bound!");
            // // }

            // size[0] = upperBound[0] - lowerBound[0];
            // size[1] = upperBound[1] - lowerBound[1];
            // size[2] = upperBound[2] - lowerBound[2];
        // // }

        // // Create the new AABB.
        // AABB aabb(lowerBound, upperBound);

        // // No need to update if the particle is still within its fattened AABB.
        // if (!alwaysReinsert && nodes[node].aabb.contains(aabb)) return false;

        // // Remove the current leaf.
        // removeLeaf(node);

        // // Fatten the new AABB.
        // // for (unsigned int i=0;i<dimension;i++)
        // // {
            // aabb.lowerBound[0] -= skinThickness * size[0];
            // aabb.upperBound[0] += skinThickness * size[0];
			
            // aabb.lowerBound[1] -= skinThickness * size[1];
            // aabb.upperBound[1] += skinThickness * size[1];
			
            // aabb.lowerBound[2] -= skinThickness * size[2];
            // aabb.upperBound[2] += skinThickness * size[2];
        // // }

        // // Assign the new AABB.
        // nodes[node].aabb = aabb;

        // // Update the surface area and centroid.
        // nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();
        // nodes[node].aabb.computeCentre();//nodes[node].aabb.centre = 

        // // Insert a new leaf node.
        // insertLeaf(node);

        // return true;
    // }

    std::vector<unsigned int> Tree::query(unsigned int particle)
    {
        // Make sure that this is a valid particle.
        if (particleMap.size() <= particle)
        {
            throw std::invalid_argument("[ERROR]: Invalid particle index!");
        }

        // Test overlap of particle AABB against all other particles.
        return query(particle, nodes[particleMap[particle]].aabb);
    }
 std::vector<unsigned int> Tree::query(const AABB& aabb)
    {
        std::vector<unsigned int> stack;
        stack.reserve(256);
        stack.push_back(root);

        std::vector<unsigned int> particles;

        while (stack.size() > 0)
        {
            unsigned int node = stack.back();
            stack.pop_back();

            if (node == NULL_NODE) continue;

				const Node& current = nodes[node];
        

            // Test for overlap between the AABBs.
            if (aabb.overlaps(current.aabb, touchIsOverlap))
            {
                // Check that we're at a leaf node.
                if (current.isLeaf())
                {
                    particles.push_back(current.particle);
                }
                else
                {
                    stack.push_back(current.left);
                    stack.push_back(current.right);
                }
            }
        }

        return particles;
    }
	
    std::vector<unsigned int> Tree::query(unsigned int particle, const AABB& aabb)
    {
        std::vector<unsigned int> stack;
        stack.reserve(256);
        stack.push_back(root);

        std::vector<unsigned int> particles;

        while (stack.size() > 0)
        {
            unsigned int node = stack.back();
            stack.pop_back();

            if (node == NULL_NODE) continue;

				const Node& current = nodes[node];
            // Copy the AABB.
            // AABB nodeAABB = current.aabb;
			
            // if (isPeriodic)
            // {
                // std::array<double, 3> separation;//(dimension);
                // std::array<double, 3> shift;//(dimension);
				
				// separation[0] = nodeAABB.centre[0] - aabb.centre[0];
				// separation[1] = nodeAABB.centre[1] - aabb.centre[1];
				// separation[2] = nodeAABB.centre[2] - aabb.centre[2];

                // bool isShifted = minimumImage(separation, shift);

                // // Shift the AABB.
                // if (isShifted)
                // {
                    // nodeAABB.lowerBound[0] += shift[0];
					// nodeAABB.upperBound[0] += shift[0];
                    // nodeAABB.lowerBound[1] += shift[1];
					// nodeAABB.upperBound[1] += shift[1];
                    // nodeAABB.lowerBound[2] += shift[2];
					// nodeAABB.upperBound[2] += shift[2];
                    
                // }
            // }

            // Test for overlap between the AABBs.
            if (aabb.overlaps(current.aabb, touchIsOverlap))
            {
                // Check that we're at a leaf node.
                if (current.isLeaf())
                {
                    // Can't interact with itself.
                    if (nodes[node].particle != particle)
                    {
                        particles.push_back(current.particle);
                    }
                }
                else
                {
                    stack.push_back(current.left);
                    stack.push_back(current.right);
                }
            }
        }

        return particles;
    }


    const AABB& Tree::getAABB(unsigned int particle)
    {
        return nodes[particleMap[particle]].aabb;
    }

    void Tree::insertLeaf(unsigned int leaf)
    {
        if (root == NULL_NODE)
        {
            root = leaf;
            nodes[root].parent = NULL_NODE;
            return;
        }

        // Find the best sibling for the node.

        AABB leafAABB = nodes[leaf].aabb;
        unsigned int index = root;

        while (!nodes[index].isLeaf())
        {
            // Extract the children of the node.
            unsigned int left  = nodes[index].left;
            unsigned int right = nodes[index].right;

            double surfaceArea = nodes[index].aabb.getSurfaceArea();

            // AABB combinedAABB;
            // combinedAABB.merge(nodes[index].aabb, leafAABB);
            // double combinedSurfaceArea = combinedAABB.getSurfaceArea();
			double combinedSurfaceArea = mergedSurfaceArea(nodes[index].aabb, leafAABB);

            // Cost of creating a new parent for this node and the new leaf.
            double cost = 2.0 * combinedSurfaceArea;

            // Minimum cost of pushing the leaf further down the tree.
            double inheritanceCost = 2.0 * (combinedSurfaceArea - surfaceArea);

            // Cost of descending to the left.
            double costLeft;
            if (nodes[left].isLeaf())
            {
                // AABB aabb;
                // aabb.merge(leafAABB, nodes[left].aabb);
                // costLeft = aabb.getSurfaceArea() + inheritanceCost;
				costLeft = mergedSurfaceArea(leafAABB, nodes[left].aabb) + inheritanceCost;
            }
            else
            {
                // AABB aabb;
                // aabb.merge(leafAABB, nodes[left].aabb);
                double oldArea = nodes[left].aabb.getSurfaceArea();
                // double newArea = aabb.getSurfaceArea();
				double newArea = mergedSurfaceArea(leafAABB, nodes[left].aabb) ;
                costLeft = (newArea - oldArea) + inheritanceCost;
            }

            // Cost of descending to the right.
            double costRight;
            if (nodes[right].isLeaf())
            {
                // AABB aabb;
                // aabb.merge(leafAABB, nodes[right].aabb);
                // costRight = aabb.getSurfaceArea() + inheritanceCost;
				costRight = mergedSurfaceArea(leafAABB, nodes[right].aabb) + inheritanceCost;
            }
            else
            {
                // AABB aabb;
                // aabb.merge(leafAABB, nodes[right].aabb);
                double oldArea = nodes[right].aabb.getSurfaceArea();
                // double newArea = aabb.getSurfaceArea();
				double newArea = mergedSurfaceArea(leafAABB, nodes[right].aabb) ;
                costRight = (newArea - oldArea) + inheritanceCost;
            }

            // Descend according to the minimum cost.
            if ((cost < costLeft) && (cost < costRight)) break;

            // Descend.
            if (costLeft < costRight) index = left;
            else                      index = right;
        }

        unsigned int sibling = index;

        // Create a new parent.
        unsigned int oldParent = nodes[sibling].parent;
        unsigned int newParent = allocateNode();
        nodes[newParent].parent = oldParent;
        nodes[newParent].aabb.merge(leafAABB, nodes[sibling].aabb);
        nodes[newParent].height = nodes[sibling].height + 1;

        // The sibling was not the root.
        if (oldParent != NULL_NODE)
        {
            if (nodes[oldParent].left == sibling) nodes[oldParent].left = newParent;
            else                                  nodes[oldParent].right = newParent;

            nodes[newParent].left = sibling;
            nodes[newParent].right = leaf;
            nodes[sibling].parent = newParent;
            nodes[leaf].parent = newParent;
        }
        // The sibling was the root.
        else
        {
            nodes[newParent].left = sibling;
            nodes[newParent].right = leaf;
            nodes[sibling].parent = newParent;
            nodes[leaf].parent = newParent;
            root = newParent;
        }

        // Walk back up the tree fixing heights and AABBs.
        index = nodes[leaf].parent;
        while (index != NULL_NODE)
        {
            index = balance(index);

            unsigned int left = nodes[index].left;
            unsigned int right = nodes[index].right;

            assert(left != NULL_NODE);
            assert(right != NULL_NODE);

            nodes[index].height = 1 + std::max(nodes[left].height, nodes[right].height);
            nodes[index].aabb.merge(nodes[left].aabb, nodes[right].aabb);

            index = nodes[index].parent;
        }
    }

    void Tree::removeLeaf(unsigned int leaf)
    {
        if (leaf == root)
        {
            root = NULL_NODE;
            return;
        }

        unsigned int parent = nodes[leaf].parent;
        unsigned int grandParent = nodes[parent].parent;
        unsigned int sibling;

        if (nodes[parent].left == leaf) sibling = nodes[parent].right;
        else                            sibling = nodes[parent].left;

        // Destroy the parent and connect the sibling to the grandparent.
        if (grandParent != NULL_NODE)
        {
            if (nodes[grandParent].left == parent) nodes[grandParent].left = sibling;
            else                                   nodes[grandParent].right = sibling;

            nodes[sibling].parent = grandParent;
            freeNode(parent);

            // Adjust ancestor bounds.
            unsigned int index = grandParent;
            while (index != NULL_NODE)
            {
                index = balance(index);

                unsigned int left = nodes[index].left;
                unsigned int right = nodes[index].right;

                nodes[index].aabb.merge(nodes[left].aabb, nodes[right].aabb);
                nodes[index].height = 1 + std::max(nodes[left].height, nodes[right].height);

                index = nodes[index].parent;
            }
        }
        else
        {
            root = sibling;
            nodes[sibling].parent = NULL_NODE;
            freeNode(parent);
        }
    }

    unsigned int Tree::balance(unsigned int node)
    {
        assert(node != NULL_NODE);

        if (nodes[node].isLeaf() || (nodes[node].height < 2))
            return node;

        unsigned int left = nodes[node].left;
        unsigned int right = nodes[node].right;

        assert(left < nodeCapacity);
        assert(right < nodeCapacity);

        int currentBalance = nodes[right].height - nodes[left].height;

        // Rotate right branch up.
        if (currentBalance > 1)
        {
            unsigned int rightLeft = nodes[right].left;
            unsigned int rightRight = nodes[right].right;

            assert(rightLeft < nodeCapacity);
            assert(rightRight < nodeCapacity);

            // Swap node and its right-hand child.
            nodes[right].left = node;
            nodes[right].parent = nodes[node].parent;
            nodes[node].parent = right;

            // The node's old parent should now point to its right-hand child.
            if (nodes[right].parent != NULL_NODE)
            {
                if (nodes[nodes[right].parent].left == node) nodes[nodes[right].parent].left = right;
                else
                {
                    assert(nodes[nodes[right].parent].right == node);
                    nodes[nodes[right].parent].right = right;
                }
            }
            else root = right;

            // Rotate.
            if (nodes[rightLeft].height > nodes[rightRight].height)
            {
                nodes[right].right = rightLeft;
                nodes[node].right = rightRight;
                nodes[rightRight].parent = node;
                nodes[node].aabb.merge(nodes[left].aabb, nodes[rightRight].aabb);
                nodes[right].aabb.merge(nodes[node].aabb, nodes[rightLeft].aabb);

                nodes[node].height = 1 + std::max(nodes[left].height, nodes[rightRight].height);
                nodes[right].height = 1 + std::max(nodes[node].height, nodes[rightLeft].height);
            }
            else
            {
                nodes[right].right = rightRight;
                nodes[node].right = rightLeft;
                nodes[rightLeft].parent = node;
                nodes[node].aabb.merge(nodes[left].aabb, nodes[rightLeft].aabb);
                nodes[right].aabb.merge(nodes[node].aabb, nodes[rightRight].aabb);

                nodes[node].height = 1 + std::max(nodes[left].height, nodes[rightLeft].height);
                nodes[right].height = 1 + std::max(nodes[node].height, nodes[rightRight].height);
            }

            return right;
        }

        // Rotate left branch up.
        if (currentBalance < -1)
        {
            unsigned int leftLeft = nodes[left].left;
            unsigned int leftRight = nodes[left].right;

            assert(leftLeft < nodeCapacity);
            assert(leftRight < nodeCapacity);

            // Swap node and its left-hand child.
            nodes[left].left = node;
            nodes[left].parent = nodes[node].parent;
            nodes[node].parent = left;

            // The node's old parent should now point to its left-hand child.
            if (nodes[left].parent != NULL_NODE)
            {
                if (nodes[nodes[left].parent].left == node) nodes[nodes[left].parent].left = left;
                else
                {
                    assert(nodes[nodes[left].parent].right == node);
                    nodes[nodes[left].parent].right = left;
                }
            }
            else root = left;

            // Rotate.
            if (nodes[leftLeft].height > nodes[leftRight].height)
            {
                nodes[left].right = leftLeft;
                nodes[node].left = leftRight;
                nodes[leftRight].parent = node;
                nodes[node].aabb.merge(nodes[right].aabb, nodes[leftRight].aabb);
                nodes[left].aabb.merge(nodes[node].aabb, nodes[leftLeft].aabb);

                nodes[node].height = 1 + std::max(nodes[right].height, nodes[leftRight].height);
                nodes[left].height = 1 + std::max(nodes[node].height, nodes[leftLeft].height);
            }
            else
            {
                nodes[left].right = leftRight;
                nodes[node].left = leftLeft;
                nodes[leftLeft].parent = node;
                nodes[node].aabb.merge(nodes[right].aabb, nodes[leftLeft].aabb);
                nodes[left].aabb.merge(nodes[node].aabb, nodes[leftRight].aabb);

                nodes[node].height = 1 + std::max(nodes[right].height, nodes[leftLeft].height);
                nodes[left].height = 1 + std::max(nodes[node].height, nodes[leftRight].height);
            }

            return left;
        }

        return node;
    }

    unsigned int Tree::computeHeight() const
    {
        return computeHeight(root);
    }

    unsigned int Tree::computeHeight(unsigned int node) const
    {
        assert(node < nodeCapacity);

        if (nodes[node].isLeaf()) return 0;

        unsigned int height1 = computeHeight(nodes[node].left);
        unsigned int height2 = computeHeight(nodes[node].right);

        return 1 + std::max(height1, height2);
    }

    unsigned int Tree::getHeight() const
    {
        if (root == NULL_NODE) return 0;
        return nodes[root].height;
    }

    unsigned int Tree::getNodeCount() const
    {
        return nodeCount;
    }

    unsigned int Tree::computeMaximumBalance() const
    {
        unsigned int maxBalance = 0;
        for (unsigned int i=0; i<nodeCapacity; i++)
        {
            if (nodes[i].height <= 1)
                continue;

            assert(nodes[i].isLeaf() == false);

            unsigned int balance = std::abs(nodes[nodes[i].left].height - nodes[nodes[i].right].height);
            maxBalance = std::max(maxBalance, balance);
        }

        return maxBalance;
    }

    double Tree::computeSurfaceAreaRatio() const
    {
        if (root == NULL_NODE) return 0.0;

        double rootArea = nodes[root].aabb.computeSurfaceArea();
        double totalArea = 0.0;

        for (unsigned int i=0; i<nodeCapacity;i++)
        {
            if (nodes[i].height < 0) continue;

            totalArea += nodes[i].aabb.computeSurfaceArea();
        }

        return totalArea / rootArea;
    }

    void Tree::validate() const
    {
#ifndef NDEBUG
        validateStructure(root);
        validateMetrics(root);

        unsigned int freeCount = 0;
        unsigned int freeIndex = freeList;

        while (freeIndex != NULL_NODE)
        {
            assert(freeIndex < nodeCapacity);
            freeIndex = nodes[freeIndex].next;
            freeCount++;
        }

        assert(getHeight() == computeHeight());
        assert((nodeCount + freeCount) == nodeCapacity);
#endif
    }

    void Tree::rebuild()
    {
        std::vector<unsigned int> nodeIndices(nodeCount);
        unsigned int count = 0;

        for (unsigned int i=0;i<nodeCapacity;i++)
        {
            // Free node.
            if (nodes[i].height < 0) continue;

            if (nodes[i].isLeaf())
            {
                nodes[i].parent = NULL_NODE;
                nodeIndices[count] = i;
                count++;
            }
            else freeNode(i);
        }

        while (count > 1)
        {
            double minCost = std::numeric_limits<double>::max();
            int iMin = -1, jMin = -1;

            for (unsigned int i=0;i<count;i++)
            {
                AABB aabbi = nodes[nodeIndices[i]].aabb;

                for (unsigned int j=i+1;j<count;j++)
                {
                    AABB aabbj = nodes[nodeIndices[j]].aabb;
                    AABB aabb;
                    aabb.merge(aabbi, aabbj);
                    double cost = aabb.getSurfaceArea();

                    if (cost < minCost)
                    {
                        iMin = i;
                        jMin = j;
                        minCost = cost;
                    }
                }
            }

            unsigned int index1 = nodeIndices[iMin];
            unsigned int index2 = nodeIndices[jMin];

            unsigned int parent = allocateNode();
            nodes[parent].left = index1;
            nodes[parent].right = index2;
            nodes[parent].height = 1 + std::max(nodes[index1].height, nodes[index2].height);
            nodes[parent].aabb.merge(nodes[index1].aabb, nodes[index2].aabb);
            nodes[parent].parent = NULL_NODE;

            nodes[index1].parent = parent;
            nodes[index2].parent = parent;

            nodeIndices[jMin] = nodeIndices[count-1];
            nodeIndices[iMin] = parent;
            count--;
        }

        root = nodeIndices[0];

        validate();
    }

    void Tree::validateStructure(unsigned int node) const
    {
        if (node == NULL_NODE) return;

        if (node == root) assert(nodes[node].parent == NULL_NODE);

        unsigned int left = nodes[node].left;
        unsigned int right = nodes[node].right;

        if (nodes[node].isLeaf())
        {
            assert(left == NULL_NODE);
            assert(right == NULL_NODE);
            assert(nodes[node].height == 0);
            return;
        }

        assert(left < nodeCapacity);
        assert(right < nodeCapacity);

        assert(nodes[left].parent == node);
        assert(nodes[right].parent == node);

        validateStructure(left);
        validateStructure(right);
    }

    void Tree::validateMetrics(unsigned int node) const
    {
        if (node == NULL_NODE) return;

        unsigned int left = nodes[node].left;
        unsigned int right = nodes[node].right;

        if (nodes[node].isLeaf())
        {
            assert(left == NULL_NODE);
            assert(right == NULL_NODE);
            assert(nodes[node].height == 0);
            return;
        }

        assert(left < nodeCapacity);
        assert(right < nodeCapacity);

        int height1 = nodes[left].height;
        int height2 = nodes[right].height;
        int height = 1 + std::max(height1, height2);
        (void)height; // Unused variable in Release build
        assert(nodes[node].height == height);

        AABB aabb;
        aabb.merge(nodes[left].aabb, nodes[right].aabb);

        // for (unsigned int i=0;i<dimension;i++)
        // {
            // assert(aabb.lowerBound[i] == nodes[node].aabb.lowerBound[i]);
            // assert(aabb.upperBound[i] == nodes[node].aabb.upperBound[i]);
        // }

        validateMetrics(left);
        validateMetrics(right);
    }

    void Tree::periodicBoundaries(std::array<double, 3>& position)
    {
        // for (unsigned int i=0;i<dimension;i++)
        // {
            if (position[0] < 0)
            {
                position[0] += boxSize[0];
            }
            else
            {
                if (position[0] >= boxSize[0])
                {
                    position[0] -= boxSize[0];
                }
            }
			
			if (position[1] < 0)
            {
                position[1] += boxSize[1];
            }
            else
            {
                if (position[1] >= boxSize[1])
                {
                    position[1] -= boxSize[1];
                }
            }
			
			if (position[2] < 0)
            {
                position[2] += boxSize[2];
            }
            else
            {
                if (position[2] >= boxSize[2])
                {
                    position[2] -= boxSize[2];
                }
            }
        // }
    }

    bool Tree::minimumImage(std::array<double, 3>& separation, std::array<double, 3>& shift)
    {
        bool isShifted = false;

        // for (unsigned int i=0;i<dimension;i++)
        // {
            if (separation[0] < negMinImage[0])
            {
                separation[0] += periodicity[0]*boxSize[0];
                shift[0] = periodicity[0]*boxSize[0];
                isShifted = true;
            }
            else
            {
                if (separation[0] >= posMinImage[0])
                {
                    separation[0] -= periodicity[0]*boxSize[0];
                    shift[0] = -periodicity[0]*boxSize[0];
                    isShifted = true;
                }
            }
			
			if (separation[1] < negMinImage[1])
            {
                separation[1] += periodicity[1]*boxSize[1];
                shift[1] = periodicity[1]*boxSize[1];
                isShifted = true;
            }
            else
            {
                if (separation[1] >= posMinImage[1])
                {
                    separation[1] -= periodicity[1]*boxSize[1];
                    shift[1] = -periodicity[1]*boxSize[1];
                    isShifted = true;
                }
            }
			
			if (separation[2] < negMinImage[2])
            {
                separation[2] += periodicity[2]*boxSize[2];
                shift[2] = periodicity[2]*boxSize[2];
                isShifted = true;
            }
            else
            {
                if (separation[2] >= posMinImage[2])
                {
                    separation[2] -= periodicity[2]*boxSize[2];
                    shift[2] = -periodicity[2]*boxSize[2];
                    isShifted = true;
                }
            }
        // }

        return isShifted;
    }
}
