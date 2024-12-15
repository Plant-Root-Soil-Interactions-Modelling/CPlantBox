#ifndef _CPLANTBOX_PLANTVISUALIZER_H
#define _CPLANTBOX_PLANTVISUALIZER_H
#pragma once

#include "mymath.h"
#include <vector>
#include <functional>
#include <string>
#include <sstream>
#include <tuple>
#include <map>
#include <memory>
#include <assert.h>
#include <optional>
#include <iostream>



namespace CPlantBox {

  // Forward declaration
  class MappedPlant;
  class Leaf;
  class Organ;
  class Plant;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


class PlantVisualiser {
  public:

public :
  PlantVisualiser();
  PlantVisualiser(const PlantVisualiser& pv);
  PlantVisualiser(std::shared_ptr<MappedPlant> plant);
  virtual ~PlantVisualiser();

  double LeafWidthScaleFactor() const { return this->leaf_width_scale_factor_; }
  void SetLeafWidthScaleFactor(double factor) { this->leaf_width_scale_factor_ = factor; }
  void SetMinimumLeafWidth(double width) {this->leaf_minimum_width_ = std::max(0.0, width);};

  void SetGeometryResolution(int resolution) { this->geometry_resolution_ = resolution; } // set the resolution of the geometry (number of cells in each direction
  void SetLeafResolution(int resolution) { this->leaf_resolution_ = resolution;}
  void SetComputeMidlineInLeaf(bool inCompute) { this->include_midline_in_leaf_ = inCompute; }

  void setPlant(const std::shared_ptr<MappedPlant>& plant);

  /* Geometry Functions for Visual Representation */

  /**
   * Compute the geometry for a single organ
   * @param organId the id of the organ
  */
  void ComputeGeometryForOrgan(int organId);
  
  /**
   * Compute the geometry for all organs of a given type
   * This merges the individual organs into one geometry!
   * @param organType the type of the organ
  */
  void ComputeGeometryForOrganType(int organType, bool clearFirst = true);
  /**
   * Compute the geometry of the whole plant
   * This merges the individual organs into one geometry!
   * These methods cannot be const as that would require returning buffers
   * individually, which would be very inefficient
  */
  void ComputeGeometry();

  /**
   * Reset the geometry buffers
  */
  void ResetGeometry();

  bool HasGeometry() const { return !geometry_.empty(); }


  std::vector<double> GetGeometry() const;
  std::vector<unsigned short> GetGeometryColors() const;
  std::vector<double> GetGeometryNormals() const;
  std::vector<double> GetGeometryTextureCoordinates() const;
  std::vector<unsigned int> GetGeometryIndices() const;
  std::vector<int> GetGeometryNodeIds() const;

    /**
   * Map Property onto vertex colors
   * This will use the complete range of the RBGA color space
   * You need to communicate the min and max values of the property seperately
   * @param property the property to map
   * @param minMax the min and max values of the property
  */
  void MapPropertyToColors(std::vector<double> property, 
                           std::pair<double, double> minMax = {0.0, 1.0});

  void SetVerbose(bool verbose) { this->verbose_ = verbose; }
  void SetAddVerticalLeafOffset(bool add) { this->add_vertical_leaf_offset_ = add; }
  void SetRightPenalty(double penalty) { this->right_penalty_ = penalty; }
  void SetShapeFunction(std::function<double(double)> shape_function) { 
    if(verbose_) std::cout << "Setting shape function" << std::endl;
    this->shape_function_ = shape_function;
  }
  void ClearShapeFunction() { this->shape_function_ = std::nullopt; }

  void SetUseStemInfluence(bool use, double radius) { 
    this->use_stem_influence_ = use;
    this->stem_influence_radius_ = radius;
  }

  void SetUseStemRadiusAsMin(bool use) { this->use_stem_radius_as_min_ = use; }

  void SetNotUseStemInfluence() { this->use_stem_influence_ = false; }

  int GetNumOrgans() const;

  std::string SelfCheck() const;

protected:
  std::shared_ptr<MappedPlant> plant_{nullptr};

  bool include_midline_in_leaf_{true};
  bool verbose_{false};
  bool add_vertical_leaf_offset_{false};
  double right_penalty_{0.1};

  double leaf_width_scale_factor_{1.0};
  double leaf_minimum_width_{0.0};
  double stem_influence_radius_{0.0};
  bool use_stem_influence_{false};
  bool use_stem_radius_as_min_{false};

  /**
   * A private method to build the attachment map for the leaf organs
   * This is used to determine the attachment point of the leaf to the stem
  */
  void BuildAttachmentMap();

  /**
   * A private method to generate a geometry for the leaf
   * Based on the existing parametrization and the template from the parameter file
   * @param leaf the leaf to generate the geometry for
   * @param petiole_zone the zone of the petiole in which no leaf surface is present
   * @param p_o the point offset (mainly internal)
   * @param c_o the cell offset (mainly internal)
   * @note using a shared_ptr as parameter might be less efficient, but only O(1) 
  */
  void GenerateLeafGeometry(std::shared_ptr<Leaf> leaf, unsigned int petiole_zone = 0, unsigned int p_o = 0, unsigned int c_o = 0);

  /**
   * A private method to generate a geometry for the stem
   * Based solely on the node-link structure
   * @param stem the stem to generate the geometry for
   * @param p_o the point offset (mainly internal)
   * @param c_o the cell offset (mainly internal)
  */
  void GenerateStemGeometry(std::shared_ptr<Organ> stem, unsigned int p_o = 0, unsigned int c_o = 0);
  
  /**
   * @brief A private method to generate a geometry for a radial leaf representation
   * @param leaf the leaf to generate the geometry for
   * @param p_o the point offset (mainly internal)
   * @param c_o the cell offset (mainly internal)
   * @note using a shared_ptr as parameter might be less efficient, but only O(1)
   */
    void GenerateRadialLeafGeometry(std::shared_ptr<Leaf> leaf, unsigned int p_o = 0, unsigned int c_o = 0);
    
  std::map<int, std::pair<int, Vector3d>> leaf_attachment_map_;

  /* Geometry buffers */
  std::vector<double> geometry_; // x,y,z coordinates
  std::vector<unsigned short> geometry_colors_; // r,g,b,a
  std::vector<double> geometry_normals_; // x,y,z coordinates
  std::vector<unsigned int> geometry_indices_; // indices for triangles
  std::vector<double> geometry_texture_coordinates_; // u,v coordinates
  std::vector<int> geometry_node_ids_; // the node ids for each vertex
  unsigned int geometry_resolution_{8}; // the resolution of the cylindric geometry
  unsigned int leaf_resolution_{20}; // the resolution of the leaf geometry
  // optional alternative shape defining function that takes [0,1] and produces [0,1]
  
  std::optional<std::function<double(double)>> shape_function_ = std::nullopt;
};

} // namespace CPlantBox


#endif 
