#pragma once

#include "oktal/geometry/Box.hpp"
#include "oktal/geometry/Vec.hpp"
#include "oktal/octree/MortonIndex.hpp"
#include <cstddef>

namespace oktal {

class OctreeGeometry {

public:
  /**
   * @brief Construct a new Octree Geometry object
   *
   * @param origin
   * @param sidelength
   */
  OctreeGeometry(const Vec3D &origin, double sidelength) noexcept
      : origin_(origin), sidelength_(sidelength) {}

  /**
   * @brief Construct a new Octree Geometry object assigned to the unitcube
   *
   */
  OctreeGeometry() : origin_(0), sidelength_(1) {};

  [[nodiscard]]
  /**
   * @brief returns the origin of the octree
   *
   * @return Vec3D
   */
  Vec3D origin() const {
    return origin_;
  }

  /**
   * @brief returns the sidelength of of the octree
   *
   * @return double
   */
  [[nodiscard]]
  double sidelength() const {
    return sidelength_;
  }

  /**
   * @brief returns the sidelength of cells on level @p level
   *
   * @param level
   * @return double
   */
  [[nodiscard]]
  double dx(size_t level) const {
    return (sidelength_ / static_cast<double>(1ULL << level));
  }

  /**
   * @brief returns a vec3D with the size of a cell on level @p level
   *
   * @param level
   * @return Vec3D
   */
  [[nodiscard]]
  Vec3D cellExtents(size_t level) const {
    const double d = dx(level);
    return Vec3D{d, d, d};
  }

  /**
   * @brief return the bottom-south-west corner of a given cell
   *
   * @param m
   * @return Vec3D
   */
  [[nodiscard]]
  Vec3D cellMinCorner(const MortonIndex &m) const;

  /**
   * @brief return the top-nord-est corner of a given cell
   *
   * @param m
   * @return Vec3D
   */
  [[nodiscard]]
  Vec3D cellMaxCorner(const MortonIndex &m) const;

  /**
   * @brief returns the volume a given cell
   *
   * @param m
   * @return double
   */
  [[nodiscard]]
  Box<double> cellBoundingBox(const MortonIndex &m) const;

  /**
   * @brief return the center a given cell
   *
   * @param m
   * @return Vec3D
   */
  [[nodiscard]]
  Vec3D cellCenter(const MortonIndex &m) const;

private:
  Vec3D origin_;

  // firstly i thought that i should use a Vec to represent the sidelength in
  // each direction but it s a cube , then side_x = side_y = side_z;
  double sidelength_;
};

} // namespace oktal