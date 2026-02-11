#include "oktal/octree/OctreeGeometry.hpp"

namespace oktal {

[[nodiscard]]
Vec3D OctreeGeometry::cellMinCorner(const MortonIndex &m) const {
  const size_t level = m.level();
  const double length = dx(level);
  auto cell_coord = m.gridCoordinates();

  return Vec3D{origin_[0] + length * static_cast<double>(cell_coord[0]),
               origin_[1] + length * static_cast<double>(cell_coord[1]),
               origin_[2] + length * static_cast<double>(cell_coord[2])};
}

[[nodiscard]]
Vec3D OctreeGeometry::cellMaxCorner(const MortonIndex &m) const {
  return Vec3D{cellMinCorner(m) + cellExtents(m.level())};
}

[[nodiscard]]
Box<double> OctreeGeometry::cellBoundingBox(const MortonIndex &m) const {
  return {cellMinCorner(m), cellMaxCorner(m)};
}

Vec3D OctreeGeometry::cellCenter(const MortonIndex &m) const {
  return Vec3D{((cellMaxCorner(m) + cellMinCorner(m)) / 2)};
}

} // namespace oktal