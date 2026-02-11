#pragma once
#include "Vec.hpp"
#include <type_traits>

namespace oktal {

/**
 * @brief Concept for scalar types (arithmetic types).
 */
template <typename T>
concept Scalar = std::is_arithmetic_v<T>; // Arithmetic implies also scalar

/**
 * @brief Concept to exclude integer types.
 */
template <typename T>
concept NoInteger = !std::is_integral_v<T>;

/**
 * @brief Represents an axis-aligned bounding box in 3D space.
 *
 * @tparam T Scalar type used for coordinates (default: double).
 */
template <Scalar T = double> class Box {
public:
  /// Type alias for 3D vector type
  using vector_type = Vec<T, 3>;

  /**
   * @brief Default constructor. Initializes min and max corners to zero.
   */
  Box() : min_corner(0), max_corner(0) {};

  /**
   * @brief Constructs a box from minimum and maximum corners.
   *
   * @param _min_corner Minimum corner of the box.
   * @param _max_corner Maximum corner of the box.
   */
  Box(vector_type _min_corner, vector_type _max_corner)
      : min_corner(_min_corner), max_corner(_max_corner) {};

  /**
   * @brief Returns a reference to the minimum corner of the box.
   *
   * Works for const and non-const objects.
   *
   * @return Reference to min_corner.
   */
  [[nodiscard]]
  auto minCorner(this auto &&box) -> decltype(auto) {
    return (box.min_corner);
  }

  /**
   * @brief Returns a reference to the maximum corner of the box.
   *
   * Works for const and non-const objects.
   *
   * @return Reference to max_corner.
   */
  [[nodiscard]]
  auto maxCorner(this auto &&box) -> decltype(auto) {
    return (box.max_corner);
  }

  /**
   * @brief Computes the center point of the box.
   *
   * @return 3D vector representing the center.
   */
  [[nodiscard]]
  vector_type center() const {
    return vector_type((min_corner + max_corner) / 2);
  }

  /**
   * @brief Computes the volume of the box.
   *
   * @return Volume as a scalar of type T.
   */
  [[nodiscard]]
  T volume() const {
    return (max_corner[0] - min_corner[0]) * (max_corner[1] - min_corner[1]) *
           (max_corner[2] - min_corner[2]);
  }

  /**
   * @brief Computes the extents (width, height, depth) of the box.
   *
   * @return 3D vector representing the extents along each axis.
   */
  [[nodiscard]]
  vector_type extents() const {
    return vector_type{max_corner[0] - min_corner[0],
                       max_corner[1] - min_corner[1],
                       max_corner[2] - min_corner[2]};
  }

private:
  vector_type min_corner; ///< Minimum corner of the box
  vector_type max_corner; ///< Maximum corner of the box
};

} // namespace oktal