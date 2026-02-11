#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <numeric>
#include <type_traits>

namespace oktal {

/**
 * @brief Concept to ensure that @p T is either an integral or floating-point
 * type.
 *
 * @tparam T Type to check.
 */
template <typename T>
concept Integral_or_FloatingPoint = std::integral<T> || std::floating_point<T>;

/**
 * @brief Represents a fixed-size algebraic vector with elements of type @p T
 * and dimension @p DIM.
 *
 * Provides element-wise arithmetic operations, magnitude computation, and type
 * conversion.
 *
 * @tparam T   Element type (integral or floating point).
 * @tparam DIM Number of dimensions (vector length).
 */
template <Integral_or_FloatingPoint T, std::size_t DIM> class Vec {
public:
  /**
   * @brief Default constructor.
   *
   * Initializes all elements to zero.
   */
  Vec() noexcept : v_{} {};

  /**
   * @brief Constructs a vector where all elements are initialized to @p
   * initValue.
   *
   * @param initValue The value assigned to all vector elements.
   */
  Vec(T initValue) noexcept : v_{} { v_.fill(initValue); }

  /**
   * @brief Constructs a vector using an initializer list.
   *
   * @details If @p init_list contains fewer than @ref DIM elements, the
   * remaining values are set to 0. If @p init_list contains more than @ref DIM
   * elements, the excess values are ignored.
   *
   * @param init_list List of values to initialize the vector.
   */
  Vec(std::initializer_list<T> init_list) noexcept : v_{} {
    const size_t copy_size = std::min(DIM, init_list.size());
    std::copy_n(init_list.begin(), copy_size, v_.begin());
  }

  /**
   * @brief Constructs a vector by copying and converting elements from another
   * vector.
   *
   * @tparam S Source element type.
   * @param other Vector to copy from.
   */
  template <typename S>
  Vec(const Vec<S, DIM> &other)
    requires std::convertible_to<S, T>
  {
    std::copy(other.begin(), other.end(), begin());
  }

  /**
   * @brief Returns the number of elements (dimensions) in the vector.
   *
   * @return Number of elements.
   */
  [[nodiscard]] std::size_t size() const { return DIM; }

  /**
   * @brief Returns a pointer to the underlying data (const).
   *
   * @return Const pointer to internal array.
   */
  [[nodiscard]] const T *data() const { return v_.data(); }

  /**
   * @brief Returns a pointer to the underlying data (mutable).
   *
   * @return Mutable pointer to internal array.
   */
  [[nodiscard]] T *data() { return v_.data(); }

  /**
   * @brief Accesses an element by index (mutable).
   *
   * @param idx Index of the element.
   * @return Reference to the element.
   */
  T &operator[](std::size_t idx) { return v_.at(idx); }

  /**
   * @brief Accesses an element by index (const).
   *
   * @param idx Index of the element.
   * @return Const reference to the element.
   */
  const T &operator[](std::size_t idx) const { return v_.at(idx); }

  /**
   * @brief Equality comparison operator.
   *
   * @param other Vector to compare against.
   * @return True if all elements are equal, false otherwise.
   */
  bool operator==(const Vec &other) const { return v_ == other.v_; }

  /**
   * @brief Inequality comparison operator.
   *
   * @param other Vector to compare against.
   * @return True if any element differs, false otherwise.
   */
  bool operator!=(const Vec &other) const { return v_ != other.v_; }

  /**
   * @brief Returns an iterator to the beginning (const).
   *
   * @return Const pointer to the first element.
   */
  [[nodiscard]] auto begin() const { return v_.cbegin(); }

  /**
   * @brief Returns an iterator to the beginning (mutable).
   *
   * @return Mutable pointer to the first element.
   */
  [[nodiscard]] auto begin() { return v_.begin(); }

  /**
   * @brief Returns an iterator to the end (const).
   *
   * @return Const pointer past the last element.
   */
  [[nodiscard]] auto end() const { return v_.cend(); }

  /**
   * @brief Returns an iterator to the end (mutable).
   *
   * @return Mutable pointer past the last element.
   */
  [[nodiscard]] auto end() { return v_.end(); }

  /**
   * @brief Unary negation operator.
   *
   * @return A new vector with all elements negated.
   */
  Vec operator-() const
    requires(!std::is_unsigned_v<T>)
  {
    Vec result{};
    std::transform(begin(), end(), result.begin(), std::negate<T>());
    return result;
  }

  /**
   * @brief Adds two vectors element-wise.
   *
   * @param other Right-hand side vector.
   * @return A new vector containing the element-wise sum.
   */
  Vec operator+(const Vec &other) const {
    Vec result{};
    std::transform(begin(), end(), other.begin(), result.begin(),
                   std::plus<T>());
    return result;
  }

  /**
   * @brief Subtracts two vectors element-wise.
   *
   * @param other Right-hand side vector.
   * @return A new vector containing the element-wise difference.
   */
  Vec operator-(const Vec &other) const {
    Vec result{};
    std::transform(begin(), end(), other.begin(), result.begin(),
                   std::minus<T>());
    return result;
  }

  /**
   * @brief Divides each element of the vector by a scalar.
   *
   * @param alpha Scalar divisor.
   * @return A new vector containing the divided elements.
   */
  Vec operator/(const T alpha) const {
    Vec result{};
    std::transform(begin(), end(), result.begin(),
                   [alpha](T vec_val) { return vec_val / alpha; });
    return result;
  }

  /**
   * @brief Multiplies each element of the vector by a scalar.
   *
   * @param vec Vector to scale.
   * @param alpha Scalar multiplier.
   * @return A new vector containing the scaled elements.
   */
  friend Vec operator*(const Vec &vec, const T alpha) {
    Vec result{};
    std::transform(vec.begin(), vec.end(), result.begin(),
                   [alpha](T vec_val) { return vec_val * alpha; });
    return result;
  }

  /**
   * @brief Multiplies each element of the vector by a scalar (commutative
   * version).
   *
   * @param alpha Scalar multiplier.
   * @param vec Vector to scale.
   * @return A new vector containing the scaled elements.
   */
  friend Vec operator*(const T alpha, const Vec &vec) {
    Vec result{};
    std::transform(vec.begin(), vec.end(), result.begin(),
                   [alpha](T vec_val) { return vec_val * alpha; });
    return result;
  }

  /**
   * @brief Adds another vector to this one in-place.
   *
   * @param other Right-hand side vector.
   * @return Reference to this vector after modification.
   */
  Vec &operator+=(const Vec &other) {
    std::transform(begin(), end(), other.begin(), begin(), std::plus<T>());
    return *this;
  }

  /**
   * @brief Subtracts another vector from this one in-place.
   *
   * @param other Right-hand side vector.
   * @return Reference to this vector after modification.
   */
  Vec &operator-=(const Vec &other) {
    std::transform(begin(), end(), other.begin(), begin(), std::minus<T>());
    return *this;
  }

  /**
   * @brief Multiplies this vector by a scalar in-place.
   *
   * @param alpha Scalar multiplier.
   * @return Reference to this vector after modification.
   */
  Vec &operator*=(const T alpha) {
    std::transform(begin(), end(), begin(),
                   [alpha](T vec_val) { return vec_val * alpha; });
    return *this;
  }

  /**
   * @brief Divides this vector by a scalar in-place.
   *
   * @param alpha Scalar divisor.
   * @return Reference to this vector after modification.
   */
  Vec &operator/=(const T alpha) {
    std::transform(begin(), end(), begin(),
                   [alpha](T vec_val) { return vec_val / alpha; });
    return *this;
  }

  /**
   * @brief Computes the magnitude (Euclidean length) of the vector.
   *
   * @return The magnitude of the vector.
   */
  [[nodiscard]] T magnitude() const
    requires std::is_floating_point_v<T>
  {
    return std::sqrt(sqrMagnitude());
  }

  /**
   * @brief Computes the squared magnitude of the vector.
   *
   * @return The squared magnitude (sum of squared elements).
   */
  [[nodiscard]] T sqrMagnitude() const {
    return std::inner_product(begin(), end(), begin(), T(0));
  }

private:
  /** @brief Internal storage array for vector elements. */
  std::array<T, DIM> v_;
};

/// @brief 3D vector using double precision floating-point values.
using Vec3D = Vec<double, 3>;

/// @brief 3D vector using single precision floating-point values.
using Vec3F = Vec<float, 3>;

} // namespace oktal
