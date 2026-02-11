#pragma once

#include "oktal/geometry/Vec.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace oktal {

using morton_bits_t = uint64_t;
using UnsignedGridCoordinates = Vec<size_t, 3>;
using SignedGridCoordinates = Vec<ptrdiff_t, 3>;

class MortonIndex {
  morton_bits_t bits;

public:
  static const size_t MAX_DEPTH = (sizeof(morton_bits_t) * 8) / 3;

  MortonIndex() noexcept : bits(1) {}
  MortonIndex(const morton_bits_t &bits_) noexcept : bits(bits_) {}

  [[nodiscard]] morton_bits_t getBits() const noexcept { return bits; }

  [[nodiscard]] static MortonIndex
  fromPath(const std::vector<morton_bits_t> &path);

  [[nodiscard]] std::vector<morton_bits_t> getPath() const;

  /* [[nodiscard]] static const std::array<size_t, MAX_DEPTH> &
  startIndices() noexcept; */

  [[nodiscard]] size_t level() const noexcept {
    return static_cast<size_t>((std::bit_width(bits) - 1) / 3);
  }

  [[nodiscard]] bool isRoot() const noexcept { return bits == 1; }

  [[nodiscard]] size_t siblingIndex() const noexcept {
    if (isRoot()) [[unlikely]] {
      return 0;
    }
    return bits & 7;
  }

  [[nodiscard]] bool isFirstSibling() const noexcept {
    return siblingIndex() == 0;
  }

  [[nodiscard]] bool isLastSibling() const noexcept {
    if (isRoot()) [[unlikely]] {
      return true;
    }
    return siblingIndex() == 7;
  }

  [[nodiscard]] MortonIndex parent() const noexcept { return {bits >> 3}; }

  [[nodiscard]] MortonIndex safeParent() const {
    if (isRoot()) [[unlikely]] {
      throw std::logic_error("Index points to root");
    }
    return parent();
  }

  [[nodiscard]] MortonIndex child(const morton_bits_t &index) const noexcept {
    return {(bits << 3) | index};
  }

  [[nodiscard]] MortonIndex safeChild(const morton_bits_t &index) const {
    if ((bits >> (sizeof(morton_bits_t) * 8 - 1)) == 1) [[unlikely]] {
      throw std::logic_error("Child would exceed maximum depth");
    }
    return child(index);
  }

  [[nodiscard]] bool operator==(const MortonIndex &other) const noexcept {
    return bits == other.bits;
  }

  [[nodiscard]] bool operator!=(const MortonIndex &other) const noexcept {
    return bits != other.bits;
  }

  [[nodiscard]] bool operator>(const MortonIndex &other) const noexcept;

  [[nodiscard]] bool operator<(const MortonIndex &other) const noexcept;

  [[nodiscard]] bool operator>=(const MortonIndex &other) const noexcept;

  [[nodiscard]] bool operator<=(const MortonIndex &other) const noexcept;

  [[nodiscard]] UnsignedGridCoordinates gridCoordinates() const;

  [[nodiscard]] static MortonIndex
  fromGridCoordinates(const size_t &refinementLevel,
                      const UnsignedGridCoordinates &coordinates);
};

} // namespace oktal