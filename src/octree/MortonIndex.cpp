#include "oktal/octree/MortonIndex.hpp"

#include <bitset>
#include <format>

namespace oktal {
/* namespace {
constexpr std::array<size_t, MortonIndex::MAX_DEPTH> getStartIndices_() {
  std::array<size_t, MortonIndex::MAX_DEPTH> mortonStartIndices{0};

  size_t offset = 1;
  for (size_t i = 1; i < mortonStartIndices.size(); ++i) {
    mortonStartIndices.at(i) = offset + mortonStartIndices.at(i - 1);
    offset = offset << 3;
  }

  return mortonStartIndices;
}

constexpr auto startIndexArray = getStartIndices_();
}; // namespace

[[nodiscard]] const std::array<size_t, MortonIndex::MAX_DEPTH> &
MortonIndex::startIndices() noexcept {
  return startIndexArray;
} */

[[nodiscard]] MortonIndex
MortonIndex::fromPath(const std::vector<morton_bits_t> &path) {
  if (path.size() > MAX_DEPTH) {
    throw std::invalid_argument(std::format(
        "The given path of length={} exceeds the maximum length {}",
        uint64_t(path.size()), uint64_t(oktal::MortonIndex::MAX_DEPTH)));
  }

  size_t shift = path.size();
  shift += shift << 1;
  morton_bits_t bits = 1 << shift;

  for (const morton_bits_t &choice : path) {

    if ((choice & 7) != choice) {
      const uint64_t index = path.size() - (shift / 3);
      throw std::invalid_argument(
          std::format("Choice {} is invalid. Bits are: {:#b}", index, choice));
    }

    shift -= 3;
    bits |= choice << shift;
  }

  return {bits};
}

[[nodiscard]] std::vector<morton_bits_t> MortonIndex::getPath() const {
  std::vector<morton_bits_t> path;
  auto shift = static_cast<size_t>(std::bit_width(bits) - 1);
  const size_t choice_count = shift / 3;
  path.resize(choice_count);

  for (morton_bits_t &choice : path) {
    shift -= 3;
    choice = 7 & (bits >> shift);
  }

  return path;
}

[[nodiscard]] bool
MortonIndex::operator>(const MortonIndex &other) const noexcept {
  const int lhv_width = std::bit_width(bits);
  const int rhv_width = std::bit_width(other.bits);
  if (operator==(other) || lhv_width >= rhv_width) {
    return false;
  }

  const int width_diff = rhv_width - lhv_width;
  const morton_bits_t ancestor = other.bits >> width_diff;

  return bits == ancestor;
}

[[nodiscard]] bool
MortonIndex::operator<(const MortonIndex &other) const noexcept {
  const int lhv_width = std::bit_width(bits);
  const int rhv_width = std::bit_width(other.bits);
  if (operator==(other) || lhv_width <= rhv_width) {
    return false;
  }

  const int width_diff = lhv_width - rhv_width;
  const morton_bits_t ancestor = bits >> width_diff;

  return other.bits == ancestor;
}

[[nodiscard]] bool
MortonIndex::operator>=(const MortonIndex &other) const noexcept {
  const int lhv_width = std::bit_width(bits);
  const int rhv_width = std::bit_width(other.bits);
  if (lhv_width > rhv_width) {
    return false;
  }

  const int width_diff = rhv_width - lhv_width;
  const morton_bits_t ancestor = other.bits >> width_diff;

  return bits == ancestor;
}

[[nodiscard]] bool
MortonIndex::operator<=(const MortonIndex &other) const noexcept {
  const int lhv_width = std::bit_width(bits);
  const int rhv_width = std::bit_width(other.bits);
  if (lhv_width < rhv_width) {
    return false;
  }

  const int width_diff = lhv_width - rhv_width;
  const morton_bits_t ancestor = bits >> width_diff;

  return other.bits == ancestor;
}

[[nodiscard]] UnsignedGridCoordinates MortonIndex::gridCoordinates() const {
  if (isRoot()) {
    return {};
  }

  const UnsignedGridCoordinates parentCoordinates = parent().gridCoordinates();
  const auto localIndex = siblingIndex();
  const UnsignedGridCoordinates localCoordinates{
      (localIndex & 1) == 1 ? 1uz : 0uz, (localIndex & 2) == 2 ? 1uz : 0uz,
      (localIndex & 4) == 4 ? 1uz : 0uz};

  return UnsignedGridCoordinates{parentCoordinates[0] << 1,
                                 parentCoordinates[1] << 1,
                                 parentCoordinates[2] << 1} +
         localCoordinates;
}

[[nodiscard]] MortonIndex
MortonIndex::fromGridCoordinates(const size_t &refinementLevel,
                                 const UnsignedGridCoordinates &coordinates) {
  if (refinementLevel == 0uz) {
    return {};
  }
  const std::array<size_t, 3> vc{coordinates[0], coordinates[1],
                                 coordinates[2]};

  const MortonIndex parent = fromGridCoordinates(
      refinementLevel - 1, {vc[0] >> 1, vc[1] >> 1, vc[2] >> 1});
  const morton_bits_t localIndex =
      ((vc[2] & 1uz) << 2) | ((vc[1] & 1uz) << 1) | (vc[0] & 1uz);

  return parent.child(localIndex);
}
}; // namespace oktal