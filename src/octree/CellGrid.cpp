#include "oktal/octree/CellOctree.hpp"
#include "oktal/octree/MortonIndex.hpp"
#include <cstddef>
#include <initializer_list>
#include <numeric>
#include <oktal/octree/CellGrid.hpp>
#include <unordered_map>

namespace std {

// This hashing function is necessary to create an unordered map for
// UnsignedGridCoordinates
template <> struct hash<oktal::UnsignedGridCoordinates> {
  size_t operator()(const oktal::UnsignedGridCoordinates &v) const noexcept {
    const size_t h1 = std::hash<size_t>{}(v[0]);
    const size_t h2 = std::hash<size_t>{}(v[1]);
    const size_t h3 = std::hash<size_t>{}(v[2]);
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
};
} // namespace std

namespace {

// Convert unsigned coordinates to signed
inline oktal::SignedGridCoordinates
toSigned(const oktal::UnsignedGridCoordinates &coords) {
  return {static_cast<std::ptrdiff_t>(coords[0]),
          static_cast<std::ptrdiff_t>(coords[1]),
          static_cast<std::ptrdiff_t>(coords[2])};
}

// Convert signed coordinates to unsigned
inline oktal::UnsignedGridCoordinates
toUnsigned(const oktal::SignedGridCoordinates &coords) {
  return {static_cast<size_t>(coords[0]), static_cast<size_t>(coords[1]),
          static_cast<size_t>(coords[2])};
}

// Check if signed coordinates are valid (non-negative)
inline bool isInvalidCoordinates(const oktal::SignedGridCoordinates &coords) {
  return (coords[0] < 0 || coords[1] < 0 || coords[2] < 0);
}

} // namespace

namespace oktal {

// -------------------------------------------------------------------------
// Cell Grid Implementation
// -------------------------------------------------------------------------
oktal::CellGridBuilder
CellGrid::create(std::shared_ptr<const CellOctree> octree) {
  return {std::move(octree)};
}

// -------------------------------------------------------------------------
// Cell Builder Implementation
// -------------------------------------------------------------------------
CellGridBuilder::CellGridBuilder(std::shared_ptr<const CellOctree> octree)
    : octree_(std::move(octree)) {}

CellGridBuilder &CellGridBuilder::levels(std::span<const std::size_t> lvls) {
  levels_.assign(lvls.begin(), lvls.end());
  return *this;
}

CellGridBuilder &
CellGridBuilder::levels(std::initializer_list<const std::size_t> lvls) {
  levels_.assign(lvls.begin(), lvls.end());
  return *this;
}

CellGridBuilder &
CellGridBuilder::neighborhood(std::span<const AdjacencyOffset> offsets) {
  adjacencyOffsets_.assign(offsets.begin(), offsets.end());
  return *this;
}

CellGridBuilder &CellGridBuilder::neighborhood(
    std::initializer_list<const AdjacencyOffset> offsets) {
  adjacencyOffsets_.assign(offsets.begin(), offsets.end());
  return *this;
}

// NOLINTNEXTLINE
CellGrid CellGridBuilder::build() {

  // DEFAULT: all levels of the octree should be taken
  if (levels_.empty()) {
    levels_.resize(octree_->getLevels().size());
    std::ranges::iota(levels_, 0);
  }

  // DEFAULT: no periodicity
  if (periodicityHandler_ == nullptr) {
    periodicityHandler_ = std::make_unique<NoPeriodicity>(NoPeriodicity());
  }

  const size_t numNonPhantomCells = octree_->numberOfNonPhantomNodes(levels_);

  std::vector<MortonIndex> mortonIndices(numNonPhantomCells);
  std::vector<size_t> streamIndexToEnum(octree_->numberOfNodes(),
                                        CellGrid::NOT_ENUMERATED);

  // Enumerate cells in Z-order via horizontalRange
  size_t idx = 0;
  for (const auto lvl : levels_) {
    for (const auto &cell : octree_->horizontalRange(lvl)) {
      mortonIndices[idx] = cell.mortonIndex();
      streamIndexToEnum[cell.streamIndex()] = idx;
      ++idx;
    }
  }
  std::vector<AdjacencyList> adjacencyLists;

  // If no neighbor offsets specified, no adjacency lists should be created
  if (!adjacencyOffsets_.empty()) {

    adjacencyLists.assign(
        adjacencyOffsets_.size(),
        AdjacencyList(mortonIndices.size(), CellGrid::NO_NEIGHBOR));

    // Cache coordinates for faster neighbor search
    // Note: there can be multiple candidates for coordinates due to multiple
    // levels, so we use multimap
    std::unordered_multimap<UnsignedGridCoordinates, size_t> coordToEnum;
    for (size_t idx = 0; idx < mortonIndices.size(); ++idx) {
      coordToEnum.emplace(mortonIndices[idx].gridCoordinates(), idx);
    }

    for (size_t offsetIdx = 0; offsetIdx < adjacencyOffsets_.size();
         ++offsetIdx) {
      const auto &currentOffset = adjacencyOffsets_.at(offsetIdx);
      for (size_t enumIdx = 0; enumIdx < mortonIndices.size(); ++enumIdx) {
        const auto mortonIdx = mortonIndices.at(enumIdx);

        const auto goalCoordsSigned =
            periodicityHandler_->getNeighborCoordinates(
                toSigned(mortonIdx.gridCoordinates()) + currentOffset,
                mortonIdx.level());

        if (isInvalidCoordinates(goalCoordsSigned)) {
          continue;
        }

        // Cast back to unsigned to compare with possible neighbors
        const auto goalCoords = toUnsigned(goalCoordsSigned);

        if (coordToEnum.contains(goalCoords)) {
          const auto possibleNeighbors = coordToEnum.equal_range(goalCoords);

          // Possible neighbors have the desired goal coordinates
          for (const auto &[coords, neighborIdx] : std::ranges::subrange(
                   possibleNeighbors.first, possibleNeighbors.second)) {
            const MortonIndex &possibleNeighbor = mortonIndices.at(neighborIdx);

            // The true neighbor is at the same level
            if (possibleNeighbor.level() == mortonIdx.level()) {
              adjacencyLists[offsetIdx][enumIdx] = neighborIdx;
              continue;
            }
          }
        }
      }
    }
  }

  return {octree_, mortonIndices, streamIndexToEnum, adjacencyOffsets_,
          adjacencyLists};
}

// -------------------------------------------------------------------------
// PeriodicityMapper Implementation
// -------------------------------------------------------------------------
SignedGridCoordinates
NoPeriodicity::getNeighborCoordinates(SignedGridCoordinates goalCoords,
                                      size_t lvl) const {
  const std::ptrdiff_t gridExtent = 1 << lvl;
  for (const auto &coord : goalCoords) {
    if (coord < 0 || coord >= gridExtent) {
      return {-1, -1, -1};
    }
  }
  return goalCoords;
}

SignedGridCoordinates
Torus::getNeighborCoordinates(SignedGridCoordinates goalCoords,
                              size_t lvl) const {
  SignedGridCoordinates resultCoords;
  const std::ptrdiff_t size = 1 << lvl; // 2^level cells per dimension

  for (size_t i = 0; i < goalCoords.size(); ++i) {
    if (periodic_.at(i)) {
      // wrap into [0, size)
      resultCoords[i] =
          static_cast<ptrdiff_t>((goalCoords[i] % size + size) % size);
    } else {
      // just cast directly
      resultCoords[i] = static_cast<ptrdiff_t>(goalCoords[i]);
    }
  }
  return resultCoords;
}
} // namespace oktal