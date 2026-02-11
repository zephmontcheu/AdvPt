#pragma once

#include "oktal/geometry/Vec.hpp"
#include "oktal/octree/CellOctree.hpp"
#include "oktal/octree/MortonIndex.hpp"
#include <algorithm>
#include <array>
#include <cstddef>
#include <format>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <memory>
#include <source_location>
#include <stdexcept>
#include <type_traits>
namespace oktal {

// -------------------------------------------------------------------------
// Type alias
// -------------------------------------------------------------------------
using AdjacencyOffset = Vec<std::ptrdiff_t, 3>;
using AdjacencyList = std::vector<size_t>;
using AdjacencyListView = std::span<const size_t>;

class CellGridBuilder;

class CellGrid {
public:
  static constexpr size_t NOT_ENUMERATED = std::numeric_limits<size_t>::max();
  static constexpr size_t NO_NEIGHBOR = NOT_ENUMERATED;

  // -------------------------------------------------------------------------
  // Nested Cell View Class
  // -------------------------------------------------------------------------
  class CellView {
  public:
    CellView(const CellGrid *cellGrid, size_t enumIdx)
        : grid_(cellGrid), enumIdx_(enumIdx) {}

    [[nodiscard]]
    size_t enumerationIndex() const {
      return enumIdx_;
    }

    [[nodiscard]]
    bool isValid() const {
      return enumIdx_ != NOT_ENUMERATED;
    }

    [[nodiscard]]
    CellView neighbor(AdjacencyOffset offset) const {
      return {grid_, grid_->neighborIndices(offset)[enumIdx_]};
    }

    [[nodiscard]]
    MortonIndex mortonIndex() const {
      return grid_->mortonIndices_.at(enumIdx_);
    }

    [[nodiscard]]
    size_t level() const {
      return mortonIndex().level();
    }

    [[nodiscard]]
    Vec3D center() const {
      // NOLINTNEXTLINE (bugprone-unchecked-optional-access)
      return grid_->octree_->getCell(grid_->mortonIndices_.at(enumIdx_))
          ->center();
    }

    [[nodiscard]]
    Box<double> boundingBox() const {
      // NOLINTNEXTLINE (bugprone-unchecked-optional-access)
      return grid_->octree_->getCell(grid_->mortonIndices_.at(enumIdx_))
          ->boundingBox();
    }

    // Implicit conversion functions
    operator size_t() const { return enumIdx_; }
    operator bool() const { return isValid(); }

  private:
    const CellGrid *grid_;
    size_t enumIdx_;
  };

  // -------------------------------------------------------------------------
  // Iterator
  // -------------------------------------------------------------------------
  class CellGridIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = CellGrid::CellView;
    using difference_type = std::ptrdiff_t;

    CellGridIterator(const CellGrid *grid, size_t idx)
        : grid_(grid), idx_(idx) {}

    CellGridIterator() = default;

    [[nodiscard]]
    value_type operator*() const {
      return {grid_, idx_};
    }

    CellGridIterator &operator++() {
      ++idx_;
      return *this;
    }

    [[nodiscard]]
    CellGridIterator operator++(int) {
      CellGridIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    [[nodiscard]]
    bool operator==(const CellGridIterator &other) const {
      return (idx_ == other.idx_ && grid_ == other.grid_);
    }

    [[nodiscard]]
    bool operator!=(const CellGridIterator &other) const {
      return !(*this == other);
    }

  private:
    const CellGrid *grid_;
    size_t idx_;
  };

  static CellGridBuilder create(std::shared_ptr<const CellOctree> octree);

  [[nodiscard]]
  size_t size() const {
    return mortonIndices_.size();
  }

  [[nodiscard]]
  const CellOctree &octree() const {
    return *octree_;
  }

  [[nodiscard]]
  std::span<const MortonIndex> mortonIndices() const {
    return mortonIndices_;
  }

  [[nodiscard]]
  std::size_t getEnumerationIndex(std::size_t streamIndex) const {
    return streamIndexToEnum_[streamIndex];
  }

  [[nodiscard]]
  std::size_t getEnumerationIndex(const CellOctree::CellView &cv) const {
    if (cv.isPhantom()) {
      return NOT_ENUMERATED;
    }
    return getEnumerationIndex(cv.streamIndex());
  }

  [[nodiscard]]
  AdjacencyListView neighborIndices(const AdjacencyOffset offset) const {

    for (size_t i = 0; i < adjacencyOffsets_.size(); ++i) {
      if (adjacencyOffsets_.at(i) == offset) {
        return adjacencyLists_.at(i);
      }
    }

    const auto sourceInfo = source_info(std::source_location::current());
    throw std::out_of_range(
        std::format("{}: Invalid adjacency offset for octree.", sourceInfo));
  }

  [[nodiscard]]
  CellOctree::CellView operator[](size_t idx) const {

    if (idx >= mortonIndices_.size()) {
      const auto srcInfo = source_info(std::source_location::current());
      throw std::out_of_range(
          std::format("{}: Cell index was out of range.", srcInfo));
    }
    // NOLINTNEXTLINE (bugprone-unchecked-optional-access)
    return octree_->getCell(mortonIndices_.at(idx)).value();
  }

  [[nodiscard]]
  bool operator==(const CellGrid &other) const {
    return octree_ == other.octree_ && mortonIndices_ == other.mortonIndices_;
  }

  // -------------------------------------------------------------------------
  // Range Interface
  // -------------------------------------------------------------------------
  [[nodiscard]]
  CellGridIterator begin() const {
    return {this, 0};
  };

  [[nodiscard]]
  CellGridIterator begin() {
    return {this, 0};
  };

  [[nodiscard]]
  CellGridIterator end() const {
    return {this, size()};
  };

  [[nodiscard]]
  CellGridIterator end() {
    return {this, size()};
  };

private:
  friend class CellGridBuilder;
  std::shared_ptr<const CellOctree> octree_;
  std::vector<MortonIndex> mortonIndices_;
  std::vector<size_t> streamIndexToEnum_;
  std::vector<double> scalarFields_;
  std::vector<Vec3D> vectorFields_;

  // We do not use hash map because the set of offsets is pretty small
  std::vector<AdjacencyOffset> adjacencyOffsets_;
  std::vector<AdjacencyList> adjacencyLists_;

  CellGrid(std::shared_ptr<const CellOctree> octree,
           std::vector<MortonIndex> mortonIndices,
           std::vector<size_t> streamIndexMapping,
           std::vector<AdjacencyOffset> offsets,
           std::vector<AdjacencyList> neighborLists)
      : octree_(std::move(octree)), mortonIndices_(std::move(mortonIndices)),
        streamIndexToEnum_(std::move(streamIndexMapping)),
        adjacencyOffsets_(std::move(offsets)),
        adjacencyLists_(std::move(neighborLists)) {}
};

// -------------------------------------------------------------------------
// Periodicity Mapper
// -------------------------------------------------------------------------
class PeriodicityMapper {
public:
  PeriodicityMapper() = default;
  virtual ~PeriodicityMapper() = default;
  PeriodicityMapper(const PeriodicityMapper &) = default;
  PeriodicityMapper(PeriodicityMapper &&) = default;
  PeriodicityMapper &operator=(const PeriodicityMapper &) = default;
  PeriodicityMapper &operator=(PeriodicityMapper &&) = default;

  [[nodiscard]]
  virtual SignedGridCoordinates
  getNeighborCoordinates(SignedGridCoordinates goalCoords,
                         size_t lvl) const = 0;
};

class NoPeriodicity : public PeriodicityMapper {
public:
  [[nodiscard]]
  SignedGridCoordinates getNeighborCoordinates(SignedGridCoordinates goalCoords,
                                               size_t lvl) const override;
};

class Torus : public PeriodicityMapper {
public:
  Torus(bool x_periodic, bool y_periodic, bool z_periodic)
      : periodic_{x_periodic, y_periodic, z_periodic} {}

  [[nodiscard]]
  SignedGridCoordinates getNeighborCoordinates(SignedGridCoordinates goalCoords,
                                               size_t lvl) const override;

private:
  std::array<bool, 3> periodic_;
};

// -------------------------------------------------------------------------
// Cell Grid Builder
// -------------------------------------------------------------------------
class CellGridBuilder {
public:
  CellGridBuilder(std::shared_ptr<const CellOctree> octree);

  [[nodiscard]]
  CellGridBuilder &levels(std::span<const std::size_t> lvls);

  [[nodiscard]]
  CellGridBuilder &levels(std::initializer_list<const std::size_t> lvls);

  [[nodiscard]]
  CellGridBuilder &neighborhood(std::span<const AdjacencyOffset> offsets);

  [[nodiscard]]
  CellGridBuilder &
  neighborhood(std::initializer_list<const AdjacencyOffset> offsets);

  template <typename T>
  [[nodiscard]]
  CellGridBuilder &periodicityMapper(T periodicityHandler) {
    static_assert(std::is_base_of_v<PeriodicityMapper, T>,
                  "T must derive from PeriodicityMapper!");
    periodicityHandler_ = std::make_unique<T>(std::move(periodicityHandler));
    return *this;
  }

  [[nodiscard]] CellGrid build();

private:
  std::shared_ptr<const CellOctree> octree_;
  std::vector<std::size_t> levels_;
  std::vector<AdjacencyOffset> adjacencyOffsets_;
  std::unique_ptr<PeriodicityMapper> periodicityHandler_;
};

} // namespace oktal