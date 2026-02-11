#pragma once
#include "oktal/geometry/Box.hpp"
#include "oktal/geometry/Vec.hpp"
#include "oktal/octree/MortonIndex.hpp"
#include "oktal/octree/OctreeGeometry.hpp"
#include <concepts>
#include <cstddef>
#include <exception>
#include <format>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

inline std::string source_info(const std::source_location &location) {
  return std::format("[File: {}, Function: {}]", location.file_name(),
                     location.function_name());
}

namespace oktal {
class OctreeCursor; // Forward declaration
template <class T>
concept OctreeIteratorPolicy =
    std::semiregular<T> && requires(const T &policy, OctreeCursor &cursor) {
      { policy.advance(cursor) } -> std::same_as<void>;
    };

template <OctreeIteratorPolicy TPolicy> class OctreeCellsRange;

// ---------------------- Depth First Search Policy ----------------------

class DfsPolicy {
public:
  DfsPolicy() = default;

  //  NOLINTNEXTLINE
  void advance(OctreeCursor &cursor) const;
};

class HorizontalPolicy {
public:
  HorizontalPolicy() = default;

  // NOLINTNEXTLINE
  void advance(OctreeCursor &cursor) const;
};

// ---------------------- Main Octree Class ----------------------

class CellOctree {
public:
  // ----- Nested Node class -----
  class Node {
  public:
    Node() : data_(0) {};
    Node(bool refined, bool phantom, std::size_t childrenIdx = 0)
        : data_((childrenIdx & CHILDREN_INDEX_MASK) | (size_t(refined) << 62) |
                (size_t(phantom) << 63)) {};
    [[nodiscard]] bool isRefined() const {
      return static_cast<bool>(data_ & REFINED_BIT_MASK);
    }
    [[nodiscard]] bool isPhantom() const {
      return static_cast<bool>(data_ & PHANTOM_BIT_MASK);
    }
    void setRefined(bool v) {
      data_ = v ? (data_ | REFINED_BIT_MASK) : (data_ & ~REFINED_BIT_MASK);
    }
    void setPhantom(bool v) {
      data_ = v ? (data_ | PHANTOM_BIT_MASK) : (data_ & ~PHANTOM_BIT_MASK);
    }
    [[nodiscard]] std::size_t childrenStartIndex() const {
      return data_ & CHILDREN_INDEX_MASK;
    }
    /**
     * @brief Set the Children Start Index object
     * @details First clears the old index bits (& ~MASK) and then assigns new
     * index via &
     * @param size
     */
    void setChildrenStartIndex(std::size_t size) {
      data_ = (data_ & ~CHILDREN_INDEX_MASK) | (size & CHILDREN_INDEX_MASK);
    }
    [[nodiscard]] std::size_t childIndex(std::size_t branch) const {
      return childrenStartIndex() + branch;
    }

  private:
    std::size_t data_;
    static constexpr std::size_t REFINED_BIT_MASK =
        std::size_t{1} << 62; // 0100 0000 ... (bit 63 is set)
    static constexpr std::size_t PHANTOM_BIT_MASK =
        std::size_t{1} << 63; // 1000 0000 ... (bit 64 is set)
    static constexpr std::size_t CHILDREN_INDEX_MASK =
        ~(PHANTOM_BIT_MASK | REFINED_BIT_MASK); // 0011 1111 ... (all bits
                                                // except the two MSB are set)
  };

  //--------Nested Class Cellview -------

  class CellView {

  public:
    CellView(const Node &m_node_, const OctreeGeometry &m_geometry_,
             const MortonIndex &m_, std::size_t m_streamIndex)
        : streamIndex_(m_streamIndex),
          node_(m_node_.isRefined(), m_node_.isPhantom(),
                m_node_.childrenStartIndex()),
          m_geometry(m_geometry_.origin(), m_geometry_.sidelength()),
          m(m_.getBits()) {}

    [[nodiscard]]
    const MortonIndex &mortonIndex() const {
      return m;
    }

    [[nodiscard]]
    bool isRoot() const {
      return m.isRoot();
    }

    [[nodiscard]]
    bool isRefined() const {
      return node_.isRefined();
    }

    [[nodiscard]]
    std::size_t level() const {
      return m.level();
    }

    [[nodiscard]]
    std::size_t streamIndex() const {
      return streamIndex_;
    }

    [[nodiscard]]
    Vec3D center() const {
      return m_geometry.cellCenter(m);
    }

    [[nodiscard]]
    Box<double> boundingBox() const {
      return m_geometry.cellBoundingBox(m);
    }

    [[nodiscard]]
    bool isPhantom() const {
      return node_.isPhantom();
    }

  private:
    std::size_t streamIndex_;
    Node node_;
    OctreeGeometry m_geometry;
    MortonIndex m;
  };

  // NOLINTNEXTLINE
  OctreeCellsRange<DfsPolicy> preOrderDepthFirstRange() const;

  // NOLINTNEXTLINE
  OctreeCellsRange<HorizontalPolicy> horizontalRange(std::size_t level) const;

private:
  // ----- Members for CellOctree -----
  std::vector<Node> nodesStream_;
  using levelStartIndex = std::size_t;
  using levelSize = std::size_t;
  std::vector<std::pair<levelStartIndex, levelSize>> levels_;
  OctreeGeometry geometry_;

public:
  CellOctree() : nodesStream_({{}}), levels_({{0, 1}}) {}

  explicit CellOctree(const OctreeGeometry &m_geometry)
      : nodesStream_({{}}), levels_({{0, 1}}), geometry_(m_geometry) {}

  CellOctree(decltype(nodesStream_) &&nodesStream, decltype(levels_) &&levels,
             const decltype(geometry_) &geometry)
      : nodesStream_(std::move(nodesStream)), levels_(std::move(levels)),
        geometry_(geometry) {}

  [[nodiscard]] std::size_t numberOfNodes() const {
    return nodesStream_.size();
  }
  [[nodiscard]] std::size_t numberOfLevels() const { return levels_.size(); }
  [[nodiscard]] std::size_t numberOfNodes(std::size_t level) const {
    if (level >= levels_.size()) {
      return 0;
    }
    return levels_.at(level).second;
  }

  [[nodiscard]] std::size_t numberOfNonPhantomNodes(std::size_t level) const {
    if (level >= levels_.size()) {
      return 0;
    }
    const auto [levelStartIdxUnsigned, levelSize] = levels_.at(level);
    const auto levelStartIdx =
        static_cast<std::ptrdiff_t>(levelStartIdxUnsigned);
    const auto levelEndIdx =
        static_cast<std::ptrdiff_t>(levelStartIdxUnsigned + levelSize);
    return std::accumulate(nodesStream_.begin() + levelStartIdx,
                           nodesStream_.begin() + levelEndIdx, size_t{0},
                           [&](size_t acc, const CellOctree::Node node) {
                             return acc +
                                    static_cast<size_t>(!node.isPhantom());
                           });
  }

  [[nodiscard]]
  std::size_t numberOfNonPhantomNodes(std::vector<size_t> &listOfLevels) const {
    size_t sum = 0;
    for (const size_t lvlNum : listOfLevels) {
      sum += numberOfNonPhantomNodes(lvlNum);
    }
    return sum;
  }

  [[nodiscard]]
  std::size_t numberOfNonPhantomNodes() const {
    size_t sum = 0;
    for (size_t i = 0; i < levels_.size(); ++i) {
      sum += numberOfNonPhantomNodes(i);
    }
    return sum;
  }

  [[nodiscard]] std::span<const Node> nodesStream() const {
    return {nodesStream_};
  }
  [[nodiscard]] std::span<const Node> nodesStream(std::size_t level) const {

    if (level >= levels_.size()) {
      return std::span<const Node>{};
    }
    return nodesStream().subspan(levels_.at(level).first,
                                 levels_.at(level).second);
  }
  [[nodiscard]] std::span<const decltype(levels_)::value_type>
  getLevels() const {
    return {levels_};
  }

  static CellOctree fromDescriptor(std::string_view descriptor);

  [[nodiscard]]
  std::optional<CellView> getCell(const MortonIndex &m) const;

  [[nodiscard]]
  bool cellExists(const MortonIndex &m) const;

  [[nodiscard]]
  std::optional<CellView> getRootCell() const;
  // Add the getter-accessor to geometry_
  [[nodiscard]]
  const OctreeGeometry &geometry() const {
    return geometry_;
  }

  [[nodiscard]] static std::shared_ptr<const CellOctree>
  createUniformGrid(OctreeGeometry geom, size_t level);

  [[nodiscard]] static std::shared_ptr<const CellOctree>
  createUniformGrid(size_t level) {
    return createUniformGrid({}, level);
  }
};

// ---------------------- Octree Cursor ----------------------

class OctreeCursor {
  const CellOctree *pOctree;
  std::vector<size_t> vPath;

  [[nodiscard]] decltype(auto) getNode(this auto &&c, const size_t &index) {
    return c.pOctree->nodesStream()[index];
  }

  [[nodiscard]] decltype(auto) currentParent(this auto &&c) {
    return c.getNode(c.vPath.at(c.vPath.size() - 2));
  }

public:
  using path_view = std::span<const decltype(vPath)::value_type>;

  OctreeCursor() : pOctree(nullptr) {}
  OctreeCursor(const CellOctree &octree_) : pOctree(&octree_), vPath{0} {}
  OctreeCursor(const CellOctree &octree_, const path_view &path_)
      : pOctree(&octree_), vPath(path_.cbegin(), path_.cend()) {}

  [[nodiscard]] const CellOctree *octree() const { return pOctree; }

  [[nodiscard]] path_view path() const { return {vPath}; }

  [[nodiscard]] bool empty() const noexcept {
    return !static_cast<bool>(pOctree);
  }

  [[nodiscard]] decltype(auto) currentNode(this auto &&c) {
    static const auto info = source_info(std::source_location::current());
    if (c.end()) {
      throw std::logic_error(
          std::format("{}: No current Node, vPath is empty", info));
    }
    const auto myLevel = c.currentLevel();
    const auto &levels = c.pOctree->getLevels();
    if (myLevel >= levels.size()) {
      throw std::out_of_range(
          std::format("{}: Current level({}) exceeds max level({})", info,
                      myLevel, levels.size() - 1));
    }
    const auto [levelStart, levelSize] = levels[myLevel];
    const auto myIndex = c.vPath.at(myLevel) - levelStart;
    if (myIndex >= levelSize) {
      throw std::out_of_range(
          std::format("{}: Current index({}) equal or greater than max "
                      "index({}) of level({})",
                      info, myIndex, levelSize, myLevel));
    }
    return c.getNode(c.currentStreamIndex());
  }

  [[nodiscard]] bool end() const noexcept { return vPath.empty(); }

  [[nodiscard]] size_t currentLevel() const noexcept {
    return vPath.size() - 1;
  }

  [[nodiscard]] const size_t &currentStreamIndex() const noexcept {
    return vPath.back();
  }

  [[nodiscard]] std::optional<CellOctree::CellView>
  currentCell() const noexcept {
    try {
      return pOctree->getCell(mortonIndex());
    } catch (...) {
      return std::nullopt;
    }
  }

  [[nodiscard]] bool firstSibling() const noexcept {
    if (vPath.size() > 1) {
      return (currentStreamIndex() & 7) == 1;
    }
    return true;
  }

  [[nodiscard]] bool lastSibling() const noexcept {
    if (vPath.size() > 1) {
      return (currentStreamIndex() & 7) == 0;
    }
    return true;
  }

  [[nodiscard]] size_t siblingIndex() const noexcept {
    if (vPath.size() > 1) {
      return currentStreamIndex() - 1;
    }
    return 0uz;
  }

  [[nodiscard]] MortonIndex mortonIndex() const
      noexcept(std::is_nothrow_convertible_v<size_t, morton_bits_t>) {
    auto bits = morton_bits_t(1);
    for (const size_t &index : std::span{vPath.begin() + 1, vPath.end()}) {
      bits = bits << 3;
      bits |= (static_cast<morton_bits_t>(index) - 1) & 7;
    }
    return {bits};
  }

  [[nodiscard]] bool operator==(const OctreeCursor &other) const noexcept {
    return vPath.size() == other.vPath.size() && pOctree == other.pOctree &&
           (end() || vPath.back() == other.vPath.back());
  }

  [[nodiscard]] bool operator!=(const OctreeCursor &other) const noexcept {
    return vPath.size() != other.vPath.size() || pOctree != other.pOctree ||
           (!end() && vPath.back() != other.vPath.back());
  }

  void ascend() {
    if (!end()) [[likely]] {
      vPath.resize(vPath.size() - 1);
    }
  }

  void descend() {
    if (!end()) [[likely]] {
      const auto &node = currentNode();
      if (node.isRefined()) {
        vPath.emplace_back(node.childrenStartIndex());
      }
    }
  }

  void descend(const size_t &childIdx) {
    static const auto info = source_info(std::source_location::current());
    if (childIdx >= 8) [[unlikely]] {
      throw std::out_of_range(std::format(
          "{}: Child Index {} exceeds the range of 0 to 7", info, childIdx));
    }
    if (!end()) [[likely]] {
      const auto &node = currentNode();
      if (node.isRefined()) {
        vPath.emplace_back(node.childIndex(childIdx));
      }
    }
  }

  void previousSibling() {
    if (!firstSibling()) {
      --vPath.back();
    }
  }

  void nextSibling() {
    if (!lastSibling()) {
      ++vPath.back();
    }
  }

  void advanceStreamIndex() {
    if (!end()) {
      const auto &streamIndex = ++vPath.back();
      const auto [levelStart, levelSize] = pOctree->getLevels()[currentLevel()];
      if (streamIndex - levelStart >= levelSize) {
        toEnd();
        return;
      }
    }
  }

  // Updates the path for the given stream index
  void updatePath(const std::size_t &streamIndex) {
    static const auto info = source_info(std::source_location::current());
    if (!end()) {
      const auto myLevel = currentLevel();
      const auto &levels = pOctree->getLevels();

      if (myLevel >= levels.size()) {
        throw std::logic_error(
            std::format("{}: Current level({}) exceeds maximum level({})", info,
                        myLevel, levels.size() - 1));
      }

      const auto [levelStart, levelSize] = levels[myLevel];
      const auto myIndex = streamIndex - levelStart;

      if (myIndex >= levelSize) {
        throw std::out_of_range(
            std::format("{}: Current index({}) equal or greater than max "
                        "index({}) of level({})",
                        info, myIndex, levelSize, myLevel));
      }

      size_t currentStreamIndex = streamIndex;
      for (size_t l = myLevel; l >= 1; --l) {
        vPath.at(l) = currentStreamIndex;

        auto i = levels[l - 1].first;
        const auto &parentLevel = pOctree->nodesStream(l - 1);
        std::optional<size_t> parentIdx = std::nullopt;

        // Iterate through all nodes of the previous level and find the
        // corresponding parent
        for (const auto &node : parentLevel) {
          if (node.isRefined()) {
            const size_t childrenStart = node.childrenStartIndex();
            // Check if currentIdx is one of this parent's children
            if (currentStreamIndex >= childrenStart &&
                currentStreamIndex < childrenStart + 8) {
              parentIdx.emplace(i);
              break;
            }
          }
          ++i;
        }

        if (!parentIdx.has_value()) {
          throw std::logic_error(
              std::format("{}: No parent value found", info));
        }

        currentStreamIndex =
            *parentIdx; // Move up to the parent for the next iteration
      }
    }
  }

  void toSibling(const size_t &siblingIdx) {
    if (vPath.size() == size_t(1)) [[unlikely]] {
      if (static_cast<bool>(siblingIdx)) {
        throw std::out_of_range(std::format(
            "Nonzero Sibling Index {} not allowed with root node", siblingIdx));
      }
      return;
    } else if (siblingIdx >= 8) [[unlikely]] {
      throw std::out_of_range(std::format(
          "Sibling Index {} exceeds the range of 0 to 7", siblingIdx));
    }

    const auto &parent = currentParent();
    vPath.back() = parent.childIndex(siblingIdx);
  }

  void toEnd() { vPath.clear(); }
};

// ---------------------- Iterator Class ----------------------
template <OctreeIteratorPolicy TPolicy, bool isConst = false>
class OctreeIterator {
public:
  // Type aliases required by ForwardIterator concept
  using value_type = std::conditional_t<isConst, const CellOctree::CellView,
                                        CellOctree::CellView>;
  using iterator_category = std::forward_iterator_tag;
  using difference_type = std::ptrdiff_t;

  OctreeIterator() = default;
  OctreeIterator(OctreeCursor cursor, const TPolicy &policy)
      : pPolicy_(&policy), cursor_(std::move(cursor)) {
    // Iterator must be on the first non-phantom node
    while (!cursor_.empty() && !cursor_.end() &&
           cursor_.currentNode().isPhantom()) {
      pPolicy_->advance(cursor_);
    }
  }
  OctreeIterator(const OctreeIterator &other) = default;
  OctreeIterator &operator=(const OctreeIterator &other) = default;
  OctreeIterator(OctreeIterator &&other) = default;
  OctreeIterator &operator=(OctreeIterator &&other) = default;
  ~OctreeIterator() = default;

  [[nodiscard]]
  value_type operator*() const {
    // NOLINTNEXTLINE
    return cursor_.currentCell().value(); // MUST always be a valid value
  }
  OctreeIterator &operator++() {
    pPolicy_->advance(cursor_);
    return *this;
  }
  OctreeIterator operator++(int) {
    OctreeIterator tmp = *this;
    ++(*this);
    return tmp;
  }
  [[nodiscard]]
  bool operator==(const OctreeIterator &other) const {
    return this->cursor_ == other.cursor_;
  }
  [[nodiscard]]
  bool operator!=(const OctreeIterator &other) const {
    // return !(*this == other);
    return this->cursor_ != other.cursor_;
  }

private:
  const TPolicy *pPolicy_;
  OctreeCursor cursor_;
};

// ---------------------- Range Class ----------------------

template <OctreeIteratorPolicy TPolicy> class OctreeCellsRange {
public:
  OctreeCellsRange(OctreeCursor start, OctreeCursor end, TPolicy policy)
      : policy_(policy), startCursor_(std::move(start)),
        endCursor_(std::move(end)) {}
  OctreeCellsRange() = default;

  [[nodiscard]]
  OctreeIterator<TPolicy, false> begin() {
    return OctreeIterator<TPolicy, false>(startCursor_, policy_);
  }
  [[nodiscard]]
  OctreeIterator<TPolicy, false> end() {
    return OctreeIterator<TPolicy, false>(endCursor_, policy_);
  }

  [[nodiscard]]
  OctreeIterator<TPolicy, true> begin() const {
    return OctreeIterator<TPolicy, true>(startCursor_, policy_);
  }
  [[nodiscard]]
  OctreeIterator<TPolicy, true> end() const {
    return OctreeIterator<TPolicy, true>(endCursor_, policy_);
  }

private:
  TPolicy policy_;
  OctreeCursor startCursor_;
  OctreeCursor endCursor_;
};

} // namespace oktal