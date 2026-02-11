#include "oktal/octree/CellOctree.hpp"
#include "oktal/octree/MortonIndex.hpp"
#include "oktal/octree/OctreeGeometry.hpp"
#include <algorithm>
#include <cstddef>
#include <optional>
#include <ranges>
#include <stdexcept>

namespace {
bool isInvalidDescriptor(const std::string_view &descriptor) {
  const std::string_view allowedChars = ".RPX";

  size_t numRefinedNodes = 0;
  size_t numTotalNodes = 0;
  for (const char c : descriptor) {

    // Skip pipes, as these are not nodes
    if (c == '|') {
      continue;
    }

    // If char not allowed -> invalid
    if (allowedChars.find(c) == std::string::npos) {
      return true;
    }

    if (c == 'X' || c == 'R') {
      numRefinedNodes++;
    }

    numTotalNodes++;
  }

  return numTotalNodes != (numRefinedNodes * 8 + 1);
}
} // anonymous namespace

namespace oktal {

CellOctree CellOctree::fromDescriptor(std::string_view descriptor) {

  if (descriptor.empty() || isInvalidDescriptor(descriptor)) {
    throw std::invalid_argument("Invalid descriptor was passed!");
  }

  CellOctree tree;
  std::vector<std::size_t> refinedNodes;
  tree.nodesStream_.pop_back(); // Remove the dummy root node
  tree.levels_.back().second--; // Also remove dummy root from counting
  for (const char c : descriptor) {
    if (c == '.' || c == 'R' || c == 'P' || c == 'X') {
      const bool refined = (c == 'R' || c == 'X');
      const bool phantom = (c == 'P' || c == 'X');
      tree.nodesStream_.emplace_back(refined, phantom);
      tree.levels_.back().second++; // increment node count for current level
      if (refined) {
        refinedNodes.emplace_back(tree.numberOfNodes() - 1);
      }

    } else if (c == '|') {

      // Start new level with all children of refined nodes from previous level
      tree.levels_.emplace_back(tree.numberOfNodes(), 0);

      size_t refinedIdx = 0;
      // Calculate childrenStartIndex (global index) for each refined node
      for (const std::size_t nodeIdx : refinedNodes) {
        Node &parentOfRefinement = tree.nodesStream_.at(nodeIdx);
        parentOfRefinement.setChildrenStartIndex(tree.numberOfNodes() - 1 +
                                                 (8 * refinedIdx++ + 1));
      }
      refinedNodes.clear();
    } else {
      throw std::invalid_argument("Invalid descriptor was passed!");
    }
  }

  return tree;
}

[[nodiscard]] std::shared_ptr<const CellOctree>
CellOctree::createUniformGrid(OctreeGeometry geom, size_t level) {
  decltype(levels_) levels;
  decltype(nodesStream_) nodes;

  levels.reserve(level + 1);
  levels.emplace_back(0, 1);
  while (levels.size() <= level) {
    auto [start, size] = levels.back();
    levels.emplace_back(start + size, size << 3);
  }

  auto [start, size] = levels.back();
  const size_t nodes_count = start + size;
  nodes.reserve(nodes_count);
  for (const auto &index : std::views::iota(size_t(1), nodes_count) |
                               std::views::stride(size_t(8))) {
    nodes.emplace_back(true, true, index);
  }
  nodes.insert(nodes.cend(), size, {});

  return std::make_shared<const CellOctree>(std::move(nodes), std::move(levels),
                                            geom);
}

[[nodiscard]]
std::optional<CellOctree::CellView>
CellOctree::getCell(const MortonIndex &m) const {
  // Root case
  if (m.isRoot()) {
    const Node &root = nodesStream_.at(0);
    if (root.isPhantom()) {
      return std::nullopt;
    }
    return CellView{root, geometry(), m, 0};
  }

  // If requested level is not present, the cell cannot exist
  if (m.level() >= numberOfLevels()) {
    return std::nullopt;
  }

  // Traverse from root following the Morton path
  std::size_t currentIdx = 0;
  const Node *current = &nodesStream_.at(0);

  const auto path = m.getPath();
  for (const auto &choice : path) {
    if (!current->isRefined()) {
      return std::nullopt;
    }
    currentIdx = current->childIndex(static_cast<std::size_t>(choice));
    current = &nodesStream_.at(currentIdx);
  }

  if (current->isPhantom()) {
    return std::nullopt;
  }
  return CellView{*current, geometry(), m, currentIdx};
}

[[nodiscard]]
bool CellOctree::cellExists(const MortonIndex &m) const {

  return getCell(m).has_value();
}

std::optional<CellOctree::CellView> CellOctree::getRootCell() const {
  const Node &root = nodesStream_.at(0);
  if (root.isPhantom()) {
    return std::nullopt;
  }

  return CellView{root, geometry_, MortonIndex(), 0};
}

// NOLINTNEXTLINE
OctreeCellsRange<DfsPolicy> CellOctree::preOrderDepthFirstRange() const {
  OctreeCursor end(*this);
  end.toEnd();
  return {{*this}, end, DfsPolicy{}};
}

// NOLINTNEXTLINE
void DfsPolicy::advance(OctreeCursor &cursor) const {

  while (true) {
    // If we are at the end or the cursor is invalid, we stop.
    if (cursor.empty() || cursor.end()) {
      return;
    }

    // If X or R, we have children -> descend to first child
    if (cursor.currentNode().isRefined()) {
      cursor.descend();
    }
    // If we do not have children, we check for the next sibling
    else if (!cursor.lastSibling()) {
      cursor.nextSibling();
    }
    // If we do not have children or siblings -> ascend and go to next sibling
    // of parent
    else {
      while (!cursor.empty() && !cursor.end()) {
        cursor.ascend(); // go to parent node
        // If root, we are finished
        if (cursor.end()) {
          return;
        }
        // If the parent has siblings, we take the next sibling, else we ascend
        // further
        if (!cursor.lastSibling()) {
          cursor.nextSibling();
          break;
        }
      }
    }

    // If the new node is not a phantom, we stop
    // If a phantom, we ignore it and continue
    if (!cursor.currentNode().isPhantom()) {
      return;
    }
  }
}

OctreeCellsRange<HorizontalPolicy>
CellOctree::horizontalRange(const std::size_t level) const {

  // If invalid level, return empty range
  if (level >= numberOfLevels()) {
    OctreeCursor end(*this);
    end.toEnd();
    return {end, end, HorizontalPolicy{}};
  }

  // Get level information
  const auto &[startIndex, size] = levels_.at(level);

  // Initialize start cursor with correct stream index
  // This should be done to keep a valid starting state
  std::vector<size_t> startPath(level + 1);
  // startPath.back() = startIndex;
  OctreeCursor start(*this, startPath);
  start.updatePath(
      startIndex); // update to the full path from the root to start index

  OctreeCursor end = {*this};
  end.toEnd();

  return {std::move(start), end, HorizontalPolicy{}};
}

// NOLINTNEXTLINE
void HorizontalPolicy::advance(OctreeCursor &cursor) const {

  if (cursor.empty() || cursor.end()) {
    return;
  }

  // Determine the original group of siblings
  const auto initialIndex = cursor.currentStreamIndex();
  const auto initialGroupIndex = (initialIndex - 1) >> 3;

  // Continue advancing the stream index until a non-phantom node is found
  // or the end of the range is reached (checked by the iterator)
  while (true) {
    cursor.advanceStreamIndex(); // Increment stream index (always, not only if
                                 // it is not last sibling)

    if (cursor.end()) {
      return;
    }

    // Stop if the new node is not a phantom.
    if (!cursor.currentNode().isPhantom()) {
      const auto nextStreamIndex = cursor.currentStreamIndex();
      const auto currentGroupIndex = (nextStreamIndex - 1) >> 3;

      // Update if initial and current node are not siblings
      if (initialGroupIndex != currentGroupIndex) {
        cursor.updatePath(nextStreamIndex);
      }

      return;
    }
  }
}

} // namespace oktal