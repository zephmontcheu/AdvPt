#include "oktal/io/VtkExport.hpp"
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <queue>
#include <tuple>

using advpt::htgfile::HyperTree;
using advpt::htgfile::SnapshotHtgFile;
using oktal::CellOctree;
using oktal::morton_bits_t;

namespace {
bool allArePhantoms(const CellOctree &octree, const CellOctree::Node &node) {
  if (!node.isPhantom()) {
    return false;
  }

  if (node.isRefined()) {
    for (const auto &child :
         octree.nodesStream().subspan(node.childrenStartIndex(), 8)) {
      if (!allArePhantoms(octree, child)) {
        return false;
      }
    }
  }

  return true;
}

auto createHtgFile(const CellOctree &octree,
                   const std::filesystem::path &filepath) {
  const auto &levelNum = octree.numberOfLevels();
  const auto &nodesStream = octree.nodesStream();
  const auto &levelInfo = octree.getLevels();
  const auto &geometry = octree.geometry();

  const auto &topLevelMin = geometry.cellMinCorner({});
  const auto &topLevelMax = geometry.cellMaxCorner({});

  const auto totalCount = nodesStream.size();
  const auto notFinestCount = totalCount - levelInfo[levelNum - 1].second;
  std::vector<morton_bits_t> levels;

  HyperTree hyperTree;
  hyperTree.xCoords = {topLevelMin[0], topLevelMax[0]};
  hyperTree.yCoords = {topLevelMin[1], topLevelMax[1]};
  hyperTree.zCoords = {topLevelMin[2], topLevelMax[2]};
  auto &nodesPerDepth = hyperTree.nodesPerDepth;
  auto &descriptor = hyperTree.descriptor;
  auto &mask = hyperTree.mask;

  /* if (levelNum == 0) {
    return SnapshotHtgFile::create(filepath, hyperTree);
  } */

  levels.reserve(totalCount);
  nodesPerDepth.reserve(levelNum);
  descriptor.resize(static_cast<size_t>(
      std::ceil(static_cast<double>(notFinestCount) / 8.0)));
  mask.resize(
      static_cast<size_t>(std::ceil(static_cast<double>(totalCount) / 8.0)));

  std::ranges::transform(
      levelInfo, std::insert_iterator(nodesPerDepth, nodesPerDepth.end()),
      std::identity(), &std::decay_t<decltype(levelInfo)>::element_type::first);

  size_t orderIndex = 0;
  size_t level = 0;
  size_t levelEndIndex = levelInfo[0].second;
  for (const auto &node : nodesStream) {
    if (orderIndex >= levelEndIndex) {
      ++level;
      levelEndIndex += levelInfo[level].second;
    }

    levels.emplace_back(level);

    if (node.isRefined() || node.isPhantom()) {
      const auto descriptorIndex = orderIndex >> 3;
      const auto bitIndex = 7 - (orderIndex & 7);

      if (node.isRefined()) {
        descriptor.at(descriptorIndex) |= (1 << bitIndex);
      }

      if (allArePhantoms(octree, node)) {
        mask.at(descriptorIndex) |= (1 << bitIndex);
      }
    }

    ++orderIndex;
  }

  return SnapshotHtgFile::create(filepath, hyperTree)
      .writeCellData("level", levels);
}
} // namespace

namespace oktal::io::vtk {
void exportOctree(const CellOctree &octree,
                  const std::filesystem::path &filepath) {
  /* static const auto& startIndices = MortonIndex::startIndices(); */
  createHtgFile(octree, filepath);
}

[[nodiscard]] CellGridExporter
exportCellGrid(const CellGrid &grid, const std::filesystem::path &filepath) {
  return {grid, createHtgFile(grid.octree(), filepath)};
}
}; // namespace oktal::io::vtk