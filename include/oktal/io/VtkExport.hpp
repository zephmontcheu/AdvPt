#pragma once

#include "advpt/htgfile/VtkHtgFile.hpp"
#include "oktal/octree/CellGrid.hpp"
#include "oktal/octree/CellOctree.hpp"

#include <filesystem>

using AccessMode = decltype(HighFive::File::ReadOnly);

template <typename F>
concept ConstructsHighFiveFile = requires(const F &f, const AccessMode &a) {
  { HighFive::File(f, a) } -> std::same_as<HighFive::File>;
};

template <ConstructsHighFiveFile F>
inline HighFive::File makeH5File(const F &filename,
                                 const AccessMode &accessMode) {
  return HighFive::File(filename, accessMode);
}

template <typename F>
inline HighFive::File makeH5File(const F &filename,
                                 const AccessMode &accessMode)
  requires(!ConstructsHighFiveFile<F> &&
           ConstructsHighFiveFile<decltype(filename.string())>)
{
  return HighFive::File(filename.string(), accessMode);
}

namespace oktal::io::vtk {
void exportOctree(const CellOctree &octree,
                  const std::filesystem::path &filepath);

class CellGridExporter {
  CellGrid const *pGrid;
  advpt::htgfile::SnapshotHtgFile file;

  CellGridExporter(const CellGrid &grid, advpt::htgfile::SnapshotHtgFile file_)
      : pGrid(&grid), file(std::move(file_)) {}

  friend CellGridExporter exportCellGrid(const CellGrid &grid,
                                         const std::filesystem::path &filepath);

public:
  template <typename T>
  CellGridExporter &writeGridVector(const std::string &name, std::vector<T> data) {
    const auto &nodes = pGrid->octree().nodesStream();
    if (data.size() < nodes.size()) {
      const auto diff = nodes.size() - data.size();
      data.insert(data.cbegin(), diff, {0});
    }
    file.writeCellData<T>(name, data);
    return *this;
  }
};

[[nodiscard]] CellGridExporter
exportCellGrid(const CellGrid &grid, const std::filesystem::path &filepath);
}; // namespace oktal::io::vtk