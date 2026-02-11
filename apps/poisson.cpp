
#include "oktal/io/VtkExport.hpp"
#include "oktal/octree/CellGrid.hpp"
#include "oktal/octree/CellOctree.hpp"
#include <cstdlib>
#include <iostream>
#include <span>
#include <stdexcept>

using oktal::CellOctree, oktal::CellGrid;
using namespace oktal;

const std::vector<Vec<std::ptrdiff_t, 3>> neighboorhood = {
    {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
namespace {

void PrintUsage(const char *name) {
  std::cerr << "Usage: " << name
            << "<refinementLevel> <max-iterations> <epsilon> <output-file> \n";
}

double eval_phi(const Vec3D &pos) {
  return std::cos(M_PI * pos[0]) * std::cos(M_PI * pos[1]) *
         std::cos(M_PI * pos[2]);
}

void initialise(const CellGrid &cells, std::vector<double> &u,
                std::vector<double> &f) {
  for (auto cell : cells) {
    const Vec3D center = cell.center();
    f[cell] = 3 * M_PI * M_PI * eval_phi(center);
    u[cell] = eval_phi(center);
  }
}

double compute_residual_norm(const CellGrid &cells,
                             const std::vector<double> &u,
                             const std::vector<double> &f,
                             std::vector<double> &residual, const double h) {
  const size_t num_cells = cells.size();
  double sum_res = 0.0;
  for (auto cell : cells) {
    double neighbor_sum = 0.0;
    for (const auto &offset : neighboorhood) {
      const size_t nbIdx = cells.neighborIndices(offset)[cell];
      if (nbIdx == CellGrid::NO_NEIGHBOR) {
        break;
      }
      neighbor_sum += u[nbIdx];
    }
    if (cell.isValid()) {
      residual[cell] = f[cell] + ((-6.0 * u[cell] + neighbor_sum) / (h * h));
    }
    sum_res += residual[cell] * residual[cell];
  }
  const double norm = std::sqrt(sum_res) / static_cast<double>(num_cells);
  return norm;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void solve_poisson(const size_t refinement_level, const double epsilon,
                   const unsigned int max_iters, std::string &output_file) {
  auto cell_octree = CellOctree::createUniformGrid(refinement_level);
  auto cell_grid =
      CellGrid::create(cell_octree).neighborhood(neighboorhood).build();
  const size_t numCells = cell_grid.size();
  const double h = cell_octree->geometry().dx(refinement_level);
  std::vector<double> f(numCells);
  std::vector<double> u(numCells, 0.0);
  std::vector<double> u_tmp(numCells, 0.0);
  std::vector<double> residual(numCells, 0.0);

  initialise(cell_grid, u, f);

  int iter = 0;
  double l2_norm = 100;

  for (unsigned int i = 0; i < max_iters; ++i) {
    for (auto cell : cell_grid) {
      bool isBoundary = false;
      double neighborSum = 0;

      for (const auto &offset : neighboorhood) {
        const size_t nbIdx = cell_grid.neighborIndices(offset)[cell];
        if (nbIdx == CellGrid::NO_NEIGHBOR) {
          isBoundary = true;
          break;
        }
        neighborSum += u[nbIdx];
      }
      if (!isBoundary) {
        u_tmp[cell] = (h * h * f[cell] + neighborSum) / 6.0;
      } else {
        u_tmp[cell] = u[cell];
      }
    }
    u.swap(u_tmp);
    l2_norm = compute_residual_norm(cell_grid, u, f, residual, h);
    if (l2_norm < epsilon) {
      break;
    }
    iter++;
  }
  std::cout << "L2 residual norm : " << l2_norm
            << " and numbers of iterations required " << iter << "\n";
  io::vtk::exportCellGrid(cell_grid, output_file)
      .writeGridVector("u", u)
      .writeGridVector("f", f)
      .writeGridVector("residual", residual);
}
} // namespace
int main(int argc, char **argv) {
  /*   try{
         std::span<char*> args(argc, argv);
 }*/
  const std::span<char *> args(argv, static_cast<size_t>(argc));
  if (args.size() != 5) {
    std::cerr << "Error : Invalid number of Arguments ! Expected 4 \n";
    PrintUsage(args[0]);
    exit(EXIT_FAILURE);
  }

  // Parse arguments
  try {
    char *endPtr = nullptr;
    const long refinement_level_ = std::strtol(args[1], &endPtr, 10);
    if (*endPtr != '\0' || refinement_level_ < 0) {
      throw std::invalid_argument(
          "refinementLevel must be a non-negative integer \n");
    }

    const auto refinement_level = static_cast<size_t>(refinement_level_);
    const long max_iter_ = std::strtol(args[2], &endPtr, 10);
    if (*endPtr != '\0' || max_iter_ < 0) {
      throw std::invalid_argument(
          "Max_iteration must be a non-negative integer \n");
    }
    const auto max_iter = static_cast<unsigned int>(max_iter_);

    const double eps = std::strtod(args[3], &endPtr);
    if (*endPtr != '\0' || eps <= 0.0) {
      throw std::invalid_argument("epsilon must be a non-negative number \n");
    }

    std::string outputFile = args[4];
    if (outputFile.empty()) {
      throw std::invalid_argument("output file cannot be empty.\n");
    }
    solve_poisson(refinement_level, eps, max_iter, outputFile);

  } catch (const std::exception &e) {
    std::cerr << "Error parsing arguments : " << e.what() << "\n";
    PrintUsage(args[0]);
    exit(EXIT_FAILURE);
  } catch (...) {
    std::cerr << "Unknown Error \n";
    std::exit(EXIT_FAILURE);
  }

  return 0;
}