#include "oktal/io/VtkExport.hpp"
#include <format>
#include <iostream>
#include <span>
#include <stdexcept>

#include <cstdlib>

using oktal::CellOctree, oktal::io::vtk::exportOctree;

namespace {
inline void convert(const std::span<char *> &args) {
  if (args.size() == 3) {
    exportOctree(CellOctree::fromDescriptor(args[2]), args[1]);
  } else {
    throw std::invalid_argument(
        std::format("Invalid number of arguments: {}", args.size()));
  }
}
} // namespace

int main(int argc, char **argv) {
  try {
    convert({argv, std::span<char *>::size_type(argc)});
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  } catch (...) {
    std::cerr << "Unknown error";
    std::exit(EXIT_FAILURE);
  }

  std::exit(EXIT_FAILURE);
  
}