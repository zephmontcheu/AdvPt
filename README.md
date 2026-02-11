# Oktal - Adaptive Octree Spatial Data Structure Library

A modern C++23 library for efficient spatial data organization and processing using adaptive octrees. Oktal provides high-performance cell-based octree structures with support for geometry operations, spatial queries, and VTK visualization.

## Features

- **Adaptive Octree Data Structure**: Efficient spatial partitioning using octrees with customizable resolution
- **Geometry Operations**: 3D bounding boxes, periodic boundary conditions, and vector operations
- **Spatial Queries**: Fast cell grid lookups and geometric intersections
- **Morton Indexing**: Space-filling curve indexing for improved cache locality
- **VTK Export**: Direct export to VTK format for visualization
- **Comprehensive Test Suite**: Extensive tests covering all major components
- **CMake Build System**: Modern build configuration with integrated testing

## Project Structure

```
oktal/
├── include/oktal/           # Public header files
│   ├── geometry/            # Geometric primitives and operations
│   │   ├── Box.hpp          # 3D axis-aligned bounding box
│   │   ├── Vec.hpp          # 3D vector type
│   │   └── PeriodicBox.hpp  # Box with periodic boundary conditions
│   ├── octree/              # Octree data structures
│   │   ├── CellOctree.hpp   # Main octree implementation
│   │   ├── CellGrid.hpp     # Grid-based cell organization
│   │   ├── MortonIndex.hpp  # Space-filling curve indexing
│   │   └── OctreeGeometry.hpp # Geometric octree utilities
│   └── io/                  # Input/Output operations
│       └── VtkExport.hpp    # VTK format export
├── src/                     # Implementation files
│   ├── geometry/            # Geometry implementation
│   ├── octree/              # Octree implementation
│   └── io/                  # I/O implementation
├── apps/                    # Example applications
│   ├── create-htgfile.cpp   # HTG file creation utility
│   └── poisson.cpp          # Poisson problem solver example
├── tests/                   # Test suite
│   ├── task1/               # Milestone 1 tests
│   ├── task2/               # Milestone 2 tests
│   ├── task3/               # Milestone 3 tests
│   ├── task4/               # Milestone 4 tests
│   └── task5/               # Milestone 5 tests
└── CMakeLists.txt           # Main build configuration
```

## Requirements

- C++23 compliant compiler (GCC 14+, Clang 17+, or MSVC 17.4+)
- CMake 3.31 or higher
- External dependency: advpt-helpers (automatically fetched during build)

## Building the Project

### Basic Build

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

### Build with Tests

Tests are enabled by default. To disable them:

```bash
cmake -DOKTAL_BUILD_TESTSUITE=OFF ..
```

### Run Tests

```bash
cd build
ctest
```

### Code Quality

Run static analysis with clang-tidy:

```bash
./clang_tidy.sh
```

## Development Status

- [x] Milestone 1: Core geometry and vector operations
- [x] Milestone 2: Box types and Morton indexing
- [x] Milestone 3: VTK export and basic octree
- [x] Milestone 4: Octree iteration and cursor operations
- [x] Milestone 5: Cell grid and Poisson solver

## Core Components

### Geometry Module
- `Box<T>`: Axis-aligned bounding box with scalar templates
- `PeriodicBox<T>`: Bounding box with periodic boundary conditions
- `Vec<T, N>`: Fixed-size N-dimensional vector type

### Octree Module
- `CellOctree<T>`: Main adaptive octree structure
- `CellGrid<T>`: Grid-based spatial partitioning
- `OctreeGeometry`: Geometric operations on octrees
- `MortonIndex`: Space-filling curve implementation for linearized indexing

### I/O Module
- `VtkExport`: Export octree and field data to VTK format for visualization

## Usage Example

```cpp
#include "oktal/octree/CellOctree.hpp"
#include "oktal/geometry/Box.hpp"

using namespace oktal;

// Create a bounding box
Box<double> domain({0, 0, 0}, {1, 1, 1});

// Create an octree with desired resolution
CellOctree<double> octree(domain, /* levels */);

// Perform spatial operations
// ... add your code here
```

## License

This is an educational project developed as part of Advanced Programming Techniques course at FAU.