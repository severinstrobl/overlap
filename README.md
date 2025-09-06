# Exact calculation of the overlap volume and area of spheres and mesh elements

[![Build Status](https://img.shields.io/github/actions/workflow/status/severinstrobl/overlap/ci.yaml?branch=main)](https://github.com/severinstrobl/overlap/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/severinstrobl/overlap/branch/main/graph/badge.svg?token=GQ2L62OXXK)](https://codecov.io/gh/severinstrobl/overlap)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=severinstrobl_overlap&metric=alert_status)](https://sonarcloud.io/summary/overall?id=severinstrobl_overlap)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://github.com/pre-commit/pre-commit)\
[![Release](https://img.shields.io/github/release/severinstrobl/overlap.svg)](https://github.com/severinstrobl/overlap/releases)
[![PyPI](https://img.shields.io/pypi/v/overlap)](https://pypi.org/project/overlap/)\
![C++](https://img.shields.io/badge/C%2B%2B-17%7C20%7C23-blue)
[![Python Version](https://img.shields.io/pypi/pyversions/overlap)](https://pypi.org/project/overlap/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](./LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1016/j.jcp.2016.02.003-blue.svg)](https://dx.doi.org/10.1016/j.jcp.2016.02.003)

Calculating the intersection or overlapping volume of a sphere and one of the
typically used mesh elements such as a tetrahedron or hexahedron is
surprisingly challenging. This header-only library implements a numerically
robust method to determine this volume.

The mathematical expressions and algorithms used in this code are described in
[S. Strobl et al.: Exact calculation of the overlap volume of spheres and mesh
elements, Journal of Computational Physics, 2016](https://dx.doi.org/10.1016/j.jcp.2016.02.003).
So if you use the code in projects resulting in any publications, please cite
this paper.

Employing the concepts and routines used for the calculation of the overlap
volume, the intersection or overlap *area* of a sphere and the facets of a mesh
element can also be calculated with this library.

## Usage

### Supported primitives

The overlap calculation directly supports these element types:

- tetrahedra (4 nodes/vertices, data type `Tetrahedron`)
- pentahedra/wedges/triangular prisms (6 nodes/vertices, data type `Wedge`)
- hexahedra (8 nodes/vertices, data type `Hexahedron`)

The elements must be convex and have to be specified as a list of three-dimensional nodes/vertices,
while the sphere (data type `Sphere`) requires a center point and the radius.

### Node ordering

The element types of the overlap library follow the node numbering conventions
of the [CFD General Notation System (CGNS)](https://cgns.github.io/) project.
Please refer to the CGNS documentation for the order of the nodes of
[hexahedral](https://cgns.github.io/standard/SIDS/convention.html#hexahedral-elements),
[tetrahedral](https://cgns.github.io/standard/SIDS/convention.html#tetrahedral-elements), and
[pentahedral/wedge-shaped](https://cgns.github.io/standard/SIDS/convention.html#pentahedral-elements)
elements of linear order, respectively. Also the ordering of the faces uses
the conventions of CGNS. This should make interfacing this library with
existing codes rather easy, often even without the need to reorder nodes.

### Dependencies

The compile-time dependencies of this code are:

- [Eigen3](http://eigen.tuxfamily.org), tested with versions 3.3.4 and 3.4.0
- [compiler supporting C++17](https://en.cppreference.com/w/cpp/compiler_support/17)

The software is currently continuously compiled and tested with the following
compilers on both x86-64 and ARM64:

| Compiler   | Versions                                                       |
| ---------- | -------------------------------------------------------------- |
| GNU G++    | 14.0.1, 13.2.0, 12.3.0, 11.4.0, 10.5.0, 9.5.0                  |
| Clang/LLVM | 18.1.3, 17.0.6, 16.0.6, 15.0.7, 14.0.0, 13.0.1, 12.0.1, 11.1.0 |

Additionally, the Intel C++ compiler starting with version 19.0 should work,
albeit this configuration is not part of the CI process. All the newer
LLVM-based oneAPI compilers are expected to work.

### C++

[![Try it online](https://img.shields.io/badge/try-online-blue.svg)](https://godbolt.org/z/jc3GfeEnd)

The library is implemented as a header-only library written in C++17. To use it
in your code, simply include the header file `include/overlap/overlap.hpp` and
make sure the **Eigen3** headers can be found by your compiler or build system.
The library exposes two relevant type aliases, namely `Scalar` for `double` and
`Vector` for `Eigen::Matrix<Scalar, 3, 1, Eigen::DontAlign>`, which are used in
the public interface for scalar and vectorial quantities, respectively. In
principle, these types can be adjusted to specific needs, yet reducing the
numerical precision of the scalar floating point type will have a significant
impact on the precision and stability of the calculations.

A minimal example calculating the overlap volume of a hexahedron with a side length
of 2 centered at the origin and a sphere with radius 1 centered at a corner of the
hexahedron could look something like this:

```cpp
using namespace overlap;

const auto vertices = std::array<Vector, 8>{{
    {-1, -1, -1},
    { 1, -1, -1},
    { 1,  1, -1},
    {-1,  1, -1},
    {-1, -1,  1},
    { 1, -1,  1},
    { 1,  1,  1},
    {-1,  1,  1}
}};

const auto hex = Hexahedron{vertices};
const auto sphere = Sphere{Vector::Constant(1), 1};

const auto volume = overlap_volume(sphere, hex);
```

This code snippet calculates the correct result (π/6) for this simple
configuration.

To obtain the overlap area of a sphere and the facets of a tetrahedron, the
function `overlap_area()` can be employed as such:

```cpp
using namespace overlap;

const auto vertices = std::array<Vector, 4>{{
    {-std::sqrt(3) / 6.0, -1.0 / 2.0, 0},
    { std::sqrt(3) / 3.0,  0.0, 0},
    {-std::sqrt(3) / 6.0,  1.0 / 2.0, 0},
    {0, 0, std::sqrt(6) / 3.0},
}};

const auto tet = Tetrahedron{vertices};
const auto sphere = Sphere{{0, 0, 1.5}, 1.25};

const auto result = overlap_area(sphere, tet);

std::cout << "surface area of sphere intersecting tetrahedron: " <<
    result.front() << "\n";

std::cout << "overlap areas per face:\n";
for(auto face_idx = 0u; face_idx < tet.faces.size(); ++face_idx) {
    std::cout << "  face #" << face_idx << ": " <<
        // the indices of the faces are NOT zero-based here!
        result[face_idx + 1] << "\n";
}

std::cout << "total surface area of tetrahedron intersecting sphere: " <<
    result.back() << "\n";
```

### Python

The Python version of the `overlap` library is available via the [Python
Package Index (PyPI)](https://pypi.org/project/overlap/), so for most
environments installation should be possible simply via `pip install overlap`.

In case no pre-built package or *wheel* is available for your system, compilation of the
wrapper code is required which in turn requires the requirements listed above
for the C++ version to be fulfilled.

The interface of Python version closely resembles the interface of the C++ version:

```python
import numpy as np
import overlap

vertices = np.array(
    (
        (-1, -np.sqrt(1.0 / 3.0), -np.sqrt(1.0 / 6.0)),
        (1, -np.sqrt(1.0 / 3.0), -np.sqrt(1.0 / 6.0)),
        (0, np.sqrt(4.0 / 3.0), -np.sqrt(1.0 / 6.0)),
        (0, 0, np.sqrt(3.0 / 2.0)),
    )
)

tet = overlap.Tetrahedron(vertices)
sphere = overlap.Sphere((0, 0, 0.5), 1)
result = overlap.overlap_volume(sphere, tet)
```

Calculation of the overlap area instead of the overlap volume is possible via
the function `overlap_area()` of the package.

## License

The `overlap` library is distributed under the MIT license, please refer to the
[LICENSE](LICENSE) file for the full license text.

This distribution uses external third-party software covered by separate
license terms. For details, please consult the corresponding license terms
of the respective package.

| Package                                            | License      |
| -------------------------------------------------- | ------------ |
| [Eigen](http://eigen.tuxfamily.org)                | MPL2         |
| [pybind11](https://github.com/pybind/pybind11)     | 3-clause BSD |
| [doctest](https://github.com/doctest/doctest)      | MIT          |
| [nanobench](https://github.com/martinus/nanobench) | MIT          |
