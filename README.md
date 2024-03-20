# Exact calculation of the overlap volume and area of spheres and mesh elements

![Build Status](https://img.shields.io/github/actions/workflow/status/severinstrobl/overlap/build.yaml?branch=master)
[![codecov](https://codecov.io/gh/severinstrobl/overlap/branch/master/graph/badge.svg?token=GQ2L62OXXK)](https://codecov.io/gh/severinstrobl/overlap)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](./COPYING)
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

# Usage

## Supported primitives

The overlap calculation directly supports these element types:

- tetrahedra (4 nodes/vertices, data type `Tetrahedron`)
- pentahedra/wedges/triangular prisms (6 nodes/vertices, data type `Wedge`)
- hexahedra (8 nodes/vertices, data type `Hexahedron`)

The elements must be convex and have to be specified as a list of three-dimensional nodes/vertices,
while the sphere (data type `Sphere`) requires a center point and the radius.

## Node ordering

The element types of the overlap library follow the node numbering conventions
of the [CFD General Notation System (CGNS)](https://cgns.github.io/) project.
Please refer to the CGNS documentation for the order of the nodes of
[hexahedral](https://cgns.github.io/CGNS_docs_current/sids/conv.html#unst_hexa),
[tetrahedral](https://cgns.github.io/CGNS_docs_current/sids/conv.html#unst_tetra), and
[pentahedral/wedge-shaped](https://cgns.github.io/CGNS_docs_current/sids/conv.html#unst_penta)
elements of linear order, respectively. Also the ordering of the faces uses
the conventions of CGNS. This should make interfacing this library with
existing codes rather easy, often even without the need to reorder nodes.

## Dependencies

The compile-time dependencies of this code are:
- [Eigen3](http://eigen.tuxfamily.org), tested with versions 3.3.4 and above
- [compiler supporting C++11](https://en.cppreference.com/w/cpp/compiler_support#cpp11)

The software is currently continuously compiled and tested with the following
compilers:

| Compiler   | Versions |
|------------|----------|
| GNU G++    | 10.3.0, 9.3.0, 8.4.0, 7.5.0 |
| Clang/LLVM | 12.0.0, 11.0.0, 10.0.0, 9.0.1, 8.0.1 |

Additionally, the Intel C++ compiler starting with version 15.0 should work,
albeit this configuration is not part of the CI process.


## C++

The library is implemented as a pure header-only library written in plain
C++11. To use it in your code, simply include the header file `overlap.hpp` and
make sure the **Eigen3** headers can be found by your compiler or build system.
The library creates two relevant type aliases, namely `scalar_t` for `double`
and `vector_t` for `Eigen::Matrix<scalar_t, 3, 1, Eigen::DontAlign>`, which are
used in the public interface for scalar and vectorial quantities, respectively.
In principle, these types can be adjusted to specific needs, yet reducing the
numerical precision of the scalar floating point type will have a significant
impact on the precision and stability of the calculations.

A minimal example calculating the overlap of a hexahedron with a side length of
2 centered at the origin and a sphere with radius 1 centered at a corner of the
hexahedron could look something like this:
```cpp
vector_t v0{-1, -1, -1};
vector_t v1{ 1, -1, -1};
vector_t v2{ 1,  1, -1};
vector_t v3{-1,  1, -1};
vector_t v4{-1, -1,  1};
vector_t v5{ 1, -1,  1};
vector_t v6{ 1,  1,  1};
vector_t v7{-1,  1,  1};

Hexahedron hex{v0, v1, v2, v3, v4, v5, v6, v7};
Sphere s{vector_t::Constant(1), 1};

scalar_t result = overlap(s, hex);
```
This code snippet calculates the correct result (&pi;/6) for this simple
configuration.

To obtain the overlap area of a sphere and the facets of a tetrahedron, the
function `overlapArea()` can be employed as such:
```cpp
vector_t v0{-std::sqrt(3) / 6.0, -1.0 / 2.0, 0};
vector_t v1{std::sqrt(3) / 3.0, 0, 0};
vector_t v2{-std::sqrt(3) / 6.0, +1.0 / 2.0, 0};
vector_t v3{0, 0, std::sqrt(6) / 3.0};

Tetrahedron tet{v0, v1, v2, v3};
Sphere s{{0, 0, 1.5}, 1.25};

auto result = overlapArea(s, tet);

std::cout << "surface area of sphere intersecting tetrahedron: " <<
    result[0] << std::endl;

std::cout << "overlap areas per face:" << std::endl;
// The indices of the faces are NOT zero-based here!
for(size_t f = 1; f < result.size() - 1; ++f)
    std::cout << "  face #" << (f - 1) << ": " << result[f] << std::endl;

std::cout << "total surface area of tetrahedron intersecting sphere: " <<
    result.back() << std::endl;
```

## Python

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

vertices = np.array((
    (-1, -np.sqrt(1./3.), -np.sqrt(1./6.)),
    (1, -np.sqrt(1./3.), -np.sqrt(1./6.)),
    (0, np.sqrt(4./3.), -np.sqrt(1./6.)),
    (0, 0, np.sqrt(3./2.))
))

tet = overlap.Tetrahedron(vertices)
sphere = overlap.Sphere((0, 0, 0.5), 1)

result = overlap.overlap(sphere, tet)
```

Calculation of the overlap area instead of the overlap volume is possible via
the function `overlap_area()` of the package.

# License

The `overlap` library is distributed under the GNU General Public
License v3, please refer to the [LICENSE](LICENSE) file for the full license
text.

This distribution bundles external third-party software covered by separate
license terms. For details please consult the corresponding license terms
included with each package in the respective subdirectory.

| Package     | License |
|-------------|----------|
| [Eigen](http://eigen.tuxfamily.org) | MPL2 |
| [Google Test](https://github.com/google/googletest) | 3-clause BSD |
| [pybind11](https://github.com/pybind/pybind11) | 3-clause BSD |
