# Exact calculation of the overlap volume and area of spheres and mesh elements

[![Build Status](https://travis-ci.org/severinstrobl/overlap.svg?branch=master)](https://travis-ci.org/severinstrobl/overlap)
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

## Dependencies

The compile-time dependencies of this code are:
- [Eigen3](http://eigen.tuxfamily.org), tested vith versions 3.3.4 and above
- [compiler supporting C++11](https://en.cppreference.com/w/cpp/compiler_support#cpp11)

The software is currently continuously compiled and tested with the following
compilers:

| Compiler   | Versions |
|------------|---------|
| GNU G++    | 8.3.0, 7.4.0, 6.4.0, 5.5.0 |
| Clang/LLVM | 8.0.1, 7.1.0, 6.0.1, 5.0.2, 4.0.1, 3.9.1 |

Additionally, the Intel C++ compiler starting with version 15.0 should work,
albeit this configuration is not part of the CI process.

## Usage

The library is implemented as a pure header-only library written in plain
C++11. To use it in your code, simply include the header file `overlap.hpp` and
make sure the Eigen3 headers can be found by your compiler or build system.  A
minimal example calculating the overlap of a hexahedron with a side length of 2
centered at the origin and a sphere with radius 1 centered at a corner of the
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
