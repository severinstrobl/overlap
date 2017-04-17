[![Build Status](https://travis-ci.org/severinstrobl/overlap.svg?branch=testing)](https://travis-ci.org/severinstrobl/overlap)

# Exact calculation of the overlap volume of spheres and mesh elements

Calculating the intersection or overlapping volume of a sphere and one of the
typically used mesh elements such as a tetrahedron or hexahedron is
surprisingly challenging. This header-only library implements a numerically
robust method to determine this volume.

The mathematical expressions and algorithms used in this code are described in
S. Strobl et al.: Exact calculation of the overlap volume of spheres and mesh
elements, Journal of Computational Physics, 2016
(http://dx.doi.org/10.1016/j.jcp.2016.02.003). So if you use the code in
projects resulting in any publications, please cite this paper.

Employing the concepts and routines used for the calculation of the overlap
volume, the intersection or overlap *area* of a sphere and the facets of a mesh
element can also be calculated with this library.

## Dependencies

The compile-time dependencies of this code are:
- Eigen3 (http://eigen.tuxfamily.org)
- C++11 compliant compiler

The software was successfully compiled and tested using the following
compilers:
- GNU G++ compiler (versions 4.8.4, 4.9.3 and 5.4.0)
- Intel C++ compiler (version 15.0.2)
- Clang C++ compiler (versions 3.6.1 and 3.9.1)

## Usage

The library is implemented as a pure header-only library written in plain
C++11. To use it in your code, simply include the header file *overlap.hpp* and
make sure the Eigen3 headers can be found by your compiler or build system.  A
minimal example calculating the overlap of a hexahedron with a side length of 2
centered at the origin and a sphere with radius 1 centered at a corner of the
hexahedron could look something like this:
```
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
This code snippet calculates the correct result (pi/6) for this simple
configuration.

To obtain the overlap area of a sphere and the facets of a tetrahedron, the
function *overlapArea* can be employed as such:
```
vector_t v0{-std::sqrt(3) / 6.0, -1.0 / 2.0, 0};
vector_t v1{std::sqrt(3) / 3.0, 0, 0};
vector_t v2{-std::sqrt(3) / 6.0, +1.0 / 2.0, 0};
vector_t v3{0, 0, std::sqrt(6) / 3.0};

Tetrahedron tet{v0, v1, v2, v3};
Sphere s{{0, 0, 1.5}, 1.25};

auto result = overlapArea(s, tet);

std::cout << "overlap areas per face:" << std::endl;
for(size_t f = 0; f < result.size() - 1; ++f)
    std::cout << "  face #" << f << ": " << result[f] << std::endl

std::cout << "total overlap area: " << result.back() << std::endl;
```
