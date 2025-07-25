# Copyright (C) 2025 Severin Strobl <git@severin-strobl.de>
#
# SPDX-License-Identifier: MIT
#
# Exact calculation of the overlap volume of spheres and mesh elements.
# http://dx.doi.org/10.1016/j.jcp.2016.02.003

"""Exact calculation of the overlap volume and area of spheres and mesh elements.

Supported primitives
--------------------
- tetrahedron (4 nodes/vertices, `overlap.Tetrahedron`)
- pentahedron/wedge (6 nodes/vertices, `overlap.Wedge`)
- hexahedron (8 nodes/vertices, `overlap.Hexahedron`)

Main functions
--------------
- `overlap.overlap_volume(sphere, element)`: Calculate the overlap volume of
    a sphere and a mesh element.
- `overlap.overlap_area(sphere, element)`: Calculate the overlap area of
    a sphere and a mesh element.
"""

from ._overlap import (
    Hexahedron,
    Sphere,
    Tetrahedron,
    Wedge,
    overlap_area,
    overlap_volume,
)
from ._version import __version__, version

__all__ = [
    "Hexahedron",
    "Sphere",
    "Tetrahedron",
    "Wedge",
    "__doc__",
    "__version__",
    "overlap_area",
    "overlap_volume",
    "version",
]
