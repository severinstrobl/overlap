# Copyright (C) 2025 Severin Strobl <git@severin-strobl.de>
#
# SPDX-License-Identifier: MIT
#
# Exact calculation of the overlap volume of spheres and mesh elements.
# http://dx.doi.org/10.1016/j.jcp.2016.02.003


import numpy as np
import overlap
import pyperf

UNIT_HEX = overlap.Hexahedron(
    np.array(
        (
            (-1, -1, -1),
            (1, -1, -1),
            (1, 1, -1),
            (-1, 1, -1),
            (-1, -1, 1),
            (1, -1, 1),
            (1, 1, 1),
            (-1, 1, 1),
        )
    )
)


UNIT_TET = overlap.Tetrahedron(
    np.array(
        (
            (-np.sqrt(3) / 6, -0.5, 0),
            (np.sqrt(3) / 3, 1.0, 0),
            (-np.sqrt(3) / 6, 0.5, 0),
            (0, 0, np.sqrt(6) / 3),
        )
    )
)

CONFIGURATIONS = {
    # sphere - hexahedron
    "hex_overlap_volume[sphere-in-hex]": (
        overlap.Sphere(np.array((0, 0, 0)), 1.0),
        UNIT_HEX,
    ),
    "hex_overlap_volume[hex-in-sphere]": (
        overlap.Sphere(np.array((0, 0, 0)), 5.0),
        UNIT_HEX,
    ),
    "hex_overlap_volume[AABB]": (overlap.Sphere(np.array((5.0, 0, 0)), 1.0), UNIT_HEX),
    # sphere - tetrahedron
    "tet_overlap_volume[sphere-in-tet]": (
        overlap.Sphere(np.array((0, 0, 0)), 0.5),
        UNIT_TET,
    ),
    "tet_overlap_volume[tet-in-sphere]": (
        overlap.Sphere(np.array((0, 0, 0)), 5.0),
        UNIT_TET,
    ),
    "tet_overlap_volume[AABB]": (overlap.Sphere(np.array((2.0, 0, 0)), 0.5), UNIT_TET),
}

runner = pyperf.Runner()

for label, (sphere, element) in CONFIGURATIONS.items():
    runner.bench_func(
        label,
        lambda sphere=sphere, element=element: overlap.overlap_volume(sphere, element),
    )
