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

CONFIGURATIONS = {
    "hex_overlap_volume[sphere-in-hex]": overlap.Sphere(np.array((0, 0, 0)), 1.0),
    "hex_overlap_volume[hex-in-sphere]": overlap.Sphere(np.array((0, 0, 0)), 5.0),
    "hex_overlap_volume[AABB]": overlap.Sphere(np.array((5, 0, 0)), 1.0),
}

runner = pyperf.Runner()

for label, sphere in CONFIGURATIONS.items():
    runner.bench_func(label, lambda sphere=sphere: overlap.overlap_volume(sphere, UNIT_HEX))
