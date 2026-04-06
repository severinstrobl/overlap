# Copyright (C) 2026 Severin Strobl <git@severin-strobl.de>
#
# SPDX-License-Identifier: MIT
#
# Exact calculation of the overlap volume of spheres and mesh elements.
# http://dx.doi.org/10.1016/j.jcp.2016.02.003


import numpy as np
import overlap
import pyperf

UNIT_TET_VERTICES = np.array(
    (
        (-np.sqrt(3) / 6, -0.5, 0),
        (np.sqrt(3) / 3, 1.0, 0),
        (-np.sqrt(3) / 6, 0.5, 0),
        (0, 0, np.sqrt(6) / 3),
    ),
    dtype=float,
)

UNIT_WEDGE_VERTICES = np.array(
    (
        (0, 0, -1),
        (1, 0, -1),
        (0, 1, -1),
        (0, 0, 1),
        (1, 0, 1),
        (0, 1, 1),
    ),
    dtype=float,
)

UNIT_HEX_VERTICES = np.array(
    (
        (-1, -1, -1),
        (1, -1, -1),
        (1, 1, -1),
        (-1, 1, -1),
        (-1, -1, 1),
        (1, -1, 1),
        (1, 1, 1),
        (-1, 1, 1),
    ),
    dtype=float,
)


runner = pyperf.Runner()

runner.bench_func(
    "construction[tet]",
    lambda: overlap.Tetrahedron(UNIT_TET_VERTICES),
)

runner.bench_func(
    "construction[wedge]",
    lambda: overlap.Wedge(UNIT_WEDGE_VERTICES),
)

runner.bench_func(
    "construction[hex]",
    lambda: overlap.Hexahedron(UNIT_HEX_VERTICES),
)
