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


seed = 79_866_982_766_580 % (2**32 - 1)
rng = np.random.default_rng(seed)


def randomized_configuration() -> tuple[overlap.Sphere, overlap.Hexahedron]:
    return overlap.Sphere(rng.uniform(-2, 2, 3), rng.uniform(0.1, 2.5)), UNIT_HEX


runner = pyperf.Runner()
runner.bench_func(
    "hex_overlap_volume[random]",
    lambda: overlap.overlap_volume(*randomized_configuration()),
)
