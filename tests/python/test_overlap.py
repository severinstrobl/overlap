# Copyright (C) 2021-2022 Severin Strobl <git@severin-strobl.de>
#
# SPDX-License-Identifier: MIT
#
# Exact calculation of the overlap volume of spheres and mesh elements.
# http://dx.doi.org/10.1016/j.jcp.2016.02.003

import numpy as np
import overlap


def test_sphere():
    sphere = overlap.Sphere((1, 1, 1), 1)
    np.testing.assert_almost_equal(sphere.center, [1, 1, 1])
    np.testing.assert_almost_equal(sphere.radius, 1.0)
    np.testing.assert_almost_equal(sphere.volume, 4.0 / 3.0 * np.pi)


class TestElements:
    def test_tetrahedron(self):
        vertices = np.array(
            (
                (-1, -np.sqrt(1.0 / 3.0), -np.sqrt(1.0 / 6.0)),
                (1, -np.sqrt(1.0 / 3.0), -np.sqrt(1.0 / 6.0)),
                (0, np.sqrt(4.0 / 3.0), -np.sqrt(1.0 / 6.0)),
                (0, 0, np.sqrt(3.0 / 2.0)),
            )
        )

        tet = overlap.Tetrahedron(vertices)

        np.testing.assert_almost_equal(tet.volume, 8.0 / (6 * np.sqrt(2)))
        np.testing.assert_almost_equal(tet.center, [0, 0, 0])
        np.testing.assert_almost_equal(tet.surface_area, 4 * np.sqrt(3))

    def test_wedge(self):
        vertices = [
            (-1, -1, -1),
            (1, -1, -1),
            (1, 1, -1),
            (-1, -1, 1),
            (1, -1, 1),
            (1, 1, 1),
        ]

        wedge = overlap.Wedge(vertices)

        np.testing.assert_almost_equal(wedge.volume, 4.0)
        np.testing.assert_almost_equal(wedge.center, [1.0 / 3.0, -1.0 / 3.0, 0])
        np.testing.assert_almost_equal(wedge.surface_area, 12.0 + 2 * np.sqrt(8))

    def test_hexahedron(self):
        vertices = [
            (-1, -1, -1),
            (1, -1, -1),
            (1, 1, -1),
            (-1, 1, -1),
            (-1, -1, 1),
            (1, -1, 1),
            (1, 1, 1),
            (-1, 1, 1),
        ]

        hexa = overlap.Hexahedron(vertices)

        np.testing.assert_almost_equal(hexa.volume, 8.0)
        np.testing.assert_almost_equal(hexa.center, [0, 0, 0])
        np.testing.assert_almost_equal(hexa.surface_area, 24.0)


class TestOverlap:
    def test_tetrahedron(self):
        tet = overlap.Tetrahedron([(0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 2)])
        sphere = overlap.Sphere((0, 0, 0), 1)

        np.testing.assert_almost_equal(overlap.overlap_volume(sphere, tet), sphere.volume / 8)

    def test_wedge(self):
        vertices = [
            (-1, -1, -1),
            (1, -1, -1),
            (1, 1, -1),
            (-1, -1, 1),
            (1, -1, 1),
            (1, 1, 1),
        ]

        wedge = overlap.Wedge(vertices)
        sphere = overlap.Sphere((1, 1, 1), 1)

        np.testing.assert_almost_equal(overlap.overlap_volume(sphere, wedge), np.pi / 12)

    def test_hexahedron(self):
        vertices = [
            (-1, -1, -1),
            (1, -1, -1),
            (1, 1, -1),
            (-1, 1, -1),
            (-1, -1, 1),
            (1, -1, 1),
            (1, 1, 1),
            (-1, 1, 1),
        ]

        hexa = overlap.Hexahedron(vertices)
        sphere = overlap.Sphere((1, 1, 1), 1)

        np.testing.assert_almost_equal(overlap.overlap_volume(sphere, hexa), np.pi / 6)


class TestOverlapArea:
    def test_wedge(self):
        vertices = [
            (-1, -1, -1),
            (1, -1, -1),
            (1, 1, -1),
            (-1, -1, 1),
            (1, -1, 1),
            (1, 1, 1),
        ]

        wedge = overlap.Wedge(vertices)
        sphere = overlap.Sphere((0, 0, 0), 1)

        np.testing.assert_almost_equal(
            overlap.overlap_area(sphere, wedge),
            [
                sphere.surface_area / 2,
                0.0,
                0.0,
                0.0,
                sphere.radius**2 * np.pi,
                0.0,
                sphere.radius**2 * np.pi,
            ],
        )

    def test_hexahedron(self):
        sphere = overlap.Sphere((1, 1, 1), 1)

        vertices = [
            (-1, -1, -1),
            (1, -1, -1),
            (1, 1, -1),
            (-1, 1, -1),
            (-1, -1, 1),
            (1, -1, 1),
            (1, 1, 1),
            (-1, 1, 1),
        ]

        hexa = overlap.Hexahedron(vertices)

        face_overlap = sphere.radius**2 * np.pi / 4

        np.testing.assert_almost_equal(
            overlap.overlap_area(sphere, hexa),
            [
                sphere.surface_area / 8,
                0.0,
                0.0,
                face_overlap,
                face_overlap,
                0,
                face_overlap,
                3 * face_overlap,
            ],
        )
