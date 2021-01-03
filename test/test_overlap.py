# Exact calculation of the overlap volume of spheres and mesh elements.
# http://dx.doi.org/10.1016/j.jcp.2016.02.003

# Copyright (C) 2021 Severin Strobl <git@severin-strobl.de>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import overlap


def test_sphere():
    s = overlap.Sphere((1, 1, 1), 1)
    np.testing.assert_almost_equal(s.center, [1, 1, 1])
    np.testing.assert_almost_equal(s.radius, 1.0)
    np.testing.assert_almost_equal(s.volume, 4./3. * np.pi)


class TestElements:
    def test_tetrahedron(self):
        vertices = np.array((
            (-1, -np.sqrt(1./3.), -np.sqrt(1./6.)),
            (1, -np.sqrt(1./3.), -np.sqrt(1./6.)),
            (0, np.sqrt(4./3.), -np.sqrt(1./6.)),
            (0, 0, np.sqrt(3./2.))
        ))
        tet = overlap.Tetrahedron(vertices)

        np.testing.assert_almost_equal(tet.volume, 8. / (6 * np.sqrt(2)))
        np.testing.assert_almost_equal(tet.center, [0, 0, 0])
        np.testing.assert_almost_equal(tet.surface_area, 4 * np.sqrt(3))

    def test_wedge(self):
        vertices = [
            (-1, -1, -1), (1, -1, -1), (1,  1, -1),
            (-1, -1,  1), (1, -1,  1), (1,  1,  1),
        ]

        wedge = overlap.Wedge(vertices)

        np.testing.assert_almost_equal(wedge.volume, 4.0)
        np.testing.assert_almost_equal(wedge.center, [1./3., -1./3., 0])
        np.testing.assert_almost_equal(wedge.surface_area,
                                       12.0 + 2 * np.sqrt(8))

    def test_hexahedron(self):
        vertices = [
            (-1, -1, -1), (1, -1, -1), (1,  1, -1), (-1,  1, -1),
            (-1, -1,  1), (1, -1,  1), (1,  1,  1), (-1,  1,  1)
        ]

        hexa = overlap.Hexahedron(vertices)

        np.testing.assert_almost_equal(hexa.volume, 8.0)
        np.testing.assert_almost_equal(hexa.center, [0, 0, 0])
        np.testing.assert_almost_equal(hexa.surface_area, 24.0)


class TestOverlap:
    def test_wedge(self):
        vertices = [
            (-1, -1, -1), (1, -1, -1), (1,  1, -1),
            (-1, -1,  1), (1, -1,  1), (1,  1,  1),
        ]

        wedge = overlap.Wedge(vertices)
        s = overlap.Sphere((1, 1, 1), 1)

        np.testing.assert_almost_equal(overlap.overlap(s, wedge), np.pi / 12)

    def test_hexahedron(self):
        vertices = [
            (-1, -1, -1), (1, -1, -1), (1,  1, -1), (-1,  1, -1),
            (-1, -1,  1), (1, -1,  1), (1,  1,  1), (-1,  1,  1)
        ]

        hexa = overlap.Hexahedron(vertices)
        s = overlap.Sphere((1, 1, 1), 1)

        np.testing.assert_almost_equal(overlap.overlap(s, hexa), np.pi / 6)


class TestOverlapArea:
    def test_wedge(self):
        vertices = [
            (-1, -1, -1), (1, -1, -1), (1,  1, -1),
            (-1, -1,  1), (1, -1,  1), (1,  1,  1),
        ]

        wedge = overlap.Wedge(vertices)
        s = overlap.Sphere((0, 0, 0), 1)

        np.testing.assert_almost_equal(overlap.overlap_area(s, wedge), [
            s.surface_area / 2,
            0., 0., 0., s.radius**2 * np.pi, 0.,
            s.radius**2 * np.pi
        ])

    def test_hexahedron(self):
        s = overlap.Sphere((1, 1, 1), 1)

        vertices = [
            (-1, -1, -1), (1, -1, -1), (1,  1, -1), (-1,  1, -1),
            (-1, -1,  1), (1, -1,  1), (1,  1,  1), (-1,  1,  1)
        ]

        hexa = overlap.Hexahedron(vertices)

        face_overlap = s.radius**2 * np.pi / 4

        np.testing.assert_almost_equal(overlap.overlap_area(s, hexa), [
            s.surface_area / 8,
            0., 0., face_overlap, face_overlap, 0, face_overlap,
            3 * face_overlap
        ])
