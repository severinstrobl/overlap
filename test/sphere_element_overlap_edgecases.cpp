/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2015-2020 Severin Strobl <severin.strobl@fau.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "gtest/gtest.h"

#include "overlap.hpp"

#include "common.hpp"

TEST(SphereElementOverlapTest, EdgeCases) {
	const std::vector<Sphere> spheres = {
		{{0, 0, 0}, 1.0},
		{{-1, 3, 9}, 10},
		{{0.5, 0.5, 0.5}, 0.50001 * std::sqrt(2.0)},

		{{-1.534427712524021, -0.6526040637766801, 3.823443102163421},
			5.459817873898927},

		{{-2.291983426015874, -3.495618444307236, 2.067917670011271},
			4.797942866073771},

		{{-0.2174878528692581, -3.076535346840716, 0.53771818665538},
			2.856370661961459},

		{{-0.4599350500370143, 4.782461234807632, -2.6269327661783},
			5.2681424105015},

		{{-0.7611917089641156, -0.8319982272779169, -0.004847761469840783},
			2.103084880441632},

		{{2.992123379449451, -0.4987719594414469, 1.44196971013958},
			4.706537474211725},

		{{15.461110118225834, -8.575216180676412, 14.487981174363448},
			21.97121086612232}
	};

	const scalar_t epsilon =
		std::sqrt(std::numeric_limits<scalar_t>::epsilon());

	for(const Sphere& sphere : spheres)
		overlap(sphere, unitHexahedron(), epsilon * sphere.volume);
}
