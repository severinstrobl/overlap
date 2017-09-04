/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2015-2017 Severin Strobl <severin.strobl@fau.de>
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

#define BOOST_TEST_MODULE detail_regularizedWedgeArea
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "overlap.hpp"

BOOST_AUTO_TEST_CASE(detail_regularizedWedgeArea) {
	scalar_t delta(2e2 * std::numeric_limits<scalar_t>::epsilon());
	using namespace detail;

	// Test distance z.
	BOOST_CHECK_CLOSE(regularizedWedgeArea(1.0, 1.0, 0.25 * M_PI), 0.0, delta);
	BOOST_CHECK_CLOSE(regularizedWedgeArea(1.0, -1.0, 0.25 * M_PI), 0.0,
		delta);

	BOOST_CHECK_CLOSE(regularizedWedgeArea(1.0, detail::tinyEpsilon,
		0.5 * M_PI), M_PI, 5 * delta);

	BOOST_CHECK_CLOSE(regularizedWedgeArea(1.0, -detail::tinyEpsilon,
		0.5 * M_PI), M_PI, 5 * delta);

	// Test angle alpha.
	BOOST_CHECK_CLOSE(regularizedWedgeArea(1.0, 0.0, 0.0), 0.0, delta);
	BOOST_CHECK_CLOSE(regularizedWedgeArea(1.0, 0.0, 0.5 * M_PI), M_PI, delta);
	BOOST_CHECK_CLOSE(regularizedWedgeArea(1.0, 0.0, 0.75 * M_PI),
		2 * M_PI - regularizedWedgeArea(1.0, 0.0, 0.25 * M_PI), delta);
}
