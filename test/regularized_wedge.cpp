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

#define BOOST_TEST_MODULE regularized_wedge
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "overlap.hpp"

BOOST_AUTO_TEST_CASE(regularized_wedge) {
	scalar_t delta(2e2 * std::numeric_limits<scalar_t>::epsilon());
	using namespace detail;

	// Test distance d.
	BOOST_CHECK_CLOSE(regularizedWedge(1.0, 1.0, 0.25 * M_PI), 0.0, delta);
	BOOST_CHECK_CLOSE(regularizedWedge(1.0, detail::tinyEpsilon, 0.25 * M_PI),
		M_PI / 6.0, 5 * delta);

	BOOST_CHECK_CLOSE(regularizedWedge(1.0, detail::tinyEpsilon, 0.5 * M_PI),
		M_PI / 3.0, 5 * delta);

	// Test angle alpha.
	BOOST_CHECK_CLOSE(regularizedWedge(1.0, 0.5, 0.0), 0.0, delta);
	BOOST_CHECK_CLOSE(regularizedWedge(1.0, 0.5, 0.5 * M_PI), 5.0 * M_PI / 48.0,
		delta);

	// Test general version for simple angles and base points very close to the
	// center.
	BOOST_CHECK_CLOSE(regularizedWedge(1.0,
		std::numeric_limits<scalar_t>::epsilon(), 0.25 * M_PI), M_PI / 6.0,
		delta);

	BOOST_CHECK_CLOSE(regularizedWedge(1.0,
		std::numeric_limits<scalar_t>::epsilon(), 0.5 * M_PI), M_PI / 3.0,
		delta);
}

BOOST_AUTO_TEST_CASE(regularized_wedge_pi_half) {
	const scalar_t delta(std::numeric_limits<scalar_t>::epsilon());
	const scalar_t alpha = 0.5 * M_PI;

	using namespace detail;

	BOOST_CHECK_CLOSE(regularizedWedge(1.0, 0.5, alpha, 0.5 * std::cos(alpha +
		0.5 * M_PI)), regularizedWedge(1.0, 0.5, alpha - delta, 0.5 *
		std::cos(alpha - delta + 0.5 * M_PI)), 5e2 * delta);

	BOOST_CHECK_CLOSE(regularizedWedge(1.0, 0.5, alpha, 0.5 * std::cos(alpha +
		0.5 * M_PI)), regularizedWedge(1.0, 0.5, alpha + delta, 0.5 *
		std::cos(alpha + delta + 0.5 * M_PI)), 5e2 * delta);

	BOOST_CHECK_CLOSE(regularizedWedge(1.0, 0.5, alpha, -0.5 * std::cos(alpha +
		0.5 * M_PI)), regularizedWedge(1.0, 0.5, alpha - delta, -0.5 *
		std::cos(alpha - delta + 0.5 * M_PI)), 5e2 * delta);

	BOOST_CHECK_CLOSE(regularizedWedge(1.0, 0.5, alpha, -0.5 * std::cos(alpha +
		0.5 * M_PI)), regularizedWedge(1.0, 0.5, alpha + delta, -0.5 *
		std::cos(alpha + delta + 0.5 * M_PI)), 5e2 * delta);
}
