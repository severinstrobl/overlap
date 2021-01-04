/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2021 Severin Strobl <git@severin-strobl.de>
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

#include <cctype>

#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

#include "overlap.hpp"

namespace py = pybind11;

template<typename... Args>
using overload_cast_ = py::detail::overload_cast_impl<Args...>;

using namespace py::literals;

template<typename Element>
void createBindings(py::module& m) {
	static_assert(std::is_same<Element, Tetrahedron>::value ||
		std::is_same<Element, Wedge>::value ||
		std::is_same<Element, Hexahedron>::value,
		"Invalid element type detected.");

  static const std::string name = std::is_same<Element, Tetrahedron>::value ?
    "Tetrahedron" : std::is_same<Element, Wedge>::value ?
    "Wedge" : "Hexahedron";

  static const std::string nameLower =
    static_cast<char>(std::tolower(static_cast<char>(name[0]))) +
    name.substr(1);

  static constexpr std::size_t nrVertices = element_trait<Element>::nrVertices;

  py::class_<Element>(m, name.c_str())
    .def(py::init<std::array<vector_t, nrVertices>>())
    .def(py::init([](py::array_t<double> vertices) {
      auto proxy = vertices.unchecked<2>();
      if(proxy.shape(0) != nrVertices || proxy.shape(1) != 3)
        throw std::invalid_argument{"invalid shape for vertex list, must be (" +
                                    std::to_string(nrVertices) + ", 3)"};

      std::array<vector_t, nrVertices> tmp{};
      for(py::ssize_t v = 0; v < proxy.shape(0); ++v)
        tmp[v] = vector_t{proxy(v, 0), proxy(v, 1), proxy(v, 2)};

      return Element(tmp);
    }))
    .def_readonly("vertices", &Element::vertices,
      "Return the vertices of the element.")
    .def_readonly("center", &Element::center,
      "Return the center point of the element.")
    .def_readonly("volume", &Element::volume,
      "Return the volume of the element.")
    .def_property_readonly("surface_area", [](const Element& elem) {
        return elem.surfaceArea();
      },
      "Return the surface area of the element.");

  m.def("overlap", overload_cast_<const Sphere&, const Element&>()
    (&overlap<Element>), "sphere"_a, py::arg("nameLower.c_str()"),
    ("Calculate the overlap volume of a sphere and a " + nameLower +
      ".").c_str());

  m.def("overlap_area", overload_cast_<const Sphere&, const Element&>()
    (&overlapArea<Element>), "sphere"_a, py::arg("nameLower.c_str()"),
    ("Calculate the overlap area of a sphere and a " + nameLower +
      ".").c_str());
}

PYBIND11_MODULE(_overlap, m) {
  m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        This originates fom CPP

        .. currentmodule:: overlap

        .. autosummary::
            :toctree: _generate
    )pbdoc";

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif

  py::class_<Sphere>(m, "Sphere")
    .def(py::init<vector_t, scalar_t>())
    .def_readonly("center", &Sphere::center,
      "Return the center point of the sphere.")
    .def_readonly("radius", &Sphere::radius,
      "Return the radius of the sphere.")
    .def_readonly("volume", &Sphere::volume,
      "Return the volume of the sphere.")
    .def_property_readonly("surface_area", [](const Sphere& s) {
        return s.surfaceArea();
      },
      "Return the surface area of the sphere.");

  createBindings<Tetrahedron>(m);
  createBindings<Wedge>(m);
  createBindings<Hexahedron>(m);
}
