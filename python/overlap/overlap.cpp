/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2021-2022 Severin Strobl <git@severin-strobl.de>
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
#include <locale>
#include <map>
#include <typeindex>
#include <typeinfo>

#include "pybind11/eigen.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "overlap/overlap.hpp"

namespace py = pybind11;

template<typename... Args>
using overload_cast_ = py::detail::overload_cast_impl<Args...>;

using namespace py::literals;

using namespace overlap;

template<typename Element,
         typename = std::enable_if_t<detail::is_element_v<Element>>>
void createBindings(py::module& m) {
  static const auto element_names = std::map<std::type_index, std::string>{
      {std::type_index(typeid(Tetrahedron)), "Tetrahedron"},
      {std::type_index(typeid(Wedge)), "Wedge"},
      {std::type_index(typeid(Hexahedron)), "Hexahedron"}};

  const auto& name = element_names.at(std::type_index(typeid(Element)));
  const auto nameLower =
      std::tolower(name.front(), std::locale{}) + name.substr(1);

  static constexpr auto nrVertices = detail::num_vertices<Element>();

  py::class_<Element>(m, name.c_str())
      .def(py::init<std::array<Vector, nrVertices>>(), py::arg("vertices"))
      .def(py::init([](py::array_t<double> vertices) {
             auto proxy = vertices.unchecked<2>();
             if (proxy.shape(0) != nrVertices || proxy.shape(1) != 3) {
               throw std::invalid_argument{
                   "invalid shape for vertex list, must be (" +
                   std::to_string(nrVertices) + ", 3)"};
             }

             std::array<Vector, nrVertices> tmp{};
             for (py::ssize_t v = 0; v < proxy.shape(0); ++v) {
               tmp[v] = Vector{proxy(v, 0), proxy(v, 1), proxy(v, 2)};
             }

             return Element{std::move(tmp)};
           }),
           py::arg("vertices"))
      .def_readonly("vertices", &Element::vertices,
                    "Return the vertices of the element.")
      .def_readonly("center", &Element::center,
                    "Return the center point of the element.")
      .def_readonly("volume", &Element::volume,
                    "Return the volume of the element.")
      .def_property_readonly(
          "surface_area",
          [](const Element& elem) { return elem.surface_area(); },
          "Return the surface area of the element.");

  m.def(
      "overlap_volume",
      overload_cast_<const Sphere&, const Element&>()(&overlap_volume<Element>),
      "sphere"_a, py::arg{nameLower.c_str()},
      ("Calculate the overlap volume of a sphere and a " + nameLower + ".")
          .c_str());

  m.def("overlap_area",
        overload_cast_<const Sphere&, const Element&>()(&overlap_area<Element>),
        "sphere"_a, py::arg{nameLower.c_str()},
        ("Calculate the overlap area of a sphere and a " + nameLower + ".")
            .c_str());
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
      .def(py::init<Vector, Scalar>(), py::arg("center"), py::arg("radius"))
      .def_readonly("center", &Sphere::center,
                    "Return the center point of the sphere.")
      .def_readonly("radius", &Sphere::radius,
                    "Return the radius of the sphere.")
      .def_readonly("volume", &Sphere::volume,
                    "Return the volume of the sphere.")
      .def_property_readonly(
          "surface_area", [](const Sphere& s) { return s.surface_area(); },
          "Return the surface area of the sphere.");

  createBindings<Tetrahedron>(m);
  createBindings<Wedge>(m);
  createBindings<Hexahedron>(m);
}
