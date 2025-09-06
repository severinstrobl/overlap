if(NOT SKBUILD)
  # CMake 3.18 should support this, but there is a bug w.r.t. multiple queries
  # with different components
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.22.0")
    set(_python_components Interpreter Development.Module)
  else()
    set(_python_components Interpreter Development)
  endif()

  find_package(
    Python 3.5
    COMPONENTS ${_python_components}
    REQUIRED
  )
endif()

if(NOT TARGET pybind11::headers)
  include(FetchContent)

  fetchcontent_declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG f5fbe867d2d26e4a0a9177a51f6e568868ad3dc8 # v3.0.1
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL
  )

  fetchcontent_makeavailable(pybind11)
endif()

set(_package_name "overlap")

pybind11_add_module(${_package_name}_python MODULE ${_package_name}.cpp)

set_target_properties(
  ${_package_name}_python PROPERTIES OUTPUT_NAME "_${_package_name}"
)

target_link_libraries(
  ${_package_name}_python PRIVATE ${_package_name}::headers Eigen3::Eigen
)

if(NOT SKBUILD)
  include(GNUInstallDirs)

  string(CONCAT python_install_dir "${CMAKE_INSTALL_LIBDIR}"
                "/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}"
                "/site-packages/${_package_name}"
  )

  install(TARGETS ${_package_name}_python
          LIBRARY DESTINATION ${python_install_dir} COMPONENT python
  )
else()
  install(TARGETS ${_package_name}_python LIBRARY DESTINATION ${_package_name})
endif()
