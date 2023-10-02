cmake_minimum_required(VERSION 3.18)
project(overlap)

enable_language(CXX)

option(OVERLAP_WITH_PYTHON "Build the Python bindings" on)
option(OVERLAP_WITH_TESTS "Build the tests" on)
option(OVERLAP_WITH_BENCHMARKS "Build the benchmarks" off)
option(OVERLAP_WITH_COVERAGE "Build with coverage support for tests" off)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
  endif()

  set(CMAKE_CXX_STANDARD_REQUIRED on)
  set(CMAKE_CXX_EXTENSIONS off)

  # prefer first detected Python verion to avoid issues with virtual envs
  set(Python_FIND_STRATEGY LOCATION)
endif()

# enable generation of Clang compilation database for IDE integration
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

# set up the dependency on Eigen3
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/third_party/eigen/cmake")
find_package(Eigen3 3.3.0 NO_MODULE QUIET)
if(NOT TARGET Eigen3::Eigen)
  message(STATUS "Failed to find Eigen3 >= 3.3.0, using internal version")

  add_library(Eigen3::Eigen IMPORTED INTERFACE)
  set_property(
    TARGET Eigen3::Eigen PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                  "${PROJECT_SOURCE_DIR}/third_party/eigen"
  )
endif()

add_library(overlap INTERFACE)
add_library(overlap::overlap ALIAS overlap)
add_library(overlap::headers ALIAS overlap)

target_compile_features(overlap INTERFACE cxx_std_17)

include(GNUInstallDirs)

target_include_directories(
  overlap INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# optional: C++ unit tests
if(OVERLAP_WITH_TESTS)
  enable_testing()
  add_subdirectory(tests/src)
endif()

# optional: C++ benchmarks
if(OVERLAP_WITH_BENCHMARKS)
  enable_testing()
  add_subdirectory(benchmarks)
endif()


# optional: build the Python bindings using pybind11
if(OVERLAP_WITH_PYTHON)
  add_subdirectory(python/overlap)
endif()
