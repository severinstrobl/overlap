cmake_minimum_required(VERSION 3.18)
project(overlap)

enable_language(CXX)

option(OVERLAP_WITH_PYTHON "Build the Python bindings" on)
option(OVERLAP_WITH_TESTS "Build the tests" on)
option(OVERLAP_WITH_BENCHMARKS "Build the benchmarks" off)
option(OVERLAP_WITH_COVERAGE "Build with coverage support for tests" off)

# cmake-lint: disable=C0103
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
  endif()

  set(CMAKE_CXX_STANDARD_REQUIRED on)
  set(CMAKE_CXX_EXTENSIONS off)

  # prefer first detected Python version to avoid issues with virtual envs
  set(Python_FIND_STRATEGY LOCATION)
endif()

# enable generation of Clang compilation database for IDE integration
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

# set up the dependency on Eigen3
include(FetchContent)

fetchcontent_declare(
  eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3147391d946bb4b6c68edd901f2add6ac1f31f8c # v3.4.0
  GIT_SHALLOW TRUE
  SOURCE_SUBDIR
  cmake # set to invalid path to avoid building and registering tests
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS
  3.4.0
  NO_MODULE
  NAMES
  Eigen3
)

fetchcontent_makeavailable(eigen3)

if(NOT TARGET Eigen3::Eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(
    Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${eigen3_SOURCE_DIR}
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
  add_subdirectory(benchmarks/cpp)
endif()

# optional: build the Python bindings using pybind11
if(OVERLAP_WITH_PYTHON)
  add_subdirectory(python)
endif()
