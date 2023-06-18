cmake_minimum_required(VERSION 3.12)

if(NOT SKBUILD OR SKBUILD_PROJECT_VERSION MATCHES "dev")
  project(overlap LANGUAGES CXX)
else()
  project(overlap VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES CXX)
endif()

option(OVERLAP_TESTS "Build the tests" on)
option(OVERLAP_COVERAGE "Build with coverage support" off)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# prefer first detected Python verion to avoid issues with virtual envs
set(Python_FIND_STRATEGY LOCATION)

# Header-only library.
include_directories(${PROJECT_SOURCE_DIR})

# Add compiler flags for C++11.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set up the Eigen3 library.
find_package(Eigen3 3.2.0 NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  message(STATUS "Failed to find Eigen3 >= 3.2.0, using internal version")

  add_library(Eigen3::Eigen IMPORTED INTERFACE)
  set_property(
    TARGET Eigen3::Eigen PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                  "${PROJECT_SOURCE_DIR}/third_party/eigen"
  )
endif()

# Add dummy interface for coverage configuration.
add_library(coverage INTERFACE)

if(OVERLAP_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "[GNU|Clang]")
  if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(FATAL_ERROR "Generation of coverage information "
                        "requires debug build, set CMAKE_BUILD_TYPE=Debug."
    )
  endif()

  message(STATUS "Enabling generation of coverage reports")
  target_compile_options(coverage INTERFACE --coverage)
  target_link_options(coverage INTERFACE --coverage)
endif()

# Activate the unit tests.
if(OVERLAP_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

# Build the Python bindings using pybind11.
option(OVERLAP_PYTHON "Build the Python bindings" ON)
if(OVERLAP_PYTHON)
  add_subdirectory(python)
endif()
