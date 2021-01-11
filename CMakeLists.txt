cmake_minimum_required(VERSION 3.10)
project(overlap)

option(OVERLAP_TESTS "Build the tests" on)
option(OVERLAP_COVERAGE "Build with coverage support" off)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Header-only library.
include_directories(${PROJECT_SOURCE_DIR})

# Add compiler flags for C++11.
if(${CMAKE_VERSION} VERSION_LESS "3.1")
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
        "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS}")
    else()
        message(FATAL_ERROR
            "Unsupported compiler detected: ${CMAKE_CXX_COMPILER_ID}")
    endif()
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# Set up the Eigen3 library.
find_package(Eigen3 3.0.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Add dummy interface for coverage configuration.
add_library(coverage INTERFACE)

if(OVERLAP_COVERAGE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(FATAL_ERROR "Generation of coverage information "
                "requires debug build, set CMAKE_BUILD_TYPE=Debug.")
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
