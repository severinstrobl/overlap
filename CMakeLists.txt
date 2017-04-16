cmake_minimum_required(VERSION 2.8)
project(overlap)

# Header-only library.
include_directories(${PROJECT_SOURCE_DIR})

# Add compiler flags for C++11.
if(${CMAKE_VERSION} VERSION_GREATER "3.1.0")
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
else()
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
        "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        add_compile_options(-std=c++11)
    else()
        message(FATAL_ERROR
            "Unsupported compiler detected: ${CMAKE_CXX_COMPILER_ID}")
    endif()
endif()

# Set up the Eigen3 library.
find_package(Eigen3 3.0.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Activate the unit tests.
enable_testing()
add_subdirectory(test)
