# Build GoogleTest shipped as a third-party dependency.
# Avoid linker issues on Windows when mixing shared and static libraries.
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
add_subdirectory("${PROJECT_SOURCE_DIR}/third_party/googletest"
                 "third_party/googletest"
                 EXCLUDE_FROM_ALL)

mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

set_target_properties(gtest PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)
set_target_properties(gmock PROPERTIES FOLDER extern)
set_target_properties(gmock_main PROPERTIES FOLDER extern)

# Load GoogleTest module and set up helper function for creating tests.
include(GoogleTest)

macro(overlap_add_test TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    target_link_libraries(${TESTNAME} PRIVATE gtest gtest_main coverage)

    gtest_discover_tests(${TESTNAME}
        WORKING_DIRECTORY ${PROJECT_DIR}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
    )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER test)
endmacro()

# List of unit tests.
set(UNIT_TESTS
    decompose_elements
    detail_clamp
    detail_regularizedWedge
    detail_regularizedWedgeArea
    elements
    normal_newell
    polygon
    sphere
    sphere_element_area_edgecases
    sphere_element_overlap
    sphere_element_overlap_edgecases
    sphere_hex_area
    sphere_tet_area
    sphere_tet_overlap_edgecases
)

# Register the individual unit tests.
foreach(UNIT_TEST ${UNIT_TESTS})
    overlap_add_test(${UNIT_TEST} "${UNIT_TEST}.cpp")
    target_link_libraries(${UNIT_TEST} PUBLIC Eigen3::Eigen)
endforeach()
