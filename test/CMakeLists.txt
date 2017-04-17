# Set up the Boost unit test library.
find_package(Boost 1.40.0 COMPONENTS unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Register the individual unit tests.
set(TEST_SOURCES
    decompose_elements.cpp
    normal_newell.cpp
    regularized_wedge.cpp
    sphere_element_overlap.cpp
    sphere_element_overlap_edgecases.cpp
    sphere_tet_overlap_edgecases.cpp
)

foreach(TEST_SOURCE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SOURCE} ${BOOST_UNIT_TEST})
    target_link_libraries(${TEST_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach(TEST_SOURCE)
