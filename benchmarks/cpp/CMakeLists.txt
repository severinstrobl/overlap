if(NOT TARGET doctest::doctest)
  include(FetchContent)

  fetchcontent_declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG 2d0a9359a60c51affe2a9bebb1be1dca47868151 # v2.5.3
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL
  )

  fetchcontent_makeavailable(doctest)
  list(APPEND CMAKE_MODULE_PATH ${doctest_SOURCE_DIR}/scripts/cmake)
endif()

if(NOT TARGET nanobench)
  include(FetchContent)

  fetchcontent_declare(
    nanobench
    GIT_REPOSITORY https://github.com/martinus/nanobench.git
    GIT_TAG b709728dba4df0356d467abfb0ba20781a8031c0 # v4.6.0
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL
  )

  fetchcontent_makeavailable(nanobench)
endif()

# list of benchmarks
set(_benchmarks details element_construction hex_overlap_volume
                tet_overlap_volume
)

# register the individual benchmarks
foreach(benchmark ${_benchmarks})
  add_executable(${benchmark} ${benchmark}.cpp)
  target_link_libraries(
    ${benchmark} PRIVATE overlap::headers Eigen3::Eigen nanobench
                         doctest::doctest
  )

  target_compile_definitions(
    ${benchmark} PRIVATE DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
  )

  add_test(NAME ${benchmark} COMMAND $<TARGET_FILE:${benchmark}>)

  set_tests_properties(${benchmark} PROPERTIES LABELS benchmark)
endforeach()
