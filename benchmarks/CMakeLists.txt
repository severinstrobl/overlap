if(NOT TARGET doctest::doctest)
  add_subdirectory(
    "${PROJECT_SOURCE_DIR}/third_party/doctest" "third_party/doctest"
    EXCLUDE_FROM_ALL
  )

  list(APPEND CMAKE_MODULE_PATH
       "${PROJECT_SOURCE_DIR}/third_party/doctest/scripts/cmake"
  )
endif()

add_subdirectory(
  "${PROJECT_SOURCE_DIR}/third_party/nanobench" "third_party/nanobench"
  EXCLUDE_FROM_ALL
)

# list of benchmarks
set(benchmarks details hex_overlap tet_overlap)

# register the individual benchmarks
foreach(benchmark ${benchmarks})
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
