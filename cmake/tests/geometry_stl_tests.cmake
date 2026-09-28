#====================================================
# Geometry STL paths
#====================================================

set(GEOMETRY_STL_TEST_OUTPUT_DIR
    "${CMAKE_BINARY_DIR}/geometry/stl"
)

set(GEOMETRY_TEST_OUTPUT_DIR
    "${CMAKE_SOURCE_DIR}/tests/geometry_tests/outputs"
)

set(GEOMETRY_TEST_SOURCE_DIR
    "${CMAKE_SOURCE_DIR}/tests/geometry_tests"
)


#====================================================
# Geometry STL test configuration
#====================================================

function(configure_geometry_stl_test TARGET_NAME)
    set_target_properties(
        ${TARGET_NAME}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
            "${GEOMETRY_STL_TEST_OUTPUT_DIR}"
    )

    target_compile_definitions(
        ${TARGET_NAME}
        PRIVATE
        "GEOMETRY_TEST_OUTPUT_DIR=\"${GEOMETRY_TEST_OUTPUT_DIR}\""
        "GEOMETRY_TEST_SOURCE_DIR=\"${GEOMETRY_TEST_SOURCE_DIR}\""
    )
endfunction()


#====================================================
# Common include directories
#====================================================

set(GEOMETRY_STL_INCLUDE_DIRS
    src/geometry
    src/geometry/stl
    src/input/stl
    src/common
    tests/geometry_tests/outputs
)


#====================================================
# test_stl_geometry
#====================================================

add_executable(test_stl_geometry
    tests/geometry_tests/stl_tests/test_stl_geometry.cpp

    src/geometry/stl/stl_geometry.cpp
    src/geometry/stl/stl_containment.cpp
    src/geometry/stl/stl_flow.cpp
    src/geometry/stl/stl_translation.cpp
    src/geometry/stl/stl_contains.cpp

    src/input/stl/stl_reader.cpp
    src/input/stl/stl_validator.cpp
)

target_include_directories(
    test_stl_geometry
    PRIVATE
    ${GEOMETRY_STL_INCLUDE_DIRS}
)

target_link_libraries(
    test_stl_geometry
    PRIVATE
    OpenMP::OpenMP_CXX
)

configure_geometry_stl_test(
    test_stl_geometry
)


#====================================================
# test_stl_containment
#====================================================

add_executable(test_stl_containment
    tests/geometry_tests/stl_tests/test_stl_containment.cpp

    src/geometry/stl/stl_geometry.cpp
    src/geometry/stl/stl_containment.cpp
    src/geometry/stl/stl_flow.cpp
    src/geometry/stl/stl_translation.cpp
    src/geometry/stl/stl_contains.cpp

    src/input/stl/stl_reader.cpp
    src/input/stl/stl_validator.cpp
)

target_include_directories(
    test_stl_containment
    PRIVATE
    ${GEOMETRY_STL_INCLUDE_DIRS}
)

target_link_libraries(
    test_stl_containment
    PRIVATE
    OpenMP::OpenMP_CXX
)

configure_geometry_stl_test(
    test_stl_containment
)


#====================================================
# test_stl_flow
#====================================================

add_executable(test_stl_flow
    tests/geometry_tests/stl_tests/test_stl_flow.cpp

    src/geometry/stl/stl_geometry.cpp
    src/geometry/stl/stl_containment.cpp
    src/geometry/stl/stl_flow.cpp
    src/geometry/stl/stl_translation.cpp
    src/geometry/stl/stl_contains.cpp

    src/input/stl/stl_reader.cpp
    src/input/stl/stl_validator.cpp
)

target_include_directories(
    test_stl_flow
    PRIVATE
    ${GEOMETRY_STL_INCLUDE_DIRS}
)

target_link_libraries(
    test_stl_flow
    PRIVATE
    OpenMP::OpenMP_CXX
)

configure_geometry_stl_test(
    test_stl_flow
)


#====================================================
# test_stl_translation
#====================================================

add_executable(test_stl_translation
    tests/geometry_tests/stl_tests/test_stl_translation.cpp

    src/geometry/stl/stl_geometry.cpp
    src/geometry/stl/stl_containment.cpp
    src/geometry/stl/stl_flow.cpp
    src/geometry/stl/stl_translation.cpp
    src/geometry/stl/stl_contains.cpp

    src/input/stl/stl_reader.cpp
    src/input/stl/stl_validator.cpp
)

target_include_directories(
    test_stl_translation
    PRIVATE
    ${GEOMETRY_STL_INCLUDE_DIRS}
)

target_link_libraries(
    test_stl_translation
    PRIVATE
    OpenMP::OpenMP_CXX
)

configure_geometry_stl_test(
    test_stl_translation
)


#====================================================
# test_stl_area_analysis
#====================================================

add_executable(test_stl_area_analysis
    tests/geometry_tests/stl_tests/test_stl_area_analysis.cpp

    src/geometry/stl/stl_geometry.cpp
    src/geometry/stl/stl_containment.cpp
    src/geometry/stl/stl_flow.cpp
    src/geometry/stl/stl_translation.cpp
    src/geometry/stl/stl_contains.cpp

    src/input/stl/stl_reader.cpp
    src/input/stl/stl_validator.cpp
)

target_include_directories(
    test_stl_area_analysis
    PRIVATE
    ${GEOMETRY_STL_INCLUDE_DIRS}
)

target_link_libraries(
    test_stl_area_analysis
    PRIVATE
    OpenMP::OpenMP_CXX
)

configure_geometry_stl_test(
    test_stl_area_analysis
)