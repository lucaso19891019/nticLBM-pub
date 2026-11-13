#====================================================
# Geometry Boolean paths
#====================================================

set(GEOMETRY_BOOLEAN_TEST_OUTPUT_DIR
    "${CMAKE_BINARY_DIR}/geometry/boolean"
)

set(GEOMETRY_TEST_OUTPUT_DIR
    "${CMAKE_SOURCE_DIR}/tests/geometry_tests/outputs"
)

set(GEOMETRY_TEST_SOURCE_DIR
    "${CMAKE_SOURCE_DIR}/tests/geometry_tests"
)


#====================================================
# Geometry Boolean test configuration
#====================================================

function(configure_geometry_boolean_test TARGET_NAME)
    set_target_properties(
        ${TARGET_NAME}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
            "${GEOMETRY_BOOLEAN_TEST_OUTPUT_DIR}"
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

set(GEOMETRY_BOOLEAN_INCLUDE_DIRS
    "${CMAKE_SOURCE_DIR}/src/geometry"
    "${CMAKE_SOURCE_DIR}/src/geometry/primitive"
    "${CMAKE_SOURCE_DIR}/src/lattice"
    "${CMAKE_SOURCE_DIR}/src/input/stl"
    "${CMAKE_SOURCE_DIR}/src/geometry/stl"
    "${CMAKE_SOURCE_DIR}/src/common"
    "${CMAKE_SOURCE_DIR}/tests/geometry_tests/outputs"
)


#====================================================
# test_boolean_annulus
#====================================================

add_executable(test_boolean_annulus
    tests/geometry_tests/boolean_tests/test_boolean_annulus.cpp
)

target_include_directories(
    test_boolean_annulus
    PRIVATE
    ${GEOMETRY_BOOLEAN_INCLUDE_DIRS}
)

target_link_libraries(
    test_boolean_annulus
    PRIVATE
    OpenMP::OpenMP_CXX
)

configure_geometry_boolean_test(
    test_boolean_annulus
)


#====================================================
# test_boolean_annular_pipe
#====================================================

add_executable(test_boolean_annular_pipe
    tests/geometry_tests/boolean_tests/test_boolean_annular_pipe.cpp
)

target_include_directories(
    test_boolean_annular_pipe
    PRIVATE
    ${GEOMETRY_BOOLEAN_INCLUDE_DIRS}
)

target_link_libraries(
    test_boolean_annular_pipe
    PRIVATE
    OpenMP::OpenMP_CXX
)

configure_geometry_boolean_test(
    test_boolean_annular_pipe
)


#====================================================
# test_boolean_branched_channel
#====================================================

add_executable(test_boolean_branched_channel
    tests/geometry_tests/boolean_tests/test_boolean_branched_channel.cpp

    src/geometry/stl/stl_geometry.cpp
    src/geometry/stl/stl_containment.cpp
    src/geometry/stl/stl_flow.cpp
    src/geometry/stl/stl_translation.cpp
    src/geometry/stl/stl_area_analysis.cpp
    src/geometry/stl/stl_boundary_features.cpp

    src/input/stl/stl_reader.cpp
    src/input/stl/stl_validator.cpp
)

target_include_directories(
    test_boolean_branched_channel
    PRIVATE
    ${GEOMETRY_BOOLEAN_INCLUDE_DIRS}
)

target_link_libraries(
    test_boolean_branched_channel
    PRIVATE
    OpenMP::OpenMP_CXX
)

configure_geometry_boolean_test(
    test_boolean_branched_channel
)


#====================================================
# test_boolean_crossed_cylinders
#====================================================

add_executable(test_boolean_crossed_cylinders
    tests/geometry_tests/boolean_tests/test_boolean_crossed_cylinders.cpp
)

target_include_directories(
    test_boolean_crossed_cylinders
    PRIVATE
    ${GEOMETRY_BOOLEAN_INCLUDE_DIRS}
)

target_link_libraries(
    test_boolean_crossed_cylinders
    PRIVATE
    OpenMP::OpenMP_CXX
)

configure_geometry_boolean_test(
    test_boolean_crossed_cylinders
)


#====================================================
# test_boolean_two_spheres
#====================================================

add_executable(test_boolean_two_spheres
    tests/geometry_tests/boolean_tests/test_boolean_two_spheres.cpp
)

target_include_directories(
    test_boolean_two_spheres
    PRIVATE
    ${GEOMETRY_BOOLEAN_INCLUDE_DIRS}
)

target_link_libraries(
    test_boolean_two_spheres
    PRIVATE
    OpenMP::OpenMP_CXX
)

configure_geometry_boolean_test(
    test_boolean_two_spheres
)
