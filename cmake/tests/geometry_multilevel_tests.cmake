#====================================================
# Geometry multilevel paths
#====================================================

set(GEOMETRY_MULTILEVEL_TEST_OUTPUT_DIR
    "${CMAKE_BINARY_DIR}/geometry/multilevel"
)

set(GEOMETRY_TEST_OUTPUT_DIR
    "${CMAKE_SOURCE_DIR}/tests/geometry_tests/outputs"
)

set(GEOMETRY_TEST_SOURCE_DIR
    "${CMAKE_SOURCE_DIR}/tests/geometry_tests"
)


#====================================================
# Geometry multilevel test configuration
#====================================================

function(configure_geometry_multilevel_test TARGET_NAME)
    set_target_properties(
        ${TARGET_NAME}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
            "${GEOMETRY_MULTILEVEL_TEST_OUTPUT_DIR}"
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

set(GEOMETRY_MULTILEVEL_INCLUDE_DIRS
    "${CMAKE_SOURCE_DIR}/src/geometry"
    "${CMAKE_SOURCE_DIR}/src/geometry/stl"
    "${CMAKE_SOURCE_DIR}/src/lattice"
    "${CMAKE_SOURCE_DIR}/src/input/stl"
    "${CMAKE_SOURCE_DIR}/src/common"
    "${CMAKE_SOURCE_DIR}/tests/geometry_tests/outputs"
)

#====================================================
# Shared STL implementation for multilevel tests
#====================================================

add_library(geometry_multilevel_test_support STATIC
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
    geometry_multilevel_test_support
    PUBLIC
    ${GEOMETRY_MULTILEVEL_INCLUDE_DIRS}
)

target_link_libraries(
    geometry_multilevel_test_support
    PUBLIC
    OpenMP::OpenMP_CXX
)


#====================================================
# test_multilevel_annulus
#====================================================

add_executable(test_multilevel_annulus
    tests/geometry_tests/multilevel_tests/test_multilevel_annulus.cpp
)

target_include_directories(
    test_multilevel_annulus
    PRIVATE
    ${GEOMETRY_MULTILEVEL_INCLUDE_DIRS}
)

target_link_libraries(
    test_multilevel_annulus
    PRIVATE
    geometry_multilevel_test_support
)

configure_geometry_multilevel_test(
    test_multilevel_annulus
)


#====================================================
# test_multilevel_spherical_obstacle
#====================================================

add_executable(test_multilevel_spherical_obstacle
    tests/geometry_tests/multilevel_tests/test_multilevel_spherical_obstacle.cpp
)

target_include_directories(
    test_multilevel_spherical_obstacle
    PRIVATE
    ${GEOMETRY_MULTILEVEL_INCLUDE_DIRS}
)

target_link_libraries(
    test_multilevel_spherical_obstacle
    PRIVATE
    geometry_multilevel_test_support
)

configure_geometry_multilevel_test(
    test_multilevel_spherical_obstacle
)


#====================================================
# test_multilevel_branched_channel
#====================================================

add_executable(test_multilevel_branched_channel
    tests/geometry_tests/multilevel_tests/test_multilevel_branched_channel.cpp
)

target_include_directories(
    test_multilevel_branched_channel
    PRIVATE
    ${GEOMETRY_MULTILEVEL_INCLUDE_DIRS}
)

target_link_libraries(
    test_multilevel_branched_channel
    PRIVATE
    geometry_multilevel_test_support
)

configure_geometry_multilevel_test(
    test_multilevel_branched_channel
)


#====================================================
# test_multilevel_topology
#====================================================

add_executable(test_multilevel_topology
    tests/geometry_tests/multilevel_tests/test_multilevel_topology.cpp
)

target_include_directories(
    test_multilevel_topology
    PRIVATE
    ${GEOMETRY_MULTILEVEL_INCLUDE_DIRS}
)

target_link_libraries(
    test_multilevel_topology
    PRIVATE
    geometry_multilevel_test_support
)

configure_geometry_multilevel_test(
    test_multilevel_topology
)


#====================================================
# test_multilevel_regions
#====================================================

add_executable(test_multilevel_regions
    tests/geometry_tests/multilevel_tests/test_multilevel_regions.cpp
)

target_include_directories(
    test_multilevel_regions
    PRIVATE
    ${GEOMETRY_MULTILEVEL_INCLUDE_DIRS}
)

target_link_libraries(
    test_multilevel_regions
    PRIVATE
    geometry_multilevel_test_support
)

configure_geometry_multilevel_test(
    test_multilevel_regions
)


#====================================================
# test_multilevel_primitive
#====================================================

add_executable(test_multilevel_primitive
    tests/geometry_tests/multilevel_tests/test_multilevel_primitive.cpp
)

target_include_directories(
    test_multilevel_primitive
    PRIVATE
    ${GEOMETRY_MULTILEVEL_INCLUDE_DIRS}
)

target_link_libraries(
    test_multilevel_primitive
    PRIVATE
    geometry_multilevel_test_support
)

configure_geometry_multilevel_test(
    test_multilevel_primitive
)


#====================================================
# test_multilevel_edits
#====================================================

add_executable(test_multilevel_edits
    tests/geometry_tests/multilevel_tests/test_multilevel_edits.cpp
)

target_include_directories(
    test_multilevel_edits
    PRIVATE
    ${GEOMETRY_MULTILEVEL_INCLUDE_DIRS}
)

target_link_libraries(
    test_multilevel_edits
    PRIVATE
    geometry_multilevel_test_support
)

configure_geometry_multilevel_test(
    test_multilevel_edits
)


#====================================================
# test_multilevel_lattice
#====================================================

add_executable(test_multilevel_lattice
    tests/geometry_tests/multilevel_tests/test_multilevel_lattice.cpp
)

target_include_directories(
    test_multilevel_lattice
    PRIVATE
    ${GEOMETRY_MULTILEVEL_INCLUDE_DIRS}
)

target_link_libraries(
    test_multilevel_lattice
    PRIVATE
    geometry_multilevel_test_support
)

configure_geometry_multilevel_test(
    test_multilevel_lattice
)


#====================================================
# test_multilevel_stl
#====================================================

add_executable(test_multilevel_stl
    tests/geometry_tests/multilevel_tests/test_multilevel_stl.cpp
)

target_include_directories(
    test_multilevel_stl
    PRIVATE
    ${GEOMETRY_MULTILEVEL_INCLUDE_DIRS}
)

target_link_libraries(
    test_multilevel_stl
    PRIVATE
    geometry_multilevel_test_support
)

configure_geometry_multilevel_test(
    test_multilevel_stl
)


#====================================================
# test_multilevel_determinism
#====================================================

add_executable(test_multilevel_determinism
    tests/geometry_tests/multilevel_tests/test_multilevel_determinism.cpp
)

target_include_directories(
    test_multilevel_determinism
    PRIVATE
    ${GEOMETRY_MULTILEVEL_INCLUDE_DIRS}
)

target_link_libraries(
    test_multilevel_determinism
    PRIVATE
    geometry_multilevel_test_support
)

configure_geometry_multilevel_test(
    test_multilevel_determinism
)


#====================================================
# Multilevel acceptance build target
#====================================================

add_custom_target(test_multilevel_acceptance
    DEPENDS
    test_multilevel_topology
    test_multilevel_regions
    test_multilevel_primitive
    test_multilevel_edits
    test_multilevel_lattice
    test_multilevel_stl
    test_multilevel_determinism
)
