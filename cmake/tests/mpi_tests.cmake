#====================================================
# MPI partition test paths
#====================================================

set(MPI_PARTITION_TEST_OUTPUT_DIR
    "${CMAKE_BINARY_DIR}/mpi/partition"
)

set(MPI_TEST_OUTPUT_DIR
    "${CMAKE_SOURCE_DIR}/tests/mpi_tests/outputs"
)

set(MPI_TEST_SOURCE_DIR
    "${CMAKE_SOURCE_DIR}/tests/mpi_tests"
)


#====================================================
# Manual partition test configuration
#====================================================

function(configure_mpi_partition_manual_test TARGET_NAME)
    target_include_directories(
        ${TARGET_NAME}
        PRIVATE
        "${CMAKE_SOURCE_DIR}/tests/mpi_tests/outputs"
    )

    target_compile_definitions(
        ${TARGET_NAME}
        PRIVATE
        "MPI_TEST_OUTPUT_DIR=\"${MPI_TEST_OUTPUT_DIR}\""
        "MPI_TEST_SOURCE_DIR=\"${MPI_TEST_SOURCE_DIR}\""
    )

    set_target_properties(
        ${TARGET_NAME}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${MPI_PARTITION_TEST_OUTPUT_DIR}"
    )
endfunction()


#====================================================
# Common include directories
#====================================================

set(MPI_PARTITION_TEST_INCLUDE_DIRS
    "${CMAKE_SOURCE_DIR}/src/mpi/partition"
    "${CMAKE_SOURCE_DIR}/src/geometry"
    "${CMAKE_SOURCE_DIR}/src/geometry/primitive"
    "${CMAKE_SOURCE_DIR}/src/geometry/multilevel"
    "${CMAKE_SOURCE_DIR}/src/geometry/stl"
    "${CMAKE_SOURCE_DIR}/src/lattice"
    "${CMAKE_SOURCE_DIR}/src/input/stl"
    "${CMAKE_SOURCE_DIR}/src/common"
)


#====================================================
# Shared partition implementation for MPI tests
#====================================================

add_library(mpi_partition_test_support STATIC
    src/mpi/partition/metis_partition.cpp
)

target_include_directories(
    mpi_partition_test_support
    PUBLIC
    ${MPI_PARTITION_TEST_INCLUDE_DIRS}
)

target_link_libraries(
    mpi_partition_test_support
    PUBLIC
    OpenMP::OpenMP_CXX
    PRIVATE
    METIS::METIS
)

set_target_properties(
    mpi_partition_test_support
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${MPI_PARTITION_TEST_OUTPUT_DIR}"
)


#====================================================
# STL implementation for manual partition tests
#====================================================

add_library(mpi_partition_stl_test_support STATIC
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
    mpi_partition_stl_test_support
    PUBLIC
    ${MPI_PARTITION_TEST_INCLUDE_DIRS}
)

target_link_libraries(
    mpi_partition_stl_test_support
    PUBLIC
    OpenMP::OpenMP_CXX
)

set_target_properties(
    mpi_partition_stl_test_support
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${MPI_PARTITION_TEST_OUTPUT_DIR}"
)


#====================================================
# test_partition_ownership
#====================================================

add_executable(test_partition_ownership
    tests/mpi_tests/partition_tests/test_partition_ownership.cpp
)

target_link_libraries(
    test_partition_ownership
    PRIVATE
    mpi_partition_test_support
)

set_target_properties(
    test_partition_ownership
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${MPI_PARTITION_TEST_OUTPUT_DIR}"
)


#====================================================
# test_partition_connectivity
#====================================================

add_executable(test_partition_connectivity
    tests/mpi_tests/partition_tests/test_partition_connectivity.cpp
)

target_link_libraries(
    test_partition_connectivity
    PRIVATE
    mpi_partition_test_support
)

set_target_properties(
    test_partition_connectivity
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${MPI_PARTITION_TEST_OUTPUT_DIR}"
)


#====================================================
# test_partition_validation
#====================================================

add_executable(test_partition_validation
    tests/mpi_tests/partition_tests/test_partition_validation.cpp
)

target_link_libraries(
    test_partition_validation
    PRIVATE
    mpi_partition_test_support
)

set_target_properties(
    test_partition_validation
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${MPI_PARTITION_TEST_OUTPUT_DIR}"
)


#====================================================
# test_partition_blocks
#====================================================

add_executable(test_partition_blocks
    tests/mpi_tests/partition_tests/test_partition_blocks.cpp
)

target_link_libraries(
    test_partition_blocks
    PRIVATE
    mpi_partition_test_support
)

set_target_properties(
    test_partition_blocks
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${MPI_PARTITION_TEST_OUTPUT_DIR}"
)


#====================================================
# test_partition_boundary_costs
#====================================================

add_executable(test_partition_boundary_costs
    tests/mpi_tests/partition_tests/test_partition_boundary_costs.cpp
)

target_link_libraries(
    test_partition_boundary_costs
    PRIVATE
    mpi_partition_test_support
)

set_target_properties(
    test_partition_boundary_costs
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${MPI_PARTITION_TEST_OUTPUT_DIR}"
)


#====================================================
# test_partition_intra_rank
#====================================================

add_executable(test_partition_intra_rank
    tests/mpi_tests/partition_tests/test_partition_intra_rank.cpp
)

target_link_libraries(
    test_partition_intra_rank
    PRIVATE
    mpi_partition_test_support
)

set_target_properties(
    test_partition_intra_rank
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${MPI_PARTITION_TEST_OUTPUT_DIR}"
)


#====================================================
# test_partition_taylor_couette (interactive)
#====================================================

add_executable(test_partition_taylor_couette
    tests/mpi_tests/partition_tests/test_partition_taylor_couette.cpp
)

target_link_libraries(
    test_partition_taylor_couette
    PRIVATE
    mpi_partition_test_support
)

configure_mpi_partition_manual_test(
    test_partition_taylor_couette
)


#====================================================
# test_partition_branched_channel (interactive)
#====================================================

add_executable(test_partition_branched_channel
    tests/mpi_tests/partition_tests/test_partition_branched_channel.cpp
)

target_link_libraries(
    test_partition_branched_channel
    PRIVATE
    mpi_partition_test_support
    mpi_partition_stl_test_support
)

configure_mpi_partition_manual_test(
    test_partition_branched_channel
)
