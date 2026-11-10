#====================================================
# MPI partition test paths
#====================================================

set(MPI_TEST_OUTPUT_DIR
    "${CMAKE_BINARY_DIR}/mpi"
)


#====================================================
# Common include directories
#====================================================

set(MPI_TEST_INCLUDE_DIRS
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}/src/geometry"
    "${CMAKE_SOURCE_DIR}/src/geometry/stl"
    "${CMAKE_SOURCE_DIR}/src/lattice"
    "${CMAKE_SOURCE_DIR}/src/input/stl"
    "${CMAKE_SOURCE_DIR}/src/common"
)


#====================================================
# Shared partition implementation for MPI tests
#====================================================

add_library(mpi_partition_test_support STATIC
    src/mpi/metis_partition.cpp
)

target_include_directories(
    mpi_partition_test_support
    PUBLIC
    ${MPI_TEST_INCLUDE_DIRS}
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
    ARCHIVE_OUTPUT_DIRECTORY "${MPI_TEST_OUTPUT_DIR}"
)


#====================================================
# test_partition_ownership
#====================================================

add_executable(test_partition_ownership
    tests/mpi_tests/test_partition_ownership.cpp
)

target_link_libraries(
    test_partition_ownership
    PRIVATE
    mpi_partition_test_support
)

set_target_properties(
    test_partition_ownership
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${MPI_TEST_OUTPUT_DIR}"
)


#====================================================
# test_partition_connectivity
#====================================================

add_executable(test_partition_connectivity
    tests/mpi_tests/test_partition_connectivity.cpp
)

target_link_libraries(
    test_partition_connectivity
    PRIVATE
    mpi_partition_test_support
)

set_target_properties(
    test_partition_connectivity
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${MPI_TEST_OUTPUT_DIR}"
)


#====================================================
# test_partition_validation
#====================================================

add_executable(test_partition_validation
    tests/mpi_tests/test_partition_validation.cpp
)

target_link_libraries(
    test_partition_validation
    PRIVATE
    mpi_partition_test_support
)

set_target_properties(
    test_partition_validation
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${MPI_TEST_OUTPUT_DIR}"
)


#====================================================
# test_partition_blocks
#====================================================

add_executable(test_partition_blocks
    tests/mpi_tests/test_partition_blocks.cpp
)

target_link_libraries(
    test_partition_blocks
    PRIVATE
    mpi_partition_test_support
)

set_target_properties(
    test_partition_blocks
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${MPI_TEST_OUTPUT_DIR}"
)
