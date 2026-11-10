cmake_minimum_required(VERSION 3.18)

project(nticLBM LANGUAGES C CXX)


#====================================================
# Dependencies
#====================================================

find_package(OpenMP REQUIRED)


#====================================================
# C++ standard
#====================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


#====================================================
# Compiler warnings
#====================================================

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(
        -Wall
        -Wextra
        -Wpedantic
    )
endif()


#====================================================
# Third-party libraries
#====================================================

include(
    cmake/3rdparties/thirdparties.cmake
)


#====================================================
# Tests
#====================================================

include(
    cmake/tests/input_config_tests.cmake
)

include(
    cmake/tests/input_stl_tests.cmake
)

include(
    cmake/tests/geometry_primitive_tests.cmake
)

include(
    cmake/tests/geometry_stl_tests.cmake
)

include(
    cmake/tests/geometry_boolean_tests.cmake
)

include(
    cmake/tests/geometry_multilevel_tests.cmake
)

include(
    cmake/tests/mpi_tests.cmake
)
