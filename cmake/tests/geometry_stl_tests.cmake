#====================================================
# Geometry STL test output
#====================================================

set(GEOMETRY_STL_TEST_OUTPUT_DIR
    "${CMAKE_BINARY_DIR}/geometry/stl"
)


function(set_geometry_stl_test_output TARGET_NAME)
    set_target_properties(
        ${TARGET_NAME}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
            "${GEOMETRY_STL_TEST_OUTPUT_DIR}"
    )
endfunction()


#====================================================
# test_stl_geometry
#====================================================

add_executable(test_stl_geometry
    tests/geometry_tests/stl_tests/test_stl_geometry.cpp

    src/geometry/stl/stl_geometry.cpp

    src/input/stl/stl_reader.cpp
    src/input/stl/stl_validator.cpp
)

target_include_directories(
    test_stl_geometry
    PRIVATE
    src/geometry
    src/geometry/stl
    src/input/stl
    src/common
)

target_link_libraries(
    test_stl_geometry
    PRIVATE
    OpenMP::OpenMP_CXX
)

set_geometry_stl_test_output(
    test_stl_geometry
)