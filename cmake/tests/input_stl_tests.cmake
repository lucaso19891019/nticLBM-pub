#====================================================
# Input STL test output
#====================================================

set(INPUT_STL_TEST_OUTPUT_DIR
    "${CMAKE_BINARY_DIR}/input/stl"
)


function(set_input_stl_test_output TARGET_NAME)
    set_target_properties(
        ${TARGET_NAME}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
            "${INPUT_STL_TEST_OUTPUT_DIR}"
    )
endfunction()


#====================================================
# test_stl_validator
#====================================================

add_executable(test_stl_validator
    tests/input_tests/stl_tests/test_stl_validator.cpp

    src/input/stl/stl_reader.cpp
    src/input/stl/stl_validator.cpp
)

target_include_directories(
    test_stl_validator
    PRIVATE
    src/input/stl
    src/common
)

target_link_libraries(
    test_stl_validator
    PRIVATE
    OpenMP::OpenMP_CXX
)

set_input_stl_test_output(
    test_stl_validator
)


#====================================================
# test_stl_normal_reconstruction
#====================================================

add_executable(test_stl_normal_reconstruction
    tests/input_tests/stl_tests/test_stl_normal_reconstruction.cpp

    src/input/stl/stl_reader.cpp
    src/input/stl/stl_validator.cpp
)

target_include_directories(
    test_stl_normal_reconstruction
    PRIVATE
    src/input/stl
    src/common
)

target_link_libraries(
    test_stl_normal_reconstruction
    PRIVATE
    OpenMP::OpenMP_CXX
)

set_input_stl_test_output(
    test_stl_normal_reconstruction
)