#====================================================
# Geometry primitive test output
#====================================================

set(GEOMETRY_PRIMITIVE_TEST_OUTPUT_DIR
    "${CMAKE_BINARY_DIR}/geometry/primitive"
)


function(set_geometry_primitive_test_output TARGET_NAME)
    set_target_properties(
        ${TARGET_NAME}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
            "${GEOMETRY_PRIMITIVE_TEST_OUTPUT_DIR}"
    )
endfunction()


#====================================================
# Common include directories
#====================================================

set(GEOMETRY_PRIMITIVE_INCLUDE_DIRS
    src/geometry
    tests/geometry_tests/outputs
    tests/geometry_tests/primitive_tests
)


#====================================================
# test_rectangle
#====================================================

add_executable(test_rectangle
    tests/geometry_tests/primitive_tests/test_rectangle.cpp
)

target_include_directories(
    test_rectangle
    PRIVATE
    ${GEOMETRY_PRIMITIVE_INCLUDE_DIRS}
)

target_link_libraries(
    test_rectangle
    PRIVATE
    OpenMP::OpenMP_CXX
)

set_geometry_primitive_test_output(
    test_rectangle
)


#====================================================
# test_circle
#====================================================

add_executable(test_circle
    tests/geometry_tests/primitive_tests/test_circle.cpp
)

target_include_directories(
    test_circle
    PRIVATE
    ${GEOMETRY_PRIMITIVE_INCLUDE_DIRS}
)

target_link_libraries(
    test_circle
    PRIVATE
    OpenMP::OpenMP_CXX
)

set_geometry_primitive_test_output(
    test_circle
)


#====================================================
# test_box
#====================================================

add_executable(test_box
    tests/geometry_tests/primitive_tests/test_box.cpp
)

target_include_directories(
    test_box
    PRIVATE
    ${GEOMETRY_PRIMITIVE_INCLUDE_DIRS}
)

target_link_libraries(
    test_box
    PRIVATE
    OpenMP::OpenMP_CXX
)

set_geometry_primitive_test_output(
    test_box
)


#====================================================
# test_sphere
#====================================================

add_executable(test_sphere
    tests/geometry_tests/primitive_tests/test_sphere.cpp
)

target_include_directories(
    test_sphere
    PRIVATE
    ${GEOMETRY_PRIMITIVE_INCLUDE_DIRS}
)

target_link_libraries(
    test_sphere
    PRIVATE
    OpenMP::OpenMP_CXX
)

set_geometry_primitive_test_output(
    test_sphere
)


#====================================================
# test_cylinder
#====================================================

add_executable(test_cylinder
    tests/geometry_tests/primitive_tests/test_cylinder.cpp
)

target_include_directories(
    test_cylinder
    PRIVATE
    ${GEOMETRY_PRIMITIVE_INCLUDE_DIRS}
)

target_link_libraries(
    test_cylinder
    PRIVATE
    OpenMP::OpenMP_CXX
)

set_geometry_primitive_test_output(
    test_cylinder
)