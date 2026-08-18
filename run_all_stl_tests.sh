#!/bin/bash

APP=./build/input/stl/test_stl_validator
STL_DIR=tests/input_tests/stl_tests/stls

$APP $STL_DIR/valid_triangle.stl valid
$APP $STL_DIR/coincident_vertices.stl invalid
$APP $STL_DIR/nearly_coincident_vertices.stl invalid
$APP $STL_DIR/collinear_vertices.stl invalid
$APP $STL_DIR/nearly_collinear_vertices.stl invalid
$APP $STL_DIR/small_valid_triangle.stl valid