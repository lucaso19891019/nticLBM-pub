#!/bin/bash

APP=./build/input/stl/test_stl_validator
STL_DIR=tests/input_tests/stl_tests/stls


#==============================================================================
# 1. Degenerate facet tests
#==============================================================================

$APP $STL_DIR/valid_triangle.stl valid
$APP $STL_DIR/coincident_vertices.stl invalid
$APP $STL_DIR/nearly_coincident_vertices.stl invalid
$APP $STL_DIR/collinear_vertices.stl invalid
$APP $STL_DIR/nearly_collinear_vertices.stl invalid
$APP $STL_DIR/small_valid_triangle.stl valid


#==============================================================================
# 2. Facet normal tests
#==============================================================================

$APP $STL_DIR/valid_normal.stl valid
$APP $STL_DIR/non_unit_normal.stl valid
$APP $STL_DIR/zero_normal.stl invalid
$APP $STL_DIR/reversed_normal.stl invalid
$APP $STL_DIR/perpendicular_normal.stl invalid
$APP $STL_DIR/slightly_tilted_normal.stl valid

./build/input/stl/test_stl_normal_reconstruction \
    $STL_DIR/non_unit_normal.stl

#==============================================================================
# 3. Duplicate facet tests
#==============================================================================

$APP $STL_DIR/different_facets.stl valid
$APP $STL_DIR/duplicate_exact.stl invalid
$APP $STL_DIR/duplicate_cyclic_order.stl invalid
$APP $STL_DIR/duplicate_reversed_winding.stl invalid
$APP $STL_DIR/duplicate_with_small_noise.stl invalid
$APP $STL_DIR/close_but_distinct_facets.stl valid
