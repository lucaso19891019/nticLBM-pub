#!/bin/bash

APP=./build/input/stl/test_stl_validator
NORMAL_APP=./build/input/stl/test_stl_normal_reconstruction

STL_DIR=tests/input_tests/stl_tests/stls


#==============================================================================
# 1. Degenerate facet tests
#==============================================================================

$APP $STL_DIR/valid_triangle.stl valid test
$APP $STL_DIR/coincident_vertices.stl invalid test
$APP $STL_DIR/nearly_coincident_vertices.stl invalid test
$APP $STL_DIR/collinear_vertices.stl invalid test
$APP $STL_DIR/nearly_collinear_vertices.stl invalid test
$APP $STL_DIR/small_valid_triangle.stl valid test


#==============================================================================
# 2. Facet normal tests
#==============================================================================

$APP $STL_DIR/valid_normal.stl valid test
$APP $STL_DIR/non_unit_normal.stl valid test
$APP $STL_DIR/zero_normal.stl invalid test
$APP $STL_DIR/reversed_normal.stl invalid test
$APP $STL_DIR/perpendicular_normal.stl invalid test
$APP $STL_DIR/slightly_tilted_normal.stl valid test

$NORMAL_APP $STL_DIR/non_unit_normal.stl


#==============================================================================
# 3. Duplicate facet tests
#==============================================================================

$APP $STL_DIR/different_facets.stl valid test
$APP $STL_DIR/duplicate_exact.stl invalid test
$APP $STL_DIR/duplicate_cyclic_order.stl invalid test
$APP $STL_DIR/duplicate_reversed_winding.stl invalid test
$APP $STL_DIR/duplicate_with_small_noise.stl invalid test
$APP $STL_DIR/close_but_distinct_facets.stl valid test

#==============================================================================
# 4. Topology tests
#==============================================================================

$APP $STL_DIR/closed_tetrahedron.stl valid full
$APP $STL_DIR/open_tetrahedron.stl invalid full
$APP $STL_DIR/non_manifold_edge.stl invalid full
$APP $STL_DIR/inconsistent_winding.stl invalid full
$APP $STL_DIR/two_closed_tetrahedra.stl valid full

#=============================================================================
# 5. Component geometry tests
#=============================================================================

$APP $STL_DIR/large_tetrahedron.stl valid full
