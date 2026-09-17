#!/bin/bash

APP=./build/input/stl/test_stl_validator
NORMAL_APP=./build/input/stl/test_stl_normal_reconstruction

STL_DIR=tests/input_tests/stl_tests/stls

#=============================================================================
# STL reader tests
#=============================================================================

$APP $STL_DIR/valid_binary.stl valid test

#==============================================================================
# 1. Degenerate facet tests
#==============================================================================

$APP $STL_DIR/valid_triangle.stl valid test
$APP $STL_DIR/coincident_vertices.stl invalid test
$APP $STL_DIR/nearly_coincident_vertices.stl invalid test
$APP $STL_DIR/collinear_vertices.stl invalid test
$APP $STL_DIR/nearly_collinear_vertices.stl invalid test
$APP $STL_DIR/small_valid_triangle.stl valid test


#=============================================================================
# 2. Facet normal reconstruction tests
#==============================================================================

$APP $STL_DIR/valid_normal.stl valid test
$APP $STL_DIR/non_unit_normal.stl valid test
$APP $STL_DIR/zero_normal.stl valid test
$APP $STL_DIR/reversed_normal.stl valid test
$APP $STL_DIR/perpendicular_normal.stl valid test
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
$APP $STL_DIR/zero_volume_closed_component.stl invalid full

#=============================================================================
# 6. Component intersections
#=============================================================================

$APP $STL_DIR/two_separated_tetrahedra.stl valid full
$APP $STL_DIR/nested_tetrahedra.stl valid full
$APP $STL_DIR/intersecting_tetrahedra.stl invalid full
$APP $STL_DIR/vertex_touching_tetrahedra.stl invalid full
$APP $STL_DIR/edge_touching_tetrahedra.stl invalid full
$APP $STL_DIR/face_touching_tetrahedra.stl invalid full
$APP $STL_DIR/near_but_not_touching_tetrahedra.stl valid full
$APP $STL_DIR/aabb_overlap_separated.stl valid full
$APP $STL_DIR/small_component_intersection.stl invalid full
$APP $STL_DIR/large_facet_intersection.stl invalid full
$APP $STL_DIR/coplanar_partial_overlap.stl invalid full
$APP $STL_DIR/coplanar_edge_crossing.stl invalid full
