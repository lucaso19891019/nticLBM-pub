#!/bin/bash

GEOMETRY_APP_DIR=./build/geometry/stl

STL_DIR=tests/geometry_tests/stl_tests/stls


#=============================================================================
# 1. STL geometry construction tests
#=============================================================================

$GEOMETRY_APP_DIR/test_stl_geometry \
    $STL_DIR/smooth_irregular_branched_channel.stl

#=============================================================================
# 2. STL component containment
#=============================================================================

$GEOMETRY_APP_DIR/test_stl_containment \
    single_root \
    $STL_DIR/single_root.stl

$GEOMETRY_APP_DIR/test_stl_containment \
    two_roots \
    $STL_DIR/two_roots.stl

$GEOMETRY_APP_DIR/test_stl_containment \
    aabb_overlap \
    $STL_DIR/aabb_overlap.stl

$GEOMETRY_APP_DIR/test_stl_containment \
    nested \
    $STL_DIR/nested.stl

$GEOMETRY_APP_DIR/test_stl_containment \
    pseudo_containment \
    $STL_DIR/pseudo_containment.stl

$GEOMETRY_APP_DIR/test_stl_containment \
    two_children \
    $STL_DIR/two_children.stl

$GEOMETRY_APP_DIR/test_stl_containment \
    two_nested_roots \
    $STL_DIR/two_nested_roots.stl

$GEOMETRY_APP_DIR/test_stl_containment \
    three_levels \
    $STL_DIR/three_levels.stl

#=============================================================================
# 3. STL flow interpreter
#=============================================================================

$GEOMETRY_APP_DIR/test_stl_flow \
    internal_root \
    $STL_DIR/single_root.stl

$GEOMETRY_APP_DIR/test_stl_flow \
    internal_nested \
    $STL_DIR/nested.stl

$GEOMETRY_APP_DIR/test_stl_flow \
    internal_nested \
    $STL_DIR/two_children.stl

$GEOMETRY_APP_DIR/test_stl_flow \
    internal_multiple_roots \
    $STL_DIR/two_roots.stl

$GEOMETRY_APP_DIR/test_stl_flow \
    external_roots \
    $STL_DIR/two_roots.stl

$GEOMETRY_APP_DIR/test_stl_flow \
    external_nested \
    $STL_DIR/two_nested_roots.stl

$GEOMETRY_APP_DIR/test_stl_flow \
    internal_to_external \
    $STL_DIR/nested.stl

#=============================================================================
# 4. STL translation tests
#=============================================================================

$GEOMETRY_APP_DIR/test_stl_translation \
    general \
    $STL_DIR/nested.stl

$GEOMETRY_APP_DIR/test_stl_translation \
    internal \
    $STL_DIR/nested.stl

$GEOMETRY_APP_DIR/test_stl_translation \
    external \
    $STL_DIR/two_nested_roots.stl

$GEOMETRY_APP_DIR/test_stl_translation \
    repeated \
    $STL_DIR/two_children.stl