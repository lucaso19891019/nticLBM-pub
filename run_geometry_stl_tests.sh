#!/bin/bash

GEOMETRY_APP_DIR=./build/geometry/stl

STL_DIR=tests/geometry_tests/stl_tests/stls


#=============================================================================
# 1. STL geometry construction tests
#=============================================================================

$GEOMETRY_APP_DIR/test_stl_geometry \
    $STL_DIR/pipe_flat_binary.stl