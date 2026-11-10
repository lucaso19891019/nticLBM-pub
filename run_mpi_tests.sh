#!/bin/bash

MPI_APP_DIR=./build/mpi


#=============================================================================
# 1. Uniform and weighted multilevel ownership
#=============================================================================

$MPI_APP_DIR/test_partition_ownership


#=============================================================================
# 2. Ghosts, peer lists, interface indices and boundary Q
#=============================================================================

$MPI_APP_DIR/test_partition_connectivity


#=============================================================================
# 3. Determinism, invalid input and empty ranks
#=============================================================================

$MPI_APP_DIR/test_partition_validation


#=============================================================================
# 4. Rectangular coverage, level frequencies and hardware capacities
#=============================================================================

$MPI_APP_DIR/test_partition_blocks
