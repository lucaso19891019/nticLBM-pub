#!/bin/bash

MPI_PARTITION_APP_DIR=./build/mpi/partition


#=============================================================================
# 1. Uniform and weighted multilevel ownership
#=============================================================================

$MPI_PARTITION_APP_DIR/test_partition_ownership


#=============================================================================
# 2. Ghosts, peer lists, interface indices and boundary Q
#=============================================================================

$MPI_PARTITION_APP_DIR/test_partition_connectivity


#=============================================================================
# 3. Determinism, invalid input and infeasible/empty domains
#=============================================================================

$MPI_PARTITION_APP_DIR/test_partition_validation


#=============================================================================
# 4. Rectangular coverage, level frequencies and hardware capacities
#=============================================================================

$MPI_PARTITION_APP_DIR/test_partition_blocks


#=============================================================================
# 5. Boundary Face membership costs and multi-Face cells
#=============================================================================

$MPI_PARTITION_APP_DIR/test_partition_boundary_costs


#=============================================================================
# 6. Single-box internal dependencies and node communication costs
#=============================================================================

$MPI_PARTITION_APP_DIR/test_partition_intra_rank
