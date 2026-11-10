# nticLBM MPI partition module

The initialization partition module is in `src/mpi/`, parallel to `geometry`
and `input`, in namespace `ntic::lbm::mpi`. It uses bundled METIS for weighted
graph partitioning. Geometry generation, MPI message transport, LBM evolution,
and AMR migration are separate responsibilities. This implementation prepares
ownership and communication metadata; it does not call `MPI_Init` or move
field data between processes yet.

## Dependencies and build layout

```text
cmake/
    tests/mpi_tests.cmake      current test targets and their shared implementation
    3rdparties/
        thirdparties.cmake    dependency entry and source/license packaging
        gklib.cmake           GKlib integration
        metis.cmake           METIS integration
        standalone.cmake.in  source-bundle rebuild entry template
3rdParties/
    METIS/
        src/                  unmodified upstream METIS 5.2.1 source
        SOURCE_VERSION.txt
    GKlib/
        src/, include/, ...   unmodified upstream GKlib source
        SOURCE_VERSION.txt
    lib/                      generated libmetis.so* and libGKlib.so*
    licenses/
    THIRD_PARTY_NOTICES.md
build/
    3rdParties/                dependency build files and generated METIS header
    mpi/                      partition test executables and support archive
```

The root CMake project enables C and C++, builds both shared libraries from
the bundled sources, then links the MPI test implementation against METIS. Internal
build files remain under `build/3rdParties/`; both dynamic libraries are written
to the shared `3rdParties/lib/` directory. Shared libraries are ignored by Git;
generated headers stay under the untracked build tree. No build-time network download is needed.
METIS uses consistent 64-bit `idx_t` and `real_t` definitions.

Commit rule: track only `3rdParties/lib/.gitkeep` to preserve the directory.
All generated files and symlinks in that directory, including `.so` version
chains, remain untracked. Do not force-add build products.

Project-authored build configuration is grouped by purpose under `cmake/`.
Current test cases use `cmake/tests/`; dependency integration uses
`cmake/3rdparties/`. There is no standalone MPI module build layer or CMake file
under `src/`. Future application/release targets can have their own purpose-based
configuration directories when needed. Upstream CMake files inside the vendored
sources are retained unmodified.

The portable build entry is the root `CMakeLists.txt`:

```bash
cmake -S . -B build
cmake --build build --parallel 16
```

CMake expresses the GKlib -> METIS -> test-target dependency order and reuses
unchanged libraries incrementally. No machine-specific shell wrapper is
required. Compiler/toolchain selection, environment modules and job count are
chosen by the user for the target machine, outside the project configuration.
MPI tests link through `mpi_partition_test_support`, defined in their test
configuration, to `METIS::METIS`. On Linux this resolves to
`3rdParties/lib/libmetis.so`, whose dependency is `3rdParties/lib/libGKlib.so`.

Preserve upstream copyright and license files. GKlib includes LGPL and BSD
components in addition to Apache-2.0. See
[third-party notices](3rdParties/THIRD_PARTY_NOTICES.md) and the per-file
upstream license inventory. Shared-library install rules also include source,
licenses and rebuild recipes; binary distributions must retain those materials.

## Partition model

Input is a finalized `geometry::MultilevelGrid<2/3, Type>`. A uniform grid is
the level-0-only case. Existing geometry CellIndex arrays are never reordered.

Atomic graph vertices are legal axis-aligned topology boxes. With
`blockLevel=0` a vertex contains one root cell and all its descendants. A larger
`blockLevel` splits covered topology nodes down to that level, while keeping
coarser real leaves indivisible. No plane crosses a real cell or separates a
coarse source from its auxiliary family. All boxes use one common integer
coordinate system at `blockLevel`, even when their physical sizes differ.

Empty topology boxes complete the spatial tiling but carry zero load and are
excluded from METIS. They are assigned afterwards by a deterministic spatial
sweep. This does not allocate ordinary Dry cell data or change METIS ownership
of nonempty boxes. Same-owner boxes are merged only across complete shared
faces with equal orthogonal extents. The result covers the snapped geometry
domain exactly, without overlaps. Each rank may own several boxes, including
disconnected ones; cuts do not have to cross the entire domain.

Same-level lattice-neighbour pairs and coarse/fine interface pairs form graph
edges. Same-level edges are a conservative cost estimate: a physical boundary
Q that interrupts an otherwise wet-to-wet neighbour link is not yet subtracted
from the graph. Coarse/fine dependencies come from `interfaces()`, rather than
reconstructing them from proximity alone.

## Work and hardware parameters

All defaults are initial estimates, not benchmarked hardware constants.

| Parameter | Default | Interpretation |
| --- | --- | --- |
| `interiorWeight`, `boundaryWeight` | 1, 1 | Cost of one real-cell update |
| `levelWeights` | empty | Unit per-update cost multiplier at every level; explicit entries tune cost while preserving the independent `2^level` step frequency |
| `auxiliaryWeight` | 0.25 | Cost per auxiliary slot per source-coarse update |
| `interfaceWeight` | 1 | Cost per coarse/fine pair per coarse update, shared equally between the two participating boxes |
| `streamingEdgeWeight`, `interfaceEdgeWeight` | 1, 1 | Relative METIS edge costs |
| `streamingBytes`, `interfaceBytes` | 8, 64 | Estimated payload per neighbour/pair event |
| `blockLevel` | 0 | Atomic graph granularity; increase when root boxes are too heavy |
| `imbalanceTolerance` | 1.03 | METIS requested load tolerance, not a guaranteed achieved bound |
| `candidates`, `seed` | 4, 0 | Deterministic candidate runs |
| `loadSlack` | 0.03 | Maximum relative compute-time increase over the best generated candidate before considering communication/block preference |
| `blockPenalty` | 0 | Secondary time penalty per maximum per-rank merged block count; fewer total blocks break equal secondary scores |
| `graphThreads`, `mergePasses` | 8, 8 | Parallel graph-construction cap and rectangle-merge pass limit |

The hardware profile contains per-rank `computeRates` and `nodeIds`, and
intra/inter-node latency and bandwidth. Empty rates mean equal capacity.
Empty node IDs use one topology-neutral link class. Default latencies are
zero and bandwidths are `1e9`; callers should replace these with consistent
measured units. Compute costs divided by compute rates and communication
latency/bandwidth costs must be expressed in compatible time units.

Real-cell work is multiplied by `2^level`, not again by the number of children.
Covered parents and ordinary Dry leaves contribute no computational work.
Auxiliary slots do not advance independently; their cost and interface costs
are charged at source-coarse frequency under this initial tunable model.

METIS receives positive integer-scaled vertex/edge weights and capacity-based
target fractions. Tiny positive weights may lose relative precision during
quantization; candidate evaluation uses original floating-point costs. Graph
construction uses OpenMP integer counts; METIS calls are serialized. The
initialization graph is centralized, not a distributed ParMETIS graph.

Candidates are first filtered by hardware-normalized maximum compute time.
Within the configured `loadSlack`, selection considers maximum communication
time plus block penalty. Traffic is aggregated by rank pair and level/phase;
message-startup cost is not charged separately to every cut edge. The estimate
assumes sequential communication cost and does not model overlap or network
contention. No global-optimality or strict block-count guarantee is claimed.
An indivisible heavy box can prevent the requested balance; inspect
`largestAtomicWork` and request smaller atomic boxes by increasing `blockLevel`
when topology permits further splitting.

## Basic use and indices

```cpp
#include "mpi/partition.hpp"

namespace part = ntic::lbm::mpi;
part::PartitionOptions options;
options.ranks = 4;
options.blockLevel = 1;
options.hardware.computeRates = {1.0, 1.0, 2.0, 2.0};
options.hardware.nodeIds = {0, 0, 1, 1};
auto plan = part::partition(grid, options);
auto local = part::localPartition(grid, plan, rank);
```

`plan.owner[globalCell]` gives the owning rank. `plan.cells` groups global
CellIndex values by rank using half-open `plan.offsets` ranges. The grouped
list does not change geometry storage. `plan.blocks` contains merged spatial
boxes; a bound `u` becomes `plan.domainMin[d] + u * plan.blockSpacing`.
`statistics` reports loads, per-level work, interface work, normalized compute
time, estimated communication time, cut cost and block counts.

Local indices address `local.globalCells`: owned entries occupy
`[0, ownedCount)`, followed by read-only ghosts. Both groups preserve geometry
CellIndex order, hence level/role/type order. `local.ranges` describes the
actual contiguous subgroups; auxiliary families are not split by physical
type. `localIndex(global)` returns an indirect local index or
`invalidCellIndex`, using a sorted local lookup rather than a global-sized
reverse array on every rank.

Ghosts conservatively include all touching real neighbours and complete
participating coarse/fine blocks. Peer send/receive lists use local indices
and match in global-cell order. They describe potential dependencies, not an
implemented time-stepping communication schedule. Local interfaces use local
indices for coarse, auxiliary and fine cells.

Local boundary data rebuild Face offsets, retain repeated Face membership,
remap every retained record to a local index, and preserve Q rows exactly.
Owned and ghost boundary records are both included. Geometry/partition plans
are snapshots: after `refine`, `derefine` or `finalize`, rebuild the plan even
if the total cell count did not change. Global CellIndex is not an AMR-stable
identity. Field migration and incremental AMR repartitioning are not yet
implemented; the ownership/local-index interfaces provide their future input.

## Automatic tests

Each `.cpp` below is one executable; shared checks are in
`tests/mpi_tests/test_partition_checks.hpp`.

| Executable under `build/mpi/` | Coverage |
| --- | --- |
| `test_partition_ownership` | Uniform and weighted multilevel ownership, auxiliary binding and unchanged geometry storage |
| `test_partition_connectivity` | 2D/3D ghosts, reciprocal peer order, local interfaces, Face/Q remapping and single-rank case |
| `test_partition_validation` | Thread-count determinism, invalid input, empty ranks and fully Dry domain |
| `test_partition_blocks` | Exact rectangular coverage, no cell cuts, default `2^level` frequencies, hardware capacities and merging |

All four targets are included by `cmake/tests/mpi_tests.cmake`; the project-root
`run_mpi_tests.sh` invokes them without interactive input. It follows existing
runner style and does not aggregate exit codes: inspect each program's result.
These are partition-metadata tests and do not require `mpirun`.

In the current development workflow, local work is source/configuration editing
only, and builds/tests run remotely. The following commands are an example
specific to the current Donau test server, not a universal build requirement:

```bash
module load metax/metax-2.33_KylinSec-3.5.3-arm64
cmake_hpcc -S . -B build
make_hpcc -C build -j 16
bash run_mpi_tests.sh
```

Remote validation on 2026-11-10 used
`metax/metax-2.33_KylinSec-3.5.3-arm64`: the complete project build succeeded,
both shared libraries were generated in `3rdParties/lib/`, and `ldd` confirmed
that the partition executable loads METIS and GKlib from that directory.
All four partition test executables passed their ten checks. This validates
partition metadata, not actual MPI transport or calibrated hardware timings.
