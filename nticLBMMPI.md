# nticLBM MPI partition module

The MPI module is in `src/mpi/`, parallel to `geometry` and `input`. Its
initialization partition functionality is in `src/mpi/partition/`, in namespace
`ntic::lbm::mpi`. It uses bundled METIS for weighted
graph partitioning. Geometry generation, MPI message transport, LBM evolution,
and AMR migration are separate responsibilities. This implementation prepares
ownership and communication metadata; it does not call `MPI_Init` or move
field data between processes yet.

## Dependencies and build layout

```text
src/
    mpi/
        partition/            partition sources; future MPI functions use sibling directories
tests/
    mpi_tests/
        partition_tests/      six automatic tests, two interactive cases and shared checks
        stls/                 MPI-owned STL fixtures
        outputs/              MPI-local VTK helpers and ignored per-case *_vtks results
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
    mpi/
        partition/            partition test executables and support archive
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
The interactive branched-channel target links its own
`mpi_partition_stl_test_support`, compiled from production STL reader,
validator and geometry sources. It does not link another module's test target.
MPI-owned STL fixtures are in `tests/mpi_tests/stls/`. The complete
`vtk_output_2d.hpp` and `vtk_output_3d.hpp` in
`tests/mpi_tests/outputs/` are maintained independently within the MPI tests.
The Taylor-Couette case includes the 2D header; the branched-channel case
includes the 3D header, which also includes its local 2D header. They provide
the existing multilevel-cell, boundary-Q and partition-block writer APIs;
the cell writer retains its additional-cell-data callback for partition
fields. Test headers and fixtures are not loaded from other test modules. `MPI_PARTITION_TEST_OUTPUT_DIR` names the build executable directory;
the `MPI_TEST_OUTPUT_DIR` compile definition names the separate generated-results
directory under `tests/mpi_tests/outputs/`. `MPI_TEST_SOURCE_DIR` points to
`tests/mpi_tests/` for loading MPI-owned fixtures. Source and generated libraries
used by the build remain inside the repository; deployment transfer directories
are not build or runtime dependencies.

Preserve upstream copyright and license files. GKlib includes LGPL and BSD
components in addition to Apache-2.0. See
[third-party notices](3rdParties/THIRD_PARTY_NOTICES.md) and the per-file
upstream license inventory. Shared-library install rules also include source,
licenses and rebuild recipes; binary distributions must retain those materials.

A recorded full vendor build emitted 244 upstream warnings: 168 unused,
60 signedness comparisons, 3 empty bodies, 9 potentially clobbered variables
and 4 absolute-value type mismatches. They have not been suppressed, and the
vendored sources remain unmodified. The clobbered group includes renumbering
state around `setjmp`/`longjmp`; this module uses zero-based numbering, but
exception paths are not fully validated. In `libmetis/parmetis.c`, the four
`abs(int)` calls on 64-bit index expressions are a potential narrowing defect
in the P-suffixed node-refinement routines, outside the current
`METIS_PartGraphKway` path. A successful build does not establish those other
upstream paths as safe; using them would require separate review and testing.
The final source-only rebuild emitted no warnings and reused the existing
dependency libraries; it did not eliminate or suppress those upstream warnings.

## Partition model

Input is a finalized `geometry::MultilevelGrid<2/3, Type>`. A uniform grid is
the level-0-only case. Existing geometry CellIndex arrays are never reordered.

Atomic graph vertices are legal axis-aligned topology boxes. With
`blockLevel=0` a vertex contains one root cell and all its descendants. A larger
`blockLevel` splits covered topology nodes down to that level, while keeping
coarser real leaves indivisible. No plane crosses a real cell or separates a
coarse source from its auxiliary family. All boxes use one common integer
coordinate system at `blockLevel`, even when their physical sizes differ.

Only positive-work atoms participate in METIS and ownership. Zero-work atoms
retain owner `-1` and an invalid block index. Each rank receives exactly one
axis-aligned rectangle/box containing all its assigned wet atoms. Dry gaps
may be padding inside that box, without allocating or charging Dry cell data.
Boxes do not overlap, every real wet/auxiliary slot is wholly contained by its
unique owner box, and full-domain Dry coverage is not required. Disconnected
wet pieces may therefore belong to the same rank inside one bounding box.
Coarse/auxiliary families and topology atoms remain indivisible.

METIS assignments guide constrained aggregation, rather than being final rank
labels: weighted label centroids order rank groups along candidate axes, and
label migration helps choose legal cuts after balancing group work against
the actual ranks' compute capacities. A separate balanced geometric seed
keeps regular central cuts visible. Every candidate is a slicing tree with
one rank per leaf; cuts may end at an ancestor cut instead of spanning the
whole domain. Leaf boxes are trimmed to their wet-atom bounding extents.

The bounded global search retains several complete feasible layouts and
rebuilds multi-rank subtrees, including the root. It varies axes, cut positions
and rank groups. Joint moves of two cut planes include ancestor/descendant
pairs: a middle rank can gain work at one boundary and lose it at another
without requiring either intermediate one-plane state to improve. Each
complete proposal is geometrically routed and its rank-pair traffic, phases,
node classes and compute times are recomputed from the full graph. The best
feasible compute time found is retained throughout the search; within its
`loadSlack` budget, maximum communication time is minimized, then compute time
and deterministic ownership break ties. A low-communication winner cannot
move the compute budget away from that historical anchor.

The search is a heuristic, not a global optimum or an exhaustive search of
all slicing trees. Non-slicing rectangular layouts are outside this search
space. `blockLevel=0` also binds all descendants of a root as an additional
granularity restriction, beyond the mandatory real coarse/auxiliary family.
Increase it when smaller topology atoms are useful. Fewer populated atoms
than ranks is rejected as infeasible; an empty domain returns an empty plan.
Failure of the bounded slicing search has a separate error and does not prove
that no general rectangular layout exists. Public `partitionGraph` callers
must supply mutually nonoverlapping populated atoms; production grid atoms
already satisfy this. Edges must join positive-work atoms with positive finite
weights, nonnegative finite bytes and supported phase bits.

Same-level lattice-neighbour pairs and coarse/fine interface pairs form graph
edges. Same-level edges are a conservative cost estimate: a physical boundary
Q that interrupts an otherwise wet-to-wet neighbour link is not yet subtracted
from the graph. Coarse/fine dependencies come from `interfaces()`, rather than
reconstructing them from proximity alone.

## Work and hardware parameters

All defaults are initial estimates, not benchmarked hardware constants.

| Parameter | Default | Interpretation |
| --- | --- | --- |
| `wetWeight` | 1 | Ordinary update cost for every real wet cell, including Boundary |
| `boundaryExtraWeight` | 1 | Uniform extra cost per stored Face/cell boundary record; zero is allowed |
| `levelWeights` | empty | Unit per-update cost multiplier at every level; explicit entries tune cost while preserving the independent `2^level` step frequency |
| `auxiliaryWeight` | 0.25 | Cost per auxiliary slot per source-coarse update |
| `interfaceWeight` | 1 | Cost per coarse/fine pair per coarse update, shared equally between the two participating boxes |
| `streamingEdgeWeight`, `interfaceEdgeWeight` | 1, 1 | Relative METIS edge costs |
| `streamingBytes`, `interfaceBytes` | 8, 64 | Estimated payload per neighbour/pair event |
| `blockLevel` | 0 | Atomic graph granularity; increase when root boxes are too heavy |
| `imbalanceTolerance` | 1.03 | METIS requested load tolerance, not a guaranteed achieved bound |
| `candidates`, `seed` | 4, 0 | Deterministic candidate runs |
| `loadSlack` | 0.03 | Relative compute budget above the best feasible compute time found during the entire search; communication is selected within that budget |
| `graphThreads` | 8 | Parallel graph-construction thread cap |
| `refinementRounds` | 8 | Global refinement round limit; zero evaluates initialization only |
| `searchWidth` | 4 | Retained feasible layouts (2 to 16), including the best compute anchor |
| `cutCandidates` | 4 | Legal split alternatives and bounded per-round subtree/pair search budget (1 to 32) |

The hardware profile contains per-rank `computeRates` and `nodeIds`, and
intra/inter-node latency and bandwidth. Empty rates mean equal capacity.
`intraRankLatency` (default 0) and `intraRankBandwidth` (default `1e9`)
remain available for local-block exchange metadata. With the single-box
constraint there is no separate same-rank block interface and its cost is zero;
ordinary cell-neighbour dependencies still exist within the box.
Empty node IDs use one topology-neutral link class. Default latencies are
zero and bandwidths are `1e9`; callers should replace these with consistent
measured units. Compute costs divided by compute rates and communication
latency/bandwidth costs must be expressed in compatible time units.

Real-cell work is multiplied by `2^level`, not again by the number of children.
Interior costs `wetWeight`. A real Boundary cell with `m` stored Face/cell
membership records costs `wetWeight + max(1, m) * boundaryExtraWeight`.
The ordinary wet update is charged once per unique cell; each stored Face
membership adds one boundary access cost, even when that record's Q links are
all invalid. Costs are not multiplied by the number of valid Q links. A real
Boundary cell with no membership records retains one extra cost as a fallback.
The total is multiplied by the cell's level cost and `2^level` update frequency.
`wetWeight` must be positive and finite; `boundaryExtraWeight` must be finite
and nonnegative. Auxiliary slots retain their separate cost model.

The uniform extra weight is a partition-stage estimate, not an implementation
of boundary conditions or their future dependency/scheduling rules. Boundary
condition types are not configured by this partition API. Any future nonlocal
boundary dependencies, including periodic pairing, must be supplied to the
connectivity and scheduling design; they cannot be inferred from this weight.
Covered parents and ordinary Dry leaves contribute no computational work.
Auxiliary slots do not advance independently; their cost and interface costs
are charged at source-coarse frequency under this initial tunable model.

METIS receives positive integer-scaled vertex/edge weights and capacity-based
target fractions. Tiny positive weights may lose relative precision during
quantization; candidate evaluation uses original floating-point costs. Graph
construction uses OpenMP integer counts; METIS calls are serialized. The
initialization graph is centralized, not a distributed ParMETIS graph.

Every candidate is evaluated with actual hardware-normalized compute and
communication costs. The load budget is relative to the best feasible layout
seen across initialization and refinement. Traffic is aggregated by rank pair
and level/phase, so startup cost is not charged to every individual edge.
Communication overlap and network contention are not modeled.

Each rank's internal cell dependencies stay inside its single box and use
the existing local indirect indices, without a separate block-exchange list
or ghost allocation. Cross-rank traffic is aggregated by rank pair and uses
same-node or cross-node parameters from `nodeIds`. The cost model does not
assume communication hiding or calibrated hardware throughput.
An indivisible heavy box can prevent the requested balance; inspect
`largestAtomicWork` and request smaller atomic boxes by increasing `blockLevel`
when topology permits further splitting.

## Basic use and indices

```cpp
#include "mpi/partition/partition.hpp"

namespace part = ntic::lbm::mpi;
part::PartitionOptions options;
options.ranks = 4;
options.blockLevel = 1;
options.hardware.computeRates = {1.0, 1.0, 2.0, 2.0};
options.hardware.nodeIds = {0, 0, 1, 1};
auto plan = part::partition(grid, options);
auto local = part::localPartition(grid, plan, rank);
```

For example, tune the uniform wet/boundary costs and local exchange model:

```cpp
options.wetWeight = 1.0;
options.boundaryExtraWeight = 0.5; // Boundary cost is 1.0 + 0.5 per Face membership.
options.hardware.intraRankLatency = 0.0;
options.hardware.intraRankBandwidth = 1.0e9;
```

`plan.owner[globalCell]` gives the owning rank. `plan.cells` groups global
CellIndex values by rank using half-open `plan.offsets` ranges. The grouped
list does not change geometry storage. `plan.blocks` contains one spatial
box per rank; a bound `u` becomes `plan.domainMin[d] + u * plan.blockSpacing`.
`statistics` reports loads, per-level work, interface work, normalized compute
time, estimated communication time, cut cost and block counts.
`plan.cellBlocks[globalCell]` references the rank box, and
`plan.blockExchanges` preserves streaming/coarse-fine graph traffic between
the boxes of different ranks. `partitionGraph` alone returns
`cellBlocks` indexed by input graph vertex; `partition(grid)` expands it to
global cells. Communication statistics include separate
`intraRankCommunicationTimes` and `interRankCommunicationTimes`; their sum is
`communicationTimes` and is used in candidate selection.

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
`local.cellBlocks` maps each owned/ghost local index to its final block.
`local.intraRankExchanges` is empty under the single-box-per-rank contract.
Internal wet-cell dependencies remain reachable through the owned local
indices, including disconnected wet pieces within the same Dry-padded box.
Ghosts and peer lists are only needed for cross-rank dependencies; no separate
same-rank block copy or message is inferred from wet-region connectedness.

Local boundary data rebuild Face offsets, retain repeated Face membership,
remap every retained record to a local index, and preserve Q rows exactly.
Owned and ghost boundary records are both included. Geometry/partition plans
are snapshots: after `refine`, `derefine` or `finalize`, rebuild the plan even
if the total cell count did not change. Global CellIndex is not an AMR-stable
identity. Field migration and incremental AMR repartitioning are not yet
implemented; the ownership/local-index interfaces provide their future input.

## Automatic tests

Each `.cpp` below is one executable; shared checks are in
`tests/mpi_tests/partition_tests/test_partition_checks.hpp`.

| Executable under `build/mpi/partition/` | Coverage |
| --- | --- |
| `test_partition_ownership` | Uniform and weighted multilevel ownership, auxiliary binding and unchanged geometry storage |
| `test_partition_connectivity` | 2D/3D ghosts, reciprocal peer order, local interfaces, Face/Q remapping and single-rank case |
| `test_partition_validation` | Thread-count determinism, invalid input, infeasible rank counts and fully Dry domain |
| `test_partition_blocks` | Unique wet-cell box coverage, default `2^level` frequencies, hardware capacities, analytic ring communication and multi-rank refinement |
| `test_partition_boundary_costs` | Analytic box load, per-membership boundary costs, multi-Face extra work without duplicate wet updates and zero surcharge |
| `test_partition_intra_rank` | Single-box L-shaped/disconnected wet regions, accessible internal dependencies, no spurious local exchange and node classes |

All six automatic targets are included by `cmake/tests/mpi_tests.cmake`; the project-root
`run_mpi_tests.sh` invokes them without interactive input. It follows existing
runner style and does not aggregate exit codes: inspect each program's result.
These are partition-metadata tests and do not require `mpirun`.
The root runner and CMake entry retain their MPI names so other MPI test groups
can be added later; their current partition executables are under
`build/mpi/partition/`. The two interactive applications below are not invoked
by the automatic runner.

## Interactive partition cases

These applications generate geometry, multilevel grids and partition metadata;
they do not advance LBM or launch MPI ranks. The annulus reads finest spacing
`hmin`, then the narrow wall-band thickness in finest-cell units. The STL case
derives `hmin` from its smallest port and reads only the wall-band cell count.
Both fix `blockLevel=0` and equal rank compute capacities.
`boundaryExtraWeight` retains its default value of 1 per stored Face record.
A bounded global partition-refinement search is used; communication hiding is not implemented. Reported hardware costs are uncalibrated model estimates, not
measured execution times.

| Executable under `build/mpi/partition/` | Recommended inputs | Preset |
| --- | --- | --- |
| `test_partition_taylor_couette` | `0.01`, then `4` | 2D concentric annulus, 4 equal ranks, node placement unspecified |
| `test_partition_branched_channel` | `2` only | Port-derived finest spacing, four topology levels, 8 equal ranks, node IDs `{0,0,0,0,1,1,1,1}` |

Run from the project root and enter the case-specific values when prompted:

```bash
OMP_NUM_THREADS=32 ./build/mpi/partition/test_partition_taylor_couette
OMP_NUM_THREADS=128 ./build/mpi/partition/test_partition_branched_channel
```

The Taylor-Couette case sets `Rin=900*hmin` and `Rout=1000*hmin`. Thus the
inner/outer radius and diameter ratios are 9:10, and the radial gap is exactly
100 finest cell widths. With the recommended inputs the radii are 9 and 10,
the three level-0/1/2 spacings are 0.04/0.02/0.01 (`maxLevel=2`), and the
level-1/level-2 wall bands are 0.08/0.04.
Both circular walls seed the bands; there are no extra rectangular refinement
zones. Scaling `hmin` scales the whole case, so the level-0 bounding topology
remains approximately 500 by 500 cells. Large band-cell inputs can fill the
gap with refinement; the program warns when the nominal wider bands meet.

The branched-channel case reads the existing
`tests/mpi_tests/stls/smooth_irregular_branched_channel.stl`, runs
full topology/normal validation, and constructs an internal-flow STL geometry.
It uses the seven formally identified planar port Faces and their actual
vertices. Each port is projected onto an orthonormal basis in its plane; the
minimum Feret width is the smallest separation between parallel supporting
lines of its convex hull. The smallest of the seven widths is `Dmin`, and
`hmin=Dmin/25` is fixed: the user cannot override spacing or the count of 25.
Area-equivalent diameters are printed for comparison, but do not determine
resolution. Bounding-box diagonals, external guide-cylinder radii and their
coverage inflation factors are not used. This metric describes the overall
port width; for a general nonconvex port it is not a local neck-width measure.

The case fixes `maxLevel=3` and `h0=8*hmin` for four topology levels. Runtime
output lists all seven measured widths, the selected Face and derived spacings;
the selected Face and numerical diameter are not hardcoded. The validated
bundled fixture gives `Dmin` approximately 4.8589440326 and `hmin` approximately
0.1943577613, in model units. With the recommended single input
`2`, the narrow band is `2*hmin`; the other bands are two and four times wider:

| Level | Spacing | Wall-band thickness |
| --- | --- | --- |
| 0 | 1.5548620904 | Base topology, no wall-band rule |
| 1 | 0.7774310452 | 1.5548620904 |
| 2 | 0.3887155226 | 0.7774310452 |
| 3 | 0.1943577613 | 0.3887155226 |

Only Other Surface (Face 0) seeds refinement; planar port caps are excluded as
wall-band seeds. Side-wall proximity and 2:1 balancing may still refine cells
around cap rims, so entire caps are not promised to remain coarse. Larger band
inputs can refine entire narrow branches; four topology levels do not guarantee
that every level has real cells in each branch. The smaller physical spacing
substantially increases cost relative to the old `hmin=1` three-level case.
The validated run started from 696080 level-0 topology cells. Retained counts,
active-level work and peak memory for the recommended setting are recorded
below; changing the band input can change those results.

The annulus caps geometry-analysis OpenMP threads at 32; the STL case caps them
at 128, both respecting the available OpenMP maximum. Graph construction is
capped at 8 threads. The annulus limits topology to 3,000,000 nodes and retained
storage to 2,000,000 slots; the STL case uses 6,000,000 nodes and 4,000,000 slots.
Both report errors if limits are exceeded. These are allocation guards, not
unconditional preallocations or a measured memory budget. The validated STL
run completed within its limits. STL analysis and VTK writing can take significant time. Same-rank, same-node and cross-node
latency/bandwidth retain their initial
uncalibrated defaults. The 2D case leaves node IDs empty and writes `node=-1`;
the 3D output uses its explicit 0/1 node assignment. These are simulated
partition profiles, not evidence of actual MPI execution on two nodes.

Each run recreates only its dedicated directory under
`tests/mpi_tests/outputs/`:

```text
test_partition_taylor_couette_vtks/
test_partition_branched_channel_vtks/
```

The directory names are checked before deletion, and links/aliases to another
case or location are rejected. Generated `*_vtks/` directories are ignored by
Git; the reusable output header remains source code. Each directory contains:

| File | Contents |
| --- | --- |
| `partition.vtk` | Real wet cells as VTK QUAD/HEXAHEDRON cells, with scalar, level, role, CellIndex, SourceNode, rank, node and BlockIndex |
| `auxiliary.vtk` | Auxiliary slots with the same fields; kept separate from real-cell output |
| `blocks.vtk` | One rectangle/box per rank with rank, node, BlockIndex and work; may include Dry padding, without requiring full Dry-domain coverage |
| `boundary_q.vtk` | Q intersection points only, retaining Face/link/CellIndex and level provenance |
| `partition_stats.txt` | Case settings and per-rank wet/auxiliary counts, per-level slot work, separate interface work, compute/local/remote cost estimates, block counts and overall imbalance |

All four VTK files use legacy VTK 3.0 `BINARY` with big-endian payloads.
Coordinates, work and Q use 64-bit doubles; connectivity/cell types and signed
fields use 32-bit integers. `CellIndex`, `SourceNode` and `BlockIndex` retain
their field names and use standard 32-bit `unsigned_int`; overflow is rejected
explicitly, as are connectivity/counts exceeding signed 32-bit limits.
`boundary_q.vtk` contains vertex cells, with no Q link segments.
`partition_stats.txt` remains ordinary text. MPI owns its output helpers;
geometry test output helpers are maintained independently.

The statistics are also printed to the console. Per-level slot work excludes
coarse/fine pair work, which is reported separately; auxiliary costs remain
charged at their source-coarse update frequency. The automatic runner still
contains only its six noninteractive applications.

Boundary records carry explicit window-local `boundaryCells` scalar indices.
Multilevel assembly maps these to grid CellIndex values before partition-local
Q mapping. Indices are not reconstructed from coordinates; retained coordinates
must exactly match the cell center. Tile centers use the shared `gridFrame`
origin and integer offsets, including through primitive, STL and Boolean
analysis. Standalone analysis retains its existing coordinate arithmetic.

Remote validation on 2026-11-13 completed the full build, all twelve checks
in the six automatic programs, and both manual applications (exit status 0).
The annulus produced four quadrant boxes with exactly equal modeled load
(`maxToMean=1`). The four-level STL case produced eight boxes, one per rank,
with `maxToMean=1.035902809`: the largest normalized load is 3.59% above the
mean. Its selected maximum compute cost, 5254988.5, is about 2.05% above the
best feasible anchor, 5149423.5, within the configured 3% load budget for
communication selection. It evaluated 645 candidates and 512 joint-cut
proposals over eight rounds; total work remains 40582869 and the geometry
cell/interface counts are unchanged. The annulus took 59.96 seconds; the STL
case took 21 minutes 55.95 seconds with 1255964 KiB maximum RSS. These are
observed initialization times, not calibrated hardware performance constants.

The manual summary reports evaluated candidate count, completed rounds,
joint-cut proposals, the retained best compute anchor and the termination
reason. Search storage is bounded by the small `searchWidth` population plus
the atom graph; per-proposal geometric routing and full graph evaluation
scale with atom/edge count. Cut enumeration sorts active atoms along each
axis, scans legal gaps and bounds retained alternatives. Recursive tree
construction also has a shared attempt budget; cut-plane search uses only populated atoms; there is no Dry-owner completion
or forced full-domain cover. Complete owner/map arrays still include unused
zero-work atoms, and routing/statistics loops skip them. No unbounded all-plane
combination search is used. No optimality guarantee
follows from reaching a round limit or a stable bounded proposal cycle.

In the current development workflow, local work is source/configuration editing
only, and builds/tests run remotely. The following commands are an example
specific to the current Donau test server, not a universal build requirement:

```bash
module load metax/metax-2.33_KylinSec-3.5.3-arm64
cmake_hpcc -S . -B build
make_hpcc -C build -j 16
bash run_mpi_tests.sh
```

The analytic boundary-cost test has 216 wet cells, 152 distinct boundary
cells and 216 Face records: ordinary wet work is charged once per cell and
the boundary extra cost once per stored Face record. Intra-rank dependencies
and local/same-node/cross-node exchange costs are covered, including
dependencies within one box containing disconnected wet pieces, with no extra local block exchange. These tests validate
partition metadata, not actual MPI transport or calibrated hardware timings.
