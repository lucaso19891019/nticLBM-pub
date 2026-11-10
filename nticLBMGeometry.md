# nticLBM Geometry

The geometry module constructs primitive, STL, and Boolean geometries and
analyzes them on a uniform Cartesian grid. Its outputs are a cell-state scalar
map, boundary-cell coordinates, boundary Face groups, and link intersection
fractions (`Q`).

All geometry types use the same analysis interface:

```cpp
geometry.interiorAreaAnalysis(result);
geometry.boundaryQAnalysis(result);

// Convenience call: perform both steps in order.
geometry.analysis(result);
```

Geometry construction and analysis do not assign velocity, pressure, wall,
inlet, or outlet boundary conditions. Face IDs identify geometric regions for
subsequent use; they are not boundary-condition types.

## Headers and namespaces

The main headers are:

| Header relative to `src/geometry` | Purpose |
| --- | --- |
| `primitive/rectangle.hpp` | Axis-aligned 2D rectangle |
| `primitive/circle.hpp` | 2D circle |
| `primitive/box.hpp` | Axis-aligned 3D box |
| `primitive/sphere.hpp` | 3D sphere |
| `primitive/cylinder.hpp` | Finite cylinder aligned with X, Y, or Z |
| `stl/stl_geometry.hpp` | Validated STL geometry and component interpretation |
| `boolean_geometry.hpp` | Boolean composition and nested analysis |
| `geometry_analysis.hpp` | Shared analysis inputs and outputs |
| `bounding_box.hpp`, `point.hpp`, `flow_type.hpp`, `axis.hpp` | Common geometry types |

Geometry types belong to `ntic::lbm::geometry`. Lattice types belong to
`ntic::lbm::lattice`; STL reading and validation belong to `ntic::lbm::stl`.
Examples below assume these aliases:

```cpp
namespace geo = ntic::lbm::geometry;
namespace lat = ntic::lbm::lattice;
```

Primitive implementations are in headers. STL analysis additionally requires
the STL geometry and input implementation files to be linked. Existing targets
in `cmake/tests/geometry_stl_tests.cmake` and
`cmake/tests/geometry_boolean_tests.cmake` show the required source files and
include directories. The current project uses C++17 and OpenMP.

## Flow type and analysis domain

`FlowType::Internal` selects the fluid region inside a primitive. Its default
analysis domain is the geometry's bounding box, `bounds`.

`FlowType::External` selects the fluid region outside a primitive, restricted
to a finite `openBox`. The open box must have positive dimensions and contain
the geometry. Pass it explicitly to avoid interactive input:

```cpp
const geo::BoundingBox openBox{{-6.0, -6.0, -6.0}, {6.0, 6.0, 6.0}};
const geo::Sphere obstacle(
    {0.0, 0.0, 0.0}, 2.0, geo::FlowType::External, &openBox);
```

The constructor copies the supplied open box. If an external geometry is
constructed without one, the constructor prompts for it. `openBox` defines
the analysis extent; its faces are not automatically added as wall, inlet,
or outlet Face groups.

For STL geometries containing several closed components, the flow region is
determined by component nesting as described below, rather than by treating
each component as an independent fluid volume.

All coordinates, dimensions, and grid spacing must use the same length unit.
Geometry constructors do not perform the configuration module's unit parsing.

## Primitive construction

Primitive constructors default to `FlowType::Internal`. Their common trailing
arguments are `FlowType` and an optional pointer to an open box.

```cpp
#include "primitive/rectangle.hpp"
#include "primitive/circle.hpp"
#include "primitive/box.hpp"
#include "primitive/sphere.hpp"
#include "primitive/cylinder.hpp"

const geo::Rectangle rectangle({0.0, 0.0}, {10.0, 6.0});
const geo::Circle circle({5.0, 3.0}, 2.0);

const geo::Box box({0.0, 0.0, 0.0}, {10.0, 6.0, 8.0});
const geo::Sphere sphere({5.0, 3.0, 4.0}, 2.0);
const geo::Cylinder cylinder(
    {5.0, 3.0, 4.0}, 2.0, 8.0, geo::Axis::Z);
```

Rectangles and boxes take minimum and maximum corners. Circles and spheres
take a center and radius. Cylinders take a center, radius, full axial length,
and axis; the two caps lie at `center[axis] +/- length / 2`.

Dimensions must be positive. Arbitrary cylinder directions and rotated boxes
or rectangles are not currently provided by these primitive types.

Each primitive exposes `bounds`, `flowType`, and `openBox`, together with its
shape parameters. `boundingBox()` returns its bounds. Prefer constructing or
translating a geometry through its supported operations rather than changing
individual parameters without updating its bounds.

### Translation

Primitive types and `STLGeometry` support `translate()`:

```cpp
geo::Circle circle({4.0, 3.0}, 2.0);
const geo::Point2D targetMin{0.0, 0.0};
circle.translate(&targetMin);
```

For internal flow, the supplied point is the new **minimum bounding-box
corner**, not the new center or a displacement. A target is required.

For external flow, call `geometry.translate()` without a target. This moves
the open-box minimum to the origin and moves the geometry by the same amount.
An explicit target is not accepted in this mode.

Translate operands before building a Boolean expression. `BooleanGeometry`
does not currently expose a translation method. Recompute analysis results
after changing a geometry's position.

## Uniform-grid analysis

Choose a lattice model compatible with the geometry dimension:

| Analysis type | Supported lattices |
| --- | --- |
| `GeometryAnalysis2D<Type>` | `D2Q9` |
| `GeometryAnalysis3D<Type>` | `D3Q15`, `D3Q19`, `D3Q27` |

A complete 3D primitive example is:

```cpp
#include "primitive/box.hpp"
#include "geometry_analysis.hpp"

const geo::Box box({0.0, 0.0, 0.0}, {6.0, 4.0, 5.0});
geo::GeometryAnalysis3D<lat::LatticeType::D3Q19> result(0.1);
box.analysis(result);
```

The constructor argument is the uniform grid spacing `h`; use a positive,
finite value. The analysis computes dimensions using the ceiling of each
domain extent divided by `h`. The final grid edge can therefore extend beyond
the requested domain maximum by less than one spacing.

For 3D, cell `(i,j,k)` has:

```text
linear index = i + nx * (j + ny * k)
center[d]    = domain.min[d] + (index[d] + 0.5) * h
```

For 2D, use `i + nx * j` and omit Z.

### Analysis settings

| Field | Meaning |
| --- | --- |
| `gridSpacing` | Uniform spacing used by this analysis |
| `domainBounds` | Optional explicit input domain; otherwise use the geometry's default domain |
| `complemented` | Temporarily analyze the complementary fluid region; default `false` |
| `domain` | Actual domain selected by the analysis |

For example, an explicit observation domain can be supplied without modifying
the geometry:

```cpp
geo::GeometryAnalysis3D<lat::LatticeType::D3Q19> result(0.2);
result.domainBounds = geo::BoundingBox{{-1.0, -1.0, -1.0}, {7.0, 5.0, 6.0}};
box.analysis(result);
```

`domainBounds` and `complemented` are analysis settings, not permanent shape
properties. Boolean analysis uses them internally to align operands and
evaluate complements. Keep the geometry, lattice, spacing, domain, and
complement setting consistent between the separate interior and Q calls.
After any change, recompute the analysis.

### Scalar map

`scalar` stores one cell-state value per cell, using the numeric values of
`CellType`:

| Value | State | Meaning |
| --- | --- | --- |
| `0` | `Dry` | Outside the selected solving region |
| `1` | `Boundary` | A solving-region cell classified near an active geometric boundary |
| `2` | `Interior` | An interior solving-region cell |

Boundary cells belong to the solving region. Cells whose centers lie exactly
on a geometric surface are treated as Dry. A Dry cell is not promoted merely
because its center or extent touches a surface. The map is a discrete
classification, not a signed-distance field.

`interiorAreaAnalysis()` produces the scalar map and boundary coordinates/group
information needed by subsequent analysis. It does not produce the complete Q
array. Call `boundaryQAnalysis()` afterwards, or use `analysis()` to obtain the
complete result.

## Boundary Q and grouped storage

For a boundary cell centered at `x`, lattice direction `e_i`, and spacing `h`,
the Q intersection is:

```text
x_wall = x + q * h * e_i
```

Valid Q values satisfy `0 < q <= 1`. `-1` means that the entry has no valid
intersection owned by this Face for this link. The rest direction is excluded.
Q is a fraction along the lattice link, not an absolute distance: the physical
distance is `q * h * |e_i|`.

The shared result uses structure-of-arrays storage:

| Field | Meaning |
| --- | --- |
| `boundaryX`, `boundaryY`, `boundaryZ` | Boundary-entry cell centers; Z exists only in 3D |
| `nBoundaries` | Number of geometric Face/edge groups, including empty groups |
| `boundaryOffsets` | Group boundaries; size is `nBoundaries + 1` |
| `q` | Entry-major Q values; size is `boundaryX.size() * (nStencils - 1)` |

Group `face` occupies the half-open entry range:

```cpp
const auto begin = result.boundaryOffsets[face];
const auto end = result.boundaryOffsets[face + 1];
```

Access a Q value using the actual lattice direction index, starting at one:

```cpp
using Lattice = decltype(result)::Lattice;
constexpr std::size_t links = Lattice::nStencils - 1;

for(std::size_t face = 0; face < result.nBoundaries; ++face)
{
    for(std::size_t entry = result.boundaryOffsets[face];
        entry < result.boundaryOffsets[face + 1]; ++entry)
    {
        for(std::size_t direction = 1; direction < Lattice::nStencils; ++direction)
        {
            const double q = result.q[entry * links + direction - 1];
            if(q == -1.0) continue;
            // Use this Face's Q intersection for this cell and direction.
        }
    }
}
```

A physical boundary cell can occur in two or more groups. Consequently,
`boundaryX.size()` counts **entries**, not unique boundary cells. An entry can
also have all Q values set to `-1`. Use the scalar map to identify physical
boundary cells, and inspect valid Q values to determine link ownership.

Each physical cell/link has a single selected Q owner. Other Face entries for
the same cell retain `-1` for that link. Empty Face groups remain present so
that group numbering does not depend on grid resolution. Equal-distance hits
use deterministic Face ordering; approximate Q-distance ties are not used to
merge unrelated surfaces.

### Primitive Face numbering

IDs are zero-based and remain the same for internal and external analysis:

| Geometry | Face IDs |
| --- | --- |
| Rectangle | `0`: bottom/Y-min, `1`: right/X-max, `2`: top/Y-max, `3`: left/X-min |
| Circle | `0`: circular boundary |
| Box | `0`: X-min, `1`: X-max, `2`: Y-min, `3`: Y-max, `4`: Z-min, `5`: Z-max |
| Sphere | `0`: spherical boundary |
| Cylinder | `0`: lateral surface, `1`: negative-axis cap, `2`: positive-axis cap |

For cylinders, negative/positive refer to the selected X, Y, or Z axis.

## STL geometry

Read and fully validate the STL before constructing `STLGeometry`:

```cpp
#include "stl/stl_geometry.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"
#include <utility>

auto data = ntic::lbm::stl::read("channel.stl");
ntic::lbm::stl::FacetTopology topology;
ntic::lbm::stl::validate(data, "full", topology);

const geo::STLGeometry channel(std::move(topology), geo::FlowType::Internal);
geo::GeometryAnalysis3D<lat::LatticeType::D3Q19> result(0.5);
channel.analysis(result);
```

For external flow, pass an explicit `BoundingBox` pointer as the third
constructor argument. Use `geometry.bounds` to access STL bounds.

Full validation requires closed, consistently wound manifold surfaces. STL
input and validation are separate from flow interpretation; the input layer
reconstructs normals from vertex winding without choosing a fluid side.
See [STL input documentation](nticLBMSTLReader.md) for input-layer details;
the example above uses the `FacetTopology` overload required by geometry.

### Initialization and component nesting

Construction computes bounds, analyzes component containment, interprets the
flow region, normalizes active component orientation, and, for internal flow,
identifies boundary features.

Current component semantics are:

| Flow type | Root components, level 0 | Child components, level 1 |
| --- | --- | --- |
| Internal | Exactly one root; fluid inside it | Active internal obstacles; fluid outside each child |
| External | One or more active root obstacles; fluid outside them | Inactive enclosed components |

Deeper nesting, level 2 or greater, is currently rejected. Multiple roots are
rejected for internal flow. Root/child relationships are available through
`containment`; active state and fluid side are available through `flow`.

Advanced inspection methods include `countRoots()`,
`pointInComponent(componentID, point)`, and
`componentInComponent(innerID, outerID)`. `pointInComponent()` tests a closed
component's geometric interior, not the final composite fluid region.

`interpretFlow()` reapplies semantics and orientation if the true flow type is
changed deliberately. This is a mutating geometry operation, unlike temporary
complement analysis. External flow also requires a valid open box, and feature
metadata should be refreshed as needed after deliberate geometry changes.

### Normal convention

After flow interpretation, **active STL boundary normals point from fluid to
solid**. Vertex winding, facet normals, and component signed volume are
updated consistently:

- an internal root's normals point outward from its enclosed volume;
- an external root obstacle's normals point inward into the solid;
- an internal child obstacle's normals point inward into that obstacle.

Inactive components are not reoriented by this normalization. Primitive and
Boolean analysis currently do not expose a shared boundary-normal array.

Temporary complementary analysis preserves the original geometry, winding,
and stored normals. Q intersection tests are two-sided. A future consumer
requiring normals oriented relative to a complementary solving region must
account for the reversed effective side rather than modifying the original
model.

### Planar Faces and feature edges

For internal STL analysis, Face `0` is the Other Surface group. Connected
planar regions receive IDs `1..nPlanarFaces`. Detection uses coplanar adjacent
facets and a fixed seed plane; a single isolated facet is not promoted to its
own planar Face. This does not approximate a curved patch as a plane by using
a user-selected angular tolerance.

```cpp
channel.printBoundaryFeatureReport();

// On a non-const geometry: move selected original planar Face IDs to Other.
// channel.identifyBoundaryFeatures({2, 5}, 1.0e-6);
```

`identifyBoundaryFeatures(excludedFaceIDs, smallFaceAreaRatio)` supports
explicit exclusions. Excluded facets move to Face `0`; remaining planar IDs
are compacted. `planarFaces` retains both `faceID` and `originalFaceID`.
Small-area classification is a warning flag, not automatic exclusion.

`boundaryFeatures` also provides per-facet Face IDs and feature-edge records.
Feature edges are metadata, not additional Q Face groups. External STL Q
analysis currently uses one group, Face `0`.

Several unrelated curved surfaces can all belong to Other Surface. Importing
externally assigned facet/Face labels is not implemented yet.

## Boolean geometry

Boolean composition combines selected fluid regions. Operands must have the
same dimension. Primitive/primitive, primitive/STL, STL/STL, and nested Boolean
operands use the same analysis interface.

| Operation | Operand flow types | Result flow type | Expression |
| --- | --- | --- | --- |
| Union | Internal, Internal | Internal | `a \| b` or `a + b` |
| Difference | Internal, Internal | Internal | `a - b` |
| Intersection | Internal, Internal | Internal | `a & b` |
| Intersection | Internal, External, in either order | Internal | `a & b` |
| Intersection | External, External | External | `a & b` |

Other flow-type/operator combinations throw. External/external intersection
requires exactly matching open-box coordinates.

### Basic examples

```cpp
#include "boolean_geometry.hpp"

// 2D annulus: internal difference.
const geo::Circle outer({0.0, 0.0}, 5.0);
const geo::Circle inner({0.0, 0.0}, 2.0);
const auto annulus = outer - inner;
geo::GeometryAnalysis2D<lat::LatticeType::D2Q9> annulusResult(0.2);
annulus.analysis(annulusResult);

// 3D annular pipe: internal/external intersection.
const geo::BoundingBox openBox{{-6.0, -6.0, -9.0}, {6.0, 6.0, 9.0}};
const geo::Cylinder tube({0.0, 0.0, 0.0}, 5.0, 12.0, geo::Axis::Z);
const geo::Cylinder core(
    {0.0, 0.0, 0.0}, 2.0, 16.0, geo::Axis::Z,
    geo::FlowType::External, &openBox);
const auto pipe = tube & core;
geo::GeometryAnalysis3D<lat::LatticeType::D3Q19> pipeResult(0.2);
pipe.analysis(pipeResult);
```

Expressions store operand geometries by value. Editing the original operand
after composition does not update the stored copy. Move operands when suitable
to avoid copying large STL data.

### Nested expressions and shared grids

Use parentheses to make composition order explicit:

```cpp
const geo::Sphere a({-1.0, 0.0, 0.0}, 3.0);
const geo::Sphere b({ 1.0, 0.0, 0.0}, 3.0);
const geo::Sphere cut({0.0, 0.0, 0.0}, 1.0);
const auto combined = (a | b) - cut;
geo::GeometryAnalysis3D<lat::LatticeType::D3Q19> result(0.2);
combined.analysis(result);
```

Each subexpression must be one of the supported modes. The analysis recursively
evaluates children on the **same spacing, domain origin, and grid dimensions**,
then merges their results. It does not merge independently sampled scalar
arrays on different grids. Existing operand result objects do not have to be
provided or cached by the caller.

Default domains are selected as follows:

- internal union: the union bounding box;
- internal difference: the left operand's bounds;
- internal/internal intersection: the left operand's bounds, which contain
  the result even for an empty intersection;
- mixed intersection: the internal operand's bounds;
- external/external intersection: the common open box; `bounds` encloses both
  underlying obstacle bounds.

`bounds` is an analysis enclosure, not necessarily the tight bounding box of
the final fluid region. An explicit `result.domainBounds` overrides the
default analysis domain.

### Scalar merge and complement

Union follows these rules: any Interior operand gives Interior; otherwise any
Boundary operand gives Boundary; otherwise the result is Dry.

Intersection follows these rules: any Dry operand gives Dry; two Interior
operands give Interior; otherwise the result is Boundary.

Difference is evaluated as intersection with a temporary complement:

```text
A - B = A intersect complement(B)
```

The complement recomputes fluid-side classification and Q analysis. It is not
implemented by numerically replacing scalar values or by changing the original
operand's `flowType`. There is no public unary `!geometry` operator. Nested
complements are propagated through the expression using De Morgan's rules.

### Boundary provenance

For a result that remains Boundary, eligible boundary information is inherited
from its operand boundaries. Left operand groups come first; right operand
groups follow with a Face-ID offset equal to the left group's count. Nested
compositions preserve this recursive numbering, including empty groups.

A shared boundary cell can therefore appear in several source groups, while
each valid link Q belongs to one selected source group. Hidden interfaces must
not contribute final valid Q values. In particular, when adding guide tubes to
a closed STL channel by union, the covered original port caps become internal
interfaces rather than the final port boundaries.

Where both operands contribute, selection depends on the Boolean operation
and link geometry, not a universal minimum/maximum rule for saved Q values.
The current union implementation can query underlying link intersections to
reject hidden hits in nested expressions. It does not need the caller to
reconstruct or supply the expression's history separately.

## OpenMP and execution

Geometry analysis uses OpenMP for suitable grid, facet, and intersection work.
The caller can choose thread count through the runtime, for example
`OMP_NUM_THREADS`. Not every stage is parallel, and timings depend on grid
resolution, STL size, and Boolean expression depth. Nested analysis currently
recomputes child results; it is not a persistent result cache.

The project's remote test environment uses:

```bash
cd /home/liugeng/nticLBM
module load metax/metax-2.33_KylinSec-3.5.3-arm64
cmake_hpcc -S . -B build
make_hpcc -C build -j 16
```

All build artifacts remain under `build/`. Use `cmake_hpcc` and `make_hpcc` in
this environment even for the current CPU/OpenMP geometry implementation.

## Manual tests and ParaView output

There are 14 interactive manual applications. Primitive, STL and Boolean
applications prompt for grid spacing. The three multilevel applications prompt
for base spacing, maximum level and wall-band thickness. Geometry dimensions
and STL filenames are specified in their respective `.cpp` sources.

| Executable under `build/geometry` | Case |
| --- | --- |
| `primitive/test_rectangle` | 2D rectangle, internal and external |
| `primitive/test_circle` | 2D circle, internal and external |
| `primitive/test_box` | 3D box, internal and external |
| `primitive/test_sphere` | 3D sphere, internal and external |
| `primitive/test_cylinder` | Axis-aligned cylinder, internal and external |
| `stl/test_stl_area_analysis` | Branched-channel STL, internal and external |
| `boolean/test_boolean_annulus` | 2D internal difference |
| `boolean/test_boolean_annular_pipe` | 3D internal/external intersection |
| `boolean/test_boolean_branched_channel` | STL union with seven cylindrical guides |
| `boolean/test_boolean_crossed_cylinders` | Internal/internal intersection |
| `boolean/test_boolean_two_spheres` | External/external intersection |
| `multilevel/test_multilevel_annulus` | 2D annulus with wall bands and disconnected volume zones |
| `multilevel/test_multilevel_spherical_obstacle` | 3D Box minus sphere with a spherical wall band |
| `multilevel/test_multilevel_branched_channel` | Branched-channel STL with an Other Surface wall band |

For example:

```bash
printf '0.1\n' | OMP_NUM_THREADS=4 ./build/geometry/primitive/test_box
printf '1.0\n' | OMP_NUM_THREADS=128 ./build/geometry/stl/test_stl_area_analysis
printf '1.0\n' | OMP_NUM_THREADS=128 ./build/geometry/boolean/test_boolean_branched_channel
printf '1\n2\n0.5\n' | OMP_NUM_THREADS=4 ./build/geometry/multilevel/test_multilevel_annulus
```

Outputs are written under `tests/geometry_tests/outputs/test_xxx_vtks/`.
Primitive and standalone STL tests use `internal_vtks/` and `external_vtks/`
subdirectories. Boolean and multilevel cases have a determined result flow
type and write directly into their own named output directories. Re-running
a manual program recreates only its own output directory.

Single-level primitive, STL and Boolean files include:

- `geometry.vtk`: cell-centered scalar map on `STRUCTURED_POINTS`;
- `bounding_box.vtk`: geometry enclosure;
- `open_box.vtk`: explicit external analysis box, where applicable;
- `boundary_q.vtk`: Q intersection points with `Q`, `FaceID`, and `LinkID`
  as Point Data.

Multilevel cases instead write quadrilateral/hexahedral real cells as an
`UNSTRUCTURED_GRID` in `geometry.vtk`, separate exchange helpers in
`auxiliary.vtk`, and level-aware Q points in `boundary_q.vtk`. Their three
named output directories and refinement definitions are listed in the
multilevel section below. Manual applications are not invoked by the automatic
test runners.

`vtk_output_boundary.hpp`, under the test output helper directory, provides
`writeVTKBoundaryPoints(path, result)` for both 2D and 3D analysis. It writes
only intersection points, with vertex cells for display, and no connecting
lines. Coincident points with different source records remain separate.
It also prints per-Face entry/link counts and the valid Q range.

For ParaView:

- use the scalar map to inspect Dry/Boundary/Interior classification;
- when using Slice with Surface With Edges, disable **Triangulate the slice**
  to retain quadrilateral cut faces where the slice geometry permits;
- Contour reconstructs an interpolated surface, not the exact original STL;
  `1` is itself a categorical Boundary value, and `1.1` can be useful for
  visualization but is not a geometric accuracy criterion;
- use a cell-data Threshold on `[1,2]` to inspect solving-region cells;
- inspect `boundary_q.vtk` or the original geometry for actual geometric
  intersections. Q points are not boundary normal vectors.

Generated `*_vtks` directories and build outputs are not source artifacts and
are kept untracked.

## Multilevel grid generation

`multilevel/multilevel_grid.hpp` adds a static cell-centered quadtree (2D) or
octree (3D) forest. Level 0 is the base grid; each refinement divides a cell
into four/eight children and halves the spacing. Geometry construction,
selection, and physical scalar/Q analysis remain separate responsibilities.

```cpp
#include "multilevel/multilevel_grid.hpp"

const geo::Box physical({0,0,0}, {12,12,12});
const geo::Sphere selection({6,6,6}, 2);

geo::MultilevelOptions options;
options.maxLevel = 3;
geo::MultilevelGrid<3> grid(physical.bounds, 1.0, options);
const auto volume = geo::RefinementRule<3>::volume(geo::volumeRegion(selection), 2);
const auto wall = geo::RefinementRule<3>::surface(
    geo::fluidRegion(physical), 3, 0.5, {0,1}); // X-min/X-max fluid-side bands.
grid.generate(physical, {volume, wall});
```

2D uses `MultilevelGrid<2>` and `RefinementRule<2>`. The default lattices are
D2Q9 and D3Q19; a 3D lattice can be supplied as the second template parameter.
All root cells cover a spacing-aligned enclosure of the requested domain.
Root and child coordinates use integers; centers and spacing are derived.

### Refinement selections

`volumeRegion(primitive)` selects the primitive's interior independently of
its `FlowType`. Region values support unrestricted `|`, `&`, `-`, and unary
`!`, so nested and disconnected selections do not inherit BooleanGeometry's
flow-type restrictions. Overlapping rules apply the highest requested target
level. Full face/edge/corner 2:1 balancing can add transition refinement.

`volumeRegion(validatedTopology)` uses even/odd containment of closed STL
components, allowing arbitrary nesting without constructing an STLGeometry.
Its surface group is Face 0. `fluidRegion(STLGeometry)` instead uses the
geometry's active components, fluid-side interpretation, and Face numbering.
`fluidRegion(BooleanGeometry)` preserves the final fluid Boolean expression.
For a Boolean volume selector, `volumeRegion()` selects the internal fluid
result, or its complement for an external result (the combined obstacle).

Surface rules take a fluid region, a target level, a physical Euclidean
thickness, and optional Face IDs. Empty Face selections mean all surfaces.
Each band is measured from the original selected geometry, not from a previous
level's stair-step envelope. A cell bounding sphere supplies conservative
distance bounds. Near curved surfaces/corners and uncertain Boolean cuts this
can refine extra cells; it avoids dropping thin bands by testing centers only.
Boolean surfaces proven hidden over the query region do not select cells.
These are conservative selectors, not an exact curved-cell clipping algorithm.
Standalone open surface meshes are not accepted as volume selectors.

STL selectors build a shared read-only triangle BVH for distance queries and
borrow the original validated topology/geometry. Those objects must outlive
the rules and must not be moved or mutated while the rules are in use.
Primitive parameters are copied into their lightweight region representation.

### Topology, cell data, and access order

`nodes()` contains indexed tree metadata: level, integer coordinate, parent,
first child, active state, and cell index. Direct children use bit order
`bx + 2*by [+ 4*bz]` and are contiguous in the topology vector. A covered parent
has `cell == invalidIndex`; it does not occupy a field-data slot.

`cells()` is one SoA index space for real leaves and auxiliary exchange cells.
It stores source-node indices, role, auxiliary-child number, and physical
`CellType`. Dry leaves remain only in `nodes()` with `cell == invalidIndex`;
they occupy no persistent cell slot. `computes(cell)` is true for every real
cell and false for auxiliary cells. No LBM distribution arrays are allocated.

Cell data are ordered by level, role, physical type (Interior/Boundary),
then Morton spatial order. `levels()` exposes contiguous ranges. Real roles
are Ordinary, FineInterface, CoarseInterface, and BothInterfaces. A real leaf
may participate in both adjacent-level interfaces without being duplicated.
Physical types are independent of interface roles. Auxiliary groups follow
real cells, with each source coarse cell's four/eight auxiliary children
contiguous in child-number order.
The existing Dry type ranges remain empty. Always check for `invalidIndex`
when resolving a topology leaf to a cell index.

`coordinate(cell)`, `center(cell)`, and `level(cell)` resolve a cell index.
`nodeBox(node)` resolves a topology index. NodeIndex and CellIndex are different
index domains even though both use `std::size_t`; do not interchange them.
Real children may occupy different role/type blocks in cell storage: access
them through `nodes()[firstChild + childID].cell`, not parent-cell arithmetic.

The index domains are:

| Index | Addresses | Access rule |
| --- | --- | --- |
| `NodeIndex` | `nodes()` topology entries | Parent and first-child relations use this domain |
| `CellIndex` | All SoA arrays in `cells()` | One shared index space across levels, roles, and physical types |
| Boundary entry index | `boundaries().cells` and the corresponding row of `q` | Several entries may reference the same `CellIndex` |
| Interface fine-list index | `interfaces().fineCells` | The stored value is a `CellIndex`; the list position is not |

Within each level, real cells are stored in this exact order:

```text
Ordinary:        Interior, Boundary
FineInterface:   Interior, Boundary
CoarseInterface: Interior, Boundary
BothInterfaces:  Interior, Boundary
Auxiliary:      source coarse cell in Morton order, then child ID
```

Each real role/type subgroup uses Morton spatial order. This is the geometry
storage order; it does not prescribe the MPI partition algorithm. Auxiliary
cells belong to the fine-sized level returned by `level(cell)`, one greater
than their source coarse node's level. Their families remain contiguous even
when the auxiliary children have different physical types; they are not split
into Interior/Boundary subgroups.

All ranges below use absolute `CellIndex` endpoints and the half-open
convention `[begin, end)`:

| Range | Meaning |
| --- | --- |
| `levels()[l].all` | Every stored cell at level `l`, including auxiliary cells |
| `levels()[l].roles[r]` | The role `r` subgroup at level `l` |
| `levels()[l].types[r][0]` | Interior cells within real role `r` |
| `levels()[l].types[r][1]` | Boundary cells within real role `r` |
| `levels()[l].types[r][2]` | Empty Dry range, retained for interface compatibility |

The role indices are Ordinary=0, FineInterface=1, CoarseInterface=2,
BothInterfaces=3, Auxiliary=4. `types` only has entries for the four real roles.
Its type indices 0/1/2 are layout positions, **not `CellType` enum values**.
An absent subgroup has `begin == end`. The independently computing portion of
one level is `[roles[Ordinary].begin, roles[Auxiliary].begin)`.

Because role takes precedence over physical type, all Boundary cells of a
level do not form one contiguous interval. Iterate the Boundary range of each
real role. Likewise, physical boundary Face membership is a separate index
list; it does not duplicate or reorder the shared cell pool.

### Multilevel boundary membership and offsets

`boundaries()` preserves the original geometry's Face IDs, including empty
Faces. For Face `f`, iterate boundary entries in
`[offsets[f], offsets[f+1])`, and resolve each entry through `cells[entry]`:

```cpp
const auto& boundary = grid.boundaries();
for (auto entry = boundary.offsets[face];
     entry < boundary.offsets[face + 1]; ++entry)
{
    const auto cell = boundary.cells[entry];
    const auto level = grid.level(cell);
    // q[entry * (nStencils - 1) + (link - 1)], for link=1..nStencils-1.
    // Physical distance along that lattice link uses grid.spacing(level).
}
```

`offsets` counts boundary records, not distinct cells and not Q scalars.
It has `nBoundaries + 1` entries, starts at zero, and ends at
`boundaries().cells.size()`. Each Q row has `nStencils - 1` entries; invalid
links are `-1`, and valid Q values are in `(0, 1]`.

Multiple Face membership is retained. For example:

```text
offsets = [0, 2, 4]
cells   = [12, 27, 12, 45]
           Face 0  Face 1
```

Cell 12 occupies one cell-data slot but appears in two Face lists. Each Face
record contains only the links assigned to that Face. Every valid
`(CellIndex, link)` has a unique Face owner; multiple membership does not
duplicate a valid link's Q.

Construction analyzes each level at its own spacing, including retained
auxiliary positions. It does not copy coarse-grid Q onto fine cells. After
Dry removal and cell sorting, every boundary record is remapped to the new
`CellIndex`; records for discarded candidates are removed. Records are then
sorted by Face and final CellIndex, and Face offsets are rebuilt from counts.
Within one Face, CellIndex order therefore also follows level order.

Currently there is one Face-offset table combining all levels and both real
and auxiliary entries. A separate `level x Face` offset table is not exposed.
Use `grid.level(cell)` and `grid.cells().isAuxiliary(cell)` to select entries.

### Coarse/fine interface representation

After refinement and balance, a full coarse-cell neighbour ring outside the
refined region is considered for CoarseInterface. Cells proven wholly outside
the physical fluid region are excluded from candidate exchange families.
After scalar analysis, only wet coarse cells with wet fine partners retain
interface roles and exchange families. Each retained interface cell keeps its real coarse slot
and gets four/eight fine-sized auxiliary slots covering that same volume.
These are an exchange representation, not actual tree children or independent
solving volume. At a flat interface their width is one coarse spacing, or two
fine spacings. FineInterface marks a conservative two-fine-cell ring inside
the real fine region. The full box neighbourhood includes face/edge/corner
contacts. Helpers are shared when a coarse cell touches several interfaces.

`interfaces().blocks` maps coarse CellIndex to its contiguous auxiliary range
and a range in `interfaces().fineCells`. This is geometric connectivity, not
a streaming/interpolation matrix or a time-stepping implementation. Real coarse
and fine indices always refer to wet cells. Complete auxiliary families can
include Dry children; these are exchange metadata, not active fluid cells.

For each `InterfaceBlock`, `coarse` directly addresses the shared cell pool;
`auxiliary.begin/end` also address that pool and cover four/eight contiguous
children in bitwise child order. In contrast, `fine.begin/end` address the
`interfaces().fineCells` list: use `fineCells[i]` to obtain the actual fine
CellIndex. A fine cell may occur in several interface blocks when it interacts
with several coarse cells; its cell-data slot is still unique. Connectivity
is rebuilt using final cell indices after sorting.

### Physical analysis, edits, and partition handoff

Physical analysis reuses `geometry.analysis()` on occupied per-level tiles;
it does not allocate the complete finest-level Cartesian box. Boundary data
retain geometric Face IDs, offsets and entry-major Q, but reference the shared
CellIndex pool rather than duplicate coordinate arrays. `boundaries()` includes
both real and auxiliary boundary entries; use cell roles to distinguish them.

`refine(nodeIndex)` and `derefine(parentNodeIndex)` edit topology and invalidate
cell data/connectivity. Derefinement currently requires direct leaf children.
Call `finalize(geometry)` afterwards to balance, analyze, sort and rebuild
indices. Balance may refine neighbours again to preserve 2:1. Repeated
`generate()` applies rules to the existing forest rather than resetting it.
Indices into the cell pool must not be retained across these rebuilds.
Retired child topology nodes remain inactive; automatic AMR compaction and
distributed migration are not implemented.

The partition module can read the cell/level ranges, node-to-cell relations,
interface blocks, and `adjacentRealCells(cell)` (face/edge/corner neighbours).
Adjacency is derived on demand instead of retaining a dense N-by-stencil table.
Topology-only Dry leaves and auxiliary entries have no independent solving volume. No partition
rank, device ownership, or communication backend is assigned here.

`maxNodes` limits topology storage; `maxCells` limits the retained cell pool.
During construction, temporary analysis candidates can include Dry leaves;
they are removed before publishing indices, connectivity and boundary data.
`analysisTileWidth` bounds each
dense analysis tile; `analysisThreads` caps simultaneous tiles to bound memory.
Selection, topology creation, balancing and interface-role detection also use
OpenMP. Allocation and canonical sorting are deterministic serial stages.

The automatic acceptance suite is split into seven `.cpp` files, each with its
own executable and no case-selection arguments. Together they cover the same
15 test points:

| Executable under `build/geometry/multilevel` | Coverage |
| --- | --- |
| `test_multilevel_topology` | Grid snapping, indexed topology, coarse/fine interface bands, topology-only Dry storage |
| `test_multilevel_regions` | Disconnected/nested regions and Boolean wall visibility |
| `test_multilevel_primitive` | Curved walls, external flow, Cylinder Face bands |
| `test_multilevel_edits` | Refine/derefine, stale views, invalid input and capacity |
| `test_multilevel_lattice` | D3Q15/19/27 topology and physical Q |
| `test_multilevel_stl` | STL Face bands and arbitrary-depth STL refinement volumes |
| `test_multilevel_determinism` | Identical 1/4-thread topology, storage, interfaces and Q |

`test_multilevel_checks.hpp` supplies shared verification functions and fixtures
included by the seven automatic tests. Every test `.cpp` has its own executable;
shared helpers live in this header. `test_multilevel_acceptance` is a build-only
aggregate target, not an executable. Use:

```bash
make_hpcc -C build -j 16 test_multilevel_acceptance
bash run_geometry_multilevel_tests.sh
./build/geometry/multilevel/test_multilevel_stl
```

Like the existing runners, run the shell script from the project root. It
invokes the seven independent automatic applications directly. Interactive
visualization is split into three separate sources and executables:

| Manual executable under `build/geometry/multilevel` | Geometry and refinement | Output directory under `tests/geometry_tests/outputs` |
| --- | --- | --- |
| `test_multilevel_annulus` | 2D annulus, both wall bands and two disconnected level-1 rectangle zones | `test_multilevel_annulus_vtks/` |
| `test_multilevel_spherical_obstacle` | 3D Box minus sphere, spherical Face 6 wall band | `test_multilevel_spherical_obstacle_vtks/` |
| `test_multilevel_branched_channel` | Internal branched-channel STL, Face 0 wall band | `test_multilevel_branched_channel_vtks/` |

Each program prompts only for base spacing, maximum level (1..4), and fluid-side
wall-band thickness. Geometry dimensions and the STL filename are fixed in the
corresponding `.cpp`. Run each independently, for example:

```bash
./build/geometry/multilevel/test_multilevel_annulus
./build/geometry/multilevel/test_multilevel_spherical_obstacle
OMP_NUM_THREADS=128 ./build/geometry/multilevel/test_multilevel_branched_channel
```

Each recreates its own output directory, leaving the other cases intact. Real
quadrilateral/hexahedral cells go to `geometry.vtk`, overlapping auxiliary
cells to `auxiliary.vtk`, and Q points to `boundary_q.vtk`. The former combined
`test_multilevel` application and `test_multilevel_vtks` output directory have
been replaced by these named cases.

## Automatic tests and shell runners

There are 13 non-interactive automatic test executables. Paths in this table
are relative to `build/geometry/`. `<case>` and `<stl_file>` are command-line
arguments, not interactive input.

| Executable | Arguments | Coverage |
| --- | --- | --- |
| `primitive/test_primitive_acceptance` | None | Primitive Face groups, Q ownership, exact shared hits, historical scalar/Q regression, repeated 1/4-thread results |
| `stl/test_stl_geometry` | `<stl_file>` | Construction, preserved topology, planar/Other Face detection, exclusions, repeated detection and empty-geometry rejection |
| `stl/test_stl_containment` | `<case> <stl_file>` | Root/child relationships and unsupported-depth rejection |
| `stl/test_stl_flow` | `<case> <stl_file>` | Internal/external component semantics, solid-side normals, repeated interpretation and flow-type switching |
| `stl/test_stl_translation` | `<case> <stl_file>` | Internal/external translation and preserved topology, normals and volume |
| `stl/test_stl_boundary_acceptance` | None | Static STL intersection fixtures, primitive/STL comparison, Face/Q contracts and deterministic output |
| `multilevel/test_multilevel_topology` | None | Grid snapping, tree indices, coverage, balance and coarse/fine interface bands |
| `multilevel/test_multilevel_regions` | None | Disconnected/nested volumes and Boolean wall visibility |
| `multilevel/test_multilevel_primitive` | None | Curved walls, external fluid domain and Cylinder Face selection |
| `multilevel/test_multilevel_edits` | None | Refine/derefine, stale views, invalid input and capacity limits |
| `multilevel/test_multilevel_lattice` | None | D3Q15/19/27 octrees, cell storage, adjacency and scalar/Q comparison |
| `multilevel/test_multilevel_stl` | None | STL Face bands and arbitrary-depth STL refinement volumes |
| `multilevel/test_multilevel_determinism` | None | Identical 1/4-thread tree, cell order, interfaces and boundary Q |

The existing STL command-line case names are:

| Program | Supported `<case>` values |
| --- | --- |
| `test_stl_containment` | `single_root`, `two_roots`, `aabb_overlap`, `pseudo_containment`, `nested`, `two_children`, `two_nested_roots`, `three_levels` |
| `test_stl_flow` | `internal_root`, `internal_nested`, `internal_multiple_roots`, `external_roots`, `external_nested`, `internal_to_external` |
| `test_stl_translation` | `internal`, `external`, `repeated` |

Match each case with a suitable static fixture under
`tests/geometry_tests/stl_tests/stls/`; the STL runner contains the established
case/file combinations. For example:

```bash
STL_DIR=tests/geometry_tests/stl_tests/stls
./build/geometry/stl/test_stl_geometry "$STL_DIR/single_root.stl"
./build/geometry/stl/test_stl_containment nested "$STL_DIR/nested.stl"
./build/geometry/stl/test_stl_flow internal_to_external "$STL_DIR/nested.stl"
./build/geometry/stl/test_stl_translation repeated "$STL_DIR/two_children.stl"
```

Cases testing unsupported configurations pass when the expected exception is
raised. In particular, physical STL `three_levels` rejection and multilevel
`stl_nested` volume-selection acceptance test different semantics.

Two geometry shell runners are provided in the project root:

| Script | Programs invoked |
| --- | --- |
| `run_geometry_stl_tests.sh` | STL construction/feature detection, containment, flow interpretation and translation with the listed fixture combinations |
| `run_geometry_multilevel_tests.sh` | Seven independent multilevel automatic executables, covering 15 test points |

The primitive and STL boundary acceptance executables are run directly; the
two scripts do not invoke them. Boolean currently has the five manual cases
listed above, rather than a separate Boolean-only automatic executable.

After loading the remote build environment described earlier, the full
automatic geometry sequence is:

```bash
cmake_hpcc -S . -B build
make_hpcc -C build -j 16

./build/geometry/primitive/test_primitive_acceptance
./build/geometry/stl/test_stl_boundary_acceptance
OMP_NUM_THREADS=128 bash run_geometry_stl_tests.sh
bash run_geometry_multilevel_tests.sh
```

Run both scripts from the project root, like the existing project runners.
Each executable returns nonzero on failure. The scripts execute sequentially
without aggregating exit codes, so inspect each program's result rather than
using only the final shell status as the suite verdict. Multilevel acceptance
applications set their own 4-thread baseline; the determinism application
explicitly compares 1 and 4 threads.

Across geometry there are 27 executable tests: 14 interactive and 13 automatic.
`test_multilevel_checks.hpp` is shared verification code, not an additional
program; `test_multilevel_acceptance` is a CMake aggregate build target.

## Current scope

Each primitive/STL/Boolean analysis still uses a single uniform grid. The
multilevel module composes sparse per-level results and provides tree/interface
metadata; LBM evolution, adaptive refinement decisions, and partitioning are
not implemented here. Separate curved STL Face labels from
external software, a unified boundary-normal result array, and general
rotated primitive geometry remain outside the current interface.
