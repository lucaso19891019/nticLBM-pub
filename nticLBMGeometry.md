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
| Union | Internal, Internal | Internal | `a | b` or `a + b` |
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

Manual applications prompt for grid spacing. Geometry dimensions and the STL
filename are specified in their `.cpp` sources.

| Executable under `build/geometry` | Case |
| --- | --- |
| `primitive/test_rectangle`, `test_circle`, `test_box`, `test_sphere`, `test_cylinder` | Each primitive, internal and external |
| `stl/test_stl_area_analysis` | Branched-channel STL, internal and external |
| `boolean/test_boolean_annulus` | 2D internal difference |
| `boolean/test_boolean_annular_pipe` | 3D internal/external intersection |
| `boolean/test_boolean_branched_channel` | STL union with seven cylindrical guides |
| `boolean/test_boolean_crossed_cylinders` | Internal/internal intersection |
| `boolean/test_boolean_two_spheres` | External/external intersection |

For example:

```bash
printf '0.1\n' | OMP_NUM_THREADS=4 ./build/geometry/primitive/test_box
printf '1.0\n' | OMP_NUM_THREADS=128 ./build/geometry/stl/test_stl_area_analysis
printf '1.0\n' | OMP_NUM_THREADS=128 ./build/geometry/boolean/test_boolean_branched_channel
```

Outputs are written under `tests/geometry_tests/outputs/test_xxx_vtks/`.
Primitive and standalone STL tests use `internal_vtks/` and `external_vtks/`
subdirectories. Boolean cases already have a determined result flow type and
write directly into their own output directories. Re-running a manual program
recreates its output directory.

Files include:

- `geometry.vtk`: cell-centered scalar map on `STRUCTURED_POINTS`;
- `bounding_box.vtk`: geometry enclosure;
- `open_box.vtk`: explicit external analysis box, where applicable;
- `boundary_q.vtk`: Q intersection points with `Q`, `FaceID`, and `LinkID`
  as Point Data.

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

## Validation applications and current scope

`primitive/test_primitive_acceptance` checks Face grouping, link Q contracts,
shared-hit behavior, historical scalar/Q regression, and repeated 1/4-thread
results. `stl/test_stl_boundary_acceptance` checks STL link intersections,
Face grouping, and deterministic output on static fixtures.

Additional STL applications cover construction/feature detection
(`test_stl_geometry`), containment (`test_stl_containment`), flow interpretation
and orientation (`test_stl_flow`), and translation (`test_stl_translation`).
These take model paths and, where applicable, a named test case; see their
usage messages and sources under `tests/geometry_tests/stl_tests/`.

The current geometry module uses a single uniform grid per analysis. Nested
Boolean expressions are supported; multilevel grid nesting/refinement is not
implemented in this analysis pipeline. Separate curved STL Face labels from
external software, a unified boundary-normal result array, and general
rotated primitive geometry remain outside the current interface.
