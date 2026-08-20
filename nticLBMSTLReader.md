# STL Input

The STL input module reads triangular surface geometry from STL files and
validates the geometry before it is passed to the geometry processing stage.

The module is designed for closed surface meshes used as simulation geometry.
Validation is performed in several stages so that malformed or topologically
invalid STL files can be detected before further geometry processing.


## Supported STL Geometry

The STL input module currently supports triangular surface meshes with the
following requirements:

- each facet must contain three valid vertices
- facet normals must be consistent with the corresponding vertex winding
- duplicate facets are not allowed
- surfaces used in full validation must be closed
- surface edges must be manifold
- neighboring facets must have consistent winding
- multiple disconnected closed surface components are supported

The STL validator does not currently perform geometry repair. Operations such
as global normal correction, component nesting, and inside/outside
classification belong to the geometry processing module.


## STL Reading

An STL file is read into an `STLData` object:

```cpp
auto data = ntic::lbm::stl::read(filename);
```

The reader extracts the triangular facets, including:

- three vertices for each facet
- the normal stored in the STL file

The stored STL normal is retained as input information and is subsequently
checked during validation.


## Validation

Validation is performed through:

```cpp
ntic::lbm::stl::STLComponents components;

ntic::lbm::stl::validate(
    data,
    mode,
    components);
```

Two validation modes are currently available:

```text
test
full
```

`test` performs the basic validation stages used for lightweight STL checks.

`full` performs the complete validation sequence, including surface topology
and component geometry analysis.


## 1. Facet Geometry

Each triangular facet is checked for basic geometric validity.

The validator detects:

- coincident vertices
- nearly coincident vertices
- collinear vertices
- nearly collinear vertices
- degenerate triangular facets

The geometric tolerance is scaled according to the overall size of the STL
geometry so that the validation is not tied to a particular coordinate scale.


## 2. Facet Normals

For each facet, a geometric normal is computed from the vertex winding:

```text
(v1 - v0) x (v2 - v0)
```

The normal stored in the STL file is checked against this geometric normal.

The validator detects:

- zero or invalid stored normals
- stored normals inconsistent with vertex winding

Non-unit stored normals are accepted when their direction is valid.

After validation, the facet normal is reconstructed from the triangle geometry
and normalized.


## 3. Duplicate Facets

Repeated triangular facets are detected geometrically.

Two facets are considered duplicates when they represent the same triangle,
independent of the ordering of their three vertices.

The duplicate-facet check therefore detects cases including:

- identical facets
- cyclic permutations of the same facet
- reversed vertex ordering
- duplicate facets containing small coordinate differences within the
  geometric tolerance

Geometrically distinct nearby facets remain valid.


## 4. Surface Topology

Full validation constructs the surface topology from the triangular facets.

The topology validation checks that:

- every surface edge is shared by exactly two facets
- no open boundary edge exists
- no non-manifold edge exists
- neighboring facets use their shared edge in opposite directions
- facet winding is consistent across each connected surface

An edge used by only one facet indicates an open surface.

An edge used by more than two facets indicates a non-manifold surface.

Disconnected closed surfaces are allowed. Each connected closed surface is
identified as an individual STL component.


## 5. Component Geometry

After topology validation, basic geometric information is computed for each
connected surface component.

The current component information includes:

- the facets belonging to the component
- the component bounding box
- the signed enclosed volume

For example, an STL containing two disconnected closed objects produces two
separate components.

The signed volume is retained as component geometry information. Its sign is
not currently used to determine or correct the global surface orientation.


## Validation Errors

Invalid STL geometry causes validation to terminate with an error describing
the detected problem.

Typical errors include:

```text
Invalid STL geometry: facet 0 contains coincident or nearly coincident vertices.

Invalid STL geometry: facet 0 has collinear or nearly collinear vertices.

Invalid STL geometry: facet 0 has a zero or invalid stored normal.

Invalid STL geometry: facet 0 has a stored normal inconsistent with its vertex winding.

Invalid STL geometry: duplicate facet detected at facet 1.

Invalid STL geometry: open boundary edge detected.

Invalid STL geometry: non-manifold edge detected.

Invalid STL geometry: inconsistent facet winding detected.
```


## Validation Modes

The validation mode controls how much of the validation pipeline is executed.

### `test`

The `test` mode is intended for lightweight validation and unit testing of the
basic STL checks.

It does not require the input surface to form a complete closed solid.

This allows individual triangles and small facet collections to be used when
testing facet geometry, normals, and duplicate detection.


### `full`

The `full` mode performs the complete STL validation pipeline.

In addition to the basic facet checks, it performs:

- surface topology validation
- closed-surface validation
- manifold validation
- winding consistency validation
- connected-component construction
- component geometry analysis

STL files intended for simulation geometry should use `full` validation.


## Current Scope

The STL input module is responsible for reading and validating STL surface
geometry.

The following operations are intentionally outside the STL input module:

- global surface orientation correction
- STL geometry repair
- component nesting analysis
- inside/outside classification
- geometry Boolean operations
- simulation-domain point generation
- lattice neighbor construction

These operations are handled by the geometry processing stages built on top of
the validated STL input.