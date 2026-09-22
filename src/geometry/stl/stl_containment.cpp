// STL component containment analysis.
//
// Responsibilities:
// - Perform broad-phase component bounding-box tests.
// - Determine candidate component containment relationships.
// - Perform geometric inside/outside tests where required.
// - Build the component containment forest.
// - Identify roots and direct children.
// - Reject unsupported nesting deeper than root + child.
//
// Performance:
// - Use inexpensive bounding-box rejection before expensive geometry tests.
// - Use OpenMP where independent geometric queries are sufficiently large.
// - Avoid unnecessary copies of STL topology data.
//
// The result describes geometric containment only.
// Flow semantics are handled separately.
//
// This file does NOT:
// - Translate STL geometry.
// - Decide internal/external fluid regions.
// - Normalize facet orientation.
// - Perform lattice discretization.
