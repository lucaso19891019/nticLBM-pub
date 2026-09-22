#pragma once

// STL component containment representation.
//
// Responsibilities:
// - Represent containment relationships among validated STL components.
// - Identify root components.
// - Identify direct child components.
// - Represent multiple independent roots as a containment forest.
// - Provide the public interface for building the containment forest.
//
// Supported containment depth:
// - Root component.
// - Direct child component.
// - Deeper nesting is considered invalid.
//
// The containment relationship is purely geometric.
// It does NOT assign fluid/solid meaning by itself.
//
// This file does NOT:
// - Read or validate STL files.
// - Translate geometry.
// - Normalize orientation.
// - Discretize geometry onto a lattice.
