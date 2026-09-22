#pragma once

// STL component orientation normalization.
//
// Responsibilities:
// - Interpret component roles using containment information and FlowType.
// - Determine which components are active for the selected flow type.
// - Determine the required orientation of active components.
// - Normalize component facet orientation when necessary.
//
// Internal flow:
// - Root components are active.
// - Direct child components are active.
// - Root orientation follows the internal-flow convention.
// - Child orientation follows the opposite convention.
//
// External flow:
// - Root components are active.
// - Components contained inside roots are inactive.
// - Root orientation follows the external-flow convention.
//
// This file does NOT:
// - Read or validate STL files.
// - Determine geometric containment.
// - Translate geometry.
// - Perform lattice discretization.
