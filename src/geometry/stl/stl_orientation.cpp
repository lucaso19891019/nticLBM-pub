// STL component orientation normalization implementation.
//
// Responsibilities:
// - Determine current component orientation.
// - Determine target orientation from FlowType and containment role.
// - Flip complete components when their orientation does not match the target.
// - Update affected topology data consistently.
// - Preserve topology connectivity and component membership.
//
// Orientation changes are performed at component level.
// The STL validator has already established consistent winding within each
// validated component, so orientation normalization must not independently
// guess or flip individual facets.
//
// Performance:
// - Use OpenMP for sufficiently large independent facet operations.
// - Avoid deep copies of topology arrays.
//
// This file does NOT:
// - Determine component containment.
// - Translate geometry.
// - Perform lattice discretization.
