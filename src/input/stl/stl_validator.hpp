#pragma once

#include "stl_reader.hpp"
#include "stl_topology.hpp"
#include <string>

namespace ntic::lbm::stl {

// Validate and normalize STL geometry.
//
// Validation is performed in several ordered stages:
//
// 1. Facet degeneracy check
// 2. Facet normal check
// 3. Duplicate facet check
// 4. Surface topology, winding, and connected-component validation
// 5. Component geometry calculation
//
// The input STLData may be modified during validation, for example when
// facet normals need to be normalized.
//
// A std::runtime_error is thrown if invalid STL geometry is detected.
void validate(
    STLData& data,
    const std::string& mode,
    FacetTopology& topology);

} // namespace ntic::lbm::stl