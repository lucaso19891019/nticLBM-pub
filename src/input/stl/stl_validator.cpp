#include "stl_validator.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace ntic::lbm::stl {

namespace {

constexpr double RELATIVE_LENGTH_TOLERANCE = 1.0e-7;
constexpr double COLLINEAR_TOLERANCE       = 1.0e-7;

//=============================================================================
// 1. Degenerate facets
//=============================================================================
//
// Validate that every STL facet defines a valid triangle.
//
// This stage will check conditions such as:
//
//   - repeated vertices within one facet
//   - three collinear vertices
//   - zero or near-zero triangle area
//
// Geometry-scale-dependent tolerances will be handled inside this stage.
//
// No topological information is required here.
//
void validateDegenerateFacets(const STLData& data)
{
    if (data.facets.empty()) {
        throw std::runtime_error(
            "Invalid STL geometry: no facets.");
    }

    //-------------------------------------------------------------------------
    // Determine the characteristic length of the geometry.
    //
    // The diagonal length of the STL bounding box is used as the global
    // reference scale for detecting nearly coincident vertices.
    //-------------------------------------------------------------------------

    STLVector minCoord = data.facets[0].vertices[0];
    STLVector maxCoord = data.facets[0].vertices[0];

    for (const auto& facet : data.facets) {
        for (const auto& vertex : facet.vertices) {
            for (std::size_t d = 0; d < 3; ++d) {
                minCoord[d] = std::min(minCoord[d], vertex[d]);
                maxCoord[d] = std::max(maxCoord[d], vertex[d]);
            }
        }
    }

    const double dx = maxCoord[0] - minCoord[0];
    const double dy = maxCoord[1] - minCoord[1];
    const double dz = maxCoord[2] - minCoord[2];

    const double scale =
        std::sqrt(dx * dx + dy * dy + dz * dz);

    if (!std::isfinite(scale) || scale <= 0.0) {
        throw std::runtime_error(
            "Invalid STL geometry: "
            "the geometry has zero or invalid extent.");
    }

    const double lengthTolerance =
        scale * RELATIVE_LENGTH_TOLERANCE;

    const double lengthToleranceSquared =
        lengthTolerance * lengthTolerance;


    //-------------------------------------------------------------------------
    // Validate every facet.
    //
    // A facet is considered degenerate if:
    //
    //   1. any of its three edges has nearly zero length, or
    //
    //   2. its three vertices are nearly collinear.
    //
    // Collinearity is measured using:
    //
    //        |e01 x e02|
    //        -----------
    //        |e01| |e02|
    //
    // which equals |sin(theta)| and is independent of triangle size.
    //-------------------------------------------------------------------------

    for (std::size_t i = 0; i < data.facets.size(); ++i) {

        const auto& v0 = data.facets[i].vertices[0];
        const auto& v1 = data.facets[i].vertices[1];
        const auto& v2 = data.facets[i].vertices[2];

        const double e01x = v1[0] - v0[0];
        const double e01y = v1[1] - v0[1];
        const double e01z = v1[2] - v0[2];

        const double e02x = v2[0] - v0[0];
        const double e02y = v2[1] - v0[1];
        const double e02z = v2[2] - v0[2];

        const double e12x = v2[0] - v1[0];
        const double e12y = v2[1] - v1[1];
        const double e12z = v2[2] - v1[2];

        const double e01Squared =
            e01x * e01x +
            e01y * e01y +
            e01z * e01z;

        const double e02Squared =
            e02x * e02x +
            e02y * e02y +
            e02z * e02z;

        const double e12Squared =
            e12x * e12x +
            e12y * e12y +
            e12z * e12z;


        //---------------------------------------------------------------------
        // Check for repeated or nearly coincident vertices.
        //---------------------------------------------------------------------

        if (e01Squared <= lengthToleranceSquared ||
            e02Squared <= lengthToleranceSquared ||
            e12Squared <= lengthToleranceSquared) {

            throw std::runtime_error(
                "Invalid STL geometry: "
                "facet " +
                std::to_string(i) +
                " contains coincident or nearly coincident vertices.");
        }


        //---------------------------------------------------------------------
        // Check for collinear or nearly collinear vertices.
        //---------------------------------------------------------------------

        const double nx =
            e01y * e02z - e01z * e02y;

        const double ny =
            e01z * e02x - e01x * e02z;

        const double nz =
            e01x * e02y - e01y * e02x;

        const double crossSquared =
            nx * nx +
            ny * ny +
            nz * nz;

        //
        // Instead of computing
        //
        //     sqrt(crossSquared) /
        //     (sqrt(e01Squared) * sqrt(e02Squared))
        //
        // compare the squared quantities directly.
        //

        const double collinearThreshold =
            COLLINEAR_TOLERANCE *
            COLLINEAR_TOLERANCE *
            e01Squared *
            e02Squared;

        if (crossSquared <= collinearThreshold) {
            throw std::runtime_error(
                "Invalid STL geometry: "
                "facet " +
                std::to_string(i) +
                " has collinear or nearly collinear vertices.");
        }
    }
}


//=============================================================================
// 2. Facet normals
//=============================================================================
//
// Validate the normal stored in each STL facet.
//
// For every facet, the geometric normal will be computed from the vertex
// winding:
//
//     (v1 - v0) x (v2 - v0)
//
// This stage will check:
//
//   - zero or invalid stored normals
//   - consistency between the stored STL normal and vertex winding
//   - normalization of the final facet normal
//
// The vertex winding remains the authoritative geometric orientation.
// The stored STL normal is treated as input information to be validated.
//
void validateFacetNormals(STLData& data)
{
    // TODO
}


//=============================================================================
// 3. Duplicate facets
//=============================================================================
//
// Detect repeated triangular facets.
//
// Two facets are considered duplicates when they represent the same
// geometric triangle, independent of the ordering of their three vertices.
//
// This stage will eventually use geometry-aware vertex comparison rather
// than relying on exact floating-point equality.
//
// Duplicate facets are rejected because they can corrupt:
//
//   - edge topology
//   - surface intersection tests
//   - inside/outside classification
//
void validateDuplicateFacets(STLData& data)
{
    // TODO
}


//=============================================================================
// 4. Surface topology, winding, and connected components
//=============================================================================
//
// Build the temporary surface-topology information required by several
// closely related validation steps.
//
// These operations are intentionally grouped together because they depend
// on the same vertex/edge/facet adjacency information.
//
// This stage will:
//
//   1. Build canonical vertex and edge representations.
//
//   2. Build facet adjacency through shared edges.
//
//   3. Validate surface closedness.
//      For a valid closed manifold surface, every edge must be shared by
//      exactly two facets.
//
//   4. Detect non-manifold topology.
//      An edge shared by more than two facets is invalid.
//
//   5. Validate local winding consistency.
//      Two facets sharing an edge must traverse that common edge in
//      opposite directions.
//
//   6. Identify connected surface components.
//
// Multiple disconnected closed components are allowed. They may represent
// separate solids, cavities, or nested geometric regions.
//
// All topology data created here is temporary validation data and is not
// part of the public STLData interface.
//
void validateTopologyWindingAndComponents(STLData& data)
{
    // TODO
}


//=============================================================================
// 5. Component nesting and orientation
//=============================================================================
//
// Validate the geometric relationship between closed connected components.
//
// This stage operates after each surface component is known to be locally
// closed, manifold, and consistently wound.
//
// It will determine:
//
//   - which components are geometrically inside other components
//   - the nesting depth of each component
//   - whether a component represents an outer shell, cavity, or nested solid
//   - whether the component orientation is consistent with that role
//
// Example:
//
//     outer shell        nesting depth 0
//         |
//         +-- cavity     nesting depth 1
//               |
//               +-- solid island
//                              nesting depth 2
//
// Orientation must alternate consistently with nesting depth so that the
// final STL geometry has an unambiguous inside/outside definition.
//
// If required, this stage may reverse the winding of an entire connected
// component and recompute its facet normals.
//
void validateNestingAndOrientation(STLData& data)
{
    // TODO
}

} // namespace


//=============================================================================
// STL validation pipeline
//=============================================================================
//
// Keep this function intentionally simple.
//
// Each call represents one logically complete validation stage. Internal
// implementation details such as tolerance computation, edge maps,
// adjacency construction, graph traversal, or orientation correction must
// remain inside the corresponding stage above.
//
void validate(STLData& data)
{
    validateDegenerateFacets(data);

    validateFacetNormals(data);

    validateDuplicateFacets(data);

    validateTopologyWindingAndComponents(data);

    validateNestingAndOrientation(data);
}

} // namespace ntic::lbm::stl