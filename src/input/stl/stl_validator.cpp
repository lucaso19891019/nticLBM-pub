#include "stl_validator.hpp"
#include "stl_topology.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

#include <array>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <queue>

namespace ntic::lbm::stl {

namespace {

constexpr double COLLINEAR_TOLERANCE       = 1.0e-7;

constexpr double RELATIVE_VERTEX_TOLERANCE = 1.0e-7;

constexpr double RELATIVE_VOLUME_TOLERANCE  = 1.0e-12;

//=============================================================================
// Spatial bin key
//=============================================================================

struct BinKey
{
    std::int64_t x;
    std::int64_t y;
    std::int64_t z;

    bool operator==(const BinKey& other) const noexcept
    {
        return
            x == other.x &&
            y == other.y &&
            z == other.z;
    }
};

struct BinKeyHash
{
    std::size_t operator()(const BinKey& key) const noexcept
    {
        const std::size_t hx =
            std::hash<std::int64_t>{}(key.x);

        const std::size_t hy =
            std::hash<std::int64_t>{}(key.y);

        const std::size_t hz =
            std::hash<std::int64_t>{}(key.z);

        return
            hx ^
            (hy << 1) ^
            (hz << 2);
    }
};

//=============================================================================
// Triangle key
//=============================================================================

struct TriangleKey
{
    std::array<std::size_t, 3> vertices;

    bool operator==(const TriangleKey& other) const noexcept
    {
        return vertices == other.vertices;
    }
};


struct TriangleKeyHash
{
    std::size_t operator()(const TriangleKey& key) const noexcept
    {
        const std::size_t h0 =
            std::hash<std::size_t>{}(key.vertices[0]);

        const std::size_t h1 =
            std::hash<std::size_t>{}(key.vertices[1]);

        const std::size_t h2 =
            std::hash<std::size_t>{}(key.vertices[2]);

        return
            h0 ^
            (h1 << 1) ^
            (h2 << 2);
    }
};

using SpatialBins =
    std::unordered_map<
        BinKey,
        std::vector<std::size_t>,
        BinKeyHash>;

struct EdgeKeyHash
{
    std::size_t operator()(const EdgeKey& key) const noexcept
    {
        const std::size_t h0 =
            std::hash<std::size_t>{}(key.v0);

        const std::size_t h1 =
            std::hash<std::size_t>{}(key.v1);

        return
            h0 ^
            (h1 << 1);
    }
};

using EdgeMap =
    std::unordered_map<
        EdgeKey,
        std::vector<EdgeUse>,
        EdgeKeyHash>;

//===============================================================================
// STL helpers
//===============================================================================

GeometryBounds computeGeometryBounds(const STLData& data)
{
    if (data.facets.empty()) {
        throw std::runtime_error(
            "Invalid STL geometry: no facets.");
    }

    GeometryBounds bounds;

    bounds.min = data.facets[0].vertices[0];
    bounds.max = data.facets[0].vertices[0];

    for (const auto& facet : data.facets) {
        for (const auto& vertex : facet.vertices) {
            for (std::size_t d = 0; d < 3; ++d) {

                bounds.min[d] =
                    std::min(bounds.min[d], vertex[d]);

                bounds.max[d] =
                    std::max(bounds.max[d], vertex[d]);
            }
        }
    }

    const double dx =
        bounds.max[0] - bounds.min[0];

    const double dy =
        bounds.max[1] - bounds.min[1];

    const double dz =
        bounds.max[2] - bounds.min[2];

    bounds.scale =
        std::sqrt(
            dx * dx +
            dy * dy +
            dz * dz);

    if (!std::isfinite(bounds.scale) ||
        bounds.scale <= 0.0) {

        throw std::runtime_error(
            "Invalid STL geometry: "
            "the geometry has zero or invalid extent.");
    }

    return bounds;
}

BinKey makeBinKey(
    const STLVector& vertex,
    const GeometryBounds& bounds,
    const double binSize)
{
    return {
        static_cast<std::int64_t>(
            std::floor(
                (vertex[0] - bounds.min[0]) /
                binSize)),

        static_cast<std::int64_t>(
            std::floor(
                (vertex[1] - bounds.min[1]) /
                binSize)),

        static_cast<std::int64_t>(
            std::floor(
                (vertex[2] - bounds.min[2]) /
                binSize))
    };
}

bool sameGeometricVertex(
    const STLVector& a,
    const STLVector& b,
    const double toleranceSquared)
{
    const double dx =
        a[0] - b[0];

    const double dy =
        a[1] - b[1];

    const double dz =
        a[2] - b[2];

    const double distanceSquared =
        dx * dx +
        dy * dy +
        dz * dz;

    return
        distanceSquared <= toleranceSquared;
}

std::size_t getOrCreateVertexID(
    const STLVector& vertex,
    const GeometryBounds& bounds,
    const double tolerance,
    const double toleranceSquared,
    std::vector<STLVector>& uniqueVertices,
    SpatialBins& bins)
{
    const BinKey baseKey =
        makeBinKey(
            vertex,
            bounds,
            tolerance);


    //-------------------------------------------------------------------------
    // Search the current spatial bin and all 26 neighboring bins.
    //
    // Vertices within the geometric tolerance may lie on opposite sides of
    // a bin boundary. Searching neighboring bins prevents such vertices from
    // being incorrectly treated as different geometric vertices.
    //-------------------------------------------------------------------------

    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {

                const BinKey neighborKey{
                    baseKey.x + dx,
                    baseKey.y + dy,
                    baseKey.z + dz
                };

                const auto binIt =
                    bins.find(neighborKey);

                if (binIt == bins.end()) {
                    continue;
                }

                for (const std::size_t vertexID :
                     binIt->second) {

                    if (sameGeometricVertex(
                            vertex,
                            uniqueVertices[vertexID],
                            toleranceSquared)) {

                        return vertexID;
                    }
                }
            }
        }
    }


    //-------------------------------------------------------------------------
    // No matching geometric vertex was found.
    // Register a new canonical vertex.
    //-------------------------------------------------------------------------

    const std::size_t newID =
        uniqueVertices.size();

    uniqueVertices.push_back(vertex);

    bins[baseKey].push_back(newID);

    return newID;
}

GeometricVertices buildGeometricVertices(
    const STLData& data,
    const GeometryBounds& bounds)
{
    const double vertexTolerance =
        bounds.scale * RELATIVE_VERTEX_TOLERANCE;

    const double vertexToleranceSquared =
        vertexTolerance * vertexTolerance;


    //-------------------------------------------------------------------------
    // Allocate the geometric vertex representation.
    //-------------------------------------------------------------------------

    GeometricVertices geometry;

    geometry.vertices.reserve(
        data.facets.size() * 3);

    geometry.facetVertexIDs.resize(
        data.facets.size());


    //-------------------------------------------------------------------------
    // Spatial bins used during vertex welding.
    //-------------------------------------------------------------------------

    SpatialBins bins;

    bins.reserve(
        data.facets.size() * 3);


    //-------------------------------------------------------------------------
    // Convert every STL vertex into a geometric vertex ID.
    //
    // Vertices within the geometry-relative tolerance are welded to the same
    // geometric vertex.
    //
    // The original vertex ordering of each facet is preserved.
    //-------------------------------------------------------------------------

    for (std::size_t i = 0;
         i < data.facets.size();
         ++i) {

        const auto& facet =
            data.facets[i];

        auto& vertexIDs =
            geometry.facetVertexIDs[i];

        for (std::size_t v = 0; v < 3; ++v) {

            vertexIDs[v] =
                getOrCreateVertexID(
                    facet.vertices[v],
                    bounds,
                    vertexTolerance,
                    vertexToleranceSquared,
                    geometry.vertices,
                    bins);
        }
    }

    return geometry;
}

EdgeKey makeEdgeKey(
    const std::size_t a,
    const std::size_t b)
{
    if (a < b) {
        return {
            a,
            b
        };
    }

    return {
        b,
        a
    };
}

bool edgeForward(
    const std::size_t a,
    const std::size_t b)
{
    return a < b;
}

void buildEdgeTopology(
    const GeometricVertices& geometry,
    EdgeMap& edgeMap)
{
    edgeMap.reserve(
        geometry.facetVertexIDs.size() * 3);


    for (std::size_t i = 0;
         i < geometry.facetVertexIDs.size();
         ++i) {


        const auto& vertices =
            geometry.facetVertexIDs[i];


        const std::array<
            std::pair<std::size_t,std::size_t>,
            3>
            edges =
            {{
                {vertices[0], vertices[1]},
                {vertices[1], vertices[2]},
                {vertices[2], vertices[0]}
            }};


        for (const auto& edge : edges) {

            const auto v0 =
                edge.first;

            const auto v1 =
                edge.second;


            const EdgeKey key =
                makeEdgeKey(
                    v0,
                    v1);


            edgeMap[key].push_back(
                {
                    i,
                    edgeForward(v0,v1)
                });
        }
    }
}

void validateEdges(
    const EdgeMap& edgeMap)
{
    for (const auto& item : edgeMap) {


        const auto& uses =
            item.second;


        //---------------------------------------------------------------------
        // Open boundary
        //---------------------------------------------------------------------

        if (uses.size() == 1) {

            throw std::runtime_error(
                "Invalid STL geometry: "
                "open boundary edge detected.");
        }


        //---------------------------------------------------------------------
        // Non-manifold edge
        //---------------------------------------------------------------------

        if (uses.size() > 2) {

            throw std::runtime_error(
                "Invalid STL geometry: "
                "non-manifold edge detected.");
        }


        //---------------------------------------------------------------------
        // Winding consistency
        //---------------------------------------------------------------------

        if (uses[0].forward ==
            uses[1].forward) {


            throw std::runtime_error(
                "Invalid STL geometry: "
                "inconsistent facet winding detected.");
        }
    }
}

void buildFacetAdjacency(
    const EdgeMap& edgeMap,
    std::vector<std::vector<std::size_t>>& adjacency)
{
    for (const auto& item : edgeMap) {

        const auto& uses =
            item.second;


        // validateEdges 已经保证 size == 2

        const auto a =
            uses[0].facetID;

        const auto b =
            uses[1].facetID;


        adjacency[a].push_back(b);
        adjacency[b].push_back(a);
    }
}

void buildComponents(
    const std::vector<std::vector<std::size_t>>& adjacency,
    STLComponents& components)
{
    const std::size_t n =
        adjacency.size();


    std::vector<bool> visited(
        n,
        false);


    for (std::size_t i = 0;
         i < n;
         ++i) {


        if (visited[i]) {
            continue;
        }


        std::vector<std::size_t>
            component;


        std::queue<std::size_t>
            queue;


        queue.push(i);

        visited[i] = true;


        while (!queue.empty()) {

            const auto current =
                queue.front();

            queue.pop();


            component.push_back(
                current);


            for (const auto next :
                 adjacency[current]) {


                if (!visited[next]) {

                    visited[next] = true;

                    queue.push(next);
                }
            }
        }


        STLComponent info;

        info.facets =
            std::move(component);


        components.push_back(
            std::move(info));
    }
}

void computeComponentBounds(
    const GeometricVertices& geometry,
    STLComponents& components)
{
    for(auto& component : components)
    {
        GeometryBounds bounds;


        bounds.min =
            geometry.vertices[
                geometry.facetVertexIDs[
                    component.facets[0]][0]
            ];

        bounds.max = bounds.min;


        for(auto facetID : component.facets)
        {
            for(auto vertexID :
                geometry.facetVertexIDs[facetID])
            {
                const auto& p =
                    geometry.vertices[vertexID];


                for(int d=0; d<3; d++)
                {
                    bounds.min[d] =
                        std::min(bounds.min[d],p[d]);

                    bounds.max[d] =
                        std::max(bounds.max[d],p[d]);
                }
            }
        }


        const double dx =
            bounds.max[0] - bounds.min[0];

        const double dy =
            bounds.max[1] - bounds.min[1];

        const double dz =
            bounds.max[2] - bounds.min[2];


        bounds.scale =
            std::sqrt(
                dx * dx +
                dy * dy +
                dz * dz);


        component.bounds =
            bounds;
    }
}

double computeSignedVolume(
    const GeometricVertices& geometry,
    const STLComponent& component)
{
    const double referenceX =
        0.5 *
        (component.bounds.min[0] +
         component.bounds.max[0]);

    const double referenceY =
        0.5 *
        (component.bounds.min[1] +
         component.bounds.max[1]);

    const double referenceZ =
        0.5 *
        (component.bounds.min[2] +
         component.bounds.max[2]);


    long double volume = 0.0L;


    for(const auto facetID : component.facets)
    {
        const auto& ids =
            geometry.facetVertexIDs[facetID];


        const auto& v0 =
            geometry.vertices[ids[0]];

        const auto& v1 =
            geometry.vertices[ids[1]];

        const auto& v2 =
            geometry.vertices[ids[2]];


        const long double x0 =
            static_cast<long double>(v0[0] - referenceX);

        const long double y0 =
            static_cast<long double>(v0[1] - referenceY);

        const long double z0 =
            static_cast<long double>(v0[2] - referenceZ);


        const long double x1 =
            static_cast<long double>(v1[0] - referenceX);

        const long double y1 =
            static_cast<long double>(v1[1] - referenceY);

        const long double z1 =
            static_cast<long double>(v1[2] - referenceZ);


        const long double x2 =
            static_cast<long double>(v2[0] - referenceX);

        const long double y2 =
            static_cast<long double>(v2[1] - referenceY);

        const long double z2 =
            static_cast<long double>(v2[2] - referenceZ);


        volume +=
            x0 *
            (y1 * z2 - z1 * y2)
            -
            y0 *
            (x1 * z2 - z1 * x2)
            +
            z0 *
            (x1 * y2 - y1 * x2);
    }


    return
        static_cast<double>(
            volume / 6.0L);
}

void validateComponentVolume(
    const STLComponent& component,
    const std::size_t componentID)
{
    const double scale =
        component.bounds.scale;


    if(!std::isfinite(scale) ||
       scale <= 0.0)
    {
        throw std::runtime_error(
            "Invalid STL geometry: "
            "component " +
            std::to_string(componentID) +
            " has zero or invalid geometric extent.");
    }


    const double volumeTolerance =
        RELATIVE_VOLUME_TOLERANCE *
        scale *
        scale *
        scale;


    if(!std::isfinite(component.signedVolume) ||
       std::abs(component.signedVolume) <= volumeTolerance)
    {
        throw std::runtime_error(
            "Invalid STL geometry: "
            "component " +
            std::to_string(componentID) +
            " has zero or nearly zero enclosed volume.");
    }
}


//////////////////////////////////////////////////////////////////////////////////////////
//=============================================================================
// 1. Degenerate facets
//=============================================================================
//
// Reject triangles that do not define a valid geometric surface element.
//
// This stage checks:
//
//   - coincident or nearly coincident vertices
//   - collinear or nearly collinear vertices
//
// The vertex tolerance is defined relative to the overall STL geometry scale.
//

void validateDegenerateFacets(
    const STLData& data,
    const GeometryBounds& bounds)
{
    const double vertexTolerance =
        bounds.scale * RELATIVE_VERTEX_TOLERANCE;

    const double vertexToleranceSquared =
        vertexTolerance * vertexTolerance;


    //-------------------------------------------------------------------------
    // Validate every facet independently.
    //
    // Error codes are stored per facet during the parallel region. Exceptions
    // are raised afterwards so that validation remains deterministic and the
    // lowest invalid facet ID is always reported.
    //
    //   0 : valid
    //   1 : coincident or nearly coincident vertices
    //   2 : collinear or nearly collinear vertices
    //-------------------------------------------------------------------------

    std::vector<unsigned char>
        errors(
            data.facets.size(),
            0);


    #pragma omp parallel for schedule(static)
    for (std::ptrdiff_t index = 0;
         index <
             static_cast<std::ptrdiff_t>(
                 data.facets.size());
         ++index) {

        const std::size_t i =
            static_cast<std::size_t>(index);


        const auto& v0 =
            data.facets[i].vertices[0];

        const auto& v1 =
            data.facets[i].vertices[1];

        const auto& v2 =
            data.facets[i].vertices[2];


        //---------------------------------------------------------------------
        // Construct the three triangle edges.
        //---------------------------------------------------------------------

        const double e01x = v1[0] - v0[0];
        const double e01y = v1[1] - v0[1];
        const double e01z = v1[2] - v0[2];

        const double e02x = v2[0] - v0[0];
        const double e02y = v2[1] - v0[1];
        const double e02z = v2[2] - v0[2];

        const double e12x = v2[0] - v1[0];
        const double e12y = v2[1] - v1[1];
        const double e12z = v2[2] - v1[2];


        //---------------------------------------------------------------------
        // Compute squared edge lengths.
        //---------------------------------------------------------------------

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
        // Check for coincident or nearly coincident vertices.
        //---------------------------------------------------------------------

        if (e01Squared <= vertexToleranceSquared ||
            e02Squared <= vertexToleranceSquared ||
            e12Squared <= vertexToleranceSquared) {

            errors[i] = 1;
            continue;
        }


        //---------------------------------------------------------------------
        // Compute the cross product:
        //
        //     (v1 - v0) x (v2 - v0)
        //---------------------------------------------------------------------

        const double nx =
            e01y * e02z -
            e01z * e02y;

        const double ny =
            e01z * e02x -
            e01x * e02z;

        const double nz =
            e01x * e02y -
            e01y * e02x;

        const double crossSquared =
            nx * nx +
            ny * ny +
            nz * nz;


        //---------------------------------------------------------------------
        // Check for collinear or nearly collinear vertices.
        //---------------------------------------------------------------------

        const double collinearThreshold =
            COLLINEAR_TOLERANCE *
            COLLINEAR_TOLERANCE *
            e01Squared *
            e02Squared;

        if (crossSquared <= collinearThreshold) {
            errors[i] = 2;
        }
    }


    //-------------------------------------------------------------------------
    // Report the first invalid facet in STL order.
    //-------------------------------------------------------------------------

    for (std::size_t i = 0;
         i < errors.size();
         ++i) {

        if (errors[i] == 1) {

            throw std::runtime_error(
                "Invalid STL geometry: "
                "facet " +
                std::to_string(i) +
                " contains coincident or nearly coincident vertices.");
        }

        if (errors[i] == 2) {

            throw std::runtime_error(
                "Invalid STL geometry: "
                "facet " +
                std::to_string(i) +
                " has collinear or nearly collinear vertices.");
        }
    }
}

//=============================================================================
// 2. Facet normal reconstruction
//=============================================================================
//
// Reconstruct the normal of each STL facet directly from its vertex winding.
//
// For every facet, the geometric normal is computed as:
//
//     (v1 - v0) x (v2 - v0)
//
// The vertex winding is the authoritative geometric orientation.
//
// The normal stored in the STL file is not used to determine facet
// orientation. It may be non-unit, zero, reversed, or otherwise unreliable.
//
// Degenerate facets have already been rejected in the previous stage, so the
// reconstructed geometric normal is guaranteed to have non-zero magnitude.
//
// The final facet normal is normalized and stored back into STLData.
//

void reconstructFacetNormals(STLData& data)
{
    #pragma omp parallel for schedule(static)
    for (std::ptrdiff_t index = 0;
         index <
             static_cast<std::ptrdiff_t>(
                 data.facets.size());
         ++index) {

        const std::size_t i =
            static_cast<std::size_t>(
                index);

        auto& facet =
            data.facets[i];


        //---------------------------------------------------------------------
        // Compute the geometric normal from the vertex winding.
        //
        //     n = (v1 - v0) x (v2 - v0)
        //
        // Degenerate facets have already been rejected by Step 1.
        //---------------------------------------------------------------------

        const auto& v0 =
            facet.vertices[0];

        const auto& v1 =
            facet.vertices[1];

        const auto& v2 =
            facet.vertices[2];


        const double e01x =
            v1[0] - v0[0];

        const double e01y =
            v1[1] - v0[1];

        const double e01z =
            v1[2] - v0[2];


        const double e02x =
            v2[0] - v0[0];

        const double e02y =
            v2[1] - v0[1];

        const double e02z =
            v2[2] - v0[2];


        const double geometricNx =
            e01y * e02z -
            e01z * e02y;

        const double geometricNy =
            e01z * e02x -
            e01x * e02z;

        const double geometricNz =
            e01x * e02y -
            e01y * e02x;


        const double geometricNorm =
            std::sqrt(
                geometricNx * geometricNx +
                geometricNy * geometricNy +
                geometricNz * geometricNz);


        //---------------------------------------------------------------------
        // Replace the stored STL normal with the normalized geometric normal.
        //---------------------------------------------------------------------

        facet.normal[0] =
            geometricNx / geometricNorm;

        facet.normal[1] =
            geometricNy / geometricNorm;

        facet.normal[2] =
            geometricNz / geometricNorm;
    }
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
// Geometric vertex IDs have already been constructed using tolerance-based
// vertex welding. This stage therefore operates only on integer vertex IDs.
//
// Duplicate facets are rejected because they can corrupt:
//
//   - edge topology
//   - surface intersection tests
//   - inside/outside classification
//

void validateDuplicateFacets(
    const GeometricVertices& geometry)
{
    //-------------------------------------------------------------------------
    // Canonical triangle registry.
    //-------------------------------------------------------------------------

    std::unordered_set<
        TriangleKey,
        TriangleKeyHash>
        triangleKeys;

    triangleKeys.reserve(
        geometry.facetVertexIDs.size());


    //-------------------------------------------------------------------------
    // Process every facet.
    //-------------------------------------------------------------------------

    for (std::size_t i = 0;
         i < geometry.facetVertexIDs.size();
         ++i) {

        // Copy the vertex IDs because the original ordering stored in
        // GeometricVertices must be preserved for later topology and
        // winding validation.

        auto vertexIDs =
            geometry.facetVertexIDs[i];


        //---------------------------------------------------------------------
        // Canonicalize the triangle.
        //
        // Sorting makes the key independent of:
        //
        //   - cyclic vertex ordering
        //   - reversed winding
        //
        // Therefore all six permutations of the same three geometric
        // vertices generate exactly the same TriangleKey.
        //---------------------------------------------------------------------

        std::sort(
            vertexIDs.begin(),
            vertexIDs.end());

        const TriangleKey key{
            vertexIDs
        };


        //---------------------------------------------------------------------
        // Detect duplicate facet.
        //---------------------------------------------------------------------

        const bool inserted =
            triangleKeys.insert(key).second;

        if (!inserted) {

            throw std::runtime_error(
                "Invalid STL geometry: "
                "duplicate facet detected at facet " +
                std::to_string(i) + ".");
        }
    }
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
// The welded geometric vertices, facet adjacency, and connected components
// are retained in FacetTopology for later geometry processing.
//
void validateTopologyWindingAndComponents(
    FacetTopology& topology)
{
    //-------------------------------------------------------------------------
    // Build edge topology
    //-------------------------------------------------------------------------

    EdgeMap edgeMap;

    buildEdgeTopology(
        topology.geometry,
        edgeMap);


    //-------------------------------------------------------------------------
    // Validate manifold and winding
    //-------------------------------------------------------------------------

    validateEdges(
        edgeMap);


    //-------------------------------------------------------------------------
    // Build facet adjacency graph
    //-------------------------------------------------------------------------

    topology.adjacency.clear();

    topology.adjacency.resize(
        topology.geometry.facetVertexIDs.size());


    buildFacetAdjacency(
        edgeMap,
        topology.adjacency);


    //-------------------------------------------------------------------------
    // Find connected components
    //-------------------------------------------------------------------------

    topology.components.clear();

    buildComponents(
        topology.adjacency,
        topology.components);
}


//=============================================================================
// 5. Component geometry
//=============================================================================
//
// Validate basic geometric properties of each closed surface component.
//
// The connected surface components are constructed in the previous topology
// validation stage.  This stage evaluates each component independently.
//
// For every component, this stage will:
//
//   - compute its geometric bounding box
//   - compute its signed enclosed volume
//   - reject components with zero or nearly zero enclosed volume
//
// The volume tolerance is defined relative to the geometric scale of each
// individual component.
//
// The signed-volume sign is retained as geometric information only. It is
// not used here to determine or correct the global surface orientation.
//
// More advanced geometric operations, including:
//
//   - global orientation analysis and correction
//   - component nesting
//   - inside/outside classification
//
// belong to the geometry module and are intentionally not performed here.
//
void validateComponentConsistency(
    FacetTopology& topology)
{
    computeComponentBounds(
        topology.geometry,
        topology.components);


    for(std::size_t i = 0;
        i < topology.components.size();
        ++i)
    {
        auto& component =
            topology.components[i];


        component.signedVolume =
            computeSignedVolume(
                topology.geometry,
                component);


        validateComponentVolume(
            component,
            i);
    }
}

} // namespace


//=============================================================================
// Public validation interface
//=============================================================================

void validate(
    STLData& data,
    const std::string& mode,
    FacetTopology& topology)
{
    //-------------------------------------------------------------------------
    // Validate mode.
    //-------------------------------------------------------------------------

    if (mode != "test" &&
        mode != "full") {

        throw std::invalid_argument(
            "Invalid STL validation mode: \"" +
            mode +
            "\". Expected \"test\" or \"full\".");
    }

    topology = FacetTopology{};

    //-------------------------------------------------------------------------
    // Compute geometry bounds once.
    //-------------------------------------------------------------------------

    const GeometryBounds bounds =
        computeGeometryBounds(data);


    //-------------------------------------------------------------------------
    // 1. Degenerate facets
    //-------------------------------------------------------------------------

    validateDegenerateFacets(
        data,
        bounds);


    //-------------------------------------------------------------------------
    // 2. Facet normal reconstruction
    //-------------------------------------------------------------------------

    reconstructFacetNormals(
        data);


    //-------------------------------------------------------------------------
    // Build the tolerance-welded geometric vertex representation.
    //-------------------------------------------------------------------------

    topology.geometry =
    buildGeometricVertices(
        data,
        bounds);


    //-------------------------------------------------------------------------
    // 3. Duplicate facets
    //-------------------------------------------------------------------------

    validateDuplicateFacets(
        topology.geometry);


    //-------------------------------------------------------------------------
    // Test mode stops after the facet-level validation stages.
    //-------------------------------------------------------------------------

    if (mode == "test") {
        return;
    }

    //-------------------------------------------------------------------------
    // 4. Topology, winding, and connected components
    //-------------------------------------------------------------------------

    validateTopologyWindingAndComponents(
        topology);


    //-------------------------------------------------------------------------
    // 5. Component geometry
    //-------------------------------------------------------------------------

    validateComponentConsistency(
        topology);
}

} // namespace ntic::lbm::stl
