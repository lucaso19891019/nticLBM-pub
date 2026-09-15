#include "stl_validator.hpp"
#include "stl_topology.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
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

// An equilateral triangle has quality = 1. Increasingly elongated triangles
// approach quality = 0.
//
// Low triangle quality is reported as a warning because highly elongated
// facets may occur in otherwise valid STL geometry.
//
// Unusually large facet area is also reported as a warning because strongly
// non-uniform facet sizes may reduce the efficiency of the later
// centroid-based spatial search.

constexpr double MIN_FACET_QUALITY =
    0.2;

constexpr double MAX_FACET_AREA_RATIO =
    4.0;

//=============================================================================
// Spatial bin key
//=============================================================================

struct BinKey
{
    std::int64_t x;
    std::int64_t y;
    std::int64_t z;

    bool operator==(
        const BinKey& other) const noexcept
    {
        return
            x == other.x &&
            y == other.y &&
            z == other.z;
    }


    bool operator<(
        const BinKey& other) const noexcept
    {
        if (x != other.x) {
            return x < other.x;
        }

        if (y != other.y) {
            return y < other.y;
        }

        return z < other.z;
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

//=============================================================================
// Facet spatial cells
//=============================================================================
//
// Flat entries used to organize facets into spatial cells.
//
// A facet may generate multiple entries because its centroid-radius search
// region can overlap multiple cells. Entries are generated independently in
// parallel and later sorted by cell key.
//
// Keeping the representation flat avoids concurrent insertion into shared
// hash tables and provides contiguous facet ranges for subsequent candidate
// pair generation.
//

struct FacetCellEntry
{
    BinKey key;

    std::size_t facetID;
};


using FacetCellEntries =
    std::vector<FacetCellEntry>;

bool facetCellEntryLess(
    const FacetCellEntry& a,
    const FacetCellEntry& b) noexcept
{
    if (a.key < b.key) {
        return true;
    }

    if (b.key < a.key) {
        return false;
    }

    return a.facetID < b.facetID;
}

//=============================================================================
// Flat edge
//=============================================================================
//
// Flat edge representation used for parallel topology construction.
//
// Each STL facet contributes exactly three FlatEdge records. The canonical
// EdgeKey identifies the undirected geometric edge, while EdgeUse retains
// the facet ID and the original traversal direction.
//

struct FlatEdge
{
    EdgeKey key;
    EdgeUse use;
};

//===============================================================================
// STL helpers
//===============================================================================

GeometryBounds computeGeometryBounds(const STLData& data)
{
    if (data.facets.empty()) {
        throw std::runtime_error(
            "Invalid STL geometry: no facets.");
    }


    //-------------------------------------------------------------------------
    // Initialize the global bounds from the first STL vertex.
    //-------------------------------------------------------------------------

    double minX =
        data.facets[0].vertices[0][0];

    double minY =
        data.facets[0].vertices[0][1];

    double minZ =
        data.facets[0].vertices[0][2];

    double maxX = minX;
    double maxY = minY;
    double maxZ = minZ;


    //-------------------------------------------------------------------------
    // Compute the global STL bounding box.
    //
    // Every facet and vertex can be processed independently. OpenMP min/max
    // reductions combine the thread-local bounds into the final global bounds.
    //-------------------------------------------------------------------------

    #pragma omp parallel for \
        reduction(min:minX,minY,minZ) \
        reduction(max:maxX,maxY,maxZ) \
        schedule(static)

    for (std::ptrdiff_t index = 0;
         index <
             static_cast<std::ptrdiff_t>(
                 data.facets.size());
         ++index) {

        const auto& facet =
            data.facets[
                static_cast<std::size_t>(index)];


        for (std::size_t v = 0;
             v < 3;
             ++v) {

            const auto& vertex =
                facet.vertices[v];


            minX =
                std::min(
                    minX,
                    vertex[0]);

            minY =
                std::min(
                    minY,
                    vertex[1]);

            minZ =
                std::min(
                    minZ,
                    vertex[2]);


            maxX =
                std::max(
                    maxX,
                    vertex[0]);

            maxY =
                std::max(
                    maxY,
                    vertex[1]);

            maxZ =
                std::max(
                    maxZ,
                    vertex[2]);
        }
    }


    //-------------------------------------------------------------------------
    // Store the final bounds.
    //-------------------------------------------------------------------------

    GeometryBounds bounds;

    bounds.min = {
        minX,
        minY,
        minZ
    };

    bounds.max = {
        maxX,
        maxY,
        maxZ
    };


    //-------------------------------------------------------------------------
    // Compute the characteristic geometry scale from the bounding-box
    // diagonal.
    //-------------------------------------------------------------------------

    const double dx =
        maxX - minX;

    const double dy =
        maxY - minY;

    const double dz =
        maxZ - minZ;


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

//=============================================================================
// Facet geometry
//=============================================================================

double computeFacetGeometry(
    const STLData& data,
    std::vector<FacetGeometry>& facetGeometry)
{
    const std::size_t facetCount =
        data.facets.size();


    facetGeometry.resize(
        facetCount);


    double totalArea =
        0.0;


    //-------------------------------------------------------------------------
    // Compute all per-facet geometric quantities in parallel.
    //
    // Each iteration writes exclusively to facetGeometry[i]. The total facet
    // area is accumulated through an OpenMP reduction.
    //-------------------------------------------------------------------------

    #pragma omp parallel for \
        reduction(+:totalArea) \
        schedule(static)

    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                facetCount);
        ++index)
    {
        const std::size_t facetID =
            static_cast<std::size_t>(
                index);


        const auto& facet =
            data.facets[facetID];

        const auto& v0 =
            facet.vertices[0];

        const auto& v1 =
            facet.vertices[1];

        const auto& v2 =
            facet.vertices[2];


        //---------------------------------------------------------------------
        // Compute facet centroid.
        //---------------------------------------------------------------------

        const STLVector centroid =
        {
            (v0[0] + v1[0] + v2[0]) / 3.0,
            (v0[1] + v1[1] + v2[1]) / 3.0,
            (v0[2] + v1[2] + v2[2]) / 3.0
        };


        //---------------------------------------------------------------------
        // Compute the three squared edge lengths.
        //---------------------------------------------------------------------

        const double e01x =
            v1[0] - v0[0];

        const double e01y =
            v1[1] - v0[1];

        const double e01z =
            v1[2] - v0[2];


        const double e12x =
            v2[0] - v1[0];

        const double e12y =
            v2[1] - v1[1];

        const double e12z =
            v2[2] - v1[2];


        const double e20x =
            v0[0] - v2[0];

        const double e20y =
            v0[1] - v2[1];

        const double e20z =
            v0[2] - v2[2];


        const double e01Squared =
            e01x * e01x +
            e01y * e01y +
            e01z * e01z;

        const double e12Squared =
            e12x * e12x +
            e12y * e12y +
            e12z * e12z;

        const double e20Squared =
            e20x * e20x +
            e20y * e20y +
            e20z * e20z;


        //---------------------------------------------------------------------
        // Compute facet area from:
        //
        //     0.5 * |(v1 - v0) x (v2 - v0)|
        //---------------------------------------------------------------------

        const double e02x =
            v2[0] - v0[0];

        const double e02y =
            v2[1] - v0[1];

        const double e02z =
            v2[2] - v0[2];


        const double nx =
            e01y * e02z -
            e01z * e02y;

        const double ny =
            e01z * e02x -
            e01x * e02z;

        const double nz =
            e01x * e02y -
            e01y * e02x;


        const double area =
            0.5 *
            std::sqrt(
                nx * nx +
                ny * ny +
                nz * nz);


        //---------------------------------------------------------------------
        // Compute dimensionless triangle quality.
        //
        // An equilateral triangle has quality = 1. The value approaches zero
        // as the triangle becomes increasingly elongated.
        //---------------------------------------------------------------------

        const double quality =
            4.0 *
            std::sqrt(3.0) *
            area /
            (e01Squared +
             e12Squared +
             e20Squared);


        //---------------------------------------------------------------------
        // Compute the maximum distance from the centroid to the three
        // vertices.
        //---------------------------------------------------------------------

        const double r0x =
            v0[0] - centroid[0];

        const double r0y =
            v0[1] - centroid[1];

        const double r0z =
            v0[2] - centroid[2];


        const double r1x =
            v1[0] - centroid[0];

        const double r1y =
            v1[1] - centroid[1];

        const double r1z =
            v1[2] - centroid[2];


        const double r2x =
            v2[0] - centroid[0];

        const double r2y =
            v2[1] - centroid[1];

        const double r2z =
            v2[2] - centroid[2];


        const double r0Squared =
            r0x * r0x +
            r0y * r0y +
            r0z * r0z;

        const double r1Squared =
            r1x * r1x +
            r1y * r1y +
            r1z * r1z;

        const double r2Squared =
            r2x * r2x +
            r2y * r2y +
            r2z * r2z;


        const double centroidRadius =
            std::sqrt(
                std::max(
                    r0Squared,
                    std::max(
                        r1Squared,
                        r2Squared)));


        //---------------------------------------------------------------------
        // Store the per-facet geometric information.
        //---------------------------------------------------------------------

        facetGeometry[facetID] =
        {
            centroid,
            area,
            quality,
            centroidRadius
        };


        totalArea +=
            area;
    }


    return
        totalArea /
        static_cast<double>(
            facetCount);
}


void checkFacetMeshQuality(
    const std::vector<FacetGeometry>& facetGeometry,
    const double averageFacetArea)
{
    const double maximumFacetArea =
        MAX_FACET_AREA_RATIO *
        averageFacetArea;


    //-------------------------------------------------------------------------
    // Inspect every facet independently.
    //
    // Low triangle quality and unusually large facet area are mesh-quality
    // diagnostics only. They do not make the STL geometry invalid.
    //
    // Both conditions are accumulated in parallel and reported as summarized
    // warnings after the parallel region.
    //-------------------------------------------------------------------------

    std::size_t lowQualityCount =
        0;

    std::size_t oversizedFacetCount =
        0;


    #pragma omp parallel for \
        reduction(+:lowQualityCount,oversizedFacetCount) \
        schedule(static)

    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                facetGeometry.size());
        ++index)
    {
        const std::size_t facetID =
            static_cast<std::size_t>(
                index);


        const auto& geometry =
            facetGeometry[facetID];


        if(geometry.quality <
           MIN_FACET_QUALITY)
        {
            ++lowQualityCount;
        }


        if(geometry.area >
           maximumFacetArea)
        {
            ++oversizedFacetCount;
        }
    }


    //-------------------------------------------------------------------------
    // Report summarized mesh-quality warnings.
    //-------------------------------------------------------------------------

    if(lowQualityCount > 0)
    {
        std::cerr
            << "Warning: STL geometry contains "
            << lowQualityCount
            << " low-quality facets "
            << "(quality < "
            << MIN_FACET_QUALITY
            << ")."
            << std::endl;
    }


    if(oversizedFacetCount > 0)
    {
        std::cerr
            << "Warning: STL geometry contains "
            << oversizedFacetCount
            << " oversized facets "
            << "(area > "
            << MAX_FACET_AREA_RATIO
            << " times the average facet area)."
            << std::endl;
    }
}

double computeFacetSpatialCellSize(
    const double averageFacetArea)
{
    constexpr double FACET_CELL_SIZE_FACTOR =
        2.0;

    if(!std::isfinite(averageFacetArea) ||
       averageFacetArea <= 0.0)
    {
        throw std::runtime_error(
            "Invalid STL geometry: "
            "invalid average facet area.");
    }


    return
        FACET_CELL_SIZE_FACTOR *
        std::sqrt(averageFacetArea);
}

void computeFacetCellRange(
    const FacetGeometry& facetGeometry,
    const GeometryBounds& bounds,
    const double cellSize,
    const double vertexTolerance,
    BinKey& minimumKey,
    BinKey& maximumKey)
{
    const double radius =
        facetGeometry.centroidRadius +
        vertexTolerance;


    minimumKey =
    {
        static_cast<std::int64_t>(
            std::floor(
                (facetGeometry.centroid[0] -
                 radius -
                 bounds.min[0]) /
                cellSize)),

        static_cast<std::int64_t>(
            std::floor(
                (facetGeometry.centroid[1] -
                 radius -
                 bounds.min[1]) /
                cellSize)),

        static_cast<std::int64_t>(
            std::floor(
                (facetGeometry.centroid[2] -
                 radius -
                 bounds.min[2]) /
                cellSize))
    };


    maximumKey =
    {
        static_cast<std::int64_t>(
            std::floor(
                (facetGeometry.centroid[0] +
                 radius -
                 bounds.min[0]) /
                cellSize)),

        static_cast<std::int64_t>(
            std::floor(
                (facetGeometry.centroid[1] +
                 radius -
                 bounds.min[1]) /
                cellSize)),

        static_cast<std::int64_t>(
            std::floor(
                (facetGeometry.centroid[2] +
                 radius -
                 bounds.min[2]) /
                cellSize))
    };
}

void buildFacetSpatialCells(
    const std::vector<FacetGeometry>& facetGeometry,
    const GeometryBounds& bounds,
    const double averageFacetArea,
    FacetCellEntries& entries)
{
    const double cellSize =
        computeFacetSpatialCellSize(
            averageFacetArea);

    const double vertexTolerance =
        bounds.scale *
        RELATIVE_VERTEX_TOLERANCE;


    //-------------------------------------------------------------------------
    // Generate facet-cell entries independently on each OpenMP thread.
    //
    // Each facet is inserted into every spatial cell overlapped by its
    // centroid-radius search region. Thread-local storage avoids synchronization
    // while entries are being generated.
    //-------------------------------------------------------------------------

    const int maxThreads =
        omp_get_max_threads();

    std::vector<FacetCellEntries>
        threadEntries(
            static_cast<std::size_t>(
                maxThreads));


    #pragma omp parallel
    {
        const int threadID =
            omp_get_thread_num();

        auto& localEntries =
            threadEntries[
                static_cast<std::size_t>(
                    threadID)];


        #pragma omp for schedule(static)
        for (std::int64_t i = 0;
             i < static_cast<std::int64_t>(
                     facetGeometry.size());
             ++i)
        {
            BinKey minimumKey;
            BinKey maximumKey;

            computeFacetCellRange(
                facetGeometry[
                    static_cast<std::size_t>(i)],
                bounds,
                cellSize,
                vertexTolerance,
                minimumKey,
                maximumKey);


            for (std::int64_t z = minimumKey.z;
                 z <= maximumKey.z;
                 ++z)
            {
                for (std::int64_t y = minimumKey.y;
                     y <= maximumKey.y;
                     ++y)
                {
                    for (std::int64_t x = minimumKey.x;
                         x <= maximumKey.x;
                         ++x)
                    {
                        localEntries.push_back(
                            {
                                {x, y, z},
                                static_cast<std::size_t>(i)
                            });
                    }
                }
            }
        }
    }


    //-------------------------------------------------------------------------
    // Concatenate thread-local entries.
    //-------------------------------------------------------------------------

    std::size_t totalEntries = 0;

    for (const auto& localEntries :
         threadEntries)
    {
        totalEntries +=
            localEntries.size();
    }


    entries.clear();

    entries.reserve(
        totalEntries);


    for (auto& localEntries :
         threadEntries)
    {
        entries.insert(
            entries.end(),
            localEntries.begin(),
            localEntries.end());
    }


    //-------------------------------------------------------------------------
    // Sort first by spatial cell and then by facet ID.
    //
    // This produces contiguous ranges for every spatial cell and makes the
    // representation independent of OpenMP execution order.
    //-------------------------------------------------------------------------

    std::sort(
        entries.begin(),
        entries.end(),
        facetCellEntryLess);
}

//=============================================================================
// Temporary facet spatial-cell coverage validation
//=============================================================================
//
// TEMPORARY:
//
// Verify that the spatial-cell range assigned to every facet fully covers
// the welding-tolerance neighborhood of all three of its vertices.
//
// This check is used only while the new facet-based candidate search is being
// validated. Remove it once the new geometric-vertex welding path has fully
// replaced the legacy spatial-bin implementation and passed all regression
// tests.
//

void validateFacetSpatialCellCoverage(
    const STLData& data,
    const std::vector<FacetGeometry>& facetGeometry,
    const GeometryBounds& bounds,
    const double averageFacetArea)
{
    const double cellSize =
        computeFacetSpatialCellSize(
            averageFacetArea);

    const double vertexTolerance =
        bounds.scale *
        RELATIVE_VERTEX_TOLERANCE;


    std::int64_t firstInvalidFacet =
        -1;

    std::int64_t firstInvalidVertex =
        -1;


    #pragma omp parallel for schedule(static)
    for (std::int64_t i = 0;
         i < static_cast<std::int64_t>(
                 data.facets.size());
         ++i)
    {
        BinKey facetMinimumKey;
        BinKey facetMaximumKey;

        computeFacetCellRange(
            facetGeometry[
                static_cast<std::size_t>(i)],
            bounds,
            cellSize,
            vertexTolerance,
            facetMinimumKey,
            facetMaximumKey);


        const auto& facet =
            data.facets[
                static_cast<std::size_t>(i)];


        for (std::size_t v = 0;
             v < 3;
             ++v)
        {
            const auto& vertex =
                facet.vertices[v];


            const BinKey vertexMinimumKey
            {
                static_cast<std::int64_t>(
                    std::floor(
                        (vertex[0] -
                         vertexTolerance -
                         bounds.min[0]) /
                        cellSize)),

                static_cast<std::int64_t>(
                    std::floor(
                        (vertex[1] -
                         vertexTolerance -
                         bounds.min[1]) /
                        cellSize)),

                static_cast<std::int64_t>(
                    std::floor(
                        (vertex[2] -
                         vertexTolerance -
                         bounds.min[2]) /
                        cellSize))
            };


            const BinKey vertexMaximumKey
            {
                static_cast<std::int64_t>(
                    std::floor(
                        (vertex[0] +
                         vertexTolerance -
                         bounds.min[0]) /
                        cellSize)),

                static_cast<std::int64_t>(
                    std::floor(
                        (vertex[1] +
                         vertexTolerance -
                         bounds.min[1]) /
                        cellSize)),

                static_cast<std::int64_t>(
                    std::floor(
                        (vertex[2] +
                         vertexTolerance -
                         bounds.min[2]) /
                        cellSize))
            };


            const bool covered =
                vertexMinimumKey.x >= facetMinimumKey.x &&
                vertexMinimumKey.y >= facetMinimumKey.y &&
                vertexMinimumKey.z >= facetMinimumKey.z &&

                vertexMaximumKey.x <= facetMaximumKey.x &&
                vertexMaximumKey.y <= facetMaximumKey.y &&
                vertexMaximumKey.z <= facetMaximumKey.z;


            if (!covered)
            {
                #pragma omp critical
                {
                    if (firstInvalidFacet < 0 ||
                        i < firstInvalidFacet ||
                        (i == firstInvalidFacet &&
                         static_cast<std::int64_t>(v) <
                             firstInvalidVertex))
                    {
                        firstInvalidFacet =
                            i;

                        firstInvalidVertex =
                            static_cast<std::int64_t>(v);
                    }
                }
            }
        }
    }


    if (firstInvalidFacet >= 0)
    {
        throw std::runtime_error(
            "Internal STL validation error: "
            "facet spatial-cell coverage does not fully cover "
            "the welding-tolerance neighborhood of facet " +
            std::to_string(firstInvalidFacet) +
            ", vertex " +
            std::to_string(firstInvalidVertex) +
            ".");
    }
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

void buildFlatEdges(
    const GeometricVertices& geometry,
    std::vector<FlatEdge>& flatEdges)
{
    const std::size_t facetCount =
        geometry.facetVertexIDs.size();


    //-------------------------------------------------------------------------
    // Allocate exactly three edge records for every triangular facet.
    //-------------------------------------------------------------------------

    flatEdges.resize(
        facetCount * 3);


    //-------------------------------------------------------------------------
    // Generate edges in parallel.
    //
    // Facet i writes only:
    //
    //     flatEdges[3*i + 0]
    //     flatEdges[3*i + 1]
    //     flatEdges[3*i + 2]
    //
    // Therefore no synchronization is required.
    //-------------------------------------------------------------------------

    #pragma omp parallel for schedule(static)

    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                facetCount);
        ++index)
    {
        const std::size_t facetID =
            static_cast<std::size_t>(
                index);


        const auto& vertices =
            geometry.facetVertexIDs[
                facetID];


        const std::array<
            std::pair<std::size_t,std::size_t>,
            3>
            edges =
            {{
                {vertices[0], vertices[1]},
                {vertices[1], vertices[2]},
                {vertices[2], vertices[0]}
            }};


        for(std::size_t e = 0;
            e < 3;
            ++e)
        {
            const std::size_t v0 =
                edges[e].first;

            const std::size_t v1 =
                edges[e].second;


            FlatEdge& flatEdge =
                flatEdges[
                    facetID * 3 + e];


            flatEdge.key =
                makeEdgeKey(
                    v0,
                    v1);


            flatEdge.use = {
                facetID,
                edgeForward(v0,v1)
            };
        }
    }
}

bool flatEdgeLess(
    const FlatEdge& a,
    const FlatEdge& b) noexcept
{
    if(a.key.v0 != b.key.v0)
    {
        return
            a.key.v0 < b.key.v0;
    }

    if(a.key.v1 != b.key.v1)
    {
        return
            a.key.v1 < b.key.v1;
    }

    if(a.use.facetID != b.use.facetID)
    {
        return
            a.use.facetID < b.use.facetID;
    }

    return
        static_cast<int>(a.use.forward) <
        static_cast<int>(b.use.forward);
}

bool sameEdgeKey(
    const EdgeKey& a,
    const EdgeKey& b) noexcept
{
    return
        a.v0 == b.v0 &&
        a.v1 == b.v1;
}

void sortFlatEdges(
    std::vector<FlatEdge>& flatEdges)
{
    std::sort(
        flatEdges.begin(),
        flatEdges.end(),
        flatEdgeLess);
}

void validateFlatEdges(
    const std::vector<FlatEdge>& flatEdges)
{
    std::size_t begin = 0;


    while(begin < flatEdges.size())
    {
        //---------------------------------------------------------------------
        // Find the complete group corresponding to one geometric edge.
        //----------------------------------------------------------------------

        std::size_t end =
            begin + 1;


        while(end < flatEdges.size() &&
              sameEdgeKey(
                  flatEdges[begin].key,
                  flatEdges[end].key))
        {
            ++end;
        }


        const std::size_t useCount =
            end - begin;


        //---------------------------------------------------------------------
        // A closed surface requires every geometric edge to be shared by
        // exactly two facets.
        //----------------------------------------------------------------------

        if(useCount == 1)
        {
            throw std::runtime_error(
                "Invalid STL geometry: "
                "open boundary edge detected.");
        }


        //---------------------------------------------------------------------
        // More than two incident facets form a non-manifold edge.
        //----------------------------------------------------------------------

        if(useCount > 2)
        {
            throw std::runtime_error(
                "Invalid STL geometry: "
                "non-manifold edge detected.");
        }


        //---------------------------------------------------------------------
        // For a consistently wound orientable surface, the two incident
        // facets must traverse their shared edge in opposite directions.
        //----------------------------------------------------------------------

        const auto& use0 =
            flatEdges[begin].use;

        const auto& use1 =
            flatEdges[begin + 1].use;


        if(use0.forward ==
           use1.forward)
        {
            throw std::runtime_error(
                "Invalid STL geometry: "
                "inconsistent facet winding detected.");
        }


        begin =
            end;
    }
}

void buildFacetAdjacency(
    const std::vector<FlatEdge>& flatEdges,
    std::vector<std::array<std::size_t,3>>& adjacency)
{
    //-------------------------------------------------------------------------
    // Count how many neighbors have been assigned to each facet.
    //
    // validateFlatEdges() has already guaranteed that every geometric edge
    // has exactly two incident facets.
    //-------------------------------------------------------------------------

    std::vector<unsigned char>
        neighborCount(
            adjacency.size(),
            0);


    std::size_t begin = 0;


    while(begin < flatEdges.size())
    {
        const std::size_t end =
            begin + 2;


        const std::size_t facetA =
            flatEdges[begin].use.facetID;

        const std::size_t facetB =
            flatEdges[begin + 1].use.facetID;


        const unsigned char slotA =
            neighborCount[facetA]++;

        const unsigned char slotB =
            neighborCount[facetB]++;


        if(slotA >= 3 ||
           slotB >= 3)
        {
            throw std::runtime_error(
                "Internal STL topology error: "
                "facet has more than three adjacent facets.");
        }


        adjacency[facetA][slotA] =
            facetB;

        adjacency[facetB][slotB] =
            facetA;


        begin =
            end;
    }


    //-------------------------------------------------------------------------
    // Every triangular facet in a valid closed manifold surface must have
    // exactly three neighboring facets.
    //-------------------------------------------------------------------------

    for(std::size_t i = 0;
        i < neighborCount.size();
        ++i)
    {
        if(neighborCount[i] != 3)
        {
            throw std::runtime_error(
                "Internal STL topology error: "
                "facet does not have exactly three adjacent facets.");
        }
    }
}

void buildComponents(
    const std::vector<std::array<std::size_t,3>>& adjacency,
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
        const auto firstFacetID =
            component.facets[0];

        const auto firstVertexID =
            geometry.facetVertexIDs[
                firstFacetID][0];

        const auto& firstPoint =
            geometry.vertices[
                firstVertexID];


        double minX = firstPoint[0];
        double minY = firstPoint[1];
        double minZ = firstPoint[2];

        double maxX = firstPoint[0];
        double maxY = firstPoint[1];
        double maxZ = firstPoint[2];


        //---------------------------------------------------------------------
        // Compute component bounds in parallel over its facets.
        //---------------------------------------------------------------------

        #pragma omp parallel for \
            reduction(min:minX,minY,minZ) \
            reduction(max:maxX,maxY,maxZ) \
            schedule(static)

        for(std::ptrdiff_t index = 0;
            index <
                static_cast<std::ptrdiff_t>(
                    component.facets.size());
            ++index)
        {
            const std::size_t facetID =
                component.facets[
                    static_cast<std::size_t>(
                        index)];


            const auto& vertexIDs =
                geometry.facetVertexIDs[
                    facetID];


            for(std::size_t v = 0;
                v < 3;
                ++v)
            {
                const auto& p =
                    geometry.vertices[
                        vertexIDs[v]];


                minX =
                    std::min(
                        minX,
                        p[0]);

                minY =
                    std::min(
                        minY,
                        p[1]);

                minZ =
                    std::min(
                        minZ,
                        p[2]);


                maxX =
                    std::max(
                        maxX,
                        p[0]);

                maxY =
                    std::max(
                        maxY,
                        p[1]);

                maxZ =
                    std::max(
                        maxZ,
                        p[2]);
            }
        }


        GeometryBounds bounds;

        bounds.min = {
            minX,
            minY,
            minZ
        };

        bounds.max = {
            maxX,
            maxY,
            maxZ
        };


        const double dx =
            maxX - minX;

        const double dy =
            maxY - minY;

        const double dz =
            maxZ - minZ;


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


    long double volume =
        0.0L;


    //---------------------------------------------------------------------
    // Accumulate signed tetrahedral volumes in parallel.
    //---------------------------------------------------------------------

    #pragma omp parallel for \
        reduction(+:volume) \
        schedule(static)

    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                component.facets.size());
        ++index)
    {
        const std::size_t facetID =
            component.facets[
                static_cast<std::size_t>(
                    index)];


        const auto& ids =
            geometry.facetVertexIDs[
                facetID];


        const auto& v0 =
            geometry.vertices[
                ids[0]];

        const auto& v1 =
            geometry.vertices[
                ids[1]];

        const auto& v2 =
            geometry.vertices[
                ids[2]];


        const long double x0 =
            static_cast<long double>(
                v0[0] - referenceX);

        const long double y0 =
            static_cast<long double>(
                v0[1] - referenceY);

        const long double z0 =
            static_cast<long double>(
                v0[2] - referenceZ);


        const long double x1 =
            static_cast<long double>(
                v1[0] - referenceX);

        const long double y1 =
            static_cast<long double>(
                v1[1] - referenceY);

        const long double z1 =
            static_cast<long double>(
                v1[2] - referenceZ);


        const long double x2 =
            static_cast<long double>(
                v2[0] - referenceX);

        const long double y2 =
            static_cast<long double>(
                v2[1] - referenceY);

        const long double z2 =
            static_cast<long double>(
                v2[2] - referenceZ);


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
// Build and validate the surface topology using a flat canonical edge
// representation.
//
// Each triangular facet contributes exactly three edge records. After sorting
// by canonical EdgeKey, all uses of the same geometric edge are contiguous.
//
// This stage will:
//
//   1. Build the flat canonical edge representation.
//
//   2. Sort geometric edge uses by EdgeKey.
//
//   3. Validate surface closedness.
//      Every geometric edge must be shared by exactly two facets.
//
//   4. Detect non-manifold topology.
//      An edge shared by more than two facets is invalid.
//
//   5. Validate local winding consistency.
//      Two facets sharing an edge must traverse that edge in opposite
//      directions.
//
//   6. Build fixed-size three-neighbor facet adjacency.
//
//   7. Identify connected surface components.
//
// Multiple connected components may be identified here. Their geometric
// admissibility is determined later by the geometry module.
//
// The welded geometric vertices, facet adjacency, and connected components
// are retained in FacetTopology for later geometry processing.
//
void validateTopologyWindingAndComponents(
    FacetTopology& topology)
{
    //-------------------------------------------------------------------------
    // Build the flat edge representation.
    //-------------------------------------------------------------------------

    std::vector<FlatEdge>
        flatEdges;

    buildFlatEdges(
        topology.geometry,
        flatEdges);


    //-------------------------------------------------------------------------
    // Sort by canonical geometric edge.
    //-------------------------------------------------------------------------

    sortFlatEdges(
        flatEdges);


    //-------------------------------------------------------------------------
    // Validate surface closedness, manifoldness, and local facet winding.
    //-------------------------------------------------------------------------

    validateFlatEdges(
        flatEdges);


    //-------------------------------------------------------------------------
    // Build fixed-size facet adjacency.
    //
    // Every valid triangular facet has exactly three neighboring facets.
    //-------------------------------------------------------------------------

    topology.adjacency.clear();

    topology.adjacency.resize(
        topology.geometry.facetVertexIDs.size());


    buildFacetAdjacency(
        flatEdges,
        topology.adjacency);


    //-------------------------------------------------------------------------
    // Find connected surface components.
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
    // Compute and validate per-facet geometric information.
    //
    // Degenerate facets have already been rejected by Step 1. The resulting
    // centroid and characteristic facet dimensions are retained for the
    // centroid-based spatial search used during geometric vertex welding.
    //-------------------------------------------------------------------------

    const double averageFacetArea =
        computeFacetGeometry(
            data,
            topology.facetGeometry);


    checkFacetMeshQuality(
        topology.facetGeometry,
        averageFacetArea);

    //-------------------------------------------------------------------------
    // Build the facet spatial search structure.
    //
    // This structure is not yet used by geometric vertex welding. It is built
    // here first so that the new parallel candidate-search path can be validated
    // independently before replacing the existing welding implementation.
    //-------------------------------------------------------------------------

    FacetCellEntries facetCellEntries;

    buildFacetSpatialCells(
        topology.facetGeometry,
        bounds,
        averageFacetArea,
        facetCellEntries);

//-------------------------------------------------------------------------
// TEMPORARY:
//
// Validate the geometric coverage used by the new facet spatial search.
// Remove this check after the new welding path has fully replaced the
// legacy implementation and passed all regression tests.
//-------------------------------------------------------------------------

validateFacetSpatialCellCoverage(
    data,
    topology.facetGeometry,
    bounds,
    averageFacetArea);

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
