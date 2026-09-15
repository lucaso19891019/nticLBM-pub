#include "stl_validator.hpp"
#include "stl_topology.hpp"
#include "parallel.hpp"

#include <omp.h>

#include <atomic>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>

#include <array>
#include <cstdint>
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

struct FacetPair
{
    std::size_t first;
    std::size_t second;

    bool operator==(
        const FacetPair& other) const noexcept
    {
        return
            first == other.first &&
            second == other.second;
    }
};


using FacetPairs =
    std::vector<FacetPair>;

//=============================================================================
// Vertex pair
//=============================================================================
//
// A raw STL vertex is identified by its flattened index:
//
//     rawVertexID = facetID * 3 + localVertexID
//
// VertexPair stores two raw STL vertices that are within the geometric
// welding tolerance. The smaller raw vertex ID is always stored first.
//

struct VertexPair
{
    std::size_t first;
    std::size_t second;

    bool operator==(
        const VertexPair& other) const noexcept
    {
        return
            first == other.first &&
            second == other.second;
    }
};


using VertexPairs =
    std::vector<VertexPair>;

using RawVertexRepresentatives =
    std::vector<std::size_t>;

bool vertexPairLess(
    const VertexPair& a,
    const VertexPair& b) noexcept
{
    if (a.first != b.first) {
        return a.first < b.first;
    }

    return a.second < b.second;
}

struct FacetCellRange
{
    std::size_t begin;
    std::size_t end;
};


using FacetCellRanges =
    std::vector<FacetCellRange>;

bool facetPairLess(
    const FacetPair& a,
    const FacetPair& b) noexcept
{
    if (a.first != b.first) {
        return a.first < b.first;
    }

    return a.second < b.second;
}

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

    std::vector<std::size_t> threadSizes(
        threadEntries.size());

    for (std::size_t i = 0;
        i < threadEntries.size();
        ++i)
    {
        threadSizes[i] =
            threadEntries[i].size();
    }


    std::vector<std::size_t> offsets;


    const std::size_t totalEntries =
        common::parallelScan(
            threadSizes,
            offsets);


    entries.resize(
        totalEntries);


    #pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0;
        i < static_cast<std::ptrdiff_t>(
                threadEntries.size());
        ++i)
    {
        const std::size_t threadID =
            static_cast<std::size_t>(i);


        std::copy(
            threadEntries[threadID].begin(),
            threadEntries[threadID].end(),
            entries.begin() +
                offsets[threadID]);
    }

    //-------------------------------------------------------------------------
    // Sort first by spatial cell and then by facet ID.
    //
    // This produces contiguous ranges for every spatial cell and makes the
    // representation independent of OpenMP execution order.
    //-------------------------------------------------------------------------

    common::parallelSort(
        entries,
        facetCellEntryLess);
}

void buildFacetCellRanges(
    const FacetCellEntries& entries,
    FacetCellRanges& ranges)
{
    ranges.clear();

    if (entries.empty()) {
        return;
    }


    std::size_t begin = 0;

    while (begin < entries.size())
    {
        std::size_t end =
            begin + 1;

        while (end < entries.size() &&
               entries[end].key ==
                   entries[begin].key)
        {
            ++end;
        }


        ranges.push_back(
            {
                begin,
                end
            });


        begin = end;
    }
}

void buildCandidateFacetPairs(
    const FacetCellEntries& entries,
    const FacetCellRanges& ranges,
    FacetPairs& pairs)
{
    const int maxThreads =
        omp_get_max_threads();

    std::vector<FacetPairs>
        threadPairs(
            static_cast<std::size_t>(
                maxThreads));


    //-------------------------------------------------------------------------
    // Generate candidate facet pairs independently for each spatial cell.
    //
    // A pair may be generated by multiple cells because facet search regions
    // overlap. Duplicate pairs are removed after all thread-local results have
    // been merged.
    //-------------------------------------------------------------------------

    #pragma omp parallel
    {
        const int threadID =
            omp_get_thread_num();

        auto& localPairs =
            threadPairs[
                static_cast<std::size_t>(
                    threadID)];


        #pragma omp for schedule(dynamic, 16)
        for (std::int64_t r = 0;
             r < static_cast<std::int64_t>(
                     ranges.size());
             ++r)
        {
            const auto& range =
                ranges[
                    static_cast<std::size_t>(r)];

            const std::size_t count =
                range.end -
                range.begin;


            if (count < 2) {
                continue;
            }


            for (std::size_t a = range.begin;
                 a + 1 < range.end;
                 ++a)
            {
                const std::size_t facetA =
                    entries[a].facetID;


                for (std::size_t b = a + 1;
                     b < range.end;
                     ++b)
                {
                    const std::size_t facetB =
                        entries[b].facetID;


                    if (facetA == facetB) {
                        continue;
                    }


                    if (facetA < facetB)
                    {
                        localPairs.push_back(
                            {
                                facetA,
                                facetB
                            });
                    }
                    else
                    {
                        localPairs.push_back(
                            {
                                facetB,
                                facetA
                            });
                    }
                }
            }
        }
    }


    //-------------------------------------------------------------------------
    // Merge thread-local candidate pairs.
    //-------------------------------------------------------------------------

    std::vector<std::size_t> threadSizes(
        threadPairs.size());


    for (std::size_t i = 0;
        i < threadPairs.size();
        ++i)
    {
        threadSizes[i] =
            threadPairs[i].size();
    }


    std::vector<std::size_t> offsets;


    const std::size_t totalPairs =
        common::parallelScan(
            threadSizes,
            offsets);


    pairs.resize(
        totalPairs);


    #pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0;
        i < static_cast<std::ptrdiff_t>(
                threadPairs.size());
        ++i)
    {
        const std::size_t threadID =
            static_cast<std::size_t>(i);


        std::copy(
            threadPairs[threadID].begin(),
            threadPairs[threadID].end(),
            pairs.begin() +
                offsets[threadID]);
    }


    //-------------------------------------------------------------------------
    // Sort and remove duplicate facet pairs.
    //-------------------------------------------------------------------------

    common::parallelSort(
        pairs,
        facetPairLess);

    pairs.erase(
        std::unique(
            pairs.begin(),
            pairs.end()),
        pairs.end());
}

bool sameGeometricVertex(
    const STLVector& a,
    const STLVector& b,
    const double toleranceSquared);

void buildCandidateVertexPairs(
    const STLData& data,
    const GeometryBounds& bounds,
    const FacetPairs& facetPairs,
    VertexPairs& vertexPairs)
{
    const double vertexTolerance =
        bounds.scale *
        RELATIVE_VERTEX_TOLERANCE;

    const double vertexToleranceSquared =
        vertexTolerance *
        vertexTolerance;


    //-------------------------------------------------------------------------
    // Generate matching raw-vertex pairs independently for every candidate
    // facet pair.
    //
    // Each candidate facet pair contributes at most 3 x 3 possible raw-vertex
    // comparisons. Only vertex pairs within the geometric welding tolerance
    // are retained.
    //
    // Thread-local storage avoids synchronization during pair generation.
    //-------------------------------------------------------------------------

    const int maxThreads =
        omp_get_max_threads();

    std::vector<VertexPairs>
        threadPairs(
            static_cast<std::size_t>(
                maxThreads));


    #pragma omp parallel
    {
        const int threadID =
            omp_get_thread_num();

        auto& localPairs =
            threadPairs[
                static_cast<std::size_t>(
                    threadID)];


        #pragma omp for schedule(static)
        for (std::int64_t pairIndex = 0;
             pairIndex <
                 static_cast<std::int64_t>(
                     facetPairs.size());
             ++pairIndex)
        {
            const FacetPair& facetPair =
                facetPairs[
                    static_cast<std::size_t>(
                        pairIndex)];

            const auto& facetA =
                data.facets[
                    facetPair.first];

            const auto& facetB =
                data.facets[
                    facetPair.second];


            for (std::size_t a = 0;
                 a < 3;
                 ++a)
            {
                const std::size_t rawVertexA =
                    facetPair.first * 3 + a;


                for (std::size_t b = 0;
                     b < 3;
                     ++b)
                {
                    if (!sameGeometricVertex(
                            facetA.vertices[a],
                            facetB.vertices[b],
                            vertexToleranceSquared))
                    {
                        continue;
                    }


                    const std::size_t rawVertexB =
                        facetPair.second * 3 + b;


                    if (rawVertexA < rawVertexB)
                    {
                        localPairs.push_back(
                            {
                                rawVertexA,
                                rawVertexB
                            });
                    }
                    else
                    {
                        localPairs.push_back(
                            {
                                rawVertexB,
                                rawVertexA
                            });
                    }
                }
            }
        }
    }


    //-------------------------------------------------------------------------
    // Merge thread-local vertex pairs.
    //-------------------------------------------------------------------------

    std::vector<std::size_t> threadSizes(
        threadPairs.size());


    for (std::size_t i = 0;
        i < threadPairs.size();
        ++i)
    {
        threadSizes[i] =
            threadPairs[i].size();
    }


    std::vector<std::size_t> offsets;


    const std::size_t totalPairs =
        common::parallelScan(
            threadSizes,
            offsets);


    vertexPairs.resize(
        totalPairs);


    #pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0;
        i < static_cast<std::ptrdiff_t>(
                threadPairs.size());
        ++i)
    {
        const std::size_t threadID =
            static_cast<std::size_t>(i);


        std::copy(
            threadPairs[threadID].begin(),
            threadPairs[threadID].end(),
            vertexPairs.begin() +
                offsets[threadID]);
    }

    //-------------------------------------------------------------------------
    // Sort and remove duplicate raw-vertex pairs.
    //-------------------------------------------------------------------------

    common::parallelSort(
        vertexPairs,
        vertexPairLess);

    vertexPairs.erase(
        std::unique(
            vertexPairs.begin(),
            vertexPairs.end()),
        vertexPairs.end());
}

void buildRawVertexRepresentatives(
    const STLData& data,
    const VertexPairs& vertexPairs,
    RawVertexRepresentatives& representatives)
{
    const std::size_t rawVertexCount =
        data.facets.size() * 3;


    //-------------------------------------------------------------------------
    // Initialize the disjoint-set forest.
    //
    // Shiloach-Vishkin style connectivity algorithm:
    // each vertex initially points to itself.
    //-------------------------------------------------------------------------

    std::vector<std::atomic<std::size_t>> parent(
        rawVertexCount);


#pragma omp parallel for
    for (std::size_t i = 0;
         i < rawVertexCount;
         ++i)
    {
        parent[i].store(
            i,
            std::memory_order_relaxed);
    }


    //-------------------------------------------------------------------------
    // Parallel Shiloach-Vishkin hooking + pointer jumping.
    //
    // Each vertex pair represents an undirected edge.
    // The algorithm repeatedly:
    //
    // 1. Hooks a larger root to a smaller root.
    // 2. Performs pointer jumping compression.
    //
    // The iteration stops when no parent changes.
    //-------------------------------------------------------------------------

    bool changed = true;


    while (changed)
    {
        changed = false;


        //-------------------------------------------------------------------------
        // Hooking phase.
        //
        // For every edge (a,b):
        //
        //     rootA = component root of a
        //     rootB = component root of b
        //
        // Connect the larger root to the smaller root.
        //
        // CAS guarantees correctness when multiple edges update the
        // same component simultaneously.
        //-------------------------------------------------------------------------

#pragma omp parallel for reduction(|| : changed)
        for (std::size_t i = 0;
             i < vertexPairs.size();
             ++i)
        {
            const VertexPair& pair =
                vertexPairs[i];


            std::size_t rootA =
                pair.first;


            while (true)
            {
                const std::size_t next =
                    parent[rootA].load(
                        std::memory_order_relaxed);

                if (next == rootA)
                {
                    break;
                }

                rootA = next;
            }


            std::size_t rootB =
                pair.second;


            while (true)
            {
                const std::size_t next =
                    parent[rootB].load(
                        std::memory_order_relaxed);

                if (next == rootB)
                {
                    break;
                }

                rootB = next;
            }


            if (rootA == rootB)
            {
                continue;
            }


            if (rootA < rootB)
            {
                std::swap(
                    rootA,
                    rootB);
            }


            std::size_t expected =
                rootA;


            if (parent[rootA].compare_exchange_strong(
                    expected,
                    rootB,
                    std::memory_order_relaxed))
            {
                changed = true;
            }
        }


        //-------------------------------------------------------------------------
        // Pointer jumping phase.
        //
        // parent[i] = parent[parent[i]]
        //
        // This flattens the forest and accelerates convergence.
        //-------------------------------------------------------------------------

#pragma omp parallel for reduction(|| : changed)
        for (std::size_t i = 0;
             i < rawVertexCount;
             ++i)
        {
            const std::size_t parentIndex =
                parent[i].load(
                    std::memory_order_relaxed);


            const std::size_t grandParent =
                parent[parentIndex].load(
                    std::memory_order_relaxed);


            if (parentIndex != grandParent)
            {
                parent[i].store(
                    grandParent,
                    std::memory_order_relaxed);

                changed = true;
            }
        }
    }


    //-------------------------------------------------------------------------
    // Store final representatives.
    //
    // After convergence every vertex points to its component root.
    //-------------------------------------------------------------------------

    representatives.resize(
        rawVertexCount);


#pragma omp parallel for
    for (std::size_t i = 0;
         i < rawVertexCount;
         ++i)
    {
        std::size_t root =
            i;


        while (true)
        {
            const std::size_t next =
                parent[root].load(
                    std::memory_order_relaxed);

            if (next == root)
            {
                break;
            }

            root = next;
        }


        representatives[i] =
            root;
    }
}

void buildGeometricVerticesFromRepresentatives(
    const STLData& data,
    const RawVertexRepresentatives& representatives,
    GeometricVertices& geometry)
{
    const std::size_t rawVertexCount =
        representatives.size();


    geometry.facetVertexIDs.resize(
        data.facets.size());


    //-------------------------------------------------------------------------
    // Mark the representative raw vertices.
    //
    // Every Union-Find root corresponds to exactly one geometric vertex.
    //-------------------------------------------------------------------------

    std::vector<std::size_t> representativeFlags(
        rawVertexCount);


    #pragma omp parallel for schedule(static)

    for (std::ptrdiff_t index = 0;
         index <
             static_cast<std::ptrdiff_t>(
                 rawVertexCount);
         ++index)
    {
        const std::size_t rawID =
            static_cast<std::size_t>(
                index);


        representativeFlags[rawID] =
            representatives[rawID] == rawID
                ? std::size_t{1}
                : std::size_t{0};
    }


    //-------------------------------------------------------------------------
    // Assign a compact geometric vertex ID to every representative.
    //-------------------------------------------------------------------------

    std::vector<std::size_t> representativeIDs;


    const std::size_t geometricVertexCount =
        common::parallelScan(
            representativeFlags,
            representativeIDs);


    geometry.vertices.resize(
        geometricVertexCount);


    //-------------------------------------------------------------------------
    // Store every unique geometric vertex at its compact geometric vertex ID.
    //-------------------------------------------------------------------------

    #pragma omp parallel for schedule(static)

    for (std::ptrdiff_t index = 0;
         index <
             static_cast<std::ptrdiff_t>(
                 rawVertexCount);
         ++index)
    {
        const std::size_t rawID =
            static_cast<std::size_t>(
                index);


        if (representativeFlags[rawID] == 0)
        {
            continue;
        }


        const std::size_t facetID =
            rawID / 3;

        const std::size_t localVertex =
            rawID % 3;


        geometry.vertices[
            representativeIDs[rawID]] =
                data.facets[facetID]
                    .vertices[localVertex];
    }


    //-------------------------------------------------------------------------
    // Convert every raw STL vertex to its geometric vertex ID.
    //-------------------------------------------------------------------------

    #pragma omp parallel for schedule(static)

    for (std::ptrdiff_t index = 0;
         index <
             static_cast<std::ptrdiff_t>(
                 data.facets.size());
         ++index)
    {
        const std::size_t facetID =
            static_cast<std::size_t>(
                index);

        auto& vertexIDs =
            geometry.facetVertexIDs[
                facetID];


        for (std::size_t v = 0;
             v < 3;
             ++v)
        {
            const std::size_t rawID =
                facetID * 3 + v;

            const std::size_t representative =
                representatives[rawID];


            vertexIDs[v] =
                representativeIDs[
                    representative];
        }
    }
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

GeometricVertices buildGeometricVertices(
    const STLData& data,
    const GeometryBounds& bounds)
{
    GeometricVertices geometry;


    std::vector<FacetGeometry> facetGeometry;

    const double averageFacetArea =
        computeFacetGeometry(
            data,
            facetGeometry);


    FacetCellEntries facetEntries;

    buildFacetSpatialCells(
        facetGeometry,
        bounds,
        averageFacetArea,
        facetEntries);


    FacetCellRanges facetRanges;

    buildFacetCellRanges(
        facetEntries,
        facetRanges);


    FacetPairs facetPairs;

    buildCandidateFacetPairs(
        facetEntries,
        facetRanges,
        facetPairs);


    VertexPairs vertexPairs;

    buildCandidateVertexPairs(
        data,
        bounds,
        facetPairs,
        vertexPairs);


    RawVertexRepresentatives representatives;

    buildRawVertexRepresentatives(
        data,
        vertexPairs,
        representatives);


    buildGeometricVerticesFromRepresentatives(
        data,
        representatives,
        geometry);


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
    common::parallelSort(
        flatEdges,
        flatEdgeLess);
}

void validateFlatEdges(
    const std::vector<FlatEdge>& flatEdges)
{
    std::vector<std::pair<std::size_t,
                          std::size_t>> groups;


    //---------------------------------------------------------------------
    // Find all geometric edge groups.
    //
    // flatEdges is already sorted by edge key.
    // Therefore all uses of one geometric edge are consecutive.
    //---------------------------------------------------------------------

    const std::size_t edgeCount =
        flatEdges.size();


    if(edgeCount == 0)
    {
        return;
    }


    //-------------------------------------------------------------------------
    // Detect the beginning of every geometric edge group.
    //
    // flatEdges is already sorted by edge key.
    // A value of 1 means that a new group starts at this position.
    //-------------------------------------------------------------------------

    std::vector<std::size_t> boundaries(
        edgeCount);


    #pragma omp parallel for
    for(std::size_t i = 0;
        i < edgeCount;
        ++i)
    {
        if(i == 0)
        {
            boundaries[i] = 1;
        }
        else
        {
            boundaries[i] =
                sameEdgeKey(
                    flatEdges[i - 1].key,
                    flatEdges[i].key)
                ?
                0
                :
                1;
        }
    }


    //-------------------------------------------------------------------------
    // Convert group boundaries into group indices.
    //
    // parallelScan is exclusive scan.
    //-------------------------------------------------------------------------

    std::vector<std::size_t> groupIds;


    common::parallelScan(
        boundaries,
        groupIds);



    const std::size_t groupCount =
        groupIds.back()
        +
        boundaries.back();



    std::vector<std::pair<std::size_t,
                        std::size_t>> groups(
        groupCount);



    #pragma omp parallel for
    for(std::size_t i = 0;
        i < edgeCount;
        ++i)
    {
        if(boundaries[i])
        {
            const std::size_t group =
                groupIds[i];


            groups[group].first =
                i;
        }
    }



    #pragma omp parallel for
    for(std::size_t i = 0;
        i < edgeCount - 1;
        ++i)
    {
        if(boundaries[i + 1])
        {
            const std::size_t group =
                groupIds[i];


            groups[group].second =
                i + 1;
        }
    }


    // Last group ends at edgeCount.

    groups[groupCount - 1].second =
        edgeCount;


    for(const auto& group : groups)
    {
        if(group.first >= group.second ||
        group.second > edgeCount)
        {
            throw std::runtime_error(
                "Invalid edge group construction.");
        }
    }


    //---------------------------------------------------------------------
    // Validate edge groups in parallel.
    //---------------------------------------------------------------------

    std::atomic<int> errorCode(
        0);


#pragma omp parallel for
    for (std::size_t group = 0;
         group < groups.size();
         ++group)
    {
        if (errorCode.load(
                std::memory_order_relaxed) != 0)
        {
            continue;
        }


        const std::size_t begin =
            groups[group].first;


        const std::size_t end =
            groups[group].second;


        const std::size_t useCount =
            end - begin;



        //-----------------------------------------------------------------
        // A closed surface requires every geometric edge to be shared by
        // exactly two facets.
        //-----------------------------------------------------------------

        if (useCount == 1)
        {
            errorCode.store(
                1,
                std::memory_order_relaxed);

            continue;
        }



        //-----------------------------------------------------------------
        // More than two incident facets form a non-manifold edge.
        //-----------------------------------------------------------------

        if (useCount > 2)
        {
            errorCode.store(
                2,
                std::memory_order_relaxed);

            continue;
        }



        //-----------------------------------------------------------------
        // Two incident facets must traverse their shared edge in opposite
        // directions.
        //-----------------------------------------------------------------

        const auto& use0 =
            flatEdges[begin].use;


        const auto& use1 =
            flatEdges[begin + 1].use;


        if (use0.forward ==
            use1.forward)
        {
            errorCode.store(
                3,
                std::memory_order_relaxed);
        }
    }



    //---------------------------------------------------------------------
    // Report validation result.
    //---------------------------------------------------------------------

    switch (errorCode.load())
    {
        case 1:
            throw std::runtime_error(
                "Invalid STL geometry: "
                "open boundary edge detected.");

        case 2:
            throw std::runtime_error(
                "Invalid STL geometry: "
                "non-manifold edge detected.");

        case 3:
            throw std::runtime_error(
                "Invalid STL geometry: "
                "inconsistent facet winding detected.");

        default:
            break;
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
    // Build candidate neighboring facet pairs from the spatial-cell structure.
    //-------------------------------------------------------------------------

    FacetCellRanges facetCellRanges;

    buildFacetCellRanges(
        facetCellEntries,
        facetCellRanges);


    FacetPairs candidateFacetPairs;

    buildCandidateFacetPairs(
        facetCellEntries,
        facetCellRanges,
        candidateFacetPairs);

    //-------------------------------------------------------------------------
    // Build matching raw-vertex pairs from the candidate facet pairs.
    //
    // This is the new parallel welding candidate path. It is currently built
    // only for validation and does not yet replace buildGeometricVertices().
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
