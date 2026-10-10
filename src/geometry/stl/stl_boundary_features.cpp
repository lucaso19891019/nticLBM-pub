
#include "stl_geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <omp.h>

namespace ntic::lbm::geometry
{
namespace
{

using FacePair = std::array<std::size_t, 2>;

double dot(const Point& a, const Point& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

Point subtract(const Point& a, const Point& b)
{
    return {
        a[0]-b[0],
        a[1]-b[1],
        a[2]-b[2]
    };
}

double norm(const Point& a)
{
    return std::sqrt(dot(a, a));
}

struct Plane
{
    Point normal{};
    double offset = 0.0;
};

struct FeatureSegment
{
    std::size_t v0 = 0;
    std::size_t v1 = 0;
    FacePair faces{};
};

// Check against a fixed reference plane.
// This prevents cumulative merging along curved surfaces.
bool facetOnPlane(
    const stl::FacetTopology& topology,
    std::size_t facetID,
    const Plane& plane,
    double distanceTolerance,
    double normalCosine)
{
    const auto& facet = topology.facetGeometry[facetID];

    if(dot(facet.normal, plane.normal) < normalCosine)
        return false;

    const auto& ids =
        topology.geometry.facetVertexIDs[facetID];

    for(std::size_t vertexID : ids)
    {
        const Point& p =
            topology.geometry.vertices[vertexID];

        if(std::abs(dot(plane.normal, p) -
                    plane.offset) > distanceTolerance)
            return false;
    }

    return true;
}

FacePair makeFacePair(std::size_t a, std::size_t b)
{
    if(a > b)
        std::swap(a, b);

    return {a, b};
}

// Number of feature segments incident to vertexID
// with the same pair of adjacent faces.
std::size_t matchingDegree(
    std::size_t vertexID,
    const FacePair& faces,
    const std::vector<FeatureSegment>& segments,
    const std::vector<std::vector<std::size_t>>& incident)
{
    std::size_t degree = 0;

    for(std::size_t segmentID : incident[vertexID])
    {
        if(segments[segmentID].faces == faces)
            ++degree;
    }

    return degree;
}

// Follow a feature curve through vertices of degree 2.
// At endpoints and junctions, the curve terminates.
void traceFeatureCurve(
    std::size_t firstSegment,
    std::size_t startVertex,
    std::size_t curveID,
    const std::vector<FeatureSegment>& segments,
    const std::vector<std::vector<std::size_t>>& incident,
    std::vector<std::size_t>& segmentCurveIDs)
{
    const std::size_t unassigned =
        std::numeric_limits<std::size_t>::max();

    const FacePair faces = segments[firstSegment].faces;

    std::size_t currentSegment = firstSegment;
    std::size_t currentVertex = startVertex;

    while(segmentCurveIDs[currentSegment] == unassigned)
    {
        segmentCurveIDs[currentSegment] = curveID;

        const FeatureSegment& segment =
            segments[currentSegment];

        const std::size_t nextVertex =
            segment.v0 == currentVertex
                ? segment.v1
                : segment.v0;

        if(matchingDegree(
               nextVertex, faces, segments, incident) != 2)
            break;

        std::size_t nextSegment = unassigned;

        for(std::size_t candidate : incident[nextVertex])
        {
            if(candidate == currentSegment)
                continue;

            if(segments[candidate].faces == faces)
            {
                nextSegment = candidate;
                break;
            }
        }

        if(nextSegment == unassigned ||
           segmentCurveIDs[nextSegment] != unassigned)
            break;

        currentVertex = nextVertex;
        currentSegment = nextSegment;
    }
}

} // namespace

void STLGeometry::identifyBoundaryFeatures()
{
    STLBoundaryFeatures result;

    const auto& mesh = topology.geometry;
    const std::size_t facetCount =
        mesh.facetVertexIDs.size();

    const std::size_t vertexCount =
        mesh.vertices.size();

    result.facetFaceIDs.assign(facetCount, 0);

    if(flowType != FlowType::Internal)
    {
        boundaryFeatures = std::move(result);
        return;
    }

    if(topology.facetGeometry.size() != facetCount ||
       topology.adjacency.size() != facetCount)
    {
        throw std::runtime_error(
            "Incomplete STL facet geometry or adjacency.");
    }

    // ---------------------------------------------------------
    // 1. Planarity tolerances
    // ---------------------------------------------------------

    const double scale = norm(subtract(
        bounds.max, bounds.min));

    if(!(scale > 0.0) || !std::isfinite(scale))
    {
        throw std::runtime_error(
            "Invalid STL geometric scale.");
    }

    // These parameters can later be made configurable.
    constexpr double relativePlaneTolerance = 1.0e-7;
    constexpr double normalAngleDegrees = 1.0;
    constexpr double relativeMinimumArea = 1.0e-8;

    constexpr double pi =
        3.14159265358979323846;

    const double distanceTolerance =
        relativePlaneTolerance * scale;

    const double normalCosine =
        std::cos(normalAngleDegrees * pi / 180.0);

    const double minimumArea =
        relativeMinimumArea * scale * scale;

    // ---------------------------------------------------------
    // 2. Mark individually valid facets in parallel
    // ---------------------------------------------------------

    std::vector<unsigned char> valid(facetCount, 0);

    #pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index < static_cast<std::ptrdiff_t>(facetCount);
        ++index)
    {
        const std::size_t id =
            static_cast<std::size_t>(index);

        const auto& f = topology.facetGeometry[id];

        const double normalLength = norm(f.normal);

        valid[id] =
            f.area > 0.0 &&
            std::isfinite(f.area) &&
            std::isfinite(normalLength) &&
            normalLength > 0.0;
    }

    // ---------------------------------------------------------
    // 3. Connected planar region growing
    // ---------------------------------------------------------

    std::vector<unsigned char> visited(facetCount, 0);
    std::vector<std::size_t> queue;
    queue.reserve(256);

    std::size_t planarCount = 0;

    for(std::size_t seed = 0; seed < facetCount; ++seed)
    {
        if(visited[seed] || !valid[seed])
            continue;

        const auto& seedFacet =
            topology.facetGeometry[seed];

        const double normalLength =
            norm(seedFacet.normal);

        Plane plane;

        for(int d = 0; d < 3; ++d)
            plane.normal[d] =
                seedFacet.normal[d] / normalLength;

        const std::size_t seedVertex =
            mesh.facetVertexIDs[seed][0];

        plane.offset = dot(
            plane.normal,
            mesh.vertices[seedVertex]);

        queue.clear();
        queue.push_back(seed);
        visited[seed] = 1;

        double regionArea = 0.0;

        for(std::size_t head = 0;
            head < queue.size();
            ++head)
        {
            const std::size_t facetID = queue[head];

            regionArea +=
                topology.facetGeometry[facetID].area;

            for(std::size_t neighbor :
                topology.adjacency[facetID])
            {
                if(neighbor >= facetCount ||
                   visited[neighbor] ||
                   !valid[neighbor])
                    continue;

                if(!facetOnPlane(
                       topology,
                       neighbor,
                       plane,
                       distanceTolerance,
                       normalCosine))
                    continue;

                visited[neighbor] = 1;
                queue.push_back(neighbor);
            }
        }

        // Reject isolated facets and tiny planar patches.
        // Rejected facets remain part of Face 0.
        if(queue.size() < 2 ||
           regionArea < minimumArea)
            continue;

        ++planarCount;

        for(std::size_t facetID : queue)
            result.facetFaceIDs[facetID] = planarCount;
    }

    result.nPlanarFaces = planarCount;

    // ---------------------------------------------------------
    // 4. Extract mesh edges between different face regions
    // ---------------------------------------------------------

    const std::size_t meshEdgeCount =
        topology.edges.size();

    std::vector<unsigned char> isFeature(
        meshEdgeCount, 0);

    #pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index < static_cast<std::ptrdiff_t>(meshEdgeCount);
        ++index)
    {
        const std::size_t id =
            static_cast<std::size_t>(index);

        const auto& edge = topology.edges[id];

        const std::size_t a =
            result.facetFaceIDs[edge.facets[0]];

        const std::size_t b =
            result.facetFaceIDs[edge.facets[1]];

        isFeature[id] = (a != b) ? 1 : 0;
    }

    std::vector<FeatureSegment> segments;
    segments.reserve(meshEdgeCount / 4);

    for(std::size_t id = 0; id < meshEdgeCount; ++id)
    {
        if(!isFeature[id])
            continue;

        const auto& edge = topology.edges[id];

        const std::size_t a =
            result.facetFaceIDs[edge.facets[0]];

        const std::size_t b =
            result.facetFaceIDs[edge.facets[1]];

        segments.push_back({
            edge.v0,
            edge.v1,
            makeFacePair(a, b)
        });
    }

    // ---------------------------------------------------------
    // 5. Build incidence only for feature segments
    // ---------------------------------------------------------

    std::vector<std::vector<std::size_t>> incident(
        vertexCount);

    for(std::size_t id = 0; id < segments.size(); ++id)
    {
        incident[segments[id].v0].push_back(id);
        incident[segments[id].v1].push_back(id);
    }

    // ---------------------------------------------------------
    // 6. Group connected segments into feature curves
    // ---------------------------------------------------------

    const std::size_t unassigned =
        std::numeric_limits<std::size_t>::max();

    std::vector<std::size_t> segmentCurveIDs(
        segments.size(), unassigned);

    std::size_t curveCount = 0;

    // First process open curves and junction-to-junction
    // curves, starting from vertices with degree != 2.
    for(std::size_t id = 0; id < segments.size(); ++id)
    {
        if(segmentCurveIDs[id] != unassigned)
            continue;

        const auto& segment = segments[id];

        const std::size_t degree0 = matchingDegree(
            segment.v0, segment.faces, segments, incident);

        const std::size_t degree1 = matchingDegree(
            segment.v1, segment.faces, segments, incident);

        if(degree0 == 2 && degree1 == 2)
            continue;

        const std::size_t startVertex =
            degree0 != 2 ? segment.v0 : segment.v1;

        traceFeatureCurve(
            id,
            startVertex,
            curveCount,
            segments,
            incident,
            segmentCurveIDs);

        ++curveCount;
    }

    // Remaining segments form closed loops, or interior
    // portions of curves whose endpoints were processed.
    for(std::size_t id = 0; id < segments.size(); ++id)
    {
        if(segmentCurveIDs[id] != unassigned)
            continue;

        traceFeatureCurve(
            id,
            segments[id].v0,
            curveCount,
            segments,
            incident,
            segmentCurveIDs);

        ++curveCount;
    }

    result.nFeatureEdges = curveCount;

    // ---------------------------------------------------------
    // 7. Store only vertex IDs and curve IDs
    // ---------------------------------------------------------

    result.edges.resize(segments.size());

    #pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index < static_cast<std::ptrdiff_t>(segments.size());
        ++index)
    {
        const std::size_t id =
            static_cast<std::size_t>(index);

        result.edges[id] = {
            segments[id].v0,
            segments[id].v1,
            segmentCurveIDs[id]
        };
    }

    boundaryFeatures = std::move(result);
}

} // namespace ntic::lbm::geometry
