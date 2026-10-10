#include "flow_type.hpp"
#include "stl_geometry.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
using ntic::lbm::geometry::BoundingBox;
using ntic::lbm::geometry::FlowType;
using ntic::lbm::geometry::Point;
using ntic::lbm::geometry::STLGeometry;
using ntic::lbm::geometry::STLBoundaryFeatures;

void require(bool condition, const std::string& message)
{
    if(!condition) throw std::runtime_error(message);
}

void requirePointEqual(const Point& a, const Point& b,
                       const std::string& message)
{
    for(std::size_t d = 0; d < 3; ++d)
        require(a[d] == b[d], message);
}

BoundingBox makeReferenceBounds(const ntic::lbm::stl::FacetTopology& topology)
{
    BoundingBox box;
    for(const auto& p : topology.geometry.vertices) box.expand(p);
    return box;
}

using EdgeKey = std::pair<std::size_t, std::size_t>;

EdgeKey edgeKey(std::size_t a, std::size_t b)
{
    if(a > b) std::swap(a, b);
    return {a, b};
}

bool nearlyEqual(double a, double b)
{
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= 1.0e-11 * scale;
}

void validateBoundaryFeatures(const STLGeometry& geometry)
{
    const auto& topology = geometry.topology;
    const auto& features = geometry.boundaryFeatures;
    const std::size_t nFacets = topology.geometry.facetVertexIDs.size();
    const std::size_t nVertices = topology.geometry.vertices.size();

    require(features.facetFaceIDs.size() == nFacets,
            "Facet Face ID array has incorrect size.");
    require(features.planarFaces.size() == features.nPlanarFaces,
            "Planar Face metadata count is incorrect.");

    std::vector<std::size_t> facetCounts(features.nPlanarFaces + 1, 0);
    std::vector<double> faceAreas(features.nPlanarFaces + 1, 0.0);
    double totalArea = 0.0;

    for(std::size_t f = 0; f < nFacets; ++f)
    {
        const std::size_t id = features.facetFaceIDs[f];
        require(id <= features.nPlanarFaces, "Facet Face ID out of range.");
        ++facetCounts[id];
        faceAreas[id] += topology.facetGeometry[f].area;
        totalArea += topology.facetGeometry[f].area;
    }

    std::size_t previousOriginalID = 0;
    for(std::size_t i = 0; i < features.planarFaces.size(); ++i)
    {
        const auto& face = features.planarFaces[i];
        require(face.faceID == i + 1, "Retained Face IDs are not contiguous.");
        require(face.originalFaceID > previousOriginalID,
                "Original Face IDs must be strictly increasing.");
        previousOriginalID = face.originalFaceID;
        require(face.facetCount >= 2, "Planar Face contains fewer than 2 facets.");
        require(facetCounts[face.faceID] == face.facetCount,
                "Planar Face facet count mismatch.");
        require(nearlyEqual(faceAreas[face.faceID], face.area),
                "Planar Face area mismatch.");
        require(std::isfinite(face.area) && face.area > 0.0,
                "Planar Face area must be positive and finite.");
    }

    // The STL topology edge lookup is O(E log E), not O(E_feature * E).
    std::map<EdgeKey, std::pair<std::size_t, std::size_t>> meshEdges;
    std::size_t expectedSegments = 0;
    for(const auto& e : topology.edges)
    {
        const EdgeKey key = edgeKey(e.v0, e.v1);
        require(meshEdges.emplace(key, std::make_pair(e.facets[0], e.facets[1])).second,
                "Duplicate geometric edge in validated STL topology.");
        if(features.facetFaceIDs[e.facets[0]] !=
           features.facetFaceIDs[e.facets[1]])
            ++expectedSegments;
    }
    require(features.edges.size() == expectedSegments,
            "Feature segment count differs from face-interface edge count.");

    std::vector<std::size_t> curveCounts(features.nFeatureEdges, 0);
    std::map<EdgeKey, bool> seen;
    for(const auto& edge : features.edges)
    {
        require(edge.v0 < nVertices && edge.v1 < nVertices && edge.v0 != edge.v1,
                "Invalid feature segment endpoints.");
        require(edge.edgeID < features.nFeatureEdges,
                "Feature curve ID out of range.");
        const EdgeKey key = edgeKey(edge.v0, edge.v1);
        const auto it = meshEdges.find(key);
        require(it != meshEdges.end(), "Feature segment is not an STL mesh edge.");
        require(seen.emplace(key, true).second,
                "Duplicate feature segment.");
        const auto facets = it->second;
        require(features.facetFaceIDs[facets.first] !=
                features.facetFaceIDs[facets.second],
                "Feature segment does not separate different Faces.");
        ++curveCounts[edge.edgeID];
    }
    for(std::size_t count : curveCounts)
        require(count > 0, "Empty feature curve group.");

    // Every segment of one curve must separate the same unordered Face pair.
    std::vector<EdgeKey> curvePairs(features.nFeatureEdges,
                                    {std::numeric_limits<std::size_t>::max(),
                                     std::numeric_limits<std::size_t>::max()});
    for(const auto& edge : features.edges)
    {
        const auto facets = meshEdges.at(edgeKey(edge.v0, edge.v1));
        const EdgeKey pair = edgeKey(features.facetFaceIDs[facets.first],
                                     features.facetFaceIDs[facets.second]);
        auto& recorded = curvePairs[edge.edgeID];
        if(recorded.first == std::numeric_limits<std::size_t>::max())
            recorded = pair;
        else
            require(recorded == pair,
                    "One feature curve crosses different Face pairs.");
    }

    std::cout << "\nBoundary feature validation passed:\n"
              << "  Total facets: " << nFacets << '\n'
              << "  Nonplanar facets: " << facetCounts[0] << '\n'
              << "  Planar Faces: " << features.nPlanarFaces << '\n'
              << "  Feature curves: " << features.nFeatureEdges << '\n'
              << "  Feature segments: " << features.edges.size() << '\n'
              << "  Total surface area: " << totalArea << '\n';
}

void requireSameFeatures(const STLBoundaryFeatures& a,
                         const STLBoundaryFeatures& b)
{
    require(a.facetFaceIDs == b.facetFaceIDs,
            "Repeated detection changed facet Face IDs.");
    require(a.nPlanarFaces == b.nPlanarFaces &&
            a.nFeatureEdges == b.nFeatureEdges,
            "Repeated detection changed feature counts.");
    require(a.edges.size() == b.edges.size(),
            "Repeated detection changed segment count.");
    for(std::size_t i = 0; i < a.edges.size(); ++i)
        require(a.edges[i].v0 == b.edges[i].v0 &&
                a.edges[i].v1 == b.edges[i].v1 &&
                a.edges[i].edgeID == b.edges[i].edgeID,
                "Repeated detection changed feature segment assignments.");
    require(a.planarFaces.size() == b.planarFaces.size(),
            "Repeated detection changed planar metadata count.");
    for(std::size_t i = 0; i < a.planarFaces.size(); ++i)
        require(a.planarFaces[i].faceID == b.planarFaces[i].faceID &&
                a.planarFaces[i].originalFaceID == b.planarFaces[i].originalFaceID &&
                a.planarFaces[i].facetCount == b.planarFaces[i].facetCount &&
                nearlyEqual(a.planarFaces[i].area, b.planarFaces[i].area),
                "Repeated detection changed planar metadata.");
}

void testExclusion(STLGeometry& geometry)
{
    const STLBoundaryFeatures baseline = geometry.boundaryFeatures;
    if(baseline.nPlanarFaces == 0)
    {
        std::cout << "Exclusion test skipped: no planar Faces detected.\n";
        return;
    }

    const std::size_t excludedOriginalID = baseline.planarFaces.front().originalFaceID;
    std::vector<unsigned char> removed(geometry.topology.geometry.facetVertexIDs.size(), 0);
    for(std::size_t i = 0; i < removed.size(); ++i)
        removed[i] = baseline.facetFaceIDs[i] == baseline.planarFaces.front().faceID;

    geometry.identifyBoundaryFeatures({excludedOriginalID});
    validateBoundaryFeatures(geometry);
    require(geometry.boundaryFeatures.nPlanarFaces + 1 == baseline.nPlanarFaces,
            "Excluding one Face did not decrease planar Face count by one.");
    for(std::size_t i = 0; i < removed.size(); ++i)
    {
        if(removed[i])
            require(geometry.boundaryFeatures.facetFaceIDs[i] == 0,
                    "Excluded Face facet was not reassigned to Face 0.");
    }
    for(const auto& face : geometry.boundaryFeatures.planarFaces)
        require(face.originalFaceID != excludedOriginalID,
                "Excluded Face still exists in retained Face metadata.");

    // A second detection must restore the exact original grouping.
    geometry.identifyBoundaryFeatures();
    validateBoundaryFeatures(geometry);
    requireSameFeatures(baseline, geometry.boundaryFeatures);
    std::cout << "Exclusion and restoration test passed (original Face "
              << excludedOriginalID << ").\n";
}

void testConstruction(const std::string& stlFile)
{
    using namespace ntic::lbm;
    stl::STLData data = stl::read(stlFile);
    stl::FacetTopology topology;
    stl::validate(data, "full", topology);
    require(!topology.geometry.vertices.empty(), "Validated STL has no vertices.");

    const BoundingBox referenceBounds = makeReferenceBounds(topology);
    const std::size_t vertexCount = topology.geometry.vertices.size();
    const std::size_t facetCount = topology.geometry.facetVertexIDs.size();
    const std::size_t facetGeometryCount = topology.facetGeometry.size();
    const std::size_t adjacencyCount = topology.adjacency.size();
    const std::size_t edgeCount = topology.edges.size();
    const std::size_t componentIDCount = topology.facetComponentIDs.size();
    const std::size_t componentCount = topology.components.size();
    const Point firstVertex = topology.geometry.vertices.front();
    const Point middleVertex = topology.geometry.vertices[vertexCount / 2];
    const Point lastVertex = topology.geometry.vertices.back();

    STLGeometry geometry(std::move(topology), FlowType::Internal);
    require(geometry.topology.geometry.vertices.size() == vertexCount,
            "Vertex count changed during construction.");
    require(geometry.topology.geometry.facetVertexIDs.size() == facetCount,
            "Facet count changed during construction.");
    require(geometry.topology.facetGeometry.size() == facetGeometryCount,
            "Facet geometry count changed during construction.");
    require(geometry.topology.adjacency.size() == adjacencyCount,
            "Adjacency count changed during construction.");
    require(geometry.topology.edges.size() == edgeCount,
            "Edge count changed during construction.");
    require(geometry.topology.facetComponentIDs.size() == componentIDCount,
            "Facet component IDs changed during construction.");
    require(geometry.topology.components.size() == componentCount,
            "Component count changed during construction.");
    requirePointEqual(geometry.topology.geometry.vertices.front(), firstVertex,
                      "First vertex changed during construction.");
    requirePointEqual(geometry.topology.geometry.vertices[vertexCount / 2], middleVertex,
                      "Middle vertex changed during construction.");
    requirePointEqual(geometry.topology.geometry.vertices.back(), lastVertex,
                      "Last vertex changed during construction.");
    requirePointEqual(geometry.bounds.min, referenceBounds.min,
                      "Minimum bound mismatch.");
    requirePointEqual(geometry.bounds.max, referenceBounds.max,
                      "Maximum bound mismatch.");
    for(const auto& vertex : geometry.topology.geometry.vertices)
        require(geometry.bounds.contains(vertex),
                "Bounding box does not contain an STL vertex.");

    geometry.printBoundaryFeatureReport();
    validateBoundaryFeatures(geometry);
    const STLBoundaryFeatures baseline = geometry.boundaryFeatures;
    geometry.identifyBoundaryFeatures();
    requireSameFeatures(baseline, geometry.boundaryFeatures);
    std::cout << "Repeated detection test passed.\n";
    testExclusion(geometry);
}

void testEmptyGeometry()
{
    ntic::lbm::stl::FacetTopology topology;
    bool caught = false;
    try { STLGeometry geometry(std::move(topology), FlowType::Internal); }
    catch(const std::runtime_error&) { caught = true; }
    require(caught, "Empty STL geometry construction did not throw.");
}
}

int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <stl_file>\n";
        return 1;
    }
    try
    {
        testConstruction(argv[1]);
        testEmptyGeometry();
        std::cout << "STL geometry and boundary feature tests passed.\n";
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr << "STL geometry tests failed: " << e.what() << '\n';
        return 1;
    }
}
