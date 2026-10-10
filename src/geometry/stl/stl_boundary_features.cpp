#include "stl_geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <omp.h>

namespace ntic::lbm::geometry {
namespace {

constexpr std::size_t invalidID = std::numeric_limits<std::size_t>::max();
using FacePair = std::array<std::size_t, 2>;

struct Vec3 {
    long double x = 0, y = 0, z = 0;
};

Vec3 diff(const Point& a, const Point& b) {
    return {static_cast<long double>(a[0]) - b[0],
            static_cast<long double>(a[1]) - b[1],
            static_cast<long double>(a[2]) - b[2]};
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}

long double dot(const Vec3& a, const Vec3& b) {
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

long double magnitude(const Vec3& v) {
    return std::sqrt(dot(v,v));
}

struct Plane {
    Point origin{};
    Vec3 normal{};
    long double normalLength = 0;
    long double edgeScale = 0;
};

Plane facetPlane(const stl::FacetTopology& topology, std::size_t facetID) {
    const auto& ids = topology.geometry.facetVertexIDs[facetID];
    const auto& v = topology.geometry.vertices;
    Plane p;
    p.origin = v[ids[0]];
    const Vec3 a = diff(v[ids[1]], p.origin);
    const Vec3 b = diff(v[ids[2]], p.origin);
    p.normal = cross(a,b);
    p.normalLength = magnitude(p.normal);
    p.edgeScale = std::max({magnitude(a), magnitude(b),
                            magnitude(diff(v[ids[2]],v[ids[1]]))});
    return p;
}

// The input vertices are double. The tolerance accounts only for arithmetic
// roundoff of cross/dot products; it is NOT a geometric/angular tolerance.
// The calculation uses translated vectors, so large absolute coordinates
// do not directly inflate the bound. Near-degenerate facets remain sensitive.
bool liesOnPlane(const stl::FacetTopology& topology,
                 std::size_t facetID, const Plane& p) {
    if(!(p.normalLength > 0) || !(p.edgeScale > 0)) return false;
    const auto& v = topology.geometry.vertices;
    const auto& ids = topology.geometry.facetVertexIDs[facetID];
    constexpr long double factor = 64.0L;
    constexpr long double eps = std::numeric_limits<double>::epsilon();
    for(std::size_t id : ids) {
        const Vec3 q = diff(v[id],p.origin);
        const long double qLength = magnitude(q);
        const long double error = factor * eps * p.edgeScale * qLength;
        // Compare unnormalized plane residual to avoid normalization error.
        if(std::abs(dot(p.normal,q)) > p.normalLength * error)
            return false;
    }
    return true;
}

// Check both directions. This is important when triangle sizes differ.
bool mutuallyCoplanar(const stl::FacetTopology& topology,
                      std::size_t a, std::size_t b,
                      const std::vector<Plane>& planes) {
    return liesOnPlane(topology,b,planes[a]) &&
           liesOnPlane(topology,a,planes[b]);
}

FacePair orderedPair(std::size_t a, std::size_t b) {
    if(a>b) std::swap(a,b);
    return {a,b};
}

struct Segment {
    std::size_t v0=0, v1=0;
    FacePair pair{};
};

void traceCurve(std::size_t first, std::size_t start,
                std::size_t curveID,
                const std::vector<Segment>& segments,
                const std::vector<std::vector<std::size_t>>& incident,
                std::vector<std::size_t>& ids) {
    const FacePair pair = segments[first].pair;
    std::size_t edge=first, vertex=start;
    while(ids[edge]==invalidID) {
        ids[edge]=curveID;
        const Segment& s=segments[edge];
        const std::size_t next=(s.v0==vertex)?s.v1:s.v0;
        // All feature curves stop at a global junction, including junctions
        // where the other incident curves have different face pairs.
        if(incident[next].size()!=2) break;
        const auto& adj=incident[next];
        const std::size_t candidate=(adj[0]==edge)?adj[1]:adj[0];
        if(segments[candidate].pair!=pair || ids[candidate]!=invalidID)
            break;
        vertex=next;
        edge=candidate;
    }
}

void buildEdges(const stl::FacetTopology& topology,
                STLBoundaryFeatures& out) {
    const std::size_t edgeCount=topology.edges.size();
    std::vector<unsigned char> marked(edgeCount,0);
    #pragma omp parallel for schedule(static)
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(edgeCount);++i) {
        const auto& e=topology.edges[static_cast<std::size_t>(i)];
        const auto a=out.facetFaceIDs[e.facets[0]];
        const auto b=out.facetFaceIDs[e.facets[1]];
        marked[static_cast<std::size_t>(i)]=(a!=b);
    }
    std::vector<Segment> segments;
    for(std::size_t i=0;i<edgeCount;++i) {
        if(!marked[i]) continue;
        const auto& e=topology.edges[i];
        segments.push_back({e.v0,e.v1,orderedPair(
            out.facetFaceIDs[e.facets[0]],out.facetFaceIDs[e.facets[1]])});
    }
    std::vector<std::vector<std::size_t>> incident(
        topology.geometry.vertices.size());
    for(std::size_t i=0;i<segments.size();++i) {
        incident[segments[i].v0].push_back(i);
        incident[segments[i].v1].push_back(i);
    }
    std::vector<std::size_t> ids(segments.size(),invalidID);
    std::size_t count=0;
    // Open curves and junction-to-junction paths first.
    for(std::size_t i=0;i<segments.size();++i) {
        if(ids[i]!=invalidID) continue;
        const auto& s=segments[i];
        const auto endpoint=[&](std::size_t v) {
            if(incident[v].size()!=2) return true;
            return segments[incident[v][0]].pair !=
                   segments[incident[v][1]].pair;
        };
        const bool a=endpoint(s.v0), b=endpoint(s.v1);
        if(!a && !b) continue;
        traceCurve(i,a?s.v0:s.v1,count++,segments,incident,ids);
    }
    // Remaining unassigned segments are closed loops.
    for(std::size_t i=0;i<segments.size();++i) {
        if(ids[i]==invalidID)
            traceCurve(i,segments[i].v0,count++,segments,incident,ids);
    }
    out.nFeatureEdges=count;
    out.edges.resize(segments.size());
    #pragma omp parallel for schedule(static)
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(segments.size());++i) {
        const std::size_t j=static_cast<std::size_t>(i);
        out.edges[j]={segments[j].v0,segments[j].v1,ids[j]};
    }
}

} // anonymous namespace

void STLGeometry::identifyBoundaryFeatures(
    const std::vector<std::size_t>& excludedFaceIDs,
    double smallFaceAreaRatio) {
    if(!std::isfinite(smallFaceAreaRatio) ||
       smallFaceAreaRatio<0 || smallFaceAreaRatio>1)
        throw std::invalid_argument("smallFaceAreaRatio must be in [0,1].");

    STLBoundaryFeatures out;
    const std::size_t count=topology.geometry.facetVertexIDs.size();
    out.facetFaceIDs.assign(count,0);
    if(flowType==FlowType::External) {
        if(!excludedFaceIDs.empty())
            throw std::invalid_argument("External flow has no planar faces to exclude.");
        boundaryFeatures=std::move(out);
        return;
    }
    if(topology.facetGeometry.size()!=count || topology.adjacency.size()!=count)
        throw std::runtime_error("Incomplete STL facet geometry or adjacency.");

    std::vector<Plane> planes(count);
    std::vector<unsigned char> valid(count,0);
    #pragma omp parallel for schedule(static)
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(count);++i) {
        const std::size_t id=static_cast<std::size_t>(i);
        planes[id]=facetPlane(topology,id);
        valid[id]=(planes[id].normalLength>0 &&
                   std::isfinite(planes[id].normalLength));
    }

    // A planar region requires at least two edge-adjacent coplanar facets.
    // Isolated curved-surface triangles never become their own planar faces.
    std::vector<unsigned char> seedCandidate(count,0);
    #pragma omp parallel for schedule(static)
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(count);++i) {
        const std::size_t id=static_cast<std::size_t>(i);
        if(!valid[id]) continue;
        for(std::size_t nb:topology.adjacency[id]) {
            if(nb<count && valid[nb] && mutuallyCoplanar(topology,id,nb,planes)) {
                seedCandidate[id]=1;
                break;
            }
        }
    }

    std::vector<unsigned char> assigned(count,0);
    std::vector<std::size_t> queue;
    queue.reserve(256);
    std::vector<STLPlanarFaceInfo> detected;
    for(std::size_t seed=0;seed<count;++seed) {
        if(assigned[seed] || !seedCandidate[seed]) continue;
        // Ensure there is an available coplanar neighbor. If a previous
        // accepted region consumed it, leave this triangle in Face 0.
        bool supported=false;
        for(std::size_t nb:topology.adjacency[seed]) {
            if(nb<count && !assigned[nb] && valid[nb] &&
               mutuallyCoplanar(topology,seed,nb,planes)) {
                supported=true;
                break;
            }
        }
        if(!supported) continue;

        queue.clear();
        queue.push_back(seed);
        assigned[seed]=1;
        double area=0;
        for(std::size_t head=0;head<queue.size();++head) {
            const std::size_t f=queue[head];
            area+=topology.facetGeometry[f].area;
            for(std::size_t nb:topology.adjacency[f]) {
                if(nb>=count || assigned[nb] || !valid[nb]) continue;
                // Every facet must match the FIXED seed plane, not just
                // the preceding facet, to avoid curvature drift.
                if(!mutuallyCoplanar(topology,seed,nb,planes)) continue;
                assigned[nb]=1;
                queue.push_back(nb);
            }
        }
        // Defensive guard: never accept a single-facet plane.
        if(queue.size()<2) {
            assigned[seed]=0;
            continue;
        }
        const std::size_t originalID=detected.size()+1;
        detected.push_back({originalID,originalID,queue.size(),area,false});
        for(std::size_t f:queue) out.facetFaceIDs[f]=originalID;
    }

    std::vector<unsigned char> excluded(detected.size()+1,0);
    for(std::size_t id:excludedFaceIDs) {
        if(id==0 || id>detected.size())
            throw std::invalid_argument("Excluded original Face ID out of range.");
        excluded[id]=1;
    }
    double totalArea=0;
    for(const auto& f:topology.facetGeometry) totalArea+=f.area;
    std::vector<std::size_t> remap(detected.size()+1,0);
    for(const auto& f:detected) {
        if(excluded[f.originalFaceID]) continue;
        const std::size_t newID=out.planarFaces.size()+1;
        remap[f.originalFaceID]=newID;
        STLPlanarFaceInfo info=f;
        info.faceID=newID;
        info.small=(totalArea>0 && info.area/totalArea<smallFaceAreaRatio);
        out.planarFaces.push_back(info);
    }
    #pragma omp parallel for schedule(static)
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(count);++i) {
        const std::size_t j=static_cast<std::size_t>(i);
        out.facetFaceIDs[j]=remap[out.facetFaceIDs[j]];
    }
    out.nPlanarFaces=out.planarFaces.size();
    buildEdges(topology,out);
    boundaryFeatures=std::move(out);
}

void STLGeometry::printBoundaryFeatureReport() const {
    const auto& f=boundaryFeatures;
    double totalArea=0;
    for(const auto& facet:topology.facetGeometry) totalArea+=facet.area;
    std::size_t nonplanar=0;
    for(std::size_t id:f.facetFaceIDs) nonplanar+=(id==0);
    std::cout << "\nSTL Boundary Feature Report\n"
              << "  Planar Faces: " << f.nPlanarFaces << '\n'
              << "  Feature Curves: " << f.nFeatureEdges << '\n'
              << "  Feature Segments: " << f.edges.size() << '\n'
              << "  Nonplanar Facets: " << nonplanar << '\n';
    for(const auto& face:f.planarFaces) {
        const double ratio=totalArea>0?face.area/totalArea:0;
        std::cout << "  Face " << face.faceID
                  << " (original " << face.originalFaceID << ")"
                  << ": facets=" << face.facetCount
                  << ", area=" << std::setprecision(12) << face.area
                  << ", areaRatio=" << ratio;
        if(face.small) std::cout << "  WARNING: small planar region; review before retaining.";
        std::cout << '\n';
    }
}

} // namespace ntic::lbm::geometry
