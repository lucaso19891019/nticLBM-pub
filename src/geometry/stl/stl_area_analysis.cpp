#include "stl_geometry.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include <omp.h>

namespace ntic::lbm::geometry
{
namespace
{

Point subtract(const Point& a, const Point& b)
{
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

double dot(const Point& a, const Point& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

Point cross(const Point& a, const Point& b)
{
    return {a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]};
}

double norm(const Point& a)
{
    return std::sqrt(dot(a,a));
}

double triangleSolidAngle(const Point& point,
                          const Point& v0,
                          const Point& v1,
                          const Point& v2)
{
    const Point a = subtract(v0,point);
    const Point b = subtract(v1,point);
    const Point c = subtract(v2,point);
    const double la = norm(a), lb = norm(b), lc = norm(c);
    const double numerator = dot(a,cross(b,c));
    const double denominator = la*lb*lc + dot(a,b)*lc
                             + dot(b,c)*la + dot(c,a)*lb;
    return 2.0*std::atan2(numerator,denominator);
}

bool pointInComponent(const STLGeometry& geometry,
                      std::size_t componentID,
                      const Point& point)
{
    const auto& topology = geometry.topology;
    const auto& facets = topology.components[componentID].facets;
    double solidAngle = 0.0;

#pragma omp parallel for reduction(+:solidAngle) schedule(static) if(!omp_in_parallel())
    for(std::int64_t i=0; i<static_cast<std::int64_t>(facets.size()); ++i)
    {
        const auto& ids = topology.geometry.facetVertexIDs[
            facets[static_cast<std::size_t>(i)]];
        const auto& vertices = topology.geometry.vertices;
        solidAngle += triangleSolidAngle(point,
                                         vertices[ids[0]],
                                         vertices[ids[1]],
                                         vertices[ids[2]]);
    }
    constexpr double pi = 3.14159265358979323846;
    return std::abs(solidAngle) > 2.0*pi;
}

bool centerIsWet(const STLGeometry& geometry, const Point& point)
{
    for(std::size_t componentID=0;
        componentID<geometry.topology.components.size(); ++componentID)
    {
        const auto& componentFlow = geometry.flow[componentID];
        if(!componentFlow.active) continue;
        const auto& bounds = geometry.topology.components[componentID].bounds;
        const bool outsideBounds =
            point[0]<bounds.min[0] || point[0]>bounds.max[0] ||
            point[1]<bounds.min[1] || point[1]>bounds.max[1] ||
            point[2]<bounds.min[2] || point[2]>bounds.max[2];
        const bool inside = !outsideBounds &&
            pointInComponent(geometry,componentID,point);
        if(componentFlow.fluidSide == FluidSide::Inside)
        {
            if(!inside) return false;
        }
        else if(inside) return false;
    }
    return true;
}

// Closed triangle/axis-aligned cube SAT. Touching is an intersection.
bool separatedOnAxisClosed(const Point& a,
                           const Point& b,
                           const Point& c,
                           const Point& axis,
                           double halfWidth)
{
    if(dot(axis,axis)==0.0) return false;
    const double pa=dot(a,axis), pb=dot(b,axis), pc=dot(c,axis);
    const double radius=halfWidth*(std::abs(axis[0])+
                                    std::abs(axis[1])+
                                    std::abs(axis[2]));
    return std::max({pa,pb,pc}) < -radius ||
           std::min({pa,pb,pc}) > radius;
}

bool triangleIntersectsCell(const Point& center,
                            double halfWidth,
                            const Point& v0,
                            const Point& v1,
                            const Point& v2)
{
    const Point a=subtract(v0,center);
    const Point b=subtract(v1,center);
    const Point c=subtract(v2,center);
    for(int axis=0; axis<3; ++axis)
    {
        if(std::max({a[axis],b[axis],c[axis]}) < -halfWidth ||
           std::min({a[axis],b[axis],c[axis]}) > halfWidth)
            return false;
    }
    const Point e0=subtract(b,a), e1=subtract(c,b), e2=subtract(a,c);
    if(separatedOnAxisClosed(a,b,c,cross(e0,e1),halfWidth))
        return false;
    const std::array<Point,3> edges{{e0,e1,e2}};
    const std::array<Point,3> axes{{
        Point{1.0,0.0,0.0},Point{0.0,1.0,0.0},Point{0.0,0.0,1.0}}};
    for(const auto& edge: edges)
        for(const auto& axis: axes)
            if(separatedOnAxisClosed(a,b,c,cross(edge,axis),halfWidth))
                return false;
    return true;
}

// Conservative center-on-surface classification for finite doubles.
// A center sufficiently close to a triangle is treated as Dry.
// This intentionally allows resolution-limited Wet nodes to be removed.
bool pointOnTriangle(const Point& p,
                     const Point& a,
                     const Point& b,
                     const Point& c)
{
    const Point ab=subtract(b,a), ac=subtract(c,a), ap=subtract(p,a);
    const Point n=cross(ab,ac);
    const double n2=dot(n,n);
    if(n2==0.0) return false;

    const double scale=std::max({norm(ab),norm(ac),norm(subtract(c,b))});
    const double tolerance=32.0*std::numeric_limits<double>::epsilon()*scale;
    if(std::abs(dot(ap,n)) > tolerance*std::sqrt(n2)) return false;

    const double d00=dot(ab,ab), d01=dot(ab,ac), d11=dot(ac,ac);
    const double d20=dot(ap,ab), d21=dot(ap,ac);
    const double denominator=d00*d11-d01*d01;
    if(denominator<=0.0) return false;

    const double u=(d11*d20-d01*d21)/denominator;
    const double v=(d00*d21-d01*d20)/denominator;
    const double baryTolerance=32.0*std::numeric_limits<double>::epsilon();
    return u>=-baryTolerance && v>=-baryTolerance &&
           u+v<=1.0+baryTolerance;
}

// Candidate grid cells whose CLOSED boxes of half-width 'halfWidth'
// may intersect the facet's axis-aligned interval [facetMin,facetMax].
bool facetCellIndexRange(double facetMin,
                         double facetMax,
                         double gridMin,
                         double spacing,
                         double halfWidth,
                         std::size_t count,
                         std::size_t& first,
                         std::size_t& last)
{
    const double lower=(facetMin-halfWidth-gridMin)/spacing-0.5;
    const double upper=(facetMax+halfWidth-gridMin)/spacing-0.5;
    // Clamp before integer conversion, also handling facets outside domain.
    if(upper<0.0 || lower>static_cast<double>(count-1)) return false;
    first=static_cast<std::size_t>(std::max(0.0,std::ceil(lower)));
    last=static_cast<std::size_t>(std::min(
        static_cast<double>(count-1),std::floor(upper)));
    return first<=last;
}

Point cellCenter(const BoundingBox& domain,
                 double spacing,
                 std::size_t i,
                 std::size_t j,
                 std::size_t k)
{
    return {domain.min[0]+(static_cast<double>(i)+0.5)*spacing,
            domain.min[1]+(static_cast<double>(j)+0.5)*spacing,
            domain.min[2]+(static_cast<double>(k)+0.5)*spacing};
}

struct FacetData
{
    Point a,b,c;
    Point min,max;
};

FacetData facetData(const STLGeometry& geometry, std::size_t facetID)
{
    const auto& ids=geometry.topology.geometry.facetVertexIDs[facetID];
    const auto& vertices=geometry.topology.geometry.vertices;
    FacetData f{vertices[ids[0]],vertices[ids[1]],vertices[ids[2]],{}, {}};
    for(int d=0; d<3; ++d)
    {
        f.min[d]=std::min({f.a[d],f.b[d],f.c[d]});
        f.max[d]=std::max({f.a[d],f.b[d],f.c[d]});
    }
    return f;
}

bool facetRange(const FacetData& f,
                const BoundingBox& domain,
                double spacing,
                double halfWidth,
                std::size_t nx,
                std::size_t ny,
                std::size_t nz,
                std::array<std::size_t,3>& lo,
                std::array<std::size_t,3>& hi)
{
    const std::array<std::size_t,3> dimensions{{nx,ny,nz}};
    for(int d=0; d<3; ++d)
        if(!facetCellIndexRange(f.min[d],f.max[d],domain.min[d],
                                spacing,halfWidth,dimensions[d],lo[d],hi[d]))
            return false;
    return true;
}

} // namespace

template <lattice::LatticeType Type>
void STLGeometry::interiorAreaAnalysis(
    GeometryAnalysis3D<Type>& analysis) const
{
    const double gridSpacing = analysis.gridSpacing;
    BoundingBox& domain = analysis.domain;
    std::size_t& nx = analysis.nx;
    std::size_t& ny = analysis.ny;
    std::size_t& nz = analysis.nz;
    std::vector<double>& scalar = analysis.scalar;
    std::vector<double>& boundaryX = analysis.boundaryX;
    std::vector<double>& boundaryY = analysis.boundaryY;
    std::vector<double>& boundaryZ = analysis.boundaryZ;

    if(!(gridSpacing>0.0) || !std::isfinite(gridSpacing))
        throw std::runtime_error("Grid spacing must be positive and finite.");
    if(flow.size()!=topology.components.size())
        throw std::runtime_error("STL flow interpretation is not available.");

    domain=(flowType==FlowType::Internal) ? bounds : openBox;
    nx=static_cast<std::size_t>(std::ceil(domain.width()/gridSpacing));
    ny=static_cast<std::size_t>(std::ceil(domain.height()/gridSpacing));
    nz=static_cast<std::size_t>(std::ceil(domain.depth()/gridSpacing));
    if(nx==0 || ny==0 || nz==0)
        throw std::runtime_error("STL analysis domain contains no grid cells.");
    if(nx>std::numeric_limits<std::size_t>::max()/ny ||
       nx*ny>std::numeric_limits<std::size_t>::max()/nz)
        throw std::runtime_error("STL grid cell count overflow.");

    const std::size_t xySize=nx*ny, cellCount=xySize*nz;
    const std::size_t facetCount=topology.geometry.facetVertexIDs.size();
    std::vector<CellType> cellTypes(cellCount,CellType::Dry);
    boundaryX.clear(); boundaryY.clear(); boundaryZ.clear();

    // Step 1: Facet-driven closed-voxel rasterization. Touch is Boundary.
    // 'initialBoundary' is immutable after this step and is the flood barrier.
    std::vector<std::atomic<std::uint8_t>> atomicBoundary(cellCount);
    std::vector<std::atomic<std::uint8_t>> atomicOnSurface(cellCount);
#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t id=0; id<static_cast<std::ptrdiff_t>(cellCount); ++id)
    {
        atomicBoundary[static_cast<std::size_t>(id)].store(0,std::memory_order_relaxed);
        atomicOnSurface[static_cast<std::size_t>(id)].store(0,std::memory_order_relaxed);
    }

#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index=0; index<static_cast<std::ptrdiff_t>(facetCount); ++index)
    {
        const std::size_t facetID=static_cast<std::size_t>(index);
        const std::size_t componentID=topology.facetComponentIDs[facetID];
        if(!flow[componentID].active) continue;
        const FacetData f=facetData(*this,facetID);
        std::array<std::size_t,3> lo{},hi{};
        if(!facetRange(f,domain,gridSpacing,0.5*gridSpacing,
                       nx,ny,nz,lo,hi)) continue;
        for(std::size_t k=lo[2]; k<=hi[2]; ++k)
            for(std::size_t j=lo[1]; j<=hi[1]; ++j)
                for(std::size_t i=lo[0]; i<=hi[0]; ++i)
                {
                    const std::size_t id=i+nx*(j+ny*k);
                    // Do not skip an already marked cell: another facet
                    // may pass exactly through its center.
                    const Point c=cellCenter(domain,gridSpacing,i,j,k);
                    if(!triangleIntersectsCell(c,0.5*gridSpacing,
                                               f.a,f.b,f.c)) continue;
                    atomicBoundary[id].store(1,std::memory_order_relaxed);
                    if(pointOnTriangle(c,f.a,f.b,f.c))
                        atomicOnSurface[id].store(1,std::memory_order_relaxed);
                }
    }

    std::vector<std::uint8_t> initialBoundary(cellCount,0);
    std::vector<std::uint8_t> onSurface(cellCount,0);
#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index=0; index<static_cast<std::ptrdiff_t>(cellCount); ++index)
    {
        const std::size_t id=static_cast<std::size_t>(index);
        initialBoundary[id]=atomicBoundary[id].load(std::memory_order_relaxed);
        onSurface[id]=atomicOnSurface[id].load(std::memory_order_relaxed);
        if(initialBoundary[id]) cellTypes[id]=CellType::Boundary;
    }

    // Step 2: Flood Fill only non-surface voxels (6-connectivity).
    std::vector<std::uint8_t> visited(cellCount,0);
    std::vector<std::size_t> queue;
    queue.reserve(cellCount);
    const std::array<std::array<int,3>,6> offsets{{
        {{-1,0,0}},{{1,0,0}},{{0,-1,0}},{{0,1,0}},{{0,0,-1}},{{0,0,1}}}};
    for(std::size_t seed=0; seed<cellCount; ++seed)
    {
        if(initialBoundary[seed] || visited[seed]) continue;
        const std::size_t sk=seed/xySize, sj=(seed%xySize)/nx, si=seed%nx;
        const bool wet=centerIsWet(*this,cellCenter(domain,gridSpacing,si,sj,sk));
        const CellType type=wet ? CellType::Interior : CellType::Dry;
        queue.clear(); queue.push_back(seed); visited[seed]=1;
        for(std::size_t head=0; head<queue.size(); ++head)
        {
            const std::size_t id=queue[head];
            cellTypes[id]=type;
            const std::size_t k=id/xySize, j=(id%xySize)/nx, i=id%nx;
            for(const auto& o: offsets)
            {
                const std::ptrdiff_t ni=static_cast<std::ptrdiff_t>(i)+o[0];
                const std::ptrdiff_t nj=static_cast<std::ptrdiff_t>(j)+o[1];
                const std::ptrdiff_t nk=static_cast<std::ptrdiff_t>(k)+o[2];
                if(ni<0 || nj<0 || nk<0 ||
                   ni>=static_cast<std::ptrdiff_t>(nx) ||
                   nj>=static_cast<std::ptrdiff_t>(ny) ||
                   nk>=static_cast<std::ptrdiff_t>(nz)) continue;
                const std::size_t neighbor=static_cast<std::size_t>(ni)+
                    nx*(static_cast<std::size_t>(nj)+ny*static_cast<std::size_t>(nk));
                if(initialBoundary[neighbor] || visited[neighbor]) continue;
                visited[neighbor]=1;
                queue.push_back(neighbor);
            }
        }
    }

    // Step 3.1: Generate unique Wet candidates from ALL initial Boundary
    // voxels, before any strict-Wet filtering. The 3x3x3 closed-voxel
    // neighborhood covers the half-width-h StencilCell of each Wet center.
    std::vector<std::atomic<std::uint8_t>> candidates(cellCount);
#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index=0; index<static_cast<std::ptrdiff_t>(cellCount); ++index)
        candidates[static_cast<std::size_t>(index)].store(0,std::memory_order_relaxed);

#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index=0; index<static_cast<std::ptrdiff_t>(cellCount); ++index)
    {
        const std::size_t id=static_cast<std::size_t>(index);
        if(!initialBoundary[id]) continue;
        const std::size_t k=id/xySize, j=(id%xySize)/nx, i=id%nx;
        const std::size_t i0=(i>0)?i-1:i, i1=std::min(i+1,nx-1);
        const std::size_t j0=(j>0)?j-1:j, j1=std::min(j+1,ny-1);
        const std::size_t k0=(k>0)?k-1:k, k1=std::min(k+1,nz-1);
        for(std::size_t kk=k0; kk<=k1; ++kk)
            for(std::size_t jj=j0; jj<=j1; ++jj)
                for(std::size_t ii=i0; ii<=i1; ++ii)
                {
                    const std::size_t neighbor=ii+nx*(jj+ny*kk);
                    if(cellTypes[neighbor]==CellType::Interior)
                        candidates[neighbor].store(1,std::memory_order_relaxed);
                }
    }

    // Step 3.2: Initial Boundary is final Boundary only for a STRICT Wet
    // center. Points exactly on any active facet are always Dry.
#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index=0; index<static_cast<std::ptrdiff_t>(cellCount); ++index)
    {
        const std::size_t id=static_cast<std::size_t>(index);
        if(!initialBoundary[id]) continue;
        if(onSurface[id])
        {
            cellTypes[id]=CellType::Dry;
            continue;
        }
        const std::size_t k=id/xySize, j=(id%xySize)/nx, i=id%nx;
        const Point c=cellCenter(domain,gridSpacing,i,j,k);
        cellTypes[id]=centerIsWet(*this,c) ? CellType::Boundary : CellType::Dry;
    }

    // Step 3.3: Facet-driven CLOSED StencilCell SAT, only for unique Wet
    // candidates. Multiple facets may hit the same candidate; atomic writes
    // only record existence of at least one hit.
    std::vector<std::atomic<std::uint8_t>> atomicHits(cellCount);
#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index=0; index<static_cast<std::ptrdiff_t>(cellCount); ++index)
        atomicHits[static_cast<std::size_t>(index)].store(0,std::memory_order_relaxed);

#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index=0; index<static_cast<std::ptrdiff_t>(facetCount); ++index)
    {
        const std::size_t facetID=static_cast<std::size_t>(index);
        const std::size_t componentID=topology.facetComponentIDs[facetID];
        if(!flow[componentID].active) continue;
        const FacetData f=facetData(*this,facetID);
        std::array<std::size_t,3> lo{},hi{};
        if(!facetRange(f,domain,gridSpacing,gridSpacing,
                       nx,ny,nz,lo,hi)) continue;
        for(std::size_t k=lo[2]; k<=hi[2]; ++k)
            for(std::size_t j=lo[1]; j<=hi[1]; ++j)
                for(std::size_t i=lo[0]; i<=hi[0]; ++i)
                {
                    const std::size_t id=i+nx*(j+ny*k);
                    if(!candidates[id].load(std::memory_order_relaxed) ||
                       atomicHits[id].load(std::memory_order_relaxed)) continue;
                    const Point c=cellCenter(domain,gridSpacing,i,j,k);
                    if(triangleIntersectsCell(c,gridSpacing,f.a,f.b,f.c))
                        atomicHits[id].store(1,std::memory_order_relaxed);
                }
    }

#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index=0; index<static_cast<std::ptrdiff_t>(cellCount); ++index)
    {
        const std::size_t id=static_cast<std::size_t>(index);
        if(candidates[id].load(std::memory_order_relaxed) &&
           atomicHits[id].load(std::memory_order_relaxed))
            cellTypes[id]=CellType::Boundary;
    }

    // Output: preserve the existing SoA boundary-coordinate and scalar API.
    boundaryX.resize(cellCount);
    boundaryY.resize(cellCount);
    boundaryZ.resize(cellCount);
    std::size_t boundaryCount=0;
#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index=0; index<static_cast<std::ptrdiff_t>(cellCount); ++index)
    {
        const std::size_t id=static_cast<std::size_t>(index);
        if(cellTypes[id]!=CellType::Boundary) continue;
        const std::size_t k=id/xySize, j=(id%xySize)/nx, i=id%nx;
        std::size_t slot=0;
#pragma omp atomic capture
        slot=boundaryCount++;
        const Point c=cellCenter(domain,gridSpacing,i,j,k);
        boundaryX[slot]=c[0]; boundaryY[slot]=c[1]; boundaryZ[slot]=c[2];
    }
    boundaryX.resize(boundaryCount);
    boundaryY.resize(boundaryCount);
    boundaryZ.resize(boundaryCount);
    scalar.resize(cellCount);
#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index=0; index<static_cast<std::ptrdiff_t>(cellCount); ++index)
    {
        const std::size_t id=static_cast<std::size_t>(index);
        scalar[id]=static_cast<double>(cellTypes[id]);
    }
}

// Explicit instantiations for every supported GeometryAnalysis3D lattice.
template void STLGeometry::interiorAreaAnalysis<lattice::LatticeType::D3Q15>(
    GeometryAnalysis3D<lattice::LatticeType::D3Q15>&) const;
template void STLGeometry::interiorAreaAnalysis<lattice::LatticeType::D3Q19>(
    GeometryAnalysis3D<lattice::LatticeType::D3Q19>&) const;
template void STLGeometry::interiorAreaAnalysis<lattice::LatticeType::D3Q27>(
    GeometryAnalysis3D<lattice::LatticeType::D3Q27>&) const;

} // namespace ntic::lbm::geometry
