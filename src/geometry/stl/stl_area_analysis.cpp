#include "stl_geometry.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
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

// Zero-geometric-tolerance predicate for Step 3.2.
// This uses double arithmetic (not an adaptive exact predicate): near-degenerate
// floating-point cases can still be affected by rounding.
bool pointOnTriangle(const Point& p,
                     const Point& a,
                     const Point& b,
                     const Point& c)
{
    const Point ab=subtract(b,a), ac=subtract(c,a), ap=subtract(p,a);
    const Point n=cross(ab,ac);
    const double n2=dot(n,n);
    if(n2==0.0) return false;

    if(dot(ap,n) != 0.0) return false;

    const double d00=dot(ab,ab), d01=dot(ab,ac), d11=dot(ac,ac);
    const double d20=dot(ap,ab), d21=dot(ap,ac);
    const double denominator=d00*d11-d01*d01;
    if(denominator<=0.0) return false;

    const double u=(d11*d20-d01*d21)/denominator;
    const double v=(d00*d21-d01*d20)/denominator;
    return u>=0.0 && v>=0.0 && u+v<=1.0;
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
    analysis.nBoundaries=0;
    analysis.boundaryOffsets.assign(1,0);
    analysis.q.clear();

    // Step 1: Facet-driven closed-voxel rasterization. Touch is Boundary.
    // The surface-center flag is evaluated without geometric tolerance,
    // because it is consumed by Step 3.2.
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
    // center. onSurface was produced by a zero-tolerance predicate in Step 1.
    // Points exactly on any active facet are always Dry.
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

namespace
{

// Two-sided, closed-edge segment/triangle intersection. Projecting along the
// link's largest component gives shared triangle edges the same arithmetic.
// Long double reduces cancellation; no approximate-distance tie tolerance is
// introduced. Coplanar links are clipped against the finite triangle edges.
bool linkTriangle(const Point& p, const Point& d, const FacetData& f, double& q)
{
    using WidePoint=std::array<long double,3>;
    const auto difference=[](const Point& a,const Point& b) {
        return WidePoint{static_cast<long double>(a[0])-b[0],
                         static_cast<long double>(a[1])-b[1],
                         static_cast<long double>(a[2])-b[2]};
    };
    const auto product=[](const WidePoint& a,const WidePoint& b) {
        return WidePoint{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],
                         a[0]*b[1]-a[1]*b[0]};
    };
    const auto scalar=[](const WidePoint& a,const WidePoint& b) {
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    };
    const WidePoint a=difference(f.a,p),b=difference(f.b,p),c=difference(f.c,p);
    const WidePoint direction{d[0],d[1],d[2]};
    const WidePoint normal=product(difference(f.b,f.a),difference(f.c,f.a));
    if(scalar(normal,normal)==0.0L) return false;
    const long double denominator=scalar(normal,direction);
    long double hit=-1.0L;
    if(denominator==0.0L)
    {
        if(scalar(normal,a)!=0.0L) return false;
        long double enter=0.0L,leave=1.0L;
        const std::array<Point,3> vertices{{f.a,f.b,f.c}};
        for(std::size_t edge=0;edge<3;++edge)
        {
            const auto e=difference(vertices[(edge+1)%3],vertices[edge]);
            const auto relative=difference(p,vertices[edge]);
            const long double value=scalar(product(e,relative),normal);
            const long double slope=scalar(product(e,direction),normal);
            if(slope==0.0L)
            {
                if(value<0.0L) return false;
            }
            else if(slope>0.0L) enter=std::max(enter,-value/slope);
            else leave=std::min(leave,-value/slope);
            if(enter>leave) return false;
        }
        hit=enter;
    }
    else
    {
        std::size_t z=0;
        for(std::size_t axis=1;axis<3;++axis)
            if(std::abs(d[axis])>std::abs(d[z])) z=axis;
        const std::size_t x=(z+1)%3,y=(x+1)%3;
        const long double sx=direction[x]/direction[z],sy=direction[y]/direction[z];
        const long double ax=a[x]-sx*a[z],ay=a[y]-sy*a[z];
        const long double bx=b[x]-sx*b[z],by=b[y]-sy*b[z];
        const long double cx=c[x]-sx*c[z],cy=c[y]-sy*c[z];
        const long double u=cx*by-cy*bx,v=ax*cy-ay*cx,w=bx*ay-by*ax;
        if((u<0.0L || v<0.0L || w<0.0L) &&
           (u>0.0L || v>0.0L || w>0.0L)) return false;
        const long double determinant=u+v+w;
        if(determinant==0.0L) return false;
        hit=(u*a[z]+v*b[z]+w*c[z])/(determinant*direction[z]);
    }
    if(!(hit>0.0L && hit<=1.0L)) return false;
    q=static_cast<double>(hit);
    return std::isfinite(q) && q>0.0 && q<=1.0;
}

template<std::size_t Links>
struct BoundaryEntry
{
    std::size_t face,cell;
    std::array<double,Links> q;
};

} // namespace

template<lattice::LatticeType Type>
void STLGeometry::boundaryQAnalysis(GeometryAnalysis3D<Type>& analysis) const
{
    using Lattice=typename GeometryAnalysis3D<Type>::Lattice;
    constexpr std::size_t links=Lattice::nStencils-1;
    const double h=analysis.gridSpacing;
    const std::size_t nx=analysis.nx,ny=analysis.ny,nz=analysis.nz;
    if(!(h>0.0) || !std::isfinite(h) || nx==0 || ny==0 || nz==0 ||
       nx>std::numeric_limits<std::size_t>::max()/ny ||
       nx*ny>std::numeric_limits<std::size_t>::max()/nz)
        throw std::invalid_argument("STL Q requires a valid interior analysis grid.");
    const std::size_t cellCount=nx*ny*nz;
    const std::size_t facetCount=topology.geometry.facetVertexIDs.size();
    const bool internal=flowType==FlowType::Internal;
    if(analysis.scalar.size()!=cellCount ||
       analysis.boundaryY.size()!=analysis.boundaryX.size() ||
       analysis.boundaryZ.size()!=analysis.boundaryX.size() ||
       flow.size()!=topology.components.size() ||
       topology.facetComponentIDs.size()!=facetCount)
        throw std::invalid_argument("Incomplete STL interior or component data.");
    if(internal && boundaryFeatures.facetFaceIDs.size()!=facetCount)
        throw std::invalid_argument("STL planar Face mapping is unavailable.");
    if(internal && boundaryFeatures.nPlanarFaces>=std::numeric_limits<std::size_t>::max()-1)
        throw std::overflow_error("STL Face count overflow.");
    const std::size_t faceCount=internal ? boundaryFeatures.nPlanarFaces+1 : 1;

    // Consume interior coordinates, canonicalizing repeated calls on already
    // Face-expanded output back to unique physical cells. Do not reclassify.
    std::vector<std::size_t> cells;
    cells.reserve(analysis.boundaryX.size());
    for(std::size_t b=0;b<analysis.boundaryX.size();++b)
    {
        const Point p{analysis.boundaryX[b],analysis.boundaryY[b],analysis.boundaryZ[b]};
        std::array<std::size_t,3> index{};
        const std::array<std::size_t,3> dimensions{{nx,ny,nz}};
        for(std::size_t axis=0;axis<3;++axis)
        {
            const double value=(p[axis]-analysis.domain.min[axis])/h-0.5;
            const double rounded=std::round(value);
            if(!std::isfinite(value) || rounded<0.0 || rounded>=dimensions[axis] ||
               std::abs(value-rounded)>64.0*std::numeric_limits<double>::epsilon()*
                                           std::max(1.0,std::abs(value)))
                throw std::invalid_argument("STL boundary coordinate is not a grid center.");
            index[axis]=static_cast<std::size_t>(rounded);
        }
        const std::size_t id=index[0]+nx*(index[1]+ny*index[2]);
        if(analysis.scalar[id]!=static_cast<double>(CellType::Boundary))
            throw std::invalid_argument("STL boundary entry is not a Boundary cell.");
        cells.push_back(id);
    }
    std::sort(cells.begin(),cells.end());
    cells.erase(std::unique(cells.begin(),cells.end()),cells.end());
    if(static_cast<std::size_t>(std::count(analysis.scalar.begin(),analysis.scalar.end(),
                              static_cast<double>(CellType::Boundary)))!=cells.size())
        throw std::invalid_argument("STL boundary coordinates do not cover interior Boundary cells.");

    // Validate malformed topology on the calling thread before entering OpenMP.
    std::vector<unsigned char> active(facetCount,0);
    for(std::size_t facet=0;facet<facetCount;++facet)
    {
        const std::size_t component=topology.facetComponentIDs[facet];
        if(component>=flow.size()) throw std::invalid_argument("Invalid STL facet component.");
        active[facet]=flow[component].active;
        if(internal && boundaryFeatures.facetFaceIDs[facet]>=faceCount)
            throw std::invalid_argument("Invalid STL facet Face ID.");
        for(std::size_t vertex:topology.geometry.facetVertexIDs[facet])
        {
            if(vertex>=topology.geometry.vertices.size())
                throw std::invalid_argument("Invalid STL facet vertex.");
            for(double coordinate:topology.geometry.vertices[vertex])
                if(!std::isfinite(coordinate))
                    throw std::invalid_argument("Non-finite STL facet vertex.");
        }
        if(active[facet])
        {
            const FacetData triangle=facetData(*this,facet);
            const Point n=cross(subtract(triangle.b,triangle.a),subtract(triangle.c,triangle.a));
            if(dot(n,n)==0.0 || !std::isfinite(dot(n,n)))
                throw std::invalid_argument("Degenerate active STL facet.");
        }
    }
    if(cells.empty())
    {
        analysis.nBoundaries=faceCount;
        analysis.boundaryOffsets.assign(faceCount+1,0);
        analysis.boundaryX.clear(); analysis.boundaryY.clear(); analysis.boundaryZ.clear();
        analysis.q.clear();
        return;
    }
    std::vector<FacetData> facets(facetCount);
#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t f=0;f<static_cast<std::ptrdiff_t>(facetCount);++f)
        if(active[static_cast<std::size_t>(f)])
            facets[static_cast<std::size_t>(f)]=facetData(*this,static_cast<std::size_t>(f));

    // Reuse interior's facet-driven range and closed half-width-h SAT. Only
    // final Boundary cells are visited. Unlike classification, keep ALL hits.
    using Candidate=std::pair<std::size_t,std::size_t>;
    const std::size_t threads=static_cast<std::size_t>(omp_get_max_threads());
    std::vector<std::vector<Candidate>> localCandidates(threads);
#pragma omp parallel
    {
        auto& local=localCandidates[static_cast<std::size_t>(omp_get_thread_num())];
#pragma omp for schedule(dynamic,32)
        for(std::ptrdiff_t f=0;f<static_cast<std::ptrdiff_t>(facetCount);++f)
        {
            const std::size_t facet=static_cast<std::size_t>(f);
            if(!active[facet]) continue;
            const auto& triangle=facets[facet];
            std::array<std::size_t,3> lo{},hi{};
            if(!facetRange(triangle,analysis.domain,h,h,nx,ny,nz,lo,hi)) continue;
            for(std::size_t k=lo[2];k<=hi[2];++k)
                for(std::size_t j=lo[1];j<=hi[1];++j)
                    for(std::size_t i=lo[0];i<=hi[0];++i)
                    {
                        const std::size_t id=i+nx*(j+ny*k);
                        if(analysis.scalar[id]!=static_cast<double>(CellType::Boundary)) continue;
                        const Point p=cellCenter(analysis.domain,h,i,j,k);
                        if(triangleIntersectsCell(p,h,triangle.a,triangle.b,triangle.c))
                            local.emplace_back(id,facet);
                    }
        }
    }
    std::vector<Candidate> candidates;
    std::size_t candidateCount=0;
    for(const auto& local:localCandidates)
    {
        if(local.size()>std::numeric_limits<std::size_t>::max()-candidateCount)
            throw std::overflow_error("STL candidate count overflow.");
        candidateCount+=local.size();
    }
    candidates.reserve(candidateCount);
    for(auto& local:localCandidates)
    {
        candidates.insert(candidates.end(),local.begin(),local.end());
        std::vector<Candidate>().swap(local);
    }
    std::sort(candidates.begin(),candidates.end());
    candidates.erase(std::unique(candidates.begin(),candidates.end()),candidates.end());
    std::vector<std::size_t> offsets(cells.size()+1,0);
    std::size_t next=0;
    for(std::size_t b=0;b<cells.size();++b)
    {
        offsets[b]=next;
        while(next<candidates.size() && candidates[next].first==cells[b]) ++next;
        if(offsets[b]==next)
            throw std::runtime_error("STL Boundary cell has no active stencil facet.");
    }
    offsets.back()=next;

    // Each worker owns every link of its cell: no atomic/critical Q updates.
    using Entry=BoundaryEntry<links>;
    std::vector<std::vector<Entry>> localEntries(threads);
#pragma omp parallel
    {
        auto& output=localEntries[static_cast<std::size_t>(omp_get_thread_num())];
#pragma omp for schedule(dynamic,32)
        for(std::ptrdiff_t index=0;index<static_cast<std::ptrdiff_t>(cells.size());++index)
        {
            const std::size_t b=static_cast<std::size_t>(index),id=cells[b];
            const Point p=cellCenter(analysis.domain,h,id%nx,(id/nx)%ny,id/(nx*ny));
            std::array<double,links> nearest;
            nearest.fill(-1.0);
            std::array<std::size_t,links> owners,winningFacets;
            owners.fill(faceCount);
            winningFacets.fill(facetCount);
            std::vector<std::size_t> faces;
            for(std::size_t c=offsets[b];c<offsets[b+1];++c)
            {
                const std::size_t facet=candidates[c].second;
                const std::size_t face=internal ? boundaryFeatures.facetFaceIDs[facet] : 0;
                faces.push_back(face);
                for(std::size_t link=1;link<Lattice::nStencils;++link)
                {
                    const Point d{h*Lattice::ex[link],h*Lattice::ey[link],h*Lattice::ez[link]};
                    double q;
                    if(!linkTriangle(p,d,facets[facet],q)) continue;
                    const std::size_t slot=link-1;
                    if(nearest[slot]<0.0 || q<nearest[slot] ||
                       (q==nearest[slot] && (face<owners[slot] ||
                        (face==owners[slot] && facet<winningFacets[slot]))))
                    {
                        nearest[slot]=q;
                        owners[slot]=face;
                        winningFacets[slot]=facet;
                    }
                }
            }
            std::sort(faces.begin(),faces.end());
            faces.erase(std::unique(faces.begin(),faces.end()),faces.end());
            for(std::size_t face:faces)
            {
                Entry entry{face,id,{}};
                entry.q.fill(-1.0);
                for(std::size_t link=0;link<links;++link)
                    if(owners[link]==face) entry.q[link]=nearest[link];
                output.push_back(entry);
            }
        }
    }
    std::vector<Entry> entries;
    for(const auto& local:localEntries)
        entries.insert(entries.end(),local.begin(),local.end());
    std::sort(entries.begin(),entries.end(),[](const Entry& a,const Entry& b) {
        return a.face<b.face || (a.face==b.face && a.cell<b.cell);
    });
    if(entries.size()>std::numeric_limits<std::size_t>::max()/links)
        throw std::overflow_error("STL Q output size overflow.");
    std::vector<double> x(entries.size()),y(entries.size()),z(entries.size());
    std::vector<double> q(entries.size()*links,-1.0);
    std::vector<std::size_t> boundaryOffsets(faceCount+1,0);
    for(const auto& entry:entries) ++boundaryOffsets[entry.face+1];
    for(std::size_t face=0;face<faceCount;++face)
        boundaryOffsets[face+1]+=boundaryOffsets[face];
#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index=0;index<static_cast<std::ptrdiff_t>(entries.size());++index)
    {
        const std::size_t b=static_cast<std::size_t>(index),id=entries[b].cell;
        const Point p=cellCenter(analysis.domain,h,id%nx,(id/nx)%ny,id/(nx*ny));
        x[b]=p[0]; y[b]=p[1]; z[b]=p[2];
        std::copy(entries[b].q.begin(),entries[b].q.end(),q.begin()+b*links);
    }
    analysis.nBoundaries=faceCount;
    analysis.boundaryOffsets.swap(boundaryOffsets);
    analysis.boundaryX.swap(x); analysis.boundaryY.swap(y); analysis.boundaryZ.swap(z);
    analysis.q.swap(q);
}

// Explicit instantiations for every supported GeometryAnalysis3D lattice.
template void STLGeometry::interiorAreaAnalysis<lattice::LatticeType::D3Q15>(
    GeometryAnalysis3D<lattice::LatticeType::D3Q15>&) const;
template void STLGeometry::interiorAreaAnalysis<lattice::LatticeType::D3Q19>(
    GeometryAnalysis3D<lattice::LatticeType::D3Q19>&) const;
template void STLGeometry::interiorAreaAnalysis<lattice::LatticeType::D3Q27>(
    GeometryAnalysis3D<lattice::LatticeType::D3Q27>&) const;

template void STLGeometry::boundaryQAnalysis<lattice::LatticeType::D3Q15>(
    GeometryAnalysis3D<lattice::LatticeType::D3Q15>&) const;
template void STLGeometry::boundaryQAnalysis<lattice::LatticeType::D3Q19>(
    GeometryAnalysis3D<lattice::LatticeType::D3Q19>&) const;
template void STLGeometry::boundaryQAnalysis<lattice::LatticeType::D3Q27>(
    GeometryAnalysis3D<lattice::LatticeType::D3Q27>&) const;

} // namespace ntic::lbm::geometry
