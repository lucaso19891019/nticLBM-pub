#pragma once

#include "box.hpp"
#include "circle.hpp"
#include "cylinder.hpp"
#include "rectangle.hpp"
#include "sphere.hpp"
#include "stl_geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ntic::lbm::geometry
{
template<typename LeftGeometry,typename RightGeometry>
struct BooleanGeometry;

namespace detail
{

enum class PointLocation { Dry, Surface, Wet };

struct LinkIntersection
{
    double q;
    std::size_t face;
};

inline PointLocation complementLocation(PointLocation location,bool complemented)
{
    if(!complemented || location==PointLocation::Surface) return location;
    return location==PointLocation::Wet ? PointLocation::Dry : PointLocation::Wet;
}

inline PointLocation flowLocation(bool inside,bool surface,FlowType flow)
{
    if(surface) return PointLocation::Surface;
    return (flow==FlowType::Internal ? inside : !inside)
        ? PointLocation::Wet : PointLocation::Dry;
}

inline std::size_t geometricFaceCount(const Rectangle&) { return 4; }
inline std::size_t geometricFaceCount(const Circle&) { return 1; }
inline std::size_t geometricFaceCount(const Box&) { return 6; }
inline std::size_t geometricFaceCount(const Sphere&) { return 1; }
inline std::size_t geometricFaceCount(const Cylinder&) { return 3; }
inline std::size_t geometricFaceCount(const STLGeometry& geometry)
{
    if(geometry.flowType==FlowType::External) return 1;
    if(geometry.boundaryFeatures.nPlanarFaces>=std::numeric_limits<std::size_t>::max()-1)
        throw std::overflow_error("STL Face count overflow.");
    if(geometry.boundaryFeatures.facetFaceIDs.size()!=geometry.topology.geometry.facetVertexIDs.size())
        throw std::invalid_argument("STL Face mapping is unavailable.");
    return geometry.boundaryFeatures.nPlanarFaces+1;
}

template<typename Geometry,typename Position>
PointLocation locateBox(const Geometry& geometry,const Position& p)
{
    bool inside=true,closed=true,surface=false;
    for(std::size_t axis=0;axis<p.size();++axis)
    {
        inside=inside && p[axis]>geometry.min[axis] && p[axis]<geometry.max[axis];
        closed=closed && p[axis]>=geometry.min[axis] && p[axis]<=geometry.max[axis];
        surface=surface || p[axis]==geometry.min[axis] || p[axis]==geometry.max[axis];
    }
    return flowLocation(inside,closed && surface,geometry.flowType);
}

inline PointLocation locatePoint(const Rectangle& geometry,const Point2D& p)
{ return locateBox(geometry,p); }
inline PointLocation locatePoint(const Box& geometry,const Point& p)
{ return locateBox(geometry,p); }

template<typename Geometry,typename Position>
PointLocation locateRound(const Geometry& geometry,const Position& p)
{
    double distance2=0.0;
    for(std::size_t axis=0;axis<p.size();++axis)
    {
        const double delta=p[axis]-geometry.center[axis];
        distance2+=delta*delta;
    }
    return flowLocation(distance2<geometry.radius*geometry.radius,
                        std::sqrt(distance2)==geometry.radius,geometry.flowType);
}

inline PointLocation locatePoint(const Circle& geometry,const Point2D& p)
{ return locateRound(geometry,p); }
inline PointLocation locatePoint(const Sphere& geometry,const Point& p)
{ return locateRound(geometry,p); }

inline std::array<std::size_t,3> cylinderDimensions(const Cylinder& geometry)
{
    if(geometry.axis==Axis::X) return {0,1,2};
    if(geometry.axis==Axis::Y) return {1,0,2};
    if(geometry.axis==Axis::Z) return {2,0,1};
    throw std::invalid_argument("Unsupported Cylinder axis.");
}

inline PointLocation locatePoint(const Cylinder& geometry,const Point& p)
{
    const auto axes=cylinderDimensions(geometry);
    const double axial=std::abs(p[axes[0]]-geometry.center[axes[0]]);
    const double u=p[axes[1]]-geometry.center[axes[1]];
    const double v=p[axes[2]]-geometry.center[axes[2]];
    const double radial2=u*u+v*v,half=0.5*geometry.length;
    const bool inside=axial<half && radial2<geometry.radius*geometry.radius;
    const bool surface=(axial==half && radial2<=geometry.radius*geometry.radius) ||
        (axial<=half && std::sqrt(radial2)==geometry.radius);
    return flowLocation(inside,surface,geometry.flowType);
}

// STL definitions share the existing triangle predicates in stl_area_analysis.cpp.
bool pointInClosedComponent(const stl::FacetTopology& topology,
                            std::size_t componentID, const Point& p);
bool fluidContainsPoint(const STLGeometry& geometry, const Point& p);
PointLocation locatePoint(const STLGeometry& geometry,const Point& p);
void appendIntersections(const STLGeometry& geometry,const Point& p,const Point& delta,
                         std::size_t faceBase,std::vector<LinkIntersection>& output);

inline void appendHit(double q,std::size_t face,std::vector<LinkIntersection>& output)
{
    if(std::isfinite(q) && q>0.0 && q<=1.0) output.push_back({q,face});
}

template<typename Geometry,typename Position>
void appendBoxIntersections(const Geometry& geometry,const Position& p,const Position& delta,
    const std::vector<std::size_t>& faceIDs,std::size_t faceBase,
    std::vector<LinkIntersection>& output)
{
    for(std::size_t axis=0;axis<p.size();++axis)
    {
        if(delta[axis]==0.0) continue;
        for(std::size_t side=0;side<2;++side)
        {
            const double plane=side==0 ? geometry.min[axis] : geometry.max[axis];
            const long double parameter=(static_cast<long double>(plane)-p[axis])/delta[axis];
            if(!(parameter>0.0L && parameter<=1.0L)) continue;
            bool finite=true;
            for(std::size_t other=0;other<p.size();++other)
            {
                if(other==axis) continue;
                const long double position=static_cast<long double>(p[other])+
                    parameter*delta[other];
                finite=finite && position>=geometry.min[other] && position<=geometry.max[other];
            }
            if(finite) appendHit(static_cast<double>(parameter),faceBase+faceIDs[2*axis+side],output);
        }
    }
}

inline void appendIntersections(const Rectangle& geometry,const Point2D& p,const Point2D& delta,
    std::size_t base,std::vector<LinkIntersection>& output)
{ appendBoxIntersections(geometry,p,delta,{3,1,0,2},base,output); }

inline void appendIntersections(const Box& geometry,const Point& p,const Point& delta,
    std::size_t base,std::vector<LinkIntersection>& output)
{ appendBoxIntersections(geometry,p,delta,{0,1,2,3,4,5},base,output); }

template<typename Consume>
void quadraticRoots(double a,double b,double c,Consume consume)
{
    if(!(a>0.0)) return;
    const double discriminant=b*b-4.0*a*c;
    if(discriminant<0.0) return;
    const double term=-0.5*(b+std::copysign(std::sqrt(discriminant),b));
    if(term==0.0) consume(-b/(2.0*a));
    else
    {
        consume(term/a);
        consume(c/term);
    }
}

template<typename Geometry,typename Position>
void appendRoundIntersections(const Geometry& geometry,const Position& p,const Position& delta,
    std::size_t base,std::vector<LinkIntersection>& output)
{
    double a=0.0,b=0.0,c=-geometry.radius*geometry.radius;
    for(std::size_t axis=0;axis<p.size();++axis)
    {
        const double relative=p[axis]-geometry.center[axis];
        a+=delta[axis]*delta[axis];
        b+=2.0*relative*delta[axis];
        c+=relative*relative;
    }
    quadraticRoots(a,b,c,[&](double q){appendHit(q,base,output);});
}

inline void appendIntersections(const Circle& geometry,const Point2D& p,const Point2D& delta,
    std::size_t base,std::vector<LinkIntersection>& output)
{ appendRoundIntersections(geometry,p,delta,base,output); }
inline void appendIntersections(const Sphere& geometry,const Point& p,const Point& delta,
    std::size_t base,std::vector<LinkIntersection>& output)
{ appendRoundIntersections(geometry,p,delta,base,output); }

inline void appendIntersections(const Cylinder& geometry,const Point& p,const Point& delta,
    std::size_t base,std::vector<LinkIntersection>& output)
{
    const auto axes=cylinderDimensions(geometry);
    const double a=p[axes[0]]-geometry.center[axes[0]],da=delta[axes[0]];
    const double u=p[axes[1]]-geometry.center[axes[1]],du=delta[axes[1]];
    const double v=p[axes[2]]-geometry.center[axes[2]],dv=delta[axes[2]];
    const double half=0.5*geometry.length,r2=geometry.radius*geometry.radius;
    quadraticRoots(du*du+dv*dv,2.0*(u*du+v*dv),u*u+v*v-r2,[&](double q) {
        if(q>0.0 && q<=1.0 && std::abs(a+q*da)<=half) appendHit(q,base,output);
    });
    if(da==0.0) return;
    for(int sign:{-1,1})
    {
        const double q=(sign*half-a)/da;
        const double ru=u+q*du,rv=v+q*dv;
        if(ru*ru+rv*rv<=r2) appendHit(q,base+(sign<0 ? 1 : 2),output);
    }
}

template<typename Geometry>
struct GeometryQuery
{
    static std::size_t faceCount(const Geometry& geometry)
    { return geometricFaceCount(geometry); }

    template<typename Position>
    static PointLocation locate(const Geometry& geometry,const Position& p,bool complemented,
        const std::vector<std::size_t>& forcedFaces,std::size_t base)
    {
        const std::size_t end=base+faceCount(geometry);
        for(std::size_t face:forcedFaces)
            if(face>=base && face<end) return PointLocation::Surface;
        return complementLocation(locatePoint(geometry,p),complemented);
    }

    template<typename Position>
    static void intersections(const Geometry& geometry,const Position& p,const Position& delta,
        std::size_t base,std::vector<LinkIntersection>& output)
    { appendIntersections(geometry,p,delta,base,output); }

    template<typename Position>
    static void surfaceFaces(const Geometry& geometry,const Position& p,bool complemented,
        const std::vector<std::size_t>& forcedFaces,std::size_t base,
        std::vector<std::size_t>& output)
    {
        if(locate(geometry,p,complemented,forcedFaces,base)!=PointLocation::Surface) return;
        const std::size_t end=base+faceCount(geometry);
        for(std::size_t face:forcedFaces)
            if(face>=base && face<end) output.push_back(face);
    }

    template<typename Analysis,typename Domain>
    static void interior(const Geometry& geometry,Analysis& analysis,const Domain& domain,bool complemented)
    {
        analysis.domainBounds=domain;
        analysis.complemented=complemented;
        geometry.interiorAreaAnalysis(analysis);
    }

    template<typename Analysis>
    static void boundary(const Geometry& geometry,Analysis& analysis,bool)
    { geometry.boundaryQAnalysis(analysis); }

    template<typename Analysis,typename Domain>
    static void analyze(const Geometry& geometry,Analysis& analysis,const Domain& domain,bool complemented)
    {
        interior(geometry,analysis,domain,complemented);
        boundary(geometry,analysis,complemented);
    }
};

template<typename LeftGeometry,typename RightGeometry>
struct GeometryQuery<BooleanGeometry<LeftGeometry,RightGeometry>>;

} // namespace detail
} // namespace ntic::lbm::geometry
