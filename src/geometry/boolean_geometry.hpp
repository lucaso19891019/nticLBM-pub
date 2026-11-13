#pragma once

#include "flow_type.hpp"
#include "geometry_query.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <omp.h>


namespace ntic::lbm::geometry
{

enum class BooleanOperation
{
    Union,
    Intersection,
    Difference
};


template <
    typename LeftGeometry,
    typename RightGeometry>
struct BooleanGeometry
{
    using BoundingBoxType =
        std::decay_t<
            decltype(
                std::declval<LeftGeometry>().bounds)>;

    using RightBoundingBoxType =
        std::decay_t<
            decltype(
                std::declval<RightGeometry>().bounds)>;

    static_assert(
        std::is_same_v<
            BoundingBoxType,
            RightBoundingBoxType>,
        "Boolean geometry operands must have "
        "the same bounding-box type.");


    LeftGeometry left;
    RightGeometry right;

    BooleanOperation operation;
    FlowType flowType;

    BoundingBoxType bounds;
    BoundingBoxType openBox;


    BooleanGeometry(
        LeftGeometry leftGeometry,
        RightGeometry rightGeometry,
        const BooleanOperation booleanOperation)
        :
        left(
            std::move(leftGeometry)),
        right(
            std::move(rightGeometry)),
        operation(
            booleanOperation),
        flowType(
            determineFlowType(
                left.flowType,
                right.flowType,
                operation)),
        bounds(
            determineBounds(
                left,
                right,
                operation)),
        openBox(
            determineOpenBox(
                left,
                right,
                operation))
    {
    }

    template<lattice::LatticeType Type>
    void analysis(GeometryAnalysis2D<Type>& result) const
    {
        static_assert(std::is_same_v<BoundingBoxType,BoundingBox2D>,
                      "A 2D analysis requires 2D Boolean operands.");
        build(result,result.domainBounds.value_or(flowType==FlowType::Internal ? bounds : openBox),
              result.complemented,true);
    }

    template<lattice::LatticeType Type>
    void analysis(GeometryAnalysis3D<Type>& result) const
    {
        static_assert(std::is_same_v<BoundingBoxType,BoundingBox>,
                      "A 3D analysis requires 3D Boolean operands.");
        build(result,result.domainBounds.value_or(flowType==FlowType::Internal ? bounds : openBox),
              result.complemented,true);
    }

    template<lattice::LatticeType Type>
    void interiorAreaAnalysis(GeometryAnalysis2D<Type>& result) const
    {
        static_assert(std::is_same_v<BoundingBoxType,BoundingBox2D>,
                      "A 2D analysis requires 2D Boolean operands.");
        build(result,result.domainBounds.value_or(flowType==FlowType::Internal ? bounds : openBox),
              result.complemented,false);
    }

    template<lattice::LatticeType Type>
    void interiorAreaAnalysis(GeometryAnalysis3D<Type>& result) const
    {
        static_assert(std::is_same_v<BoundingBoxType,BoundingBox>,
                      "A 3D analysis requires 3D Boolean operands.");
        build(result,result.domainBounds.value_or(flowType==FlowType::Internal ? bounds : openBox),
              result.complemented,false);
    }

    template<lattice::LatticeType Type>
    void boundaryQAnalysis(GeometryAnalysis2D<Type>& result) const
    { completeBoundary(result,result.complemented); }

    template<lattice::LatticeType Type>
    void boundaryQAnalysis(GeometryAnalysis3D<Type>& result) const
    { completeBoundary(result,result.complemented); }

private:

    template<typename Geometry>
    friend struct detail::GeometryQuery;

    using LeftQuery=detail::GeometryQuery<LeftGeometry>;
    using RightQuery=detail::GeometryQuery<RightGeometry>;
    using Location=detail::PointLocation;
    using Intersection=detail::LinkIntersection;
    static constexpr std::size_t dimension=
        std::tuple_size<decltype(std::declval<BoundingBoxType>().min)>::value;

    std::size_t geometricFaceCount() const
    {
        const std::size_t countLeft=LeftQuery::faceCount(left);
        const std::size_t countRight=RightQuery::faceCount(right);
        if(countRight>=std::numeric_limits<std::size_t>::max()-countLeft)
            throw std::overflow_error("Boolean Face count overflow.");
        return countLeft+countRight;
    }

    // Difference is intersection with the entire right operand's complement.
    // Complementing a Boolean recursively applies De Morgan's laws; no saved
    // scalar values are numerically inverted and no member flowType is changed.
    bool isUnion(bool complemented) const
    { return (operation==BooleanOperation::Union) != complemented; }

    bool rightComplement(bool complemented) const
    { return complemented != (operation==BooleanOperation::Difference); }

    template<typename Position>
    Location queryPoint(const Position& p,bool complemented,
        const std::vector<std::size_t>& forcedFaces,std::size_t base) const
    {
        const Location a=LeftQuery::locate(left,p,complemented,forcedFaces,base);
        const Location b=RightQuery::locate(right,p,rightComplement(complemented),
                                           forcedFaces,base+LeftQuery::faceCount(left));
        if(isUnion(complemented))
        {
            if(a==Location::Wet || b==Location::Wet) return Location::Wet;
            if(a==Location::Surface || b==Location::Surface) return Location::Surface;
            return Location::Dry;
        }
        if(a==Location::Dry || b==Location::Dry) return Location::Dry;
        if(a==Location::Surface || b==Location::Surface) return Location::Surface;
        return Location::Wet;
    }

    template<typename Position>
    void queryLinks(const Position& p,const Position& delta,std::size_t base,
                    std::vector<Intersection>& output) const
    {
        LeftQuery::intersections(left,p,delta,base,output);
        RightQuery::intersections(right,p,delta,base+LeftQuery::faceCount(left),output);
    }

    template<typename Position>
    void querySurfaceFaces(const Position& p,bool complemented,
        const std::vector<std::size_t>& forcedFaces,std::size_t base,
        std::vector<std::size_t>& output) const
    {
        const std::size_t rightBase=base+LeftQuery::faceCount(left);
        const Location a=LeftQuery::locate(left,p,complemented,forcedFaces,base);
        const Location b=RightQuery::locate(right,p,rightComplement(complemented),forcedFaces,rightBase);
        if(isUnion(complemented))
        {
            if(a==Location::Wet || b==Location::Wet) return;
        }
        else if(a==Location::Dry || b==Location::Dry) return;
        if(a==Location::Surface)
            LeftQuery::surfaceFaces(left,p,complemented,forcedFaces,base,output);
        if(b==Location::Surface)
            RightQuery::surfaceFaces(right,p,rightComplement(complemented),forcedFaces,rightBase,output);
    }

    static double mergeState(double a,double b,bool unite)
    {
        const double dry=static_cast<double>(CellType::Dry);
        const double boundary=static_cast<double>(CellType::Boundary);
        const double interior=static_cast<double>(CellType::Interior);
        if(unite)
        {
            if(a==interior || b==interior) return interior;
            if(a==boundary || b==boundary) return boundary;
            return dry;
        }
        if(a==dry || b==dry) return dry;
        if(a==interior && b==interior) return interior;
        return boundary;
    }

    static std::size_t validateDomain(const BoundingBoxType& domain,double h)
    {
        if(!(h>0.0) || !std::isfinite(h))
            throw std::invalid_argument("Grid spacing must be finite and positive.");
        std::size_t count=1;
        for(std::size_t axis=0;axis<dimension;++axis)
        {
            const double width=domain.max[axis]-domain.min[axis];
            const double cells=std::ceil(width/h);
            if(!std::isfinite(domain.min[axis]) || !std::isfinite(domain.max[axis]) ||
               !(width>0.0) || !std::isfinite(cells) || cells<1.0 ||
               cells>=static_cast<double>(std::numeric_limits<std::ptrdiff_t>::max()))
                throw std::invalid_argument("Invalid Boolean analysis domain.");
            const std::size_t size=static_cast<std::size_t>(cells);
            if(count>static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max())/size)
                throw std::overflow_error("Boolean grid size overflow.");
            count*=size;
        }
        return count;
    }

    template<typename Analysis>
    static std::size_t entryCell(const Analysis& data,std::size_t b)
    {
        const auto id=data.boundaryCells.at(b);
        if(id>=data.scalar.size()) throw std::logic_error("Invalid child Boundary index.");
        return id;
    }

    template<typename Analysis>
    static auto cellCenter(const Analysis& data,std::size_t id)
    {
        auto p=data.domain.min;
        p[0]=data.centerCoordinate(0,id%data.nx);
        p[1]=data.centerCoordinate(1,(id/data.nx)%data.ny);
        if constexpr(dimension==3)
            p[2]=data.centerCoordinate(2,id/(data.nx*data.ny));
        return p;
    }

    template<typename Analysis>
    void completeBoundary(Analysis& result,bool complemented) const
    {
        result.validateBoundaryCells();
        Analysis complete(result.gridSpacing);
        complete.domainBounds=result.domainBounds;
        complete.gridFrame=result.gridFrame;
        complete.complemented=complemented;
        build(complete,result.domainBounds.value_or(flowType==FlowType::Internal ? bounds : openBox),
              complemented,true);
        bool valid=result.nx==complete.nx && result.ny==complete.ny &&
                   result.scalar==complete.scalar && result.domain.min==complete.domain.min &&
                   result.domain.max==complete.domain.max;
        if constexpr(dimension==3) valid=valid && result.nz==complete.nz;
        if(!valid) throw std::invalid_argument("Boolean Q requires its matching interior analysis.");
        result=std::move(complete);
    }

    template<typename Analysis>
    void build(Analysis& result,const BoundingBoxType& domain,bool complemented,bool withQ) const
    {
        using Lattice=typename Analysis::Lattice;
        constexpr std::size_t links=Lattice::nStencils-1;
        const std::size_t count=validateDomain(domain,result.gridSpacing);
        const std::size_t leftFaces=LeftQuery::faceCount(left),faces=geometricFaceCount();
        Analysis a(result.gridSpacing),b(result.gridSpacing);
        a.gridFrame=result.gridFrame; b.gridFrame=result.gridFrame;
        // Child results live through the whole merge. Nested Boolean children
        // recursively complete their own results using this very same grid.
        // Face memberships are inherited from complete child outputs, including
        // stencil-only entries. A standalone interior call leaves its own Q
        // empty; analysis() performs the complete merge in a single traversal.
        LeftQuery::analyze(left,a,domain,complemented);
        RightQuery::analyze(right,b,domain,rightComplement(complemented));
        if(a.scalar.size()!=count || b.scalar.size()!=count || a.nx!=b.nx || a.ny!=b.ny)
            throw std::logic_error("Boolean operands did not use the common grid.");
        if constexpr(dimension==3)
            if(a.nz!=b.nz) throw std::logic_error("Boolean operand depths differ.");
        const bool unite=isUnion(complemented);
        Analysis merged(result.gridSpacing);
        merged.domainBounds=result.domainBounds;
        merged.gridFrame=result.gridFrame;
        merged.complemented=complemented;
        merged.domain=domain; merged.nx=a.nx; merged.ny=a.ny;
        if constexpr(dimension==3) merged.nz=a.nz;
        merged.validateGridFrame();
        merged.scalar.resize(count);
#pragma omp parallel for schedule(static)
        for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(count);++i)
        {
            const std::size_t id=static_cast<std::size_t>(i);
            merged.scalar[id]=mergeState(a.scalar[id],b.scalar[id],unite);
        }

        struct SourceLink { std::size_t cell,link,face; double q; };
        using Membership=std::pair<std::size_t,std::size_t>;
        std::vector<Membership> memberships;
        std::vector<SourceLink> sourceLinks;
        const auto collect=[&](const Analysis& child,std::size_t base,std::size_t expectedFaces) {
            if(child.nBoundaries!=expectedFaces || child.boundaryOffsets.size()!=expectedFaces+1 ||
               child.boundaryOffsets.front()!=0 || child.boundaryOffsets.back()!=child.boundaryX.size() ||
               !std::is_sorted(child.boundaryOffsets.begin(),child.boundaryOffsets.end()) ||
               child.boundaryY.size()!=child.boundaryX.size() ||
               child.q.size()!=child.boundaryX.size()*links)
                throw std::logic_error("Invalid child Face/Q output.");
            if constexpr(dimension==3)
                if(child.boundaryZ.size()!=child.boundaryX.size())
                    throw std::logic_error("Invalid child Z coordinate output.");
            child.validateBoundaryCells();
            for(std::size_t face=0;face<expectedFaces;++face)
                for(std::size_t entry=child.boundaryOffsets[face];entry<child.boundaryOffsets[face+1];++entry)
                {
                    const std::size_t id=entryCell(child,entry);
                    if(merged.scalar[id]!=static_cast<double>(CellType::Boundary)) continue;
                    memberships.emplace_back(id,base+face);
                    if(withQ)
                        for(std::size_t link=0;link<links;++link)
                        {
                            const double q=child.q[entry*links+link];
                            if(q>0.0) sourceLinks.push_back({id,link,base+face,q});
                        }
                }
        };
        collect(a,0,leftFaces);
        collect(b,leftFaces,faces-leftFaces);
        std::sort(memberships.begin(),memberships.end());
        memberships.erase(std::unique(memberships.begin(),memberships.end()),memberships.end());
        std::sort(sourceLinks.begin(),sourceLinks.end(),[](const SourceLink& x,const SourceLink& y) {
            return std::tie(x.cell,x.link,x.face)<std::tie(y.cell,y.link,y.face);
        });
        std::vector<std::size_t> cells,memberOffsets{0},linkOffsets{0};
        std::size_t nextMember=0,nextLink=0;
        for(std::size_t id=0;id<count;++id)
        {
            if(merged.scalar[id]!=static_cast<double>(CellType::Boundary)) continue;
            cells.push_back(id);
            while(nextMember<memberships.size() && memberships[nextMember].first==id) ++nextMember;
            while(nextLink<sourceLinks.size() && sourceLinks[nextLink].cell==id) ++nextLink;
            if(nextMember==memberOffsets.back())
                throw std::logic_error("Boolean Boundary cell has no source Face.");
            memberOffsets.push_back(nextMember); linkOffsets.push_back(nextLink);
        }
        struct Entry { std::size_t cell,face; std::array<double,links> q; };
        std::vector<std::vector<Entry>> perCell(cells.size());
#pragma omp parallel for schedule(dynamic,32)
        for(std::ptrdiff_t index=0;index<static_cast<std::ptrdiff_t>(cells.size());++index)
        {
            const std::size_t cell=static_cast<std::size_t>(index),id=cells[cell];
            std::array<double,links> nearest;
            std::array<std::size_t,links> owner;
            nearest.fill(-1.0); owner.fill(faces);
            const bool bothBoundary=a.scalar[id]==static_cast<double>(CellType::Boundary) &&
                                    b.scalar[id]==static_cast<double>(CellType::Boundary);
            if(withQ && unite && bothBoundary)
            {
                const auto p=cellCenter(merged,id);
                for(std::size_t direction=1;direction<Lattice::nStencils;++direction)
                {
                    auto delta=p;
                    delta[0]=merged.gridSpacing*Lattice::ex[direction];
                    delta[1]=merged.gridSpacing*Lattice::ey[direction];
                    if constexpr(dimension==3) delta[2]=merged.gridSpacing*Lattice::ez[direction];
                    std::vector<Intersection> hits;
                    queryLinks(p,delta,0,hits);
                    std::sort(hits.begin(),hits.end(),[](const Intersection& x,const Intersection& y) {
                        return x.q<y.q || (x.q==y.q && x.face<y.face);
                    });
                    for(std::size_t begin=0;begin<hits.size();)
                    {
                        std::size_t end=begin+1;
                        std::vector<std::size_t> forced{hits[begin].face};
                        while(end<hits.size() && hits[end].q==hits[begin].q)
                            forced.push_back(hits[end++].face);
                        auto point=p;
                        for(std::size_t axis=0;axis<dimension;++axis)
                            point[axis]+=hits[begin].q*delta[axis];
                        if(queryPoint(point,complemented,forced,0)==Location::Surface)
                        {
                            std::vector<std::size_t> exposed;
                            querySurfaceFaces(point,complemented,forced,0,exposed);
                            if(!exposed.empty())
                            {
                                nearest[direction-1]=hits[begin].q;
                                owner[direction-1]=*std::min_element(exposed.begin(),exposed.end());
                                break;
                            }
                        }
                        begin=end;
                    }
                }
            }
            else if(withQ)
            {
                // Intersection leaves the result at the first child boundary.
                // A one-source union inherits that source's link unchanged.
                for(std::size_t i=linkOffsets[cell];i<linkOffsets[cell+1];++i)
                {
                    const auto& candidate=sourceLinks[i];
                    const std::size_t link=candidate.link;
                    if(nearest[link]<0.0 || candidate.q<nearest[link] ||
                       (candidate.q==nearest[link] && candidate.face<owner[link]))
                    { nearest[link]=candidate.q; owner[link]=candidate.face; }
                }
            }
            std::vector<std::size_t> cellFaces;
            for(std::size_t i=memberOffsets[cell];i<memberOffsets[cell+1];++i)
                cellFaces.push_back(memberships[i].second);
            if(withQ)
                for(std::size_t link=0;link<links;++link)
                    if(nearest[link]>0.0) cellFaces.push_back(owner[link]);
            std::sort(cellFaces.begin(),cellFaces.end());
            cellFaces.erase(std::unique(cellFaces.begin(),cellFaces.end()),cellFaces.end());
            for(std::size_t face:cellFaces)
            {
                Entry entry{id,face,{}};
                entry.q.fill(-1.0);
                if(withQ)
                    for(std::size_t link=0;link<links;++link)
                        if(owner[link]==face) entry.q[link]=nearest[link];
                perCell[cell].push_back(entry);
            }
        }
        std::vector<Entry> entries;
        for(const auto& group:perCell) entries.insert(entries.end(),group.begin(),group.end());
        std::sort(entries.begin(),entries.end(),[](const Entry& x,const Entry& y) {
            return std::tie(x.face,x.cell)<std::tie(y.face,y.cell);
        });
        if(entries.size()>std::numeric_limits<std::size_t>::max()/links)
            throw std::overflow_error("Boolean Q output size overflow.");
        merged.nBoundaries=faces;
        merged.boundaryOffsets.assign(faces+1,0);
        merged.boundaryCells.resize(entries.size());
        merged.boundaryX.resize(entries.size()); merged.boundaryY.resize(entries.size());
        if constexpr(dimension==3) merged.boundaryZ.resize(entries.size());
        if(withQ) merged.q.assign(entries.size()*links,-1.0);
        for(const auto& entry:entries) ++merged.boundaryOffsets[entry.face+1];
        for(std::size_t face=0;face<faces;++face) merged.boundaryOffsets[face+1]+=merged.boundaryOffsets[face];
#pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index=0;index<static_cast<std::ptrdiff_t>(entries.size());++index)
        {
            const std::size_t b=static_cast<std::size_t>(index);
            const auto p=cellCenter(merged,entries[b].cell);
            merged.boundaryCells[b]=entries[b].cell;
            merged.boundaryX[b]=p[0]; merged.boundaryY[b]=p[1];
            if constexpr(dimension==3) merged.boundaryZ[b]=p[2];
            if(withQ) std::copy(entries[b].q.begin(),entries[b].q.end(),merged.q.begin()+b*links);
        }
        result=std::move(merged);
    }

    static FlowType determineFlowType(
        const FlowType leftFlowType,
        const FlowType rightFlowType,
        const BooleanOperation booleanOperation)
    {
        if(booleanOperation ==
           BooleanOperation::Union)
        {
            if(leftFlowType !=
                   FlowType::Internal ||
               rightFlowType !=
                   FlowType::Internal)
            {
                throw std::invalid_argument(
                    "Geometry union requires two "
                    "internal-flow geometries.");
            }

            return FlowType::Internal;
        }


        if(booleanOperation ==
           BooleanOperation::Difference)
        {
            if(leftFlowType !=
                   FlowType::Internal ||
               rightFlowType !=
                   FlowType::Internal)
            {
                throw std::invalid_argument(
                    "Geometry difference requires two "
                    "internal-flow geometries.");
            }

            return FlowType::Internal;
        }


        if(booleanOperation ==
           BooleanOperation::Intersection)
        {
            if(leftFlowType == FlowType::Internal &&
               rightFlowType == FlowType::Internal)
            {
                return FlowType::Internal;
            }

            if(leftFlowType ==
                   FlowType::External &&
               rightFlowType ==
                   FlowType::External)
            {
                return FlowType::External;
            }

            if(leftFlowType ==
                   FlowType::Internal &&
               rightFlowType ==
                   FlowType::External)
            {
                return FlowType::Internal;
            }

            if(leftFlowType ==
                   FlowType::External &&
               rightFlowType ==
                   FlowType::Internal)
            {
                return FlowType::Internal;
            }

            throw std::invalid_argument(
                "Geometry intersection requires "
                "internal-internal, external-external or "
                "internal-external geometries.");
        }


        throw std::invalid_argument(
            "Unsupported geometry Boolean operation.");
    }


    static BoundingBoxType unionBounds(
        const BoundingBoxType& leftBounds,
        const BoundingBoxType& rightBounds)
    {
        BoundingBoxType result =
            leftBounds;

        for(std::size_t d = 0;
            d < result.min.size();
            ++d)
        {
            result.min[d] =
                std::min(
                    leftBounds.min[d],
                    rightBounds.min[d]);

            result.max[d] =
                std::max(
                    leftBounds.max[d],
                    rightBounds.max[d]);
        }

        return result;
    }


    static bool sameBounds(
        const BoundingBoxType& leftBounds,
        const BoundingBoxType& rightBounds)
    {
        for(std::size_t d = 0;
            d < leftBounds.min.size();
            ++d)
        {
            if(leftBounds.min[d] !=
                   rightBounds.min[d] ||
               leftBounds.max[d] !=
                   rightBounds.max[d])
            {
                return false;
            }
        }

        return true;
    }


    static BoundingBoxType determineBounds(
        const LeftGeometry& leftGeometry,
        const RightGeometry& rightGeometry,
        const BooleanOperation booleanOperation)
    {
        if(booleanOperation ==
           BooleanOperation::Union)
        {
            return unionBounds(
                leftGeometry.bounds,
                rightGeometry.bounds);
        }


        if(booleanOperation ==
           BooleanOperation::Difference)
        {
            return leftGeometry.bounds;
        }


        if(booleanOperation ==
           BooleanOperation::Intersection)
        {
            if(leftGeometry.flowType == FlowType::Internal &&
               rightGeometry.flowType == FlowType::Internal)
            {
                // The result is contained in the left internal domain. Keep a
                // valid analysis domain even when the intersection is empty.
                return leftGeometry.bounds;
            }

            if(leftGeometry.flowType ==
                   FlowType::External &&
               rightGeometry.flowType ==
                   FlowType::External)
            {
                return unionBounds(
                    leftGeometry.bounds,
                    rightGeometry.bounds);
            }

            if(leftGeometry.flowType ==
                   FlowType::Internal &&
               rightGeometry.flowType ==
                   FlowType::External)
            {
                return leftGeometry.bounds;
            }

            if(leftGeometry.flowType ==
                   FlowType::External &&
               rightGeometry.flowType ==
                   FlowType::Internal)
            {
                return rightGeometry.bounds;
            }
        }


        throw std::invalid_argument(
            "Cannot determine Boolean geometry bounds.");
    }


    static BoundingBoxType determineOpenBox(
        const LeftGeometry& leftGeometry,
        const RightGeometry& rightGeometry,
        const BooleanOperation booleanOperation)
    {
        if(booleanOperation ==
               BooleanOperation::Intersection &&
           leftGeometry.flowType ==
               FlowType::External &&
           rightGeometry.flowType ==
               FlowType::External)
        {
            if(!sameBounds(
                   leftGeometry.openBox,
                   rightGeometry.openBox))
            {
                throw std::invalid_argument(
                    "External geometry intersection "
                    "requires identical open boxes.");
            }

            return leftGeometry.openBox;
        }

        return BoundingBoxType{};
    }
};


namespace detail
{

template<typename LeftGeometry,typename RightGeometry>
struct GeometryQuery<BooleanGeometry<LeftGeometry,RightGeometry>>
{
    using Geometry=BooleanGeometry<LeftGeometry,RightGeometry>;

    static std::size_t faceCount(const Geometry& geometry)
    { return geometry.geometricFaceCount(); }

    template<typename Position>
    static PointLocation locate(const Geometry& geometry,const Position& p,bool complemented,
        const std::vector<std::size_t>& forcedFaces,std::size_t base)
    { return geometry.queryPoint(p,complemented,forcedFaces,base); }

    template<typename Position>
    static void intersections(const Geometry& geometry,const Position& p,const Position& delta,
        std::size_t base,std::vector<LinkIntersection>& output)
    { geometry.queryLinks(p,delta,base,output); }

    template<typename Position>
    static void surfaceFaces(const Geometry& geometry,const Position& p,bool complemented,
        const std::vector<std::size_t>& forcedFaces,std::size_t base,
        std::vector<std::size_t>& output)
    { geometry.querySurfaceFaces(p,complemented,forcedFaces,base,output); }

    template<typename Analysis,typename Domain>
    static void interior(const Geometry& geometry,Analysis& result,const Domain& domain,bool complemented)
    {
        result.domainBounds=domain;
        result.complemented=complemented;
        geometry.build(result,domain,complemented,false);
    }

    template<typename Analysis>
    static void boundary(const Geometry& geometry,Analysis& result,bool complemented)
    { geometry.completeBoundary(result,complemented); }

    template<typename Analysis,typename Domain>
    static void analyze(const Geometry& geometry,Analysis& result,const Domain& domain,bool complemented)
    {
        result.domainBounds=domain;
        result.complemented=complemented;
        geometry.build(result,domain,complemented,true);
    }
};

// Keep Boolean operators out of unrelated arithmetic and iterator overloads.
template<typename Geometry, typename = void>
struct IsGeometry : std::false_type {};

template<typename Geometry>
struct IsGeometry<Geometry, std::void_t<
    decltype(std::declval<Geometry>().bounds),
    decltype(std::declval<Geometry>().flowType)>> : std::true_type {};

} // namespace detail

template <
    typename LeftGeometry,
    typename RightGeometry,
    std::enable_if_t<
        detail::IsGeometry<LeftGeometry>::value &&
        detail::IsGeometry<RightGeometry>::value, int> = 0>
auto operator|(
    LeftGeometry left,
    RightGeometry right)
{
    return BooleanGeometry<
        LeftGeometry,
        RightGeometry>(
            std::move(left),
            std::move(right),
            BooleanOperation::Union);
}


template <
    typename LeftGeometry,
    typename RightGeometry,
    std::enable_if_t<
        detail::IsGeometry<LeftGeometry>::value &&
        detail::IsGeometry<RightGeometry>::value, int> = 0>
auto operator+(
    LeftGeometry left,
    RightGeometry right)
{
    return BooleanGeometry<
        LeftGeometry,
        RightGeometry>(
            std::move(left),
            std::move(right),
            BooleanOperation::Union);
}


template <
    typename LeftGeometry,
    typename RightGeometry,
    std::enable_if_t<
        detail::IsGeometry<LeftGeometry>::value &&
        detail::IsGeometry<RightGeometry>::value, int> = 0>
auto operator&(
    LeftGeometry left,
    RightGeometry right)
{
    return BooleanGeometry<
        LeftGeometry,
        RightGeometry>(
            std::move(left),
            std::move(right),
            BooleanOperation::Intersection);
}


template <
    typename LeftGeometry,
    typename RightGeometry,
    std::enable_if_t<
        detail::IsGeometry<LeftGeometry>::value &&
        detail::IsGeometry<RightGeometry>::value, int> = 0>
auto operator-(
    LeftGeometry left,
    RightGeometry right)
{
    return BooleanGeometry<
        LeftGeometry,
        RightGeometry>(
            std::move(left),
            std::move(right),
            BooleanOperation::Difference);
}

} // namespace ntic::lbm::geometry
