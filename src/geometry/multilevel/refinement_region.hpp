#pragma once

#include "boolean_geometry.hpp"

#include <functional>
#include <memory>
#include <numeric>

namespace ntic::lbm::geometry
{

template<std::size_t D>
using Position = std::array<double,D>;

template<std::size_t D>
using RegionBox = std::conditional_t<D == 2, BoundingBox2D, BoundingBox>;

enum class BoxRelation { Outside, Mixed, Inside };

namespace detail
{

template<std::size_t D>
Position<D> center(const RegionBox<D>& box)
{
    Position<D> result{};
    for(std::size_t d=0; d<D; ++d) result[d]=box.min[d]+0.5*(box.max[d]-box.min[d]);
    return result;
}

template<std::size_t D>
double radius(const RegionBox<D>& box)
{
    double squared=0;
    for(std::size_t d=0; d<D; ++d) squared+=0.25*std::pow(box.max[d]-box.min[d],2);
    return std::sqrt(squared);
}

template<std::size_t D>
RegionBox<D> expand(RegionBox<D> box, double thickness)
{
    for(std::size_t d=0; d<D; ++d) { box.min[d]-=thickness; box.max[d]+=thickness; }
    return box;
}

template<std::size_t D>
double boxDistance(const Position<D>& p, const RegionBox<D>& box)
{
    double squared=0;
    for(std::size_t d=0; d<D; ++d)
    {
        const double delta=std::max({box.min[d]-p[d],0.0,p[d]-box.max[d]});
        squared+=delta*delta;
    }
    return std::sqrt(squared);
}

template<std::size_t D>
double signedBoxDistance(const Position<D>& p, const RegionBox<D>& box)
{
    const double outside=boxDistance<D>(p,box);
    if(outside>0) return outside;
    double inside=std::numeric_limits<double>::infinity();
    for(std::size_t d=0; d<D; ++d) inside=std::min({inside,p[d]-box.min[d],box.max[d]-p[d]});
    return -inside;
}

inline BoxRelation invert(BoxRelation state)
{
    return state==BoxRelation::Inside ? BoxRelation::Outside :
           state==BoxRelation::Outside ? BoxRelation::Inside : BoxRelation::Mixed;
}

template<std::size_t D>
struct Node
{
    std::size_t faces=0;
    virtual ~Node() = default;
    virtual BoxRelation classify(const RegionBox<D>& box) const = 0;
    virtual bool nearSurface(const RegionBox<D>& box, double thickness,
        const std::vector<std::size_t>& selected, std::size_t base) const = 0;
};

template<std::size_t D>
struct Leaf final : Node<D>
{
    std::function<double(const Position<D>&)> signedDistance;
    std::function<double(const Position<D>&,const std::vector<std::size_t>&,std::size_t)> distance;
    std::function<BoxRelation(const RegionBox<D>&)> boxClassification;

    BoxRelation classify(const RegionBox<D>& box) const override
    {
        if(boxClassification) return boxClassification(box);
        const double value=signedDistance(center<D>(box)), r=radius<D>(box);
        return value>r ? BoxRelation::Outside : value<-r ? BoxRelation::Inside : BoxRelation::Mixed;
    }

    bool nearSurface(const RegionBox<D>& box, double thickness,
        const std::vector<std::size_t>& selected, std::size_t base) const override
    {
        // A bounding sphere gives a conservative Euclidean distance bound.
        // False positives near corners are permitted; thin surfaces cannot be
        // lost through cell-center sampling of the requested band.
        return distance(center<D>(box),selected,base)<=thickness+radius<D>(box);
    }
};

template<std::size_t D>
struct Composite final : Node<D>
{
    enum class Operation { Union, Intersection, Complement };
    using Ptr=std::shared_ptr<const Node<D>>;
    Ptr left,right;
    Operation operation;

    Composite(Ptr a, Ptr b, Operation op) : left(std::move(a)), right(std::move(b)), operation(op)
    { this->faces=left->faces+(right ? right->faces : 0); }

    BoxRelation classify(const RegionBox<D>& box) const override
    {
        const auto a=left->classify(box);
        if(operation==Operation::Complement) return invert(a);
        const auto b=right->classify(box);
        if(operation==Operation::Union)
        {
            if(a==BoxRelation::Inside || b==BoxRelation::Inside) return BoxRelation::Inside;
            if(a==BoxRelation::Outside && b==BoxRelation::Outside) return BoxRelation::Outside;
        }
        else
        {
            if(a==BoxRelation::Outside || b==BoxRelation::Outside) return BoxRelation::Outside;
            if(a==BoxRelation::Inside && b==BoxRelation::Inside) return BoxRelation::Inside;
        }
        return BoxRelation::Mixed;
    }

    bool nearSurface(const RegionBox<D>& box, double thickness,
        const std::vector<std::size_t>& selected, std::size_t base) const override
    {
        if(operation==Operation::Complement) return left->nearSurface(box,thickness,selected,base);
        const auto expanded=expand<D>(box,thickness);
        const auto a=left->classify(expanded), b=right->classify(expanded);
        const bool unite=operation==Operation::Union;
        const bool leftExposed=unite ? b!=BoxRelation::Inside : b!=BoxRelation::Outside;
        const bool rightExposed=unite ? a!=BoxRelation::Inside : a!=BoxRelation::Outside;
        return (leftExposed && left->nearSurface(box,thickness,selected,base)) ||
               (rightExposed && right->nearSurface(box,thickness,selected,base+left->faces));
    }
};

inline bool selectedFace(std::size_t face, const std::vector<std::size_t>& selected)
{ return selected.empty() || std::find(selected.begin(),selected.end(),face)!=selected.end(); }

template<std::size_t D>
std::shared_ptr<Leaf<D>> boxLeaf(RegionBox<D> bounds, std::array<std::size_t,2*D> ids)
{
    auto result=std::make_shared<Leaf<D>>(); result->faces=2*D;
    result->signedDistance=[bounds](const Position<D>& p) { return signedBoxDistance<D>(p,bounds); };
    result->boxClassification=[bounds](const RegionBox<D>& cell) {
        bool inside=true;
        for(std::size_t d=0; d<D; ++d)
        {
            // Boundary-only contact has no volume to refine.
            if(cell.max[d]<=bounds.min[d] || cell.min[d]>=bounds.max[d]) return BoxRelation::Outside;
            inside=inside && cell.min[d]>=bounds.min[d] && cell.max[d]<=bounds.max[d];
        }
        return inside ? BoxRelation::Inside : BoxRelation::Mixed;
    };
    result->distance=[bounds,ids](const Position<D>& p,const auto& selected,std::size_t base) {
        double distance=std::numeric_limits<double>::infinity();
        for(std::size_t axis=0; axis<D; ++axis)
            for(std::size_t side=0; side<2; ++side)
            {
                if(!selectedFace(base+ids[2*axis+side],selected)) continue;
                auto face=bounds;
                face.min[axis]=face.max[axis]=side==0 ? bounds.min[axis] : bounds.max[axis];
                distance=std::min(distance,boxDistance<D>(p,face));
            }
        return distance;
    };
    return result;
}

template<std::size_t D>
std::shared_ptr<Leaf<D>> roundLeaf(Position<D> center, double radius)
{
    auto result=std::make_shared<Leaf<D>>(); result->faces=1;
    result->signedDistance=[center,radius](const Position<D>& p) {
        double squared=0; for(std::size_t d=0; d<D; ++d) squared+=std::pow(p[d]-center[d],2);
        return std::sqrt(squared)-radius;
    };
    const auto sdf=result->signedDistance;
    result->distance=[sdf](const Position<D>& p,const auto& selected,std::size_t base) {
        return selectedFace(base,selected) ? std::abs(sdf(p)) : std::numeric_limits<double>::infinity();
    };
    return result;
}

// A read-only BVH over the existing validated topology. No vertex/facet copy.
// It is built once per STL selector and shared by all OpenMP selection queries.
class STLDistance
{
    struct Branch
    {
        BoundingBox box;
        std::size_t first=0,count=0,left=0,right=0;
    };
    const stl::FacetTopology* topology_;
    std::vector<std::size_t> facets_,faces_;
    std::vector<Branch> tree_;

    static Point subtract(const Point& a,const Point& b)
    { return {a[0]-b[0],a[1]-b[1],a[2]-b[2]}; }
    static double dot(const Point& a,const Point& b)
    { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
    static double segmentDistance(const Point& p,const Point& a,const Point& b)
    {
        const auto ab=subtract(b,a),ap=subtract(p,a);
        const double length=dot(ab,ab);
        const double t=length>0 ? std::clamp(dot(ap,ab)/length,0.0,1.0) : 0.0;
        Point delta{}; for(std::size_t d=0; d<3; ++d) delta[d]=ap[d]-t*ab[d];
        return std::sqrt(dot(delta,delta));
    }
    double triangleDistance(const Point& p,std::size_t facet) const
    {
        const auto& ids=topology_->geometry.facetVertexIDs[facet];
        const auto& a=topology_->geometry.vertices[ids[0]];
        const auto& b=topology_->geometry.vertices[ids[1]];
        const auto& c=topology_->geometry.vertices[ids[2]];
        const auto ab=subtract(b,a),ac=subtract(c,a),ap=subtract(p,a);
        const double aa=dot(ab,ab),bb=dot(ab,ac),cc=dot(ac,ac);
        const double pa=dot(ap,ab),pc=dot(ap,ac),det=aa*cc-bb*bb;
        if(det>0)
        {
            const double u=(pa*cc-pc*bb)/det,v=(pc*aa-pa*bb)/det;
            if(u>=0 && v>=0 && u+v<=1)
            {
                Point delta{}; for(std::size_t d=0; d<3; ++d) delta[d]=ap[d]-u*ab[d]-v*ac[d];
                return std::sqrt(dot(delta,delta));
            }
        }
        return std::min({segmentDistance(p,a,b),segmentDistance(p,b,c),segmentDistance(p,c,a)});
    }
    std::size_t build(std::size_t first,std::size_t count)
    {
        const auto index=tree_.size(); tree_.emplace_back();
        BoundingBox box;
        for(std::size_t i=first; i<first+count; ++i)
            for(const auto v:topology_->geometry.facetVertexIDs[facets_[i]]) box.expand(topology_->geometry.vertices[v]);
        tree_[index].box=box;
        if(count<=8) { tree_[index].first=first; tree_[index].count=count; return index; }
        std::size_t axis=0;
        for(std::size_t d=1; d<3; ++d)
            if(box.max[d]-box.min[d]>box.max[axis]-box.min[axis]) axis=d;
        const auto coordinate=[&](std::size_t facet) {
            double value=0; for(auto v:topology_->geometry.facetVertexIDs[facet]) value+=topology_->geometry.vertices[v][axis]/3.0;
            return value;
        };
        const auto middle=first+count/2;
        std::nth_element(facets_.begin()+first,facets_.begin()+middle,facets_.begin()+first+count,
            [&](auto a,auto b) { const auto x=coordinate(a),y=coordinate(b); return x<y || (x==y && a<b); });
        const auto left=build(first,count/2),right=build(middle,count-count/2);
        tree_[index].left=left; tree_[index].right=right;
        return index;
    }
    void nearest(std::size_t index,const Point& p,const std::vector<std::size_t>& selected,
                 std::size_t base,double& best) const
    {
        const auto& node=tree_[index];
        if(boxDistance<3>(p,node.box)>best) return;
        if(node.count)
        {
            for(std::size_t i=node.first; i<node.first+node.count; ++i)
            {
                const auto facet=facets_[i];
                if(selectedFace(base+faces_[facet],selected)) best=std::min(best,triangleDistance(p,facet));
            }
            return;
        }
        const bool leftFirst=boxDistance<3>(p,tree_[node.left].box)<=boxDistance<3>(p,tree_[node.right].box);
        nearest(leftFirst ? node.left : node.right,p,selected,base,best);
        nearest(leftFirst ? node.right : node.left,p,selected,base,best);
    }
public:
    // topology must outlive this selector. Active/Face arrays are compact metadata.
    STLDistance(const stl::FacetTopology& topology,std::vector<std::size_t> facets,
                std::vector<std::size_t> faces)
        : topology_(&topology),facets_(std::move(facets)),faces_(std::move(faces))
    {
        if(facets_.empty() || faces_.size()!=topology.geometry.facetVertexIDs.size())
            throw std::invalid_argument("STL refinement selector requires validated nonempty topology.");
        tree_.reserve(2*facets_.size()); build(0,facets_.size());
    }
    double distance(const Point& p,const std::vector<std::size_t>& selected={},std::size_t base=0) const
    {
        double result=std::numeric_limits<double>::infinity(); nearest(0,p,selected,base,result); return result;
    }
};

} // namespace detail

template<std::size_t D>
class RefinementRegion
{
    using Node=detail::Node<D>;
    std::shared_ptr<const Node> node_;
public:
    explicit RefinementRegion(std::shared_ptr<const Node> node) : node_(std::move(node))
    { if(!node_) throw std::invalid_argument("Null refinement region."); }
    std::size_t faceCount() const { return node_->faces; }
    BoxRelation classify(const RegionBox<D>& box) const { return node_->classify(box); }
    bool nearSurface(const RegionBox<D>& box,double thickness,const std::vector<std::size_t>& faces={}) const
    { return node_->nearSurface(box,thickness,faces,0); }
    friend RefinementRegion operator!(const RefinementRegion& a)
    {
        using Composite=detail::Composite<D>;
        return RefinementRegion(std::make_shared<Composite>(a.node_,nullptr,Composite::Operation::Complement));
    }
    friend RefinementRegion operator|(const RefinementRegion& a,const RefinementRegion& b)
    {
        using Composite=detail::Composite<D>;
        return RefinementRegion(std::make_shared<Composite>(a.node_,b.node_,Composite::Operation::Union));
    }
    friend RefinementRegion operator&(const RefinementRegion& a,const RefinementRegion& b)
    {
        using Composite=detail::Composite<D>;
        return RefinementRegion(std::make_shared<Composite>(a.node_,b.node_,Composite::Operation::Intersection));
    }
    friend RefinementRegion operator-(const RefinementRegion& a,const RefinementRegion& b) { return a & !b; }
};

// Volume selectors deliberately ignore primitive FlowType: select the shape's interior.
inline RefinementRegion<2> volumeRegion(const Rectangle& g)
{ return RefinementRegion<2>(detail::boxLeaf<2>(g.bounds,{3,1,0,2})); }
inline RefinementRegion<3> volumeRegion(const Box& g)
{ return RefinementRegion<3>(detail::boxLeaf<3>(g.bounds,{0,1,2,3,4,5})); }
inline RefinementRegion<2> volumeRegion(const Circle& g)
{ return RefinementRegion<2>(detail::roundLeaf<2>(g.center,g.radius)); }
inline RefinementRegion<3> volumeRegion(const Sphere& g)
{ return RefinementRegion<3>(detail::roundLeaf<3>(g.center,g.radius)); }
inline RefinementRegion<3> volumeRegion(const Cylinder& g)
{
    auto node=std::make_shared<detail::Leaf<3>>(); node->faces=3;
    const auto axes=detail::cylinderDimensions(g);
    const auto center=g.center; const auto radius=g.radius,half=0.5*g.length;
    node->signedDistance=[axes,center,radius,half](const Point& p) {
        const double r=std::hypot(p[axes[1]]-center[axes[1]],p[axes[2]]-center[axes[2]])-radius;
        const double z=std::abs(p[axes[0]]-center[axes[0]])-half;
        return std::hypot(std::max(r,0.0),std::max(z,0.0))+std::min(std::max(r,z),0.0);
    };
    node->distance=[axes,center,radius,half](const Point& p,const auto& selected,std::size_t base) {
        const double r=std::hypot(p[axes[1]]-center[axes[1]],p[axes[2]]-center[axes[2]]);
        const double z=p[axes[0]]-center[axes[0]];
        double best=std::numeric_limits<double>::infinity();
        if(detail::selectedFace(base,selected)) best=std::hypot(r-radius,std::max(std::abs(z)-half,0.0));
        for(int sign:{-1,1})
            if(detail::selectedFace(base+(sign<0 ? 1 : 2),selected))
                best=std::min(best,std::hypot(z-sign*half,std::max(r-radius,0.0)));
        return best;
    };
    return RefinementRegion<3>(std::move(node));
}

// Raw validated STL selectors have independent even/odd solid semantics. They
// support arbitrary component nesting; no STLGeometry flow interpretation occurs.
inline RefinementRegion<3> volumeRegion(const stl::FacetTopology& topology)
{
    std::vector<std::size_t> facets(topology.geometry.facetVertexIDs.size());
    std::iota(facets.begin(),facets.end(),0);
    auto distance=std::make_shared<detail::STLDistance>(topology,std::move(facets),
        std::vector<std::size_t>(topology.geometry.facetVertexIDs.size(),0));
    auto node=std::make_shared<detail::Leaf<3>>(); node->faces=1;
    node->signedDistance=[distance,&topology](const Point& p) {
        const double d=distance->distance(p); if(d==0) return 0.0;
        bool inside=false;
        for(std::size_t c=0; c<topology.components.size(); ++c)
        {
            const auto& box=topology.components[c].bounds;
            bool contained=true;
            for(std::size_t axis=0; axis<3; ++axis) contained=contained && p[axis]>=box.min[axis] && p[axis]<=box.max[axis];
            if(contained && detail::pointInClosedComponent(topology,c,p)) inside=!inside;
        }
        return inside ? -d : d;
    };
    node->distance=[distance](const Point& p,const auto& selected,std::size_t base) { return distance->distance(p,selected,base); };
    return RefinementRegion<3>(std::move(node));
}

inline RefinementRegion<3> fluidRegion(const STLGeometry& g)
{
    std::vector<std::size_t> facets;
    const auto count=g.topology.geometry.facetVertexIDs.size();
    std::vector<std::size_t> faces(count,0);
    for(std::size_t f=0; f<count; ++f)
        if(g.flow.at(g.topology.facetComponentIDs[f]).active)
        {
            facets.push_back(f);
            if(g.flowType==FlowType::Internal) faces[f]=g.boundaryFeatures.facetFaceIDs.at(f);
        }
    auto distance=std::make_shared<detail::STLDistance>(g.topology,std::move(facets),std::move(faces));
    auto node=std::make_shared<detail::Leaf<3>>(); node->faces=detail::GeometryQuery<STLGeometry>::faceCount(g);
    node->signedDistance=[distance,&g](const Point& p) {
        const double d=distance->distance(p); if(d==0) return 0.0;
        // The nearest-surface query already excluded on-surface points. Reuse
        // STL fluid semantics without scanning every triangle a second time.
        return detail::fluidContainsPoint(g,p) ? -d : d;
    };
    node->distance=[distance](const Point& p,const auto& selected,std::size_t base) { return distance->distance(p,selected,base); };
    return RefinementRegion<3>(std::move(node));
}

template<typename Geometry>
auto fluidRegion(const Geometry& g)
{
    auto region=volumeRegion(g);
    return g.flowType==FlowType::Internal ? region : !region;
}

template<typename A,typename B>
auto fluidRegion(const BooleanGeometry<A,B>& g)
{
    const auto a=fluidRegion(g.left),b=fluidRegion(g.right);
    if(g.operation==BooleanOperation::Union) return a | b;
    if(g.operation==BooleanOperation::Difference) return a - b;
    return a & b;
}

template<typename A,typename B>
auto volumeRegion(const BooleanGeometry<A,B>& g)
{
    auto region=fluidRegion(g);
    return g.flowType==FlowType::Internal ? region : !region;
}

template<std::size_t D>
struct RefinementRule
{
    RefinementRegion<D> region;
    unsigned targetLevel;
    bool surfaceBand=false;
    double thickness=0;
    std::vector<std::size_t> faces;

    static RefinementRule volume(RefinementRegion<D> region,unsigned level)
    { return {std::move(region),level,false,0,{}}; }
    static RefinementRule surface(RefinementRegion<D> region,unsigned level,double thickness,
                                  std::vector<std::size_t> faces={})
    {
        if(!(thickness>=0) || !std::isfinite(thickness))
            throw std::invalid_argument("Surface refinement thickness must be finite and nonnegative.");
        for(auto face:faces) if(face>=region.faceCount()) throw std::invalid_argument("Invalid refinement Face ID.");
        return {std::move(region),level,true,thickness,std::move(faces)};
    }
    bool intersects(const RegionBox<D>& box) const
    {
        if(!surfaceBand) return region.classify(box)!=BoxRelation::Outside;
        // Restrict the band to the selected fluid region; never flip wall normals.
        return region.classify(box)!=BoxRelation::Outside && region.nearSurface(box,thickness,faces);
    }
};

} // namespace ntic::lbm::geometry
