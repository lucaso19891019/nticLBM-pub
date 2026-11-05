#include "box.hpp"
#include "stl_geometry.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace ntic::lbm;
using namespace ntic::lbm::geometry;

namespace
{

void require(bool condition,const std::string& message)
{
    if(!condition) throw std::runtime_error(message);
}

std::filesystem::path stlFile(const char* filename)
{
    return std::filesystem::path(GEOMETRY_TEST_SOURCE_DIR)/"stl_tests"/"stls"/filename;
}

STLGeometry makeGeometry(const char* filename,FlowType flow,const BoundingBox& domain)
{
    auto data=stl::read(stlFile(filename));
    stl::FacetTopology topology;
    stl::validate(data,"full",topology);
    return STLGeometry(std::move(topology),flow,&domain);
}

template<lattice::LatticeType Type>
std::vector<double> canonical(const GeometryAnalysis3D<Type>& analysis)
{
    constexpr std::size_t links=GeometryAnalysis3D<Type>::Lattice::nStencils-1;
    std::vector<double> result(analysis.scalar.size()*links,-1.0);
    const std::size_t entries=analysis.boundaryX.size();
    require(analysis.nBoundaries>0 && analysis.boundaryOffsets.size()==analysis.nBoundaries+1,
            "Invalid Face count/offset size");
    require(analysis.boundaryOffsets.front()==0 && analysis.boundaryOffsets.back()==entries &&
            analysis.boundaryY.size()==entries && analysis.boundaryZ.size()==entries &&
            analysis.q.size()==entries*links,"Invalid coordinate/Q sizes");
    std::vector<bool> represented(analysis.scalar.size(),false);
    for(std::size_t face=0;face<analysis.nBoundaries;++face)
    {
        const std::size_t begin=analysis.boundaryOffsets[face],end=analysis.boundaryOffsets[face+1];
        require(begin<=end,"Offsets must be monotone");
        std::size_t previous=0;
        for(std::size_t b=begin;b<end;++b)
        {
            const auto i=static_cast<std::size_t>(std::llround(
                (analysis.boundaryX[b]-analysis.domain.min[0])/analysis.gridSpacing-0.5));
            const auto j=static_cast<std::size_t>(std::llround(
                (analysis.boundaryY[b]-analysis.domain.min[1])/analysis.gridSpacing-0.5));
            const auto k=static_cast<std::size_t>(std::llround(
                (analysis.boundaryZ[b]-analysis.domain.min[2])/analysis.gridSpacing-0.5));
            require(i<analysis.nx && j<analysis.ny && k<analysis.nz,"Entry outside grid");
            const std::size_t id=i+analysis.nx*(j+analysis.ny*k);
            require(b==begin || previous<id,"Non-unique or unstable (Face, Cell)");
            previous=id;
            require(analysis.scalar[id]==static_cast<double>(CellType::Boundary),
                    "Entry is not an interior-classified Boundary cell");
            represented[id]=true;
            for(std::size_t link=0;link<links;++link)
            {
                const double q=analysis.q[b*links+link];
                require(q==-1.0 || (std::isfinite(q) && q>0.0 && q<=1.0),"Invalid Q range");
                if(q==-1.0) continue;
                double& value=result[id*links+link];
                require(value==-1.0,"Duplicate physical (Cell, Link) Q");
                value=q;
            }
        }
    }
    for(std::size_t id=0;id<represented.size();++id)
        require(represented[id]==(analysis.scalar[id]==static_cast<double>(CellType::Boundary)),
                "A physical Boundary cell was lost");
    return result;
}

template<lattice::LatticeType Type>
void repeat(const STLGeometry& geometry,const GeometryAnalysis3D<Type>& baseline)
{
    auto analysis=baseline;
    // Standalone Q is idempotent even after coordinates have been Face-expanded.
    geometry.boundaryQAnalysis(analysis);
    for(int threads:{1,4,4})
    {
        require(analysis.scalar==baseline.scalar && analysis.nBoundaries==baseline.nBoundaries &&
                analysis.boundaryOffsets==baseline.boundaryOffsets &&
                analysis.boundaryX==baseline.boundaryX && analysis.boundaryY==baseline.boundaryY &&
                analysis.boundaryZ==baseline.boundaryZ && analysis.q==baseline.q,
                "Repeated/threaded STL output differs");
        omp_set_num_threads(threads);
        geometry.analysis(analysis);
    }
    require(analysis.scalar==baseline.scalar && analysis.boundaryOffsets==baseline.boundaryOffsets &&
            analysis.boundaryX==baseline.boundaryX && analysis.boundaryY==baseline.boundaryY &&
            analysis.boundaryZ==baseline.boundaryZ && analysis.q==baseline.q,
            "Final threaded STL output differs");
}

template<lattice::LatticeType Type>
void compareBox(bool external,bool nested,bool exclude)
{
    const BoundingBox domain{{-1.0,-1.0,-1.0},{3.0,3.0,3.0}};
    const FlowType flow=external ? FlowType::External : FlowType::Internal;
    STLGeometry geometry=makeGeometry(nested ? "boundary_nested_box.stl" : "boundary_box.stl",flow,domain);
    if(exclude) geometry.identifyBoundaryFeatures({1});
    if(nested)
        require(std::count_if(geometry.flow.begin(),geometry.flow.end(),
                    [](const auto& component){return !component.active;})==1,
                "External nested fixture must contain an inactive component");
    for(double h:{0.25,0.5,1.0,4.0})
    {
        omp_set_num_threads(1);
        GeometryAnalysis3D<Type> analysis(h),primitive(h);
        geometry.interiorAreaAnalysis(analysis);
        const auto scalar=analysis.scalar;
        geometry.boundaryQAnalysis(analysis);
        require(analysis.scalar==scalar,"STL Q changed CellType");
        const Box box{{0.0,0.0,0.0},{2.0,2.0,2.0},flow,&domain};
        box.analysis(primitive);
        require(scalar==primitive.scalar,"STL/Primitive Box classification differs");
        const auto actual=canonical(analysis),expected=canonical(primitive);
        require(actual.size()==expected.size(),"Canonical grid mismatch");
        for(std::size_t link=0;link<actual.size();++link)
            require(std::abs(actual[link]-expected[link])<=1.0e-12,
                    "STL Q differs from analytic Box Q");
        require(analysis.nBoundaries==(external ? 1 : (exclude ? 6 : 7)),
                "Incorrect STL Other/planar Face count");
        if(!external && !exclude)
            require(analysis.boundaryOffsets[0]==analysis.boundaryOffsets[1],
                    "Box Other Surface should be empty");
        if(exclude && std::count(scalar.begin(),scalar.end(),
                                static_cast<double>(CellType::Boundary))>0)
            require(analysis.boundaryOffsets[1]>0,"Excluded plane must belong to Face 0");
        repeat(geometry,analysis);
    }
}

// Manufactured interior output isolates Q edge cases without relying on flood fill.
template<lattice::LatticeType Type>
GeometryAnalysis3D<Type> probe(STLGeometry geometry,const Point& p,
    const char* filename,const std::vector<std::size_t>& faces)
{
    auto data=stl::read(stlFile(filename));
    stl::FacetTopology topology;
    // Open triangle fixtures exercise Q only, so use facet-level validation.
    stl::validate(data,"test",topology);
    require(topology.geometry.facetVertexIDs.size()==faces.size(),"Fixture Face map mismatch");
    topology.components.resize(1);
    topology.facetComponentIDs.assign(faces.size(),0);
    geometry.topology=std::move(topology);
    geometry.flow.resize(1);
    geometry.flow[0].active=true;
    geometry.flowType=FlowType::Internal;
    geometry.boundaryFeatures.nPlanarFaces=*std::max_element(faces.begin(),faces.end());
    geometry.boundaryFeatures.facetFaceIDs=faces;
    GeometryAnalysis3D<Type> analysis(0.5);
    analysis.nx=analysis.ny=analysis.nz=1;
    analysis.domain={{p[0]-0.25,p[1]-0.25,p[2]-0.25},{p[0]+0.25,p[1]+0.25,p[2]+0.25}};
    analysis.scalar={static_cast<double>(CellType::Boundary)};
    analysis.boundaryX={p[0]}; analysis.boundaryY={p[1]}; analysis.boundaryZ={p[2]};
    geometry.boundaryQAnalysis(analysis);
    canonical(analysis);
    const auto before=analysis;
    // Reverse facet iteration while preserving Face IDs and geometry.
    std::reverse(geometry.topology.geometry.facetVertexIDs.begin(),
                 geometry.topology.geometry.facetVertexIDs.end());
    std::reverse(geometry.boundaryFeatures.facetFaceIDs.begin(),geometry.boundaryFeatures.facetFaceIDs.end());
    omp_set_num_threads(4);
    geometry.boundaryQAnalysis(analysis);
    require(analysis.boundaryOffsets==before.boundaryOffsets && analysis.q==before.q &&
            analysis.boundaryX==before.boundaryX && analysis.boundaryY==before.boundaryY &&
            analysis.boundaryZ==before.boundaryZ,"Facet traversal changed Q ownership");
    return analysis;
}

template<lattice::LatticeType Type>
void assertHit(const GeometryAnalysis3D<Type>& analysis,std::size_t face,double expected)
{
    using Lattice=typename GeometryAnalysis3D<Type>::Lattice;
    std::size_t direction=0;
    for(std::size_t link=1;link<Lattice::nStencils;++link)
        if(Lattice::ex[link]==1 && Lattice::ey[link]==0 && Lattice::ez[link]==0) direction=link;
    require(direction!=0,"Missing axial lattice direction");
    std::size_t hits=0;
    for(std::size_t f=0;f<analysis.nBoundaries;++f)
        for(std::size_t b=analysis.boundaryOffsets[f];b<analysis.boundaryOffsets[f+1];++b)
        {
            const double q=analysis.q[b*(Lattice::nStencils-1)+direction-1];
            if(q<0.0) continue;
            ++hits;
            require(f==face && std::abs(q-expected)<=1.0e-12,"Incorrect nearest/shared hit");
        }
    require(hits==(expected<0.0 ? 0u : 1u),"Incorrect axial hit count");
}

template<lattice::LatticeType Type>
void intersectionCases()
{
    const BoundingBox domain{{-1.0,-1.0,-1.0},{3.0,3.0,3.0}};
    const auto geometry=makeGeometry("boundary_box.stl",FlowType::Internal,domain);
    assertHit(probe<Type>(geometry,{0.0,0.0,0.0},"boundary_nearest_planes.stl",{1,2}),2,0.5);
    assertHit(probe<Type>(geometry,{0.0,0.0,0.0},"boundary_shared_edge.stl",{2,1}),1,0.5);
    assertHit(probe<Type>(geometry,{0.0,0.0,0.0},"boundary_endpoint.stl",{0}),0,1.0);
    assertHit(probe<Type>(geometry,{0.0,0.0,0.0},"boundary_behind.stl",{0}),0,-1.0);
    assertHit(probe<Type>(geometry,{-0.25,0.0,0.0},"boundary_coplanar.stl",{0}),0,0.5);
    const auto stencilOnly=probe<Type>(geometry,{-0.25,0.0,0.125},"boundary_coplanar.stl",{0});
    assertHit(stencilOnly,0,-1.0);
    require(!stencilOnly.boundaryX.empty() &&
            std::all_of(stencilOnly.q.begin(),stencilOnly.q.end(),[](double q){return q==-1.0;}),
            "Stencil-only Boundary entry was lost or assigned a false Q");
    // Shared diagonal of two triangles, with hits on either side of the seam.
    const double epsilon=std::ldexp(1.0,-36);
    for(double y:{-epsilon,0.0,epsilon})
        assertHit(probe<Type>(geometry,{0.0,y,0.0},"boundary_shared_edge.stl",{1,2}),y<0.0 ? 2 : 1,0.5);
    assertHit(probe<Type>(geometry,{0.0,1.0,1.0},"boundary_shared_edge.stl",{1,2}),1,0.5);
    bool rejected=false;
    auto degenerate=stl::read(stlFile("boundary_degenerate.stl"));
    require(degenerate.size()==1,"Degenerate fixture must contain one facet");
    stl::FacetTopology rejectedTopology;
    try { stl::validate(degenerate,"test",rejectedTopology); }
    catch(const std::runtime_error&) { rejected=true; }
    require(rejected,"Degenerate STL fixture was not rejected by validation");
}

void branchedChannel()
{
    auto data=stl::read(stlFile("smooth_irregular_branched_channel.stl"));
    stl::FacetTopology topology;
    stl::validate(data,"full",topology);
    const STLGeometry geometry(std::move(topology),FlowType::Internal);
    require(geometry.boundaryFeatures.nPlanarFaces==7,"Branched channel should have seven planes");
    const double h=std::max({geometry.bounds.width(),geometry.bounds.height(),geometry.bounds.depth()})/12.0;
    GeometryAnalysis3D<lattice::LatticeType::D3Q19> analysis(h);
    omp_set_num_threads(1);
    geometry.interiorAreaAnalysis(analysis);
    const auto scalar=analysis.scalar;
    geometry.boundaryQAnalysis(analysis);
    require(analysis.scalar==scalar && analysis.nBoundaries==8,"Branched channel contract failed");
    canonical(analysis);
    require(std::any_of(analysis.q.begin(),analysis.q.end(),[](double q){return q>0.0;}),
            "Branched channel has no valid Q");
    require(analysis.boundaryOffsets[1]>0,"Branched channel lost Other Surface");
    repeat(geometry,analysis);
}

template<lattice::LatticeType Type>
void checkModel()
{
    compareBox<Type>(false,false,false);
    compareBox<Type>(false,false,true);
    compareBox<Type>(true,false,false);
    compareBox<Type>(true,true,false);
    intersectionCases<Type>();
}

} // namespace

int main()
{
    omp_set_dynamic(0);
    try
    {
        checkModel<lattice::LatticeType::D3Q15>();
        checkModel<lattice::LatticeType::D3Q19>();
        checkModel<lattice::LatticeType::D3Q27>();
        branchedChannel();
        std::cout<<"STL boundary acceptance passed\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"STL boundary acceptance failed: "<<error.what()<<'\n';
        return 1;
    }
}
