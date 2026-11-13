#include "test_multilevel_checks.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"

#include <cmath>
#include <limits>
#include <filesystem>
#include <functional>
#include <iostream>
#include <stdexcept>

namespace geo=ntic::lbm::geometry;
namespace lat=ntic::lbm::lattice;

namespace
{
using geo::test::require;

template<std::size_t D,lat::LatticeType Type,typename Geometry>
void compareGeometry(const geo::MultilevelGrid<D,Type>& grid,const Geometry& geometry)
{
    using Analysis=std::conditional_t<D==2,geo::GeometryAnalysis2D<Type>,geo::GeometryAnalysis3D<Type>>;
    geo::test::checkGeometry(grid,[&](Analysis& result) { geometry.analysis(result); });
}

void stlWall()
{
    const auto directory=std::filesystem::path(GEOMETRY_TEST_SOURCE_DIR)/"stl_tests"/"stls";
    auto data=ntic::lbm::stl::read((directory/"boundary_box.stl").string());
    ntic::lbm::stl::FacetTopology topologyData;
    ntic::lbm::stl::validate(data,"full",topologyData);
    {
        const auto raw=geo::volumeRegion(topologyData);
        require(raw.classify(geo::BoundingBox{{0.5,0.5,0.5},{0.75,0.75,0.75}})==geo::BoxRelation::Inside,"Raw STL volume selector failed");
    }
    geo::STLGeometry geometry(std::move(topologyData),geo::FlowType::Internal);
    geo::MultilevelOptions options; options.maxLevel=2; options.analysisTileWidth=2;
    geo::MultilevelGrid<3> grid(geometry.bounds,0.5,options);
    grid.generate(geometry,{geo::RefinementRule<3>::surface(geo::fluidRegion(geometry),1,0.15,{1})});
    geo::test::checkTopology(grid); compareGeometry(grid,geometry);
}

void stlNonbinaryGridCenters()
{
    const auto directory=std::filesystem::path(GEOMETRY_TEST_SOURCE_DIR)/"stl_tests"/"stls";
    auto data=ntic::lbm::stl::read((directory/"boundary_box.stl").string());
    ntic::lbm::stl::FacetTopology topologyData;
    ntic::lbm::stl::validate(data,"full",topologyData);
    geo::STLGeometry geometry(std::move(topologyData),geo::FlowType::Internal);
    const geo::Point origin{-75.0,-63.837493896484375,-61.046390533447266};
    geometry.translate(&origin);
    const double h=4.858944032554682/25.0;
    geo::GeometryAnalysis3D<lat::LatticeType::D3Q19> analysis(h);
    geometry.interiorAreaAnalysis(analysis);
    require(!analysis.boundaryX.empty(),"Nonbinary-spacing fixture has no boundary cells");
    analysis.validateBoundaryCells();
    auto interiorCells=analysis.boundaryCells;
    std::sort(interiorCells.begin(),interiorCells.end());
    const auto scalar=analysis.scalar;
    geometry.boundaryQAnalysis(analysis);
    require(analysis.scalar==scalar,"STL Q changed scalar classification");
    auto expandedCells=analysis.boundaryCells;
    std::sort(expandedCells.begin(),expandedCells.end());
    expandedCells.erase(std::unique(expandedCells.begin(),expandedCells.end()),expandedCells.end());
    require(expandedCells==interiorCells,"Face expansion changed physical boundary cell indices");
    const auto baseline=analysis;
    geometry.boundaryQAnalysis(analysis);
    require(analysis.boundaryCells==baseline.boundaryCells && analysis.q==baseline.q,
            "Repeated Q changed boundary indices or links");
    const auto rejects=[&](auto invalid,const char* message) {
        bool rejected=false;
        try { geometry.boundaryQAnalysis(invalid); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,message);
    };
    auto invalid=analysis;
    invalid.boundaryCells.front()=invalid.scalar.size();
    rejects(invalid,"Out-of-range boundary index was accepted");
    const auto other=std::find_if(analysis.boundaryCells.begin(),analysis.boundaryCells.end(),
        [&](auto id) { return id!=analysis.boundaryCells.front(); });
    require(other!=analysis.boundaryCells.end(),"Index fixture has no second Boundary cell");
    invalid=analysis; invalid.boundaryCells.front()=*other;
    rejects(invalid,"A different valid boundary index with mismatched coordinates was accepted");
    invalid=analysis; invalid.boundaryCells.pop_back();
    rejects(invalid,"Mismatched boundary index count was accepted");
    invalid=analysis;
    invalid.boundaryX.front()=std::nextafter(invalid.boundaryX.front(),std::numeric_limits<double>::infinity());
    rejects(invalid,"A one-ULP coordinate change was accepted");
    invalid=analysis; invalid.boundaryX.front()+=0.125*h;
    rejects(invalid,"A genuine one-eighth-cell coordinate offset was accepted");

    // Exercise the real multilevel tile origin/coordinate round trip with a
    // tiny 12-facet fixture, rather than rerunning the large channel model.
    geo::MultilevelOptions options;
    options.maxLevel=3; options.analysisTileWidth=2; options.analysisThreads=4;
    geo::MultilevelGrid<3> grid(geometry.bounds,8*h,options);
    grid.generate(geometry,{geo::RefinementRule<3>::surface(geo::fluidRegion(geometry),3,h)});
    require(grid.levels().size()==4 && grid.levels()[3].all.size()>0,"Nonbinary-spacing case did not reach level 3");
    geo::test::checkTopology(grid); compareGeometry(grid,geometry);

    // The same anchor must propagate through primitive and nested Boolean
    // analyses, including the temporary complement used by subtraction.
    const geo::Box box(geometry.bounds.min,geometry.bounds.max);
    const auto check=[&](const auto& shape) {
        geo::MultilevelGrid<3> other(shape.bounds,8*h,options);
        other.generate(shape,{geo::RefinementRule<3>::surface(geo::fluidRegion(shape),3,h)});
        geo::test::checkTopology(other); compareGeometry(other,shape);
    };
    check(box);
    const geo::Point center{origin[0]+1,origin[1]+1,origin[2]+1};
    const auto nested=(box-geo::Sphere(center,0.3)) & box;
    check(nested);
}

void stlNested()
{
    const auto directory=std::filesystem::path(GEOMETRY_TEST_SOURCE_DIR)/"stl_tests"/"stls";
    geo::MultilevelOptions options; options.maxLevel=2; options.analysisTileWidth=2;
    auto nestedData=ntic::lbm::stl::read((directory/"three_levels.stl").string());
    ntic::lbm::stl::FacetTopology nested;
    ntic::lbm::stl::validate(nestedData,"full",nested);
    const auto deep=geo::volumeRegion(nested); // No STLGeometry two-level flow restriction.
    const geo::Box enclosing({-12,-12,-12},{12,12,12});
    geo::MultilevelGrid<3> selectorGrid(enclosing.bounds,4,options);
    selectorGrid.generate(enclosing,{geo::RefinementRule<3>::volume(deep,1)});
    geo::test::checkTopology(selectorGrid);
}

} // namespace

int main()
{
    try
    {
        omp_set_dynamic(0); omp_set_num_threads(4);
        stlWall();
        std::cout << "[PASS] stl_wall\n";
        stlNonbinaryGridCenters();
        std::cout << "[PASS] stl_nonbinary_grid_centers\n";
        stlNested();
        std::cout << "[PASS] stl_nested\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "[FAIL] test_multilevel_stl: " << error.what() << '\n';
        return 1;
    }
}
