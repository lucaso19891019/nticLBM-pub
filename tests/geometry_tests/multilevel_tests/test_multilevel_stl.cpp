#include "test_multilevel_checks.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"

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
