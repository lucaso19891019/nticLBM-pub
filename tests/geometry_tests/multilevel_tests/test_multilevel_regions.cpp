#include "multilevel/multilevel_grid.hpp"

#include <functional>
#include <iostream>
#include <stdexcept>

namespace geo=ntic::lbm::geometry;
namespace lat=ntic::lbm::lattice;

// Shared verification is implemented in multilevel_test_checks.cpp.
namespace ntic::lbm::geometry::test
{
void require(bool condition,const char* message);
geo::MultilevelGrid<2> makeRegionGrid(int threads);
void checkTopology(const geo::MultilevelGrid<2>& grid);
void checkGeometry(
    const geo::MultilevelGrid<2>& grid,
    const std::function<void(geo::GeometryAnalysis2D<lat::LatticeType::D2Q9>&)>& analyze);
} // namespace ntic::lbm::geometry::test

namespace
{
using geo::test::require;

template<std::size_t D,lat::LatticeType Type,typename Geometry>
void compareGeometry(const geo::MultilevelGrid<D,Type>& grid,const Geometry& geometry)
{
    using Analysis=std::conditional_t<D==2,geo::GeometryAnalysis2D<Type>,geo::GeometryAnalysis3D<Type>>;
    geo::test::checkGeometry(grid,[&](Analysis& result) { geometry.analysis(result); });
}

void regionsAndNesting()
{
    const auto physical=geo::Rectangle({0,0},{5,8}) | geo::Rectangle({3,0},{8,8});
    const auto grid=geo::test::makeRegionGrid(4);
    geo::test::checkTopology(grid); compareGeometry(grid,physical);
    bool third=false,coarseGap=false;
    for(std::size_t c=0; c<grid.cells().size(); ++c)
    {
        if(grid.cells().isAuxiliary(c)) continue;
        third=third || grid.level(c)==3;
        const auto p=grid.center(c); coarseGap=coarseGap || (grid.level(c)==0 && p[0]>3 && p[0]<5);
    }
    require(third && coarseGap,"Nested/disconnected refinement was not preserved");
}

void booleanWall()
{
    const geo::BoundingBox2D domain{{0,0},{8,8}};
    const auto physical=geo::Rectangle({0,0},{5,8}) | geo::Rectangle({3,0},{8,8});
    geo::MultilevelOptions options; options.maxLevel=3; options.analysisTileWidth=4;
    const auto region=geo::fluidRegion(geo::Rectangle({0,0},{4,4}) | geo::Rectangle({2,0},{6,4}));
    require(!region.nearSurface(geo::BoundingBox2D{{3.9,1.5},{4.1,2.5}},0.1,{1}),"Hidden Boolean Face selected for wall refinement");
    const auto exposed=geo::RefinementRule<2>::surface(geo::fluidRegion(physical),2,0.3,{0});
    geo::MultilevelGrid<2> wall(domain,1,options); wall.generate(physical,{exposed}); geo::test::checkTopology(wall);
    compareGeometry(wall,physical);
}

} // namespace

int main()
{
    try
    {
        omp_set_dynamic(0); omp_set_num_threads(4);
        regionsAndNesting();
        std::cout << "[PASS] regions\n";
        booleanWall();
        std::cout << "[PASS] boolean_wall\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "[FAIL] test_multilevel_regions: " << error.what() << '\n';
        return 1;
    }
}
