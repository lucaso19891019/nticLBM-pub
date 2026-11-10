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
void checkTopology(const geo::MultilevelGrid<2>& grid);
void checkGeometry(
    const geo::MultilevelGrid<2>& grid,
    const std::function<void(geo::GeometryAnalysis2D<lat::LatticeType::D2Q9>&)>& analyze);
void checkTopology(const geo::MultilevelGrid<3>& grid);
void checkGeometry(
    const geo::MultilevelGrid<3>& grid,
    const std::function<void(geo::GeometryAnalysis3D<lat::LatticeType::D3Q19>&)>& analyze);
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

void curvedWall()
{
    const auto annulus=geo::Circle({0,0},4)-geo::Circle({0,0},1);
    geo::MultilevelOptions options; options.maxLevel=2; options.analysisTileWidth=3;
    geo::MultilevelGrid<2> grid(annulus.bounds,1,options);
    grid.generate(annulus,{geo::RefinementRule<2>::surface(geo::fluidRegion(annulus),2,0.4)});
    geo::test::checkTopology(grid); compareGeometry(grid,annulus);
}

void externalFlow()
{
    geo::MultilevelOptions options; options.maxLevel=2; options.analysisTileWidth=3;
    const geo::BoundingBox open{{-2,-2,-2},{2,2,2}};
    const geo::Sphere obstacle({0,0,0},0.75,geo::FlowType::External,&open);
    const geo::Cylinder zone({0,0,0},1.2,3,geo::Axis::Z);
    geo::MultilevelGrid<3> external(open,1,options);
    external.generate(obstacle,{geo::RefinementRule<3>::volume(geo::volumeRegion(zone),1)});
    geo::test::checkTopology(external); compareGeometry(external,obstacle);
}

void cylinderFaces()
{
    const geo::Cylinder zone({0,0,0},1.2,3,geo::Axis::Z);
    const auto caps=geo::RefinementRule<3>::surface(geo::fluidRegion(zone),1,0.1,{1,2});
    require(caps.intersects(geo::BoundingBox{{-0.1,-0.1,1.4},{0.1,0.1,1.5}}),"Cylinder cap selector missed cap");
    require(!caps.intersects(geo::BoundingBox{{-0.1,-0.1,-0.1},{0.1,0.1,0.1}}),"Cylinder cap selector selected deep interior");
}

} // namespace

int main()
{
    try
    {
        omp_set_dynamic(0); omp_set_num_threads(4);
        curvedWall();
        std::cout << "[PASS] curved_wall\n";
        externalFlow();
        std::cout << "[PASS] external_flow\n";
        cylinderFaces();
        std::cout << "[PASS] cylinder_faces\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "[FAIL] test_multilevel_primitive: " << error.what() << '\n';
        return 1;
    }
}
