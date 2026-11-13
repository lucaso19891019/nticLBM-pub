#include "test_multilevel_checks.hpp"

#include <functional>
#include <iostream>
#include <stdexcept>

namespace geo=ntic::lbm::geometry;
namespace lat=ntic::lbm::lattice;

namespace
{


template<std::size_t D,lat::LatticeType Type,typename Geometry>
void compareGeometry(const geo::MultilevelGrid<D,Type>& grid,const Geometry& geometry)
{
    using Analysis=std::conditional_t<D==2,geo::GeometryAnalysis2D<Type>,geo::GeometryAnalysis3D<Type>>;
    checkGeometry(grid,[&](Analysis& result) { geometry.analysis(result); });
}

template<lat::LatticeType Type>
void threeDimensions()
{
    const geo::BoundingBox domain{{0,0,0},{4,4,4}};
    const geo::Box physical(domain.min,domain.max);
    const geo::Sphere selection({2,2,2},0.4);
    geo::MultilevelOptions options; options.maxLevel=2; options.analysisTileWidth=2;
    geo::MultilevelGrid<3,Type> grid(domain,1,options);
    grid.generate(physical,{geo::RefinementRule<3>::volume(geo::volumeRegion(selection),2)});
    checkTopology(grid); compareGeometry(grid,physical);
}

} // namespace

int main()
{
    try
    {
        omp_set_dynamic(0); omp_set_num_threads(4);
        threeDimensions<lat::LatticeType::D3Q15>();
        std::cout << "[PASS] d3q15\n";
        threeDimensions<lat::LatticeType::D3Q19>();
        std::cout << "[PASS] d3q19\n";
        threeDimensions<lat::LatticeType::D3Q27>();
        std::cout << "[PASS] d3q27\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "[FAIL] test_multilevel_lattice: " << error.what() << '\n';
        return 1;
    }
}
