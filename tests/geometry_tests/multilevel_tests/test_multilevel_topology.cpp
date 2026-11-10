#include "test_multilevel_checks.hpp"

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

void slab()
{
    const geo::BoundingBox2D domain{{0,0},{8,8}};
    const geo::Rectangle physical(domain.min,domain.max),selection({2.3,0},{6.7,8});
    geo::MultilevelOptions options; options.maxLevel=3; options.analysisTileWidth=2;
    const auto rule=geo::RefinementRule<2>::volume(geo::volumeRegion(selection),1);
    geo::MultilevelGrid<2> grid(domain,1,options); grid.generate(physical,{rule});
    geo::test::checkTopology(grid); compareGeometry(grid,physical);
    require(grid.cells().size()==248,"Unexpected slab storage count");
    require(grid.interfaces().blocks.size()==16,"Slab does not have one coarse layer on each side");
    std::size_t fine=0;
    for(std::size_t c=0; c<grid.cells().size(); ++c)
    {
        const auto p=grid.coordinate(c);
        if(grid.cells().isAuxiliary(c)) require(p[0]==2 || p[0]==3 || p[0]==14 || p[0]==15,"Auxiliary band exceeds two fine layers");
        else if(grid.cells().role[c]&geo::FineInterface)
        { ++fine; require(p[0]==4 || p[0]==5 || p[0]==12 || p[0]==13,"Wrong two-layer fine interface band"); }
    }
    require(fine==64,"Wrong fine interface cell count");
}

void sparseDryStorage()
{
    const geo::BoundingBox2D domain{{-4,-4},{4,4}};
    const geo::Circle physical({0,0},2);
    geo::MultilevelOptions options; options.maxLevel=2; options.analysisTileWidth=2;
    geo::MultilevelGrid<2> grid(domain,1,options);
    grid.generate(physical,{geo::RefinementRule<2>::surface(geo::fluidRegion(physical),2,0.5)});
    geo::test::checkTopology(grid); compareGeometry(grid,physical);
    std::size_t dry=0;
    for(const auto& node:grid.nodes()) dry+=node.leaf() && node.cell==geo::invalidIndex;
    require(dry>0,"Dry topology-only leaves were not exercised");

    // maxCells applies to retained slots, not temporary dry candidates.
    options.maxCells=16;
    geo::MultilevelGrid<2> coarse(domain,1,options); coarse.generate(physical);
    geo::test::checkTopology(coarse); compareGeometry(coarse,physical);
    require(coarse.cells().size()<=16,"Dry candidates consume the cell capacity");
    const geo::Circle absent({20,20},1);
    coarse.finalize(absent);
    geo::test::checkTopology(coarse); compareGeometry(coarse,absent);
    require(coarse.cells().size()==0 && coarse.interfaces().blocks.empty() && coarse.boundaries().cells.empty(),
            "Fully dry domain retains cell data or stale connectivity");
    coarse.finalize(physical);
    geo::test::checkTopology(coarse); compareGeometry(coarse,physical);
}

} // namespace

int main()
{
    try
    {
        omp_set_dynamic(0); omp_set_num_threads(4);
        slab();
        std::cout << "[PASS] slab\n";
        sparseDryStorage();
        std::cout << "[PASS] sparse_dry_storage\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "[FAIL] test_multilevel_topology: " << error.what() << '\n';
        return 1;
    }
}
