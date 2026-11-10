#include "test_multilevel_checks.hpp"

#include <functional>
#include <iostream>
#include <stdexcept>

namespace geo=ntic::lbm::geometry;
namespace lat=ntic::lbm::lattice;

namespace
{
using geo::test::require;

template<typename Function>
void rejects(Function function)
{
    bool rejected=false;
    try { function(); } catch(const std::exception&) { rejected=true; }
    require(rejected,"Invalid input was accepted");
}

void topologyEdits()
{
    const geo::BoundingBox2D domain{{0,0},{4,4}}; const geo::Rectangle physical(domain.min,domain.max);
    geo::MultilevelOptions options; options.maxLevel=2; options.analysisTileWidth=2;
    geo::MultilevelGrid<2> grid(domain,1,options); grid.generate(physical);
    grid.refine(0); rejects([&] { grid.cells(); }); grid.finalize(physical); geo::test::checkTopology(grid);
    require(grid.nodes()[0].cell==geo::invalidIndex,"Refined parent still has data");
    grid.derefine(0); grid.finalize(physical); geo::test::checkTopology(grid);
    require(grid.cells().size()==16 && grid.interfaces().blocks.empty(),"derefine failed to restore coarse leaves");
}

void invalidInput()
{
    const geo::BoundingBox2D domain{{0,0},{4,4}}; const geo::Rectangle physical(domain.min,domain.max);
    geo::MultilevelOptions options; options.maxLevel=2; options.analysisTileWidth=2;
    geo::MultilevelGrid<2> grid(domain,1,options); grid.generate(physical);
    rejects([&] { grid.derefine(0); });
    rejects([&] { geo::MultilevelGrid<2> invalid(domain,0); });
    rejects([&] { auto invalid=geo::RefinementRule<2>::surface(geo::fluidRegion(physical),1,0.1,{4}); (void)invalid; });
    options.maxNodes=16; geo::MultilevelGrid<2> limited(domain,1,options);
    rejects([&] { limited.refine(0); });
    rejects([&] { grid.generate(physical,{geo::RefinementRule<2>::volume(geo::volumeRegion(physical),3)}); });
}

} // namespace

int main()
{
    try
    {
        omp_set_dynamic(0); omp_set_num_threads(4);
        topologyEdits();
        std::cout << "[PASS] edits\n";
        invalidInput();
        std::cout << "[PASS] invalid_input\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "[FAIL] test_multilevel_edits: " << error.what() << '\n';
        return 1;
    }
}
