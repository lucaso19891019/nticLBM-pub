#include "test_partition_checks.hpp"
#include <iostream>

namespace geo=ntic::lbm::geometry;
namespace part=ntic::lbm::mpi;

int main()
{
    try
    {
        omp_set_num_threads(4);
        const auto annulus=part::test::annulus();
        part::PartitionOptions options; options.ranks=4;
        part::test::checkLocalPartitions(annulus,part::partition(annulus,options));
        std::cout << "[PASS] 2d_ghosts_interfaces_q\n";
        const auto physical=geo::Box({0,0,0},{4,4,4})-geo::Sphere({2,2,2},1);
        geo::MultilevelOptions gridOptions; gridOptions.maxLevel=1;
        geo::MultilevelGrid<3> grid(physical.bounds,1,gridOptions);
        grid.generate(physical,{geo::RefinementRule<3>::surface(geo::fluidRegion(physical),1,0.3,{6})});
        const auto plan=part::partition(grid,options);
        part::test::checkOwnership(grid,plan); part::test::checkLocalPartitions(grid,plan);
        std::cout << "[PASS] 3d_ghosts_interfaces_q\n";
        options.ranks=1;
        const auto single=part::localPartition(grid,part::partition(grid,options),0);
        part::test::require(single.peers.empty() && single.ownedCount==grid.cells().size() && single.globalCells.size()==single.ownedCount,
                            "Single rank has unnecessary ghosts or peers");
        std::cout << "[PASS] single_rank\n";
        return 0;
    }
    catch(const std::exception& error) { std::cerr << "[FAIL] connectivity: " << error.what() << '\n'; return 1; }
}
