#include "test_partition_checks.hpp"
#include <iostream>

namespace geo=ntic::lbm::geometry;
namespace part=ntic::lbm::mpi;

int main()
{
    try
    {
        omp_set_num_threads(4);
        // At h=1 this box has 6^3 wet cells. Its 4^3 interior cells
        // leave 152 boundary cells, although the six Face lists have 216 entries.
        const geo::Box physical({0,0,0},{6,6,6});
        geo::MultilevelOptions gridOptions; gridOptions.maxLevel=0;
        geo::MultilevelGrid<3> grid(physical.bounds,1,gridOptions); grid.generate(physical);
        constexpr std::size_t wetCells=216,interiorCells=64,boundaryCells=152,faceRecords=216;
        requirePartition(grid.cells().size()==wetCells,"Boundary fixture wet-cell count differs");
        const auto countedBoundary=std::count(grid.cells().cellType.begin(),grid.cells().cellType.end(),geo::CellType::Boundary);
        const auto countedInterior=std::count(grid.cells().cellType.begin(),grid.cells().cellType.end(),geo::CellType::Interior);
        requirePartition(static_cast<std::size_t>(countedBoundary)==boundaryCells &&
                            static_cast<std::size_t>(countedInterior)==interiorCells,"Boundary fixture layer counts differ");
        const auto& boundary=grid.boundaries();
        requirePartition(boundary.nBoundaries==6 && boundary.cells.size()==faceRecords,"Boundary fixture Face counts differ");
        std::vector<std::size_t> memberships(wetCells,0);
        for(auto cell:boundary.cells) ++memberships[cell];
        requirePartition(std::any_of(memberships.begin(),memberships.end(),[](auto n) { return n>1; }),
                            "Boundary fixture has no multi-Face cell");

        part::PartitionOptions options; options.wetWeight=2; options.boundaryExtraWeight=3;
        const auto weighted=part::partition(grid,options);
        // Charge ordinary work for 216 unique wet cells and extra work for
        // all 216 Face memberships, including repeated boundary-cell indices.
        requirePartition(weighted.statistics.totalWeight==2*wetCells+3*faceRecords,
                            "Boundary work is not ordinary wet work plus stored Face access costs");
        requirePartition(weighted.statistics.totalWeight>2*wetCells+3*boundaryCells,
                            "Multi-Face cells did not receive multiple boundary surcharges");
        checkLocalPartitions(grid,weighted);

        options.boundaryExtraWeight=0;
        const auto noExtra=part::partition(grid,options);
        requirePartition(noExtra.statistics.totalWeight==2*wetCells,"Zero surcharge changed ordinary wet work");
        requirePartition(weighted.statistics.totalWeight-noExtra.statistics.totalWeight==3*faceRecords,
                            "Boundary surcharge does not follow all Face memberships");
        options.wetWeight=4; options.boundaryExtraWeight=3;
        const auto moreWet=part::partition(grid,options);
        requirePartition(moreWet.statistics.totalWeight-weighted.statistics.totalWeight==2*wetCells,
                            "Ordinary wet work was omitted or duplicated at boundary cells");
        std::cout << "[PASS] boundary_face_membership_costs\n";
        return 0;
    }
    catch(const std::exception& error) { std::cerr << "[FAIL] boundary costs: " << error.what() << '\n'; return 1; }
}
