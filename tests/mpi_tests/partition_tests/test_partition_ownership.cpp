#include "test_partition_checks.hpp"
#include <iostream>

namespace geo=ntic::lbm::geometry;
namespace part=ntic::lbm::mpi;

int main()
{
    try
    {
        omp_set_num_threads(4);
        const geo::Rectangle physical({0,0},{8,8});
        geo::MultilevelOptions gridOptions; gridOptions.maxLevel=0;
        geo::MultilevelGrid<2> uniform(physical.bounds,1,gridOptions); uniform.generate(physical);
        part::PartitionOptions options; options.ranks=4; options.boundaryExtraWeight=0;
        const auto equal=part::partition(uniform,options);
        checkOwnership(uniform,equal);
        for(auto count:equal.statistics.realCells) requirePartition(count==16,"Equal-cost uniform partition is unbalanced");
        requirePartition(equal.statistics.maxToMean==1,"Uniform weight statistics incorrect");
        std::cout << "[PASS] uniform_ownership\n";

        const auto grid=makePartitionAnnulus();
        const auto sourceNodes=grid.cells().sourceNode;
        options.wetWeight=2; options.boundaryExtraWeight=3;
        options.auxiliaryWeight=0.5; options.interfaceWeight=0; options.levelWeights={1,2,4};
        const auto weighted=part::partition(grid,options);
        checkOwnership(grid,weighted);
        std::vector<double> expected(options.ranks,0);
        std::vector<std::size_t> memberships(grid.cells().size(),0);
        for(auto cell:grid.boundaries().cells) ++memberships[cell];
        double boundaryExtra=0;
        std::size_t boundaryCells=0;
        for(std::size_t c=0; c<grid.cells().size(); ++c)
        {
            const bool boundary=!grid.cells().isAuxiliary(c) && grid.cells().cellType[c]==geo::CellType::Boundary;
            const double extra=boundary ? 3*std::max(std::size_t{1},memberships[c]) : 0;
            const double base=grid.cells().isAuxiliary(c) ? 0.5 : 2+extra;
            const auto level=grid.level(c)-(grid.cells().isAuxiliary(c) ? 1 : 0);
            expected[weighted.owner[c]]+=base*std::ldexp(1.0,static_cast<int>(level))*options.levelWeights[level];
            if(boundary)
            {
                ++boundaryCells;
                boundaryExtra+=extra*std::ldexp(1.0,static_cast<int>(level))*options.levelWeights[level];
            }
        }
        requirePartition(expected==weighted.statistics.rankWeights,"Weighted load accounting differs");
        requirePartition(boundaryCells>0,"Boundary load case contains no boundary cells");
        auto noExtra=options; noExtra.boundaryExtraWeight=0;
        const auto baseline=part::partition(grid,noExtra);
        requirePartition(weighted.statistics.totalWeight-baseline.statistics.totalWeight==boundaryExtra,
                            "Boundary extra work does not follow stored Face memberships");
        requirePartition(grid.cells().sourceNode==sourceNodes,"Partition reordered geometry storage");
        std::cout << "[PASS] weighted_multilevel_ownership\n";
        return 0;
    }
    catch(const std::exception& error) { std::cerr << "[FAIL] ownership: " << error.what() << '\n'; return 1; }
}
