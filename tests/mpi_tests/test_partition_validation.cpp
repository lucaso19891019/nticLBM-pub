#include "test_partition_checks.hpp"
#include <iostream>

namespace geo=ntic::lbm::geometry;
namespace part=ntic::lbm::mpi;

int main()
{
    try
    {
        omp_set_num_threads(4);
        const auto grid=part::test::annulus(); part::PartitionOptions options; options.ranks=4;
        options.graphThreads=1; const auto serial=part::partition(grid,options);
        options.graphThreads=4; const auto parallel=part::partition(grid,options);
        part::test::require(serial.owner==parallel.owner && serial.cells==parallel.cells && serial.offsets==parallel.offsets &&
                            serial.statistics.rankWeights==parallel.statistics.rankWeights,"Nondeterministic partition");
        std::cout << "[PASS] deterministic_partition\n";
        part::test::rejects([&] { auto bad=options; bad.ranks=0; part::partition(grid,bad); });
        part::test::rejects([&] { auto bad=options; bad.boundaryWeight=0; part::partition(grid,bad); });
        part::test::rejects([&] { auto bad=options; bad.auxiliaryWeight=std::numeric_limits<double>::infinity(); part::partition(grid,bad); });
        part::test::rejects([&] { auto bad=options; bad.levelWeights={1}; part::partition(grid,bad); });
        part::test::rejects([&] { auto bad=options; bad.hardware.computeRates={1,2}; part::partition(grid,bad); });
        part::test::rejects([&] { auto bad=options; bad.hardware.computeRates={1,1,0,1}; part::partition(grid,bad); });
        part::test::rejects([&] { part::localPartition(grid,serial,4); });
        part::test::rejects([&] { auto bad=serial; bad.owner.pop_back(); part::localPartition(grid,bad,0); });
        part::test::rejects([&] { auto bad=serial; bad.cells[0]=bad.cells[1]; part::localPartition(grid,bad,0); });
        std::cout << "[PASS] invalid_input\n";
        const geo::Rectangle physical({0,0},{1,1});
        geo::MultilevelOptions gridOptions; gridOptions.maxLevel=0;
        geo::MultilevelGrid<2> small(physical.bounds,1,gridOptions); small.generate(physical);
        auto sparse=part::partition(small,options);
        part::test::checkOwnership(small,sparse); part::test::checkLocalPartitions(small,sparse);
        const geo::Circle absent({10,10},1); small.finalize(absent);
        auto empty=part::partition(small,options);
        part::test::require(empty.statistics.totalWeight==0 && empty.statistics.maxToMean==0,"Empty weight statistics invalid");
        part::test::checkLocalPartitions(small,empty);
        std::cout << "[PASS] empty_ranks_and_empty_domain\n";
        return 0;
    }
    catch(const std::exception& error) { std::cerr << "[FAIL] validation: " << error.what() << '\n'; return 1; }
}
