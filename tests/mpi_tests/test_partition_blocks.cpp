#include "test_partition_checks.hpp"
#include <iostream>

namespace geo=ntic::lbm::geometry;
namespace part=ntic::lbm::mpi;

namespace
{
template<typename Grid>
void checkBlocks(const Grid& grid,const part::PartitionPlan& plan)
{
    constexpr auto dimensions=Grid::dimension;
    long double volume=0,expected=1;
    for(std::size_t d=0; d<dimensions; ++d)
        expected*=(grid.domain().max[d]-grid.domain().min[d])/plan.blockSpacing;
    for(std::size_t b=0; b<plan.blocks.size(); ++b)
    {
        const auto& block=plan.blocks[b]; long double size=1;
        part::test::require(block.owner>=0 && block.owner<plan.ranks(),"Invalid spatial owner");
        for(std::size_t d=0; d<dimensions; ++d)
        {
            part::test::require(block.upper[d]>block.lower[d],"Degenerate spatial block");
            part::test::require(block.upper[d]<=(grid.domain().max[d]-grid.domain().min[d])/plan.blockSpacing,
                                "Spatial block exceeds domain");
            size*=block.upper[d]-block.lower[d];
        }
        volume+=size;
        for(std::size_t a=0; a<b; ++a)
        {
            bool separated=false;
            for(std::size_t d=0; d<dimensions; ++d)
                separated=separated || block.upper[d]<=plan.blocks[a].lower[d] || plan.blocks[a].upper[d]<=block.lower[d];
            part::test::require(separated,"Spatial blocks overlap");
        }
    }
    part::test::require(volume==expected,"Spatial tiling omits part of the domain");
    for(std::size_t c=0; c<grid.cells().size(); ++c)
    {
        const auto p=grid.coordinate(c); const auto level=grid.level(c);
        const auto common=std::max(level,plan.blockLevel); std::size_t matches=0;
        for(const auto& block:plan.blocks)
        {
            bool contains=true;
            for(std::size_t d=0; d<dimensions; ++d)
            {
                const auto low=static_cast<std::uint64_t>(p[d])<<(common-level);
                const auto high=(static_cast<std::uint64_t>(p[d])+1)<<(common-level);
                contains=contains && (static_cast<std::uint64_t>(block.lower[d])<<(common-plan.blockLevel))<=low &&
                    (static_cast<std::uint64_t>(block.upper[d])<<(common-plan.blockLevel))>=high;
            }
            if(contains) { ++matches; part::test::require(block.owner==plan.owner[c],"Cell owner differs from spatial block"); }
        }
        part::test::require(matches==1,"A cell is cut by a plane or has no unique spatial block");
    }
}
}

int main()
{
    try
    {
        omp_set_num_threads(4);
        const auto grid=part::test::annulus();
        part::PartitionOptions options; options.ranks=4; options.blockLevel=1;
        const auto plan=part::partition(grid,options);
        checkBlocks(grid,plan); part::test::checkOwnership(grid,plan); part::test::checkLocalPartitions(grid,plan);
        auto explicitFrequency=options; explicitFrequency.levelWeights={1,1,1};
        const auto explicitPlan=part::partition(grid,explicitFrequency);
        part::test::require(plan.owner==explicitPlan.owner && plan.statistics.totalWeight==explicitPlan.statistics.totalWeight,
                            "Default level frequencies are not 2^level");
        std::cout << "[PASS] rectangular_tiling_and_level_frequency\n";

        const geo::Rectangle physical({0,0},{16,16});
        geo::MultilevelOptions gridOptions; gridOptions.maxLevel=0;
        geo::MultilevelGrid<2> uniform(physical.bounds,1,gridOptions); uniform.generate(physical);
        options.ranks=2; options.blockLevel=0;
        options.hardware.computeRates={1,2}; options.hardware.nodeIds={0,1};
        options.hardware.interNodeLatency=0.01; options.hardware.interNodeBandwidth=1000;
        const auto hardware=part::partition(uniform,options); checkBlocks(uniform,hardware);
        part::test::require(hardware.statistics.rankWeights[1]>hardware.statistics.rankWeights[0],"Faster rank not assigned more work");
        part::test::require(hardware.statistics.maxToMean<1.15,"Hardware-normalized load is unbalanced");
        part::test::require(hardware.statistics.communicationTimes[0]>0,"Cross-node communication cost not evaluated");
        part::test::require(hardware.blocks.size()<256,"Adjacent rectangles were not merged");
        std::cout << "[PASS] hardware_capacity_and_block_merging\n";
        return 0;
    }
    catch(const std::exception& error) { std::cerr << "[FAIL] blocks: " << error.what() << '\n'; return 1; }
}
