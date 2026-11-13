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
    requirePartition(plan.blocks.size()==static_cast<std::size_t>(plan.ranks()),"Expected one box per rank");
    std::vector<unsigned> owners(plan.ranks(),0);
    for(std::size_t b=0; b<plan.blocks.size(); ++b)
    {
        const auto& block=plan.blocks[b];
        requirePartition(block.owner>=0 && block.owner<plan.ranks(),"Invalid spatial owner");
        requirePartition(owners[block.owner]++==0,"Rank owns multiple boxes");
        for(std::size_t d=0; d<dimensions; ++d)
        {
            requirePartition(block.upper[d]>block.lower[d],"Degenerate spatial block");
            requirePartition(block.upper[d]<=(grid.domain().max[d]-grid.domain().min[d])/plan.blockSpacing,
                                "Spatial block exceeds domain");
        }
        for(std::size_t a=0; a<b; ++a)
        {
            bool separated=false;
            for(std::size_t d=0; d<dimensions; ++d)
                separated=separated || block.upper[d]<=plan.blocks[a].lower[d] || plan.blocks[a].upper[d]<=block.lower[d];
            requirePartition(separated,"Spatial blocks overlap");
        }
    }
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
            if(contains) { ++matches; requirePartition(block.owner==plan.owner[c],"Cell owner differs from spatial block"); }
        }
        requirePartition(matches==1,"A cell is cut by a plane or has no unique spatial block");
    }
}

part::PartitionGraph squareRing()
{
    part::PartitionGraph graph;graph.dimension=2;
    std::array<std::size_t,64> ids;ids.fill(part::invalidCellIndex);
    for(unsigned y=0;y<8;++y) for(unsigned x=0;x<8;++x)
    {
        if(x>=2 && x<6 && y>=2 && y<6) continue;
        ids[x+8*y]=graph.blocks.size();
        graph.blocks.push_back({{x,y,0},{x+1,y+1,1},0,1});
    }
    for(unsigned y=0;y<8;++y) for(unsigned x=0;x<8;++x)
    {
        const auto a=ids[x+8*y];if(a==part::invalidCellIndex) continue;
        if(x+1<8 && ids[x+1+8*y]!=part::invalidCellIndex) graph.edges.push_back({a,ids[x+1+8*y],1,8,1});
        if(y+1<8 && ids[x+8*(y+1)]!=part::invalidCellIndex) graph.edges.push_back({a,ids[x+8*(y+1)],1,8,1});
    }
    return graph;
}

void checkGlobalSearch()
{
    auto ring=squareRing();part::PartitionOptions options;options.ranks=4;
    options.hardware.intraNodeBandwidth=8;
    const auto quadrants=part::partitionGraph(ring,options);
    for(const auto& box:quadrants.blocks)
        for(unsigned d=0;d<2;++d)
            requirePartition((box.lower[d]==0 || box.lower[d]==4) && box.upper[d]==box.lower[d]+4,
                                "Balanced square ring did not select central quadrants");
    requirePartition(quadrants.statistics.cutWeight==8,"Square-ring cut edge count differs");
    for(int rank=0;rank<4;++rank)
        requirePartition(quadrants.statistics.rankWeights[rank]==12 && quadrants.statistics.communicationTimes[rank]==4,
                            "Independent square-ring load/communication calculation differs");

    // Root/child cuts share the middle rank: moving both produces a
    // three-rank one-in/one-out compensation, unlike disjoint subtree moves.
    part::PartitionGraph chain;chain.dimension=2;
    for(unsigned x=0;x<9;++x) chain.blocks.push_back({{x,0,0},{x+1,1,1},0,1});
    for(std::size_t x=0;x<8;++x) chain.edges.push_back({x,x+1,1,8,1});
    auto chainOptions=options;chainOptions.ranks=3;
    const auto compensated=part::partitionGraph(chain,chainOptions);
    requirePartition(compensated.statistics.jointCutCandidates>0 && compensated.statistics.maxComputeTime==3,
                        "Three-rank ancestor/child compensation was not evaluated or harmed load balance");

    // METIS weights prefer four compact clusters. Physical byte costs instead
    // favour four horizontal rows: the two costly links then become local.
    part::PartitionGraph weighted;weighted.dimension=2;
    for(unsigned y=0;y<4;++y) for(unsigned x=0;x<4;++x)
        weighted.blocks.push_back({{x,y,0},{x+1,y+1,1},0,(x+y)%2 ? 1.5 : 0.5});
    for(unsigned y=0;y<4;++y) for(unsigned x=0;x<4;++x)
    {
        const std::size_t a=x+4*y;
        if(x<3) weighted.edges.push_back({a,a+1,x/2==(x+1)/2 ? 1000.0 : 1.0,
                                          x==1 && (y==0 || y==3) ? 1.0e6 : 0.0,1});
        if(y<3) weighted.edges.push_back({a,a+4,y/2==(y+1)/2 ? 1000.0 : 1.0,0,1});
    }
    options.candidates=1;options.refinementRounds=0;
    const auto initial=part::partitionGraph(weighted,options);
    options.refinementRounds=8;
    const auto refined=part::partitionGraph(weighted,options);
    requirePartition(refined.statistics.maxComputeTime<=initial.statistics.maxComputeTime*(1+options.loadSlack),
                        "Global refinement exceeded the initial load budget");
    const auto maxCommunication=[](const auto& plan)
    {return *std::max_element(plan.statistics.communicationTimes.begin(),plan.statistics.communicationTimes.end());};
    requirePartition(maxCommunication(refined)<maxCommunication(initial),"Grouped search did not improve physical interface cost");
    std::set<int> changed;
    for(std::size_t i=0;i<initial.owner.size();++i) if(initial.owner[i]!=refined.owner[i])
    {changed.insert(initial.owner[i]);changed.insert(refined.owner[i]);}
    requirePartition(changed.size()>=3 && refined.statistics.evaluatedCandidates>initial.statistics.evaluatedCandidates,
                        "Refinement did not exercise a multi-rank compensated change");
}

}

int main()
{
    try
    {
        omp_set_num_threads(4);
        const auto grid=makePartitionAnnulus();
        part::PartitionOptions options; options.ranks=4; options.blockLevel=1;
        const auto plan=part::partition(grid,options);
        checkBlocks(grid,plan); checkOwnership(grid,plan); checkLocalPartitions(grid,plan);
        auto explicitFrequency=options; explicitFrequency.levelWeights={1,1,1};
        const auto explicitPlan=part::partition(grid,explicitFrequency);
        requirePartition(plan.owner==explicitPlan.owner && plan.statistics.totalWeight==explicitPlan.statistics.totalWeight,
                            "Default level frequencies are not 2^level");
        std::cout << "[PASS] rectangular_wet_cover_and_level_frequency\n";

        const geo::Rectangle physical({0,0},{16,16});
        geo::MultilevelOptions gridOptions; gridOptions.maxLevel=0;
        geo::MultilevelGrid<2> uniform(physical.bounds,1,gridOptions); uniform.generate(physical);
        options.ranks=2; options.blockLevel=0;
        options.hardware.computeRates={1,2}; options.hardware.nodeIds={0,1};
        options.hardware.interNodeLatency=0.01; options.hardware.interNodeBandwidth=1000;
        const auto hardware=part::partition(uniform,options); checkBlocks(uniform,hardware);
        requirePartition(hardware.statistics.rankWeights[1]>hardware.statistics.rankWeights[0],"Faster rank not assigned more work");
        requirePartition(hardware.statistics.maxToMean<1.15,"Hardware-normalized load is unbalanced");
        requirePartition(hardware.statistics.communicationTimes[0]>0,"Cross-node communication cost not evaluated");
        requirePartition(hardware.blocks.size()==2,"Hardware partition must have one box per rank");
        checkGlobalSearch();
        std::cout << "[PASS] hardware_capacity_and_single_boxes\n";
        return 0;
    }
    catch(const std::exception& error) { std::cerr << "[FAIL] blocks: " << error.what() << '\n'; return 1; }
}
