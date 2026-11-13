#include "test_partition_checks.hpp"
#include <iostream>

namespace geo=ntic::lbm::geometry;
namespace part=ntic::lbm::mpi;

int main()
{
    try
    {
        // L-shaped and disconnected wet regions use one bounding box; Dry
        // padding does not allocate fields or create artificial local messages.
        part::PartitionGraph graph;graph.dimension=2;
        graph.blocks={{{0,0,0},{1,1,1},0,1},{{1,0,0},{2,1,1},0,1},{{0,1,0},{1,2,1},0,1}};
        graph.edges={{0,1,1,8,1},{0,2,1,8,1}};
        graph.blocks.push_back({{20,20,0},{21,21,1},0,0});
        part::PartitionOptions options;options.hardware.intraRankLatency=2;options.hardware.intraRankBandwidth=8;
        const auto plan=part::partitionGraph(graph,options);
        part::test::require(plan.owner.back()==-1 && plan.cellBlocks.back()==part::invalidCellIndex && plan.statistics.totalWeight==3,
                            "Zero-work background atom acquired ownership or load");
        part::test::require(plan.blocks.size()==1 && plan.blockExchanges.empty(),"L-shaped wet region did not use one padded box");
        part::test::require(plan.blocks[0].lower[0]==0 && plan.blocks[0].lower[1]==0 &&
                            plan.blocks[0].upper[0]==2 && plan.blocks[0].upper[1]==2,"L-shaped box bounds differ");
        part::test::require(plan.statistics.communicationTimes[0]==0 && plan.statistics.maxEstimatedTime==3,
                            "Single-box internal dependencies were charged as block exchange");
        part::PartitionGraph disconnected;disconnected.dimension=2;
        disconnected.blocks={{{0,0,0},{1,1,1},0,1},{{3,0,0},{4,1,1},0,1}};
        disconnected.edges={{0,1,1,16,2}};
        const auto nonlocal=part::partitionGraph(disconnected,options);
        part::test::require(nonlocal.blocks.size()==1 && nonlocal.blocks[0].upper[0]==4 &&
                            nonlocal.blockExchanges.empty() && nonlocal.statistics.communicationTimes[0]==0,
                            "Disconnected wet regions lost their single padded box");
        const auto physical=geo::Rectangle({0,0},{4,2})+geo::Rectangle({0,2},{2,4});
        geo::MultilevelOptions gridOptions;gridOptions.maxLevel=0;
        geo::MultilevelGrid<2> grid(physical.bounds,1,gridOptions);grid.generate(physical);
        const auto actual=part::partition(grid,options);
        const auto local=part::localPartition(grid,actual,0);
        part::test::require(local.intraRankExchanges.empty() && local.peers.empty() && local.globalCells.size()==local.ownedCount,
                            "Single rank allocated cross-block exchange or ghosts");
        for(std::size_t c=0;c<local.ownedCount;++c)
            for(auto neighbour:grid.adjacentRealCells(local.globalCells[c]))
                part::test::require(local.localIndex(neighbour)<local.ownedCount,"Internal cell dependency is inaccessible");
        part::test::checkLocalPartitions(grid,actual);
        // A remote dependency uses same-node or cross-node hardware parameters instead.
        options.ranks=2; options.hardware.nodeIds={0,0}; options.hardware.intraNodeLatency=4; options.hardware.intraNodeBandwidth=8;
        const auto sameNode=part::partitionGraph(graph,options);
        options.hardware.nodeIds={0,1}; options.hardware.interNodeLatency=8; options.hardware.interNodeBandwidth=4;
        const auto crossNode=part::partitionGraph(graph,options);
        part::test::require(sameNode.statistics.interRankCommunicationTimes[0]>0 &&
                            crossNode.statistics.interRankCommunicationTimes[0]>sameNode.statistics.interRankCommunicationTimes[0],
                            "Node topology cost not applied");
        std::cout << "[PASS] single_box_dependencies_and_node_costs\n";
        return 0;
    }
    catch(const std::exception& error) { std::cerr << "[FAIL] intra-rank: " << error.what() << '\n'; return 1; }
}
