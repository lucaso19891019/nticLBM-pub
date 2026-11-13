#pragma once

#include "partition.hpp"
#include <set>
#include <stdexcept>

namespace ntic::lbm::mpi::test
{
inline void require(bool condition,const char* message)
{ if(!condition) throw std::runtime_error(message); }

template<typename Function>
void rejects(Function function)
{
    bool rejected=false;
    try { function(); } catch(const std::exception&) { rejected=true; }
    require(rejected,"Invalid partition input was accepted");
}

inline geometry::MultilevelGrid<2> annulus()
{
    const auto physical=geometry::Circle({0,0},4)-geometry::Circle({0,0},1);
    geometry::MultilevelOptions options; options.maxLevel=2;
    geometry::MultilevelGrid<2> grid(physical.bounds,1,options);
    grid.generate(physical,{geometry::RefinementRule<2>::surface(geometry::fluidRegion(physical),2,0.3)});
    return grid;
}

template<typename Grid>
void checkOwnership(const Grid& grid,const PartitionPlan& plan)
{
    require(plan.owner.size()==grid.cells().size(),"Incomplete owner map");
    ntic::lbm::mpi::detail::validatePlan(grid,plan);
    if(!plan.owner.empty())
    {
        require(plan.blocks.size()==static_cast<std::size_t>(plan.ranks()),"Rank does not own exactly one box");
        for(auto count:plan.statistics.blocksPerRank) require(count==1,"Invalid per-rank box count");
    }
    std::vector<unsigned> seen(grid.cells().size(),0);
    for(int rank=0; rank<plan.ranks(); ++rank)
        for(auto i=plan.offsets[rank]; i<plan.offsets[rank+1]; ++i)
        {
            const auto c=plan.cells[i];
            require(c<seen.size() && seen[c]++==0 && plan.owner[c]==rank,"Duplicated/misowned cell");
            if(grid.cells().isAuxiliary(c))
                require(plan.owner[grid.nodes()[grid.cells().sourceNode[c]].cell]==rank,"Auxiliary family crosses ownership");
        }
    require(std::all_of(seen.begin(),seen.end(),[](auto count) { return count==1; }),"Omitted cell");
}

template<typename Grid>
void checkLocalPartitions(const Grid& grid,const PartitionPlan& plan)
{
    std::vector<LocalPartition> locals;
    for(int rank=0; rank<plan.ranks(); ++rank) locals.push_back(localPartition(grid,plan,rank));
    for(const auto& local:locals)
    {
        std::set<GlobalCellIndex> seen;
        for(std::size_t c=0; c<local.globalCells.size(); ++c)
        {
            const auto global=local.globalCells[c];
            require(seen.insert(global).second,"Duplicate local cell");
            require(local.localIndex(global)==c,"Global/local lookup mismatch");
            require((c<local.ownedCount)==(plan.owner[global]==local.rank),"Owned/ghost ordering mismatch");
            if(c<local.ownedCount && !grid.cells().isAuxiliary(global))
                for(auto neighbour:grid.adjacentRealCells(global))
                    require(local.localIndex(neighbour)!=invalidCellIndex,"Missing real neighbour ghost");
        }
        require(local.localIndex(invalidCellIndex)==invalidCellIndex,"Missing lookup did not return invalid index");
        require(local.cellBlocks.size()==local.globalCells.size(),"Missing local block indices");
        std::set<std::pair<GlobalCellIndex,GlobalCellIndex>> internal;
        for(const auto& exchange:local.intraRankExchanges)
        {
            require(exchange.a<exchange.b && exchange.b<plan.blocks.size(),"Invalid internal block pair");
            require(plan.blocks[exchange.a].owner==local.rank && plan.blocks[exchange.b].owner==local.rank,
                    "Internal block exchange crosses rank ownership");
            for(const auto& pair:exchange.dependencies)
            {
                require(pair.first<local.ownedCount && pair.second<local.ownedCount,"Internal exchange refers to a ghost");
                require(local.cellBlocks[pair.first]==exchange.a && local.cellBlocks[pair.second]==exchange.b,
                        "Internal dependency block indices differ");
                const auto a=local.globalCells[pair.first],b=local.globalCells[pair.second];
                require(internal.insert(std::minmax(a,b)).second,"Duplicate internal dependency");
            }
        }
        for(std::size_t c=0; c<local.ownedCount; ++c)
            if(!grid.cells().isAuxiliary(local.globalCells[c]))
                for(auto neighbour:grid.adjacentRealCells(local.globalCells[c]))
                    if(plan.owner[neighbour]==local.rank && plan.cellBlocks[neighbour]!=local.cellBlocks[c])
                        require(internal.count(std::minmax(local.globalCells[c],neighbour))!=0,"Missing intra-rank neighbour dependency");
        for(const auto& peer:local.peers)
        {
            const auto& remote=locals[peer.rank];
            const auto counterpart=std::find_if(remote.peers.begin(),remote.peers.end(),[&](const auto& p) { return p.rank==local.rank; });
            require(counterpart!=remote.peers.end(),"Missing reverse peer");
            require(peer.send.size()==counterpart->receive.size() && peer.receive.size()==counterpart->send.size(),"Peer sizes differ");
            for(std::size_t i=0; i<peer.send.size(); ++i)
            {
                require(peer.send[i]<local.ownedCount,"Sending a ghost as owner");
                require(local.globalCells[peer.send[i]]==remote.globalCells[counterpart->receive[i]],"Send/receive cell order differs");
            }
            for(std::size_t i=0; i<peer.receive.size(); ++i)
            {
                require(peer.receive[i]>=local.ownedCount,"Receive target is not a ghost");
                require(local.globalCells[peer.receive[i]]==remote.globalCells[counterpart->send[i]],"Receive/send cell order differs");
            }
        }
        for(const auto& block:local.interfaces)
        {
            require(block.coarse<local.globalCells.size(),"Missing local coarse interface");
            require(block.auxiliary.size()==(std::size_t{1}<<Grid::dimension),"Incomplete local auxiliary family");
            for(auto c:block.auxiliary) require(c<local.globalCells.size(),"Missing local auxiliary cell");
            for(auto c:block.fine) require(c<local.globalCells.size(),"Missing local fine interface");
        }
        const auto& global=grid.boundaries(); const auto& boundary=local.boundaries;
        constexpr auto links=lattice::Model<Grid::latticeType>::nStencils-1;
        require(boundary.offsets.size()==global.offsets.size(),"Local Face IDs changed");
        require(boundary.q.size()==boundary.cells.size()*links,"Local Q dimensions changed");
        for(std::size_t face=0; face<global.nBoundaries; ++face)
        {
            auto entry=boundary.offsets[face];
            for(auto i=global.offsets[face]; i<global.offsets[face+1]; ++i)
            {
                const auto cell=local.localIndex(global.cells[i]);
                if(cell==invalidCellIndex) continue;
                require(entry<boundary.offsets[face+1] && boundary.cells[entry]==cell,"Local boundary membership differs");
                for(std::size_t q=0; q<links; ++q) require(boundary.q[entry*links+q]==global.q[i*links+q],"Partition changed Q");
                ++entry;
            }
            require(entry==boundary.offsets[face+1],"Extra local boundary entry");
        }
    }
}
} // namespace ntic::lbm::mpi::test
