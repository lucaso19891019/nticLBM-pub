#pragma once

#include "partition_types.hpp"
#include "geometry/multilevel/multilevel_grid.hpp"

#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>
#include <utility>
#include <exception>

namespace ntic::lbm::mpi
{

namespace detail
{
inline void positiveWeight(double weight)
{
    if(!(weight>0.0) || !std::isfinite(weight))
        throw std::invalid_argument("Partition weights must be positive and finite.");
}

template<typename Grid>
GlobalCellIndex auxiliarySource(const Grid& grid,GlobalCellIndex cell)
{
    const auto source=grid.nodes().at(grid.cells().sourceNode.at(cell)).cell;
    if(source>=grid.cells().size() || grid.cells().isAuxiliary(source) || !grid.cells().computes(source))
        throw std::logic_error("Auxiliary cell has no wet real coarse source.");
    return source;
}

template<typename Grid>
void validatePlan(const Grid& grid,const PartitionPlan& plan)
{
    if(plan.ranks()<1 || plan.owner.size()!=grid.cells().size() || plan.cells.size()!=plan.owner.size() ||
       plan.cellBlocks.size()!=plan.owner.size() ||
       plan.offsets.front()!=0 || plan.offsets.back()!=plan.owner.size())
        throw std::invalid_argument("Partition plan does not match the finalized grid.");
    for(std::size_t d=0;d<Grid::dimension;++d)
        if(plan.domainMin[d]!=grid.domain().min[d]) throw std::invalid_argument("Partition origin does not match the grid.");
    if(plan.blockLevel>30 || plan.dimension!=Grid::dimension ||
       plan.blockSpacing!=grid.spacing(plan.blockLevel))
        throw std::invalid_argument("Partition coordinate frame does not match the grid.");
    if(plan.owner.empty())
    {
        if(!plan.blocks.empty()) throw std::invalid_argument("Empty grid has spatial partition blocks.");
    }
    else if(plan.blocks.size()!=static_cast<std::size_t>(plan.ranks()))
        throw std::invalid_argument("Each rank must own exactly one spatial box.");
    std::vector<bool> rankBoxes(plan.ranks(),false);
    for(std::size_t b=0;b<plan.blocks.size();++b)
    {
        const auto& box=plan.blocks[b];
        if(box.owner<0 || box.owner>=plan.ranks() || rankBoxes[box.owner])
            throw std::invalid_argument("Repeated or invalid spatial box owner.");
        rankBoxes[box.owner]=true;
        for(std::size_t d=0;d<Grid::dimension;++d)
            if(box.lower[d]>=box.upper[d]) throw std::invalid_argument("Degenerate partition box.");
        for(std::size_t a=0;a<b;++a)
        {
            bool separated=false;
            for(std::size_t d=0;d<Grid::dimension;++d)
                separated=separated || box.upper[d]<=plan.blocks[a].lower[d] || plan.blocks[a].upper[d]<=box.lower[d];
            if(!separated) throw std::invalid_argument("Partition boxes overlap.");
        }
    }
    std::vector<std::uint8_t> seen(plan.owner.size(),0);
    for(int rank=0; rank<plan.ranks(); ++rank)
    {
        if(plan.offsets[rank]>plan.offsets[rank+1] || plan.offsets[rank+1]>plan.cells.size())
            throw std::invalid_argument("Invalid partition offsets.");
        for(auto i=plan.offsets[rank]; i<plan.offsets[rank+1]; ++i)
        {
            const auto c=plan.cells[i];
            if(c>=plan.owner.size() || seen[c]++ || plan.owner[c]!=rank)
                throw std::invalid_argument("Partition ownership is inconsistent.");
            if(plan.cellBlocks[c]>=plan.blocks.size() || plan.blocks[plan.cellBlocks[c]].owner!=rank)
                throw std::invalid_argument("Partition cell/block mapping is inconsistent.");
            const auto coordinate=grid.coordinate(c);const auto level=grid.level(c);
            const auto common=std::max(level,plan.blockLevel);
            const auto& box=plan.blocks[plan.cellBlocks[c]];
            for(std::size_t d=0;d<Grid::dimension;++d)
            {
                const auto lower=static_cast<std::uint64_t>(coordinate[d])<<(common-level);
                const auto upper=(static_cast<std::uint64_t>(coordinate[d])+1)<<(common-level);
                if((static_cast<std::uint64_t>(box.lower[d])<<(common-plan.blockLevel))>lower ||
                   (static_cast<std::uint64_t>(box.upper[d])<<(common-plan.blockLevel))<upper)
                    throw std::invalid_argument("Stored cell is not contained by its rank box.");
            }
        }
    }
    for(GlobalCellIndex c=0; c<plan.owner.size(); ++c)
        if(grid.cells().isAuxiliary(c) && plan.owner[c]!=plan.owner[auxiliarySource(grid,c)])
            throw std::invalid_argument("Auxiliary family is split from its coarse source.");
}
} // namespace detail

namespace detail
{
inline double frequency(unsigned level,const PartitionOptions&)
{
    return std::ldexp(1.0,static_cast<int>(level));
}

inline double levelCost(unsigned level,const PartitionOptions& options)
{
    if(options.levelWeights.empty()) return 1;
    if(level>=options.levelWeights.size()) throw std::invalid_argument("Level cost weights are incomplete.");
    return options.levelWeights[level];
}

inline void nonnegative(double value)
{
    if(value<0 || !std::isfinite(value)) throw std::invalid_argument("Partition cost must be finite and nonnegative.");
}

inline void validateOptions(const PartitionOptions& o)
{
    if(o.ranks<1 || o.blockLevel>30 || o.candidates<1 || o.candidates>64 || o.seed<0 ||
       o.seed>std::numeric_limits<int>::max()-o.candidates || o.graphThreads<1 || o.graphThreads>1024 || o.refinementRounds<0 || o.refinementRounds>256 ||
       o.searchWidth<2 || o.searchWidth>16 || o.cutCandidates<1 || o.cutCandidates>32 ||
       !(o.imbalanceTolerance>1) || !std::isfinite(o.imbalanceTolerance))
        throw std::invalid_argument("Invalid partition options.");
    if(o.ranks>std::numeric_limits<int>::max()/2 ||
       static_cast<std::size_t>(o.ranks)>std::numeric_limits<std::size_t>::max()/(static_cast<std::size_t>(o.cutCandidates)*8))
        throw std::length_error("Partition search budget overflow.");
    positiveWeight(o.wetWeight); nonnegative(o.boundaryExtraWeight);
    nonnegative(o.auxiliaryWeight); nonnegative(o.interfaceWeight);
    positiveWeight(o.streamingEdgeWeight); positiveWeight(o.interfaceEdgeWeight);
    nonnegative(o.streamingBytes); nonnegative(o.interfaceBytes); nonnegative(o.loadSlack);
    for(auto w:o.levelWeights) positiveWeight(w);
    const auto& h=o.hardware;
    if((!h.computeRates.empty() && h.computeRates.size()!=static_cast<std::size_t>(o.ranks)) ||
       (!h.nodeIds.empty() && h.nodeIds.size()!=static_cast<std::size_t>(o.ranks)))
        throw std::invalid_argument("Hardware profile must describe every rank.");
    for(auto w:h.computeRates) positiveWeight(w);
    for(auto n:h.nodeIds) if(n<0) throw std::invalid_argument("Negative hardware node ID.");
    nonnegative(h.intraNodeLatency); nonnegative(h.interNodeLatency);
    positiveWeight(h.intraNodeBandwidth); positiveWeight(h.interNodeBandwidth);
    nonnegative(h.intraRankLatency); positiveWeight(h.intraRankBandwidth);
}
} // namespace detail

// Axis-aligned topology boxes contain their own leaf cells and descendants.
// Their weighted graph, rather than a space-filling curve, is given to METIS.
// Zero-work atoms are excluded from ownership. METIS labels guide constrained
// aggregation into one box per rank; Dry padding has no coverage requirement.
template<std::size_t D,lattice::LatticeType Type>
PartitionPlan partition(const geometry::MultilevelGrid<D,Type>& grid,const PartitionOptions& options={})
{
    detail::validateOptions(options);
    const auto& cells=grid.cells();
    // Each stored Face/cell record represents boundary data access, even if
    // all of its Q links are invalid. Auxiliary slots use their separate cost.
    std::vector<std::uint32_t> boundaryMemberships(cells.size(),0);
    for(auto c:grid.boundaries().cells)
    {
        if(c>=cells.size()) throw std::logic_error("Boundary record references an invalid cell.");
        if(cells.isAuxiliary(c) || cells.cellType[c]!=geometry::CellType::Boundary) continue;
        auto& count=boundaryMemberships[c];
        if(count==std::numeric_limits<std::uint32_t>::max()) throw std::length_error("Boundary membership count overflow.");
        ++count;
    }
    if(!(grid.spacing(options.blockLevel)>0) || !std::isfinite(grid.spacing(options.blockLevel)))
        throw std::invalid_argument("Atomic-block coordinate spacing is not representable.");
    const auto& nodes=grid.nodes();
    std::vector<std::size_t> nodeBlock(nodes.size(),invalidCellIndex);
    PartitionGraph graph; graph.dimension=D;
    for(std::size_t n=0; n<nodes.size(); ++n)
    {
        const auto& node=nodes[n];
        if(!node.active || node.level>options.blockLevel || (!node.leaf() && node.level!=options.blockLevel)) continue;
        PartitionBlock block;
        for(std::size_t d=0; d<3; ++d)
        {
            const unsigned shift=options.blockLevel-node.level;
            const auto lower=d<D ? static_cast<std::uint64_t>(node.coordinate[d])<<shift : 0;
            const auto upper=d<D ? (static_cast<std::uint64_t>(node.coordinate[d])+1)<<shift : 1;
            if(upper>std::numeric_limits<std::uint32_t>::max()) throw std::length_error("Atomic block coordinate capacity exceeded.");
            block.lower[d]=static_cast<std::uint32_t>(lower); block.upper[d]=static_cast<std::uint32_t>(upper);
        }
        nodeBlock[n]=graph.blocks.size(); graph.blocks.push_back(block);
    }
    const auto blockCount=graph.blocks.size();
    const auto blockOf=[&](GlobalCellIndex cell) {
        auto node=cells.sourceNode[cell];
        while(nodes[node].level>options.blockLevel) node=nodes[node].parent;
        return nodeBlock[node];
    };
    const auto cellWork=[&](GlobalCellIndex c) {
        const auto level=grid.level(c);
        const bool auxiliary=cells.isAuxiliary(c);
        const double boundaryExtra=!auxiliary && cells.cellType[c]==geometry::CellType::Boundary ?
            options.boundaryExtraWeight*std::max(std::uint32_t{1},boundaryMemberships[c]) : 0.0;
        const double base=auxiliary ? options.auxiliaryWeight : options.wetWeight+boundaryExtra;
        const auto costLevel=level-(auxiliary ? 1 : 0);
        const double result=base*detail::frequency(costLevel,options)*detail::levelCost(costLevel,options);
        if(base>0) detail::positiveWeight(result); else detail::nonnegative(result);
        return result;
    };
    unsigned maxLevel=0;
    for(GlobalCellIndex c=0; c<cells.size(); ++c)
    {
        maxLevel=std::max(maxLevel,grid.level(c));
        if(!cells.isAuxiliary(c) && !cells.computes(c)) throw std::logic_error("Dry real cell in partition input.");
        if(cells.isAuxiliary(c)) detail::auxiliarySource(grid,c);
        auto& work=graph.blocks[blockOf(c)].work; work+=cellWork(c); detail::nonnegative(work);
    }
    if(!options.levelWeights.empty() && options.levelWeights.size()<=maxLevel)
        throw std::invalid_argument("Level cost weights do not cover the grid.");
    struct EdgeUnits
    {
        struct Count { unsigned level; std::uint64_t stream=0,interface=0; };
        std::vector<Count> levels;
        Count& at(unsigned level)
        {
            for(auto& count:levels) if(count.level==level) return count;
            levels.push_back({level,0,0}); return levels.back();
        }
    };
    using EdgeMap=std::map<std::pair<std::size_t,std::size_t>,EdgeUnits>;
    const auto add=[](std::uint64_t& to,std::uint64_t value) {
        if(to>std::numeric_limits<std::uint64_t>::max()-value) throw std::overflow_error("Graph edge count overflow.");
        to+=value;
    };
    const int threads=std::min(options.graphThreads,omp_get_max_threads());
    std::vector<EdgeMap> edges(static_cast<std::size_t>(threads));
    std::exception_ptr error;
    // Fixed integer counts make parallel graph aggregation independent of
    // floating-point accumulation order and OpenMP scheduling.
#pragma omp parallel for schedule(static) num_threads(threads)
    for(std::ptrdiff_t id=0; id<static_cast<std::ptrdiff_t>(cells.size()); ++id)
    {
        try
        {
            const auto c=static_cast<GlobalCellIndex>(id);
            if(cells.isAuxiliary(c)) continue;
            const auto level=grid.level(c); const auto a=blockOf(c); const auto p=grid.coordinate(c);
            using Lattice=lattice::Model<Type>;
            for(std::size_t link=1; link<Lattice::nStencils; ++link)
            {
                typename geometry::MultilevelGrid<D,Type>::SignedCoordinate neighbour{};
                neighbour[0]=static_cast<std::int64_t>(p[0])+Lattice::ex[link];
                neighbour[1]=static_cast<std::int64_t>(p[1])+Lattice::ey[link];
                if constexpr(D==3) neighbour[2]=static_cast<std::int64_t>(p[2])+Lattice::ez[link];
                const auto n=grid.coveringNode(level,neighbour);
                if(n==geometry::invalidIndex) continue;
                const auto& node=grid.nodes()[n];
                if(!node.leaf() || node.level!=level || node.cell==geometry::invalidIndex || node.cell<=c) continue;
                const auto b=blockOf(node.cell);
                if(a!=b) add(edges[omp_get_thread_num()][std::minmax(a,b)].at(level).stream,1);
            }
        }
        catch(...)
        {
#pragma omp critical(ntic_partition_graph_error)
            { if(!error) error=std::current_exception(); }
        }
    }
    if(error) std::rethrow_exception(error);
    EdgeMap merged;
    for(auto& thread:edges)
        for(const auto& item:thread)
            for(const auto& count:item.second.levels) add(merged[item.first].at(count.level).stream,count.stream);
    edges.clear();
    std::vector<double> interfaceWork(blockCount,0);
    const auto& interfaces=grid.interfaces();
    for(const auto& block:interfaces.blocks)
    {
        const auto a=blockOf(block.coarse); const auto level=grid.level(block.coarse);
        const double work=options.interfaceWeight*detail::frequency(level,options)*detail::levelCost(level,options);
        detail::nonnegative(work);
        for(auto i=block.fine.begin; i<block.fine.end; ++i)
        {
            const auto b=blockOf(interfaces.fineCells[i]);
            graph.blocks[a].work+=0.5*work; graph.blocks[b].work+=0.5*work;
            interfaceWork[a]+=0.5*work; interfaceWork[b]+=0.5*work;
            if(a!=b) add(merged[std::minmax(a,b)].at(level).interface,1);
        }
    }
    for(auto& item:merged)
    {
        PartitionEdge edge; edge.a=item.first.first; edge.b=item.first.second;
        std::sort(item.second.levels.begin(),item.second.levels.end(),[](const auto& a,const auto& b) { return a.level<b.level; });
        for(const auto& count:item.second.levels)
        {
            const auto l=count.level;
            const double frequency=detail::frequency(l,options);
            const double stream=count.stream*frequency,coupling=count.interface*frequency;
            edge.weight+=stream*options.streamingEdgeWeight+coupling*options.interfaceEdgeWeight;
            edge.bytes+=stream*options.streamingBytes+coupling*options.interfaceBytes;
            if(count.stream) edge.phases|=std::uint64_t{1}<<l;
            if(count.interface) edge.phases|=std::uint64_t{1}<<(l+32);
        }
        graph.edges.push_back(edge);
    }
    auto plan=partitionGraph(graph,options); const auto blockOwners=std::move(plan.owner);
    const auto atomToBlock=std::move(plan.cellBlocks);
    plan.dimension=D; plan.blockLevel=options.blockLevel; plan.blockSpacing=grid.spacing(options.blockLevel);
    for(std::size_t d=0; d<D; ++d) plan.domainMin[d]=grid.domain().min[d];
    plan.owner.resize(cells.size()); plan.offsets.assign(static_cast<std::size_t>(options.ranks)+1,0);
    plan.cellBlocks.resize(cells.size());
    plan.statistics.realCells.assign(options.ranks,0); plan.statistics.auxiliaryCells.assign(options.ranks,0);
    plan.statistics.interfaceWork.assign(options.ranks,0);
    plan.statistics.workPerLevel.assign(options.ranks,std::vector<double>(maxLevel+1,0));
    for(std::size_t b=0; b<blockCount; ++b)
        if(blockOwners[b]>=0) plan.statistics.interfaceWork[blockOwners[b]]+=interfaceWork[b];
    for(GlobalCellIndex c=0; c<cells.size(); ++c)
    {
        const int rank=blockOwners[blockOf(c)]; plan.owner[c]=rank;
        plan.cellBlocks[c]=atomToBlock[blockOf(c)];
        ++plan.offsets[static_cast<std::size_t>(rank)+1];
        plan.statistics.workPerLevel[rank][grid.level(c)]+=cellWork(c);
        if(cells.isAuxiliary(c)) ++plan.statistics.auxiliaryCells[rank]; else ++plan.statistics.realCells[rank];
    }
    std::partial_sum(plan.offsets.begin(),plan.offsets.end(),plan.offsets.begin());
    plan.cells.resize(cells.size()); auto cursor=plan.offsets;
    for(GlobalCellIndex c=0; c<cells.size(); ++c) plan.cells[cursor[plan.owner[c]]++]=c;
    return plan;
}

// Metadata preparation for one rank. Ghosts conservatively include full
// touching-cell adjacency and complete coarse/fine interface blocks. These
// lists are not a lattice-specific streaming or interpolation operator.
template<std::size_t D,lattice::LatticeType Type>
LocalPartition localPartition(const geometry::MultilevelGrid<D,Type>& grid,const PartitionPlan& plan,int rank)
{
    detail::validatePlan(grid,plan);
    if(rank<0 || rank>=plan.ranks()) throw std::invalid_argument("Local rank is out of range.");
    LocalPartition local; local.rank=rank;
    local.globalCells.assign(plan.cells.begin()+plan.offsets[rank],plan.cells.begin()+plan.offsets[rank+1]);
    local.ownedCount=local.globalCells.size();
    struct Exchange { std::vector<GlobalCellIndex> send,receive; };
    std::map<int,Exchange> exchanges;
    std::map<std::pair<std::size_t,std::size_t>,std::vector<std::pair<GlobalCellIndex,GlobalCellIndex>>> internal;
    const auto dependency=[&](GlobalCellIndex owned,GlobalCellIndex remote) {
        if(plan.owner[remote]!=rank)
        { auto& peer=exchanges[plan.owner[remote]]; peer.send.push_back(owned); peer.receive.push_back(remote); }
        else
        {
            const auto a=plan.cellBlocks[owned],b=plan.cellBlocks[remote];
            if(a<b) internal[{a,b}].push_back({owned,remote});
            else if(b<a) internal[{b,a}].push_back({remote,owned});
        }
    };
    for(auto c:local.globalCells)
        if(!grid.cells().isAuxiliary(c))
            for(auto neighbour:grid.adjacentRealCells(c)) dependency(c,neighbour);
    std::vector<std::size_t> localBlocks;
    const auto& connectivity=grid.interfaces();
    for(std::size_t b=0; b<connectivity.blocks.size(); ++b)
    {
        const auto& block=connectivity.blocks[b];
        std::vector<GlobalCellIndex> participants{block.coarse};
        for(auto c=block.auxiliary.begin; c<block.auxiliary.end; ++c) participants.push_back(c);
        for(auto i=block.fine.begin; i<block.fine.end; ++i) participants.push_back(connectivity.fineCells[i]);
        if(std::none_of(participants.begin(),participants.end(),[&](auto c) { return plan.owner[c]==rank; })) continue;
        localBlocks.push_back(b);
        // Each participating rank receives the entire block. Auxiliary cells
        // stay owned by the coarse rank, including physically Dry children.
        for(auto owned:participants) if(plan.owner[owned]==rank)
            for(auto remote:participants) dependency(owned,remote);
    }
    const auto unique=[](auto& values) {
        std::sort(values.begin(),values.end()); values.erase(std::unique(values.begin(),values.end()),values.end());
    };
    std::vector<GlobalCellIndex> ghosts;
    for(auto& item:exchanges)
    {
        unique(item.second.send); unique(item.second.receive);
        ghosts.insert(ghosts.end(),item.second.receive.begin(),item.second.receive.end());
    }
    unique(ghosts);
    local.globalCells.insert(local.globalCells.end(),ghosts.begin(),ghosts.end());
    local.lookupLocal.resize(local.globalCells.size()); std::iota(local.lookupLocal.begin(),local.lookupLocal.end(),0);
    std::sort(local.lookupLocal.begin(),local.lookupLocal.end(),[&](auto a,auto b) { return local.globalCells[a]<local.globalCells[b]; });
    for(auto c:local.lookupLocal) local.lookupGlobal.push_back(local.globalCells[c]);
    for(auto c:local.globalCells) local.cellBlocks.push_back(plan.cellBlocks[c]);
    for(auto& item:internal)
    {
        unique(item.second);
        LocalBlockExchange exchange; exchange.a=item.first.first; exchange.b=item.first.second;
        for(const auto& pair:item.second)
            exchange.dependencies.push_back({local.localIndex(pair.first),local.localIndex(pair.second)});
        local.intraRankExchanges.push_back(std::move(exchange));
    }
    for(std::size_t c=0; c<local.globalCells.size(); ++c)
    {
        const auto global=local.globalCells[c]; const auto& cells=grid.cells();
        const auto level=grid.level(global); const auto role=cells.role[global];
        const int type=cells.isAuxiliary(global) ? -1 : cells.cellType[global]==geometry::CellType::Interior ? 0 : 1;
        const bool owned=c<local.ownedCount;
        if(local.ranges.empty() || local.ranges.back().level!=level || local.ranges.back().role!=role ||
           local.ranges.back().typeOrder!=type || local.ranges.back().owned!=owned)
            local.ranges.push_back({level,role,type,owned,c,c+1});
        else local.ranges.back().end=c+1;
    }
    for(const auto& item:exchanges)
    {
        PeerExchange peer; peer.rank=item.first;
        for(auto c:item.second.send) peer.send.push_back(local.localIndex(c));
        for(auto c:item.second.receive) peer.receive.push_back(local.localIndex(c));
        local.peers.push_back(std::move(peer));
    }
    for(auto b:localBlocks)
    {
        const auto& block=connectivity.blocks[b]; LocalInterfaceBlock result;
        result.coarse=local.localIndex(block.coarse);
        for(auto c=block.auxiliary.begin; c<block.auxiliary.end; ++c) result.auxiliary.push_back(local.localIndex(c));
        for(auto i=block.fine.begin; i<block.fine.end; ++i) result.fine.push_back(local.localIndex(connectivity.fineCells[i]));
        local.interfaces.push_back(std::move(result));
    }
    const auto& boundary=grid.boundaries();
    constexpr auto links=geometry::MultilevelGrid<D,Type>::latticeType;
    constexpr std::size_t qCount=lattice::Model<links>::nStencils-1;
    local.boundaries.offsets.push_back(0);
    for(std::size_t face=0; face<boundary.nBoundaries; ++face)
    {
        for(auto i=boundary.offsets[face]; i<boundary.offsets[face+1]; ++i)
        {
            const auto cell=local.localIndex(boundary.cells[i]);
            if(cell==invalidCellIndex) continue;
            local.boundaries.cells.push_back(cell);
            local.boundaries.q.insert(local.boundaries.q.end(),boundary.q.begin()+i*qCount,boundary.q.begin()+(i+1)*qCount);
        }
        local.boundaries.offsets.push_back(local.boundaries.cells.size());
    }
    return local;
}

} // namespace ntic::lbm::mpi
