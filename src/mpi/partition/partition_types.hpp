#pragma once

#include <cstddef>
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ntic::lbm::mpi
{

using GlobalCellIndex=std::size_t;
using LocalCellIndex=std::size_t;
inline constexpr std::size_t invalidCellIndex=std::numeric_limits<std::size_t>::max();

struct HardwareProfile
{
    std::vector<double> computeRates; // Empty: all ranks process one work unit per time unit.
    std::vector<int> nodeIds;         // Empty: topology-neutral link parameters below.
    double intraNodeLatency=0.0,interNodeLatency=0.0;
    double intraNodeBandwidth=1.0e9,interNodeBandwidth=1.0e9;
    double intraRankLatency=0.0,intraRankBandwidth=1.0e9;
};

struct PartitionOptions
{
    int ranks=1;
    double wetWeight=1.0;           // Ordinary update cost for every real wet cell, including Boundary.
    double boundaryExtraWeight=1.0; // Per stored Face/cell record; real Boundary cells without records pay once.
    double auxiliaryWeight=0.25; // Per auxiliary slot, per source-coarse step.
    double interfaceWeight=1.0;  // Per coarse/fine pair, per coarse step.
    // Empty: unit per-update cost at all levels. Entries multiply cost only;
    // the physical 2^level step frequency is always applied independently.
    std::vector<double> levelWeights;
    unsigned blockLevel=0; // Split topology into atoms down to this level, where legal.
    double streamingEdgeWeight=1.0,interfaceEdgeWeight=1.0;
    double streamingBytes=8.0,interfaceBytes=64.0;
    double imbalanceTolerance=1.03;
    double loadSlack=0.03;       // Relative to the best feasible compute time found across the search.
    int candidates=4,seed=0,graphThreads=8;
    int refinementRounds=8,searchWidth=4,cutCandidates=4;
    HardwareProfile hardware;
};

// Half-open boxes in blockLevel integer coordinates, shared by every mesh level.
// Each working rank owns exactly one box. Dry padding need not cover the
// domain; every populated atom belongs wholly to one box. Geometry retains
// its own CellIndex order, independently of block ownership.
struct PartitionBlock
{
    std::array<std::uint32_t,3> lower{},upper{};
    int owner=0;
    double work=0.0;
};

struct PartitionEdge
{
    std::size_t a=0,b=0;
    double weight=0.0,bytes=0.0;
    // Streaming phases use bits 0..30, coarse/fine phases use bits 32..62.
    std::uint64_t phases=0;
};

// Public graph input requires mutually nonoverlapping positive-work atoms.
// Zero-work background atoms are unowned and cannot be graph edge endpoints.
struct PartitionGraph
{
    unsigned dimension=0;
    std::vector<PartitionBlock> blocks;
    std::vector<PartitionEdge> edges;
};

// Final rank-box IDs; one box per rank implies cross-rank pairs only.
struct BlockExchange
{
    std::size_t a=0,b=0;
    double bytes=0.0;
    std::uint64_t phases=0;
};

struct PartitionStatistics
{
    std::vector<double> rankWeights;
    std::vector<std::size_t> realCells,auxiliaryCells;
    double totalWeight=0.0;
    double maxToMean=0.0;
    std::vector<double> computeTimes,communicationTimes,interfaceWork;
    std::vector<double> intraRankCommunicationTimes,interRankCommunicationTimes;
    std::vector<std::size_t> blocksPerRank;
    std::vector<std::vector<double>> workPerLevel;
    double maxComputeTime=0.0,maxEstimatedTime=0.0,cutWeight=0.0;
    double largestAtomicWork=0.0;
    std::size_t graphVertices=0,graphEdges=0;
    std::size_t evaluatedCandidates=0,jointCutCandidates=0;
    double bestFeasibleComputeTime=0;
    unsigned refinementRounds=0;
    bool searchExhausted=false; // Bounded proposal cycle retained no new beam ownership.
};

// A snapshot of geometry CellIndex ownership. Geometry indices and arrays are
// never reordered. Rebuild this snapshot after geometry topology is finalized
// again; CellIndex is not a persistent identity across refine/derefine.
struct PartitionPlan
{
    std::vector<int> owner;
    std::vector<std::size_t> offsets;
    std::vector<GlobalCellIndex> cells;
    unsigned dimension=0;
    std::array<double,3> domainMin{};
    unsigned blockLevel=0;
    double blockSpacing=0.0;
    std::vector<PartitionBlock> blocks;
    std::vector<std::size_t> cellBlocks; // Per global cell; partitionGraph alone returns per graph vertex.
    std::vector<BlockExchange> blockExchanges;
    PartitionStatistics statistics;

    int ranks() const { return offsets.empty() ? 0 : static_cast<int>(offsets.size()-1); }
};

// Implemented in metis_partition.cpp; no METIS types escape the public API.
PartitionPlan partitionGraph(const PartitionGraph& graph,const PartitionOptions& options);

struct PeerExchange
{
    int rank=0;
    // Lists on either side are ordered by GlobalCellIndex, so corresponding
    // send/receive lists agree without depending on local storage order.
    std::vector<LocalCellIndex> send,receive;
};

struct LocalBoundaryData
{
    std::vector<std::size_t> offsets;
    std::vector<LocalCellIndex> cells;
    std::vector<double> q;
};

struct LocalInterfaceBlock
{
    LocalCellIndex coarse=invalidCellIndex;
    std::vector<LocalCellIndex> auxiliary,fine;
};

struct LocalCellRange
{
    unsigned level=0;
    std::uint8_t role=0;
    int typeOrder=0; // Interior=0, Boundary=1, auxiliary family=-1.
    bool owned=true;
    std::size_t begin=0,end=0;
};

struct LocalBlockExchange
{
    std::size_t a=0,b=0; // Final plan.blocks indices.
    // Conservative paired dependencies, oriented from block a to block b.
    // Both indices address owned cells; no local ghost allocation is needed.
    std::vector<std::pair<LocalCellIndex,LocalCellIndex>> dependencies;
};

// Only the rank-local metadata are assembled here; no distribution functions
// or MPI messages are allocated. Owned cells precede read-only ghost cells.
struct LocalPartition
{
    int rank=0;
    std::size_t ownedCount=0;
    std::vector<GlobalCellIndex> globalCells;
    std::vector<PeerExchange> peers;
    std::vector<std::size_t> cellBlocks; // Parallel to globalCells, including ghosts.
    std::vector<LocalBlockExchange> intraRankExchanges;
    LocalBoundaryData boundaries;
    std::vector<LocalInterfaceBlock> interfaces;
    std::vector<LocalCellRange> ranges;
    // Sorted global keys with matching local indices; O(log N) indirect lookup
    // avoids a global-sized reverse map on every rank.
    std::vector<GlobalCellIndex> lookupGlobal;
    std::vector<LocalCellIndex> lookupLocal;

    LocalCellIndex localIndex(GlobalCellIndex global) const
    {
        const auto it=std::lower_bound(lookupGlobal.begin(),lookupGlobal.end(),global);
        return it==lookupGlobal.end() || *it!=global ? invalidCellIndex : lookupLocal[it-lookupGlobal.begin()];
    }
};

} // namespace ntic::lbm::mpi
