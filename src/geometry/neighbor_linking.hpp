#pragma once

#include "multilevel_grid.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <exception>
#include <limits>
#include <stdexcept>
#include <vector>
#include <omp.h>

namespace ntic::lbm::geometry
{

// A snapshot of storage addresses, not ownership or distribution storage.
// Existing CellPool IDs remain unchanged; Ghost endpoints follow that prefix.
template<std::size_t D,lattice::LatticeType Type>
struct NeighborLinks
{
    using Traits=lattice::Twist<Type>;
    using Coordinate=std::array<std::int64_t,D>;
    struct Ghost
    {
        unsigned level=0;
        Coordinate coordinate{};
        bool operator<(const Ghost& other) const
        { return level<other.level || (level==other.level && coordinate<other.coordinate); }
        bool operator==(const Ghost& other) const
        { return level==other.level && coordinate==other.coordinate; }
    };
    std::size_t originalCount=0;
    std::vector<std::array<CellIndex,Traits::nAnchors>> anchors;
    std::vector<Ghost> ghosts;
    // Global address ranges for Ghosts, by level; original IDs are not reordered.
    std::vector<IndexRange> ghostLevels;
    std::size_t size() const { return originalCount+ghosts.size(); }
    bool isGhost(CellIndex cell) const { return cell>=originalCount && cell<size(); }
    CellIndex anchor(CellIndex source,std::size_t direction) const
    {
        const auto slot=Traits::directionAnchor.at(direction);
        if(source>=originalCount) throw std::out_of_range("Ghost endpoints have no neighbor table.");
        return slot==0 ? source : anchors.at(source)[slot-1];
    }
};

template<std::size_t D,lattice::LatticeType Type>
NeighborLinks<D,Type> linkNeighbors(const MultilevelGrid<D,Type>& grid)
{
    using Result=NeighborLinks<D,Type>;
    using Ghost=typename Result::Ghost;
    using Traits=typename Result::Traits;
    static_assert(D==Traits::dimension,"Lattice and grid dimensions must match.");
    // cells() rejects unfinalized grids, including an invalidated previous snapshot.
    const auto& cells=grid.cells();
    Result result; result.originalCount=cells.size();
    if(cells.size()>static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()))
        throw std::length_error("Neighbor source count exceeds parallel index range.");
    result.anchors.resize(cells.size());
    result.ghostLevels.resize(grid.levels().size());
    if(cells.size()==0) return result;
    std::vector<CellIndex> auxiliaryBase(grid.nodes().size(),invalidIndex);
    for(const auto& block:grid.interfaces().blocks)
        auxiliaryBase[cells.sourceNode[block.coarse]]=block.auxiliary.begin;
    auto endpoint=[&](CellIndex c,std::size_t a) {
        Ghost key; key.level=grid.level(c);
        const auto coordinate=grid.coordinate(c);
        for(std::size_t d=0;d<D;++d)
            key.coordinate[d]=static_cast<std::int64_t>(coordinate[d])+((Traits::anchorMasks[a]>>d)&1U);
        return key;
    };
    auto existing=[&](const Ghost& key) {
        const auto n=grid.coveringNode(key.level,key.coordinate);
        if(n==invalidIndex) return invalidIndex;
        const auto& node=grid.nodes()[n];
        if(node.level==key.level) return node.cell;
        if(node.level+1!=key.level || auxiliaryBase[n]==invalidIndex) return invalidIndex;
        std::size_t child=0;
        for(std::size_t d=0;d<D;++d) child|=static_cast<std::size_t>(key.coordinate[d]&1)<<d;
        return auxiliaryBase[n]+child;
    };
    const int threads=omp_in_parallel() ? 1 : std::max(1,std::min(omp_get_max_threads(),
        static_cast<int>(std::min<std::size_t>(cells.size(),std::numeric_limits<int>::max()))));
    std::vector<std::vector<Ghost>> missing(static_cast<std::size_t>(threads));
    std::exception_ptr error;
#pragma omp parallel for schedule(static) num_threads(threads) if(!omp_in_parallel())
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(cells.size());++i)
    {
        try
        {
            const auto c=static_cast<CellIndex>(i);
            auto& local=missing[static_cast<std::size_t>(omp_get_thread_num())];
            for(std::size_t a=0;a<Traits::nAnchors;++a)
            {
                const auto key=endpoint(c,a);
                result.anchors[c][a]=existing(key);
                if(result.anchors[c][a]==invalidIndex) local.push_back(key);
            }
        }
        catch(...)
        {
#pragma omp critical(ntic_neighbor_link_error)
            { if(!error) error=std::current_exception(); }
        }
    }
    if(error) std::rethrow_exception(error);
    for(auto& local:missing)
    {
        std::sort(local.begin(),local.end());
        local.erase(std::unique(local.begin(),local.end()),local.end());
        result.ghosts.insert(result.ghosts.end(),local.begin(),local.end());
        std::vector<Ghost>().swap(local);
    }
    std::sort(result.ghosts.begin(),result.ghosts.end());
    result.ghosts.erase(std::unique(result.ghosts.begin(),result.ghosts.end()),result.ghosts.end());
    if(result.ghosts.size()>invalidIndex-result.originalCount)
        throw std::length_error("Neighbor address space exceeds CellIndex.");
    std::size_t cursor=0;
    for(unsigned level=0;level<result.ghostLevels.size();++level)
    {
        const auto begin=cursor;
        while(cursor<result.ghosts.size() && result.ghosts[cursor].level==level) ++cursor;
        result.ghostLevels[level]={result.originalCount+begin,result.originalCount+cursor};
    }
#pragma omp parallel for schedule(static) num_threads(threads) if(!omp_in_parallel())
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(cells.size());++i)
    {
        const auto c=static_cast<CellIndex>(i);
        for(std::size_t a=0;a<Traits::nAnchors;++a)
            if(result.anchors[c][a]==invalidIndex)
            {
                const auto key=endpoint(c,a);
                result.anchors[c][a]=result.originalCount+static_cast<std::size_t>(
                    std::lower_bound(result.ghosts.begin(),result.ghosts.end(),key)-result.ghosts.begin());
            }
    }
    return result;
}

} // namespace ntic::lbm::geometry
