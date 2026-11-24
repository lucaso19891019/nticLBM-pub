#pragma once

#include "partition.hpp"
#include "neighbor_linking.hpp"

#include <cstdint>
#include <map>
#include <numeric>
#include <tuple>
#include <omp.h>

namespace ntic::lbm::mpi
{

struct LocalIndexSpan { std::size_t offset=0,length=0; };
enum class LayoutGroup : std::uint8_t
{
    All,Owned,Computes,Auxiliary,Ghost,Halo,Send,Receive,ReceiveDependent,
    Interior,Boundary,FineInterface,CoarseInterface
};
struct LocalLayoutGroup
{
    unsigned level=0;
    LayoutGroup kind=LayoutGroup::All;
    LocalIndexSpan ranges,points; // Views into selectionRanges/selectionPoints.
};
struct LocalLayoutCell
{
    unsigned level=0;
    int owner=-1;
    GlobalCellIndex global=invalidCellIndex;
    geometry::NodeIndex globalSourceNode=geometry::invalidIndex;
    std::size_t localSourceNode=invalidCellIndex;
    std::uint8_t role=geometry::Ordinary;
    geometry::CellType type=geometry::CellType::Dry;
    bool owned=false,ghost=false,halo=false,send=false,receiveDependent=false;
};
struct LocalTreeCell
{
    geometry::NodeIndex globalNode=geometry::invalidIndex;
    LocalCellIndex realCell=invalidCellIndex;
};
struct LocalLayoutOptions { bool tryMorton=true,tryTypeFirst=true; };
struct LocalLayoutScore
{
    std::uint64_t communicationRanges=0,interfaceRanges=0,boundaryRanges=0,anchorDistance=0;
    auto key() const { return std::tie(communicationRanges,interfaceRanges,boundaryRanges,anchorDistance); }
    bool operator<(const LocalLayoutScore& other) const { return key()<other.key(); }
    bool operator==(const LocalLayoutScore& other) const { return key()==other.key(); }
};

template<std::size_t D>
struct LayoutFrame
{
    std::array<double,D> origin{};
    double coarseSpacing=1;
    unsigned levelCount=0;
};

template<std::size_t D>
struct RankTreeNode
{
    geometry::NodeIndex globalNode=geometry::invalidIndex,parentGlobal=geometry::invalidIndex,firstGlobalChild=geometry::invalidIndex;
    GlobalCellIndex globalRealCell=invalidCellIndex;
    unsigned level=0;
    std::array<std::uint32_t,D> coordinate{};
    bool active=false,leaf=false;
    std::array<geometry::NodeIndex,std::size_t{1}<<D> globalChildren{};
};

template<std::size_t D>
struct LocalLayoutNode : RankTreeNode<D>
{
    std::size_t parent=invalidCellIndex;
    LocalCellIndex realCell=invalidCellIndex;
    std::array<std::size_t,std::size_t{1}<<D> children{};
};

struct RankLayoutPeer
{
    int rank=0;
    std::vector<GlobalCellIndex> send,receive;
};
struct RankLayoutInterface
{
    GlobalCellIndex coarse=invalidCellIndex;
    std::vector<GlobalCellIndex> auxiliary,fine;
};
struct RankLayoutBlockExchange
{
    std::size_t a=0,b=0;
    std::vector<std::pair<GlobalCellIndex,GlobalCellIndex>> dependencies;
};

template<std::size_t D,lattice::LatticeType Type>
struct RankLayoutInput
{
    using Traits=lattice::Twist<Type>;
    int rank=0,rankCount=0;
    LayoutFrame<D> frame;
    // Sorted global storage IDs, independent of every rank-zero array.
    std::vector<LocalLayoutCell> cells;
    std::vector<std::array<std::int64_t,D>> coordinates;
    std::vector<std::uint8_t> auxiliaryChild;
    std::vector<GlobalCellIndex> sourceIDs;
    std::vector<std::array<GlobalCellIndex,Traits::nAnchors>> anchors;
    std::vector<RankLayoutPeer> peers;
    LocalBoundaryData boundaries; // cells are GLOBAL storage IDs in this input.
    std::vector<RankLayoutInterface> interfaces;
    std::vector<RankLayoutBlockExchange> intraRankExchanges;
    std::vector<RankTreeNode<D>> nodes; // Source nodes and their ancestor closure.
};

template<std::size_t D,lattice::LatticeType Type>
struct LocalLayout
{
    using Traits=lattice::Twist<Type>;
    int rank=0;
    LayoutFrame<D> frame;
    std::vector<LocalLayoutCell> cells;
    std::vector<std::array<std::int64_t,D>> coordinates;
    std::vector<LocalLayoutNode<D>> nodes;
    std::vector<GlobalCellIndex> lookupGlobal;
    std::vector<LocalCellIndex> lookupLocal;
    // Only owned original real/auxiliary sources have rows; endpoints have none.
    std::vector<LocalCellIndex> sourceIDs;
    std::vector<std::array<LocalCellIndex,Traits::nAnchors>> anchors;
    std::vector<PeerExchange> peers;
    LocalBoundaryData boundaries;
    std::vector<LocalIndexSpan> boundaryFaces;
    std::vector<LocalInterfaceBlock> interfaces;
    std::vector<LocalBlockExchange> intraRankExchanges;
    std::vector<LocalIndexSpan> levelRanges;
    std::vector<LocalLayoutGroup> groups;
    std::vector<LocalIndexSpan> selectionRanges;
    std::vector<LocalCellIndex> selectionPoints;
    std::vector<LocalTreeCell> treeCells;
    LocalLayoutScore score;
    unsigned candidates=0,selectedCandidate=0;
    std::array<double,D> center(LocalCellIndex cell) const
    {
        const auto& coordinate=coordinates.at(cell);
        const double h=std::ldexp(frame.coarseSpacing,-static_cast<int>(cells.at(cell).level));
        std::array<double,D> result{};
        for(std::size_t d=0;d<D;++d) result[d]=std::fma(static_cast<double>(coordinate[d])+0.5,h,frame.origin[d]);
        return result;
    }
    LocalCellIndex anchor(LocalCellIndex source,std::size_t direction) const
    {
        const auto row=std::lower_bound(sourceIDs.begin(),sourceIDs.end(),source);
        if(row==sourceIDs.end() || *row!=source) throw std::out_of_range("Only owned original sources have neighbor rows.");
        const auto slot=Traits::directionAnchor.at(direction);
        return slot==0 ? source : anchors[static_cast<std::size_t>(row-sourceIDs.begin())][slot-1];
    }
    LocalCellIndex localIndex(GlobalCellIndex global) const
    {
        const auto it=std::lower_bound(lookupGlobal.begin(),lookupGlobal.end(),global);
        return it==lookupGlobal.end() || *it!=global ? invalidCellIndex : lookupLocal[it-lookupGlobal.begin()];
    }
};

namespace detail
{
// Full 32 bits per coordinate; words are most-significant first.
template<std::size_t D>
constexpr std::array<std::uint32_t,3> layoutMortonKey(const std::array<std::uint32_t,D>& coordinate)
{
    std::array<std::uint32_t,3> result{};
    for(unsigned bit=0;bit<32;++bit) for(std::size_t d=0;d<D;++d)
    {
        const auto position=D*bit+d;
        result[2-position/32]|=((coordinate[d]>>bit)&1U)<<(position%32);
    }
    return result;
}
// Stable LSD radix passes preserve atom order for equal 96-bit keys.
// Each thread owns a contiguous input chunk and its histogram/scatter cursors.
template<typename Digit>
void layoutRadixOrder(std::vector<std::size_t>& order,unsigned passes,Digit digitAt)
{
    if(order.size()<2) return;
    const int threads=omp_in_parallel() || order.size()<4096 ? 1 : std::max(1,std::min(omp_get_max_threads(),
        static_cast<int>(std::min<std::size_t>(order.size(),std::numeric_limits<int>::max()))));
    std::vector<std::size_t> temporary(order.size());
    std::vector<std::array<std::size_t,256>> histogram(static_cast<std::size_t>(threads)),offsets(histogram.size());
    for(unsigned pass=0;pass<passes;++pass)
    {
        auto digit=[&](std::size_t id) { return digitAt(id,pass); };
        if(threads==1)
        {
            auto& counts=histogram[0]; counts.fill(0);
            for(auto id:order) ++counts[digit(id)];
            std::size_t cursor=0;
            for(unsigned b=0;b<256;++b) { offsets[0][b]=cursor; cursor+=counts[b]; }
            for(auto id:order) temporary[offsets[0][digit(id)]++]=id;
        }
        else
        {
#pragma omp parallel num_threads(threads)
            {
                const auto thread=static_cast<std::size_t>(omp_get_thread_num());
                const auto team=static_cast<std::size_t>(omp_get_num_threads());
                const auto quotient=order.size()/team,remainder=order.size()%team;
                const auto begin=thread*quotient+std::min(thread,remainder);
                const auto end=begin+quotient+(thread<remainder ? 1 : 0);
                histogram[thread].fill(0);
                for(auto i=begin;i<end;++i) ++histogram[thread][digit(order[i])];
#pragma omp barrier
#pragma omp single
                {
                    std::size_t cursor=0;
                    for(unsigned b=0;b<256;++b) for(std::size_t t=0;t<team;++t)
                    { offsets[t][b]=cursor; cursor+=histogram[t][b]; }
                }
                for(auto i=begin;i<end;++i) temporary[offsets[thread][digit(order[i])]++]=order[i];
            }
        }
        order.swap(temporary);
    }
}

template<typename Atoms>
void layoutMortonOrder(std::vector<std::size_t>& order,const Atoms& atoms)
{
    std::uint32_t high=0,middle=0,low=0;
#pragma omp parallel for schedule(static) reduction(|:high,middle,low) if(!omp_in_parallel() && order.size()>=4096)
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(order.size());++i)
    {
        const auto& key=atoms[order[static_cast<std::size_t>(i)]].morton;
        high|=key[0]; middle|=key[1]; low|=key[2];
    }
    const std::array<std::uint32_t,3> bits{high,middle,low};
    unsigned passes=12;
    while(passes && ((bits[2-(passes-1)/4]>>(((passes-1)%4)*8))&255U)==0) --passes;
    layoutRadixOrder(order,passes,[&](std::size_t id,unsigned pass) {
        return (atoms[id].morton[2-pass/4]>>((pass%4)*8))&255U;
    });
}

inline std::size_t layoutRangeCount(std::vector<LocalCellIndex> values,const std::vector<std::size_t>& levelStarts)
{
    if(values.empty()) return 0;
    std::sort(values.begin(),values.end()); values.erase(std::unique(values.begin(),values.end()),values.end());
    std::size_t count=1;
    for(std::size_t i=1;i<values.size();++i)
        count+=values[i]!=values[i-1]+1 || std::binary_search(levelStarts.begin(),levelStarts.end(),values[i]);
    return count;
}
} // namespace detail

// Borrowed snapshots must outlive this builder and stay unchanged during use.
// Global source anchors are scanned during preparation, not once per rank.
template<std::size_t D,lattice::LatticeType Type>
class LocalLayoutBuilder
{
    using Grid=geometry::MultilevelGrid<D,Type>;
    using Links=geometry::NeighborLinks<D,Type>;
    using Traits=lattice::Twist<Type>;
    using Result=LocalLayout<D,Type>;
    struct Exchange { std::vector<GlobalCellIndex> send,receive; };
    const Grid& grid_;
    const Links& links_;
    const PartitionPlan& plan_;
    std::vector<int> ghostOwners_;
    std::vector<std::map<int,Exchange>> extra_;
    std::vector<std::uint8_t> sourceReceive_,sourceSend_;
    std::vector<GlobalCellIndex> auxiliaryBase_;
    static void unique(std::vector<GlobalCellIndex>& values)
    { std::sort(values.begin(),values.end()); values.erase(std::unique(values.begin(),values.end()),values.end()); }
    int owner(GlobalCellIndex id) const
    { return id<links_.originalCount ? plan_.owner[id] : ghostOwners_[id-links_.originalCount]; }
    unsigned level(GlobalCellIndex id) const
    { return id<links_.originalCount ? grid_.level(id) : links_.ghosts[id-links_.originalCount].level; }
    std::array<std::int64_t,D> coordinate(GlobalCellIndex id) const
    {
        if(id>=links_.originalCount) return links_.ghosts[id-links_.originalCount].coordinate;
        const auto p=grid_.coordinate(id); std::array<std::int64_t,D> result{};
        std::copy(p.begin(),p.end(),result.begin()); return result;
    }
    void request(int rank,GlobalCellIndex id)
    {
        const auto remote=owner(id); if(remote==rank) return;
        extra_[rank][remote].receive.push_back(id); extra_[remote][rank].send.push_back(id);
    }
public:
    LocalLayoutBuilder(const Grid& grid,const Links& links,const PartitionPlan& plan)
        :grid_(grid),links_(links),plan_(plan)
    {
        detail::validatePlan(grid,plan);
        if(links.originalCount!=grid.cells().size() || links.anchors.size()!=links.originalCount)
            throw std::invalid_argument("Neighbor snapshot does not match the finalized grid.");
        if(links.originalCount>static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()))
            throw std::length_error("Global layout source count exceeds parallel index range.");
        ghostOwners_.assign(links.ghosts.size(),plan.ranks());
        const int threads=omp_in_parallel() ? 1 : std::max(1,std::min(omp_get_max_threads(),
            static_cast<int>(std::min<std::size_t>(links.originalCount,std::numeric_limits<int>::max()))));
        using Request=std::pair<GlobalCellIndex,int>;
        std::vector<std::vector<Request>> requests(static_cast<std::size_t>(threads));
        std::exception_ptr error;
#pragma omp parallel for schedule(static) num_threads(threads) if(!omp_in_parallel())
        for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(links.originalCount);++i)
        {
            try
            {
                const auto c=static_cast<GlobalCellIndex>(i);
                auto& local=requests[static_cast<std::size_t>(omp_get_thread_num())];
                for(auto id:links.anchors[c])
                {
                    if(id>=links.size()) throw std::invalid_argument("Invalid global storage anchor.");
                    if(id>=links.originalCount) local.push_back({id,plan.owner[c]});
                }
            }
            catch(...)
            {
#pragma omp critical(ntic_layout_prepare_error)
                { if(!error) error=std::current_exception(); }
            }
        }
        if(error) std::rethrow_exception(error);
        auto compact=[](auto& local) {
            std::sort(local.begin(),local.end()); local.erase(std::unique(local.begin(),local.end()),local.end());
        };
        for(auto& local:requests)
        {
            compact(local);
            for(const auto& entry:local)
                ghostOwners_[entry.first-links.originalCount]=std::min(ghostOwners_[entry.first-links.originalCount],entry.second);
            local.clear();
        }
        for(auto rank:ghostOwners_) if(rank==plan.ranks())
            throw std::invalid_argument("Unreferenced Ghost in neighbor snapshot.");
        extra_.resize(static_cast<std::size_t>(plan.ranks()));
        sourceReceive_.assign(links.originalCount,0); sourceSend_.assign(links.originalCount,0);
        auxiliaryBase_.assign(grid.nodes().size(),invalidCellIndex);
        for(const auto& block:grid.interfaces().blocks)
            auxiliaryBase_[grid.cells().sourceNode[block.coarse]]=block.auxiliary.begin;
#pragma omp parallel for schedule(static) num_threads(threads) if(!omp_in_parallel())
        for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(links.originalCount);++i)
        {
            try
            {
                const auto c=static_cast<GlobalCellIndex>(i);
                auto& local=requests[static_cast<std::size_t>(omp_get_thread_num())];
                for(auto id:links.anchors[c]) if(owner(id)!=plan.owner[c])
                {
                    sourceReceive_[c]=1;
                    if(id<links.originalCount && grid.cells().isAuxiliary(id))
                    {
                        const auto begin=auxiliaryBase_[grid.cells().sourceNode[id]];
                        if(begin==invalidCellIndex) throw std::logic_error("Auxiliary storage lacks its interface family.");
                        for(std::size_t child=0;child<(std::size_t{1}<<D);++child) local.push_back({begin+child,plan.owner[c]});
                    }
                    else local.push_back({id,plan.owner[c]});
                }
            }
            catch(...)
            {
#pragma omp critical(ntic_layout_prepare_error)
                { if(!error) error=std::current_exception(); }
            }
        }
        if(error) std::rethrow_exception(error);
        for(auto& local:requests) { compact(local); for(const auto& entry:local) request(entry.second,entry.first); }
        std::vector<std::uint8_t> requested(links.size(),0);
        for(auto& rank:extra_) for(auto& peer:rank)
        { unique(peer.second.send); unique(peer.second.receive); for(auto id:peer.second.send) requested[id]=1; }
#pragma omp parallel for schedule(static) num_threads(threads) if(!omp_in_parallel())
        for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(links.originalCount);++i)
        {
            const auto c=static_cast<GlobalCellIndex>(i);
            sourceSend_[c]=requested[c];
            for(auto id:links.anchors[c]) if(owner(id)==plan.owner[c] && requested[id]) sourceSend_[c]=1;
        }
        std::vector<GlobalCellIndex>().swap(auxiliaryBase_);
    }
    RankLayoutInput<D,Type> prepare(int rank) const;
    Result build(int rank,const LocalLayoutOptions& options={}) const;
};

template<std::size_t D,lattice::LatticeType Type>
auto LocalLayoutBuilder<D,Type>::prepare(int rank) const -> RankLayoutInput<D,Type>
{
    const auto base=localPartition(grid_,plan_,rank);
    std::map<int,Exchange> exchanges;
    std::vector<GlobalCellIndex> ids=base.globalCells;
    for(const auto& peer:base.peers)
    {
        auto& e=exchanges[peer.rank];
        for(auto c:peer.send) e.send.push_back(base.globalCells[c]);
        for(auto c:peer.receive) e.receive.push_back(base.globalCells[c]);
    }
    for(const auto& peer:extra_[rank])
    {
        auto& e=exchanges[peer.first];
        e.send.insert(e.send.end(),peer.second.send.begin(),peer.second.send.end());
        e.receive.insert(e.receive.end(),peer.second.receive.begin(),peer.second.receive.end());
    }
    // Only this rank's sources are traversed here. All ranks together visit N*K.
    for(auto i=plan_.offsets[rank];i<plan_.offsets[rank+1];++i)
        for(auto id:links_.anchors[plan_.cells[i]]) ids.push_back(id);
    std::vector<GlobalCellIndex> sent,dependent;
    for(auto& peer:exchanges)
    {
        unique(peer.second.send); unique(peer.second.receive);
        ids.insert(ids.end(),peer.second.send.begin(),peer.second.send.end());
        ids.insert(ids.end(),peer.second.receive.begin(),peer.second.receive.end());
        sent.insert(sent.end(),peer.second.send.begin(),peer.second.send.end());
    }
    // Legacy send cells conservatively identify owned adjacency/interface work
    // with a remote dependency; peer.receive itself contains the remote Halo.
    for(const auto& peer:base.peers)
        for(auto c:peer.send) dependent.push_back(base.globalCells[c]);
    unique(sent); unique(dependent); unique(ids);
    if(ids.size()>static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()))
        throw std::length_error("Local layout exceeds parallel index range.");
    std::vector<LocalLayoutCell> records(ids.size());
    const auto& pool=grid_.cells();
#pragma omp parallel for schedule(static) if(!omp_in_parallel())
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(ids.size());++i)
    {
        const auto index=static_cast<std::size_t>(i),id=ids[index]; auto& cell=records[index];
        cell.global=id; cell.level=level(id); cell.owner=owner(id); cell.owned=cell.owner==rank;
        cell.ghost=id>=links_.originalCount; cell.halo=!cell.owned;
        cell.send=std::binary_search(sent.begin(),sent.end(),id);
        if(!cell.ghost)
        {
            cell.globalSourceNode=pool.sourceNode[id]; cell.role=pool.role[id]; cell.type=pool.cellType[id];
            if(cell.owned)
            {
                cell.receiveDependent=sourceReceive_[id] || std::binary_search(dependent.begin(),dependent.end(),id);
                cell.send=cell.send || sourceSend_[id];
                // A source can access an owned endpoint requested by another rank
                // even when the source itself is absent from peer.send.
                for(auto endpoint:links_.anchors[id])
                    if(owner(endpoint)==rank && std::binary_search(sent.begin(),sent.end(),endpoint)) cell.send=true;
            }
        }
    }
    RankLayoutInput<D,Type> input; input.rank=rank; input.rankCount=plan_.ranks();
    input.frame.origin=grid_.domain().min; input.frame.coarseSpacing=grid_.spacing(0);
    input.frame.levelCount=static_cast<unsigned>(grid_.levels().size());
    input.cells=std::move(records); input.coordinates.resize(ids.size()); input.auxiliaryChild.resize(ids.size(),255);
#pragma omp parallel for schedule(static) if(!omp_in_parallel())
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(ids.size());++i)
    {
        const auto c=static_cast<std::size_t>(i); input.coordinates[c]=coordinate(ids[c]);
        if(ids[c]<links_.originalCount) input.auxiliaryChild[c]=pool.auxiliaryChild[ids[c]];
    }
    for(const auto& cell:input.cells) if(cell.owned && !cell.ghost)
    { input.sourceIDs.push_back(cell.global); input.anchors.push_back(links_.anchors[cell.global]); }
    for(auto& peer:exchanges) input.peers.push_back({peer.first,std::move(peer.second.send),std::move(peer.second.receive)});
    for(const auto& block:base.interfaces)
    {
        RankLayoutInterface b; b.coarse=base.globalCells[block.coarse];
        for(auto c:block.auxiliary) b.auxiliary.push_back(base.globalCells[c]);
        for(auto c:block.fine) b.fine.push_back(base.globalCells[c]);
        input.interfaces.push_back(std::move(b));
    }
    for(const auto& exchange:base.intraRankExchanges)
    {
        RankLayoutBlockExchange e; e.a=exchange.a; e.b=exchange.b;
        for(const auto& pair:exchange.dependencies) e.dependencies.push_back({base.globalCells[pair.first],base.globalCells[pair.second]});
        input.intraRankExchanges.push_back(std::move(e));
    }
    const auto& boundary=grid_.boundaries(); constexpr auto qCount=lattice::Model<Type>::nStencils-1;
    input.boundaries.offsets.push_back(0);
    for(std::size_t face=0;face<boundary.nBoundaries;++face)
    {
        for(auto i=boundary.offsets[face];i<boundary.offsets[face+1];++i)
            if(std::binary_search(ids.begin(),ids.end(),boundary.cells[i]))
            {
                input.boundaries.cells.push_back(boundary.cells[i]);
                input.boundaries.q.insert(input.boundaries.q.end(),boundary.q.begin()+i*qCount,boundary.q.begin()+(i+1)*qCount);
            }
        input.boundaries.offsets.push_back(input.boundaries.cells.size());
    }
    std::vector<geometry::NodeIndex> nodeIDs;
    for(const auto& cell:input.cells) if(!cell.ghost) nodeIDs.push_back(cell.globalSourceNode);
    std::sort(nodeIDs.begin(),nodeIDs.end()); nodeIDs.erase(std::unique(nodeIDs.begin(),nodeIDs.end()),nodeIDs.end());
    const auto sources=nodeIDs;
    for(auto n:sources) for(auto parent=grid_.nodes()[n].parent;parent!=geometry::invalidIndex;parent=grid_.nodes()[parent].parent)
        nodeIDs.push_back(parent);
    std::sort(nodeIDs.begin(),nodeIDs.end()); nodeIDs.erase(std::unique(nodeIDs.begin(),nodeIDs.end()),nodeIDs.end());
    for(auto n:nodeIDs)
    {
        const auto& source=grid_.nodes()[n]; RankTreeNode<D> node;
        node.globalNode=n; node.parentGlobal=source.parent; node.firstGlobalChild=source.firstChild;
        node.globalRealCell=source.cell; node.level=source.level; node.coordinate=source.coordinate;
        node.active=source.active; node.leaf=source.leaf(); node.globalChildren.fill(geometry::invalidIndex);
        if(source.firstChild!=geometry::invalidIndex)
            for(std::size_t child=0;child<(std::size_t{1}<<D);++child) node.globalChildren[child]=source.firstChild+child;
        input.nodes.push_back(node);
    }
    return input;
}

template<std::size_t D,lattice::LatticeType Type>
void validateRankLayoutInput(const RankLayoutInput<D,Type>& input)
{
    const auto fail=[] { throw std::invalid_argument("Invalid rank-local layout input."); };
    if(input.rank<0 || input.rank>=input.rankCount || input.frame.levelCount==0 || input.frame.levelCount>31 ||
       !(input.frame.coarseSpacing>0) || !std::isfinite(input.frame.coarseSpacing)) fail();
    for(auto x:input.frame.origin) if(!std::isfinite(x)) fail();
    const auto n=input.cells.size();
    if(input.coordinates.size()!=n || input.auxiliaryChild.size()!=n || input.sourceIDs.size()!=input.anchors.size() ||
       n>static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max())) fail();
    std::vector<GlobalCellIndex> ids; ids.reserve(n);
    for(const auto& cell:input.cells)
    {
        if(cell.global==invalidCellIndex || (!ids.empty() && ids.back()>=cell.global) || cell.level>=input.frame.levelCount ||
           cell.owner<0 || cell.owner>=input.rankCount || cell.owned!=(cell.owner==input.rank) || cell.halo==cell.owned ||
           cell.role>geometry::Auxiliary || static_cast<unsigned>(cell.type)>static_cast<unsigned>(geometry::CellType::Interior) ||
           cell.localSourceNode!=invalidCellIndex) fail();
        ids.push_back(cell.global);
    }
    auto position=[&](GlobalCellIndex id) {
        const auto it=std::lower_bound(ids.begin(),ids.end(),id);
        if(it==ids.end() || *it!=id) fail();
        return static_cast<std::size_t>(it-ids.begin());
    };
    std::vector<GlobalCellIndex> expectedSources;
    for(std::size_t i=0;i<n;++i)
    {
        const auto& cell=input.cells[i];
        if(cell.ghost)
        {
            if(cell.globalSourceNode!=geometry::invalidIndex || input.auxiliaryChild[i]!=255 ||
               cell.role!=geometry::Ordinary || cell.type!=geometry::CellType::Dry || cell.receiveDependent) fail();
        }
        else
        {
            if(cell.globalSourceNode==geometry::invalidIndex) fail();
            if(cell.owned) expectedSources.push_back(cell.global);
            if(cell.role==geometry::Auxiliary && input.auxiliaryChild[i]>=(std::size_t{1}<<D)) fail();
        }
    }
    if(input.sourceIDs!=expectedSources) fail();
    for(std::size_t row=0;row<input.sourceIDs.size();++row)
    {
        const auto source=position(input.sourceIDs[row]);
        for(std::size_t a=0;a<lattice::Twist<Type>::nAnchors;++a)
        {
            const auto endpoint=position(input.anchors[row][a]);
            if(input.cells[endpoint].level!=input.cells[source].level) fail();
            for(std::size_t d=0;d<D;++d)
            {
                const auto x=input.coordinates[source][d];
                const auto delta=(lattice::Twist<Type>::anchorMasks[a]>>d)&1U;
                if(x>std::numeric_limits<std::int64_t>::max()-delta || input.coordinates[endpoint][d]!=x+delta) fail();
            }
        }
    }
    int previous=-1;
    for(const auto& peer:input.peers)
    {
        if(peer.rank<=previous || peer.rank<0 || peer.rank>=input.rankCount || peer.rank==input.rank) fail();
        previous=peer.rank;
        auto check=[&](const auto& values,bool send) {
            GlobalCellIndex before=invalidCellIndex;
            for(auto id:values)
            {
                if(before!=invalidCellIndex && before>=id) fail();
                before=id;
                const auto& cell=input.cells[position(id)];
                if(cell.owner!=(send ? input.rank : peer.rank)) fail();
            }
        };
        check(peer.send,true); check(peer.receive,false);
    }
    constexpr auto qCount=lattice::Model<Type>::nStencils-1;
    const auto& boundary=input.boundaries;
    if(boundary.offsets.empty() || boundary.offsets.front()!=0 || boundary.offsets.back()!=boundary.cells.size() ||
       boundary.cells.size()>std::numeric_limits<std::size_t>::max()/qCount || boundary.q.size()!=boundary.cells.size()*qCount) fail();
    for(std::size_t i=1;i<boundary.offsets.size();++i) if(boundary.offsets[i]<boundary.offsets[i-1]) fail();
    for(auto id:boundary.cells) if(input.cells[position(id)].ghost) fail();
    for(const auto& block:input.interfaces)
    {
        (void)position(block.coarse);
        if(block.auxiliary.size()!=(std::size_t{1}<<D)) fail();
        for(auto id:block.auxiliary) if(input.cells[position(id)].role!=geometry::Auxiliary) fail();
        for(auto id:block.fine) (void)position(id);
    }
    for(const auto& exchange:input.intraRankExchanges)
        for(const auto& pair:exchange.dependencies) { (void)position(pair.first); (void)position(pair.second); }
    std::vector<geometry::NodeIndex> nodes;
    for(const auto& node:input.nodes)
    {
        if(node.globalNode==geometry::invalidIndex || (!nodes.empty() && nodes.back()>=node.globalNode) || node.level>=input.frame.levelCount) fail();
        nodes.push_back(node.globalNode);
    }
    auto nodePosition=[&](geometry::NodeIndex id) {
        const auto it=std::lower_bound(nodes.begin(),nodes.end(),id);
        if(it==nodes.end() || *it!=id) fail();
        return static_cast<std::size_t>(it-nodes.begin());
    };
    for(const auto& node:input.nodes)
    {
        if(!node.active || node.leaf!=(node.firstGlobalChild==geometry::invalidIndex)) fail();
        if(!node.leaf && node.firstGlobalChild>geometry::invalidIndex-((std::size_t{1}<<D)-1)) fail();
        for(std::size_t child=0;child<(std::size_t{1}<<D);++child)
            if(node.globalChildren[child]!=(node.leaf ? geometry::invalidIndex : node.firstGlobalChild+child)) fail();
        if(node.parentGlobal!=geometry::invalidIndex)
        {
            const auto& parent=input.nodes[nodePosition(node.parentGlobal)];
            if(parent.level+1!=node.level || parent.leaf) fail();
            std::size_t child=0;
            for(std::size_t d=0;d<D;++d)
            {
                if(node.coordinate[d]/2!=parent.coordinate[d]) fail();
                child|=static_cast<std::size_t>(node.coordinate[d]&1U)<<d;
            }
            if(parent.globalChildren[child]!=node.globalNode) fail();
        }
        else if(node.level!=0) fail();
        for(auto child:node.globalChildren) if(child!=geometry::invalidIndex && std::binary_search(nodes.begin(),nodes.end(),child))
            if(input.nodes[nodePosition(child)].parentGlobal!=node.globalNode) fail();
    }
    for(std::size_t i=0;i<input.cells.size();++i) if(!input.cells[i].ghost)
    {
        const auto& cell=input.cells[i];const auto& source=input.nodes[nodePosition(cell.globalSourceNode)];
        if(!source.leaf) fail();
        const bool auxiliary=cell.role==geometry::Auxiliary;
        if(cell.level!=source.level+(auxiliary ? 1U : 0U)) fail();
        if(!auxiliary && source.globalRealCell!=cell.global) fail();
        for(std::size_t d=0;d<D;++d)
        {
            const auto expected=auxiliary ? 2*static_cast<std::int64_t>(source.coordinate[d])+((input.auxiliaryChild[i]>>d)&1U) : source.coordinate[d];
            if(input.coordinates[i][d]!=expected) fail();
        }
    }
}

template<std::size_t D,lattice::LatticeType Type>
LocalLayout<D,Type> buildLocalLayout(const RankLayoutInput<D,Type>& input,const LocalLayoutOptions& options={})
{
    validateRankLayoutInput(input);
    using Traits=lattice::Twist<Type>; using Result=LocalLayout<D,Type>;
    const auto rank=input.rank; const auto& records=input.cells;
    std::vector<GlobalCellIndex> ids; for(const auto& cell:records) ids.push_back(cell.global);
    // Resolve sparse global addresses once; candidate scores only permute indices.
    std::vector<std::size_t> sourcePositions;
    for(std::size_t i=0;i<records.size();++i) if(records[i].owned && !records[i].ghost) sourcePositions.push_back(i);
    std::vector<std::array<std::size_t,Traits::nAnchors>> endpointPositions(sourcePositions.size());
#pragma omp parallel for schedule(static) if(!omp_in_parallel())
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(sourcePositions.size());++i)
    {
        const auto row=static_cast<std::size_t>(i),c=row;
        for(std::size_t a=0;a<Traits::nAnchors;++a)
            endpointPositions[row][a]=static_cast<std::size_t>(std::lower_bound(ids.begin(),ids.end(),input.anchors[c][a])-ids.begin());
    }
    const auto& globalBoundary=input.boundaries;
    std::vector<std::vector<std::pair<std::size_t,std::size_t>>> boundaryRows((globalBoundary.offsets.size()-1));
    for(std::size_t face=0;face<(globalBoundary.offsets.size()-1);++face)
        for(auto entry=globalBoundary.offsets[face];entry<globalBoundary.offsets[face+1];++entry)
        {
            const auto it=std::lower_bound(ids.begin(),ids.end(),globalBoundary.cells[entry]);
            if(it!=ids.end() && *it==globalBoundary.cells[entry])
                boundaryRows[face].push_back({static_cast<std::size_t>(it-ids.begin()),entry});
        }
    struct Atom
    {
        std::size_t first=0,count=1;
        unsigned level=0,storage=0,communication=2,role=0,type=0;
        std::array<std::uint32_t,3> morton{};
    };
    std::vector<Atom> atoms; atoms.reserve(ids.size());
    bool mortonAvailable=true;
    for(std::size_t i=0;i<records.size();)
    {
        Atom atom; atom.first=i; const auto& first=records[i]; atom.level=first.level;
        atom.storage=first.halo ? 3U : first.ghost ? 2U : first.role==geometry::Auxiliary ? 1U : 0U;
        if(!first.ghost && first.role==geometry::Auxiliary)
        {
            atom.count=std::size_t{1}<<D;
            if(input.auxiliaryChild[i]!=0 || atom.count>records.size()-i)
                throw std::logic_error("Local auxiliary family is incomplete.");
            for(std::size_t child=0;child<atom.count;++child)
                if(records[i+child].global!=first.global+child || records[i+child].globalSourceNode!=first.globalSourceNode ||
                   input.auxiliaryChild[i+child]!=child)
                    throw std::logic_error("Local auxiliary family is incomplete or unordered.");
        }
        for(std::size_t child=0;child<atom.count;++child)
        {
            const auto& cell=records[i+child];
            atom.communication=std::min(atom.communication,cell.halo || cell.receiveDependent ? 0U : cell.send ? 1U : 2U);
            atom.role|=cell.role;
            atom.type=std::max(atom.type,cell.type==geometry::CellType::Boundary ? 1U : cell.type==geometry::CellType::Dry ? 2U : 0U);
        }
        std::array<std::uint32_t,D> p{}; const auto c=input.coordinates[i];
        for(std::size_t d=0;d<D;++d)
        {
            if(c[d]<0 || static_cast<std::uint64_t>(c[d])>std::numeric_limits<std::uint32_t>::max()) mortonAvailable=false;
            else p[d]=static_cast<std::uint32_t>(c[d]);
        }
        atom.morton=detail::layoutMortonKey(p); atoms.push_back(atom); i+=atom.count;
    }
    auto bucket=[](const Atom& a,bool typeFirst) {
        const auto first=typeFirst ? a.type : a.role,second=typeFirst ? a.role : a.type;
        return ((((static_cast<std::uint64_t>(a.level)*4+a.storage)*3+a.communication)*8+first)*8+second);
    };
    std::vector<std::size_t> levelStarts,levelCounts(input.frame.levelCount,0);
    for(const auto& record:records) ++levelCounts.at(record.level);
    std::size_t levelOffset=0;
    for(auto count:levelCounts) { if(count) levelStarts.push_back(levelOffset); levelOffset+=count; }
    std::vector<std::size_t> bestOrder;
    LocalLayoutScore bestScore;
    unsigned candidates=0,selected=0;
    bool hasBest=false;
    auto evaluate=[&](const std::vector<std::size_t>& order) {
        std::vector<LocalCellIndex> inverse(ids.size());
        std::vector<std::size_t> offsets(order.size()+1,0);
        for(std::size_t i=0;i<order.size();++i) offsets[i+1]=offsets[i]+atoms[order[i]].count;
#pragma omp parallel for schedule(static) if(!omp_in_parallel())
        for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(order.size());++i)
        {
            const auto index=static_cast<std::size_t>(i),a=order[index];
            for(std::size_t child=0;child<atoms[a].count;++child) inverse[atoms[a].first+child]=offsets[index]+child;
        }
        auto local=[&](GlobalCellIndex id) {
            const auto it=std::lower_bound(ids.begin(),ids.end(),id);
            if(it==ids.end() || *it!=id) throw std::logic_error("Local dependency endpoint is missing.");
            return inverse[static_cast<std::size_t>(it-ids.begin())];
        };
        LocalLayoutScore score;
        for(const auto& peer:input.peers)
        {
            std::vector<LocalCellIndex> send,receive;
            for(auto c:peer.send) send.push_back(local(c));
            for(auto c:peer.receive) receive.push_back(local(c));
            score.communicationRanges+=detail::layoutRangeCount(std::move(send),levelStarts)+detail::layoutRangeCount(std::move(receive),levelStarts);
        }
        for(const auto& block:input.interfaces)
        {
            std::vector<LocalCellIndex> members{local(block.coarse)};
            for(auto c:block.auxiliary) members.push_back(local(c));
            for(auto c:block.fine) members.push_back(local(c));
            score.interfaceRanges+=detail::layoutRangeCount(std::move(members),levelStarts);
        }
        for(const auto& face:boundaryRows)
        {
            std::vector<LocalCellIndex> members;
            for(const auto& entry:face) members.push_back(inverse[entry.first]);
            score.boundaryRanges+=detail::layoutRangeCount(std::move(members),levelStarts);
        }
        std::vector<LocalCellIndex> sendWork,receiveWork;
        for(std::size_t i=0;i<records.size();++i) if(records[i].owned && !records[i].ghost)
        {
            if(records[i].send) sendWork.push_back(inverse[i]);
            if(records[i].receiveDependent) receiveWork.push_back(inverse[i]);
        }
        score.communicationRanges+=detail::layoutRangeCount(std::move(sendWork),levelStarts)+detail::layoutRangeCount(std::move(receiveWork),levelStarts);
        const auto sourceCount=sourcePositions.size();
        if(ids.size()>1 && sourceCount>std::numeric_limits<std::uint64_t>::max()/Traits::nAnchors/(ids.size()-1))
            throw std::length_error("Local anchor-distance score exceeds uint64 range.");
        std::uint64_t distanceSum=0;
#pragma omp parallel for schedule(static) reduction(+:distanceSum) if(!omp_in_parallel())
        for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(sourceCount);++i)
        {
            const auto row=static_cast<std::size_t>(i),source=inverse[sourcePositions[row]];
            for(auto endpoint:endpointPositions[row])
            {
                const auto target=inverse[endpoint];
                distanceSum+=source>target ? source-target : target-source;
            }
        }
        score.anchorDistance=distanceSum;
        if(!hasBest || score<bestScore) { bestOrder=order; bestScore=score; selected=candidates; hasBest=true; }
        ++candidates;
    };
    for(unsigned typeFirst=0;typeFirst<(options.tryTypeFirst ? 2U : 1U);++typeFirst)
    {
        std::uint64_t bucketBits=0;
#pragma omp parallel for schedule(static) reduction(|:bucketBits) if(!omp_in_parallel())
        for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(atoms.size());++i)
            bucketBits|=bucket(atoms[static_cast<std::size_t>(i)],typeFirst!=0);
        unsigned bucketPasses=0; for(auto bits=bucketBits;bits;bits>>=8) ++bucketPasses;
        for(unsigned morton=0;morton<(options.tryMorton && mortonAvailable ? 2U : 1U);++morton)
        {
            std::vector<std::size_t> order(atoms.size()); std::iota(order.begin(),order.end(),0);
            if(morton) detail::layoutMortonOrder(order,atoms);
            detail::layoutRadixOrder(order,bucketPasses,[&](std::size_t a,unsigned pass) {
                return static_cast<unsigned>((bucket(atoms[a],typeFirst!=0)>>(pass*8))&255U);
            });
            evaluate(order);
        }
    }
    Result result; result.rank=rank; result.frame=input.frame; result.score=bestScore; result.candidates=candidates; result.selectedCandidate=selected;
    result.lookupGlobal=ids; result.lookupLocal.resize(ids.size()); result.cells.resize(ids.size()); result.coordinates.resize(ids.size());
    std::vector<std::size_t> finalOffsets(bestOrder.size()+1,0);
    for(std::size_t i=0;i<bestOrder.size();++i) finalOffsets[i+1]=finalOffsets[i]+atoms[bestOrder[i]].count;
#pragma omp parallel for schedule(static) if(!omp_in_parallel())
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(bestOrder.size());++i)
    {
        const auto index=static_cast<std::size_t>(i),a=bestOrder[index];
        for(std::size_t child=0;child<atoms[a].count;++child)
        {
            const auto source=atoms[a].first+child,target=finalOffsets[index]+child;
            result.lookupLocal[source]=target; result.cells[target]=records[source]; result.coordinates[target]=input.coordinates[source];
        }
    }
    auto local=[&](GlobalCellIndex id) {
        const auto value=result.localIndex(id);
        if(value==invalidCellIndex) throw std::logic_error("Final local dependency endpoint is missing.");
        return value;
    };
    for(const auto& peer:input.peers)
    {
        PeerExchange e; e.rank=peer.rank;
        for(auto c:peer.send) e.send.push_back(local(c));
        for(auto c:peer.receive) e.receive.push_back(local(c));
        result.peers.push_back(std::move(e));
    }
    std::vector<std::pair<LocalCellIndex,std::size_t>> sourceRows;
    sourceRows.reserve(sourcePositions.size());
    for(std::size_t row=0;row<sourcePositions.size();++row) sourceRows.push_back({result.lookupLocal[sourcePositions[row]],row});
    std::sort(sourceRows.begin(),sourceRows.end());
    result.sourceIDs.resize(sourceRows.size()); result.anchors.resize(sourceRows.size());
#pragma omp parallel for schedule(static) if(!omp_in_parallel())
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(sourceRows.size());++i)
    {
        const auto row=static_cast<std::size_t>(i),source=sourceRows[row].second;
        result.sourceIDs[row]=sourceRows[row].first;
        for(std::size_t a=0;a<Traits::nAnchors;++a) result.anchors[row][a]=result.lookupLocal[endpointPositions[source][a]];
    }
    for(const auto& block:input.interfaces)
    {
        LocalInterfaceBlock b; b.coarse=local(block.coarse);
        for(auto c:block.auxiliary) b.auxiliary.push_back(local(c));
        for(auto c:block.fine) b.fine.push_back(local(c));
        result.interfaces.push_back(std::move(b));
    }
    for(const auto& exchange:input.intraRankExchanges)
    {
        LocalBlockExchange e; e.a=exchange.a; e.b=exchange.b;
        for(const auto& pair:exchange.dependencies) e.dependencies.push_back({local(pair.first),local(pair.second)});
        result.intraRankExchanges.push_back(std::move(e));
    }
    constexpr auto qCount=lattice::Model<Type>::nStencils-1;
    result.boundaries.offsets.push_back(0);
    std::vector<std::size_t> qRows;
    std::size_t boundaryCount=0; for(const auto& face:boundaryRows) boundaryCount+=face.size();
    qRows.reserve(boundaryCount); result.boundaries.cells.reserve(boundaryCount);
    for(const auto& face:boundaryRows)
    {
        std::vector<std::pair<LocalCellIndex,std::size_t>> entries;
        for(const auto& entry:face) entries.push_back({result.lookupLocal[entry.first],entry.second});
        std::sort(entries.begin(),entries.end());
        const auto begin=result.boundaries.cells.size();
        for(const auto& entry:entries) { result.boundaries.cells.push_back(entry.first); qRows.push_back(entry.second); }
        result.boundaryFaces.push_back({begin,entries.size()});
        result.boundaries.offsets.push_back(result.boundaries.cells.size());
    }
    if(qRows.size()>std::numeric_limits<std::size_t>::max()/qCount)
        throw std::length_error("Local boundary Q count exceeds size range.");
    result.boundaries.q.resize(qRows.size()*qCount);
#pragma omp parallel for schedule(static) if(!omp_in_parallel())
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(qRows.size());++i)
    {
        const auto row=static_cast<std::size_t>(i);
        std::copy_n(globalBoundary.q.begin()+qRows[row]*qCount,qCount,result.boundaries.q.begin()+row*qCount);
    }
    auto selectedBy=[](const LocalLayoutCell& cell,LayoutGroup kind) {
        switch(kind)
        {
            case LayoutGroup::All:return true;
            case LayoutGroup::Owned:return cell.owned;
            case LayoutGroup::Computes:return cell.owned && !cell.ghost && cell.role!=geometry::Auxiliary && cell.type!=geometry::CellType::Dry;
            case LayoutGroup::Auxiliary:return !cell.ghost && cell.role==geometry::Auxiliary;
            case LayoutGroup::Ghost:return cell.ghost;
            case LayoutGroup::Halo:return cell.halo;
            case LayoutGroup::Send:return cell.send;
            case LayoutGroup::Receive:return cell.halo;
            case LayoutGroup::ReceiveDependent:return cell.receiveDependent;
            case LayoutGroup::Interior:return !cell.ghost && cell.type==geometry::CellType::Interior;
            case LayoutGroup::Boundary:return !cell.ghost && cell.type==geometry::CellType::Boundary;
            case LayoutGroup::FineInterface:return !cell.ghost && (cell.role&geometry::FineInterface)!=0 && cell.role!=geometry::Auxiliary;
            case LayoutGroup::CoarseInterface:return !cell.ghost && (cell.role&geometry::CoarseInterface)!=0 && cell.role!=geometry::Auxiliary;
        }
        return false;
    };
    result.levelRanges.resize(input.frame.levelCount);
    std::size_t levelCursor=0;
    for(std::size_t l=0;l<result.levelRanges.size();++l)
    {
        const auto begin=levelCursor;
        while(levelCursor<result.cells.size() && result.cells[levelCursor].level==l) ++levelCursor;
        result.levelRanges[l]={begin,levelCursor-begin};
    }
    for(std::size_t first=0;first<result.cells.size();)
    {
        auto last=first+1; const auto l=result.cells[first].level;
        while(last<result.cells.size() && result.cells[last].level==l) ++last;
        for(unsigned k=0;k<=static_cast<unsigned>(LayoutGroup::CoarseInterface);++k)
        {
            LocalLayoutGroup group; group.level=l; group.kind=static_cast<LayoutGroup>(k);
            group.ranges.offset=result.selectionRanges.size(); group.points.offset=result.selectionPoints.size();
            for(auto i=first;i<last;)
            {
                if(!selectedBy(result.cells[i],group.kind)) { ++i; continue; }
                const auto begin=i++;
                while(i<last && selectedBy(result.cells[i],group.kind)) ++i;
                if(i-begin==1) result.selectionPoints.push_back(begin);
                else result.selectionRanges.push_back({begin,i-begin});
            }
            group.ranges.length=result.selectionRanges.size()-group.ranges.offset;
            group.points.length=result.selectionPoints.size()-group.points.offset;
            result.groups.push_back(group);
        }
        first=last;
    }
    std::vector<geometry::NodeIndex> nodeKeys;
    for(const auto& node:input.nodes) nodeKeys.push_back(node.globalNode);
    auto nodeIndex=[&](geometry::NodeIndex global) {
        const auto it=std::lower_bound(nodeKeys.begin(),nodeKeys.end(),global);
        return it==nodeKeys.end() || *it!=global ? invalidCellIndex : static_cast<std::size_t>(it-nodeKeys.begin());
    };
    for(const auto& source:input.nodes)
    {
        LocalLayoutNode<D> node; static_cast<RankTreeNode<D>&>(node)=source;
        node.parent=nodeIndex(source.parentGlobal); node.realCell=result.localIndex(source.globalRealCell);
        for(std::size_t child=0;child<(std::size_t{1}<<D);++child) node.children[child]=nodeIndex(source.globalChildren[child]);
        result.nodes.push_back(node); result.treeCells.push_back({source.globalNode,node.realCell});
    }
    for(auto& cell:result.cells) if(!cell.ghost) cell.localSourceNode=nodeIndex(cell.globalSourceNode);
    return result;
}

template<std::size_t D,lattice::LatticeType Type>
auto LocalLayoutBuilder<D,Type>::build(int rank,const LocalLayoutOptions& options) const -> Result
{ return buildLocalLayout(prepare(rank),options); }

} // namespace ntic::lbm::mpi
