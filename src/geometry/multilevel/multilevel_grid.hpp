#pragma once

#include "refinement_region.hpp"

#include <atomic>
#include <cstdint>
#include <exception>
#include <tuple>
#include <omp.h>

namespace ntic::lbm::geometry
{

using NodeIndex = std::size_t;
using CellIndex = std::size_t;
inline constexpr std::size_t invalidIndex = std::numeric_limits<std::size_t>::max();

struct IndexRange
{
    std::size_t begin=0,end=0;
    std::size_t size() const { return end-begin; }
};

enum CellRole : std::uint8_t
{
    Ordinary=0,
    FineInterface=1,   // Real leaf, communicates with its coarser level.
    CoarseInterface=2, // Real leaf, communicates with its finer level.
    BothInterfaces=3,
    Auxiliary=4       // Fine-size exchange representation of a real coarse leaf.
};

struct MultilevelOptions
{
    unsigned maxLevel=6;
    std::size_t maxNodes=5000000;
    std::size_t maxCells=10000000;
    // Bound dense geometry analysis memory by sparse, occupied tiles rather
    // than allocating the entire finest Cartesian bounding box.
    unsigned analysisTileWidth=16;
    unsigned analysisThreads=8;
};

template<std::size_t D>
struct TreeNode
{
    std::array<std::uint32_t,D> coordinate{};
    NodeIndex parent=invalidIndex,firstChild=invalidIndex;
    CellIndex cell=invalidIndex;
    unsigned level=0;
    bool active=true;
    bool leaf() const { return active && firstChild==invalidIndex; }
};

// Only real leaves and auxiliary exchange cells have CellIndex slots. Covered
// parents remain in nodes, with cell == invalidIndex and no field-data slot.
struct CellPool
{
    std::vector<NodeIndex> sourceNode;
    std::vector<std::uint8_t> role,auxiliaryChild;
    std::vector<CellType> cellType;
    std::size_t size() const { return sourceNode.size(); }
    bool isAuxiliary(CellIndex c) const { return role.at(c)==Auxiliary; }
    bool computes(CellIndex c) const { return !isAuxiliary(c) && cellType.at(c)!=CellType::Dry; }
};

struct LevelRanges
{
    IndexRange all;
    // Ordinary, FineInterface, CoarseInterface, BothInterfaces, Auxiliary.
    std::array<IndexRange,5> roles{};
    // Real-role subdivisions: Interior, Boundary, Dry. Each computes range
    // is types[r][0].begin .. types[r][1].end, without duplication of cells.
    std::array<std::array<IndexRange,3>,4> types{};
};

struct InterfaceBlock
{
    CellIndex coarse=invalidIndex;
    IndexRange auxiliary; // Four/eight contiguous fine-size slots.
    IndexRange fine;      // Range into InterfaceConnectivity::fineCells.
};

struct InterfaceConnectivity
{
    std::vector<InterfaceBlock> blocks;
    std::vector<CellIndex> fineCells;
};

struct GridBoundaryData
{
    std::size_t nBoundaries=0;
    std::vector<std::size_t> offsets{0};
    std::vector<CellIndex> cells;
    std::vector<double> q; // entry-major; same convention as GeometryAnalysis.
};

template<std::size_t D, lattice::LatticeType Type =
    D==2 ? lattice::LatticeType::D2Q9 : lattice::LatticeType::D3Q19>
class MultilevelGrid
{
public:
    static constexpr std::size_t dimension=D;
    static constexpr lattice::LatticeType latticeType=Type;
    using Coordinate=std::array<std::uint32_t,D>;
    using SignedCoordinate=std::array<std::int64_t,D>;
private:
    static_assert(D==2 || D==3,"Multilevel grids support 2D and 3D.");
    using Analysis=std::conditional_t<D==2,GeometryAnalysis2D<Type>,GeometryAnalysis3D<Type>>;
    using Lattice=typename Analysis::Lattice;
    static constexpr std::size_t childrenPerCell=std::size_t{1}<<D;
    static constexpr std::size_t links=Lattice::nStencils-1;
    struct BoundaryEntry { std::size_t face; CellIndex cell; std::array<double,links> q; };

    RegionBox<D> domain_;
    double coarseSpacing_;
    MultilevelOptions options_;
    Coordinate rootDimensions_{};
    std::vector<TreeNode<D>> nodes_;
    CellPool cells_;
    std::vector<LevelRanges> levels_;
    InterfaceConnectivity interfaces_;
    GridBoundaryData boundaries_;
    bool finalized_=false;

    static bool mortonLess(const Coordinate& a,const Coordinate& b)
    {
        // Compare interleaved bits without a fixed-width packed Morton key.
        std::uint32_t difference=0;
        for(std::size_t d=0; d<D; ++d) difference|=a[d]^b[d];
        if(difference==0) return false;
        unsigned bit=0; while(difference>>1U) { difference>>=1U; ++bit; }
        for(std::size_t axis=D; axis-->0;)
        {
            const auto x=(a[axis]>>bit)&1U,y=(b[axis]>>bit)&1U;
            if(x!=y) return x<y;
        }
        return false;
    }
    void invalidate()
    {
        finalized_=false; cells_=CellPool{}; levels_.clear();
        interfaces_=InterfaceConnectivity{}; boundaries_=GridBoundaryData{};
        for(auto& node:nodes_) node.cell=invalidIndex;
    }
    void requireFinalized() const
    { if(!finalized_) throw std::logic_error("Finalize the grid after topology edits before accessing cell data."); }

    void split(const std::vector<NodeIndex>& parents)
    {
        if(parents.empty()) return;
        if(parents.size()>(options_.maxNodes-nodes_.size())/childrenPerCell)
            throw std::length_error("Multilevel topology exceeds maxNodes.");
        for(auto p:parents)
            if(!nodes_.at(p).leaf() || nodes_[p].level>=options_.maxLevel)
                throw std::invalid_argument("Only leaves below maxLevel can be refined.");
        const auto first=nodes_.size();
        nodes_.resize(first+parents.size()*childrenPerCell);
#pragma omp parallel for schedule(static)
        for(std::ptrdiff_t i=0; i<static_cast<std::ptrdiff_t>(parents.size()); ++i)
        {
            const auto p=parents[static_cast<std::size_t>(i)];
            const auto begin=first+static_cast<std::size_t>(i)*childrenPerCell;
            const auto coordinate=nodes_[p].coordinate;
            const unsigned level=nodes_[p].level+1;
            nodes_[p].firstChild=begin;
            for(std::size_t child=0; child<childrenPerCell; ++child)
            {
                auto& n=nodes_[begin+child]; n.parent=p; n.level=level;
                for(std::size_t d=0; d<D; ++d) n.coordinate[d]=2*coordinate[d]+((child>>d)&1U);
            }
        }
    }
    template<typename Function>
    static void offsets(int width,Function consume)
    {
        SignedCoordinate delta{};
        for(int x=-width; x<=width; ++x)
            for(int y=-width; y<=width; ++y)
            {
                delta[0]=x; delta[1]=y;
                if constexpr(D==2) consume(delta);
                else for(int z=-width; z<=width; ++z) { delta[2]=z; consume(delta); }
            }
    }
    NodeIndex neighbourAt(const TreeNode<D>& node,const SignedCoordinate& offset) const
    {
        SignedCoordinate coordinate{};
        for(std::size_t d=0; d<D; ++d) coordinate[d]=static_cast<std::int64_t>(node.coordinate[d])+offset[d];
        return coveringNode(node.level,coordinate);
    }
    void balance()
    {
        // Full face/edge/corner 2:1 balance, independent of a particular lattice.
        for(;;)
        {
            // One byte per node avoids collecting many duplicate requests from
            // fine neighbours. Atomic stores only mark coarse target nodes.
            std::vector<std::uint8_t> marked(nodes_.size(),0);
#pragma omp parallel for schedule(static)
            for(std::ptrdiff_t i=0; i<static_cast<std::ptrdiff_t>(nodes_.size()); ++i)
            {
                const auto& node=nodes_[static_cast<std::size_t>(i)];
                if(!node.leaf() || node.level<2) continue;
                offsets(1,[&](const auto& delta) {
                    const auto neighbour=neighbourAt(node,delta);
                    if(neighbour!=invalidIndex && nodes_[neighbour].level+1<node.level)
                    {
#pragma omp atomic write
                        marked[neighbour]=1;
                    }
                });
            }
            std::vector<NodeIndex> marks;
            for(NodeIndex i=0; i<marked.size(); ++i) if(marked[i]) marks.push_back(i);
            if(marks.empty()) break;
            split(marks);
        }
    }
    std::vector<std::uint8_t> interfaceRoles(const RefinementRegion<D>& physical) const
    {
        std::vector<std::uint8_t> roles(nodes_.size(),Ordinary);
        std::exception_ptr error;
#pragma omp parallel for schedule(static)
        for(std::ptrdiff_t i=0; i<static_cast<std::ptrdiff_t>(nodes_.size()); ++i)
        {
            const auto& node=nodes_[static_cast<std::size_t>(i)];
            if(!node.leaf()) continue;
            bool coarse=false;
            offsets(1,[&](const auto& delta) {
                const auto other=neighbourAt(node,delta);
                coarse=coarse || (other!=invalidIndex && nodes_[other].level==node.level && !nodes_[other].leaf());
            });
            if(coarse)
            {
                try
                {
                    // Do not allocate exchange families in space proven wholly
                    // outside the fluid. A dry center alone is NOT sufficient
                    // to discard a cut coarse cell or its wet auxiliary cells.
                    if(physical.classify(nodeBox(static_cast<std::size_t>(i)))!=BoxRelation::Outside)
                        roles[static_cast<std::size_t>(i)]=CoarseInterface;
                }
                catch(...)
                {
#pragma omp critical(ntic_multilevel_interface_error)
                    { if(!error) error=std::current_exception(); }
                }
            }
        }
        if(error) std::rethrow_exception(error);
        // Barrier between the passes: only the immutable coarse flags are read.
        const auto coarse=roles;
#pragma omp parallel for schedule(static)
        for(std::ptrdiff_t i=0; i<static_cast<std::ptrdiff_t>(nodes_.size()); ++i)
        {
            const auto& node=nodes_[static_cast<std::size_t>(i)];
            if(!node.leaf() || node.level==0) continue;
            bool fine=false;
            offsets(2,[&](const auto& delta) {
                const auto other=neighbourAt(node,delta);
                fine=fine || (other!=invalidIndex && nodes_[other].level+1==node.level && coarse[other]==CoarseInterface);
            });
            if(fine) roles[static_cast<std::size_t>(i)]|=FineInterface;
        }
        return roles;
    }
    void appendCell(NodeIndex source,std::uint8_t role,std::uint8_t child=0)
    {
        if(cells_.size()>=options_.maxCells) throw std::length_error("Multilevel cell pool exceeds maxCells.");
        cells_.sourceNode.push_back(source); cells_.role.push_back(role);
        cells_.auxiliaryChild.push_back(child); cells_.cellType.push_back(CellType::Dry);
    }

    template<typename Geometry>
    std::vector<BoundaryEntry> analyzeTiles(const Geometry& geometry)
    {
        struct Tile { unsigned level; Coordinate coordinate; IndexRange entries; };
        std::vector<CellIndex> order(cells_.size()); std::iota(order.begin(),order.end(),0);
        const auto width=options_.analysisTileWidth;
        const auto tileCoordinate=[&](CellIndex c) {
            auto key=coordinate(c); for(auto& value:key) value/=width; return key;
        };
        std::sort(order.begin(),order.end(),[&](auto a,auto b) {
            if(level(a)!=level(b)) return level(a)<level(b);
            const auto x=tileCoordinate(a),y=tileCoordinate(b);
            return x<y || (x==y && coordinate(a)<coordinate(b));
        });
        std::vector<Tile> tiles;
        for(std::size_t first=0; first<order.size();)
        {
            const auto l=level(order[first]); const auto key=tileCoordinate(order[first]);
            auto last=first+1;
            while(last<order.size() && level(order[last])==l && tileCoordinate(order[last])==key) ++last;
            tiles.push_back({l,key,{first,last}}); first=last;
        }
        const int threads=std::max(1,std::min(omp_get_max_threads(),static_cast<int>(options_.analysisThreads)));
        std::vector<std::vector<BoundaryEntry>> output(static_cast<std::size_t>(threads));
        std::exception_ptr error; std::atomic<bool> failed{false};
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        for(std::ptrdiff_t t=0; t<static_cast<std::ptrdiff_t>(tiles.size()); ++t)
        {
            if(failed.load(std::memory_order_relaxed)) continue;
            try
            {
                const auto& tile=tiles[static_cast<std::size_t>(t)]; const double h=spacing(tile.level);
                Coordinate start{},last{};
                for(std::size_t d=0; d<D; ++d)
                {
                    start[d]=tile.coordinate[d]*width;
                    last[d]=std::min<std::uint64_t>(static_cast<std::uint64_t>(start[d])+width,
                        static_cast<std::uint64_t>(rootDimensions_[d])<<tile.level);
                }
                RegionBox<D> box;
                for(std::size_t d=0; d<D; ++d)
                {
                    box.min[d]=domain_.min[d]+static_cast<double>(start[d])*h;
                    box.max[d]=domain_.min[d]+static_cast<double>(last[d])*h;
                }
                Analysis analysis(h); analysis.domainBounds=box;
                geometry.analysis(analysis); // Reuse primitive/STL/Boolean scalar, Face and Q contracts.
                if(analysis.nBoundaries!=detail::GeometryQuery<Geometry>::faceCount(geometry))
                    throw std::logic_error("Geometry analysis returned inconsistent Face groups.");
                std::vector<CellIndex> map(analysis.scalar.size(),invalidIndex);
                const auto linear=[&](const Coordinate& p) {
                    std::size_t id=p[0]-start[0]+analysis.nx*(p[1]-start[1]);
                    if constexpr(D==3) id+=analysis.nx*analysis.ny*(p[2]-start[2]);
                    return id;
                };
                for(std::size_t i=tile.entries.begin; i<tile.entries.end; ++i)
                {
                    const auto c=order[i],id=linear(coordinate(c));
                    if(id>=map.size() || map[id]!=invalidIndex) throw std::logic_error("Invalid or duplicate grid cell.");
                    map[id]=c; cells_.cellType[c]=static_cast<CellType>(static_cast<int>(analysis.scalar[id]));
                }
                auto& entries=output[static_cast<std::size_t>(omp_get_thread_num())];
                for(std::size_t face=0; face<analysis.nBoundaries; ++face)
                    for(std::size_t entry=analysis.boundaryOffsets[face]; entry<analysis.boundaryOffsets[face+1]; ++entry)
                    {
                        Coordinate p{};
                        p[0]=start[0]+static_cast<std::uint32_t>(std::llround((analysis.boundaryX[entry]-box.min[0])/h-0.5));
                        p[1]=start[1]+static_cast<std::uint32_t>(std::llround((analysis.boundaryY[entry]-box.min[1])/h-0.5));
                        if constexpr(D==3) p[2]=start[2]+static_cast<std::uint32_t>(std::llround((analysis.boundaryZ[entry]-box.min[2])/h-0.5));
                        const auto id=linear(p);
                        if(id>=map.size() || map[id]==invalidIndex) continue;
                        BoundaryEntry record{face,map[id],{}};
                        std::copy_n(analysis.q.begin()+entry*links,links,record.q.begin());
                        entries.push_back(record);
                    }
            }
            catch(...)
            {
                failed.store(true,std::memory_order_relaxed);
#pragma omp critical(ntic_multilevel_analysis_error)
                { if(!error) error=std::current_exception(); }
            }
        }
        if(error) std::rethrow_exception(error);
        std::vector<BoundaryEntry> result;
        std::size_t count=0; for(const auto& list:output) count+=list.size(); result.reserve(count);
        for(auto& list:output) result.insert(result.end(),std::make_move_iterator(list.begin()),std::make_move_iterator(list.end()));
        return result;
    }
    static unsigned typeOrder(CellType type)
    { return type==CellType::Interior ? 0 : type==CellType::Boundary ? 1 : 2; }

    void sortCells(std::vector<BoundaryEntry>& boundaryEntries)
    {
        std::vector<CellIndex> order(cells_.size()),newIndex(cells_.size()); std::iota(order.begin(),order.end(),0);
        std::sort(order.begin(),order.end(),[&](auto a,auto b) {
            if(level(a)!=level(b)) return level(a)<level(b);
            if(cells_.role[a]!=cells_.role[b]) return cells_.role[a]<cells_.role[b];
            if(!cells_.isAuxiliary(a) && cells_.cellType[a]!=cells_.cellType[b])
                return typeOrder(cells_.cellType[a])<typeOrder(cells_.cellType[b]);
            if(cells_.isAuxiliary(a))
            {
                const auto& x=nodes_[cells_.sourceNode[a]].coordinate;
                const auto& y=nodes_[cells_.sourceNode[b]].coordinate;
                if(x!=y) return mortonLess(x,y);
                return cells_.auxiliaryChild[a]<cells_.auxiliaryChild[b];
            }
            return mortonLess(coordinate(a),coordinate(b));
        });
        CellPool sorted;
        sorted.sourceNode.reserve(order.size()); sorted.role.reserve(order.size());
        sorted.auxiliaryChild.reserve(order.size()); sorted.cellType.reserve(order.size());
        for(std::size_t i=0; i<order.size(); ++i)
        {
            const auto old=order[i]; newIndex[old]=i;
            sorted.sourceNode.push_back(cells_.sourceNode[old]); sorted.role.push_back(cells_.role[old]);
            sorted.auxiliaryChild.push_back(cells_.auxiliaryChild[old]); sorted.cellType.push_back(cells_.cellType[old]);
        }
        cells_=std::move(sorted);
        for(std::size_t c=0; c<cells_.size(); ++c)
            if(!cells_.isAuxiliary(c)) nodes_[cells_.sourceNode[c]].cell=c;
        for(auto& record:boundaryEntries) record.cell=newIndex[record.cell];
        unsigned maxLevel=0; for(std::size_t c=0; c<cells_.size(); ++c) maxLevel=std::max(maxLevel,level(c));
        levels_.resize(maxLevel+1);
        std::size_t cursor=0;
        for(unsigned l=0; l<=maxLevel; ++l)
        {
            auto& ranges=levels_[l]; ranges.all.begin=cursor;
            for(unsigned role=0; role<5; ++role)
            {
                ranges.roles[role].begin=cursor;
                if(role<4)
                    for(unsigned type=0; type<3; ++type)
                    {
                        auto& range=ranges.types[role][type]; range.begin=cursor;
                        while(cursor<cells_.size() && level(cursor)==l && cells_.role[cursor]==role &&
                              typeOrder(cells_.cellType[cursor])==type) ++cursor;
                        range.end=cursor;
                    }
                else while(cursor<cells_.size() && level(cursor)==l && cells_.role[cursor]==Auxiliary) ++cursor;
                ranges.roles[role].end=cursor;
            }
            ranges.all.end=cursor;
        }
        if(cursor!=cells_.size()) throw std::logic_error("Cell range construction failed.");
    }
    void buildInterfaces()
    {
        std::vector<std::pair<NodeIndex,CellIndex>> fine;
        for(CellIndex c=0; c<cells_.size(); ++c)
        {
            if(cells_.isAuxiliary(c) || !(cells_.role[c]&FineInterface)) continue;
            const auto& node=nodes_[cells_.sourceNode[c]];
            std::vector<NodeIndex> sources;
            offsets(2,[&](const auto& delta) {
                const auto other=neighbourAt(node,delta);
                if(other!=invalidIndex && nodes_[other].level+1==node.level && nodes_[other].cell!=invalidIndex &&
                   (cells_.role[nodes_[other].cell]&CoarseInterface)) sources.push_back(other);
            });
            std::sort(sources.begin(),sources.end()); sources.erase(std::unique(sources.begin(),sources.end()),sources.end());
            for(auto source:sources) fine.emplace_back(source,c);
        }
        std::sort(fine.begin(),fine.end()); fine.erase(std::unique(fine.begin(),fine.end()),fine.end());
        for(const auto& ranges:levels_)
            for(auto first=ranges.roles[Auxiliary].begin; first<ranges.roles[Auxiliary].end; first+=childrenPerCell)
            {
                const auto source=cells_.sourceNode[first];
                const auto begin=std::lower_bound(fine.begin(),fine.end(),std::make_pair(source,CellIndex{0}));
                const auto end=std::upper_bound(fine.begin(),fine.end(),std::make_pair(source,invalidIndex));
                const auto offset=interfaces_.fineCells.size();
                for(auto it=begin; it!=end; ++it) interfaces_.fineCells.push_back(it->second);
                interfaces_.blocks.push_back({nodes_[source].cell,{first,first+childrenPerCell},
                    {offset,interfaces_.fineCells.size()}});
            }
    }
    void flattenBoundary(std::vector<BoundaryEntry>& entries,std::size_t faces)
    {
        std::sort(entries.begin(),entries.end(),[](const auto& a,const auto& b) {
            return std::tie(a.face,a.cell)<std::tie(b.face,b.cell);
        });
        boundaries_.nBoundaries=faces; boundaries_.offsets.assign(faces+1,0);
        if(entries.size()>std::numeric_limits<std::size_t>::max()/links) throw std::length_error("Boundary Q size overflow.");
        boundaries_.cells.resize(entries.size()); boundaries_.q.resize(entries.size()*links);
        for(std::size_t i=0; i<entries.size(); ++i)
        {
            boundaries_.cells[i]=entries[i].cell; ++boundaries_.offsets[entries[i].face+1];
            std::copy(entries[i].q.begin(),entries[i].q.end(),boundaries_.q.begin()+i*links);
        }
        for(std::size_t f=0; f<faces; ++f) boundaries_.offsets[f+1]+=boundaries_.offsets[f];
    }
    void collectTouching(NodeIndex current,NodeIndex source,std::vector<CellIndex>& output) const
    {
        const auto& a=nodes_[current]; const auto& b=nodes_[source];
        const unsigned common=std::max(a.level,b.level);
        bool boundary=false;
        for(std::size_t d=0; d<D; ++d)
        {
            const auto amin=static_cast<std::uint64_t>(a.coordinate[d])<<(common-a.level);
            const auto amax=amin+(std::uint64_t{1}<<(common-a.level));
            const auto bmin=static_cast<std::uint64_t>(b.coordinate[d])<<(common-b.level);
            const auto bmax=bmin+(std::uint64_t{1}<<(common-b.level));
            if(amax<bmin || bmax<amin) return;
            boundary=boundary || amax==bmin || bmax==amin;
        }
        if(a.leaf()) { if(boundary && a.cell!=invalidIndex) output.push_back(a.cell); return; }
        for(std::size_t child=0; child<childrenPerCell; ++child) collectTouching(a.firstChild+child,source,output);
    }
public:
    MultilevelGrid(const RegionBox<D>& domain,double spacing,MultilevelOptions options={})
        : domain_(domain),coarseSpacing_(spacing),options_(options)
    {
        if(!(spacing>0) || !std::isfinite(spacing) || options.maxLevel>30 || options.analysisTileWidth==0 ||
           options.analysisTileWidth>64 || options.analysisThreads==0 || options.analysisThreads>1024 || options.maxNodes==0 || options.maxCells==0)
            throw std::invalid_argument("Invalid multilevel grid options.");
        std::size_t count=1;
        for(std::size_t d=0; d<D; ++d)
        {
            const double extent=domain.max[d]-domain.min[d],n=std::ceil(extent/spacing);
            if(!std::isfinite(domain.min[d]) || !std::isfinite(domain.max[d]) || !(extent>0) ||
               !std::isfinite(n) || n<1 || n>static_cast<double>(std::numeric_limits<std::uint32_t>::max()>>options.maxLevel))
                throw std::invalid_argument("Grid domain/level exceeds integer coordinate capacity.");
            rootDimensions_[d]=static_cast<std::uint32_t>(n);
            if(count>options.maxNodes/rootDimensions_[d]) throw std::length_error("Root grid exceeds maxNodes.");
            count*=rootDimensions_[d];
            // Use the snapped outer grid domain consistently on every level.
            domain_.max[d]=domain_.min[d]+rootDimensions_[d]*spacing;
            const double finest=std::ldexp(spacing,-static_cast<int>(options.maxLevel));
            if(!(finest>0) || domain_.min[d]+0.5*finest==domain_.min[d] || domain_.max[d]-0.5*finest==domain_.max[d])
                throw std::invalid_argument("Grid coordinates cannot represent the finest cell centers.");
        }
        nodes_.resize(count);
        for(std::size_t i=0; i<count; ++i)
        {
            auto remaining=i;
            for(std::size_t d=0; d<D; ++d) { nodes_[i].coordinate[d]=remaining%rootDimensions_[d]; remaining/=rootDimensions_[d]; }
        }
    }
    const RegionBox<D>& domain() const { return domain_; }
    double spacing(unsigned level) const { return std::ldexp(coarseSpacing_,-static_cast<int>(level)); }
    const std::vector<TreeNode<D>>& nodes() const { return nodes_; }
    const CellPool& cells() const { requireFinalized(); return cells_; }
    const std::vector<LevelRanges>& levels() const { requireFinalized(); return levels_; }
    const InterfaceConnectivity& interfaces() const { requireFinalized(); return interfaces_; }
    const GridBoundaryData& boundaries() const { requireFinalized(); return boundaries_; }
    bool finalized() const { return finalized_; }
    unsigned level(CellIndex cell) const
    { return nodes_.at(cells_.sourceNode.at(cell)).level+(cells_.role[cell]==Auxiliary ? 1 : 0); }
    Coordinate coordinate(CellIndex cell) const
    {
        auto p=nodes_.at(cells_.sourceNode.at(cell)).coordinate;
        if(cells_.role[cell]==Auxiliary)
            for(std::size_t d=0; d<D; ++d) p[d]=2*p[d]+((cells_.auxiliaryChild[cell]>>d)&1U);
        return p;
    }
    Position<D> center(CellIndex cell) const
    {
        const auto p=coordinate(cell); const auto h=spacing(level(cell)); Position<D> result{};
        for(std::size_t d=0; d<D; ++d) result[d]=domain_.min[d]+(p[d]+0.5)*h;
        return result;
    }
    RegionBox<D> nodeBox(NodeIndex index) const
    {
        const auto& node=nodes_.at(index); const auto h=spacing(node.level); RegionBox<D> result;
        for(std::size_t d=0; d<D; ++d)
        {
            result.min[d]=domain_.min[d]+node.coordinate[d]*h;
            result.max[d]=domain_.min[d]+(static_cast<double>(node.coordinate[d])+1)*h;
        }
        return result;
    }
    // Return the node at the requested level, or its covering coarser leaf.
    NodeIndex coveringNode(unsigned level,const SignedCoordinate& coordinate) const
    {
        if(level>options_.maxLevel) return invalidIndex;
        NodeIndex root=0; std::size_t stride=1;
        for(std::size_t d=0; d<D; ++d)
        {
            if(coordinate[d]<0 || static_cast<std::uint64_t>(coordinate[d])>=(static_cast<std::uint64_t>(rootDimensions_[d])<<level)) return invalidIndex;
            root+=static_cast<std::size_t>(coordinate[d]>>level)*stride; stride*=rootDimensions_[d];
        }
        auto index=root;
        for(unsigned l=0; l<level && !nodes_[index].leaf(); ++l)
        {
            std::size_t child=0;
            for(std::size_t d=0; d<D; ++d) child|=((coordinate[d]>>(level-l-1))&1)<<d;
            index=nodes_[index].firstChild+child;
        }
        return index;
    }
    void refine(NodeIndex leaf)
    {
        if(!nodes_.at(leaf).leaf() || nodes_[leaf].level>=options_.maxLevel)
            throw std::invalid_argument("Invalid refinement request.");
        invalidate(); split({leaf});
    }
    void derefine(NodeIndex parent)
    {
        const auto& node=nodes_.at(parent);
        if(!node.active || node.leaf()) throw std::invalid_argument("derefine requires an active parent.");
        const auto first=node.firstChild;
        for(std::size_t child=0; child<childrenPerCell; ++child)
            if(!nodes_[first+child].leaf()) throw std::invalid_argument("derefine requires direct leaf children.");
        invalidate();
        for(std::size_t child=0; child<childrenPerCell; ++child) nodes_[first+child].active=false;
        nodes_[parent].firstChild=invalidIndex;
    }
    template<typename Geometry>
    void generate(const Geometry& geometry,const std::vector<RefinementRule<D>>& rules={})
    {
        for(const auto& rule:rules)
        {
            if(rule.targetLevel>options_.maxLevel) throw std::invalid_argument("Refinement target exceeds maxLevel.");
            if(rule.surfaceBand && (!(rule.thickness>=0) || !std::isfinite(rule.thickness)))
                throw std::invalid_argument("Invalid surface-band thickness.");
            for(auto face:rule.faces) if(face>=rule.region.faceCount()) throw std::invalid_argument("Invalid refinement Face ID.");
        }
        invalidate();
        const auto physical=fluidRegion(geometry);
        // Geometry selector failures must leave no finalized/partially valid result.
        for(;;)
        {
            std::vector<std::uint8_t> marked(nodes_.size(),0);
            std::exception_ptr error; std::atomic<bool> failed{false};
#pragma omp parallel for schedule(dynamic,64)
            for(std::ptrdiff_t i=0; i<static_cast<std::ptrdiff_t>(nodes_.size()); ++i)
            {
                const auto id=static_cast<std::size_t>(i); const auto& node=nodes_[id];
                if(!node.leaf() || node.level>=options_.maxLevel || failed.load(std::memory_order_relaxed)) continue;
                try
                {
                    const auto box=nodeBox(id);
                    bool requested=false;
                    for(const auto& rule:rules)
                        if(rule.targetLevel>node.level && rule.intersects(box)) { requested=true; break; }
                    if(requested && physical.classify(box)!=BoxRelation::Outside) marked[id]=1;
                }
                catch(...)
                {
                    failed.store(true,std::memory_order_relaxed);
#pragma omp critical(ntic_multilevel_selection_error)
                    { if(!error) error=std::current_exception(); }
                }
            }
            if(error) std::rethrow_exception(error);
            std::vector<NodeIndex> parents;
            for(NodeIndex i=0; i<marked.size(); ++i) if(marked[i]) parents.push_back(i);
            if(parents.empty()) break;
            split(parents);
        }
        finalize(geometry);
    }
    template<typename Geometry>
    void finalize(const Geometry& geometry)
    {
        invalidate(); balance();
        const auto roles=interfaceRoles(fluidRegion(geometry));
        for(NodeIndex n=0; n<nodes_.size(); ++n)
            if(nodes_[n].leaf())
            {
                appendCell(n,roles[n]);
                if(roles[n]&CoarseInterface)
                    for(std::size_t child=0; child<childrenPerCell; ++child) appendCell(n,Auxiliary,static_cast<std::uint8_t>(child));
            }
        try
        {
            auto entries=analyzeTiles(geometry);
            sortCells(entries); buildInterfaces();
            flattenBoundary(entries,detail::GeometryQuery<Geometry>::faceCount(geometry));
            finalized_=true;
        }
        catch(...) { invalidate(); throw; }
    }
    // Partition module consumes these read-only views; no rank/device ownership
    // is embedded in geometric indices. Adjacency is generated on demand rather
    // than storing a dense N*27 table. Auxiliary cells are represented by blocks.
    std::vector<CellIndex> adjacentRealCells(CellIndex cell) const
    {
        requireFinalized();
        if(cells_.isAuxiliary(cell)) throw std::invalid_argument("Request adjacency for a real leaf.");
        const auto source=cells_.sourceNode.at(cell); const auto& node=nodes_[source];
        std::vector<CellIndex> result;
        offsets(1,[&](const auto& delta) {
            const auto other=neighbourAt(node,delta);
            if(other!=invalidIndex && other!=source) collectTouching(other,source,result);
        });
        std::sort(result.begin(),result.end()); result.erase(std::unique(result.begin(),result.end()),result.end());
        return result;
    }
};

} // namespace ntic::lbm::geometry
