#include "partition_mpi.hpp"
#include "test_partition_checks.hpp"

#include <iostream>
#include <set>
#include <memory>
#include <optional>

namespace geo=ntic::lbm::geometry;
namespace part=ntic::lbm::mpi;
namespace lat=ntic::lbm::lattice;

namespace
{
constexpr auto highMorton=part::detail::layoutMortonKey(std::array<std::uint32_t,3>{0,0,1U<<31});
static_assert(highMorton[0]==(1U<<31) && highMorton[1]==0 && highMorton[2]==0);
static_assert(part::detail::layoutMortonKey(std::array<std::uint32_t,3>{1U<<22,0,0})[0]!=0);

void radixOrder()
{
    struct Atom { std::array<std::uint32_t,3> morton{}; };
    std::vector<Atom> atoms(5003);
    for(std::size_t i=0;i<atoms.size();++i)
        atoms[i].morton=part::detail::layoutMortonKey(std::array<std::uint32_t,3>{
            static_cast<std::uint32_t>((i%113)*7919),static_cast<std::uint32_t>((i%7)<<29),
            static_cast<std::uint32_t>((i%5)<<30)});
    // Deliberate duplicate high-bit keys exercise stable tie order.
    atoms[4000]=atoms[3]; atoms[4999]=atoms[3];
    std::vector<std::size_t> expected(atoms.size()); std::iota(expected.begin(),expected.end(),0);
    std::stable_sort(expected.begin(),expected.end(),[&](auto a,auto b) { return atoms[a].morton<atoms[b].morton; });
    auto serial=expected; std::iota(serial.begin(),serial.end(),0);
    omp_set_num_threads(1); part::detail::layoutMortonOrder(serial,atoms);
    auto parallel=serial; std::iota(parallel.begin(),parallel.end(),0);
    omp_set_num_threads(4); part::detail::layoutMortonOrder(parallel,atoms);
    requirePartition(serial==expected && parallel==expected,"96-bit radix order differs from independent stable sort");
    auto nested=serial; std::iota(nested.begin(),nested.end(),0);
#pragma omp parallel num_threads(4)
    {
#pragma omp single
        part::detail::layoutMortonOrder(nested,atoms);
    }
    requirePartition(nested==expected,"Nested radix order differs");
}

template<typename Layout>
std::vector<part::LocalCellIndex> selection(const Layout& layout,const part::LocalLayoutGroup& group)
{
    std::vector<part::LocalCellIndex> result;
    for(auto i=group.ranges.offset;i<group.ranges.offset+group.ranges.length;++i)
    {
        const auto range=layout.selectionRanges.at(i);
        requirePartition(range.length>1,"Singleton range should be an explicit point");
        for(auto c=range.offset;c<range.offset+range.length;++c) result.push_back(c);
    }
    for(auto i=group.points.offset;i<group.points.offset+group.points.length;++i) result.push_back(layout.selectionPoints.at(i));
    std::sort(result.begin(),result.end());
    requirePartition(std::adjacent_find(result.begin(),result.end())==result.end(),"A logical selection contains duplicates");
    return result;
}

template<std::size_t D,lat::LatticeType Type>
void checkLayout(const geo::MultilevelGrid<D,Type>& grid,const geo::NeighborLinks<D,Type>& links,
                 const part::PartitionPlan& plan,const part::LocalLayout<D,Type>& layout)
{
    std::set<std::size_t> originals; unsigned lastLevel=0;
    for(std::size_t c=0;c<layout.cells.size();++c)
    {
        const auto& cell=layout.cells[c];
        requirePartition(c==0 || cell.level>=lastLevel,"Levels are not contiguous"); lastLevel=cell.level;
        requirePartition(layout.localIndex(cell.global)==c && originals.insert(cell.global).second,"Local lookup is not one-to-one");
        requirePartition(cell.owned==(cell.owner==layout.rank) && cell.halo!=cell.owned && cell.ghost==(cell.global>=links.originalCount),"Storage classifications are inconsistent");
        if(!cell.ghost)
        {
            requirePartition(cell.owner==plan.owner[cell.global] && cell.globalSourceNode==grid.cells().sourceNode[cell.global],"Physical owner/source node changed");
            if(cell.role==geo::Auxiliary && grid.cells().auxiliaryChild[cell.global]==0)
                for(std::size_t child=0;child<(std::size_t{1}<<D);++child)
                    requirePartition(c+child<layout.cells.size() && layout.cells[c+child].global==cell.global+child,
                                     "Auxiliary family is split or its child order changed");
        }
        else
        {
            int owner=plan.ranks();
            for(std::size_t source=0;source<links.originalCount;++source)
                if(std::find(links.anchors[source].begin(),links.anchors[source].end(),cell.global)!=links.anchors[source].end())
                    owner=std::min(owner,plan.owner[source]);
            requirePartition(cell.owner==owner,"Ghost metadata owner is not the minimum referring rank");
        }
    }
    std::size_t covered=0;
    for(std::size_t level=0;level<layout.levelRanges.size();++level)
    {
        const auto span=layout.levelRanges[level];
        requirePartition(span.offset==covered,"Level ranges overlap or leave gaps");
        for(auto c=span.offset;c<span.offset+span.length;++c) requirePartition(layout.cells.at(c).level==level,"Wrong level span");
        covered+=span.length;
    }
    requirePartition(covered==layout.cells.size(),"Level ranges do not cover storage");
    std::size_t ownedSources=0;
    for(auto owner:plan.owner) ownedSources+=owner==layout.rank;
    requirePartition(layout.sourceIDs.size()==ownedSources && layout.anchors.size()==ownedSources,"Neighbor rows include an endpoint or omit an owned source");
    std::uint64_t anchorDistance=0;
    for(std::size_t row=0;row<layout.sourceIDs.size();++row)
    {
        const auto local=layout.sourceIDs[row],global=layout.cells[local].global;
        for(std::size_t a=0;a<lat::Twist<Type>::nAnchors;++a)
        {
            const auto endpoint=layout.anchors[row][a];
            requirePartition(layout.cells.at(endpoint).global==links.anchors[global][a],"Anchor address changed after local sorting");
            anchorDistance+=local>endpoint ? local-endpoint : endpoint-local;
        }
        for(std::size_t q=0;q<lat::Model<Type>::nStencils;++q)
        {
            const auto a=lat::Twist<Type>::directionAnchor[q];
            requirePartition(layout.anchor(local,q)==(a==0 ? local : layout.anchors[row][a-1]),"Local direction lookup is incorrect");
        }
    }
    requirePartition(anchorDistance==layout.score.anchorDistance,"Reported anchor-locality score differs from final addresses");
    for(const auto& group:layout.groups)
    {
        const auto values=selection(layout,group);
        for(auto c:values) requirePartition(layout.cells.at(c).level==group.level,"Selection crossed a level boundary");
        std::vector<std::size_t> expected;
        for(std::size_t c=0;c<layout.cells.size();++c)
        {
            const auto& cell=layout.cells[c]; if(cell.level!=group.level) continue;
            bool member=false;
            switch(group.kind)
            {
                case part::LayoutGroup::All:member=true;break;
                case part::LayoutGroup::Owned:member=cell.owned;break;
                case part::LayoutGroup::Computes:member=cell.owned && !cell.ghost && grid.cells().computes(cell.global);break;
                case part::LayoutGroup::Auxiliary:member=!cell.ghost && grid.cells().isAuxiliary(cell.global);break;
                case part::LayoutGroup::Ghost:member=cell.ghost;break;
                case part::LayoutGroup::Halo:case part::LayoutGroup::Receive:member=cell.halo;break;
                case part::LayoutGroup::Send:member=cell.send;break;
                case part::LayoutGroup::ReceiveDependent:member=cell.receiveDependent;break;
                case part::LayoutGroup::Interior:member=!cell.ghost && cell.type==geo::CellType::Interior;break;
                case part::LayoutGroup::Boundary:member=!cell.ghost && cell.type==geo::CellType::Boundary;break;
                case part::LayoutGroup::FineInterface:member=!cell.ghost && cell.role!=geo::Auxiliary && (cell.role&geo::FineInterface);break;
                case part::LayoutGroup::CoarseInterface:member=!cell.ghost && cell.role!=geo::Auxiliary && (cell.role&geo::CoarseInterface);break;
            }
            if(member) expected.push_back(c);
        }
        requirePartition(values==expected,"Range/point selection lost or added members");
    }
    const auto& boundary=grid.boundaries(); constexpr auto qCount=lat::Model<Type>::nStencils-1;
    std::size_t expectedRecords=0;
    for(auto global:boundary.cells) expectedRecords+=layout.localIndex(global)!=part::invalidCellIndex;
    requirePartition(expectedRecords==layout.boundaries.cells.size(),"Visible Face/Q records were omitted");
    for(std::size_t face=0;face<boundary.nBoundaries;++face)
        for(auto row=layout.boundaryFaces[face].offset;row<layout.boundaryFaces[face].offset+layout.boundaryFaces[face].length;++row)
        {
            const auto global=layout.cells[layout.boundaries.cells[row]].global;
            bool found=false;
            for(auto entry=boundary.offsets[face];entry<boundary.offsets[face+1];++entry) if(boundary.cells[entry]==global)
            {
                found=true;
                requirePartition(std::equal(boundary.q.begin()+entry*qCount,boundary.q.begin()+(entry+1)*qCount,
                                           layout.boundaries.q.begin()+row*qCount),"Face/Q record changed during remapping");
            }
            requirePartition(found,"Boundary record has the wrong Face/cell");
        }
    for(const auto& node:layout.treeCells)
        requirePartition(node.realCell==layout.localIndex(grid.nodes()[node.globalNode].cell),"Sparse tree view points to an auxiliary slot");
    for(const auto& block:layout.interfaces)
    {
        const auto coarse=layout.cells[block.coarse].global;
        const auto match=std::find_if(grid.interfaces().blocks.begin(),grid.interfaces().blocks.end(),[&](const auto& b) {return b.coarse==coarse;});
        requirePartition(match!=grid.interfaces().blocks.end() && block.auxiliary.size()==(std::size_t{1}<<D),"Local interface lost its family");
        for(std::size_t child=0;child<block.auxiliary.size();++child)
            requirePartition(layout.cells[block.auxiliary[child]].global==match->auxiliary.begin+child,"Interface child mapping changed");
        requirePartition(block.fine.size()==match->fine.size(),"Interface fine list changed");
        for(std::size_t i=0;i<block.fine.size();++i)
            requirePartition(layout.cells[block.fine[i]].global==grid.interfaces().fineCells[match->fine.begin+i],"Interface fine mapping changed");
    }
}

template<std::size_t D,lat::LatticeType Type>
std::vector<unsigned char> signature(const part::LocalLayout<D,Type>& layout)
{
    part::detail::LayoutWriter w;w.word(layout.rank);w.word(layout.candidates);w.word(layout.selectedCandidate);
    w.word(layout.score.communicationRanges);w.word(layout.score.interfaceRanges);w.word(layout.score.boundaryRanges);w.word(layout.score.anchorDistance);
    for(auto x:layout.frame.origin) w.real(x);
    w.real(layout.frame.coarseSpacing);w.word(layout.frame.levelCount);
    w.word(layout.cells.size());
    for(std::size_t i=0;i<layout.cells.size();++i)
    {
        const auto& c=layout.cells[i];w.id(c.global);w.id(c.globalSourceNode);w.id(c.localSourceNode);w.word(c.level);w.word(c.owner);
        w.word(c.role);w.word(static_cast<unsigned>(c.type));w.word(c.owned);w.word(c.ghost);w.word(c.halo);w.word(c.send);w.word(c.receiveDependent);
        for(auto x:layout.coordinates[i]) w.integer(x);
    }
    w.ids(layout.lookupGlobal);w.ids(layout.lookupLocal);w.ids(layout.sourceIDs);
    for(const auto& row:layout.anchors) for(auto id:row) w.id(id);
    w.word(layout.peers.size());for(const auto& peer:layout.peers) {w.word(peer.rank);w.ids(peer.send);w.ids(peer.receive);}
    w.ids(layout.boundaries.offsets);w.ids(layout.boundaries.cells);w.word(layout.boundaries.q.size());for(auto q:layout.boundaries.q) w.real(q);
    w.word(layout.interfaces.size());for(const auto& b:layout.interfaces) {w.id(b.coarse);w.ids(b.auxiliary);w.ids(b.fine);}
    w.word(layout.intraRankExchanges.size());for(const auto& e:layout.intraRankExchanges)
    {w.id(e.a);w.id(e.b);w.word(e.dependencies.size());for(const auto& pair:e.dependencies) {w.id(pair.first);w.id(pair.second);}}
    auto spans=[&](const auto& values) {w.word(values.size());for(const auto& span:values) {w.id(span.offset);w.id(span.length);}};
    spans(layout.levelRanges);spans(layout.boundaryFaces);spans(layout.selectionRanges);w.ids(layout.selectionPoints);
    w.word(layout.groups.size());for(const auto& g:layout.groups)
    {w.word(g.level);w.word(static_cast<unsigned>(g.kind));w.id(g.ranges.offset);w.id(g.ranges.length);w.id(g.points.offset);w.id(g.points.length);}
    w.word(layout.nodes.size());for(const auto& n:layout.nodes)
    {
        w.id(n.globalNode);w.id(n.parentGlobal);w.id(n.firstGlobalChild);w.id(n.globalRealCell);w.word(n.level);
        for(auto x:n.coordinate) w.word(x);
        w.word(n.active);w.word(n.leaf);
        for(auto child:n.globalChildren) w.id(child);
        w.id(n.parent);w.id(n.realCell);for(auto child:n.children) w.id(child);
    }
    w.word(layout.treeCells.size());for(const auto& n:layout.treeCells) {w.id(n.globalNode);w.id(n.realCell);}
    return std::move(w.bytes);
}

template<std::size_t D,lat::LatticeType Type,typename MakeGrid>
void distributedFixture(MPI_Comm comm,int root,MakeGrid makeGrid,bool expectSharedGhost=false)
{
    int rank=0,size=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    using Grid=geo::MultilevelGrid<D,Type>;using Links=geo::NeighborLinks<D,Type>;using Builder=part::LocalLayoutBuilder<D,Type>;
    std::unique_ptr<Grid> grid;std::optional<Links> links;std::optional<part::PartitionPlan> plan;std::unique_ptr<Builder> builder;
    std::string rootError;
    if(rank==root)
    {
        try
        {
            grid=std::make_unique<Grid>(makeGrid());links=geo::linkNeighbors(*grid);
            part::PartitionOptions options;options.ranks=size;plan=part::partition(*grid,options);
            builder=std::make_unique<Builder>(*grid,*links,*plan);
        }
        catch(const std::exception& e) {rootError=e.what();}
    }
    const auto input=part::distributeLayoutInput<D,Type>(builder.get(),comm,root,rootError);
    part::detail::LayoutError error;
    try
    {
        const auto packet=part::encodeLayoutInput(input);
        requirePartition(packet==part::encodeLayoutInput(part::decodeLayoutInput<D,Type>(packet)),"Rank packet round-trip changed fields");
        auto damaged=packet;damaged.at(23)^=1;bool rejected=false;
        try {(void)part::decodeLayoutInput<D,Type>(damaged);}catch(const std::exception&) {rejected=true;}
        requirePartition(rejected,"Wrong packet dimension was accepted");
        damaged=packet;damaged.pop_back();rejected=false;
        try {(void)part::decodeLayoutInput<D,Type>(damaged);}catch(const std::exception&) {rejected=true;}
        requirePartition(rejected,"Truncated packet was accepted");
    }
    catch(const std::exception& e) {error=e.what();}
    part::detail::layoutErrorGate(comm,error);
    omp_set_num_threads(4);const auto layout=part::initializeLocalLayout(input,comm);
    try
    {
        omp_set_num_threads(1);const auto one=part::buildLocalLayout(input);
        omp_set_num_threads(4);const auto repeat=part::buildLocalLayout(input);
        requirePartition(signature(one)==signature(layout) && signature(repeat)==signature(layout),"Independent local build depends on thread count/history");
        part::LocalLayoutOptions baseline;baseline.tryMorton=false;baseline.tryTypeFirst=false;
        requirePartition(!(part::buildLocalLayout(input,baseline).score<layout.score),"Local sorting discarded a better baseline");
        for(const auto& cell:layout.cells) if(!cell.ghost)
            requirePartition(cell.localSourceNode<layout.nodes.size() && layout.nodes[cell.localSourceNode].globalNode==cell.globalSourceNode,
                             "Cell source node is not rank-local");
        for(std::size_t i=0;i<layout.nodes.size();++i)
        {
            const auto& n=layout.nodes[i];
            if(n.parent!=part::invalidCellIndex) requirePartition(layout.nodes[n.parent].globalNode==n.parentGlobal,"Compact parent link is wrong");
            for(std::size_t child=0;child<(std::size_t{1}<<D);++child)
                if(n.children[child]!=part::invalidCellIndex) requirePartition(layout.nodes[n.children[child]].globalNode==n.globalChildren[child],"Compact child link is wrong");
        }
        if(size==1) {requirePartition(layout.peers.empty(),"Single rank has peers");for(const auto& c:layout.cells) requirePartition(!c.halo,"Single rank has Halo");}
    }
    catch(const std::exception& e) {error=e.what();}
    part::detail::layoutErrorGate(comm,error);
    int shared=0;for(const auto& cell:layout.cells) shared=shared || (cell.ghost && cell.halo);
    MPI_Allreduce(MPI_IN_PLACE,&shared,1,MPI_INT,MPI_MAX,comm);
    if(expectSharedGhost && size>1 && !shared) error="Shared cross-rank Ghost fixture was not exercised.";
    part::detail::layoutErrorGate(comm,error);
    // Only the test root constructs serial references, never the production transport.
    for(int destination=0;destination<size;++destination)
    {
        std::vector<unsigned char> expected;
        if(rank==root)
        {
            try {const auto reference=builder->build(destination);checkLayout(*grid,*links,*plan,reference);expected=signature(reference);}
            catch(const std::exception& e) {error=e.what();}
        }
        part::detail::layoutErrorGate(comm,error);
        std::uint64_t bytes=expected.size();MPI_Bcast(&bytes,1,MPI_UINT64_T,root,comm);
        if(rank!=root) expected.resize(static_cast<std::size_t>(bytes));
        requirePartition(bytes<=static_cast<std::uint64_t>(std::numeric_limits<int>::max()),"Test reference packet is unexpectedly large");
        MPI_Bcast(expected.data(),static_cast<int>(bytes),MPI_BYTE,root,comm);
        if(rank==destination && expected!=signature(layout)) error="Distributed local layout differs from exact serial reference.";
        part::detail::layoutErrorGate(comm,error);
    }
}

template<lat::LatticeType Type>
void threeDimensions(MPI_Comm comm)
{
    distributedFixture<3,Type>(comm,0,[] {
        const geo::Box physical({0,0,0},{4,4,4});const geo::Sphere selection({2,2,2},0.4);
        geo::MultilevelOptions options;options.maxLevel=1;geo::MultilevelGrid<3,Type> grid(physical.bounds,1,options);
        grid.generate(physical,{geo::RefinementRule<3>::volume(geo::volumeRegion(selection),1)});return grid;
    });
}

} // namespace

int main(int argc,char** argv)
{
    int provided=0;if(MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided)!=MPI_SUCCESS) return 1;
    int mainThread=0;MPI_Is_thread_main(&mainThread);
    if(provided<MPI_THREAD_FUNNELED || !mainThread) MPI_Abort(MPI_COMM_WORLD,1);
    int rank=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);int result=0;
    try
    {
        omp_set_dynamic(0);omp_set_num_threads(4);radixOrder();
        auto square=[] {const geo::Rectangle shape({0,0},{4,4});geo::MultilevelGrid<2> grid(shape.bounds,1);grid.generate(shape);return grid;};
        distributedFixture<2,lat::LatticeType::D2Q9>(MPI_COMM_WORLD,0,square,true);
        threeDimensions<lat::LatticeType::D3Q15>(MPI_COMM_WORLD);threeDimensions<lat::LatticeType::D3Q19>(MPI_COMM_WORLD);threeDimensions<lat::LatticeType::D3Q27>(MPI_COMM_WORLD);
        distributedFixture<2,lat::LatticeType::D2Q9>(MPI_COMM_WORLD,0,[] {
            const geo::Rectangle shape({0,0},{4,4});const geo::Circle absent({20,20},1);geo::MultilevelGrid<2> grid(shape.bounds,1);grid.generate(absent);return grid;
        });
        bool rejected=false;
        try {(void)part::distributeLayoutInput<2,lat::LatticeType::D2Q9>(nullptr,MPI_COMM_WORLD,0,rank==0 ? "Intentional root preparation failure" : "");}
        catch(const std::exception&) {rejected=true;}
        int rejectedCount=rejected ? 1 : 0,size=0;MPI_Comm_size(MPI_COMM_WORLD,&size);MPI_Allreduce(MPI_IN_PLACE,&rejectedCount,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        requirePartition(rejectedCount==size,"Root preparation failure did not reach every rank");
        MPI_Comm split=MPI_COMM_NULL;MPI_Comm_split(MPI_COMM_WORLD,rank%2,-rank,&split);int splitSize=0;MPI_Comm_size(split,&splitSize);
        distributedFixture<2,lat::LatticeType::D2Q9>(split,splitSize-1,square);MPI_Comm_free(&split);
        if(rank==0) std::cout << "[PASS] local_layout_anchors_and_remapping\n[PASS] local_layout_ranges_and_radix\n";
    }
    catch(const std::exception& e) {std::cerr << "[FAIL] rank " << rank << " layout: " << e.what() << '\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
    MPI_Allreduce(MPI_IN_PLACE,&result,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);MPI_Finalize();return result;
}
