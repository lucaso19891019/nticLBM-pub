#include "partition.hpp"
#include <metis.h>

#include <cmath>
#include <map>
#include <mutex>
#include <numeric>
#include <string>
#include <tuple>

namespace ntic::lbm::mpi
{
namespace
{
idx_t checkedIndex(std::size_t value)
{
    if(value>static_cast<std::size_t>(std::numeric_limits<idx_t>::max())) throw std::length_error("METIS index capacity exceeded.");
    return static_cast<idx_t>(value);
}

std::vector<idx_t> integerWeights(const std::vector<double>& values)
{
    double largest=0; long double sum=0;
    for(auto value:values)
    {
        if(!(value>0) || !std::isfinite(value)) throw std::invalid_argument("Invalid METIS graph weight.");
        largest=std::max(largest,value); sum+=value;
    }
    if(values.empty()) return {};
    // Reserve headroom for METIS's integer sums and refinement arithmetic.
    const long double budget=std::numeric_limits<idx_t>::max()/64;
    if(values.size()>budget/2) throw std::length_error("Graph exceeds METIS weight capacity.");
    const long double scale=std::min(1.0e6L/largest,(budget-values.size())/sum);
    std::vector<idx_t> result; result.reserve(values.size());
    for(auto value:values) result.push_back(std::max<idx_t>(1,static_cast<idx_t>(std::floor(value*scale+0.5L))));
    return result;
}

std::vector<PartitionBlock> mergeBlocks(std::vector<PartitionBlock> blocks,unsigned dimension,int passes)
{
    // Merge only equal cross sections with a complete shared face. Every merge
    // preserves an exact rectangular tiling, including empty spatial boxes.
    for(int pass=0; pass<passes; ++pass)
    {
        const auto initial=blocks.size();
        for(unsigned axis=0; axis<dimension; ++axis)
        {
            const auto sameSection=[&](const auto& a,const auto& b) {
                if(a.owner!=b.owner) return false;
                for(unsigned d=0; d<dimension; ++d)
                    if(d!=axis && (a.lower[d]!=b.lower[d] || a.upper[d]!=b.upper[d])) return false;
                return true;
            };
            std::sort(blocks.begin(),blocks.end(),[&](const auto& a,const auto& b) {
                if(a.owner!=b.owner) return a.owner<b.owner;
                for(unsigned d=0; d<dimension; ++d) if(d!=axis)
                {
                    if(a.lower[d]!=b.lower[d]) return a.lower[d]<b.lower[d];
                    if(a.upper[d]!=b.upper[d]) return a.upper[d]<b.upper[d];
                }
                return a.lower[axis]<b.lower[axis];
            });
            std::size_t count=0;
            for(std::size_t i=0; i<blocks.size(); ++i)
            {
                if(count && sameSection(blocks[count-1],blocks[i]) && blocks[count-1].upper[axis]==blocks[i].lower[axis])
                { blocks[count-1].upper[axis]=blocks[i].upper[axis]; blocks[count-1].work+=blocks[i].work; }
                else { if(count!=i) blocks[count]=blocks[i]; ++count; }
            }
            blocks.resize(count);
        }
        if(blocks.size()==initial) break;
    }
    std::sort(blocks.begin(),blocks.end(),[](const auto& a,const auto& b) {
        return std::tie(a.owner,a.lower,a.upper)<std::tie(b.owner,b.lower,b.upper);
    });
    return blocks;
}

void fillEmptyOwners(const PartitionGraph& graph,std::vector<int>& owners)
{
    // Empty topology boxes carry neither graph weight nor communication.
    // Assign them in lexicographic spatial order to complete the tiling; this
    // never changes METIS ownership of a nonempty block.
    std::vector<std::size_t> order(owners.size()); std::iota(order.begin(),order.end(),0);
    std::sort(order.begin(),order.end(),[&](auto a,auto b) { return graph.blocks[a].lower<graph.blocks[b].lower; });
    int previous=-1;
    for(auto b:order) { if(owners[b]>=0) previous=owners[b]; else if(previous>=0) owners[b]=previous; }
    int next=0;
    for(auto it=order.rbegin(); it!=order.rend(); ++it)
    { if(owners[*it]>=0) next=owners[*it]; else owners[*it]=next; }
}

PartitionStatistics evaluate(const PartitionGraph& graph,const std::vector<int>& owner,const PartitionOptions& options,
                             const std::vector<PartitionBlock>& merged)
{
    PartitionStatistics stats;
    stats.rankWeights.assign(options.ranks,0); stats.computeTimes.assign(options.ranks,0);
    stats.communicationTimes.assign(options.ranks,0); stats.blocksPerRank.assign(options.ranks,0);
    for(std::size_t b=0; b<graph.blocks.size(); ++b)
    {
        const double work=graph.blocks[b].work;
        stats.rankWeights[owner[b]]+=work; stats.totalWeight+=work;
        if(work>0) ++stats.graphVertices;
        stats.largestAtomicWork=std::max(stats.largestAtomicWork,work);
    }
    const auto& hardware=options.hardware;
    double totalRate=0;
    for(int rank=0; rank<options.ranks; ++rank)
    {
        const double rate=hardware.computeRates.empty() ? 1 : hardware.computeRates[rank];
        totalRate+=rate; stats.computeTimes[rank]=stats.rankWeights[rank]/rate;
        stats.maxComputeTime=std::max(stats.maxComputeTime,stats.computeTimes[rank]);
    }
    if(!std::isfinite(stats.totalWeight) || !std::isfinite(totalRate) || !std::isfinite(stats.maxComputeTime))
        throw std::overflow_error("Partition load/time overflow.");
    if(stats.totalWeight>0) stats.maxToMean=stats.maxComputeTime/(stats.totalWeight/totalRate);
    struct Traffic { double bytes=0; std::uint64_t phases=0; };
    std::map<std::pair<int,int>,Traffic> traffic;
    for(const auto& edge:graph.edges)
    {
        ++stats.graphEdges;
        const int a=owner[edge.a],b=owner[edge.b];
        if(a==b) continue;
        stats.cutWeight+=edge.weight;
        auto& item=traffic[std::minmax(a,b)]; item.bytes+=edge.bytes; item.phases|=edge.phases;
    }
    for(const auto& entry:traffic)
    {
        const auto a=entry.first.first,b=entry.first.second;
        const bool same=hardware.nodeIds.empty() || hardware.nodeIds[a]==hardware.nodeIds[b];
        const auto latency=same ? hardware.intraNodeLatency : hardware.interNodeLatency;
        const auto bandwidth=same ? hardware.intraNodeBandwidth : hardware.interNodeBandwidth;
        double messages=0;
        for(unsigned l=0; l<31; ++l)
        {
            const double frequency=std::ldexp(1.0,static_cast<int>(l));
            if(entry.second.phases&(std::uint64_t{1}<<l)) messages+=frequency;
            if(entry.second.phases&(std::uint64_t{1}<<(l+32))) messages+=frequency;
        }
        const double time=messages*latency+entry.second.bytes/bandwidth;
        stats.communicationTimes[a]+=time; stats.communicationTimes[b]+=time;
    }
    for(const auto& block:merged) ++stats.blocksPerRank[block.owner];
    for(int rank=0; rank<options.ranks; ++rank)
    {
        const double time=stats.computeTimes[rank]+stats.communicationTimes[rank];
        if(!std::isfinite(time)) throw std::overflow_error("Communication time overflow.");
        stats.maxEstimatedTime=std::max(stats.maxEstimatedTime,time);
    }
    return stats;
}
} // namespace

PartitionPlan partitionGraph(const PartitionGraph& graph,const PartitionOptions& options)
{
    detail::validateOptions(options);
    if(options.ranks<1 || graph.dimension<2 || graph.dimension>3 || options.candidates<1)
        throw std::invalid_argument("Invalid block graph configuration.");
    std::vector<std::size_t> vertices,map(graph.blocks.size(),invalidCellIndex);
    std::vector<double> work;
    for(std::size_t b=0; b<graph.blocks.size(); ++b)
    {
        if(graph.blocks[b].work<0 || !std::isfinite(graph.blocks[b].work)) throw std::invalid_argument("Invalid block work.");
        if(graph.blocks[b].work>0) { map[b]=vertices.size(); vertices.push_back(b); work.push_back(graph.blocks[b].work); }
    }
    std::vector<int> ranks(options.ranks); std::iota(ranks.begin(),ranks.end(),0);
    const auto rate=[&](int rank) { return options.hardware.computeRates.empty() ? 1.0 : options.hardware.computeRates.at(rank); };
    // If there are fewer indivisible work blocks than ranks, use the fastest
    // ranks and report the unavoidable empty ranks rather than splitting cells.
    std::stable_sort(ranks.begin(),ranks.end(),[&](int a,int b) { return rate(a)>rate(b); });
    ranks.resize(std::min(ranks.size(),vertices.size()));
    std::sort(ranks.begin(),ranks.end());
    const auto n=vertices.size();
    std::vector<idx_t> xadj(n+1,0);
    for(const auto& edge:graph.edges)
    {
        if(edge.a>=map.size() || edge.b>=map.size() || edge.a==edge.b || map[edge.a]==invalidCellIndex || map[edge.b]==invalidCellIndex)
            throw std::invalid_argument("Invalid block-graph edge.");
        ++xadj[map[edge.a]+1]; ++xadj[map[edge.b]+1];
    }
    for(std::size_t i=1; i<xadj.size(); ++i)
    {
        if(xadj[i]>std::numeric_limits<idx_t>::max()-xadj[i-1]) throw std::length_error("METIS adjacency overflow.");
        xadj[i]+=xadj[i-1];
    }
    const auto entries=static_cast<std::size_t>(xadj.back());
    std::vector<idx_t> adjacency(std::max<std::size_t>(1,entries),0); // Non-null for isolated graphs.
    std::vector<double> edgeValues(entries); auto cursor=xadj;
    for(const auto& edge:graph.edges)
    {
        const auto a=map[edge.a],b=map[edge.b];
        const auto i=static_cast<std::size_t>(cursor[a]++),j=static_cast<std::size_t>(cursor[b]++);
        adjacency[i]=checkedIndex(b); adjacency[j]=checkedIndex(a); edgeValues[i]=edgeValues[j]=edge.weight;
    }
    auto vertexWeights=integerWeights(work),edgeWeights=integerWeights(edgeValues);
    std::vector<real_t> targets;
    double totalRate=0; for(auto rank:ranks) totalRate+=rate(rank);
    if(!ranks.empty() && (!(totalRate>0) || !std::isfinite(totalRate))) throw std::invalid_argument("Invalid rank capacities.");
    for(auto rank:ranks)
    {
        const auto fraction=rate(rank)/totalRate;
        if(!(fraction>0)) throw std::invalid_argument("Rank capacity ratio is not representable.");
        targets.push_back(fraction);
    }
    struct Candidate { std::vector<int> owners; std::vector<PartitionBlock> boxes; PartitionStatistics stats; };
    std::vector<Candidate> candidates;
    const int trials=ranks.size()>1 ? options.candidates : 1;
    for(int trial=0; trial<trials; ++trial)
    {
        Candidate candidate; candidate.owners.assign(graph.blocks.size(),-1);
        std::vector<idx_t> assignment(n,0);
        if(ranks.size()>1)
        {
            idx_t count=checkedIndex(n),constraints=1,parts=checkedIndex(ranks.size()),objective=0;
            idx_t settings[METIS_NOPTIONS]; METIS_SetDefaultOptions(settings);
            settings[METIS_OPTION_SEED]=options.seed+trial;
            settings[METIS_OPTION_NUMBERING]=0; settings[METIS_OPTION_CONTIG]=0;
            settings[METIS_OPTION_OBJTYPE]=METIS_OBJTYPE_CUT;
            settings[METIS_OPTION_MINCONN]=1; settings[METIS_OPTION_NCUTS]=1;
            real_t tolerance=options.imbalanceTolerance;
            // METIS owns mutable allocator/RNG state. Serialize invocations,
            // while graph construction remains parallel in the caller.
            static std::mutex metisMutex;
            std::lock_guard<std::mutex> lock(metisMutex);
            const auto status=METIS_PartGraphKway(&count,&constraints,xadj.data(),adjacency.data(),vertexWeights.data(),
                nullptr,edgeWeights.empty() ? nullptr : edgeWeights.data(),&parts,targets.data(),&tolerance,settings,&objective,assignment.data());
            if(status!=METIS_OK) throw std::runtime_error("METIS_PartGraphKway failed with status "+std::to_string(status));
        }
        for(std::size_t v=0; v<n; ++v)
        {
            if(assignment[v]<0 || static_cast<std::size_t>(assignment[v])>=ranks.size()) throw std::runtime_error("METIS returned an invalid partition.");
            candidate.owners[vertices[v]]=ranks[static_cast<std::size_t>(assignment[v])];
        }
        fillEmptyOwners(graph,candidate.owners);
        auto boxes=graph.blocks;
        for(std::size_t b=0; b<boxes.size(); ++b) boxes[b].owner=candidate.owners[b];
        candidate.boxes=mergeBlocks(std::move(boxes),graph.dimension,options.mergePasses);
        candidate.stats=evaluate(graph,candidate.owners,options,candidate.boxes);
        candidates.push_back(std::move(candidate));
    }
    double bestCompute=std::numeric_limits<double>::infinity();
    for(const auto& candidate:candidates) bestCompute=std::min(bestCompute,candidate.stats.maxComputeTime);
    std::size_t best=0; double bestSecondary=std::numeric_limits<double>::infinity();
    std::size_t bestBlocks=std::numeric_limits<std::size_t>::max();
    for(std::size_t c=0; c<candidates.size(); ++c)
    {
        const auto& candidate=candidates[c];
        if(candidate.stats.maxComputeTime>bestCompute*(1+options.loadSlack)) continue;
        const auto maxBlocks=*std::max_element(candidate.stats.blocksPerRank.begin(),candidate.stats.blocksPerRank.end());
        const double communication=*std::max_element(candidate.stats.communicationTimes.begin(),candidate.stats.communicationTimes.end());
        const double secondary=communication+options.blockPenalty*maxBlocks;
        if(secondary<bestSecondary || (secondary==bestSecondary && candidate.boxes.size()<bestBlocks))
        { best=c; bestSecondary=secondary; bestBlocks=candidate.boxes.size(); }
    }
    PartitionPlan result;
    result.owner=std::move(candidates[best].owners); result.blocks=std::move(candidates[best].boxes);
    result.statistics=std::move(candidates[best].stats);
    return result;
}
} // namespace ntic::lbm::mpi
