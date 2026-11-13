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

std::vector<BlockExchange> blockExchanges(const PartitionGraph& graph,const std::vector<std::size_t>& atomToBlock)
{
    std::map<std::pair<std::size_t,std::size_t>,BlockExchange> exchanges;
    for(const auto& edge:graph.edges)
    {
        const auto a=atomToBlock[edge.a],b=atomToBlock[edge.b];
        if(a==b) continue;
        const auto key=std::minmax(a,b); auto& item=exchanges[key];
        item.a=key.first; item.b=key.second; item.bytes+=edge.bytes; item.phases|=edge.phases;
    }
    std::vector<BlockExchange> result; result.reserve(exchanges.size());
    for(const auto& item:exchanges) result.push_back(item.second);
    return result;
}

double exchangeTime(double bytes,std::uint64_t phases,double latency,double bandwidth)
{
    double messages=0;
    for(unsigned l=0; l<31; ++l)
    {
        const double frequency=std::ldexp(1.0,static_cast<int>(l));
        if(phases&(std::uint64_t{1}<<l)) messages+=frequency;
        if(phases&(std::uint64_t{1}<<(l+32))) messages+=frequency;
    }
    return messages*latency+bytes/bandwidth;
}

PartitionStatistics evaluate(const PartitionGraph& graph,const std::vector<int>& owner,const PartitionOptions& options,
                             const std::vector<PartitionBlock>& merged,const std::vector<BlockExchange>& exchanges)
{
    PartitionStatistics stats;
    stats.rankWeights.assign(options.ranks,0); stats.computeTimes.assign(options.ranks,0);
    stats.communicationTimes.assign(options.ranks,0); stats.blocksPerRank.assign(options.ranks,0);
    stats.intraRankCommunicationTimes.assign(options.ranks,0); stats.interRankCommunicationTimes.assign(options.ranks,0);
    for(std::size_t b=0; b<graph.blocks.size(); ++b)
    {
        const double work=graph.blocks[b].work;
        if(work==0) continue;
        if(owner[b]<0 || owner[b]>=options.ranks) throw std::logic_error("Unowned populated atom.");
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
        if(a!=b) stats.cutWeight+=edge.weight;
    }
    for(const auto& exchange:exchanges)
    {
        const int a=merged[exchange.a].owner,b=merged[exchange.b].owner;
        if(a==b)
        {
            stats.intraRankCommunicationTimes[a]+=exchangeTime(exchange.bytes,exchange.phases,
                                                              hardware.intraRankLatency,hardware.intraRankBandwidth);
        }
        else
        {
            auto& item=traffic[std::minmax(a,b)]; item.bytes+=exchange.bytes; item.phases|=exchange.phases;
        }
    }
    for(const auto& entry:traffic)
    {
        const auto a=entry.first.first,b=entry.first.second;
        const bool same=hardware.nodeIds.empty() || hardware.nodeIds[a]==hardware.nodeIds[b];
        const auto latency=same ? hardware.intraNodeLatency : hardware.interNodeLatency;
        const auto bandwidth=same ? hardware.intraNodeBandwidth : hardware.interNodeBandwidth;
        const double time=exchangeTime(entry.second.bytes,entry.second.phases,latency,bandwidth);
        stats.interRankCommunicationTimes[a]+=time; stats.interRankCommunicationTimes[b]+=time;
    }
    for(const auto& block:merged) ++stats.blocksPerRank[block.owner];
    for(int rank=0; rank<options.ranks; ++rank)
    {
        stats.communicationTimes[rank]=stats.intraRankCommunicationTimes[rank]+stats.interRankCommunicationTimes[rank];
        const double time=stats.computeTimes[rank]+stats.communicationTimes[rank];
        if(!std::isfinite(time)) throw std::overflow_error("Communication time overflow.");
        stats.maxEstimatedTime=std::max(stats.maxEstimatedTime,time);
    }
    return stats;
}
// Rectangular aggregation of the weighted METIS atom assignment.
struct CutNode
{
    unsigned axis=0;
    std::uint32_t plane=0;
    int left=-1,right=-1,rank=-1;
};

struct Candidate
{
    bool geometrySeed=false;
    std::vector<CutNode> tree;
    std::vector<int> labels;
    std::vector<int> owners;
    std::vector<PartitionBlock> boxes;
    std::vector<std::size_t> atomToBlock;
    std::vector<BlockExchange> exchanges;
    PartitionStatistics stats;
};

double rankRate(int rank,const PartitionOptions& options)
{
    return options.hardware.computeRates.empty() ? 1.0 : options.hardware.computeRates[rank];
}

// Only populated atoms constrain the cover. Cuts may span arbitrary Dry gaps,
// but cannot intersect an indivisible coarse/auxiliary family atom.
bool route(const PartitionGraph& graph,Candidate& candidate,int ranks)
{
    candidate.owners.assign(graph.blocks.size(),-1);
    candidate.atomToBlock.assign(graph.blocks.size(),invalidCellIndex);
    candidate.boxes.assign(ranks,{});
    std::vector<bool> occupied(ranks,false);
    for(int rank=0;rank<ranks;++rank) candidate.boxes[rank].owner=rank;
    for(std::size_t atom=0;atom<graph.blocks.size();++atom)
    {
        const auto& input=graph.blocks[atom];
        if(input.work==0) continue;
        int node=0;
        while(candidate.tree[node].rank<0)
        {
            const auto& cut=candidate.tree[node];
            if(input.upper[cut.axis]<=cut.plane) node=cut.left;
            else if(input.lower[cut.axis]>=cut.plane) node=cut.right;
            else return false;
        }
        const int rank=candidate.tree[node].rank;
        auto& box=candidate.boxes[rank];
        if(!occupied[rank]) { box.lower=input.lower;box.upper=input.upper;occupied[rank]=true; }
        else for(unsigned d=0;d<graph.dimension;++d)
        {
            box.lower[d]=std::min(box.lower[d],input.lower[d]);
            box.upper[d]=std::max(box.upper[d],input.upper[d]);
        }
        box.work+=input.work;
        candidate.owners[atom]=rank;candidate.atomToBlock[atom]=static_cast<std::size_t>(rank);
    }
    return std::all_of(occupied.begin(),occupied.end(),[](bool value){return value;});
}

void score(const PartitionGraph& graph,Candidate& candidate,const PartitionOptions& options)
{
    candidate.exchanges=blockExchanges(graph,candidate.atomToBlock);
    candidate.stats=evaluate(graph,candidate.owners,options,candidate.boxes,candidate.exchanges);
}

struct CutChoice
{
    unsigned axis=0;
    std::uint32_t plane=0;
    std::size_t leftRanks=0;
    double balance=0,migration=0,shape=0;
    std::vector<int> order;
};

bool cutOrder(const CutChoice& a,const CutChoice& b)
{
    return std::tie(a.balance,a.migration,a.shape,a.axis,a.plane,a.leftRanks,a.order)<
           std::tie(b.balance,b.migration,b.shape,b.axis,b.plane,b.leftRanks,b.order);
}

std::vector<CutChoice> cutChoices(const PartitionGraph& graph,const std::vector<std::size_t>& atoms,
                                const std::vector<int>& ranks,const std::vector<int>& labels,
                                const PartitionOptions& options,unsigned variant,bool neutral)
{
    std::vector<CutChoice> choices;
    std::array<std::uint32_t,3> low{},high{};
    low=graph.blocks[atoms.front()].lower;high=graph.blocks[atoms.front()].upper;
    for(auto atom:atoms) for(unsigned d=0;d<graph.dimension;++d)
    {low[d]=std::min(low[d],graph.blocks[atom].lower[d]);high[d]=std::max(high[d],graph.blocks[atom].upper[d]);}
    double total=0;
    for(auto atom:atoms) total+=graph.blocks[atom].work;

    for(unsigned axis=0;axis<graph.dimension;++axis)
    {
        std::vector<double> centers(options.ranks,0),weights(options.ranks,0);
        for(auto atom:atoms)
        {
            const auto& box=graph.blocks[atom];const int rank=labels[atom];
            centers[rank]+=box.work*(0.5*box.lower[axis]+0.5*box.upper[axis]);weights[rank]+=box.work;
        }
        auto order=ranks;
        std::stable_sort(order.begin(),order.end(),[&](int a,int b)
        {
            const double ca=weights[a]>0 ? centers[a]/weights[a] : 0;
            const double cb=weights[b]>0 ? centers[b]/weights[b] : 0;
            return ca==cb ? a<b : ca<cb;
        });
        if(neutral) std::sort(order.begin(),order.end());
        if(variant%3==2) std::reverse(order.begin(),order.end());
        std::vector<std::size_t> sorted=atoms;
        std::sort(sorted.begin(),sorted.end(),[&](auto a,auto b)
        { return std::tie(graph.blocks[a].lower[axis],a)<std::tie(graph.blocks[b].lower[axis],b); });
        std::vector<std::uint32_t> suffix(sorted.size());
        suffix.back()=graph.blocks[sorted.back()].lower[axis];
        for(std::size_t i=sorted.size()-1;i>0;--i)
            suffix[i-1]=std::min(suffix[i],graph.blocks[sorted[i-1]].lower[axis]);
        std::vector<std::size_t> rankPosition(options.ranks,ranks.size());
        for(std::size_t i=0;i<order.size();++i) rankPosition[order[i]]=i;
        std::vector<double> prefixLabel(ranks.size(),0),suffixRate(ranks.size()+1,0);
        for(std::size_t k=ranks.size();k>0;--k) suffixRate[k-1]=suffixRate[k]+rankRate(order[k-1],options);
        double prefix=0;std::uint32_t upper=0;
        for(std::size_t i=0;i+1<sorted.size();++i)
        {
            const auto atom=sorted[i];const auto& box=graph.blocks[atom];
            prefix+=box.work;upper=std::max(upper,box.upper[axis]);
            const auto position=rankPosition[labels[atom]];
            if(position<ranks.size()) prefixLabel[position]+=box.work;
            if(upper>suffix[i+1]) continue;
            double leftRate=0,leftLabels=0,totalLeftLabels=0;
            for(std::size_t k=1;k<ranks.size();++k)
            {
                leftRate+=rankRate(order[k-1],options);
                leftLabels+=prefixLabel[k-1];totalLeftLabels+=weights[order[k-1]];
                if(i+1<k || sorted.size()-i-1<ranks.size()-k) continue;
                // Bound rank-group enumeration for very large process counts.
                if(ranks.size()>8 && k!=ranks.size()/2 && k!=(ranks.size()+1)/2) continue;
                CutChoice choice;choice.axis=axis;choice.plane=upper;choice.leftRanks=k;choice.order=order;
                choice.balance=std::max(prefix/leftRate,(total-prefix)/suffixRate[k]);
                choice.migration=neutral ? 0 : (prefix-leftLabels)+(totalLeftLabels-leftLabels);
                // Balanced rank counts favour shallower trees when load ties.
                choice.shape=0;
                for(bool left:{true,false})
                {
                    double shortest=std::numeric_limits<double>::infinity(),longest=0;
                    for(unsigned d=0;d<graph.dimension;++d)
                    {
                        const double length=d==axis ? (left ? upper-low[d] : high[d]-upper) : high[d]-low[d];
                        shortest=std::min(shortest,length);longest=std::max(longest,length);
                    }
                    choice.shape=std::max(choice.shape,longest/shortest);
                }
                auto place=std::lower_bound(choices.begin(),choices.end(),choice,cutOrder);
                if(place!=choices.end() || choices.size()<static_cast<std::size_t>(options.cutCandidates))
                {
                    choices.insert(place,std::move(choice));
                    if(choices.size()>static_cast<std::size_t>(options.cutCandidates)) choices.pop_back();
                }
            }
        }
    }
    std::sort(choices.begin(),choices.end(),cutOrder);
    // Select diverse legal planes, rather than retaining many copies of one cut.
    std::vector<CutChoice> selected;
    for(auto& choice:choices)
    {
        bool duplicate=false;
        for(const auto& item:selected)
            if(item.axis==choice.axis && item.plane==choice.plane && item.leftRanks==choice.leftRanks && item.order==choice.order)
                duplicate=true;
        if(!duplicate) selected.push_back(std::move(choice));
        if(selected.size()==static_cast<std::size_t>(options.cutCandidates)) break;
    }
    return selected;
}

bool buildTree(const PartitionGraph& graph,const std::vector<std::size_t>& atoms,const std::vector<int>& ranks,
               const std::vector<int>& labels,const PartitionOptions& options,unsigned variant,bool neutral,
               std::vector<CutNode>& tree,int& result,std::size_t& attempts)
{
    if(atoms.size()<ranks.size() || attempts==0) return false;
    --attempts;
    const auto original=tree.size();result=static_cast<int>(original);tree.push_back({});
    if(ranks.size()==1) { tree[result].rank=ranks.front();return true; }
    const auto choices=cutChoices(graph,atoms,ranks,labels,options,variant,neutral);
    if(choices.empty()) { tree.resize(original);return false; }
    // Limited backtracking keeps initialization predictable while admitting
    // different split axes/rank-group sizes when an early cut is infeasible.
    for(std::size_t attempt=0;attempt<choices.size();++attempt)
    {
        tree.resize(original+1);
        const auto& cut=choices[(attempt+variant)%choices.size()];
        std::vector<std::size_t> left,right;
        for(auto atom:atoms)
            (graph.blocks[atom].upper[cut.axis]<=cut.plane ? left : right).push_back(atom);
        std::vector<int> leftRanks(cut.order.begin(),cut.order.begin()+cut.leftRanks);
        std::vector<int> rightRanks(cut.order.begin()+cut.leftRanks,cut.order.end());
        int a=-1,b=-1;
        if(buildTree(graph,left,leftRanks,labels,options,variant,neutral,tree,a,attempts) &&
           buildTree(graph,right,rightRanks,labels,options,variant,neutral,tree,b,attempts))
        { tree[result]={cut.axis,cut.plane,a,b,-1};return true; }
    }
    tree.resize(original);return false;
}

void subtreeRanks(const std::vector<CutNode>& tree,int node,std::vector<int>& ranks)
{
    const auto& cut=tree[node];
    if(cut.rank>=0) ranks.push_back(cut.rank);
    else {subtreeRanks(tree,cut.left,ranks);subtreeRanks(tree,cut.right,ranks);}
}

void compactTree(const std::vector<CutNode>& input,int node,std::vector<CutNode>& output,int& result)
{
    result=static_cast<int>(output.size());output.push_back(input[node]);
    if(input[node].rank<0)
    {
        int left=-1,right=-1;
        compactTree(input,input[node].left,output,left);compactTree(input,input[node].right,output,right);
        output[result].left=left;output[result].right=right;
    }
}

std::vector<std::size_t> atomsInSubtree(const Candidate& candidate,int node,const std::vector<std::size_t>& vertices,int ranks)
{
    std::vector<int> selected;subtreeRanks(candidate.tree,node,selected);
    std::vector<bool> member(ranks,false);for(auto rank:selected) member[rank]=true;
    std::vector<std::size_t> atoms;
    for(auto atom:vertices) if(member[candidate.owners[atom]]) atoms.push_back(atom);
    return atoms;
}

double communication(const Candidate& candidate)
{
    return *std::max_element(candidate.stats.communicationTimes.begin(),candidate.stats.communicationTimes.end());
}

bool candidateOrder(const Candidate& a,const Candidate& b,double bestCompute,const PartitionOptions& options)
{
    const bool eligibleA=a.stats.maxComputeTime<=bestCompute*(1+options.loadSlack);
    const bool eligibleB=b.stats.maxComputeTime<=bestCompute*(1+options.loadSlack);
    if(eligibleA!=eligibleB) return eligibleA;
    if(eligibleA && communication(a)!=communication(b)) return communication(a)<communication(b);
    if(a.stats.maxComputeTime!=b.stats.maxComputeTime) return a.stats.maxComputeTime<b.stats.maxComputeTime;
    return a.owners<b.owners;
}

// Recut an entire multi-rank subtree. Atom labels remain the METIS proposal;
// a candidate is admitted only after complete geometric routing and exact
// full-graph communication evaluation, never from its local cut score alone.
bool rebuild(const PartitionGraph& graph,const Candidate& source,int node,const std::vector<std::size_t>& vertices,
             const std::vector<int>& labels,const PartitionOptions& options,unsigned variant,Candidate& target)
{
    target.tree=source.tree;target.labels=source.labels;target.geometrySeed=source.geometrySeed;
    auto atoms=atomsInSubtree(source,node,vertices,options.ranks);
    std::vector<int> ranks;subtreeRanks(source.tree,node,ranks);
    int replacement=-1;
    std::size_t attempts=static_cast<std::size_t>(options.ranks)*options.cutCandidates*8;
    if(!buildTree(graph,atoms,ranks,labels,options,variant,source.geometrySeed,target.tree,replacement,attempts)) return false;
    target.tree[node]=target.tree[replacement];
    std::vector<CutNode> compact;int root=-1;compactTree(target.tree,0,compact,root);target.tree=std::move(compact);
    if(!route(graph,target,options.ranks)) return false;
    score(graph,target,options);return true;
}

// Neighbouring legal positions for an existing cut, preserving its child
// rank sets. Two independent cuts move simultaneously before evaluation.
std::vector<std::uint32_t> adjacentPlanes(const PartitionGraph& graph,const Candidate& source,int node,
                                        const std::vector<std::size_t>& vertices,int ranks)
{
    auto atoms=atomsInSubtree(source,node,vertices,ranks);const auto& cut=source.tree[node];
    std::sort(atoms.begin(),atoms.end(),[&](auto a,auto b)
    { return std::tie(graph.blocks[a].lower[cut.axis],a)<std::tie(graph.blocks[b].lower[cut.axis],b); });
    std::vector<std::uint32_t> suffix(atoms.size());
    suffix.back()=graph.blocks[atoms.back()].lower[cut.axis];
    for(std::size_t i=atoms.size()-1;i>0;--i)
        suffix[i-1]=std::min(suffix[i],graph.blocks[atoms[i-1]].lower[cut.axis]);
    std::uint32_t upper=0,before=cut.plane,after=cut.plane;
    for(std::size_t i=0;i+1<atoms.size();++i)
    {
        upper=std::max(upper,graph.blocks[atoms[i]].upper[cut.axis]);
        if(upper>suffix[i+1]) continue;
        if(upper<cut.plane) before=upper;
        else if(upper>cut.plane && after==cut.plane) after=upper;
    }
    std::vector<std::uint32_t> result;
    if(before!=cut.plane) result.push_back(before);
    if(after!=cut.plane) result.push_back(after);
    return result;
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
    for(const auto& edge:graph.edges)
    {
        if(edge.a>=graph.blocks.size() || edge.b>=graph.blocks.size() || edge.a==edge.b ||
           graph.blocks[edge.a].work==0 || graph.blocks[edge.b].work==0)
            throw std::invalid_argument("Graph edges must connect populated distinct atoms.");
        detail::positiveWeight(edge.weight);detail::nonnegative(edge.bytes);
        if(edge.phases&((std::uint64_t{1}<<31)|(std::uint64_t{1}<<63)))
            throw std::invalid_argument("Graph edge uses an unsupported update phase.");
    }
    if(vertices.empty())
    {
        PartitionPlan empty;empty.dimension=graph.dimension;empty.owner.assign(graph.blocks.size(),-1);
        empty.cellBlocks.assign(graph.blocks.size(),invalidCellIndex);
        empty.statistics=evaluate(graph,empty.owner,options,{},{});return empty;
    }
    if(vertices.size()<ranks.size())
        throw std::invalid_argument("Fewer populated indivisible atoms than requested ranks.");
    for(auto atom:vertices)
        for(unsigned d=0;d<graph.dimension;++d)
            if(graph.blocks[atom].upper[d]<=graph.blocks[atom].lower[d])
                throw std::invalid_argument("Degenerate populated atom.");
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
    std::vector<Candidate> beam;
    std::size_t evaluated=0,jointCandidates=0;
    double bestCompute=std::numeric_limits<double>::infinity();
    const auto retain=[&](Candidate candidate)
    {
        ++evaluated;bestCompute=std::min(bestCompute,candidate.stats.maxComputeTime);
        for(const auto& item:beam) if(item.owners==candidate.owners) return false;
        beam.push_back(std::move(candidate));
        const auto lowest=std::min_element(beam.begin(),beam.end(),[](const auto& a,const auto& b)
        {return a.stats.maxComputeTime<b.stats.maxComputeTime;});
        const auto loadOwners=lowest->owners;
        std::sort(beam.begin(),beam.end(),[&](const auto& a,const auto& b)
        {return candidateOrder(a,b,bestCompute,options);});
        if(beam.size()>static_cast<std::size_t>(options.searchWidth))
        {
            auto load=std::find_if(beam.begin(),beam.end(),[&](const auto& item){return item.owners==loadOwners;});
            if(static_cast<std::size_t>(load-beam.begin())>=static_cast<std::size_t>(options.searchWidth))
                std::swap(beam[options.searchWidth-1],*load);
            beam.resize(options.searchWidth);
        }
        return true;
    };
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
        candidate.labels=candidate.owners;
        int root=-1;
        std::size_t attempts=static_cast<std::size_t>(options.ranks)*options.cutCandidates*8;
        if(!buildTree(graph,vertices,ranks,candidate.labels,options,static_cast<unsigned>(trial),false,candidate.tree,root,attempts)) continue;
        if(!route(graph,candidate,options.ranks)) continue;
        score(graph,candidate,options);retain(std::move(candidate));
    }
    // An independent balanced geometry seed provides regular central cuts
    // even when every METIS label centroid favours a strip arrangement.
    {
        Candidate neutral;neutral.geometrySeed=true;neutral.labels.assign(graph.blocks.size(),ranks.front());
        int root=-1;std::size_t attempts=static_cast<std::size_t>(options.ranks)*options.cutCandidates*8;
        if(buildTree(graph,vertices,ranks,neutral.labels,options,0,true,neutral.tree,root,attempts) &&
           route(graph,neutral,options.ranks))
        {score(graph,neutral,options);retain(std::move(neutral));}
    }
    if(beam.empty())
        throw std::runtime_error("Bounded slicing-tree search found no feasible one-box-per-rank layout.");
    unsigned completed=0,stalled=0;bool exhausted=false;
    for(int round=0;round<options.refinementRounds;++round)
    {
        const auto sources=beam;
        for(const auto& source:sources)
        {
            std::vector<int> nodes;
            for(std::size_t node=0;node<source.tree.size();++node)
                if(source.tree[node].rank<0) nodes.push_back(static_cast<int>(node));
            // Root recuts involve every rank. Other subtrees rotate through the
            // budget rather than revisiting only the heaviest pair.
            const auto budget=std::min(nodes.size(),static_cast<std::size_t>(options.cutCandidates));
            for(std::size_t step=0;step<budget;++step)
            {
                const int node=step==0 ? 0 : nodes[(static_cast<std::size_t>(round)+step)%nodes.size()];
                Candidate proposal;
                if(rebuild(graph,source,node,vertices,source.labels,options,
                           static_cast<unsigned>(round+step+1),proposal))
                    retain(std::move(proposal));
            }
            // Joint positive/negative load compensation: modify two
            // interfaces, including ancestor/descendant pairs sharing ranks, in a single feasible proposal. No intermediate
            // one-cut improvement is required or used as an acceptance gate.
            std::size_t pairs=0;
            for(std::size_t i=0;i<nodes.size() && pairs<static_cast<std::size_t>(options.cutCandidates);++i)
                for(std::size_t j=i+1;j<nodes.size() && pairs<static_cast<std::size_t>(options.cutCandidates);++j)
                {
                    ++pairs;
                    const auto a=adjacentPlanes(graph,source,nodes[i],vertices,options.ranks);
                    const auto b=adjacentPlanes(graph,source,nodes[j],vertices,options.ranks);
                    for(auto first:a) for(auto second:b)
                    {
                        Candidate proposal;proposal.tree=source.tree;proposal.labels=source.labels;proposal.geometrySeed=source.geometrySeed;
                        proposal.tree[nodes[i]].plane=first;proposal.tree[nodes[j]].plane=second;
                        if(!route(graph,proposal,options.ranks)) continue;
                        score(graph,proposal,options);
                        if(proposal.owners!=source.owners) ++jointCandidates;
                        retain(std::move(proposal));
                    }
                }
        }
        ++completed;
        bool changed=sources.size()!=beam.size();
        for(const auto& candidate:beam)
            if(std::none_of(sources.begin(),sources.end(),[&](const auto& item){return item.owners==candidate.owners;})) changed=true;
        stalled=changed ? 0 : stalled+1;
        if(stalled>=static_cast<unsigned>(3*options.cutCandidates)) {exhausted=true;break;}

    }
    std::sort(beam.begin(),beam.end(),[&](const auto& a,const auto& b)
    {return candidateOrder(a,b,bestCompute,options);});
    auto& best=beam.front();
    PartitionPlan result;result.dimension=graph.dimension;
    result.owner=std::move(best.owners);result.blocks=std::move(best.boxes);
    result.cellBlocks=std::move(best.atomToBlock);
    result.blockExchanges=std::move(best.exchanges);result.statistics=std::move(best.stats);
    result.statistics.evaluatedCandidates=evaluated;result.statistics.refinementRounds=completed;
    result.statistics.jointCutCandidates=jointCandidates;result.statistics.bestFeasibleComputeTime=bestCompute;
    result.statistics.searchExhausted=exhausted;
    return result;
}
} // namespace ntic::lbm::mpi
