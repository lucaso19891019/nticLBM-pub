#pragma once

#include "local_layout.hpp"
#include <mpi.h>

#include <cstring>
#include <string>
#include <string_view>

namespace ntic::lbm::mpi
{
namespace detail
{
class LayoutError
{
    std::array<char,4096> text_{};
public:
    LayoutError& operator=(std::string_view text)
    {
        text_.fill(0);std::copy_n(text.data(),std::min(text.size(),text_.size()-1),text_.data());return *this;
    }
    bool empty() const { return text_[0]==0; }
    operator std::string_view() const { return std::string_view(text_.data()); }
};

class LayoutWriter
{
    static_assert(sizeof(std::size_t)<=sizeof(std::uint64_t),"Layout wire IDs require size_t no wider than uint64.");
public:
    std::vector<unsigned char> bytes;
    void word(std::uint64_t value)
    { for(unsigned i=0;i<8;++i) bytes.push_back(static_cast<unsigned char>(value>>(56-8*i))); }
    void id(std::size_t value)
    { word(value==invalidCellIndex ? std::numeric_limits<std::uint64_t>::max() : static_cast<std::uint64_t>(value)); }
    void real(double value)
    { std::uint64_t bits=0; static_assert(sizeof(value)==sizeof(bits) && std::numeric_limits<double>::is_iec559,"Layout wire doubles require IEEE binary64."); std::memcpy(&bits,&value,8); word(bits); }
    void integer(std::int64_t value)
    { std::uint64_t bits=0; std::memcpy(&bits,&value,8); word(bits); }
    void ids(const std::vector<GlobalCellIndex>& values)
    { word(values.size()); for(auto value:values) id(value); }
};
class LayoutReader
{
    const std::vector<unsigned char>& bytes_;
    std::size_t cursor_=0;
public:
    explicit LayoutReader(const std::vector<unsigned char>& bytes):bytes_(bytes) {}
    std::size_t remaining() const { return bytes_.size()-cursor_; }
    std::uint64_t word()
    {
        if(remaining()<8) throw std::invalid_argument("Truncated layout packet.");
        std::uint64_t value=0;
        for(unsigned i=0;i<8;++i) value=(value<<8)|bytes_[cursor_++];
        return value;
    }
    std::size_t count(std::size_t minimumBytes=8)
    {
        const auto value=word();
        if(value>std::numeric_limits<std::size_t>::max() || value>static_cast<std::uint64_t>(std::numeric_limits<std::ptrdiff_t>::max()) ||
           (minimumBytes && value>remaining()/minimumBytes)) throw std::length_error("Invalid layout packet count.");
        return static_cast<std::size_t>(value);
    }
    std::size_t id()
    {
        const auto value=word();
        if(value==std::numeric_limits<std::uint64_t>::max()) return invalidCellIndex;
        if(value>=static_cast<std::uint64_t>(invalidCellIndex)) throw std::overflow_error("Layout global ID is not representable.");
        return static_cast<std::size_t>(value);
    }
    unsigned unsignedValue()
    {
        const auto value=word(); if(value>std::numeric_limits<unsigned>::max()) throw std::overflow_error("Layout unsigned overflow.");
        return static_cast<unsigned>(value);
    }
    int rank()
    {
        const auto value=word(); if(value>static_cast<std::uint64_t>(std::numeric_limits<int>::max())) throw std::overflow_error("Layout rank overflow.");
        return static_cast<int>(value);
    }
    bool flag()
    { const auto value=word(); if(value>1) throw std::invalid_argument("Invalid layout boolean."); return value!=0; }
    std::uint8_t byte()
    { const auto value=word(); if(value>255) throw std::overflow_error("Layout byte overflow."); return static_cast<std::uint8_t>(value); }
    double real()
    { const auto bits=word(); double value=0; std::memcpy(&value,&bits,8); return value; }
    std::int64_t integer()
    { const auto bits=word(); std::int64_t value=0; std::memcpy(&value,&bits,8); return value; }
    std::vector<GlobalCellIndex> ids()
    {
        std::vector<GlobalCellIndex> values(count());
        for(auto& value:values) value=id();
        return values;
    }
};

inline void layoutMpiCheck(int status,MPI_Comm comm)
{
    if(status!=MPI_SUCCESS) { MPI_Abort(comm,status); std::terminate(); }
}
class LayoutComm
{
    MPI_Comm comm_=MPI_COMM_NULL;
public:
    explicit LayoutComm(MPI_Comm parent)
    {
        layoutMpiCheck(MPI_Comm_dup(parent,&comm_),parent);
        layoutMpiCheck(MPI_Comm_set_errhandler(comm_,MPI_ERRORS_RETURN),comm_);
    }
    LayoutComm(const LayoutComm&)=delete;
    LayoutComm& operator=(const LayoutComm&)=delete;
    ~LayoutComm()
    {
        if(comm_!=MPI_COMM_NULL) {const auto original=comm_;const int status=MPI_Comm_free(&comm_);layoutMpiCheck(status,original);}
    }
    operator MPI_Comm() const { return comm_; }
};

// Fixed-size error broadcast requires no allocation before collective progress.
inline void layoutErrorGate(MPI_Comm comm,std::string_view error)
{
    int rank=0; layoutMpiCheck(MPI_Comm_rank(comm,&rank),comm);
    int first=error.empty() ? std::numeric_limits<int>::max() : rank;
    layoutMpiCheck(MPI_Allreduce(MPI_IN_PLACE,&first,1,MPI_INT,MPI_MIN,comm),comm);
    if(first==std::numeric_limits<int>::max()) return;
    std::array<char,4096> message{};
    if(rank==first) std::copy_n(error.data(),std::min(error.size(),message.size()-1),message.data());
    layoutMpiCheck(MPI_Bcast(message.data(),static_cast<int>(message.size()),MPI_CHAR,first,comm),comm);
    throw std::runtime_error(message.data());
}
inline void layoutSendBytes(const std::vector<unsigned char>& bytes,int rank,int tag,MPI_Comm comm)
{
    for(std::size_t offset=0;offset<bytes.size();)
    {
        const auto count=std::min<std::size_t>(bytes.size()-offset,std::numeric_limits<int>::max());
        layoutMpiCheck(MPI_Send(bytes.data()+offset,static_cast<int>(count),MPI_BYTE,rank,tag,comm),comm); offset+=count;
    }
}
inline void layoutReceiveBytes(std::vector<unsigned char>& bytes,int rank,int tag,MPI_Comm comm)
{
    for(std::size_t offset=0;offset<bytes.size();)
    {
        const auto count=std::min<std::size_t>(bytes.size()-offset,std::numeric_limits<int>::max());
        layoutMpiCheck(MPI_Recv(bytes.data()+offset,static_cast<int>(count),MPI_BYTE,rank,tag,comm,MPI_STATUS_IGNORE),comm); offset+=count;
    }
}
} // namespace detail

template<std::size_t D,lattice::LatticeType Type>
std::vector<unsigned char> encodeLayoutInput(const RankLayoutInput<D,Type>& input)
{
    validateRankLayoutInput(input);
    detail::LayoutWriter w; w.word(0x4e5449434c415931ULL); w.word(1); w.word(D); w.word(static_cast<unsigned>(Type));
    w.word(input.rank); w.word(input.rankCount); w.word(input.frame.levelCount); w.real(input.frame.coarseSpacing);
    for(auto x:input.frame.origin) w.real(x);
    w.word(input.cells.size());
    for(std::size_t i=0;i<input.cells.size();++i)
    {
        const auto& c=input.cells[i]; w.id(c.global); w.id(c.globalSourceNode); w.word(c.level); w.word(c.owner);
        w.word(c.role); w.word(static_cast<unsigned>(c.type)); w.word(c.owned); w.word(c.ghost); w.word(c.halo);
        w.word(c.send); w.word(c.receiveDependent); w.word(input.auxiliaryChild[i]);
        for(auto x:input.coordinates[i]) w.integer(x);
    }
    w.word(input.sourceIDs.size());
    for(std::size_t i=0;i<input.sourceIDs.size();++i)
    { w.id(input.sourceIDs[i]); for(auto id:input.anchors[i]) w.id(id); }
    w.word(input.peers.size()); for(const auto& p:input.peers) {w.word(p.rank);w.ids(p.send);w.ids(p.receive);}
    w.ids(input.boundaries.offsets); w.ids(input.boundaries.cells); w.word(input.boundaries.q.size()); for(auto q:input.boundaries.q) w.real(q);
    w.word(input.interfaces.size()); for(const auto& b:input.interfaces) {w.id(b.coarse);w.ids(b.auxiliary);w.ids(b.fine);}
    w.word(input.intraRankExchanges.size());
    for(const auto& e:input.intraRankExchanges)
    { w.id(e.a);w.id(e.b);w.word(e.dependencies.size());for(const auto& p:e.dependencies) {w.id(p.first);w.id(p.second);} }
    w.word(input.nodes.size());
    for(const auto& n:input.nodes)
    {
        w.id(n.globalNode);w.id(n.parentGlobal);w.id(n.firstGlobalChild);w.id(n.globalRealCell);w.word(n.level);
        for(auto x:n.coordinate) w.word(x);
        w.word(n.active);w.word(n.leaf);for(auto child:n.globalChildren) w.id(child);
    }
    return std::move(w.bytes);
}

template<std::size_t D,lattice::LatticeType Type>
RankLayoutInput<D,Type> decodeLayoutInput(const std::vector<unsigned char>& bytes)
{
    detail::LayoutReader r(bytes);
    if(r.word()!=0x4e5449434c415931ULL || r.word()!=1 || r.word()!=D || r.word()!=static_cast<unsigned>(Type))
        throw std::invalid_argument("Layout packet version/dimension/lattice mismatch.");
    RankLayoutInput<D,Type> input; input.rank=r.rank();input.rankCount=r.rank();input.frame.levelCount=r.unsignedValue();input.frame.coarseSpacing=r.real();
    for(auto& x:input.frame.origin) x=r.real();
    const auto cells=r.count((12+D)*8); input.cells.resize(cells);input.coordinates.resize(cells);input.auxiliaryChild.resize(cells);
    for(std::size_t i=0;i<cells;++i)
    {
        auto& c=input.cells[i];c.global=r.id();c.globalSourceNode=r.id();c.level=r.unsignedValue();c.owner=r.rank();
        c.role=r.byte();const auto type=r.unsignedValue();if(type>2) throw std::invalid_argument("Invalid layout cell type.");
        c.type=static_cast<geometry::CellType>(type);c.owned=r.flag();c.ghost=r.flag();c.halo=r.flag();c.send=r.flag();c.receiveDependent=r.flag();
        input.auxiliaryChild[i]=r.byte();for(auto& x:input.coordinates[i]) x=r.integer();
    }
    const auto rows=r.count((lattice::Twist<Type>::nAnchors+1)*8);input.sourceIDs.resize(rows);input.anchors.resize(rows);
    for(std::size_t i=0;i<rows;++i) {input.sourceIDs[i]=r.id();for(auto& id:input.anchors[i]) id=r.id();}
    input.peers.resize(r.count(24));for(auto& p:input.peers) {p.rank=r.rank();p.send=r.ids();p.receive=r.ids();}
    input.boundaries.offsets=r.ids();input.boundaries.cells=r.ids();input.boundaries.q.resize(r.count());for(auto& q:input.boundaries.q) q=r.real();
    input.interfaces.resize(r.count(24));for(auto& b:input.interfaces) {b.coarse=r.id();b.auxiliary=r.ids();b.fine=r.ids();}
    input.intraRankExchanges.resize(r.count(24));
    for(auto& e:input.intraRankExchanges) {e.a=r.id();e.b=r.id();e.dependencies.resize(r.count(16));for(auto& p:e.dependencies) {p.first=r.id();p.second=r.id();}}
    input.nodes.resize(r.count((7+D+(std::size_t{1}<<D))*8));
    for(auto& n:input.nodes)
    {
        n.globalNode=r.id();n.parentGlobal=r.id();n.firstGlobalChild=r.id();n.globalRealCell=r.id();n.level=r.unsignedValue();
        for(auto& x:n.coordinate) x=r.unsignedValue();
        n.active=r.flag();n.leaf=r.flag();
        for(auto& child:n.globalChildren) child=r.id();
    }
    if(r.remaining()) throw std::invalid_argument("Trailing layout packet data.");
    validateRankLayoutInput(input); return input;
}

template<std::size_t D,lattice::LatticeType Type>
RankLayoutInput<D,Type> distributeLayoutInput(const LocalLayoutBuilder<D,Type>* rootBuilder,MPI_Comm communicator,
                                             int root=0,const std::string& rootError={})
{
    detail::LayoutComm privateComm(communicator); MPI_Comm comm=privateComm;
    int rank=0,size=0;detail::layoutMpiCheck(MPI_Comm_rank(comm,&rank),comm);detail::layoutMpiCheck(MPI_Comm_size(comm,&size),comm);
    detail::LayoutError error; if(root<0 || root>=size) error="Layout distribution root is out of range.";
    int minimumRoot=root,maximumRoot=root;
    detail::layoutMpiCheck(MPI_Allreduce(MPI_IN_PLACE,&minimumRoot,1,MPI_INT,MPI_MIN,comm),comm);
    detail::layoutMpiCheck(MPI_Allreduce(MPI_IN_PLACE,&maximumRoot,1,MPI_INT,MPI_MAX,comm),comm);
    if(minimumRoot!=maximumRoot) error="MPI ranks disagree on layout root.";
    detail::layoutErrorGate(comm,error); RankLayoutInput<D,Type> input;
    if(rank==root)
    {
        error=rootError;
        if(error.empty() && !rootBuilder) error="Root layout builder is missing.";
        for(int destination=0;destination<size;++destination)
        {
            std::vector<unsigned char> packet;
            if(error.empty())
            {
                try
                {
                    auto prepared=rootBuilder->prepare(destination);
                    if(prepared.rankCount!=size) throw std::invalid_argument("Partition rank count differs from MPI communicator size.");
                    if(destination==root) input=std::move(prepared); else packet=encodeLayoutInput(prepared);
                }
                catch(const std::exception& e) {error=e.what();}
                catch(...) {error="Unknown root layout preparation error.";}
            }
            if(destination==root) continue;
            std::array<std::uint64_t,2> header{error.empty() ? 0U : 1U,static_cast<std::uint64_t>(packet.size())};
            detail::layoutMpiCheck(MPI_Send(header.data(),2,MPI_UINT64_T,destination,0,comm),comm);
            if(header[0]==0)
            {
                int ready=0;detail::layoutMpiCheck(MPI_Recv(&ready,1,MPI_INT,destination,1,comm,MPI_STATUS_IGNORE),comm);
                if(ready) detail::layoutSendBytes(packet,destination,2,comm);
            }
        }
    }
    else
    {
        std::array<std::uint64_t,2> header{};
        detail::layoutMpiCheck(MPI_Recv(header.data(),2,MPI_UINT64_T,root,0,comm,MPI_STATUS_IGNORE),comm);
        if(header[0]==0)
        {
            std::vector<unsigned char> packet;int ready=1;
            try
            {
                if(header[1]>packet.max_size()) throw std::length_error("Layout packet is too large for this rank.");
                packet.resize(static_cast<std::size_t>(header[1]));
            }
            catch(const std::exception& e) {ready=0;error=e.what();}
            catch(...) {ready=0;error="Layout packet allocation failed.";}
            detail::layoutMpiCheck(MPI_Send(&ready,1,MPI_INT,root,1,comm),comm);
            if(ready)
            {
                detail::layoutReceiveBytes(packet,root,2,comm);
                try {input=decodeLayoutInput<D,Type>(packet);if(input.rank!=rank || input.rankCount!=size) throw std::invalid_argument("Layout packet addressed to a different rank.");}
                catch(const std::exception& e) {error=e.what();}
                catch(...) {error="Unknown layout decode error.";}
            }
        }
    }
    detail::layoutErrorGate(comm,error);return input;
}

namespace detail
{
template<std::size_t D,lattice::LatticeType Type>
void checkLayoutPeers(const LocalLayout<D,Type>& layout,MPI_Comm comm)
{
    int rank=0,size=0;layoutMpiCheck(MPI_Comm_rank(comm,&rank),comm);layoutMpiCheck(MPI_Comm_size(comm,&size),comm);
    LayoutError error;
    std::vector<std::uint64_t> sendCounts,receiveCounts;
    std::vector<const PeerExchange*> peers;
    try
    {
        if(layout.rank!=rank) throw std::invalid_argument("Local layout rank differs from MPI rank.");
        sendCounts.assign(static_cast<std::size_t>(size),0);receiveCounts.resize(sendCounts.size());peers.assign(sendCounts.size(),nullptr);
        for(const auto& peer:layout.peers)
        {
            if(peer.rank<0 || peer.rank>=size || peer.rank==rank || peers[peer.rank]) throw std::invalid_argument("Invalid MPI layout peer.");
            peers[peer.rank]=&peer;sendCounts[peer.rank]=peer.send.size();
            auto check=[&](const auto& indices) {
                GlobalCellIndex previous=invalidCellIndex;
                for(auto index:indices)
                {
                    const auto global=layout.cells.at(index).global;
                    if(previous!=invalidCellIndex && previous>=global) throw std::invalid_argument("Peer storage IDs are not canonical.");
                    previous=global;
                }
            };
            check(peer.send);check(peer.receive);
        }
    }
    catch(const std::exception& e) {error=e.what();}
    catch(...) {error="MPI peer metadata allocation failed.";}
    layoutErrorGate(comm,error);
    layoutMpiCheck(MPI_Alltoall(sendCounts.data(),1,MPI_UINT64_T,receiveCounts.data(),1,MPI_UINT64_T,comm),comm);
    std::vector<std::uint64_t> outgoing,incoming;
    try
    {
        std::uint64_t maxSend=0,maxReceive=0;
        for(int peer=0;peer<size;++peer)
        {
            const auto expected=peers[peer] ? peers[peer]->receive.size() : 0;
            if(receiveCounts[peer]!=expected) throw std::invalid_argument("MPI peer storage counts do not match.");
            maxSend=std::max(maxSend,sendCounts[peer]);maxReceive=std::max(maxReceive,receiveCounts[peer]);
        }
        if(maxSend>outgoing.max_size() || maxReceive>incoming.max_size()) throw std::length_error("MPI peer storage list is too large.");
        outgoing.resize(static_cast<std::size_t>(maxSend));incoming.resize(static_cast<std::size_t>(maxReceive));
    }
    catch(const std::exception& e) {error=e.what();}
    catch(...) {error="MPI peer buffer allocation failed.";}
    layoutErrorGate(comm,error);
    // Identical lexicographic rank-pair schedule avoids unmatched blocking calls.
    // IDs, not hashes, are compared, independently of local memory order.
    for(int a=0;a<size;++a) for(int b=a+1;b<size;++b)
    {
        if(rank!=a && rank!=b) continue;
        const int peer=rank==a ? b : a;
        const auto sending=static_cast<std::size_t>(sendCounts[peer]),receiving=static_cast<std::size_t>(receiveCounts[peer]);
        if(sending==0 && receiving==0) continue;
        for(std::size_t i=0;i<sending;++i) outgoing[i]=layout.cells[peers[peer]->send[i]].global;
        const auto total=std::max(sending,receiving);
        for(std::size_t offset=0;offset<total;)
        {
            const auto remaining=total-offset;
            const auto advance=std::min<std::size_t>(remaining,std::numeric_limits<int>::max());
            const auto send=offset<sending ? std::min(advance,sending-offset) : 0;
            const auto receive=offset<receiving ? std::min(advance,receiving-offset) : 0;
            layoutMpiCheck(MPI_Sendrecv(send ? outgoing.data()+offset : nullptr,static_cast<int>(send),MPI_UINT64_T,peer,3,
                receive ? incoming.data()+offset : nullptr,static_cast<int>(receive),MPI_UINT64_T,peer,3,comm,MPI_STATUS_IGNORE),comm);
            offset+=advance;
        }
        for(std::size_t i=0;i<receiving;++i)
            if(incoming[i]!=layout.cells[peers[peer]->receive[i]].global) error="MPI peer global storage IDs do not match.";
    }
    layoutErrorGate(comm,error);
}
} // namespace detail

template<std::size_t D,lattice::LatticeType Type>
LocalLayout<D,Type> initializeLocalLayout(const RankLayoutInput<D,Type>& input,MPI_Comm communicator,
                                         const LocalLayoutOptions& options={})
{
    detail::LayoutComm privateComm(communicator);MPI_Comm comm=privateComm;
    int rank=0,size=0;detail::layoutMpiCheck(MPI_Comm_rank(comm,&rank),comm);detail::layoutMpiCheck(MPI_Comm_size(comm,&size),comm);
    detail::LayoutError error;LocalLayout<D,Type> layout;
    try
    {
        if(input.rank!=rank || input.rankCount!=size) throw std::invalid_argument("Rank input differs from its MPI communicator.");
        layout=buildLocalLayout(input,options);
    }
    catch(const std::exception& e) {error=e.what();}
    catch(...) {error="Unknown local layout construction error.";}
    detail::layoutErrorGate(comm,error);detail::checkLayoutPeers(layout,comm);return layout;
}

} // namespace ntic::lbm::mpi
