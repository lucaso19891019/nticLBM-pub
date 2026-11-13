#pragma once

#include "partition_types.hpp"
#include "lattice_model.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <ostream>
#include <stdexcept>
#include <system_error>
#include <type_traits>
#include <vector>
#include <functional>
#include <string>

namespace lat=ntic::lbm::lattice;
namespace part=ntic::lbm::mpi;


inline void recreateOutputDirectory(
    const std::filesystem::path& path)
{
    std::error_code error;


    const bool exists =
        std::filesystem::exists(
            path,
            error);

    if(error)
    {
        throw std::runtime_error(
            "Failed to check output directory: " +
            path.string());
    }


    if(exists)
    {
        std::filesystem::remove_all(
            path,
            error);

        if(error)
        {
            throw std::runtime_error(
                "Failed to remove output directory: " +
                path.string());
        }
    }


    std::filesystem::create_directories(
        path,
        error);

    if(error)
    {
        throw std::runtime_error(
            "Failed to create output directory: " +
            path.string());
    }
}


// Legacy VTK binary arrays use big-endian fixed-width values.
inline std::size_t vtkProduct(std::size_t a,std::size_t b)
{
    if(b && a>std::numeric_limits<std::size_t>::max()/b)
        throw std::overflow_error("VTK element count overflow.");
    return a*b;
}

inline std::size_t vtkCount(std::size_t count)
{
    if(count>static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
        throw std::overflow_error("Legacy VTK count exceeds int32.");
    return count;
}

template<typename T,typename U>
T vtkInteger(U value)
{
    static_assert(std::is_integral_v<T> && std::is_integral_v<U>);
    if constexpr(std::is_signed_v<U>)
        if(value<0)
        {
            if constexpr(std::is_unsigned_v<T>)
                throw std::overflow_error("Negative VTK unsigned index.");
            else if(static_cast<std::intmax_t>(value)<std::numeric_limits<T>::min())
                throw std::overflow_error("VTK integer field is out of range.");
            return static_cast<T>(value);
        }
    if(static_cast<std::uintmax_t>(value)>static_cast<std::uintmax_t>(std::numeric_limits<T>::max()))
        throw std::overflow_error("VTK integer field is out of range.");
    return static_cast<T>(value);
}

inline std::ofstream vtkOpen(const std::filesystem::path& path)
{
    std::ofstream output(path,std::ios::binary);
    if(!output) throw std::runtime_error("Cannot open " + path.string());
    output.imbue(std::locale::classic());
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    return output;
}

inline void vtkFinish(std::ofstream& output,const std::filesystem::path& path)
{
    output.flush();
    if(!output) throw std::runtime_error("Failed to write " + path.string());
}

template<typename T,typename Getter>
void vtkArray(std::ostream& output,std::size_t count,Getter value)
{
    static_assert(std::is_same_v<T,double> || std::is_same_v<T,std::int32_t> || std::is_same_v<T,std::uint32_t>);
    static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559);
    const std::uint16_t one=1;
    const bool littleEndian=*reinterpret_cast<const unsigned char*>(&one)==1;
    std::array<T,8192> buffer;
    for(std::size_t first=0; first<count;)
    {
        const auto size=std::min(buffer.size(),count-first);
        for(std::size_t i=0; i<size; ++i)
        {
            if constexpr(std::is_integral_v<T>) buffer[i]=vtkInteger<T>(value(first+i));
            else buffer[i]=value(first+i);
        }
        auto* bytes=reinterpret_cast<unsigned char*>(buffer.data());
        if(littleEndian)
            for(std::size_t i=0; i<size; ++i) std::reverse(bytes+i*sizeof(T),bytes+(i+1)*sizeof(T));
        output.write(reinterpret_cast<const char*>(buffer.data()),static_cast<std::streamsize>(size*sizeof(T)));
        first+=size;
    }
    output.put('\n');
    if(!output) throw std::runtime_error("Binary VTK payload write failed.");
}

template<typename T,typename Getter>
void vtkScalar(std::ostream& output,const char* name,std::size_t count,Getter value)
{
    const char* type=std::is_same_v<T,double> ? "double" : std::is_same_v<T,std::uint32_t> ? "unsigned_int" : "int";
    output << "SCALARS " << name << ' ' << type << " 1\nLOOKUP_TABLE default\n";
    vtkArray<T>(output,count,value);
}


// Optional fields use the same selected CellIndex order as the base fields.
using MultilevelCellDataWriter=
    std::function<void(std::ostream&,const std::vector<std::size_t>&)>;


template<typename Grid>
void vtkWriteMultilevelCells(const std::filesystem::path& path,
    const Grid& grid,bool auxiliary=false,
    const MultilevelCellDataWriter& writeAdditionalCellData={})
{
    constexpr auto D=Grid::dimension;
    const auto& cells=grid.cells();
    std::vector<std::size_t> selected;
    for(std::size_t c=0; c<cells.size(); ++c) if(cells.isAuxiliary(c)==auxiliary) selected.push_back(c);
    constexpr std::size_t corners=std::size_t{1}<<D;
    constexpr std::array<unsigned,8> order{{0,1,3,2,4,5,7,6}};
    const auto pointCount=vtkCount(vtkProduct(selected.size(),corners));
    const auto connectivityCount=vtkCount(vtkProduct(selected.size(),corners+1));
    auto out=vtkOpen(path);
    out << "# vtk DataFile Version 3.0\nnticLBM multilevel cells\nBINARY\n"
        << "DATASET UNSTRUCTURED_GRID\nPOINTS " << pointCount << " double\n";
    vtkArray<double>(out,vtkProduct(pointCount,3),[&](std::size_t i) {
        const auto cell=selected[i/(3*corners)];
        const auto corner=(i/3)%corners,axis=i%3;
        const auto p=grid.coordinate(cell); const double h=grid.spacing(grid.level(cell));
        return axis<D ? grid.domain().min[axis]+(p[axis]+((order[corner]>>axis)&1U))*h : 0;
    });
    out << "CELLS " << selected.size() << ' ' << connectivityCount << '\n';
    vtkArray<std::int32_t>(out,connectivityCount,[&](std::size_t i) {
        const auto offset=i%(corners+1);
        return vtkInteger<std::int32_t>(offset==0 ? corners : (i/(corners+1))*corners+offset-1);
    });
    out << "CELL_TYPES " << selected.size() << '\n';
    vtkArray<std::int32_t>(out,selected.size(),[](std::size_t) { return D==2 ? 9 : 12; });
    out << "CELL_DATA " << selected.size() << '\n';
    vtkScalar<std::int32_t>(out,"scalar",selected.size(),[&](std::size_t i) { return static_cast<int>(cells.cellType[selected[i]]); });
    vtkScalar<std::int32_t>(out,"level",selected.size(),[&](std::size_t i) { return grid.level(selected[i]); });
    vtkScalar<std::int32_t>(out,"role",selected.size(),[&](std::size_t i) { return static_cast<unsigned>(cells.role[selected[i]]); });
    vtkScalar<std::uint32_t>(out,"CellIndex",selected.size(),[&](std::size_t i) { return selected[i]; });
    vtkScalar<std::uint32_t>(out,"SourceNode",selected.size(),[&](std::size_t i) { return cells.sourceNode[selected[i]]; });
    if(writeAdditionalCellData) writeAdditionalCellData(out,selected);
    vtkFinish(out,path);
}

template<typename Grid>
void vtkWriteMultilevelBoundaryPoints(const std::filesystem::path& path,
    const Grid& grid)
{
    constexpr auto D=Grid::dimension;
    using Lattice=lat::Model<Grid::latticeType>;
    constexpr auto links=Lattice::nStencils-1;
    const auto& boundary=grid.boundaries();
    struct Hit { std::size_t entry,face,direction; };
    std::vector<Hit> hits;
    for(std::size_t face=0; face<boundary.nBoundaries; ++face)
        for(auto e=boundary.offsets[face]; e<boundary.offsets[face+1]; ++e)
            for(std::size_t i=1; i<=links; ++i)
            {
                const double q=boundary.q[e*links+i-1];
                if(q==-1) continue;
                if(!(q>0 && q<=1)) throw std::runtime_error("Invalid multilevel boundary Q.");
                hits.push_back({e,face,i});
            }
    vtkCount(hits.size());
    const auto connectivityCount=vtkCount(vtkProduct(hits.size(),2));
    auto out=vtkOpen(path);
    out << "# vtk DataFile Version 3.0\nMultilevel Q intersections\nBINARY\n"
        << "DATASET POLYDATA\nPOINTS " << hits.size() << " double\n";
    vtkArray<double>(out,vtkProduct(hits.size(),3),[&](std::size_t i) {
        const auto& hit=hits[i/3]; const auto axis=i%3;
        const auto c=boundary.cells[hit.entry]; auto p=grid.center(c);
        const double distance=boundary.q[hit.entry*links+hit.direction-1]*grid.spacing(grid.level(c));
        p[0]+=distance*Lattice::ex[hit.direction]; p[1]+=distance*Lattice::ey[hit.direction];
        if constexpr(D==3) p[2]+=distance*Lattice::ez[hit.direction];
        return axis<D ? p[axis] : 0;
    });
    out << "VERTICES " << hits.size() << ' ' << connectivityCount << '\n';
    vtkArray<std::int32_t>(out,connectivityCount,[](std::size_t i) { return vtkInteger<std::int32_t>(i%2==0 ? 1 : i/2); });
    out << "POINT_DATA " << hits.size() << '\n';
    vtkScalar<double>(out,"Q",hits.size(),[&](std::size_t i) { return boundary.q[hits[i].entry*links+hits[i].direction-1]; });
    vtkScalar<std::int32_t>(out,"FaceID",hits.size(),[&](std::size_t i) { return hits[i].face; });
    vtkScalar<std::int32_t>(out,"LinkID",hits.size(),[&](std::size_t i) { return hits[i].direction; });
    vtkScalar<std::uint32_t>(out,"CellIndex",hits.size(),[&](std::size_t i) { return boundary.cells[hits[i].entry]; });
    vtkScalar<std::int32_t>(out,"level",hits.size(),[&](std::size_t i) { return grid.level(boundary.cells[hits[i].entry]); });
    vtkScalar<std::int32_t>(out,"role",hits.size(),[&](std::size_t i) { return static_cast<unsigned>(grid.cells().role[boundary.cells[hits[i].entry]]); });
    vtkFinish(out,path);
}


template<typename Grid,std::enable_if_t<Grid::dimension==2,int> =0>
void writeVTKMultilevelCells(const std::filesystem::path& path,
    const Grid& grid,bool auxiliary=false,
    const MultilevelCellDataWriter& writeAdditionalCellData={})
{ vtkWriteMultilevelCells(path,grid,auxiliary,writeAdditionalCellData); }

template<typename Grid,std::enable_if_t<Grid::dimension==2,int> =0>
void writeVTKMultilevelBoundaryPoints(const std::filesystem::path& path,
    const Grid& grid)
{ vtkWriteMultilevelBoundaryPoints(path,grid); }


inline int vtkPartitionNode(const part::HardwareProfile& hardware,int rank)
{
    return hardware.nodeIds.empty() ? -1 : hardware.nodeIds.at(static_cast<std::size_t>(rank));
}

inline std::filesystem::path vtkPreparePartitionOutputDirectory(
    const std::filesystem::path& root,const std::string& caseName)
{
    if(caseName!="test_partition_taylor_couette_vtks" &&
       caseName!="test_partition_branched_channel_vtks")
        throw std::invalid_argument("Unknown manual partition output directory.");
    const auto outputRoot=std::filesystem::weakly_canonical(root);
    const auto directory=outputRoot/caseName;
    // Refuse links or aliases to another output location before removing data.
    if(directory.parent_path()!=outputRoot || std::filesystem::weakly_canonical(directory)!=directory)
        throw std::runtime_error("Manual partition output resolves outside its dedicated directory.");
    recreateOutputDirectory(directory);
    return directory;
}

template<typename Grid>
void vtkWritePartitionCells(const std::filesystem::path& path,
    const Grid& grid,const part::PartitionPlan& plan,
    const part::HardwareProfile& hardware,bool auxiliary=false)
{
    if(plan.owner.size()!=grid.cells().size() || plan.cellBlocks.size()!=grid.cells().size())
        throw std::invalid_argument("Partition output does not match the grid.");
    vtkWriteMultilevelCells(path,grid,auxiliary,
        [&](std::ostream& out,const std::vector<std::size_t>& selected) {
            vtkScalar<std::int32_t>(out,"rank",selected.size(),[&](std::size_t i) { return plan.owner[selected[i]]; });
            vtkScalar<std::int32_t>(out,"node",selected.size(),[&](std::size_t i) { return vtkPartitionNode(hardware,plan.owner[selected[i]]); });
            vtkScalar<std::uint32_t>(out,"BlockIndex",selected.size(),[&](std::size_t i) { return plan.cellBlocks[selected[i]]; });
        });
}

inline void vtkWritePartitionBlocks(const std::filesystem::path& path,
    const part::PartitionPlan& plan,const part::HardwareProfile& hardware)
{
    if(plan.dimension!=2 && plan.dimension!=3)
        throw std::invalid_argument("Partition block output requires two or three dimensions.");
    const std::size_t corners=std::size_t{1}<<plan.dimension;
    constexpr std::array<unsigned,8> order{{0,1,3,2,4,5,7,6}};
    const auto pointCount=vtkCount(vtkProduct(plan.blocks.size(),corners));
    const auto connectivityCount=vtkCount(vtkProduct(plan.blocks.size(),corners+1));
    auto out=vtkOpen(path);
    out << "# vtk DataFile Version 3.0\nnticLBM rank partition boxes\nBINARY\n"
        << "DATASET UNSTRUCTURED_GRID\nPOINTS " << pointCount << " double\n";
    vtkArray<double>(out,vtkProduct(pointCount,3),[&](std::size_t i) {
        const auto& block=plan.blocks[i/(3*corners)];
        const auto corner=(i/3)%corners,axis=i%3;
        const auto index=((order[corner]>>axis)&1U)!=0 ? block.upper[axis] : block.lower[axis];
        return axis<plan.dimension ? plan.domainMin[axis]+index*plan.blockSpacing : 0;
    });
    out << "CELLS " << plan.blocks.size() << ' ' << connectivityCount << '\n';
    vtkArray<std::int32_t>(out,connectivityCount,[&](std::size_t i) {
        const auto offset=i%(corners+1);
        return vtkInteger<std::int32_t>(offset==0 ? corners : (i/(corners+1))*corners+offset-1);
    });
    out << "CELL_TYPES " << plan.blocks.size() << '\n';
    vtkArray<std::int32_t>(out,plan.blocks.size(),[&](std::size_t) { return plan.dimension==2 ? 9 : 12; });
    out << "CELL_DATA " << plan.blocks.size() << '\n';
    vtkScalar<std::int32_t>(out,"rank",plan.blocks.size(),[&](std::size_t i) { return plan.blocks[i].owner; });
    vtkScalar<std::int32_t>(out,"node",plan.blocks.size(),[&](std::size_t i) { return vtkPartitionNode(hardware,plan.blocks[i].owner); });
    vtkScalar<std::uint32_t>(out,"BlockIndex",plan.blocks.size(),[](std::size_t i) { return i; });
    vtkScalar<double>(out,"work",plan.blocks.size(),[&](std::size_t i) { return plan.blocks[i].work; });
    vtkFinish(out,path);
}

template<typename Grid>
void vtkWritePartitionSummary(std::ostream& out,const Grid& grid,
    const part::PartitionPlan& plan,const part::PartitionOptions& options,const std::string& description,unsigned analysisThreads)
{
    const auto& stats=plan.statistics;
    out << std::setprecision(10) << description << '\n'
        << "Partition metadata only: one box per rank with bounded global refinement; no MPI transport, LBM evolution or communication hiding.\n"
        << "Hardware costs are uncalibrated model units, not measured seconds.\n"
        << "Node -1 means unspecified; it does not assert a physical node assignment.\n"
        << "Ranks: " << options.ranks << ", blockLevel: " << options.blockLevel
        << ", candidates: " << options.candidates << ", boundary extra per Face record: " << options.boundaryExtraWeight << '\n'
        << "Local exchange latency/bandwidth: " << options.hardware.intraRankLatency << " / " << options.hardware.intraRankBandwidth << '\n'
        << "Same-node exchange latency/bandwidth: " << options.hardware.intraNodeLatency << " / " << options.hardware.intraNodeBandwidth << '\n'
        << "Cross-node exchange latency/bandwidth: " << options.hardware.interNodeLatency << " / " << options.hardware.interNodeBandwidth << '\n'
        << "Analysis thread cap: " << analysisThreads << ", graph thread cap: " << options.graphThreads << '\n'
        << "Tree nodes: " << grid.nodes().size() << ", retained slots: " << grid.cells().size()
        << ", coarse/fine interface blocks: " << grid.interfaces().blocks.size() << '\n'
        << "Work per level below counts slot work; coarse/fine pair work is reported separately.\n";
    for(int rank=0; rank<options.ranks; ++rank)
    {
        const double rate=options.hardware.computeRates.empty() ? 1 : options.hardware.computeRates[rank];
        out << "Rank " << rank << ": node=" << vtkPartitionNode(options.hardware,rank)
            << ", computeRate=" << rate << ", wet=" << stats.realCells[rank]
            << ", auxiliary=" << stats.auxiliaryCells[rank] << ", blocks=" << stats.blocksPerRank[rank]
            << ", work=" << stats.rankWeights[rank] << ", interfaceWork=" << stats.interfaceWork[rank]
            << ", computeTime=" << stats.computeTimes[rank]
            << ", localExchangeTime=" << stats.intraRankCommunicationTimes[rank]
            << ", remoteExchangeTime=" << stats.interRankCommunicationTimes[rank] << '\n';
        for(std::size_t level=0; level<stats.workPerLevel[rank].size(); ++level)
            out << "  level " << level << ": h=" << grid.spacing(static_cast<unsigned>(level))
                << ", work=" << stats.workPerLevel[rank][level] << '\n';
    }
    out << "Overall max-to-mean compute ratio: " << stats.maxToMean
        << "\nMaximum compute time: " << stats.maxComputeTime
        << "\nMaximum estimated compute + exchange time: " << stats.maxEstimatedTime
        << "\nGraph vertices/edges: " << stats.graphVertices << " / " << stats.graphEdges
        << "\nLargest atomic work: " << stats.largestAtomicWork
        << "\nBest feasible compute time: " << stats.bestFeasibleComputeTime
        << "\nEvaluated candidates: " << stats.evaluatedCandidates
        << "\nRefinement rounds: " << stats.refinementRounds
        << "\nJoint-cut candidates: " << stats.jointCutCandidates
        << "\nSearch termination: " << (stats.searchExhausted ? "unchanged beam proposal cycle" : "round budget")
        << "\nRank boxes: " << plan.blocks.size() << '\n';
}

template<typename Grid>
void vtkWritePartitionOutputs(const std::filesystem::path& root,const std::string& caseName,
    const Grid& grid,const part::PartitionPlan& plan,
    const part::PartitionOptions& options,const std::string& description,unsigned analysisThreads)
{
    const auto directory=vtkPreparePartitionOutputDirectory(root,caseName);
    vtkWritePartitionCells(directory/"partition.vtk",grid,plan,options.hardware);
    vtkWritePartitionCells(directory/"auxiliary.vtk",grid,plan,options.hardware,true);
    vtkWritePartitionBlocks(directory/"blocks.vtk",plan,options.hardware);
    vtkWriteMultilevelBoundaryPoints(directory/"boundary_q.vtk",grid);
    vtkWritePartitionSummary(std::cout,grid,plan,options,description,analysisThreads);
    std::ofstream report(directory/"partition_stats.txt");
    if(!report) throw std::runtime_error("Cannot open partition statistics output.");
    vtkWritePartitionSummary(report,grid,plan,options,description,analysisThreads);
    if(!report) throw std::runtime_error("Failed to write partition statistics output.");
    std::cout << "Output directory: " << directory << '\n';
}


template<typename Grid,std::enable_if_t<Grid::dimension==2,int> =0>
void writeVTKPartitionCells(const std::filesystem::path& path,
    const Grid& grid,const part::PartitionPlan& plan,
    const part::HardwareProfile& hardware,bool auxiliary=false)
{ vtkWritePartitionCells(path,grid,plan,hardware,auxiliary); }

template<typename Grid,std::enable_if_t<Grid::dimension==2,int> =0>
void writePartitionSummary(std::ostream& out,const Grid& grid,
    const part::PartitionPlan& plan,const part::PartitionOptions& options,const std::string& description,unsigned analysisThreads)
{ vtkWritePartitionSummary(out,grid,plan,options,description,analysisThreads); }

template<typename Grid,std::enable_if_t<Grid::dimension==2,int> =0>
void writePartitionOutputs(const std::filesystem::path& root,const std::string& caseName,
    const Grid& grid,const part::PartitionPlan& plan,
    const part::PartitionOptions& options,const std::string& description,unsigned analysisThreads)
{ vtkWritePartitionOutputs(root,caseName,grid,plan,options,description,analysisThreads); }

// PartitionPlan stores its dimension at runtime; preserve its direct writer API.
inline void writeVTKPartitionBlocks(const std::filesystem::path& path,
    const part::PartitionPlan& plan,const part::HardwareProfile& hardware)
{ vtkWritePartitionBlocks(path,plan,hardware); }
