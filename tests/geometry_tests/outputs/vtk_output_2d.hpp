#pragma once

#include "bounding_box.hpp"
#include "geometry_analysis.hpp"
#include "point.hpp"

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

namespace lat=ntic::lbm::lattice;
namespace geo=ntic::lbm::geometry;


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


//=============================================================================
// 2D scalar field
//=============================================================================

inline void writeVTK2DScalar(
    const std::filesystem::path& path,
    const std::size_t nx,
    const std::size_t ny,
    const geo::Point2D& origin,
    const double spacing,
    const std::vector<double>& scalar)
{
    if(!(spacing>0.0) || !std::isfinite(spacing))
        throw std::invalid_argument("Grid spacing must be positive and finite.");
    if(nx==0 || ny==0)
        throw std::invalid_argument("VTK grid dimensions must be positive.");
    const auto pointCount=vtkCount(vtkProduct(nx,ny));
    if(scalar.size()!=pointCount)
        throw std::invalid_argument("2D scalar size does not match the VTK grid dimensions.");
    auto output=vtkOpen(path);
    output << "# vtk DataFile Version 3.0\nnticLBM 2D scalar field\nBINARY\n"
           << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << nx << ' ' << ny << " 1\n"
           << "ORIGIN " << origin[0] << ' ' << origin[1] << " 0\n"
           << "SPACING " << spacing << ' ' << spacing << " 1\n"
           << "POINT_DATA " << pointCount << '\n';
    vtkScalar<double>(output,"scalar",pointCount,[&](std::size_t i) { return scalar[i]; });
    vtkFinish(output,path);
}

//=============================================================================
// 2D cell scalar field
//=============================================================================

inline void writeVTK2DCellScalar(
    const std::filesystem::path& path,
    const std::size_t nx,
    const std::size_t ny,
    const geo::Point2D& origin,
    const double spacing,
    const std::vector<double>& scalar)
{
    if(!(spacing>0.0) || !std::isfinite(spacing))
        throw std::invalid_argument("Grid spacing must be positive and finite.");
    if(nx==0 || ny==0)
        throw std::invalid_argument("VTK cell dimensions must be positive.");
    vtkCount(nx); vtkCount(ny);
    const auto cellCount=vtkCount(vtkProduct(nx,ny));
    vtkCount(nx+1); vtkCount(ny+1);
    vtkCount(vtkProduct(nx+1,ny+1));
    if(scalar.size()!=cellCount)
        throw std::invalid_argument("2D scalar size does not match the VTK cell dimensions.");
    auto output=vtkOpen(path);
    output << "# vtk DataFile Version 3.0\nnticLBM 2D cell scalar field\nBINARY\n"
           << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << nx+1 << ' ' << ny+1 << " 1\n"
           << "ORIGIN " << origin[0] << ' ' << origin[1] << " 0\n"
           << "SPACING " << spacing << ' ' << spacing << " 1\n"
           << "CELL_DATA " << cellCount << '\n';
    vtkScalar<double>(output,"scalar",cellCount,[&](std::size_t i) { return scalar[i]; });
    vtkFinish(output,path);
}

//=============================================================================
// 2D bounding box
//=============================================================================

inline void writeVTK2DBoundingBox(
    const std::filesystem::path& path,
    const geo::BoundingBox2D& box)
{
    const std::array<std::array<double,3>,4> points{{
        {{box.min[0],box.min[1],0}},{{box.max[0],box.min[1],0}},
        {{box.max[0],box.max[1],0}},{{box.min[0],box.max[1],0}}}};
    const std::array<std::int32_t,12> lines{{2,0,1,2,1,2,2,2,3,2,3,0}};
    auto output=vtkOpen(path);
    output << "# vtk DataFile Version 3.0\nnticLBM 2D bounding box\nBINARY\n"
           << "DATASET POLYDATA\nPOINTS 4 double\n";
    vtkArray<double>(output,12,[&](std::size_t i) { return points[i/3][i%3]; });
    output << "LINES 4 12\n";
    vtkArray<std::int32_t>(output,lines.size(),[&](std::size_t i) { return lines[i]; });
    vtkFinish(output,path);
}

//=============================================================================
// Boundary Q points
//=============================================================================


// Output only Q intersection points. Coincident points with different
// source Face/link records remain separate, preserving provenance.
template<typename Analysis, typename Position, typename Direction>
void vtkWriteBoundaryPoints(
    const std::filesystem::path& path,
    const Analysis& analysis,
    Position position,
    Direction direction)
{
    using Lattice = typename Analysis::Lattice;
    constexpr std::size_t links = Lattice::nStencils - 1;
    struct IntersectionPoint
    {
        std::array<double,3> position;
        double q;
        std::size_t face, link;
    };
    if(analysis.q.size() != vtkProduct(analysis.boundaryX.size(),links) ||
       analysis.boundaryOffsets.size() != analysis.nBoundaries + 1)
        throw std::runtime_error("Incomplete boundary analysis for VTK output.");

    std::vector<IntersectionPoint> points;
    double minimumQ = std::numeric_limits<double>::infinity();
    double maximumQ = 0.0;
    for(std::size_t face = 0; face < analysis.nBoundaries; ++face)
    {
        const auto first = points.size();
        for(std::size_t entry = analysis.boundaryOffsets[face];
            entry < analysis.boundaryOffsets[face + 1]; ++entry)
        {
            const auto start = position(entry);
            for(std::size_t link = 0; link < links; ++link)
            {
                const double q = analysis.q[entry * links + link];
                if(q == -1.0) continue;
                if(!(q > 0.0 && q <= 1.0))
                    throw std::runtime_error("Invalid Q value in boundary output.");
                minimumQ = std::min(minimumQ, q);
                maximumQ = std::max(maximumQ, q);
                auto end = start;
                const auto delta = direction(link + 1);
                for(std::size_t d = 0; d < 3; ++d)
                    end[d] += q * analysis.gridSpacing * delta[d];
                points.push_back({end, q, face, link + 1});
            }
        }
        std::cout << "Face " << face << ": "
                  << analysis.boundaryOffsets[face + 1] - analysis.boundaryOffsets[face]
                  << " entries, " << points.size() - first << " valid Q links\n";
    }

    if(!points.empty())
        std::cout << "Q range: [" << minimumQ << ", " << maximumQ << "]\n";

    vtkCount(points.size());
    const auto connectivityCount=vtkCount(vtkProduct(points.size(),2));
    auto output=vtkOpen(path);
    output << "# vtk DataFile Version 3.0\nBoundary Q intersections\nBINARY\nDATASET POLYDATA\n"
           << "POINTS " << points.size() << " double\n";
    vtkArray<double>(output,vtkProduct(points.size(),3),[&](std::size_t i) {
        return points[i/3].position[i%3];
    });
    output << "VERTICES " << points.size() << ' ' << connectivityCount << '\n';
    vtkArray<std::int32_t>(output,connectivityCount,[](std::size_t i) {
        return vtkInteger<std::int32_t>(i%2==0 ? 1 : i/2);
    });
    output << "POINT_DATA " << points.size() << '\n';
    vtkScalar<double>(output,"Q",points.size(),[&](std::size_t i) { return points[i].q; });
    vtkScalar<std::int32_t>(output,"FaceID",points.size(),[&](std::size_t i) { return points[i].face; });
    vtkScalar<std::int32_t>(output,"LinkID",points.size(),[&](std::size_t i) { return points[i].link; });
    vtkFinish(output,path);
}


template<lat::LatticeType Type>
void writeVTKBoundaryPoints(
    const std::filesystem::path& path, const geo::GeometryAnalysis2D<Type>& analysis)
{
    using Lattice = typename geo::GeometryAnalysis2D<Type>::Lattice;
    vtkWriteBoundaryPoints(path, analysis,
        [&](std::size_t entry) { return std::array<double,3>{
            analysis.boundaryX[entry], analysis.boundaryY[entry], 0.0}; },
        [](std::size_t link) { return std::array<int,3>{
            Lattice::ex[link], Lattice::ey[link], 0}; });
}


//=============================================================================
// Multilevel cell and Q output
//=============================================================================


template<typename Grid>
void vtkWriteMultilevelCells(const std::filesystem::path& path,
    const Grid& grid,bool auxiliary=false)
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
    const Grid& grid,bool auxiliary=false)
{ vtkWriteMultilevelCells(path,grid,auxiliary); }

template<typename Grid,std::enable_if_t<Grid::dimension==2,int> =0>
void writeVTKMultilevelBoundaryPoints(const std::filesystem::path& path,
    const Grid& grid)
{ vtkWriteMultilevelBoundaryPoints(path,grid); }
