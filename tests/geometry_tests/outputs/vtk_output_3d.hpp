#pragma once

#include "vtk_output_2d.hpp"

namespace ntic::lbm::geometry::test
{

template<typename FluidFunction>
void writeVTK3DGeometry(
    const std::filesystem::path& path,
    const BoundingBox& domainBounds,
    const double spacing,
    FluidFunction isFluid)
{
    if(!(spacing>0.0) || !std::isfinite(spacing))
        throw std::invalid_argument("Grid spacing must be positive and finite.");
    const auto dimension=[&](double width) {
        const double cells=std::floor(width/spacing);
        if(!(cells>=0.0) || !std::isfinite(cells) ||
           cells>=static_cast<double>(std::numeric_limits<std::int32_t>::max()))
            throw std::invalid_argument("Invalid VTK geometry dimension.");
        return static_cast<std::size_t>(cells)+1;
    };
    const auto nx=dimension(domainBounds.width());
    const auto ny=dimension(domainBounds.height());
    const auto nz=dimension(domainBounds.depth());
    const auto xySize=detail::vtkProduct(nx,ny);
    const auto pointCount=detail::vtkCount(detail::vtkProduct(xySize,nz));
    std::vector<std::int32_t> fluid(pointCount,0);
#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index=0; index<static_cast<std::ptrdiff_t>(pointCount); ++index)
    {
        const auto id=static_cast<std::size_t>(index);
        const Point point{
            domainBounds.min[0]+static_cast<double>(id%nx)*spacing,
            domainBounds.min[1]+static_cast<double>((id/nx)%ny)*spacing,
            domainBounds.min[2]+static_cast<double>(id/xySize)*spacing};
        fluid[id]=isFluid(point) ? 1 : 0;
    }
    auto output=detail::vtkOpen(path);
    output << "# vtk DataFile Version 3.0\nnticLBM 3D geometry test\nBINARY\n"
           << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << nx << ' ' << ny << ' ' << nz << '\n'
           << "ORIGIN " << domainBounds.min[0] << ' ' << domainBounds.min[1] << ' ' << domainBounds.min[2] << '\n'
           << "SPACING " << spacing << ' ' << spacing << ' ' << spacing << '\n'
           << "POINT_DATA " << pointCount << '\n';
    detail::vtkScalar<std::int32_t>(output,"fluid",pointCount,[&](std::size_t i) { return fluid[i]; });
    detail::vtkFinish(output,path);
}

//=============================================================================
// 3D cell scalar field
//=============================================================================

inline void writeVTK3DCellScalar(
    const std::filesystem::path& path,
    const std::size_t nx,
    const std::size_t ny,
    const std::size_t nz,
    const Point& origin,
    const double spacing,
    const std::vector<double>& scalar)
{
    if(!(spacing>0.0) || !std::isfinite(spacing))
        throw std::invalid_argument("Grid spacing must be positive and finite.");
    if(nx==0 || ny==0 || nz==0)
        throw std::invalid_argument("VTK cell dimensions must be positive.");
    detail::vtkCount(nx); detail::vtkCount(ny); detail::vtkCount(nz);
    detail::vtkCount(nx+1); detail::vtkCount(ny+1); detail::vtkCount(nz+1);
    detail::vtkCount(detail::vtkProduct(detail::vtkProduct(nx+1,ny+1),nz+1));
    const auto cellCount=detail::vtkCount(detail::vtkProduct(detail::vtkProduct(nx,ny),nz));
    if(scalar.size()!=cellCount)
        throw std::invalid_argument("3D scalar size does not match the VTK cell dimensions.");
    auto output=detail::vtkOpen(path);
    output << "# vtk DataFile Version 3.0\nnticLBM 3D cell scalar field\nBINARY\n"
           << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << nx+1 << ' ' << ny+1 << ' ' << nz+1 << '\n'
           << "ORIGIN " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << '\n'
           << "SPACING " << spacing << ' ' << spacing << ' ' << spacing << '\n'
           << "CELL_DATA " << cellCount << '\n';
    detail::vtkScalar<double>(output,"scalar",cellCount,[&](std::size_t i) { return scalar[i]; });
    detail::vtkFinish(output,path);
}

inline void writeVTK3DBoundingBox(
    const std::filesystem::path& path,
    const BoundingBox& box)
{
    constexpr std::array<unsigned,8> order{{0,1,3,2,4,5,7,6}};
    constexpr std::array<std::int32_t,36> lines{{
        2,0,1,2,1,2,2,2,3,2,3,0,2,4,5,2,5,6,2,6,7,2,7,4,2,0,4,2,1,5,2,2,6,2,3,7}};
    auto output=detail::vtkOpen(path);
    output << "# vtk DataFile Version 3.0\nnticLBM 3D bounding box\nBINARY\n"
           << "DATASET POLYDATA\nPOINTS 8 double\n";
    detail::vtkArray<double>(output,24,[&](std::size_t i) {
        const auto axis=i%3;
        return ((order[i/3]>>axis)&1U)!=0 ? box.max[axis] : box.min[axis];
    });
    output << "LINES 12 36\n";
    detail::vtkArray<std::int32_t>(output,lines.size(),[&](std::size_t i) { return lines[i]; });
    detail::vtkFinish(output,path);
}

//=============================================================================
// Boundary Q points
//=============================================================================

template<lattice::LatticeType Type>
void writeVTKBoundaryPoints(
    const std::filesystem::path& path, const GeometryAnalysis3D<Type>& analysis)
{
    using Lattice = typename GeometryAnalysis3D<Type>::Lattice;
    detail::writeVTKBoundaryPoints(path, analysis,
        [&](std::size_t entry) { return std::array<double,3>{
            analysis.boundaryX[entry], analysis.boundaryY[entry], analysis.boundaryZ[entry]}; },
        [](std::size_t link) { return std::array<int,3>{
            Lattice::ex[link], Lattice::ey[link], Lattice::ez[link]}; });
}

//=============================================================================
// Multilevel cell and Q output
//=============================================================================

template<lattice::LatticeType Type>
void writeVTKMultilevelCells(const std::filesystem::path& path,
    const MultilevelGrid<3,Type>& grid,bool auxiliary=false)
{ detail::writeVTKMultilevelCells(path,grid,auxiliary); }

template<lattice::LatticeType Type>
void writeVTKMultilevelBoundaryPoints(const std::filesystem::path& path,
    const MultilevelGrid<3,Type>& grid)
{ detail::writeVTKMultilevelBoundaryPoints(path,grid); }

} // namespace ntic::lbm::geometry::test
