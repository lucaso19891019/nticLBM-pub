#pragma once

#include "multilevel/multilevel_grid.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>

namespace ntic::lbm::geometry::test
{

template<std::size_t D,lattice::LatticeType Type>
void writeVTKMultilevelCells(const std::filesystem::path& path,
    const MultilevelGrid<D,Type>& grid,bool auxiliary=false)
{
    const auto& cells=grid.cells();
    std::vector<std::size_t> selected;
    for(std::size_t c=0; c<cells.size(); ++c) if(cells.isAuxiliary(c)==auxiliary) selected.push_back(c);
    constexpr std::size_t corners=std::size_t{1}<<D;
    // VTK quad/hexahedron ordering differs from bitwise child ordering.
    constexpr std::array<unsigned,8> order{0,1,3,2,4,5,7,6};
    std::ofstream out(path);
    if(!out) throw std::runtime_error("Cannot open " + path.string());
    out << std::setprecision(17) << "# vtk DataFile Version 3.0\nnticLBM multilevel cells\nASCII\n"
        << "DATASET UNSTRUCTURED_GRID\nPOINTS " << selected.size()*corners << " double\n";
    for(auto cell:selected)
    {
        const auto p=grid.coordinate(cell); const double h=grid.spacing(grid.level(cell));
        for(std::size_t corner=0; corner<corners; ++corner)
        {
            for(std::size_t d=0; d<3; ++d)
            {
                const double value=d<D ? grid.domain().min[d]+(p[d]+((order[corner]>>d)&1U))*h : 0;
                out << value << (d==2 ? '\n' : ' ');
            }
        }
    }
    out << "CELLS " << selected.size() << ' ' << selected.size()*(corners+1) << '\n';
    for(std::size_t cell=0; cell<selected.size(); ++cell)
    {
        out << corners; for(std::size_t corner=0; corner<corners; ++corner) out << ' ' << cell*corners+corner;
        out << '\n';
    }
    out << "CELL_TYPES " << selected.size() << '\n';
    for(std::size_t i=0; i<selected.size(); ++i) out << (D==2 ? 9 : 12) << '\n';
    out << "CELL_DATA " << selected.size() << "\nSCALARS scalar int 1\nLOOKUP_TABLE default\n";
    for(auto c:selected) out << static_cast<int>(cells.cellType[c]) << '\n';
    out << "SCALARS level int 1\nLOOKUP_TABLE default\n";
    for(auto c:selected) out << grid.level(c) << '\n';
    out << "SCALARS role int 1\nLOOKUP_TABLE default\n";
    for(auto c:selected) out << static_cast<unsigned>(cells.role[c]) << '\n';
    out << "SCALARS CellIndex unsigned_long 1\nLOOKUP_TABLE default\n";
    for(auto c:selected) out << c << '\n';
    out << "SCALARS SourceNode unsigned_long 1\nLOOKUP_TABLE default\n";
    for(auto c:selected) out << cells.sourceNode[c] << '\n';
    if(!out) throw std::runtime_error("Failed to write " + path.string());
}

template<std::size_t D,lattice::LatticeType Type>
void writeVTKMultilevelBoundaryPoints(const std::filesystem::path& path,
    const MultilevelGrid<D,Type>& grid)
{
    using Lattice=lattice::Model<Type>;
    constexpr auto links=Lattice::nStencils-1;
    const auto& boundary=grid.boundaries();
    struct Hit { std::size_t entry,face,direction; };
    std::vector<Hit> hits;
    for(std::size_t face=0; face<boundary.nBoundaries; ++face)
        for(auto e=boundary.offsets[face]; e<boundary.offsets[face+1]; ++e)
            for(std::size_t i=1; i<=links; ++i) if(boundary.q[e*links+i-1]>0) hits.push_back({e,face,i});
    std::ofstream out(path);
    if(!out) throw std::runtime_error("Cannot open " + path.string());
    out << std::setprecision(17) << "# vtk DataFile Version 3.0\nMultilevel Q intersections\nASCII\n"
        << "DATASET POLYDATA\nPOINTS " << hits.size() << " double\n";
    for(const auto& hit:hits)
    {
        const auto c=boundary.cells[hit.entry]; auto p=grid.center(c);
        const double distance=boundary.q[hit.entry*links+hit.direction-1]*grid.spacing(grid.level(c));
        p[0]+=distance*Lattice::ex[hit.direction]; p[1]+=distance*Lattice::ey[hit.direction];
        if constexpr(D==3) p[2]+=distance*Lattice::ez[hit.direction];
        out << p[0] << ' ' << p[1] << ' ';
        if constexpr(D==3) out << p[2]; else out << 0;
        out << '\n';
    }
    out << "VERTICES " << hits.size() << ' ' << 2*hits.size() << '\n';
    for(std::size_t i=0; i<hits.size(); ++i) out << "1 " << i << '\n';
    out << "POINT_DATA " << hits.size() << "\nSCALARS Q double 1\nLOOKUP_TABLE default\n";
    for(const auto& hit:hits) out << boundary.q[hit.entry*links+hit.direction-1] << '\n';
    out << "SCALARS FaceID int 1\nLOOKUP_TABLE default\n";
    for(const auto& hit:hits) out << hit.face << '\n';
    out << "SCALARS LinkID int 1\nLOOKUP_TABLE default\n";
    for(const auto& hit:hits) out << hit.direction << '\n';
    out << "SCALARS CellIndex unsigned_long 1\nLOOKUP_TABLE default\n";
    for(const auto& hit:hits) out << boundary.cells[hit.entry] << '\n';
    out << "SCALARS level int 1\nLOOKUP_TABLE default\n";
    for(const auto& hit:hits) out << grid.level(boundary.cells[hit.entry]) << '\n';
    out << "SCALARS role int 1\nLOOKUP_TABLE default\n";
    for(const auto& hit:hits) out << static_cast<unsigned>(grid.cells().role[boundary.cells[hit.entry]]) << '\n';
    if(!out) throw std::runtime_error("Failed to write " + path.string());
}

} // namespace ntic::lbm::geometry::test
