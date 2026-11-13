#pragma once

#include "vtk_output_2d.hpp"

namespace ntic::lbm::geometry::test
{

template<lattice::LatticeType Type>
void writeVTKMultilevelCells(const std::filesystem::path& path,
    const MultilevelGrid<3,Type>& grid,bool auxiliary=false,
    const MultilevelCellDataWriter& writeAdditionalCellData={})
{ detail::writeVTKMultilevelCells(path,grid,auxiliary,writeAdditionalCellData); }

template<lattice::LatticeType Type>
void writeVTKMultilevelBoundaryPoints(const std::filesystem::path& path,
    const MultilevelGrid<3,Type>& grid)
{ detail::writeVTKMultilevelBoundaryPoints(path,grid); }

} // namespace ntic::lbm::geometry::test

namespace ntic::lbm::mpi::test
{

template<lattice::LatticeType Type>
void writeVTKPartitionCells(const std::filesystem::path& path,
    const geometry::MultilevelGrid<3,Type>& grid,const PartitionPlan& plan,
    const HardwareProfile& hardware,bool auxiliary=false)
{ detail::writeVTKPartitionCells(path,grid,plan,hardware,auxiliary); }

template<lattice::LatticeType Type>
void writePartitionSummary(std::ostream& out,const geometry::MultilevelGrid<3,Type>& grid,
    const PartitionPlan& plan,const PartitionOptions& options,const std::string& description,unsigned analysisThreads)
{ detail::writePartitionSummary(out,grid,plan,options,description,analysisThreads); }

template<lattice::LatticeType Type>
void writePartitionOutputs(const std::filesystem::path& root,const std::string& caseName,
    const geometry::MultilevelGrid<3,Type>& grid,const PartitionPlan& plan,
    const PartitionOptions& options,const std::string& description,unsigned analysisThreads)
{ detail::writePartitionOutputs(root,caseName,grid,plan,options,description,analysisThreads); }

} // namespace ntic::lbm::mpi::test
