#pragma once

#include "vtk_output_2d.hpp"


template<typename Grid,std::enable_if_t<Grid::dimension==3,int> =0>
void writeVTKMultilevelCells(const std::filesystem::path& path,
    const Grid& grid,bool auxiliary=false,
    const MultilevelCellDataWriter& writeAdditionalCellData={})
{ vtkWriteMultilevelCells(path,grid,auxiliary,writeAdditionalCellData); }

template<typename Grid,std::enable_if_t<Grid::dimension==3,int> =0>
void writeVTKMultilevelBoundaryPoints(const std::filesystem::path& path,
    const Grid& grid)
{ vtkWriteMultilevelBoundaryPoints(path,grid); }


template<typename Grid,std::enable_if_t<Grid::dimension==3,int> =0>
void writeVTKPartitionCells(const std::filesystem::path& path,
    const Grid& grid,const part::PartitionPlan& plan,
    const part::HardwareProfile& hardware,bool auxiliary=false)
{ vtkWritePartitionCells(path,grid,plan,hardware,auxiliary); }

template<typename Grid,std::enable_if_t<Grid::dimension==3,int> =0>
void writePartitionSummary(std::ostream& out,const Grid& grid,
    const part::PartitionPlan& plan,const part::PartitionOptions& options,const std::string& description,unsigned analysisThreads)
{ vtkWritePartitionSummary(out,grid,plan,options,description,analysisThreads); }

template<typename Grid,std::enable_if_t<Grid::dimension==3,int> =0>
void writePartitionOutputs(const std::filesystem::path& root,const std::string& caseName,
    const Grid& grid,const part::PartitionPlan& plan,
    const part::PartitionOptions& options,const std::string& description,unsigned analysisThreads)
{ vtkWritePartitionOutputs(root,caseName,grid,plan,options,description,analysisThreads); }
