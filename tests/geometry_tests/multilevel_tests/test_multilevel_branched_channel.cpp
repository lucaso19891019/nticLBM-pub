#include "multilevel_grid.hpp"
#include "vtk_output_3d.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"

#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace geo=ntic::lbm::geometry;

int main()
{
    try
    {
        std::cout << "Multilevel STL branched channel\n\nInternal branched-channel STL. Refine the fluid-side band of Face 0 (Other Surface).\n\n";
        double spacing=0.0,thickness=0.0;
        unsigned maxLevel=0;
        std::cout << "Base grid spacing: "; std::cin >> spacing;
        std::cout << "Maximum level (1..4): "; std::cin >> maxLevel;
        std::cout << "Fluid-side wall refinement thickness: "; std::cin >> thickness;
        if(!std::cin || !(spacing>0.0) || !std::isfinite(spacing) || maxLevel<1 || maxLevel>4 ||
           !(thickness>=0.0) || !std::isfinite(thickness))
            throw std::invalid_argument("Invalid manual test parameters.");

        const auto file=std::filesystem::path(GEOMETRY_TEST_SOURCE_DIR) /
            "stl_tests" / "stls" / "smooth_irregular_branched_channel.stl";
        std::cout << "STL file: " << file << '\n';
        auto data=ntic::lbm::stl::read(file.string());
        ntic::lbm::stl::FacetTopology topology;
        ntic::lbm::stl::validate(data,"full",topology);
        const geo::STLGeometry geometry(std::move(topology),geo::FlowType::Internal);
        const std::vector<geo::RefinementRule<3>> rules{
            geo::RefinementRule<3>::surface(geo::fluidRegion(geometry),maxLevel,thickness,{0})};
        geo::MultilevelOptions options;
        options.maxLevel=maxLevel;
        geo::MultilevelGrid<3> grid(geometry.bounds,spacing,options);
        grid.generate(geometry,rules);

        const auto directory=std::filesystem::path(GEOMETRY_TEST_OUTPUT_DIR) /
            "test_multilevel_branched_channel_vtks";
        geo::test::recreateOutputDirectory(directory);
        geo::test::writeVTKMultilevelCells(directory / "geometry.vtk",grid);
        geo::test::writeVTKMultilevelCells(directory / "auxiliary.vtk",grid,true);
        geo::test::writeVTKMultilevelBoundaryPoints(directory / "boundary_q.vtk",grid);

        for(unsigned level=0; level<grid.levels().size(); ++level)
        {
            const auto& ranges=grid.levels()[level];
            std::size_t computes=0;
            for(auto cell=ranges.all.begin; cell<ranges.all.end; ++cell)
                computes+=grid.cells().computes(cell);
            std::cout << "Level " << level << ": h=" << grid.spacing(level) << ", "
                      << ranges.all.size() << " slots, " << computes << " real wet cells, "
                      << ranges.roles[geo::Auxiliary].size() << " auxiliary slots\n";
        }
        std::cout << "Tree nodes: " << grid.nodes().size()
                  << "\nInterface blocks: " << grid.interfaces().blocks.size()
                  << "\nVTK output written to:\n" << directory << '\n';
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
