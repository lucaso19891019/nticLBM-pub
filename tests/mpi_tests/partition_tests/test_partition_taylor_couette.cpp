#include "partition.hpp"
#include "vtk_output_2d.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace geo=ntic::lbm::geometry;
namespace part=ntic::lbm::mpi;

int main()
{
    try
    {
        std::cout << "Taylor-Couette annulus: geometry and partition only\n"
                  << "Three levels (0,1,2), four equal-capacity ranks. No LBM evolution or MPI transport.\n"
                  << "Inner/outer radius ratio 9:10; the gap spans 100 finest cells.\n\n";
        double finestSpacing=0,bandCells=0;
        std::cout << "Finest grid spacing (recommended 0.01): "; std::cin >> finestSpacing;
        std::cout << "Narrow wall-band thickness in finest cells (recommended 4): "; std::cin >> bandCells;
        const double innerRadius=900*finestSpacing,outerRadius=1000*finestSpacing;
        const double coarseSpacing=4*finestSpacing,narrowBand=bandCells*finestSpacing,wideBand=2*narrowBand;
        if(!std::cin || !(finestSpacing>0) || !std::isfinite(outerRadius) ||
           !(bandCells>0) || !std::isfinite(bandCells) || !(narrowBand>0) || !std::isfinite(wideBand))
            throw std::invalid_argument("Spacing and wall-band cell count must be finite and positive.");
        if(4*bandCells>=100)
            std::cout << "Warning: the two wide bands span the gap; little or no coarse interior may remain.\n";

        const int threads=std::max(1,std::min(32,omp_get_max_threads()));
        omp_set_num_threads(threads);
        const auto geometry=geo::Circle({0,0},outerRadius)-geo::Circle({0,0},innerRadius);
        const auto walls=geo::fluidRegion(geometry);
        const std::vector<geo::RefinementRule<2>> rules{
            geo::RefinementRule<2>::surface(walls,1,wideBand),
            geo::RefinementRule<2>::surface(walls,2,narrowBand)};
        geo::MultilevelOptions gridOptions;
        gridOptions.maxLevel=2; gridOptions.analysisThreads=static_cast<unsigned>(threads);
        gridOptions.maxNodes=3000000; gridOptions.maxCells=2000000;
        std::cout << "Generating the concentric annulus (about 500 x 500 level-0 topology cells).\n"
                  << "Limits: " << gridOptions.maxNodes << " tree nodes, " << gridOptions.maxCells << " retained slots.\n"
                  << std::flush;
        geo::MultilevelGrid<2> grid(geometry.bounds,coarseSpacing,gridOptions);
        grid.generate(geometry,rules);

        part::PartitionOptions options;
        options.ranks=4; options.blockLevel=0; options.graphThreads=std::min(8,threads);
        options.hardware.computeRates.assign(4,1.0);
        // Node IDs remain unspecified; no physical node placement is implied.
        std::cout << "Partitioning with METIS...\n" << std::flush;
        const auto plan=part::partition(grid,options);
        std::ostringstream description;
        description << std::setprecision(17) << "Taylor-Couette concentric annulus, no additional refinement rectangles.\n"
                    << "Rin=" << innerRadius << ", Rout=" << outerRadius << ", radius/diameter ratio=9:10.\n"
                    << "hmin=" << finestSpacing << ", h0=" << coarseSpacing << ", gap=100*hmin.\n"
                    << "Level-1 wall band=" << wideBand << ", level-2 wall band=" << narrowBand
                    << " (" << bandCells << " finest cells).\n"
                    << "Refinement selectors are conservative; transition balancing can add nearby cells.";
        part::test::writePartitionOutputs(MPI_TEST_OUTPUT_DIR,"test_partition_taylor_couette_vtks",
            grid,plan,options,description.str(),gridOptions.analysisThreads);
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
