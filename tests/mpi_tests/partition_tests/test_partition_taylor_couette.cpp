#include "partition.hpp"
#include "partition_mpi.hpp"
#include "vtk_output_2d.hpp"

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <optional>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace geo=ntic::lbm::geometry;
namespace part=ntic::lbm::mpi;

int main(int argc,char** argv)
{
    int provided=0;
    if(MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided)!=MPI_SUCCESS) return 1;
    int mainThread=0;MPI_Is_thread_main(&mainThread);
    if(provided<MPI_THREAD_FUNNELED || !mainThread) MPI_Abort(MPI_COMM_WORLD,1);
    int rank=0,size=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);int result=0;
    try
    {
        using Grid=geo::MultilevelGrid<2>;
        using Links=geo::NeighborLinks<2,Grid::latticeType>;
        using Builder=part::LocalLayoutBuilder<2,Grid::latticeType>;
        MPI_Comm shared=MPI_COMM_NULL;MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&shared);
        int localRanks=1;MPI_Comm_size(shared,&localRanks);MPI_Comm_free(&shared);
        const int localThreads=std::getenv("OMP_NUM_THREADS") ? std::max(1,omp_get_max_threads()) :
            std::max(1,omp_get_num_procs()/localRanks);
        const int globalThreads=std::max(1,std::min(32,omp_get_num_procs()));
        omp_set_dynamic(0);
        std::unique_ptr<Grid> globalGrid;std::optional<Links> globalLinks;std::optional<part::PartitionPlan> globalPlan;
        std::unique_ptr<Builder> builder;part::PartitionOptions options;
        std::string rootError,rootDescription;unsigned rootAnalysisThreads=1;
        if(rank==0)
        {
            try
            {
                if(size!=4) throw std::invalid_argument("This manual case requires mpiexec -n 4.");
                omp_set_num_threads(globalThreads);

                std::cout << "Taylor-Couette annulus: geometry and partition only\n"
                          << "Three levels (0,1,2), four equal-capacity ranks. No LBM evolution or PDF transport.\n"
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

                const int threads=globalThreads;
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
                globalGrid=std::make_unique<Grid>(geometry.bounds,coarseSpacing,gridOptions);
                auto& grid=*globalGrid;
                grid.generate(geometry,rules);
                globalLinks=geo::linkNeighbors(grid);


                options.ranks=4; options.blockLevel=0; options.graphThreads=std::min(8,threads);
                options.hardware.computeRates.assign(4,1.0);
                // Node IDs remain unspecified; no physical node placement is implied.
                std::cout << "Partitioning with METIS...\n" << std::flush;
                globalPlan=part::partition(grid,options);
                const auto& plan=*globalPlan;
                std::ostringstream description;
                description << std::setprecision(17) << "Taylor-Couette concentric annulus, no additional refinement rectangles.\n"
                            << "Rin=" << innerRadius << ", Rout=" << outerRadius << ", radius/diameter ratio=9:10.\n"
                            << "hmin=" << finestSpacing << ", h0=" << coarseSpacing << ", gap=100*hmin.\n"
                            << "Level-1 wall band=" << wideBand << ", level-2 wall band=" << narrowBand
                            << " (" << bandCells << " finest cells).\n"
                            << "Refinement selectors are conservative; transition balancing can add nearby cells.";
                rootDescription=description.str(); rootAnalysisThreads=gridOptions.analysisThreads;
                builder=std::make_unique<Builder>(grid,*globalLinks,plan);

            }
            catch(const std::exception& e) {rootError=e.what();}
            catch(...) {rootError="Unknown root geometry/partition error.";}
        }
        omp_set_num_threads(localThreads);
        const auto input=part::distributeLayoutInput<2,Grid::latticeType>(builder.get(),MPI_COMM_WORLD,0,rootError);
        const auto layout=part::initializeLocalLayout(input,MPI_COMM_WORLD);
        std::array<std::uint64_t,10> summary{};
        summary[0]=layout.cells.size();summary[1]=layout.sourceIDs.size();
        for(const auto& cell:layout.cells) {summary[2]+=cell.ghost;summary[3]+=cell.halo;}
        summary[4]=layout.candidates;summary[5]=layout.selectedCandidate;summary[6]=layout.score.communicationRanges;
        summary[7]=layout.score.interfaceRanges;summary[8]=layout.score.boundaryRanges;summary[9]=layout.score.anchorDistance;
        std::vector<std::uint64_t> summaries;part::detail::LayoutError summaryError;
        try {if(rank==0) summaries.resize(static_cast<std::size_t>(size)*summary.size());}
        catch(const std::exception& e) {summaryError=e.what();}
        part::detail::layoutErrorGate(MPI_COMM_WORLD,summaryError);
        part::detail::layoutMpiCheck(MPI_Gather(summary.data(),static_cast<int>(summary.size()),MPI_UINT64_T,
            rank==0 ? summaries.data() : nullptr,static_cast<int>(summary.size()),MPI_UINT64_T,0,MPI_COMM_WORLD),MPI_COMM_WORLD);
        if(rank==0)
        {
            std::cout << "Actual MPI ranks=" << size << "; rank-local layouts were built on their respective ranks.\n"
                      << "Node IDs in the cost model are modeled placement, not physical host detection.\n";
            for(int peer=0;peer<size;++peer)
            {
                const auto offset=static_cast<std::size_t>(peer)*summary.size();
                std::cout << "Local layout rank " << peer << ": cells=" << summaries[offset] << ", sources=" << summaries[offset+1]
                          << ", Ghost=" << summaries[offset+2] << ", Halo=" << summaries[offset+3]
                          << ", candidates=" << summaries[offset+4] << ", selected=" << summaries[offset+5]
                          << ", score=(" << summaries[offset+6] << ',' << summaries[offset+7] << ','
                          << summaries[offset+8] << ',' << summaries[offset+9] << ")\n";
            }
            writePartitionOutputs(MPI_TEST_OUTPUT_DIR,"test_partition_taylor_couette_vtks",*globalGrid,*globalPlan,options,rootDescription,rootAnalysisThreads);
        }
    }
    catch(const std::exception& e) {std::cerr << "Error on rank " << rank << ": " << e.what() << '\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
    MPI_Allreduce(MPI_IN_PLACE,&result,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);MPI_Finalize();return result;
}
