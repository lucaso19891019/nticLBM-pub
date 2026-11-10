#include "multilevel/multilevel_grid.hpp"

#include <functional>
#include <iostream>
#include <stdexcept>

namespace geo=ntic::lbm::geometry;
namespace lat=ntic::lbm::lattice;

// Shared verification is implemented in multilevel_test_checks.cpp.
namespace ntic::lbm::geometry::test
{
void require(bool condition,const char* message);
geo::MultilevelGrid<2> makeRegionGrid(int threads);
void checkTopology(const geo::MultilevelGrid<2>& grid);
} // namespace ntic::lbm::geometry::test

namespace
{
using geo::test::require;

void determinism()
{
    const auto serial=geo::test::makeRegionGrid(1),parallel=geo::test::makeRegionGrid(4);
    geo::test::checkTopology(serial); geo::test::checkTopology(parallel);
    require(serial.cells().sourceNode==parallel.cells().sourceNode &&
            serial.cells().auxiliaryChild==parallel.cells().auxiliaryChild &&
            serial.cells().role==parallel.cells().role && serial.cells().cellType==parallel.cells().cellType &&
            serial.boundaries().offsets==parallel.boundaries().offsets &&
            serial.boundaries().q==parallel.boundaries().q && serial.boundaries().cells==parallel.boundaries().cells,
            "1/4-thread cell/Face/Q output differs");
    require(serial.nodes().size()==parallel.nodes().size(),"1/4-thread topology size differs");
    for(std::size_t i=0; i<serial.nodes().size(); ++i)
    {
        const auto& a=serial.nodes()[i]; const auto& b=parallel.nodes()[i];
        require(a.coordinate==b.coordinate && a.parent==b.parent && a.firstChild==b.firstChild &&
                a.cell==b.cell && a.level==b.level && a.active==b.active,"1/4-thread topology differs");
    }
    require(serial.interfaces().fineCells==parallel.interfaces().fineCells &&
            serial.interfaces().blocks.size()==parallel.interfaces().blocks.size(),"1/4-thread interface connectivity differs");
    for(std::size_t i=0; i<serial.interfaces().blocks.size(); ++i)
    {
        const auto& a=serial.interfaces().blocks[i]; const auto& b=parallel.interfaces().blocks[i];
        require(a.coarse==b.coarse && a.auxiliary.begin==b.auxiliary.begin && a.auxiliary.end==b.auxiliary.end &&
                a.fine.begin==b.fine.begin && a.fine.end==b.fine.end,"1/4-thread interface ranges differ");
    }
}

} // namespace

int main()
{
    try
    {
        omp_set_dynamic(0); omp_set_num_threads(4);
        determinism();
        std::cout << "[PASS] determinism\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "[FAIL] test_multilevel_determinism: " << error.what() << '\n';
        return 1;
    }
}
