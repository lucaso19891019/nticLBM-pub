#include "rectangle.hpp"
#include "vtk_output.hpp"
#include "vtk_output_2d.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace ntic::lbm::geometry;
using namespace ntic::lbm::geometry::test;

namespace
{

std::vector<double> buildScalar(
    const std::vector<CellType>& cellTypes)
{
    std::vector<double> scalar(
        cellTypes.size(),
        0.0);

    #pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                cellTypes.size());
        ++index)
    {
        const std::size_t cellIndex =
            static_cast<std::size_t>(
                index);

        switch(cellTypes[cellIndex])
        {
            case CellType::Dry:
                scalar[cellIndex] =
                    0.0;
                break;

            case CellType::Boundary:
                scalar[cellIndex] =
                    1.0;
                break;

            case CellType::Interior:
                scalar[cellIndex] =
                    2.0;
                break;
        }
    }

    return scalar;
}

void writeAnalysis(
    const Rectangle& geometry,
    const BoundingBox2D& domainBounds,
    const double spacing,
    const std::filesystem::path& directory)
{
    const std::size_t nx =
        static_cast<std::size_t>(
            std::ceil(
                domainBounds.width() /
                spacing));

    const std::size_t ny =
        static_cast<std::size_t>(
            std::ceil(
                domainBounds.height() /
                spacing));

    std::vector<CellType> cellTypes;
    std::vector<double> boundaryX;
    std::vector<double> boundaryY;

    geometry.interiorAreaAnalysis(
        spacing,
        cellTypes,
        boundaryX,
        boundaryY);

    const std::vector<double> scalar =
        buildScalar(
            cellTypes);

    const Point2D firstCellCenter{
        domainBounds.min[0] +
            0.5 * spacing,
        domainBounds.min[1] +
            0.5 * spacing
    };

    writeVTK2DScalar(
        directory /
            "geometry.vtk",
        nx,
        ny,
        firstCellCenter,
        spacing,
        scalar);

    writeVTK2DBoundingBox(
        directory /
            "bounding_box.vtk",
        geometry.boundingBox());
}

} // namespace

int main()
{
    try
    {
        std::cout
            << "========================================\n"
            << "Rectangle Geometry Test\n"
            << "========================================\n\n";

        const BoundingBox2D openBox{
            {-2.0, -1.0},
            {10.0, 9.0}
        };
