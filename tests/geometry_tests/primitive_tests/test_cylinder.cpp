#include "cylinder.hpp"
#include "vtk_output.hpp"
#include "vtk_output_3d.hpp"

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
    const Cylinder& geometry,
    const BoundingBox& domainBounds,
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

    const std::size_t nz =
        static_cast<std::size_t>(
            std::ceil(
                domainBounds.depth() /
                spacing));

    std::vector<CellType> cellTypes;
    std::vector<double> boundaryX;
    std::vector<double> boundaryY;
    std::vector<double> boundaryZ;

    geometry.interiorAreaAnalysis(
        spacing,
        cellTypes,
        boundaryX,
        boundaryY,
        boundaryZ);

    const std::vector<double> scalar =
        buildScalar(
            cellTypes);

    const Point firstCellCenter{
        domainBounds.min[0] +
            0.5 * spacing,
        domainBounds.min[1] +
            0.5 * spacing,
        domainBounds.min[2] +
            0.5 * spacing
    };

    writeVTK3DScalar(
        directory /
            "geometry.vtk",
        nx,
        ny,
        nz,
        firstCellCenter,
        spacing,
        scalar);

    writeVTK3DBoundingBox(
        directory /
            "bounding_box.vtk",
        geometry.boundingBox());
}

} // namespace

int main()
{
    try