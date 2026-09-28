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
    const Rectangle& rectangle,
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

    rectangle.interiorAreaAnalysis(
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
        rectangle.boundingBox());
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

        const Rectangle internalRectangle{
            {1.0, 2.0},
            {7.0, 6.0},
            FlowType::Internal
        };

        const Rectangle externalRectangle{
            {1.0, 2.0},
            {7.0, 6.0},
            FlowType::External,
            &openBox
        };

        const BoundingBox2D boundingBox =
            internalRectangle.boundingBox();

        std::cout
            << "Rectangle:\n"
            << "  min = ("
            << internalRectangle.min[0] << ", "
            << internalRectangle.min[1] << ")\n"
            << "  max = ("
            << internalRectangle.max[0] << ", "
            << internalRectangle.max[1] << ")\n\n";

        std::cout
            << "Bounding box:\n"
            << "  min = ("
            << boundingBox.min[0] << ", "
            << boundingBox.min[1] << ")\n"
            << "  max = ("
            << boundingBox.max[0] << ", "
            << boundingBox.max[1] << ")\n\n";

        std::cout
            << "Open box:\n"
            << "  min = ("
            << openBox.min[0] << ", "
            << openBox.min[1] << ")\n"
            << "  max = ("
            << openBox.max[0] << ", "
            << openBox.max[1] << ")\n\n";

        double spacing =
            0.0;

        std::cout
            << "Enter grid spacing: ";

        std::cin
            >> spacing;

        if(!std::cin)
        {
            throw std::runtime_error(
                "Failed to read grid spacing.");
        }

        if(spacing <= 0.0)
        {
            throw std::runtime_error(
                "Grid spacing must be positive.");
        }

        const std::filesystem::path outputDirectory =
            std::filesystem::path(
                GEOMETRY_TEST_OUTPUT_DIR) /
            "test_rectangle_vtks";

        recreateOutputDirectory(
            outputDirectory);

        const std::filesystem::path internalDirectory =
            outputDirectory /
            "internal_vtks";

        const std::filesystem::path externalDirectory =
            outputDirectory /
            "external_vtks";

        std::filesystem::create_directories(
            internalDirectory);

        std::filesystem::create_directories(
            externalDirectory);

        writeAnalysis(
            internalRectangle,
            boundingBox,
            spacing,
            internalDirectory);

        writeAnalysis(
            externalRectangle,
            openBox,
            spacing,
            externalDirectory);

        writeVTK2DBoundingBox(
            externalDirectory /
                "open_box.vtk",
            openBox);

        std::cout
            << "\nVTK output written to:\n"
            << outputDirectory
            << "\n";

        return 0;
    }
    catch(const std::exception& exception)
    {
        std::cerr
            << "Error: "
            << exception.what()
            << "\n";

        return 1;
    }
}