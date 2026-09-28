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
    const Cylinder& cylinder,
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

    cylinder.interiorAreaAnalysis(
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
        cylinder.boundingBox());
}

} // namespace


int main()
{
    try
    {
        std::cout
            << "========================================\n"
            << "Cylinder Geometry Test\n"
            << "========================================\n\n";

        const BoundingBox openBox{
            {-1.0, -2.0, -1.0},
            {9.0, 8.0, 11.0}
        };

        const Cylinder internalCylinder{
            {4.0, 3.0, 5.0},
            2.0,
            6.0,
            Axis::Z,
            FlowType::Internal
        };

        const Cylinder externalCylinder{
            {4.0, 3.0, 5.0},
            2.0,
            6.0,
            Axis::Z,
            FlowType::External,
            &openBox
        };

        const BoundingBox boundingBox =
            internalCylinder.boundingBox();

        std::cout
            << "Cylinder:\n"
            << "  center = ("
            << internalCylinder.center[0] << ", "
            << internalCylinder.center[1] << ", "
            << internalCylinder.center[2] << ")\n"
            << "  radius = "
            << internalCylinder.radius
            << "\n"
            << "  length = "
            << internalCylinder.length
            << "\n"
            << "  axis = Z\n\n";

        std::cout
            << "Bounding box:\n"
            << "  min = ("
            << boundingBox.min[0] << ", "
            << boundingBox.min[1] << ", "
            << boundingBox.min[2] << ")\n"
            << "  max = ("
            << boundingBox.max[0] << ", "
            << boundingBox.max[1] << ", "
            << boundingBox.max[2] << ")\n\n";

        std::cout
            << "Open box:\n"
            << "  min = ("
            << openBox.min[0] << ", "
            << openBox.min[1] << ", "
            << openBox.min[2] << ")\n"
            << "  max = ("
            << openBox.max[0] << ", "
            << openBox.max[1] << ", "
            << openBox.max[2] << ")\n\n";

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
            "test_cylinder_vtks";

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
            internalCylinder,
            boundingBox,
            spacing,
            internalDirectory);

        writeAnalysis(
            externalCylinder,
            openBox,
            spacing,
            externalDirectory);

        writeVTK3DBoundingBox(
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