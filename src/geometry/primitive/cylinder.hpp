#pragma once

#include "axis.hpp"
#include "bounding_box.hpp"
#include "cell_type.hpp"
#include "flow_type.hpp"
#include "geometry_analysis.hpp"
#include "point.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace ntic::lbm::geometry
{

struct Cylinder
{
    Point center;
    double radius;
    double length;
    Axis axis;

    BoundingBox bounds;
    BoundingBox openBox;

    FlowType flowType;


    Cylinder(
        const Point& inputCenter,
        const double inputRadius,
        const double inputLength,
        const Axis inputAxis,
        const FlowType inputFlowType = FlowType::Internal,
        const BoundingBox* inputOpenBox = nullptr)
        : center(inputCenter),
          radius(inputRadius),
          length(inputLength),
          axis(inputAxis),
          flowType(inputFlowType)
    {
        if(radius <= 0.0)
        {
            throw std::invalid_argument(
                "Cylinder radius must be positive.");
        }

        if(length <= 0.0)
        {
            throw std::invalid_argument(
                "Cylinder length must be positive.");
        }

        bounds =
            makeBounds();

        if(flowType == FlowType::External)
        {
            if(inputOpenBox != nullptr)
            {
                validateOpenBox(
                    *inputOpenBox);

                openBox =
                    *inputOpenBox;
            }
            else
            {
                openBox =
                    readOpenBox();
            }
        }
    }


    [[nodiscard]]
    BoundingBox boundingBox() const noexcept
    {
        return bounds;
    }


    void translate(
        const Point* targetPoint = nullptr)
    {
        Point displacement{};

        if(flowType == FlowType::Internal)
        {
            if(targetPoint == nullptr)
            {
                throw std::invalid_argument(
                    "Internal Cylinder translation requires a target point.");
            }

            displacement = {
                (*targetPoint)[0] - bounds.min[0],
                (*targetPoint)[1] - bounds.min[1],
                (*targetPoint)[2] - bounds.min[2]
            };
        }
        else if(flowType == FlowType::External)
        {
            if(targetPoint != nullptr)
            {
                throw std::invalid_argument(
                    "External Cylinder translation does not accept a target point.");
            }

            displacement = {
                -openBox.min[0],
                -openBox.min[1],
                -openBox.min[2]
            };
        }
        else
        {
            throw std::invalid_argument(
                "Unsupported Cylinder flow type.");
        }

        for(std::size_t d = 0;
            d < 3;
            ++d)
        {
            center[d] +=
                displacement[d];

            bounds.min[d] +=
                displacement[d];

            bounds.max[d] +=
                displacement[d];

            if(flowType ==
               FlowType::External)
            {
                openBox.min[d] +=
                    displacement[d];

                openBox.max[d] +=
                    displacement[d];
            }
        }
    }


    template <lattice::LatticeType Type>
    void analysis(GeometryAnalysis3D<Type>& analysis) const
    {
        interiorAreaAnalysis(analysis);
        boundaryQAnalysis(analysis);
    }


    template <lattice::LatticeType Type>
    void boundaryQAnalysis(
        GeometryAnalysis3D<Type>& analysis) const
    {
        using Lattice =
            typename GeometryAnalysis3D<Type>::Lattice;
    
        constexpr std::size_t nLinks =
            Lattice::nStencils - 1;
    
        const std::size_t boundaryCount =
            analysis.boundaryX.size();
    
        if(analysis.boundaryY.size() != boundaryCount ||
           analysis.boundaryZ.size() != boundaryCount)
        {
            throw std::invalid_argument(
                "Boundary coordinate sizes do not match.");
        }
    
        const double h = analysis.gridSpacing;
    
        if(!std::isfinite(h) || h <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be finite and positive.");
        }
    
        analysis.q.assign(
            boundaryCount * nLinks, -1.0);
    
        // Identify axial and radial coordinate dimensions.
        std::size_t axialDimension = 2;
        std::size_t radialDimension0 = 0;
        std::size_t radialDimension1 = 1;
    
        if(axis == Axis::X)
        {
            axialDimension = 0;
            radialDimension0 = 1;
            radialDimension1 = 2;
        }
        else if(axis == Axis::Y)
        {
            axialDimension = 1;
            radialDimension0 = 0;
            radialDimension1 = 2;
        }
        else if(axis != Axis::Z)
        {
            throw std::invalid_argument(
                "Unsupported Cylinder axis.");
        }
    
        const double halfLength = 0.5 * length;
        const double radiusSquared = radius * radius;
    
        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index < static_cast<std::ptrdiff_t>(boundaryCount);
            ++index)
        {
            const std::size_t b =
                static_cast<std::size_t>(index);
    
            // Position relative to the Cylinder center.
            const double p[3] = {
                analysis.boundaryX[b] - center[0],
                analysis.boundaryY[b] - center[1],
                analysis.boundaryZ[b] - center[2]
            };
    
            const double pa = p[axialDimension];
            const double pr0 = p[radialDimension0];
            const double pr1 = p[radialDimension1];
    
            for(std::size_t i = 1;
                i < Lattice::nStencils;
                ++i)
            {
                const double d[3] = {
                    h * static_cast<double>(Lattice::ex[i]),
                    h * static_cast<double>(Lattice::ey[i]),
                    h * static_cast<double>(Lattice::ez[i])
                };
    
                const double da = d[axialDimension];
                const double dr0 = d[radialDimension0];
                const double dr1 = d[radialDimension1];
    
                double qValue = -1.0;
    
                // Part 1: intersection with the lateral surface.
                //
                // (pr0 + q*dr0)^2 +
                // (pr1 + q*dr1)^2 = radius^2
    
                const double A =
                    dr0 * dr0 + dr1 * dr1;
    
                const double B =
                    2.0 * (pr0 * dr0 + pr1 * dr1);
    
                const double C =
                    pr0 * pr0 + pr1 * pr1 -
                    radiusSquared;
    
                if(A > 0.0)
                {
                    const double discriminant =
                        B * B - 4.0 * A * C;
    
                    if(discriminant >= 0.0)
                    {
                        const double sqrtD =
                            std::sqrt(discriminant);
    
                        const double rootTerm =
                            -0.5 *
                            (B + std::copysign(sqrtD, B));
    
                        double root1;
                        double root2;
    
                        if(rootTerm == 0.0)
                        {
                            root1 = -B / (2.0 * A);
                            root2 = root1;
                        }
                        else
                        {
                            root1 = rootTerm / A;
                            root2 = C / rootTerm;
                        }
    
                        // A lateral intersection must lie
                        // between the two end caps.
                        if(root1 > 0.0 && root1 <= 1.0)
                        {
                            const double axial =
                                pa + root1 * da;
    
                            if(std::abs(axial) <= halfLength)
                            {
                                qValue = root1;
                            }
                        }
    
                        if(root2 > 0.0 && root2 <= 1.0)
                        {
                            const double axial =
                                pa + root2 * da;
    
                            if(std::abs(axial) <= halfLength &&
                               (qValue < 0.0 ||
                                root2 < qValue))
                            {
                                qValue = root2;
                            }
                        }
                    }
                }
    
                // Part 2: intersection with the two end caps.
                //
                // pa + q*da = +/- halfLength
    
                if(da != 0.0)
                {
                    for(int sign = -1; sign <= 1; sign += 2)
                    {
                        const double capPosition =
                            static_cast<double>(sign) *
                            halfLength;
    
                        const double root =
                            (capPosition - pa) / da;
    
                        if(root <= 0.0 || root > 1.0)
                        {
                            continue;
                        }
    
                        const double radial0 =
                            pr0 + root * dr0;
    
                        const double radial1 =
                            pr1 + root * dr1;
    
                        const double radialSquared =
                            radial0 * radial0 +
                            radial1 * radial1;
    
                        if(radialSquared <= radiusSquared &&
                           (qValue < 0.0 ||
                            root < qValue))
                        {
                            qValue = root;
                        }
                    }
                }
    
                analysis.q[b * nLinks + (i - 1)] =
                    qValue;
            }
        }
    }


    template <lattice::LatticeType Type>
    void interiorAreaAnalysis(GeometryAnalysis3D<Type>& analysis) const
    {
        const double h = analysis.gridSpacing;
        if(!std::isfinite(h) || h <= 0.0)
        {
            throw std::invalid_argument("Grid spacing must be finite and positive.");
        }

        analysis.domain = flowType == FlowType::Internal ? bounds : openBox;
        const BoundingBox& domain = analysis.domain;
        analysis.nx = static_cast<std::size_t>(std::ceil(domain.width() / h));
        analysis.ny = static_cast<std::size_t>(std::ceil(domain.height() / h));
        analysis.nz = static_cast<std::size_t>(std::ceil(domain.depth() / h));

        const std::size_t nx = analysis.nx;
        const std::size_t ny = analysis.ny;
        const std::size_t nz = analysis.nz;
        const std::size_t cellCount = nx * ny * nz;
        std::vector<CellType> cellTypes(cellCount, CellType::Dry);

        std::vector<std::size_t> cellBoundaryIDs(cellCount, 0);

        // StencilCell half-width is h; its half-diagonal is sqrt(3)*h.
        const double stencilHalfDiagonal = std::sqrt(3.0) * h;

        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index < static_cast<std::ptrdiff_t>(cellCount); ++index)
        {
            const std::size_t id = static_cast<std::size_t>(index);
            const std::size_t i = id % nx;
            const std::size_t j = (id / nx) % ny;
            const std::size_t k = id / (nx * ny);
            const Point centerPoint{
                domain.min[0] + (static_cast<double>(i) + 0.5) * h,
                domain.min[1] + (static_cast<double>(j) + 0.5) * h,
                domain.min[2] + (static_cast<double>(k) + 0.5) * h
            };

            const bool inside = pointInside(centerPoint);
            const bool wet = flowType == FlowType::Internal ? inside : !inside;
            const double surfaceDistance = distanceToBoundary(centerPoint);

            // The center must be strictly wet and must not lie on the surface.
            if(!wet || surfaceDistance == 0.0)
            {
                cellTypes[id] = CellType::Dry;
                continue;
            }
            if(surfaceDistance > stencilHalfDiagonal)
            {
                cellTypes[id] = CellType::Interior;
                continue;
            }

            const Point stencilMin{
                centerPoint[0] - h, centerPoint[1] - h, centerPoint[2] - h
            };
            const Point stencilMax{
                centerPoint[0] + h, centerPoint[1] + h, centerPoint[2] + h
            };

            bool intersects = false;
            std::size_t boundaryID = 0;
            
            boundaryIntersection(
                stencilMin,
                stencilMax,
                centerPoint,
                intersects,
                boundaryID);
            
            cellTypes[id] =
                intersects ? CellType::Boundary : CellType::Interior;
            
            if(intersects)
            {
                cellBoundaryIDs[id] = boundaryID;
            }
        }

        // Compact Boundary coordinates and geometric element IDs.
        // Fixed-size blocks preserve deterministic cell-index order.
        constexpr std::size_t blockSize = 4096;
        
        const std::size_t blockCount =
            cellCount / blockSize +
            (cellCount % blockSize != 0 ? 1 : 0);
        
        std::vector<std::size_t> blockOffsets(blockCount + 1, 0);
        
        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t block = 0;
            block < static_cast<std::ptrdiff_t>(blockCount);
            ++block)
        {
            const std::size_t b = static_cast<std::size_t>(block);
            const std::size_t begin = b * blockSize;
            const std::size_t end =
                std::min(begin + blockSize, cellCount);
        
            std::size_t count = 0;
        
            for(std::size_t id = begin; id < end; ++id)
            {
                if(cellTypes[id] == CellType::Boundary)
                {
                    ++count;
                }
            }
        
            blockOffsets[b + 1] = count;
        }
        
        for(std::size_t b = 0; b < blockCount; ++b)
        {
            blockOffsets[b + 1] += blockOffsets[b];
        }
        
        const std::size_t boundaryCount =
            blockOffsets[blockCount];
        
        analysis.boundaryX.resize(boundaryCount);
        analysis.boundaryY.resize(boundaryCount);
        analysis.boundaryZ.resize(boundaryCount);
        
        std::vector<std::size_t> boundaryIDs(boundaryCount, 0);
        
        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t block = 0;
            block < static_cast<std::ptrdiff_t>(blockCount);
            ++block)
        {
            const std::size_t b = static_cast<std::size_t>(block);
            const std::size_t begin = b * blockSize;
            const std::size_t end =
                std::min(begin + blockSize, cellCount);
        
            std::size_t slot = blockOffsets[b];
        
            for(std::size_t id = begin; id < end; ++id)
            {
                if(cellTypes[id] != CellType::Boundary)
                {
                    continue;
                }
        
                const std::size_t i = id % nx;
                const std::size_t j = (id / nx) % ny;
                const std::size_t k = id / (nx * ny);
        
                analysis.boundaryX[slot] =
                    domain.min[0] +
                    (static_cast<double>(i) + 0.5) * h;
        
                analysis.boundaryY[slot] =
                    domain.min[1] +
                    (static_cast<double>(j) + 0.5) * h;
        
                analysis.boundaryZ[slot] =
                    domain.min[2] +
                    (static_cast<double>(k) + 0.5) * h;
        
                boundaryIDs[slot] = cellBoundaryIDs[id];
        
                ++slot;
            }
        }
        
        // Fixed geometric element count:
        analysis.nBoundaries = 5; 
        
        analysis.boundaryOffsets.assign(
            analysis.nBoundaries + 1, 0);
        
        // Histogram.
        for(const std::size_t id : boundaryIDs)
        {
            ++analysis.boundaryOffsets[id + 1];
        }
        
        // Prefix sum.
        for(std::size_t k = 0; k < analysis.nBoundaries; ++k)
        {
            analysis.boundaryOffsets[k + 1] +=
                analysis.boundaryOffsets[k];
        }
        
        // Reorder coordinates by boundary element ID.
        // Within each group, preserve the original cell-index order.
        std::vector<std::size_t> next = analysis.boundaryOffsets;
        
        std::vector<double> groupedX(boundaryCount);
        std::vector<double> groupedY(boundaryCount);
        std::vector<double> groupedZ(boundaryCount);
        
        for(std::size_t b = 0; b < boundaryCount; ++b)
        {
            const std::size_t slot = next[boundaryIDs[b]]++;
        
            groupedX[slot] = analysis.boundaryX[b];
            groupedY[slot] = analysis.boundaryY[b];
            groupedZ[slot] = analysis.boundaryZ[b];
        }
        
        analysis.boundaryX.swap(groupedX);
        analysis.boundaryY.swap(groupedY);
        analysis.boundaryZ.swap(groupedZ);

        analysis.scalar.resize(cellCount);
        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index < static_cast<std::ptrdiff_t>(cellCount); ++index)
        {
            const std::size_t id = static_cast<std::size_t>(index);
            analysis.scalar[id] = static_cast<double>(static_cast<int>(cellTypes[id]));
        }
    }


private:


    [[nodiscard]]
    BoundingBox makeBounds() const noexcept
    {
        const double halfLength =
            0.5 * length;

        switch(axis)
        {
            case Axis::X:
                return {
                    {
                        center[0] - halfLength,
                        center[1] - radius,
                        center[2] - radius
                    },
                    {
                        center[0] + halfLength,
                        center[1] + radius,
                        center[2] + radius
                    }
                };

            case Axis::Y:
                return {
                    {
                        center[0] - radius,
                        center[1] - halfLength,
                        center[2] - radius
                    },
                    {
                        center[0] + radius,
                        center[1] + halfLength,
                        center[2] + radius
                    }
                };

            case Axis::Z:
                return {
                    {
                        center[0] - radius,
                        center[1] - radius,
                        center[2] - halfLength
                    },
                    {
                        center[0] + radius,
                        center[1] + radius,
                        center[2] + halfLength
                    }
                };
        }

        return {};
    }


    void coordinates(
        const Point& point,
        double& axial,
        double& radial0,
        double& radial1) const noexcept
    {
        switch(axis)
        {
            case Axis::X:
                axial =
                    point[0] - center[0];

                radial0 =
                    point[1] - center[1];

                radial1 =
                    point[2] - center[2];

                return;

            case Axis::Y:
                axial =
                    point[1] - center[1];

                radial0 =
                    point[0] - center[0];

                radial1 =
                    point[2] - center[2];

                return;

            case Axis::Z:
                axial =
                    point[2] - center[2];

                radial0 =
                    point[0] - center[0];

                radial1 =
                    point[1] - center[1];

                return;
        }
    }


    [[nodiscard]]
    bool pointInside(
        const Point& point) const noexcept
    {
        double axial = 0.0;
        double radial0 = 0.0;
        double radial1 = 0.0;

        coordinates(
            point,
            axial,
            radial0,
            radial1);

        return
            std::abs(axial) <
                0.5 * length &&
            radial0 * radial0 +
            radial1 * radial1 <
                radius * radius;
    }


    [[nodiscard]]
    double distanceToBoundary(
        const Point& point) const noexcept
    {
        double axial = 0.0;
        double radial0 = 0.0;
        double radial1 = 0.0;

        coordinates(
            point,
            axial,
            radial0,
            radial1);

        const double radialDistance =
            std::sqrt(
                radial0 * radial0 +
                radial1 * radial1);

        const double qRadial =
            radialDistance - radius;

        const double qAxial =
            std::abs(axial) -
            0.5 * length;

        const double outsideRadial =
            std::max(
                qRadial,
                0.0);

        const double outsideAxial =
            std::max(
                qAxial,
                0.0);

        const double outsideDistance =
            std::sqrt(
                outsideRadial *
                    outsideRadial +
                outsideAxial *
                    outsideAxial);

        const double insideDistance =
            std::min(
                std::max(
                    qRadial,
                    qAxial),
                0.0);

        return
            std::abs(
                outsideDistance +
                insideDistance);
    }


    // Boundary element indices:
    // 0: Lateral cylindrical surface
    // 1: Negative-axis end cap
    // 2: Positive-axis end cap
    // 3: Negative-axis circular rim
    // 4: Positive-axis circular rim
    //
    // Priority: circular rim (1D) before surface (2D).
    // Same dimension: nearest finite feature, then lowest ID.
    void boundaryIntersection(
        const Point& cellMin,
        const Point& cellMax,
        const Point& point,
        bool& intersects,
        std::size_t& boundaryID) const noexcept
    {
        std::size_t a = 2, u = 0, v = 1;
    
        if(axis == Axis::X)
        {
            a = 0; u = 1; v = 2;
        }
        else if(axis == Axis::Y)
        {
            a = 1; u = 0; v = 2;
        }
    
        const double axialMin = center[a] - 0.5 * length;
        const double axialMax = center[a] + 0.5 * length;
        const double radiusSquared = radius * radius;
    
        const double closestU =
            std::clamp(center[u], cellMin[u], cellMax[u]);
        const double closestV =
            std::clamp(center[v], cellMin[v], cellMax[v]);
    
        const double du = closestU - center[u];
        const double dv = closestV - center[v];
    
        const double minRadialSquared = du * du + dv * dv;
    
        const double farU = std::max(
            std::abs(cellMin[u] - center[u]),
            std::abs(cellMax[u] - center[u]));
    
        const double farV = std::max(
            std::abs(cellMin[v] - center[v]),
            std::abs(cellMax[v] - center[v]));
    
        const double maxRadialSquared =
            farU * farU + farV * farV;
    
        const bool axialOverlap =
            cellMax[a] >= axialMin &&
            cellMin[a] <= axialMax;
    
        const bool diskOverlap =
            minRadialSquared <= radiusSquared;
    
        const bool circleOverlap =
            minRadialSquared <= radiusSquared &&
            maxRadialSquared >= radiusSquared;
    
        const bool negativePlane =
            cellMin[a] <= axialMin &&
            cellMax[a] >= axialMin;
    
        const bool positivePlane =
            cellMin[a] <= axialMax &&
            cellMax[a] >= axialMax;
    
        const double radialDistance = std::hypot(
            point[u] - center[u],
            point[v] - center[v]);
    
        const double radialDelta = radialDistance - radius;
    
        intersects = false;
        boundaryID = 0;
    
        int bestDimension = 3;
        double bestDistance =
            std::numeric_limits<double>::infinity();
    
        const auto consider = [&](std::size_t id,
                                  int dimension,
                                  bool hit,
                                  double distance)
        {
            if(!hit)
            {
                return;
            }
    
            intersects = true;
    
            if(dimension < bestDimension ||
               (dimension == bestDimension &&
                (distance < bestDistance ||
                 (distance == bestDistance && id < boundaryID))))
            {
                bestDimension = dimension;
                bestDistance = distance;
                boundaryID = id;
            }
        };
    
        const double axialClamped =
            std::clamp(point[a], axialMin, axialMax);
    
        const double lateralAxialDelta =
            point[a] - axialClamped;
    
        const double lateralDistanceSquared =
            radialDelta * radialDelta +
            lateralAxialDelta * lateralAxialDelta;
    
        const double outsideRadial =
            std::max(radialDelta, 0.0);
    
        const double negativeAxialDelta =
            point[a] - axialMin;
    
        const double positiveAxialDelta =
            point[a] - axialMax;
    
        // Finite lateral cylindrical surface.
        consider(
            0, 2,
            axialOverlap && circleOverlap,
            lateralDistanceSquared);
    
        // Finite circular end caps.
        consider(
            1, 2,
            negativePlane && diskOverlap,
            negativeAxialDelta * negativeAxialDelta +
                outsideRadial * outsideRadial);
    
        consider(
            2, 2,
            positivePlane && diskOverlap,
            positiveAxialDelta * positiveAxialDelta +
                outsideRadial * outsideRadial);
    
        // Circular rims: lower-dimensional priority.
        consider(
            3, 1,
            negativePlane && circleOverlap,
            negativeAxialDelta * negativeAxialDelta +
                radialDelta * radialDelta);
    
        consider(
            4, 1,
            positivePlane && circleOverlap,
            positiveAxialDelta * positiveAxialDelta +
                radialDelta * radialDelta);
    }

    void validateOpenBox(
        const BoundingBox& box) const
    {
        if(box.width() <= 0.0 ||
           box.height() <= 0.0 ||
           box.depth() <= 0.0)
        {
            throw std::invalid_argument(
                "Open box dimensions must be positive.");
        }

        if(!box.contains(bounds))
        {
            throw std::invalid_argument(
                "Open box must contain the Cylinder bounds.");
        }
    }


    [[nodiscard]]
    BoundingBox readOpenBox() const
    {
        std::cout
            << "Cylinder bounds:\n"
            << "  min = ("
            << bounds.min[0] << ", "
            << bounds.min[1] << ", "
            << bounds.min[2] << ")\n"
            << "  max = ("
            << bounds.max[0] << ", "
            << bounds.max[1] << ", "
            << bounds.max[2] << ")\n";

        Point origin{};
        Point size{};

        std::cout
            << "Enter open box origin x y z: ";

        if(!(std::cin >>
             origin[0] >>
             origin[1] >>
             origin[2]))
        {
            throw std::runtime_error(
                "Failed to read open box origin.");
        }

        std::cout
            << "Enter open box size x y z: ";

        if(!(std::cin >>
             size[0] >>
             size[1] >>
             size[2]))
        {
            throw std::runtime_error(
                "Failed to read open box size.");
        }

        const BoundingBox box{
            origin,
            {
                origin[0] + size[0],
                origin[1] + size[1],
                origin[2] + size[2]
            }
        };

        validateOpenBox(box);

        return box;
    }
};

} // namespace ntic::lbm::geometry
