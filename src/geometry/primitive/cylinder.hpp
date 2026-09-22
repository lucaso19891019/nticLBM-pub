#pragma once

#include "point.hpp"
#include "axis.hpp"
#include "bounding_box.hpp"

#include "flow_type.hpp"

#include <cmath>
#include <stdexcept>

namespace ntic::lbm::geometry
{
//=============================================================================
// Cylinder
//=============================================================================
//
// Axis-aligned three-dimensional cylindrical primitive geometry.
//
// The cylinder is defined by:
//
//   - center
//   - radius
//   - length
//   - axis (X, Y, or Z)
//
// The center is located at the geometric center of the cylinder.
//
// Boundary points are considered inside.
//
//=============================================================================

struct Cylinder
{
    Point center;
    double radius;
    double length;
    Axis axis;


    bool contains(
        const Point& point,
        const FlowType flowType,
        const BoundingBox* openBox = nullptr) const
    {
        double axialDistance =
            0.0;

        double radialDistanceSquared =
            0.0;


        switch(axis)
        {
            case Axis::X:
            {
                axialDistance =
                    point[0] - center[0];

                const double dy =
                    point[1] - center[1];

                const double dz =
                    point[2] - center[2];

                radialDistanceSquared =
                    dy * dy +
                    dz * dz;

                break;
            }


            case Axis::Y:
            {
                axialDistance =
                    point[1] - center[1];

                const double dx =
                    point[0] - center[0];

                const double dz =
                    point[2] - center[2];

                radialDistanceSquared =
                    dx * dx +
                    dz * dz;

                break;
            }


            case Axis::Z:
            {
                axialDistance =
                    point[2] - center[2];

                const double dx =
                    point[0] - center[0];

                const double dy =
                    point[1] - center[1];

                radialDistanceSquared =
                    dx * dx +
                    dy * dy;

                break;
            }
        }


        const double halfLength =
            0.5 * length;

        const double radiusSquared =
            radius * radius;

        const double absoluteAxialDistance =
            std::abs(
                axialDistance);


        const bool strictlyInside =
            absoluteAxialDistance <
                halfLength &&
            radialDistanceSquared <
                radiusSquared;


        if(flowType ==
        FlowType::Internal)
        {
            return
                strictlyInside;
        }


        if(openBox == nullptr)
        {
            throw std::invalid_argument(
                "External flow requires an open box.");
        }


        const bool onEndSurface =
            absoluteAxialDistance ==
                halfLength &&
            radialDistanceSquared <=
                radiusSquared;


        const bool onSideSurface =
            absoluteAxialDistance <=
                halfLength &&
            radialDistanceSquared ==
                radiusSquared;


        const bool strictlyOutside =
            !strictlyInside &&
            !onEndSurface &&
            !onSideSurface;


        return
            openBox->strictlyContains(point) &&
            strictlyOutside;
    }


    BoundingBox boundingBox() const noexcept
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
};

} // namespace ntic::lbm::geometry