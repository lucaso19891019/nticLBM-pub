#pragma once

#include "bounding_box.hpp"
#include "flow_type.hpp"
#include "point.hpp"

#include <stdexcept>


namespace ntic::lbm::geometry
{

struct Rectangle
{
    Point min{};
    Point max{};


    bool contains(
        const Point& point,
        const FlowType flowType,
        const BoundingBox* openBox = nullptr) const
    {
        const bool strictlyInside =
            point[0] > min[0] &&
            point[0] < max[0] &&

            point[1] > min[1] &&
            point[1] < max[1];


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


        const bool onBoundary =
            (
                point[0] >= min[0] &&
                point[0] <= max[0] &&

                point[1] >= min[1] &&
                point[1] <= max[1]
            ) &&
            (
                point[0] == min[0] ||
                point[0] == max[0] ||

                point[1] == min[1] ||
                point[1] == max[1]
            );


        const bool strictlyOutside =
            !strictlyInside &&
            !onBoundary;


        const bool strictlyInsideOpenBox =
            point[0] > openBox->min[0] &&
            point[0] < openBox->max[0] &&

            point[1] > openBox->min[1] &&
            point[1] < openBox->max[1];


        return
            strictlyInsideOpenBox &&
            strictlyOutside;
    }


    BoundingBox boundingBox() const noexcept
    {
        return BoundingBox{
            {
                min[0],
                min[1],
                0.0
            },
            {
                max[0],
                max[1],
                0.0
            }
        };
    }
};

} // namespace ntic::lbm::geometry