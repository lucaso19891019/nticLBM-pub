#pragma once

#include "point.hpp"

#include <algorithm>
#include <limits>


namespace ntic::lbm::geometry
{

struct BoundingBox
{
    Point min{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()
    };

    Point max{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()
    };


    double width() const noexcept
    {
        return
            max[0] - min[0];
    }


    double height() const noexcept
    {
        return
            max[1] - min[1];
    }


    double depth() const noexcept
    {
        return
            max[2] - min[2];
    }


    bool contains(
        const Point& point) const noexcept
    {
        return
            point[0] >= min[0] &&
            point[0] <= max[0] &&

            point[1] >= min[1] &&
            point[1] <= max[1] &&

            point[2] >= min[2] &&
            point[2] <= max[2];
    }


    bool strictlyContains(
        const Point& point) const noexcept
    {
        return
            point[0] > min[0] &&
            point[0] < max[0] &&

            point[1] > min[1] &&
            point[1] < max[1] &&

            point[2] > min[2] &&
            point[2] < max[2];
    }


    bool contains(
        const BoundingBox& box) const noexcept
    {
        return
            box.min[0] >= min[0] &&
            box.max[0] <= max[0] &&

            box.min[1] >= min[1] &&
            box.max[1] <= max[1] &&

            box.min[2] >= min[2] &&
            box.max[2] <= max[2];
    }


    void expand(
        const Point& point) noexcept
    {
        min[0] =
            std::min(
                min[0],
                point[0]);

        min[1] =
            std::min(
                min[1],
                point[1]);

        min[2] =
            std::min(
                min[2],
                point[2]);


        max[0] =
            std::max(
                max[0],
                point[0]);

        max[1] =
            std::max(
                max[1],
                point[1]);

        max[2] =
            std::max(
                max[2],
                point[2]);
    }
};

} // namespace ntic::lbm::geometry