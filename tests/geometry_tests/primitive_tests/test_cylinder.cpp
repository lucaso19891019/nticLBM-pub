#include "cylinder.hpp"

#include <cassert>
#include <cmath>
#include <iostream>


using namespace ntic::lbm::geometry;


namespace
{

bool nearlyEqual(
    const double a,
    const double b,
    const double tolerance = 1.0e-12)
{
    return
        std::abs(a - b) <= tolerance;
}

} // namespace


int main()
{
    {
        const Cylinder cylinder{
            {1.0, 2.0, 3.0},
            2.0,
            6.0,
            Axis::X
        };


        assert(cylinder.contains(
            Point{1.0, 2.0, 3.0}));

        assert(cylinder.contains(
            Point{4.0, 2.0, 3.0}));

        assert(cylinder.contains(
            Point{-2.0, 2.0, 3.0}));

        assert(cylinder.contains(
            Point{1.0, 4.0, 3.0}));


        assert(!cylinder.contains(
            Point{4.1, 2.0, 3.0}));

        assert(!cylinder.contains(
            Point{1.0, 4.1, 3.0}));


        const BoundingBox bounds =
            cylinder.boundingBox();


        assert(nearlyEqual(
            bounds.min[0],
            -2.0));

        assert(nearlyEqual(
            bounds.min[1],
            0.0));

        assert(nearlyEqual(
            bounds.min[2],
            1.0));

        assert(nearlyEqual(
            bounds.max[0],
            4.0));

        assert(nearlyEqual(
            bounds.max[1],
            4.0));

        assert(nearlyEqual(
            bounds.max[2],
            5.0));
    }


    {
        const Cylinder cylinder{
            {1.0, 2.0, 3.0},
            2.0,
            6.0,
            Axis::Y
        };


        assert(cylinder.contains(
            Point{1.0, 2.0, 3.0}));

        assert(cylinder.contains(
            Point{1.0, 5.0, 3.0}));

        assert(cylinder.contains(
            Point{1.0, -1.0, 3.0}));

        assert(cylinder.contains(
            Point{3.0, 2.0, 3.0}));


        assert(!cylinder.contains(
            Point{1.0, 5.1, 3.0}));

        assert(!cylinder.contains(
            Point{3.1, 2.0, 3.0]));


        const BoundingBox bounds =
            cylinder.boundingBox();


        assert(nearlyEqual(
            bounds.min[0],
            -1.0));

        assert(nearlyEqual(
            bounds.min[1],
            -1.0));

        assert(nearlyEqual(
            bounds.min[2],
            1.0));

        assert(nearlyEqual(
            bounds.max[0],
            3.0));

        assert(nearlyEqual(
            bounds.max[1],
            5.0));

        assert(nearlyEqual(
            bounds.max[2],
            5.0));
    }


    {
        const Cylinder cylinder{
            {1.0, 2.0, 3.0},
            2.0,
            6.0,
            Axis::Z
        };


        assert(cylinder.contains(
            Point{1.0, 2.0, 3.0}));

        assert(cylinder.contains(
            Point{1.0, 2.0, 6.0}));

        assert(cylinder.contains(
            Point{1.0, 2.0, 0.0}));

        assert(cylinder.contains(
            Point{3.0, 2.0, 3.0]));


        assert(!cylinder.contains(
            Point{1.0, 2.0, 6.1}));

        assert(!cylinder.contains(
            Point{3.1, 2.0, 3.0]));


        const BoundingBox bounds =
            cylinder.boundingBox();


        assert(nearlyEqual(
            bounds.min[0],
            -1.0));

        assert(nearlyEqual(
            bounds.min[1],
            0.0));

        assert(nearlyEqual(
            bounds.min[2],
            0.0));

        assert(nearlyEqual(
            bounds.max[0],
            3.0));

        assert(nearlyEqual(
            bounds.max[1],
            4.0));

        assert(nearlyEqual(
            bounds.max[2],
            6.0));
    }


    std::cout
        << "Cylinder tests passed.\n";

    return 0;
}