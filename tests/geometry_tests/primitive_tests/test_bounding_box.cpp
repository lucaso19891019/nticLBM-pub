#include "bounding_box.hpp"

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
        const BoundingBox box{
            {1.0, 2.0, 3.0},
            {5.0, 8.0, 13.0}
        };

        assert(nearlyEqual(
            box.width(),
            4.0));

        assert(nearlyEqual(
            box.height(),
            6.0));

        assert(nearlyEqual(
            box.depth(),
            10.0));
    }


    {
        const BoundingBox box{
            {1.0, 2.0, 3.0},
            {5.0, 8.0, 13.0}
        };

        assert(box.contains(
            Point{2.0, 4.0, 6.0}));

        assert(box.contains(
            Point{1.0, 2.0, 3.0}));

        assert(box.contains(
            Point{5.0, 8.0, 13.0}));

        assert(!box.contains(
            Point{0.0, 4.0, 6.0}));

        assert(!box.contains(
            Point{2.0, 9.0, 6.0}));
    }


    {
        const BoundingBox outer{
            {0.0, 0.0, 0.0},
            {10.0, 20.0, 30.0}
        };

        const BoundingBox inner{
            {2.0, 3.0, 4.0},
            {8.0, 15.0, 25.0}
        };

        const BoundingBox touching{
            {0.0, 5.0, 5.0},
            {10.0, 15.0, 20.0}
        };

        const BoundingBox outside{
            {-1.0, 3.0, 4.0},
            {8.0, 15.0, 25.0}
        };

        assert(outer.contains(
            inner));

        assert(outer.contains(
            touching));

        assert(!outer.contains(
            outside));
    }


    {
        BoundingBox box;

        box.expand(
            Point{10.0, 20.0, 30.0});

        assert(nearlyEqual(
            box.min[0],
            10.0));

        assert(nearlyEqual(
            box.min[1],
            20.0));

        assert(nearlyEqual(
            box.min[2],
            30.0));

        assert(nearlyEqual(
            box.max[0],
            10.0));

        assert(nearlyEqual(
            box.max[1],
            20.0));

        assert(nearlyEqual(
            box.max[2],
            30.0));


        box.expand(
            Point{15.0, 18.0, 40.0});

        assert(nearlyEqual(
            box.min[0],
            10.0));

        assert(nearlyEqual(
            box.min[1],
            18.0));

        assert(nearlyEqual(
            box.min[2],
            30.0));

        assert(nearlyEqual(
            box.max[0],
            15.0));

        assert(nearlyEqual(
            box.max[1],
            20.0));

        assert(nearlyEqual(
            box.max[2],
            40.0));
    }


    {
        BoundingBox box;

        box.expand(
            Point{-10.0, -20.0, -30.0});

        assert(nearlyEqual(
            box.min[0],
            -10.0));

        assert(nearlyEqual(
            box.min[1],
            -20.0));

        assert(nearlyEqual(
            box.min[2],
            -30.0));

        assert(nearlyEqual(
            box.max[0],
            -10.0));

        assert(nearlyEqual(
            box.max[1],
            -20.0));

        assert(nearlyEqual(
            box.max[2],
            -30.0));
    }


    std::cout
        << "BoundingBox tests passed.\n";

    return 0;
}