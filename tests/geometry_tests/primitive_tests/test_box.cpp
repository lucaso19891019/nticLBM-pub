#include "box.hpp"

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
    const Box box{
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

    assert(!box.contains(
        Point{2.0, 4.0, 14.0}));


    const BoundingBox bounds =
        box.boundingBox();


    assert(nearlyEqual(
        bounds.min[0],
        1.0));

    assert(nearlyEqual(
        bounds.min[1],
        2.0));

    assert(nearlyEqual(
        bounds.min[2],
        3.0));

    assert(nearlyEqual(
        bounds.max[0],
        5.0));

    assert(nearlyEqual(
        bounds.max[1],
        8.0));

    assert(nearlyEqual(
        bounds.max[2],
        13.0));


    std::cout
        << "Box tests passed.\n";

    return 0;
}