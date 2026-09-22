#include "circle.hpp"

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
    const Circle circle{
        {1.0, 2.0, 0.0},
        2.0
    };


    assert(circle.contains(
        Point{1.0, 2.0, 0.0}));

    assert(circle.contains(
        Point{3.0, 2.0, 0.0}));

    assert(circle.contains(
        Point{1.0, 4.0, 0.0}));


    assert(!circle.contains(
        Point{3.1, 2.0, 0.0}));

    assert(!circle.contains(
        Point{1.0, 4.1, 0.0}));


    // Circle is a two-dimensional XY primitive.
    // The z coordinate does not participate in containment.
    assert(circle.contains(
        Point{1.0, 2.0, 10.0}));

    assert(circle.contains(
        Point{3.0, 2.0, -10.0}));


    const BoundingBox bounds =
        circle.boundingBox();


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
        0.0));


    std::cout
        << "Circle tests passed.\n";

    return 0;
}