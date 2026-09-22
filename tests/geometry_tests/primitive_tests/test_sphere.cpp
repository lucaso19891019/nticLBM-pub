#include "sphere.hpp"

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
    const Sphere sphere{
        {1.0, 2.0, 3.0},
        2.0
    };


    assert(sphere.contains(
        Point{1.0, 2.0, 3.0}));

    assert(sphere.contains(
        Point{3.0, 2.0, 3.0}));

    assert(sphere.contains(
        Point{1.0, 4.0, 3.0}));

    assert(sphere.contains(
        Point{1.0, 2.0, 5.0}));


    assert(!sphere.contains(
        Point{3.1, 2.0, 3.0}));

    assert(!sphere.contains(
        Point{1.0, 4.1, 3.0}));

    assert(!sphere.contains(
        Point{1.0, 2.0, 5.1}));


    const BoundingBox bounds =
        sphere.boundingBox();


    assert(nearlyEqual(
        bounds.min[0],
        -1.0));

    assert(nearlyEqual(
        bounds.min[1],
        0.0));

    assert(nearlyEqual(
        bounds.min[2],
        1.0));

    assert(nearlyEqual(
        bounds.max[0],
        3.0));

    assert(nearlyEqual(
        bounds.max[1],
        4.0));

    assert(nearlyEqual(
        bounds.max[2],
        5.0));


    std::cout
        << "Sphere tests passed.\n";

    return 0;
}