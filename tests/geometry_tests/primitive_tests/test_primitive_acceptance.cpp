#include "box.hpp"
#include "circle.hpp"
#include "cylinder.hpp"
#include "rectangle.hpp"
#include "sphere.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iostream>
#include <omp.h>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace single_boundary_suite
{
using namespace ntic::lbm;
using namespace ntic::lbm::geometry;

namespace
{

void require(bool condition, const char* message)
{
    if(!condition) throw std::runtime_error(message);
}

template<bool ThreeD, typename Geometry, typename Analysis>
void checkSingleBoundary(const Geometry& geometry, Analysis& analysis,
                         bool expectEmpty)
{
    geometry.interiorAreaAnalysis(analysis);
    const auto scalar = analysis.scalar;
    const auto x = analysis.boundaryX;
    const auto y = analysis.boundaryY;
    const auto offsets = analysis.boundaryOffsets;
    std::vector<double> z;
    if constexpr(ThreeD) z = analysis.boundaryZ;
    geometry.boundaryQAnalysis(analysis);

    using Lattice = typename Analysis::Lattice;
    constexpr std::size_t nLinks = Lattice::nStencils - 1;
    const std::size_t count = analysis.boundaryX.size();
    require(analysis.nBoundaries == 1, "Expected one geometric boundary");
    require(analysis.boundaryOffsets == std::vector<std::size_t>{0, count},
            "Invalid single-boundary offsets");
    require(analysis.boundaryY.size() == count, "Invalid Y coordinate count");
    if constexpr(ThreeD)
    {
        require(analysis.boundaryZ.size() == count, "Invalid Z coordinate count");
        require(analysis.boundaryZ == z, "Q analysis changed Z coordinates");
    }
    require(analysis.q.size() == count * nLinks, "Invalid Q block count");
    require((count == 0) == expectEmpty, "Unexpected empty boundary");
    require(analysis.scalar == scalar && analysis.boundaryX == x &&
            analysis.boundaryY == y && analysis.boundaryOffsets == offsets,
            "Q analysis changed classification or grouping");
    require(static_cast<std::size_t>(std::count(
                scalar.begin(), scalar.end(),
                static_cast<double>(CellType::Boundary))) == count,
            "Boundary entries do not cover the classified cells");

    std::vector<bool> seen(scalar.size(), false);
    std::size_t validLinks = 0;
    for(std::size_t b = 0; b < count; ++b)
    {
        const double h = analysis.gridSpacing;
        const auto i = static_cast<std::size_t>(std::llround(
            (analysis.boundaryX[b] - analysis.domain.min[0]) / h - 0.5));
        const auto j = static_cast<std::size_t>(std::llround(
            (analysis.boundaryY[b] - analysis.domain.min[1]) / h - 0.5));
        std::size_t id = i + analysis.nx * j;
        if constexpr(ThreeD)
        {
            const auto k = static_cast<std::size_t>(std::llround(
                (analysis.boundaryZ[b] - analysis.domain.min[2]) / h - 0.5));
            id += analysis.nx * analysis.ny * k;
        }
        require(id < scalar.size(), "Boundary coordinate outside grid");
        require(!seen[id], "Duplicate entry within boundary ID 0");
        seen[id] = true;
        require(scalar[id] == static_cast<double>(CellType::Boundary),
                "Entry is not a Boundary cell");

        for(std::size_t direction = 1; direction < Lattice::nStencils; ++direction)
        {
            const double q = analysis.q[b * nLinks + direction - 1];
            require(q == -1.0 || (std::isfinite(q) && q > 0.0 && q <= 1.0),
                    "Invalid Q range");
            if(q == -1.0) continue;
            ++validLinks;
            const double dx = analysis.boundaryX[b] - geometry.center[0] +
                              q * h * Lattice::ex[direction];
            const double dy = analysis.boundaryY[b] - geometry.center[1] +
                              q * h * Lattice::ey[direction];
            double squaredRadius = dx * dx + dy * dy;
            if constexpr(ThreeD)
            {
                const double dz = analysis.boundaryZ[b] - geometry.center[2] +
                                  q * h * Lattice::ez[direction];
                squaredRadius += dz * dz;
            }
            require(std::abs(squaredRadius - geometry.radius * geometry.radius)
                        <= 1.0e-11 * geometry.radius * geometry.radius,
                    "Q block does not intersect the circular/spherical surface");
        }
    }
    require(expectEmpty || validLinks > 0, "Expected valid boundary links");
}

template<lattice::LatticeType Type>
void checkSpheres()
{
    const BoundingBox openBox{{-2.0, -2.0, -2.0}, {2.0, 2.0, 2.0}};
    for(const auto flow : {FlowType::Internal, FlowType::External})
    {
        const Sphere sphere{{0.0, 0.0, 0.0}, 1.0, flow, &openBox};
        for(const double h : {0.25, 0.4, 4.0})
        {
            GeometryAnalysis3D<Type> analysis(h);
            checkSingleBoundary<true>(sphere, analysis, h == 4.0);
        }
    }
}

} // namespace

int run()
{
    try
    {
        const BoundingBox2D openBox{{-2.0, -2.0}, {2.0, 2.0}};
        for(const auto flow : {FlowType::Internal, FlowType::External})
        {
            const Circle circle{{0.0, 0.0}, 1.0, flow, &openBox};
            for(const double h : {0.25, 0.4, 4.0})
            {
                GeometryAnalysis2D<lattice::LatticeType::D2Q9> analysis(h);
                checkSingleBoundary<false>(circle, analysis, h == 4.0);
            }
        }
        checkSpheres<lattice::LatticeType::D3Q15>();
        checkSpheres<lattice::LatticeType::D3Q19>();
        checkSpheres<lattice::LatticeType::D3Q27>();
        std::cout << "Single-boundary contract checks passed\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
} // namespace single_boundary_suite

namespace rectangle_suite
{
using namespace ntic::lbm;
using namespace ntic::lbm::geometry;
using Analysis = GeometryAnalysis2D<lattice::LatticeType::D2Q9>;
using Lattice = Analysis::Lattice;

namespace
{

void require(bool condition, const char* message)
{
    if(!condition) throw std::runtime_error(message);
}

// Independent finite-segment oracle: intersect each edge separately.
std::pair<double, std::size_t> nearestLink(
    const Rectangle& rectangle, double x, double y, double h, std::size_t link)
{
    double nearest = -1.0;
    std::size_t owner = 4;
    for(std::size_t edge = 0; edge < 4; ++edge)
    {
        const bool vertical = edge == 1 || edge == 3;
        const std::size_t axis = vertical ? 0 : 1;
        const std::size_t other = 1 - axis;
        const double p[2] = {x, y};
        const double d[2] = {h * Lattice::ex[link], h * Lattice::ey[link]};
        if(d[axis] == 0.0) continue;
        const double plane = (edge == 0 || edge == 3)
            ? rectangle.min[axis] : rectangle.max[axis];
        const double q = (plane - p[axis]) / d[axis];
        if(q <= 0.0 || q > 1.0) continue;
        const double crossing = p[other] + q * d[other];
        if(crossing < rectangle.min[other] || crossing > rectangle.max[other])
            continue;
        if(nearest < 0.0 || q < nearest)
        {
            nearest = q;
            owner = edge;
        }
    }
    return {nearest, owner};
}

void checkCase(const Rectangle& rectangle, double h,
               bool& sawDuplicateCell, bool& sawEmpty, bool& sawCornerTie,
               bool& sawBoundaryWithoutQ)
{
    Analysis analysis(h);
    omp_set_num_threads(1);
    rectangle.interiorAreaAnalysis(analysis);
    const auto originalScalar = analysis.scalar;
    rectangle.boundaryQAnalysis(analysis);
    require(analysis.scalar == originalScalar, "Q changed CellType");
    const std::size_t count = analysis.boundaryX.size();
    require(analysis.nBoundaries == 4, "Rectangle must have four edges");
    require(analysis.boundaryOffsets.size() == 5 &&
            analysis.boundaryOffsets.front() == 0 &&
            analysis.boundaryOffsets.back() == count,
            "Invalid edge offsets");
    require(analysis.boundaryY.size() == count && analysis.q.size() == count * 8,
            "Invalid coordinate/Q lengths");

    const std::size_t cellCount = analysis.nx * analysis.ny;
    require(analysis.scalar.size() == cellCount, "Invalid scalar size");
    std::vector<std::array<bool, 4>> entries(cellCount);
    std::vector<std::array<int, 8>> validCount(cellCount);
    std::vector<std::array<double, 8>> actualQ(cellCount);
    std::vector<std::array<std::size_t, 8>> actualOwner(cellCount);
    for(auto& values : actualQ) values.fill(-1.0);

    for(std::size_t edge = 0; edge < 4; ++edge)
    {
        const std::size_t begin = analysis.boundaryOffsets[edge];
        const std::size_t end = analysis.boundaryOffsets[edge + 1];
        require(begin <= end, "Offsets must be monotone");
        sawEmpty = sawEmpty || begin == end;
        std::size_t previous = 0;
        for(std::size_t b = begin; b < end; ++b)
        {
            const double x = analysis.boundaryX[b];
            const double y = analysis.boundaryY[b];
            const auto i = static_cast<std::size_t>(std::llround(
                (x - analysis.domain.min[0]) / h - 0.5));
            const auto j = static_cast<std::size_t>(std::llround(
                (y - analysis.domain.min[1]) / h - 0.5));
            require(i < analysis.nx && j < analysis.ny, "Entry outside grid");
            const std::size_t id = i + analysis.nx * j;
            require(b == begin || previous < id, "Edge order is not stable/unique");
            previous = id;
            entries[id][edge] = true;
            for(std::size_t link = 1; link < Lattice::nStencils; ++link)
            {
                const double q = analysis.q[b * 8 + link - 1];
                require(q == -1.0 || (std::isfinite(q) && q > 0.0 && q <= 1.0),
                        "Invalid Q range");
                if(q == -1.0) continue;
                require(++validCount[id][link - 1] == 1, "Duplicate (Cell, Link) Q");
                actualQ[id][link - 1] = q;
                actualOwner[id][link - 1] = edge;
            }
        }
    }

    for(std::size_t id = 0; id < cellCount; ++id)
    {
        const double x = analysis.domain.min[0] + (id % analysis.nx + 0.5) * h;
        const double y = analysis.domain.min[1] + (id / analysis.nx + 0.5) * h;
        const bool inside = x > rectangle.min[0] && x < rectangle.max[0] &&
                            y > rectangle.min[1] && y < rectangle.max[1];
        const bool onSurface =
            ((x == rectangle.min[0] || x == rectangle.max[0]) &&
             y >= rectangle.min[1] && y <= rectangle.max[1]) ||
            ((y == rectangle.min[1] || y == rectangle.max[1]) &&
             x >= rectangle.min[0] && x <= rectangle.max[0]);
        const bool wet = !onSurface &&
            (rectangle.flowType == FlowType::Internal ? inside : !inside);
        const bool overlapX = x + h >= rectangle.min[0] && x - h <= rectangle.max[0];
        const bool overlapY = y + h >= rectangle.min[1] && y - h <= rectangle.max[1];
        const std::array<bool, 4> expected{
            wet && overlapX && std::abs(y - rectangle.min[1]) <= h,
            wet && overlapY && std::abs(x - rectangle.max[0]) <= h,
            wet && overlapX && std::abs(y - rectangle.max[1]) <= h,
            wet && overlapY && std::abs(x - rectangle.min[0]) <= h
        };
        std::size_t memberships = 0;
        for(const bool member : expected) if(member) ++memberships;
        const CellType expectedType = !wet ? CellType::Dry :
            (memberships != 0 ? CellType::Boundary : CellType::Interior);
        require(analysis.scalar[id] == static_cast<double>(expectedType),
                "Cell classification differs from closed-stencil semantics");
        require(entries[id] == expected, "Incorrect finite-edge membership");
        sawDuplicateCell = sawDuplicateCell || memberships > 1;
        bool hasQ = false;
        for(std::size_t link = 1; link < Lattice::nStencils; ++link)
        {
            if(expectedType != CellType::Boundary) continue;
            const auto hit = nearestLink(rectangle, x, y, h, link);
            const double actual = actualQ[id][link - 1];
            require(std::abs(actual - hit.first) <= 1.0e-12,
                    "Q differs from nearest finite-edge intersection");
            if(hit.first < 0.0) continue;
            hasQ = true;
            require(actualOwner[id][link - 1] == hit.second, "Incorrect Q Edge ID");
            const double ix = x + hit.first * h * Lattice::ex[link];
            const double iy = y + hit.first * h * Lattice::ey[link];
            sawCornerTie = sawCornerTie ||
                ((ix == rectangle.min[0] || ix == rectangle.max[0]) &&
                 (iy == rectangle.min[1] || iy == rectangle.max[1]));
        }
        sawBoundaryWithoutQ = sawBoundaryWithoutQ ||
            (expectedType == CellType::Boundary && !hasQ);
    }

    // Reuse the object and compare complete output across thread counts.
    const Analysis baseline = analysis;
    omp_set_num_threads(4);
    rectangle.analysis(analysis);
    require(analysis.scalar == baseline.scalar &&
            analysis.boundaryX == baseline.boundaryX &&
            analysis.boundaryY == baseline.boundaryY &&
            analysis.boundaryOffsets == baseline.boundaryOffsets &&
            analysis.q == baseline.q, "Thread count changed Rectangle output");
}

} // namespace

int run()
{
    try
    {
        bool sawDuplicateCell = false, sawEmpty = false, sawCornerTie = false;
        bool sawBoundaryWithoutQ = false;
        const BoundingBox2D openBox{{-1.25, -2.0}, {3.75, 3.0}};
        for(const auto flow : {FlowType::Internal, FlowType::External})
        {
            const Rectangle rectangle{{0.25, -0.5}, {2.25, 1.5}, flow, &openBox};
            // Binary-exact fixtures include endpoint/tangent corner hits,
            // surface-centered Dry cells, empty edges, and coarse grids.
            for(const double h : {0.125, 0.25, 0.5, 1.0, 4.0})
                checkCase(rectangle, h, sawDuplicateCell, sawEmpty, sawCornerTie,
                          sawBoundaryWithoutQ);
        }
        require(sawDuplicateCell && sawEmpty && sawCornerTie,
                "Fixtures did not exercise duplicate cells, empty edges, and corners");
        // A stencil can touch the surface without a discrete link hitting it.
        const BoundingBox2D thinOpenBox{{-0.75, -0.75}, {1.25, 1.25}};
        const Rectangle thin{{0.0625, 0.0625}, {0.125, 0.125},
                             FlowType::External, &thinOpenBox};
        checkCase(thin, 0.5, sawDuplicateCell, sawEmpty, sawCornerTie,
                  sawBoundaryWithoutQ);
        require(sawBoundaryWithoutQ, "Missing stencil-only Boundary fixture");
        std::cout << "Rectangle boundary contract checks passed\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
} // namespace rectangle_suite

namespace box_cylinder_suite
{
using namespace ntic::lbm;
using namespace ntic::lbm::geometry;

namespace
{

void require(bool condition, const char* message)
{
    if(!condition) throw std::runtime_error(message);
}

struct Hit
{
    long double q = -1.0L;
    std::size_t face = 0;
    bool tied = false;

    void consider(long double candidate, std::size_t id)
    {
        if(candidate <= 0.0L || candidate > 1.0L) return;
        if(q < 0.0L || candidate < q)
        {
            q = candidate;
            face = id;
            tied = false;
        }
        else if(candidate == q && id != face)
        {
            tied = true;
            if(id < face) face = id;
        }
    }
};

struct Reference
{
    bool wet = false;
    std::vector<bool> faces;
};

Reference classify(const Box& box, const Point& p, double h)
{
    bool inside = true, closedInside = true, surface = false;
    for(std::size_t d = 0; d < 3; ++d)
    {
        inside = inside && p[d] > box.min[d] && p[d] < box.max[d];
        closedInside = closedInside && p[d] >= box.min[d] && p[d] <= box.max[d];
        surface = surface || p[d] == box.min[d] || p[d] == box.max[d];
    }
    Reference result;
    result.wet = !(closedInside && surface) &&
        (box.flowType == FlowType::Internal ? inside : !inside);
    result.faces.resize(6, false);
    for(std::size_t face = 0; face < 6; ++face)
    {
        const std::size_t axis = face / 2;
        const double plane = face % 2 == 0 ? box.min[axis] : box.max[axis];
        bool hit = result.wet && std::abs(p[axis] - plane) <= h;
        for(std::size_t d = 0; d < 3; ++d)
        {
            if(d != axis)
                hit = hit && p[d] + h >= box.min[d] && p[d] - h <= box.max[d];
        }
        result.faces[face] = hit;
    }
    return result;
}

Hit intersect(const Box& box, const Point& p, const Point& delta)
{
    Hit best;
    for(std::size_t face = 0; face < 6; ++face)
    {
        const std::size_t axis = face / 2;
        if(delta[axis] == 0.0) continue;
        const long double plane = face % 2 == 0 ? box.min[axis] : box.max[axis];
        const long double q = (plane - p[axis]) / delta[axis];
        bool finiteFace = true;
        for(std::size_t d = 0; d < 3; ++d)
        {
            if(d == axis) continue;
            const long double crossing = p[d] + q * delta[d];
            finiteFace = finiteFace && crossing >= box.min[d] && crossing <= box.max[d];
        }
        if(finiteFace) best.consider(q, face);
    }
    return best;
}

std::array<std::size_t, 3> dimensions(const Cylinder& cylinder)
{
    if(cylinder.axis == Axis::X) return {0, 1, 2};
    if(cylinder.axis == Axis::Y) return {1, 0, 2};
    return {2, 0, 1};
}

Reference classify(const Cylinder& cylinder, const Point& p, double h)
{
    const auto dims = dimensions(cylinder);
    const double a = p[dims[0]] - cylinder.center[dims[0]];
    const double u = p[dims[1]] - cylinder.center[dims[1]];
    const double v = p[dims[2]] - cylinder.center[dims[2]];
    const double halfLength = cylinder.length / 2.0;
    const double r2 = cylinder.radius * cylinder.radius;
    const double radial2 = u * u + v * v;
    const bool inside = std::abs(a) < halfLength && radial2 < r2;
    const bool onSurface =
        (std::abs(a) == halfLength && radial2 <= r2) ||
        (std::abs(a) <= halfLength && radial2 == r2);
    Reference result;
    result.wet = !onSurface &&
        (cylinder.flowType == FlowType::Internal ? inside : !inside);
    // Radial square's minimum and maximum distance to the axis.
    const double nearU = std::max(std::abs(u) - h, 0.0);
    const double nearV = std::max(std::abs(v) - h, 0.0);
    const double min2 = nearU * nearU + nearV * nearV;
    const double max2 = (std::abs(u) + h) * (std::abs(u) + h) +
                        (std::abs(v) + h) * (std::abs(v) + h);
    result.faces = {
        result.wet && std::abs(a) <= halfLength + h && min2 <= r2 && max2 >= r2,
        result.wet && std::abs(a + halfLength) <= h && min2 <= r2,
        result.wet && std::abs(a - halfLength) <= h && min2 <= r2
    };
    return result;
}

Hit intersect(const Cylinder& cylinder, const Point& point, const Point& delta)
{
    const auto dims = dimensions(cylinder);
    const long double a = point[dims[0]] - cylinder.center[dims[0]];
    const long double u = point[dims[1]] - cylinder.center[dims[1]];
    const long double v = point[dims[2]] - cylinder.center[dims[2]];
    const long double da = delta[dims[0]], du = delta[dims[1]], dv = delta[dims[2]];
    const long double halfLength = cylinder.length / 2.0L;
    const long double r2 = static_cast<long double>(cylinder.radius) * cylinder.radius;
    Hit best;
    // Long-double direct roots are independent of the production stable-root formula.
    const long double A = du * du + dv * dv;
    const long double B = 2.0L * (u * du + v * dv);
    const long double C = u * u + v * v - r2;
    const long double discriminant = B * B - 4.0L * A * C;
    if(A > 0.0L && discriminant >= 0.0L)
    {
        for(const long double sign : {-1.0L, 1.0L})
        {
            const long double root = (-B + sign * std::sqrt(discriminant)) / (2.0L * A);
            if(std::abs(a + root * da) <= halfLength) best.consider(root, 0);
        }
    }
    if(da != 0.0L)
    {
        for(int sign : {-1, 1})
        {
            const long double root = (sign * halfLength - a) / da;
            const long double ru = u + root * du, rv = v + root * dv;
            if(ru * ru + rv * rv <= r2) best.consider(root, sign < 0 ? 1 : 2);
        }
    }
    return best;
}

struct Coverage
{
    bool duplicate = false;
    bool empty = false;
    bool tie = false;
    bool stencilOnly = false;
};

template<lattice::LatticeType Type, typename Geometry>
void checkCase(const Geometry& geometry, double h, std::size_t nFaces, Coverage& coverage)
{
    using Analysis = GeometryAnalysis3D<Type>;
    using Lattice = typename Analysis::Lattice;
    constexpr std::size_t nLinks = Lattice::nStencils - 1;
    Analysis analysis(h);
    omp_set_num_threads(1);
    geometry.interiorAreaAnalysis(analysis);
    const auto scalar = analysis.scalar;
    geometry.boundaryQAnalysis(analysis);
    require(scalar == analysis.scalar, "Q changed classification");
    const std::size_t count = analysis.boundaryX.size();
    require(analysis.nBoundaries == nFaces && analysis.boundaryOffsets.size() == nFaces + 1,
            "Incorrect Face count");
    require(analysis.boundaryOffsets.front() == 0 && analysis.boundaryOffsets.back() == count,
            "Invalid offsets endpoints");
    require(analysis.boundaryY.size() == count && analysis.boundaryZ.size() == count &&
            analysis.q.size() == count * nLinks, "Invalid coordinate/Q size");
    const std::size_t cellCount = analysis.nx * analysis.ny * analysis.nz;
    require(scalar.size() == cellCount, "Invalid scalar size");
    std::vector<std::vector<bool>> membership(cellCount, std::vector<bool>(nFaces, false));
    std::vector<std::array<double, nLinks>> qValues(cellCount);
    std::vector<std::array<std::size_t, nLinks>> owners(cellCount);
    for(auto& values : qValues) values.fill(-1.0);

    for(std::size_t face = 0; face < nFaces; ++face)
    {
        const std::size_t begin = analysis.boundaryOffsets[face];
        const std::size_t end = analysis.boundaryOffsets[face + 1];
        require(begin <= end, "Offsets are not monotone");
        coverage.empty = coverage.empty || begin == end;
        std::size_t previous = 0;
        for(std::size_t b = begin; b < end; ++b)
        {
            const auto i = static_cast<std::size_t>(std::llround(
                (analysis.boundaryX[b] - analysis.domain.min[0]) / h - 0.5));
            const auto j = static_cast<std::size_t>(std::llround(
                (analysis.boundaryY[b] - analysis.domain.min[1]) / h - 0.5));
            const auto k = static_cast<std::size_t>(std::llround(
                (analysis.boundaryZ[b] - analysis.domain.min[2]) / h - 0.5));
            require(i < analysis.nx && j < analysis.ny && k < analysis.nz, "Entry outside grid");
            const std::size_t id = i + analysis.nx * (j + analysis.ny * k);
            require(b == begin || previous < id, "Duplicate or unordered (Face, Cell)");
            previous = id;
            membership[id][face] = true;
            for(std::size_t link = 0; link < nLinks; ++link)
            {
                const double q = analysis.q[b * nLinks + link];
                require(q == -1.0 || (std::isfinite(q) && q > 0.0 && q <= 1.0), "Invalid Q");
                if(q == -1.0) continue;
                require(qValues[id][link] == -1.0, "Duplicate (Cell, Link) Q");
                qValues[id][link] = q;
                owners[id][link] = face;
            }
        }
    }

    for(std::size_t id = 0; id < cellCount; ++id)
    {
        const Point p{
            analysis.domain.min[0] + (id % analysis.nx + 0.5) * h,
            analysis.domain.min[1] + ((id / analysis.nx) % analysis.ny + 0.5) * h,
            analysis.domain.min[2] + (id / (analysis.nx * analysis.ny) + 0.5) * h
        };
        const Reference expected = classify(geometry, p, h);
        std::size_t nMemberships = 0;
        for(bool member : expected.faces) if(member) ++nMemberships;
        const CellType type = !expected.wet ? CellType::Dry :
            (nMemberships ? CellType::Boundary : CellType::Interior);
        require(scalar[id] == static_cast<double>(type), "CellType differs from reference");
        require(membership[id] == expected.faces, "Incorrect finite Face membership");
        coverage.duplicate = coverage.duplicate || nMemberships > 1;
        if(type != CellType::Boundary) continue;
        bool hasQ = false;
        for(std::size_t link = 1; link < Lattice::nStencils; ++link)
        {
            const Point delta{h * Lattice::ex[link], h * Lattice::ey[link], h * Lattice::ez[link]};
            const Hit hit = intersect(geometry, p, delta);
            require(std::abs(static_cast<long double>(qValues[id][link - 1]) - hit.q) <= 1.0e-12L,
                    "Q differs from nearest finite-surface intersection");
            if(hit.q < 0.0L) continue;
            hasQ = true;
            require(owners[id][link - 1] == hit.face, "Incorrect nearest Face/tie ownership");
            coverage.tie = coverage.tie || hit.tied;
        }
        coverage.stencilOnly = coverage.stencilOnly || !hasQ;
    }

    const Analysis baseline = analysis;
    omp_set_num_threads(4);
    geometry.analysis(analysis);
    require(analysis.nBoundaries == baseline.nBoundaries && analysis.scalar == baseline.scalar &&
            analysis.boundaryOffsets == baseline.boundaryOffsets &&
            analysis.boundaryX == baseline.boundaryX && analysis.boundaryY == baseline.boundaryY &&
            analysis.boundaryZ == baseline.boundaryZ && analysis.q == baseline.q,
            "Thread count changed output");
}

template<lattice::LatticeType Type>
void checkModel()
{
    Coverage boxCoverage, cylinderCoverage;
    const BoundingBox openBox{{-2.0, -2.0, -2.0}, {3.0, 3.0, 3.0}};
    for(const auto flow : {FlowType::Internal, FlowType::External})
    {
        const Box box{{0.0, 0.0, 0.0}, {2.0, 2.0, 2.0}, flow, &openBox};
        for(double h : {0.25, 0.5, 1.0, 4.0, 8.0})
            checkCase<Type>(box, h, 6, boxCoverage);
        for(const auto axis : {Axis::X, Axis::Y, Axis::Z})
        {
            const Cylinder cylinder{{0.25, 0.25, 0.25}, 1.25, 2.5, axis, flow, &openBox};
            for(double h : {0.25, 0.5, 1.0, 4.0, 8.0})
                checkCase<Type>(cylinder, h, 3, cylinderCoverage);
        }
    }
    require(boxCoverage.duplicate && boxCoverage.empty && boxCoverage.tie,
            "Missing Box edge/corner or empty Face coverage");
    require(cylinderCoverage.duplicate && cylinderCoverage.empty,
            "Missing Cylinder rim or empty Face coverage");
    // D3Q19/27 have the face diagonals used by these exact rim-tie fixtures.
    if constexpr(Type != lattice::LatticeType::D3Q15)
        require(cylinderCoverage.tie, "Missing exact Cylinder rim-tie coverage");

    const Box thinBox{{0.0625, 0.0625, 0.0625}, {0.125, 0.125, 0.125},
                      FlowType::External, &openBox};
    checkCase<Type>(thinBox, 0.5, 6, boxCoverage);
    const Cylinder thinCylinder{{0.0625, 0.0625, 0.0625}, 0.03125, 0.0625,
                                Axis::Z, FlowType::External, &openBox};
    checkCase<Type>(thinCylinder, 0.5, 3, cylinderCoverage);
    require(boxCoverage.stencilOnly && cylinderCoverage.stencilOnly,
            "Missing stencil-only Boundary coverage");
}

} // namespace

int run()
{
    try
    {
        checkModel<lattice::LatticeType::D3Q15>();
        checkModel<lattice::LatticeType::D3Q19>();
        checkModel<lattice::LatticeType::D3Q27>();
        std::cout << "Box and Cylinder boundary contract checks passed\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
} // namespace box_cylinder_suite

namespace acceptance_audit
{
// Historical validation data from pre-refactor commit
// 26819215034af3f1aeb16273c4b72b0186427c58, generated on Donau/Metax.
// Fingerprints cover CellType and physical (Cell, Link) Q, independent of grouping.
// The historical Rectangle scalar buffer was preallocated to avoid its old
// missing-resize bug; classification and Q formulas were unchanged.
constexpr std::array<std::array<std::uint64_t, 2>, 306> primitiveBaseline{{
    {UINT64_C(17529772238581266725), UINT64_C(4319174414620639781)},
    {UINT64_C(6005294088549850407), UINT64_C(8048700912394160109)},
    {UINT64_C(16105219819455793831), UINT64_C(4300002631480217651)},
    {UINT64_C(4593983902668796263), UINT64_C(5275651730362915497)},
    {UINT64_C(12261239528159138021), UINT64_C(6080320474771328421)},
    {UINT64_C(10028924121191548612), UINT64_C(7639070646424335441)},
    {UINT64_C(16412019938526469761), UINT64_C(12529288243611296293)},
    {UINT64_C(190511049365602690), UINT64_C(5355070838892192445)},
    {UINT64_C(9808139124141015426), UINT64_C(10683476381634026869)},
    {UINT64_C(14524487114347692418), UINT64_C(790106436341850501)},
    {UINT64_C(4906859039572279682), UINT64_C(9389505303891832101)},
    {UINT64_C(9534658432582202753), UINT64_C(5855745380169270069)},
    {UINT64_C(13692813073196568256), UINT64_C(4478744923719171035)},
    {UINT64_C(2098813954852614976), UINT64_C(4559092065829810822)},
    {UINT64_C(5727161853241798720), UINT64_C(11570389910607401874)},
    {UINT64_C(1516871628827865152), UINT64_C(12432055392918947522)},
    {UINT64_C(16056160626730960708), UINT64_C(178259736477288837)},
    {UINT64_C(16056160626730960708), UINT64_C(178259736477288837)},
    {UINT64_C(18288476033698550117), UINT64_C(16895400119619834044)},
    {UINT64_C(18288476033698550117), UINT64_C(14975771004508513196)},
    {UINT64_C(18288476033698550117), UINT64_C(5677929284519584360)},
    {UINT64_C(16412019938526469761), UINT64_C(16804616312541776677)},
    {UINT64_C(190511049365602690), UINT64_C(3392250210262164133)},
    {UINT64_C(9808139124141015426), UINT64_C(15481334951995106293)},
    {UINT64_C(14524487114347692418), UINT64_C(3385195753324462965)},
    {UINT64_C(4906859039572279682), UINT64_C(8810206873614456869)},
    {UINT64_C(9534658432582202753), UINT64_C(13666135582079685517)},
    {UINT64_C(13692813073196568256), UINT64_C(14268824917809511564)},
    {UINT64_C(2098813954852614976), UINT64_C(455374634889707937)},
    {UINT64_C(5727161853241798720), UINT64_C(6995411927532299253)},
    {UINT64_C(1516871628827865152), UINT64_C(13886442423011866509)},
    {UINT64_C(16056160626730960708), UINT64_C(3389907508410574021)},
    {UINT64_C(16056160626730960708), UINT64_C(3389907508410574021)},
    {UINT64_C(18288476033698550117), UINT64_C(7618192059014386001)},
    {UINT64_C(18288476033698550117), UINT64_C(2592998699419060689)},
    {UINT64_C(18288476033698550117), UINT64_C(8953675158504359057)},
    {UINT64_C(16412019938526469761), UINT64_C(13133158651486399013)},
    {UINT64_C(190511049365602690), UINT64_C(8833730871417054165)},
    {UINT64_C(9808139124141015426), UINT64_C(8815223170328241781)},
    {UINT64_C(14524487114347692418), UINT64_C(2736294282833495957)},
    {UINT64_C(4906859039572279682), UINT64_C(7209592730590210021)},
    {UINT64_C(9534658432582202753), UINT64_C(3944927049747119983)},
    {UINT64_C(13692813073196568256), UINT64_C(3022064510511092580)},
    {UINT64_C(2098813954852614976), UINT64_C(9175407157885834066)},
    {UINT64_C(5727161853241798720), UINT64_C(1416450631847663754)},
    {UINT64_C(1516871628827865152), UINT64_C(15438394364018957958)},
    {UINT64_C(16056160626730960708), UINT64_C(13205604908187629381)},
    {UINT64_C(16056160626730960708), UINT64_C(13205604908187629381)},
    {UINT64_C(18288476033698550117), UINT64_C(5039524975278124720)},
    {UINT64_C(18288476033698550117), UINT64_C(10920431327725718144)},
    {UINT64_C(18288476033698550117), UINT64_C(14294074852666228920)},
    {UINT64_C(6118361367620604389), UINT64_C(6793616165656203301)},
    {UINT64_C(6200362462167956967), UINT64_C(12605157513019236901)},
    {UINT64_C(1103322887634160101), UINT64_C(7426809150745912625)},
    {UINT64_C(2514334355018534629), UINT64_C(13747833417089453694)},
    {UINT64_C(18353432366564669156), UINT64_C(14873418896780126869)},
    {UINT64_C(18353432366564669156), UINT64_C(728185502825205492)},
    {UINT64_C(7973526579148381967), UINT64_C(10724813643074296869)},
    {UINT64_C(2574289612368350285), UINT64_C(8605976230149345421)},
    {UINT64_C(12913388442691800141), UINT64_C(6915666455450334453)},
    {UINT64_C(14635107636427887181), UINT64_C(8868437321807184597)},
    {UINT64_C(4296008806104437325), UINT64_C(17703086099892501493)},
    {UINT64_C(5480526077471321355), UINT64_C(10898457597808254921)},
    {UINT64_C(1025508559937555433), UINT64_C(15721145385470726077)},
    {UINT64_C(7840475474084902027), UINT64_C(7381596977730045835)},
    {UINT64_C(6090586613237727883), UINT64_C(7963292140317423311)},
    {UINT64_C(11560964049533452939), UINT64_C(15857334094968305011)},
    {UINT64_C(7772328307559942662), UINT64_C(16086766315172914597)},
    {UINT64_C(7772328307559942662), UINT64_C(18304076119469168598)},
    {UINT64_C(7772328307559942662), UINT64_C(5339724641159779088)},
    {UINT64_C(7772328307559942662), UINT64_C(6368301153868989460)},
    {UINT64_C(7772328307559942662), UINT64_C(6883198673966640652)},
    {UINT64_C(7973526579148381967), UINT64_C(2877094519536596005)},
    {UINT64_C(2574289612368350285), UINT64_C(15275405273028215581)},
    {UINT64_C(12913388442691800141), UINT64_C(458232027823803781)},
    {UINT64_C(14635107636427887181), UINT64_C(13266614496316230157)},
    {UINT64_C(4296008806104437325), UINT64_C(15578313636299799117)},
    {UINT64_C(5480526077471321355), UINT64_C(16592183088108073277)},
    {UINT64_C(1025508559937555433), UINT64_C(12810374562225286336)},
    {UINT64_C(7840475474084902027), UINT64_C(15012954339648264717)},
    {UINT64_C(6090586613237727883), UINT64_C(6942857930343463169)},
    {UINT64_C(11560964049533452939), UINT64_C(2793614193456337453)},
    {UINT64_C(7772328307559942662), UINT64_C(7381888576150164005)},
    {UINT64_C(7772328307559942662), UINT64_C(1089840772317338633)},
    {UINT64_C(7772328307559942662), UINT64_C(13630794864739916240)},
    {UINT64_C(7772328307559942662), UINT64_C(8452227174923795460)},
    {UINT64_C(7772328307559942662), UINT64_C(16422466617307186896)},
    {UINT64_C(7973526579148381967), UINT64_C(14648205755822065445)},
    {UINT64_C(2574289612368350285), UINT64_C(17902755944630931597)},
    {UINT64_C(12913388442691800141), UINT64_C(6384911256136576101)},
    {UINT64_C(14635107636427887181), UINT64_C(13997255305543674941)},
    {UINT64_C(4296008806104437325), UINT64_C(13971631135704189637)},
    {UINT64_C(5480526077471321355), UINT64_C(4047563972814338569)},
    {UINT64_C(1025508559937555433), UINT64_C(5194050551578385491)},
    {UINT64_C(7840475474084902027), UINT64_C(16292894234588561095)},
    {UINT64_C(6090586613237727883), UINT64_C(10501438993531793323)},
    {UINT64_C(11560964049533452939), UINT64_C(16078310437440349051)},
    {UINT64_C(7772328307559942662), UINT64_C(8142422980555827093)},
    {UINT64_C(7772328307559942662), UINT64_C(14195512183508244246)},
    {UINT64_C(7772328307559942662), UINT64_C(14123632840540065233)},
    {UINT64_C(7772328307559942662), UINT64_C(2610000272319599469)},
    {UINT64_C(7772328307559942662), UINT64_C(11407892658350104145)},
    {UINT64_C(17529772238581266725), UINT64_C(4319174414620639781)},
    {UINT64_C(6005294088549850407), UINT64_C(8048700912394160109)},
    {UINT64_C(16105219819455793831), UINT64_C(4300002631480217651)},
    {UINT64_C(4593983902668796263), UINT64_C(5275651730362915497)},
    {UINT64_C(12261239528159138021), UINT64_C(6080320474771328421)},
    {UINT64_C(10028924121191548612), UINT64_C(7639070646424335441)},
    {UINT64_C(16412019938526469761), UINT64_C(12529288243611296293)},
    {UINT64_C(190511049365602690), UINT64_C(5355070838892192445)},
    {UINT64_C(9808139124141015426), UINT64_C(10683476381634026869)},
    {UINT64_C(14524487114347692418), UINT64_C(790106436341850501)},
    {UINT64_C(4906859039572279682), UINT64_C(9389505303891832101)},
    {UINT64_C(9534658432582202753), UINT64_C(5855745380169270069)},
    {UINT64_C(13692813073196568256), UINT64_C(4478744923719171035)},
    {UINT64_C(2098813954852614976), UINT64_C(4559092065829810822)},
    {UINT64_C(5727161853241798720), UINT64_C(11570389910607401874)},
    {UINT64_C(1516871628827865152), UINT64_C(12432055392918947522)},
    {UINT64_C(16056160626730960708), UINT64_C(178259736477288837)},
    {UINT64_C(16056160626730960708), UINT64_C(178259736477288837)},
    {UINT64_C(18288476033698550117), UINT64_C(16895400119619834044)},
    {UINT64_C(18288476033698550117), UINT64_C(14975771004508513196)},
    {UINT64_C(18288476033698550117), UINT64_C(5677929284519584360)},
    {UINT64_C(16412019938526469761), UINT64_C(16804616312541776677)},
    {UINT64_C(190511049365602690), UINT64_C(3392250210262164133)},
    {UINT64_C(9808139124141015426), UINT64_C(15481334951995106293)},
    {UINT64_C(14524487114347692418), UINT64_C(3385195753324462965)},
    {UINT64_C(4906859039572279682), UINT64_C(8810206873614456869)},
    {UINT64_C(9534658432582202753), UINT64_C(13666135582079685517)},
    {UINT64_C(13692813073196568256), UINT64_C(14268824917809511564)},
    {UINT64_C(2098813954852614976), UINT64_C(455374634889707937)},
    {UINT64_C(5727161853241798720), UINT64_C(6995411927532299253)},
    {UINT64_C(1516871628827865152), UINT64_C(13886442423011866509)},
    {UINT64_C(16056160626730960708), UINT64_C(3389907508410574021)},
    {UINT64_C(16056160626730960708), UINT64_C(3389907508410574021)},
    {UINT64_C(18288476033698550117), UINT64_C(7618192059014386001)},
    {UINT64_C(18288476033698550117), UINT64_C(2592998699419060689)},
    {UINT64_C(18288476033698550117), UINT64_C(8953675158504359057)},
    {UINT64_C(16412019938526469761), UINT64_C(13133158651486399013)},
    {UINT64_C(190511049365602690), UINT64_C(8833730871417054165)},
    {UINT64_C(9808139124141015426), UINT64_C(8815223170328241781)},
    {UINT64_C(14524487114347692418), UINT64_C(2736294282833495957)},
    {UINT64_C(4906859039572279682), UINT64_C(7209592730590210021)},
    {UINT64_C(9534658432582202753), UINT64_C(3944927049747119983)},
    {UINT64_C(13692813073196568256), UINT64_C(3022064510511092580)},
    {UINT64_C(2098813954852614976), UINT64_C(9175407157885834066)},
    {UINT64_C(5727161853241798720), UINT64_C(1416450631847663754)},
    {UINT64_C(1516871628827865152), UINT64_C(15438394364018957958)},
    {UINT64_C(16056160626730960708), UINT64_C(13205604908187629381)},
    {UINT64_C(16056160626730960708), UINT64_C(13205604908187629381)},
    {UINT64_C(18288476033698550117), UINT64_C(5039524975278124720)},
    {UINT64_C(18288476033698550117), UINT64_C(10920431327725718144)},
    {UINT64_C(18288476033698550117), UINT64_C(14294074852666228920)},
    {UINT64_C(6118361367620604389), UINT64_C(5686763850236306269)},
    {UINT64_C(6200362462167956967), UINT64_C(9308562019327599457)},
    {UINT64_C(1103322887634160101), UINT64_C(2495938664264873149)},
    {UINT64_C(2514334355018534629), UINT64_C(7022098130113497353)},
    {UINT64_C(18353432366564669156), UINT64_C(12068665635328776853)},
    {UINT64_C(18353432366564669156), UINT64_C(11761776663144745878)},
    {UINT64_C(7973526579148381967), UINT64_C(18637546818426613)},
    {UINT64_C(2574289612368350285), UINT64_C(14512821094547214861)},
    {UINT64_C(12913388442691800141), UINT64_C(8336256873765065021)},
    {UINT64_C(14635107636427887181), UINT64_C(8008496248792891785)},
    {UINT64_C(4296008806104437325), UINT64_C(18203409950178298385)},
    {UINT64_C(5480526077471321355), UINT64_C(7395928214074335681)},
    {UINT64_C(1025508559937555433), UINT64_C(13942266533849210910)},
    {UINT64_C(7840475474084902027), UINT64_C(18438630782429201586)},
    {UINT64_C(6090586613237727883), UINT64_C(1555585555534952770)},
    {UINT64_C(11560964049533452939), UINT64_C(403432038223259738)},
    {UINT64_C(7772328307559942662), UINT64_C(6222858710957859109)},
    {UINT64_C(7772328307559942662), UINT64_C(11240928467316906623)},
    {UINT64_C(7772328307559942662), UINT64_C(7374878176824610192)},
    {UINT64_C(7772328307559942662), UINT64_C(13098316280003358934)},
    {UINT64_C(7772328307559942662), UINT64_C(794725453854245630)},
    {UINT64_C(7973526579148381967), UINT64_C(3623360142711291301)},
    {UINT64_C(2574289612368350285), UINT64_C(10089316379041623285)},
    {UINT64_C(12913388442691800141), UINT64_C(17719302089510356365)},
    {UINT64_C(14635107636427887181), UINT64_C(4640688663459797617)},
    {UINT64_C(4296008806104437325), UINT64_C(223948707857011865)},
    {UINT64_C(5480526077471321355), UINT64_C(11207822778065945269)},
    {UINT64_C(1025508559937555433), UINT64_C(16456806555875394888)},
    {UINT64_C(7840475474084902027), UINT64_C(23081477164562885)},
    {UINT64_C(6090586613237727883), UINT64_C(18043302670704068880)},
    {UINT64_C(11560964049533452939), UINT64_C(312975905786305420)},
    {UINT64_C(7772328307559942662), UINT64_C(9395448374323094181)},
    {UINT64_C(7772328307559942662), UINT64_C(2656349953407302501)},
    {UINT64_C(7772328307559942662), UINT64_C(2547672251331875920)},
    {UINT64_C(7772328307559942662), UINT64_C(1486547923951691632)},
    {UINT64_C(7772328307559942662), UINT64_C(1352849872122568040)},
    {UINT64_C(7973526579148381967), UINT64_C(6385577496219341973)},
    {UINT64_C(2574289612368350285), UINT64_C(17603013259568581269)},
    {UINT64_C(12913388442691800141), UINT64_C(7317795772112164285)},
    {UINT64_C(14635107636427887181), UINT64_C(13582218787879714649)},
    {UINT64_C(4296008806104437325), UINT64_C(8744812649567398393)},
    {UINT64_C(5480526077471321355), UINT64_C(10460349299943950705)},
    {UINT64_C(1025508559937555433), UINT64_C(5258193340624821431)},
    {UINT64_C(7840475474084902027), UINT64_C(8237717772653832370)},
    {UINT64_C(6090586613237727883), UINT64_C(8320256125309131195)},
    {UINT64_C(11560964049533452939), UINT64_C(9274969340286476823)},
    {UINT64_C(7772328307559942662), UINT64_C(11556589936170044949)},
    {UINT64_C(7772328307559942662), UINT64_C(18431797050028854427)},
    {UINT64_C(7772328307559942662), UINT64_C(8933904411293669457)},
    {UINT64_C(7772328307559942662), UINT64_C(5442482372789353649)},
    {UINT64_C(7772328307559942662), UINT64_C(8247452009781323545)},
    {UINT64_C(17529772238581266725), UINT64_C(4319174414620639781)},
    {UINT64_C(6005294088549850407), UINT64_C(8048700912394160109)},
    {UINT64_C(16105219819455793831), UINT64_C(4300002631480217651)},
    {UINT64_C(4593983902668796263), UINT64_C(5275651730362915497)},
    {UINT64_C(12261239528159138021), UINT64_C(6080320474771328421)},
    {UINT64_C(10028924121191548612), UINT64_C(7639070646424335441)},
    {UINT64_C(16412019938526469761), UINT64_C(12529288243611296293)},
    {UINT64_C(190511049365602690), UINT64_C(5355070838892192445)},
    {UINT64_C(9808139124141015426), UINT64_C(10683476381634026869)},
    {UINT64_C(14524487114347692418), UINT64_C(790106436341850501)},
    {UINT64_C(4906859039572279682), UINT64_C(9389505303891832101)},
    {UINT64_C(9534658432582202753), UINT64_C(5855745380169270069)},
    {UINT64_C(13692813073196568256), UINT64_C(4478744923719171035)},
    {UINT64_C(2098813954852614976), UINT64_C(4559092065829810822)},
    {UINT64_C(5727161853241798720), UINT64_C(11570389910607401874)},
    {UINT64_C(1516871628827865152), UINT64_C(12432055392918947522)},
    {UINT64_C(16056160626730960708), UINT64_C(178259736477288837)},
    {UINT64_C(16056160626730960708), UINT64_C(178259736477288837)},
    {UINT64_C(18288476033698550117), UINT64_C(16895400119619834044)},
    {UINT64_C(18288476033698550117), UINT64_C(14975771004508513196)},
    {UINT64_C(18288476033698550117), UINT64_C(5677929284519584360)},
    {UINT64_C(16412019938526469761), UINT64_C(16804616312541776677)},
    {UINT64_C(190511049365602690), UINT64_C(3392250210262164133)},
    {UINT64_C(9808139124141015426), UINT64_C(15481334951995106293)},
    {UINT64_C(14524487114347692418), UINT64_C(3385195753324462965)},
    {UINT64_C(4906859039572279682), UINT64_C(8810206873614456869)},
    {UINT64_C(9534658432582202753), UINT64_C(13666135582079685517)},
    {UINT64_C(13692813073196568256), UINT64_C(14268824917809511564)},
    {UINT64_C(2098813954852614976), UINT64_C(455374634889707937)},
    {UINT64_C(5727161853241798720), UINT64_C(6995411927532299253)},
    {UINT64_C(1516871628827865152), UINT64_C(13886442423011866509)},
    {UINT64_C(16056160626730960708), UINT64_C(3389907508410574021)},
    {UINT64_C(16056160626730960708), UINT64_C(3389907508410574021)},
    {UINT64_C(18288476033698550117), UINT64_C(7618192059014386001)},
    {UINT64_C(18288476033698550117), UINT64_C(2592998699419060689)},
    {UINT64_C(18288476033698550117), UINT64_C(8953675158504359057)},
    {UINT64_C(16412019938526469761), UINT64_C(13133158651486399013)},
    {UINT64_C(190511049365602690), UINT64_C(8833730871417054165)},
    {UINT64_C(9808139124141015426), UINT64_C(8815223170328241781)},
    {UINT64_C(14524487114347692418), UINT64_C(2736294282833495957)},
    {UINT64_C(4906859039572279682), UINT64_C(7209592730590210021)},
    {UINT64_C(9534658432582202753), UINT64_C(3944927049747119983)},
    {UINT64_C(13692813073196568256), UINT64_C(3022064510511092580)},
    {UINT64_C(2098813954852614976), UINT64_C(9175407157885834066)},
    {UINT64_C(5727161853241798720), UINT64_C(1416450631847663754)},
    {UINT64_C(1516871628827865152), UINT64_C(15438394364018957958)},
    {UINT64_C(16056160626730960708), UINT64_C(13205604908187629381)},
    {UINT64_C(16056160626730960708), UINT64_C(13205604908187629381)},
    {UINT64_C(18288476033698550117), UINT64_C(5039524975278124720)},
    {UINT64_C(18288476033698550117), UINT64_C(10920431327725718144)},
    {UINT64_C(18288476033698550117), UINT64_C(14294074852666228920)},
    {UINT64_C(6118361367620604389), UINT64_C(15316090990454498013)},
    {UINT64_C(6200362462167956967), UINT64_C(323346772523691577)},
    {UINT64_C(1103322887634160101), UINT64_C(18372777657055449081)},
    {UINT64_C(2514334355018534629), UINT64_C(11696109613017217016)},
    {UINT64_C(16619548789707927269), UINT64_C(18289371085434518465)},
    {UINT64_C(18353432366564669156), UINT64_C(7644459062064234130)},
    {UINT64_C(7973526579148381967), UINT64_C(12216562495307517941)},
    {UINT64_C(2574289612368350285), UINT64_C(2649050479276245165)},
    {UINT64_C(12913388442691800141), UINT64_C(18252953203326120573)},
    {UINT64_C(14635107636427887181), UINT64_C(2103061363339106233)},
    {UINT64_C(4296008806104437325), UINT64_C(13979528039607224993)},
    {UINT64_C(5480526077471321355), UINT64_C(5242125516557665001)},
    {UINT64_C(1025508559937555433), UINT64_C(12877691611464030631)},
    {UINT64_C(7840475474084902027), UINT64_C(10986135676882349751)},
    {UINT64_C(6090586613237727883), UINT64_C(2822173326669380445)},
    {UINT64_C(11560964049533452939), UINT64_C(10704312881701829737)},
    {UINT64_C(2318667767331600903), UINT64_C(7737282333605346085)},
    {UINT64_C(7772328307559942662), UINT64_C(10478502942685121073)},
    {UINT64_C(7772328307559942662), UINT64_C(7104734646058285735)},
    {UINT64_C(7772328307559942662), UINT64_C(9779291802088990706)},
    {UINT64_C(7772328307559942662), UINT64_C(17237424088601802514)},
    {UINT64_C(7973526579148381967), UINT64_C(3426888703566113317)},
    {UINT64_C(2574289612368350285), UINT64_C(4467486569087384645)},
    {UINT64_C(12913388442691800141), UINT64_C(12875547261353827773)},
    {UINT64_C(14635107636427887181), UINT64_C(17687672987865575969)},
    {UINT64_C(4296008806104437325), UINT64_C(14177984405988245041)},
    {UINT64_C(5480526077471321355), UINT64_C(6085151586173160781)},
    {UINT64_C(1025508559937555433), UINT64_C(6243233335779126609)},
    {UINT64_C(7840475474084902027), UINT64_C(10901810581920758613)},
    {UINT64_C(6090586613237727883), UINT64_C(9065935095112014545)},
    {UINT64_C(11560964049533452939), UINT64_C(13136395680976341433)},
    {UINT64_C(2318667767331600903), UINT64_C(1552800623604058384)},
    {UINT64_C(7772328307559942662), UINT64_C(16307246027150490149)},
    {UINT64_C(7772328307559942662), UINT64_C(5868386494180115063)},
    {UINT64_C(7772328307559942662), UINT64_C(5474075761079556139)},
    {UINT64_C(7772328307559942662), UINT64_C(8847349230907604991)},
    {UINT64_C(7973526579148381967), UINT64_C(15497643636461889557)},
    {UINT64_C(2574289612368350285), UINT64_C(7005413394709144205)},
    {UINT64_C(12913388442691800141), UINT64_C(17344703416369517917)},
    {UINT64_C(14635107636427887181), UINT64_C(3646675969123731529)},
    {UINT64_C(4296008806104437325), UINT64_C(11226898609948075025)},
    {UINT64_C(5480526077471321355), UINT64_C(17777757937478391401)},
    {UINT64_C(1025508559937555433), UINT64_C(17592772885626733857)},
    {UINT64_C(7840475474084902027), UINT64_C(15071043382335160035)},
    {UINT64_C(6090586613237727883), UINT64_C(16086228777408884381)},
    {UINT64_C(11560964049533452939), UINT64_C(17443398384089126361)},
    {UINT64_C(2318667767331600903), UINT64_C(711754513150529121)},
    {UINT64_C(7772328307559942662), UINT64_C(15673842636445652578)},
    {UINT64_C(7772328307559942662), UINT64_C(7014822769712745898)},
    {UINT64_C(7772328307559942662), UINT64_C(16998303407575202741)},
    {UINT64_C(7772328307559942662), UINT64_C(398774488559822105)},
}};


using namespace ntic::lbm;
using namespace ntic::lbm::geometry;
namespace shapes = ntic::lbm::geometry;

void require(bool condition, const char* message)
{
    if(!condition) throw std::runtime_error(message);
}

void hashWord(std::uint64_t& hash, std::uint64_t word)
{
    for(unsigned byte = 0; byte < 8; ++byte)
    {
        hash ^= (word >> (8 * byte)) & 255;
        hash *= UINT64_C(1099511628211);
    }
}

struct Snapshot
{
    std::uint64_t scalar = UINT64_C(14695981039346656037);
    std::uint64_t links = UINT64_C(14695981039346656037);
};

template<bool ThreeD, typename Analysis>
Snapshot snapshot(const Analysis& analysis)
{
    using Lattice = typename Analysis::Lattice;
    constexpr std::size_t nLinks = Lattice::nStencils - 1;
    require(analysis.boundaryX.size() == analysis.boundaryY.size(), "Mismatched XY sizes");
    if constexpr(ThreeD)
        require(analysis.boundaryX.size() == analysis.boundaryZ.size(), "Mismatched Z size");
    require(analysis.q.size() == analysis.boundaryX.size() * nLinks, "Mismatched Q size");
    std::vector<double> canonical(analysis.scalar.size() * nLinks, -1.0);
    for(std::size_t b = 0; b < analysis.boundaryX.size(); ++b)
    {
        const auto i = static_cast<std::size_t>(std::llround(
            (analysis.boundaryX[b] - analysis.domain.min[0]) / analysis.gridSpacing - 0.5));
        const auto j = static_cast<std::size_t>(std::llround(
            (analysis.boundaryY[b] - analysis.domain.min[1]) / analysis.gridSpacing - 0.5));
        require(i < analysis.nx && j < analysis.ny, "Invalid XY cell index");
        std::size_t id = i + analysis.nx * j;
        if constexpr(ThreeD)
        {
            const auto k = static_cast<std::size_t>(std::llround(
                (analysis.boundaryZ[b] - analysis.domain.min[2]) / analysis.gridSpacing - 0.5));
            require(k < analysis.nz, "Invalid Z cell index");
            id += analysis.nx * analysis.ny * k;
        }
        require(analysis.scalar[id] == static_cast<double>(CellType::Boundary),
                "Boundary entry refers to a non-Boundary cell");
        for(std::size_t link = 0; link < nLinks; ++link)
        {
            const double q = analysis.q[b * nLinks + link];
            require(q == -1.0 || (std::isfinite(q) && q > 0.0 && q <= 1.0), "Invalid Q range");
            if(q < 0.0) continue;
            double& slot = canonical[id * nLinks + link];
            require(slot == -1.0, "Duplicate valid physical (Cell, Link)");
            slot = q;
        }
    }
    Snapshot result;
    hashWord(result.scalar, analysis.nx);
    hashWord(result.scalar, analysis.ny);
    if constexpr(ThreeD) hashWord(result.scalar, analysis.nz);
    for(double value : analysis.scalar)
        hashWord(result.scalar, static_cast<std::uint64_t>(value));
    for(double value : canonical)
    {
        std::uint64_t bits;
        static_assert(sizeof(bits) == sizeof(value), "Requires 64-bit double");
        std::memcpy(&bits, &value, sizeof(bits));
        hashWord(result.links, bits);
    }
    return result;
}

template<bool ThreeD, typename Geometry, typename Analysis>
void repeat(const Geometry& geometry, Analysis& analysis)
{
    const Analysis baseline = analysis;
    for(int threads : {1, 4, 4})
    {
        omp_set_num_threads(threads);
        geometry.analysis(analysis);
        require(analysis.nBoundaries == baseline.nBoundaries &&
                analysis.scalar == baseline.scalar && analysis.boundaryX == baseline.boundaryX &&
                analysis.boundaryY == baseline.boundaryY &&
                analysis.boundaryOffsets == baseline.boundaryOffsets && analysis.q == baseline.q,
                "Repeated/threaded output differs");
        if constexpr(ThreeD)
            require(analysis.boundaryZ == baseline.boundaryZ, "Repeated Z output differs");
    }
}

std::size_t caseNumber = 0;

template<bool ThreeD, lattice::LatticeType Type, typename Geometry>
void baselineCase(const Geometry& geometry, double h, const char* name)
{
    using Analysis = std::conditional_t<ThreeD, GeometryAnalysis3D<Type>, GeometryAnalysis2D<Type>>;
    Analysis analysis(h);
    omp_set_num_threads(1);
    geometry.analysis(analysis);
    const Snapshot result = snapshot<ThreeD>(analysis);
    if(caseNumber >= primitiveBaseline.size() || result.scalar != primitiveBaseline[caseNumber][0] ||
       result.links != primitiveBaseline[caseNumber][1])
    {
        throw std::runtime_error(std::string("Historical regression failed: case ") +
            std::to_string(caseNumber) + " " + name + " h=" + std::to_string(h) +
            " scalar=" + std::to_string(result.scalar) + " links=" + std::to_string(result.links));
    }
    repeat<ThreeD>(geometry, analysis);
    ++caseNumber;
}

template<lattice::LatticeType Type>
void baseline3D(double shift, FlowType flow, const BoundingBox& openBox)
{
    const Point center{0.25 + shift, 0.25, 0.25};
    const shapes::Box box{{shift, 0.0, 0.0}, {2.0 + shift, 2.0, 2.0}, flow, &openBox};
    const shapes::Sphere sphere{center, 1.25, flow, &openBox};
    for(double h : {0.5, 0.37, 4.0})
    {
        baselineCase<true, Type>(box, h, "Box");
        baselineCase<true, Type>(sphere, h, "Sphere");
        for(auto axis : {Axis::X, Axis::Y, Axis::Z})
        {
            const shapes::Cylinder cylinder{center, 1.25, 2.5, axis, flow, &openBox};
            baselineCase<true, Type>(cylinder, h, "Cylinder");
        }
    }
}

void historicalRegression()
{
    const BoundingBox2D open2D{{-2.0, -2.0}, {3.0, 3.0}};
    const BoundingBox open3D{{-2.0, -2.0, -2.0}, {3.0, 3.0, 3.0}};
    const double perturbation = std::ldexp(1.0, -36);
    for(double shift : {0.0, -perturbation, perturbation})
    {
        for(auto flow : {FlowType::Internal, FlowType::External})
        {
            const shapes::Rectangle rectangle{{shift, 0.0}, {2.0 + shift, 2.0}, flow, &open2D};
            const shapes::Circle circle{{0.25 + shift, 0.25}, 1.25, flow, &open2D};
            for(double h : {0.5, 0.37, 4.0})
            {
                baselineCase<false, lattice::LatticeType::D2Q9>(rectangle, h, "Rectangle");
                baselineCase<false, lattice::LatticeType::D2Q9>(circle, h, "Circle");
            }
            baseline3D<lattice::LatticeType::D3Q15>(shift, flow, open3D);
            baseline3D<lattice::LatticeType::D3Q19>(shift, flow, open3D);
            baseline3D<lattice::LatticeType::D3Q27>(shift, flow, open3D);
        }
    }
    require(caseNumber == primitiveBaseline.size(), "Historical case count differs");
    std::cout << "Historical CellType/Q regression and repeated 1/4-thread output: "
              << caseNumber << " cases passed\n";
}

template<bool ThreeD, lattice::LatticeType Type, typename Geometry, typename Position>
void probe(const Geometry& geometry, const Position& p, const Position& delta,
           std::size_t nFaces, std::size_t expectedFace)
{
    using Analysis = std::conditional_t<ThreeD, GeometryAnalysis3D<Type>, GeometryAnalysis2D<Type>>;
    using Lattice = typename Analysis::Lattice;
    Analysis analysis(0.5);
    analysis.nBoundaries = nFaces;
    analysis.boundaryOffsets.resize(nFaces + 1);
    for(std::size_t face = 0; face <= nFaces; ++face) analysis.boundaryOffsets[face] = face;
    analysis.boundaryX.assign(nFaces, p[0]);
    analysis.boundaryY.assign(nFaces, p[1]);
    if constexpr(ThreeD) analysis.boundaryZ.assign(nFaces, p[2]);
    geometry.boundaryQAnalysis(analysis);
    std::size_t direction = 0;
    for(std::size_t i = 1; i < Lattice::nStencils; ++i)
    {
        bool match = Lattice::ex[i] == delta[0] && Lattice::ey[i] == delta[1];
        if constexpr(ThreeD) match = match && Lattice::ez[i] == delta[2];
        if(match) direction = i;
    }
    require(direction != 0, "Probe direction missing");
    for(std::size_t face = 0; face < nFaces; ++face)
    {
        const double q = analysis.q[face * (Lattice::nStencils - 1) + direction - 1];
        require(face == expectedFace ? q > 0.0 && q <= 1.0 : q == -1.0,
                "Exact/perturbed shared-hit ownership failed");
    }
}

void sharedHitPerturbations()
{
    const Rectangle rectangle{{0.0, 0.0}, {2.0, 2.0}};
    const Box box{{0.0, 0.0, 0.0}, {2.0, 2.0, 2.0}};
    const Cylinder cylinder{{0.0, 0.0, 0.0}, 1.25, 2.5, Axis::Z};
    const double epsilon = std::ldexp(1.0, -36);
    for(double offset : {-epsilon, 0.0, epsilon})
    {
        probe<false, lattice::LatticeType::D2Q9>(rectangle,
            Point2D{0.25, 0.25 + offset}, Point2D{-1.0, -1.0}, 4, offset > 0.0 ? 3 : 0);
        probe<true, lattice::LatticeType::D3Q27>(box,
            Point{0.25, 0.25 + offset, 1.0}, Point{-1.0, -1.0, -1.0}, 6, offset < 0.0 ? 2 : 0);
        probe<true, lattice::LatticeType::D3Q19>(cylinder,
            Point{0.5, 1.0, 1.0 + offset}, Point{1.0, 0.0, 1.0}, 3, offset > 0.0 ? 2 : 0);
    }
    std::cout << "Exact shared-hit and two-sided perturbation checks passed\n";
}

void multiBlockCircle()
{
    const BoundingBox2D domain{{-2.0, -2.0}, {3.0, 3.0}};
    const Circle circle{{0.25, 0.25}, 1.25, FlowType::External, &domain};
    GeometryAnalysis2D<lattice::LatticeType::D2Q9> analysis(0.05);
    omp_set_num_threads(1);
    single_boundary_suite::checkSingleBoundary<false>(circle, analysis, false);
    require(analysis.scalar.size() > 2 * 4096, "Circle fixture must span multiple blocks");
    repeat<false>(circle, analysis);
    std::cout << "Multi-block Circle compaction checks passed\n";
}
} // namespace acceptance_audit

int main()
{
    omp_set_dynamic(0);
    try
    {
        if(single_boundary_suite::run() != 0 || rectangle_suite::run() != 0 ||
           box_cylinder_suite::run() != 0) return 1;
        acceptance_audit::sharedHitPerturbations();
        acceptance_audit::multiBlockCircle();
        acceptance_audit::historicalRegression();
        std::cout << "Primitive stage acceptance passed\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "Primitive stage acceptance failed: " << error.what() << '\n';
        return 1;
    }
}
