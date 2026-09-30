#include "stl_geometry.hpp"

#include <algorithm>
#include <atomic>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <omp.h>


namespace ntic::lbm::geometry
{

namespace
{

//=============================================================================
// Vector operations
//=============================================================================

Point subtract(
    const Point& a,
    const Point& b)
{
    return
    {
        a[0] - b[0],
        a[1] - b[1],
        a[2] - b[2]
    };
}


double dot(
    const Point& a,
    const Point& b)
{
    return
        a[0] * b[0] +
        a[1] * b[1] +
        a[2] * b[2];
}


Point cross(
    const Point& a,
    const Point& b)
{
    return
    {
        a[1] * b[2] -
            a[2] * b[1],

        a[2] * b[0] -
            a[0] * b[2],

        a[0] * b[1] -
            a[1] * b[0]
    };
}


double norm(
    const Point& a)
{
    return
        std::sqrt(
            dot(
                a,
                a));
}


//=============================================================================
// Triangle solid angle
//=============================================================================

double triangleSolidAngle(
    const Point& point,
    const Point& vertex0,
    const Point& vertex1,
    const Point& vertex2)
{
    const Point a =
        subtract(
            vertex0,
            point);

    const Point b =
        subtract(
            vertex1,
            point);

    const Point c =
        subtract(
            vertex2,
            point);


    const double lengthA =
        norm(
            a);

    const double lengthB =
        norm(
            b);

    const double lengthC =
        norm(
            c);


    const double numerator =
        dot(
            a,
            cross(
                b,
                c));


    const double denominator =
        lengthA *
            lengthB *
            lengthC +
        dot(
            a,
            b) *
            lengthC +
        dot(
            b,
            c) *
            lengthA +
        dot(
            c,
            a) *
            lengthB;


    return
        2.0 *
        std::atan2(
            numerator,
            denominator);
}


//=============================================================================
// Point in component
//=============================================================================

bool pointInComponent(
    const STLGeometry& geometry,
    const std::size_t componentID,
    const Point& point)
{
    const auto& topology =
        geometry.topology;

    const auto& component =
        topology.components[
            componentID];


    double solidAngle =
        0.0;


#pragma omp parallel for reduction(+:solidAngle) schedule(static) \
    if(!omp_in_parallel())
    for(std::int64_t localFacetID = 0;
        localFacetID <
            static_cast<std::int64_t>(
                component.facets.size());
        ++localFacetID)
    {
        const std::size_t facetID =
            component.facets[
                static_cast<std::size_t>(
                    localFacetID)];


        const auto& vertexIDs =
            topology.geometry.facetVertexIDs[
                facetID];


        const Point& vertex0 =
            topology.geometry.vertices[
                vertexIDs[0]];

        const Point& vertex1 =
            topology.geometry.vertices[
                vertexIDs[1]];

        const Point& vertex2 =
            topology.geometry.vertices[
                vertexIDs[2]];


        solidAngle +=
            triangleSolidAngle(
                point,
                vertex0,
                vertex1,
                vertex2);
    }


    constexpr double pi =
        3.14159265358979323846;


    return
        std::abs(
            solidAngle) >
        2.0 * pi;
}


//=============================================================================
// Wet-side test
//=============================================================================

bool centerIsWet(
    const STLGeometry& geometry,
    const Point& point)
{
    const std::size_t componentCount =
        geometry.topology.components.size();


    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        const STLComponentFlow& componentFlow =
            geometry.flow[
                componentID];


        if(!componentFlow.active)
        {
            continue;
        }


        const auto& component =
            geometry.topology.components[
                componentID];


        const bool outsideBounds =
            point[0] <
                component.bounds.min[0] ||
            point[0] >
                component.bounds.max[0] ||
            point[1] <
                component.bounds.min[1] ||
            point[1] >
                component.bounds.max[1] ||
            point[2] <
                component.bounds.min[2] ||
            point[2] >
                component.bounds.max[2];


        const bool inside =
            outsideBounds
                ? false
                : pointInComponent(
                    geometry,
                    componentID,
                    point);


        if(componentFlow.fluidSide ==
           FluidSide::Inside)
        {
            if(!inside)
            {
                return false;
            }
        }
        else
        {
            if(inside)
            {
                return false;
            }
        }
    }


    return true;
}


//=============================================================================
// Separating-axis tests
//=============================================================================

bool separatedOnAxisClosed(
    const Point& vertex0,
    const Point& vertex1,
    const Point& vertex2,
    const Point& axis,
    const double halfGridSpacing)
{
    const double axisLengthSquared =
        dot(
            axis,
            axis);


    if(axisLengthSquared ==
       0.0)
    {
        return false;
    }


    const double projection0 =
        dot(
            vertex0,
            axis);

    const double projection1 =
        dot(
            vertex1,
            axis);

    const double projection2 =
        dot(
            vertex2,
            axis);


    const double triangleMin =
        std::min(
            projection0,
            std::min(
                projection1,
                projection2));

    const double triangleMax =
        std::max(
            projection0,
            std::max(
                projection1,
                projection2));


    const double boxRadius =
        halfGridSpacing *
        (
            std::abs(
                axis[0]) +
            std::abs(
                axis[1]) +
            std::abs(
                axis[2])
        );


    return
        triangleMax <
            -boxRadius ||
        triangleMin >
            boxRadius;
}


bool separatedOnAxisOpen(
    const Point& vertex0,
    const Point& vertex1,
    const Point& vertex2,
    const Point& axis,
    const double halfGridSpacing)
{
    const double axisLengthSquared =
        dot(
            axis,
            axis);


    if(axisLengthSquared ==
       0.0)
    {
        return false;
    }


    const double projection0 =
        dot(
            vertex0,
            axis);

    const double projection1 =
        dot(
            vertex1,
            axis);

    const double projection2 =
        dot(
            vertex2,
            axis);


    const double triangleMin =
        std::min(
            projection0,
            std::min(
                projection1,
                projection2));

    const double triangleMax =
        std::max(
            projection0,
            std::max(
                projection1,
                projection2));


    const double boxRadius =
        halfGridSpacing *
        (
            std::abs(
                axis[0]) +
            std::abs(
                axis[1]) +
            std::abs(
                axis[2])
        );


    return
        triangleMax <=
            -boxRadius ||
        triangleMin >=
            boxRadius;
}


//=============================================================================
// Triangle / closed-cell intersection
//=============================================================================

bool triangleIntersectsCell(
    const Point& point,
    const double gridSpacing,
    const Point& triangleVertex0,
    const Point& triangleVertex1,
    const Point& triangleVertex2)
{
    const double halfGridSpacing =
        0.5 *
        gridSpacing;


    const Point vertex0 =
        subtract(
            triangleVertex0,
            point);

    const Point vertex1 =
        subtract(
            triangleVertex1,
            point);

    const Point vertex2 =
        subtract(
            triangleVertex2,
            point);


    const double triangleMinX =
        std::min(
            vertex0[0],
            std::min(
                vertex1[0],
                vertex2[0]));

    const double triangleMaxX =
        std::max(
            vertex0[0],
            std::max(
                vertex1[0],
                vertex2[0]));


    if(triangleMaxX <
           -halfGridSpacing ||
       triangleMinX >
           halfGridSpacing)
    {
        return false;
    }


    const double triangleMinY =
        std::min(
            vertex0[1],
            std::min(
                vertex1[1],
                vertex2[1]));

    const double triangleMaxY =
        std::max(
            vertex0[1],
            std::max(
                vertex1[1],
                vertex2[1]));


    if(triangleMaxY <
           -halfGridSpacing ||
       triangleMinY >
           halfGridSpacing)
    {
        return false;
    }


    const double triangleMinZ =
        std::min(
            vertex0[2],
            std::min(
                vertex1[2],
                vertex2[2]));

    const double triangleMaxZ =
        std::max(
            vertex0[2],
            std::max(
                vertex1[2],
                vertex2[2]));


    if(triangleMaxZ <
           -halfGridSpacing ||
       triangleMinZ >
           halfGridSpacing)
    {
        return false;
    }


    const Point edge0 =
        subtract(
            vertex1,
            vertex0);

    const Point edge1 =
        subtract(
            vertex2,
            vertex1);

    const Point edge2 =
        subtract(
            vertex0,
            vertex2);


    const Point triangleNormal =
        cross(
            edge0,
            edge1);


    if(separatedOnAxisClosed(
           vertex0,
           vertex1,
           vertex2,
           triangleNormal,
           halfGridSpacing))
    {
        return false;
    }


    const std::array<Point,3> edges =
    {{
        edge0,
        edge1,
        edge2
    }};


    const std::array<Point,3> cellAxes =
    {{
        {
            1.0,
            0.0,
            0.0
        },
        {
            0.0,
            1.0,
            0.0
        },
        {
            0.0,
            0.0,
            1.0
        }
    }};


    for(const Point& edge :
        edges)
    {
        for(const Point& cellAxis :
            cellAxes)
        {
            const Point axis =
                cross(
                    edge,
                    cellAxis);


            if(separatedOnAxisClosed(
                   vertex0,
                   vertex1,
                   vertex2,
                   axis,
                   halfGridSpacing))
            {
                return false;
            }
        }
    }


    return true;
}


//=============================================================================
// Triangle / open-cell-interior intersection
//=============================================================================

bool triangleIntersectsCellInterior(
    const Point& point,
    const double gridSpacing,
    const Point& triangleVertex0,
    const Point& triangleVertex1,
    const Point& triangleVertex2)
{
    const double halfGridSpacing =
        0.5 *
        gridSpacing;


    const Point vertex0 =
        subtract(
            triangleVertex0,
            point);

    const Point vertex1 =
        subtract(
            triangleVertex1,
            point);

    const Point vertex2 =
        subtract(
            triangleVertex2,
            point);


    const double triangleMinX =
        std::min(
            vertex0[0],
            std::min(
                vertex1[0],
                vertex2[0]));

    const double triangleMaxX =
        std::max(
            vertex0[0],
            std::max(
                vertex1[0],
                vertex2[0]));


    if(triangleMaxX <=
           -halfGridSpacing ||
       triangleMinX >=
           halfGridSpacing)
    {
        return false;
    }


    const double triangleMinY =
        std::min(
            vertex0[1],
            std::min(
                vertex1[1],
                vertex2[1]));

    const double triangleMaxY =
        std::max(
            vertex0[1],
            std::max(
                vertex1[1],
                vertex2[1]));


    if(triangleMaxY <=
           -halfGridSpacing ||
       triangleMinY >=
           halfGridSpacing)
    {
        return false;
    }


    const double triangleMinZ =
        std::min(
            vertex0[2],
            std::min(
                vertex1[2],
                vertex2[2]));

    const double triangleMaxZ =
        std::max(
            vertex0[2],
            std::max(
                vertex1[2],
                vertex2[2]));


    if(triangleMaxZ <=
           -halfGridSpacing ||
       triangleMinZ >=
           halfGridSpacing)
    {
        return false;
    }


    const Point edge0 =
        subtract(
            vertex1,
            vertex0);

    const Point edge1 =
        subtract(
            vertex2,
            vertex1);

    const Point edge2 =
        subtract(
            vertex0,
            vertex2);


    const Point triangleNormal =
        cross(
            edge0,
            edge1);


    if(separatedOnAxisOpen(
           vertex0,
           vertex1,
           vertex2,
           triangleNormal,
           halfGridSpacing))
    {
        return false;
    }


    const std::array<Point,3> edges =
    {{
        edge0,
        edge1,
        edge2
    }};


    const std::array<Point,3> cellAxes =
    {{
        {
            1.0,
            0.0,
            0.0
        },
        {
            0.0,
            1.0,
            0.0
        },
        {
            0.0,
            0.0,
            1.0
        }
    }};


    for(const Point& edge :
        edges)
    {
        for(const Point& cellAxis :
            cellAxes)
        {
            const Point axis =
                cross(
                    edge,
                    cellAxis);


            if(separatedOnAxisOpen(
                   vertex0,
                   vertex1,
                   vertex2,
                   axis,
                   halfGridSpacing))
            {
                return false;
            }
        }
    }


    return true;
}


//=============================================================================
// Surface-cell state
//=============================================================================

enum class SurfaceCellState : std::uint8_t
{
    None = 0,
    Touch = 1,
    Boundary = 2
};


void upgradeSurfaceCellState(
    std::atomic<std::uint8_t>& state,
    const SurfaceCellState requestedState)
{
    const std::uint8_t requested =
        static_cast<std::uint8_t>(
            requestedState);

    std::uint8_t current =
        state.load(
            std::memory_order_relaxed);


    while(current < requested &&
          !state.compare_exchange_weak(
              current,
              requested,
              std::memory_order_relaxed,
              std::memory_order_relaxed))
    {
    }
}


bool facetCellIndexRange(
    const double facetMin,
    const double facetMax,
    const double gridMin,
    const double gridSpacing,
    const std::size_t cellCount,
    std::size_t& firstCell,
    std::size_t& lastCell)
{
    const double relativeMin =
        (facetMin - gridMin) /
        gridSpacing;

    const double relativeMax =
        (facetMax - gridMin) /
        gridSpacing;


    std::ptrdiff_t first =
        static_cast<std::ptrdiff_t>(
            std::ceil(
                relativeMin)) -
        1;

    std::ptrdiff_t last =
        static_cast<std::ptrdiff_t>(
            std::floor(
                relativeMax));


    if(last < 0 ||
       first >=
           static_cast<std::ptrdiff_t>(
               cellCount))
    {
        return false;
    }


    first =
        std::max<std::ptrdiff_t>(
            first,
            0);

    last =
        std::min<std::ptrdiff_t>(
            last,
            static_cast<std::ptrdiff_t>(
                cellCount) -
            1);


    if(first > last)
    {
        return false;
    }


    firstCell =
        static_cast<std::size_t>(
            first);

    lastCell =
        static_cast<std::size_t>(
            last);


    return true;
}


Point cellCenter(
    const BoundingBox& domain,
    const double gridSpacing,
    const std::size_t i,
    const std::size_t j,
    const std::size_t k)
{
    return
    {
        domain.min[0] +
            (static_cast<double>(i) + 0.5) *
            gridSpacing,

        domain.min[1] +
            (static_cast<double>(j) + 0.5) *
            gridSpacing,

        domain.min[2] +
            (static_cast<double>(k) + 0.5) *
            gridSpacing
    };
}

} // namespace


//=============================================================================
// Interior area analysis
//=============================================================================

void STLGeometry::interiorAreaAnalysis(
    const double gridSpacing,
    BoundingBox& domain,
    std::size_t& nx,
    std::size_t& ny,
    std::size_t& nz,
    std::vector<double>& scalar,
    std::vector<double>& boundaryX,
    std::vector<double>& boundaryY,
    std::vector<double>& boundaryZ) const
{
    if(gridSpacing <=
       0.0)
    {
        throw std::runtime_error(
            "Grid spacing must be positive.");
    }


    if(flow.size() !=
       topology.components.size())
    {
        throw std::runtime_error(
            "STL flow interpretation is not available.");
    }


    domain =
        flowType == FlowType::Internal
            ? bounds
            : openBox;


    nx =
        static_cast<std::size_t>(
            std::ceil(
                domain.width() /
                gridSpacing));

    ny =
        static_cast<std::size_t>(
            std::ceil(
                domain.height() /
                gridSpacing));

    nz =
        static_cast<std::size_t>(
            std::ceil(
                domain.depth() /
                gridSpacing));


    if(nx == 0 ||
       ny == 0 ||
       nz == 0)
    {
        throw std::runtime_error(
            "STL analysis domain contains no grid cells.");
    }


    const std::size_t xySize =
        nx * ny;

    const std::size_t cellCount =
        xySize * nz;


    std::vector<CellType> cellTypes;


    cellTypes.assign(
        cellCount,
        CellType::Dry);

    boundaryX.clear();
    boundaryY.clear();
    boundaryZ.clear();


    //--------------------------------------------------------------------------
    // Step 1:
    //
    // Rasterize active STL facets directly onto the structured grid.
    //
    // For each facet, its exact axis-aligned bounding box is converted
    // algebraically to the range of grid cells whose closed boxes can touch
    // that facet. Only those cells are tested.
    //
    // SurfaceCellState::Touch:
    //     at least one active facet intersects the closed cell, but no active
    //     facet has yet been found to enter the open cell interior.
    //
    // SurfaceCellState::Boundary:
    //     at least one active facet enters the open cell interior. This is a
    //     definite final boundary cell.
    //
    // Multiple facets may affect the same cell. The state therefore changes
    // monotonically:
    //
    //     None -> Touch -> Boundary
    //--------------------------------------------------------------------------

    std::vector<std::atomic<std::uint8_t>> atomicSurfaceStates(
        cellCount);

#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                cellCount);
        ++index)
    {
        atomicSurfaceStates[
            static_cast<std::size_t>(
                index)].store(
                    static_cast<std::uint8_t>(
                        SurfaceCellState::None),
                    std::memory_order_relaxed);
    }


    const std::size_t facetCount =
        topology.geometry.facetVertexIDs.size();


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                facetCount);
        ++index)
    {
        const std::size_t facetID =
            static_cast<std::size_t>(
                index);

        const std::size_t componentID =
            topology.facetComponentIDs[
                facetID];


        if(!flow[
                componentID].active)
        {
            continue;
        }


        const auto& vertexIDs =
            topology.geometry.facetVertexIDs[
                facetID];

        const Point& vertex0 =
            topology.geometry.vertices[
                vertexIDs[0]];

        const Point& vertex1 =
            topology.geometry.vertices[
                vertexIDs[1]];

        const Point& vertex2 =
            topology.geometry.vertices[
                vertexIDs[2]];


        const double facetMinX =
            std::min(
                vertex0[0],
                std::min(
                    vertex1[0],
                    vertex2[0]));

        const double facetMaxX =
            std::max(
                vertex0[0],
                std::max(
                    vertex1[0],
                    vertex2[0]));

        const double facetMinY =
            std::min(
                vertex0[1],
                std::min(
                    vertex1[1],
                    vertex2[1]));

        const double facetMaxY =
            std::max(
                vertex0[1],
                std::max(
                    vertex1[1],
                    vertex2[1]));

        const double facetMinZ =
            std::min(
                vertex0[2],
                std::min(
                    vertex1[2],
                    vertex2[2]));

        const double facetMaxZ =
            std::max(
                vertex0[2],
                std::max(
                    vertex1[2],
                    vertex2[2]));


        std::size_t firstI =
            0;

        std::size_t lastI =
            0;

        std::size_t firstJ =
            0;

        std::size_t lastJ =
            0;

        std::size_t firstK =
            0;

        std::size_t lastK =
            0;


        if(!facetCellIndexRange(
               facetMinX,
               facetMaxX,
               domain.min[0],
               gridSpacing,
               nx,
               firstI,
               lastI) ||
           !facetCellIndexRange(
               facetMinY,
               facetMaxY,
               domain.min[1],
               gridSpacing,
               ny,
               firstJ,
               lastJ) ||
           !facetCellIndexRange(
               facetMinZ,
               facetMaxZ,
               domain.min[2],
               gridSpacing,
               nz,
               firstK,
               lastK))
        {
            continue;
        }


        for(std::size_t k = firstK;
            k <= lastK;
            ++k)
        {
            for(std::size_t j = firstJ;
                j <= lastJ;
                ++j)
            {
                for(std::size_t i = firstI;
                    i <= lastI;
                    ++i)
                {
                    const std::size_t cellID =
                        i +
                        nx *
                        (
                            j +
                            ny * k
                        );


                    if(atomicSurfaceStates[
                           cellID].load(
                               std::memory_order_relaxed) ==
                       static_cast<std::uint8_t>(
                           SurfaceCellState::Boundary))
                    {
                        continue;
                    }


                    const Point point =
                        cellCenter(
                            domain,
                            gridSpacing,
                            i,
                            j,
                            k);


                    if(!triangleIntersectsCell(
                           point,
                           gridSpacing,
                           vertex0,
                           vertex1,
                           vertex2))
                    {
                        continue;
                    }


                    if(triangleIntersectsCellInterior(
                           point,
                           gridSpacing,
                           vertex0,
                           vertex1,
                           vertex2))
                    {
                        upgradeSurfaceCellState(
                            atomicSurfaceStates[
                                cellID],
                            SurfaceCellState::Boundary);
                    }
                    else
                    {
                        upgradeSurfaceCellState(
                            atomicSurfaceStates[
                                cellID],
                            SurfaceCellState::Touch);
                    }
                }
            }
        }
    }


    std::vector<SurfaceCellState> surfaceStates(
        cellCount,
        SurfaceCellState::None);


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                cellCount);
        ++index)
    {
        const std::size_t cellID =
            static_cast<std::size_t>(
                index);

        surfaceStates[cellID] =
            static_cast<SurfaceCellState>(
                atomicSurfaceStates[
                    cellID].load(
                        std::memory_order_relaxed));


        if(surfaceStates[cellID] ==
           SurfaceCellState::Boundary)
        {
            cellTypes[cellID] =
                CellType::Boundary;
        }
    }


    //--------------------------------------------------------------------------
    // Step 2:
    //
    // Classify all non-Boundary cells by 6-neighbor connected regions.
    //
    // Touch cells participate in region analysis because the STL does not
    // enter their open cell interior. Definite Boundary cells do not.
    //
    // A conservative rule is used for Touch-Touch neighbors: they are not
    // connected directly. If an STL facet lies exactly on their shared grid
    // face, the two cell interiors can be on opposite physical sides even
    // though both cells are only Touch cells. Refusing direct Touch-Touch
    // propagation can split one physical region into several numerical
    // regions, but that does not change classification because every numerical
    // region is independently classified by one strict centerIsWet() query.
    //--------------------------------------------------------------------------

    std::vector<std::uint8_t> visited(
        cellCount,
        0);

    std::vector<std::size_t> queue;

    queue.reserve(
        cellCount);


    const std::array<std::array<int,3>,6> neighborOffsets =
    {{
        {{-1,  0,  0}},
        {{ 1,  0,  0}},
        {{ 0, -1,  0}},
        {{ 0,  1,  0}},
        {{ 0,  0, -1}},
        {{ 0,  0,  1}}
    }};


    for(std::size_t seedID = 0;
        seedID < cellCount;
        ++seedID)
    {
        if(surfaceStates[seedID] ==
           SurfaceCellState::Boundary)
        {
            continue;
        }

        if(visited[seedID] !=
           0)
        {
            continue;
        }


        const std::size_t seedK =
            seedID /
            xySize;

        const std::size_t seedRemainder =
            seedID %
            xySize;

        const std::size_t seedJ =
            seedRemainder /
            nx;

        const std::size_t seedI =
            seedRemainder %
            nx;


        const Point seedPoint =
            cellCenter(
                domain,
                gridSpacing,
                seedI,
                seedJ,
                seedK);


        const bool wet =
            centerIsWet(
                *this,
                seedPoint);

        const CellType noneType =
            wet
                ? CellType::Interior
                : CellType::Dry;

        const CellType touchType =
            wet
                ? CellType::Boundary
                : CellType::Dry;


        queue.clear();

        queue.push_back(
            seedID);

        visited[seedID] =
            1;


        std::size_t head =
            0;


        while(head <
              queue.size())
        {
            const std::size_t cellID =
                queue[
                    head++];


            cellTypes[cellID] =
                surfaceStates[cellID] ==
                    SurfaceCellState::Touch
                    ? touchType
                    : noneType;


            const std::size_t k =
                cellID /
                xySize;

            const std::size_t remainder =
                cellID %
                xySize;

            const std::size_t j =
                remainder /
                nx;

            const std::size_t i =
                remainder %
                nx;


            for(const auto& offset :
                neighborOffsets)
            {
                const std::ptrdiff_t ni =
                    static_cast<std::ptrdiff_t>(
                        i) +
                    offset[0];

                const std::ptrdiff_t nj =
                    static_cast<std::ptrdiff_t>(
                        j) +
                    offset[1];

                const std::ptrdiff_t nk =
                    static_cast<std::ptrdiff_t>(
                        k) +
                    offset[2];


                if(ni < 0 ||
                   nj < 0 ||
                   nk < 0 ||
                   ni >=
                       static_cast<std::ptrdiff_t>(
                           nx) ||
                   nj >=
                       static_cast<std::ptrdiff_t>(
                           ny) ||
                   nk >=
                       static_cast<std::ptrdiff_t>(
                           nz))
                {
                    continue;
                }


                const std::size_t neighborID =
                    static_cast<std::size_t>(
                        ni) +
                    nx *
                    (
                        static_cast<std::size_t>(
                            nj) +
                        ny *
                        static_cast<std::size_t>(
                            nk)
                    );


                if(visited[neighborID] !=
                   0)
                {
                    continue;
                }

                if(surfaceStates[neighborID] ==
                   SurfaceCellState::Boundary)
                {
                    continue;
                }


                if(surfaceStates[cellID] ==
                       SurfaceCellState::Touch &&
                   surfaceStates[neighborID] ==
                       SurfaceCellState::Touch)
                {
                    continue;
                }


                visited[neighborID] =
                    1;

                queue.push_back(
                    neighborID);
            }
        }
    }


    //--------------------------------------------------------------------------
    // Step 3:
    //
    // Compact final boundary-cell coordinates into SoA vectors.
    //--------------------------------------------------------------------------

    boundaryX.resize(
        cellCount);

    boundaryY.resize(
        cellCount);

    boundaryZ.resize(
        cellCount);


    std::size_t boundaryCount =
        0;


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                cellCount);
        ++index)
    {
        const std::size_t cellID =
            static_cast<std::size_t>(
                index);


        if(cellTypes[cellID] !=
           CellType::Boundary)
        {
            continue;
        }


        const std::size_t k =
            cellID /
            xySize;

        const std::size_t remainder =
            cellID %
            xySize;

        const std::size_t j =
            remainder /
            nx;

        const std::size_t i =
            remainder %
            nx;


        std::size_t slot =
            0;

#pragma omp atomic capture
        slot = boundaryCount++;


        const Point point =
            cellCenter(
                domain,
                gridSpacing,
                i,
                j,
                k);


        boundaryX[slot] =
            point[0];

        boundaryY[slot] =
            point[1];

        boundaryZ[slot] =
            point[2];
    }


    boundaryX.resize(
        boundaryCount);

    boundaryY.resize(
        boundaryCount);

    boundaryZ.resize(
        boundaryCount);


    scalar.resize(
        cellCount);


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                cellCount);
        ++index)
    {
        const std::size_t cellID =
            static_cast<std::size_t>(
                index);

        scalar[cellID] =
            static_cast<double>(
                cellTypes[cellID]);
    }
}

} // namespace ntic::lbm::geometry
