#pragma once

#include "flow_type.hpp"

#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <utility>


namespace ntic::lbm::geometry
{

enum class BooleanOperation
{
    Union,
    Intersection,
    Difference
};


template <
    typename LeftGeometry,
    typename RightGeometry>
struct BooleanGeometry
{
    using BoundingBoxType =
        std::decay_t<
            decltype(
                std::declval<LeftGeometry>().bounds)>;

    using RightBoundingBoxType =
        std::decay_t<
            decltype(
                std::declval<RightGeometry>().bounds)>;

    static_assert(
        std::is_same_v<
            BoundingBoxType,
            RightBoundingBoxType>,
        "Boolean geometry operands must have "
        "the same bounding-box type.");


    LeftGeometry left;
    RightGeometry right;

    BooleanOperation operation;
    FlowType flowType;

    BoundingBoxType bounds;
    BoundingBoxType openBox;


    BooleanGeometry(
        LeftGeometry leftGeometry,
        RightGeometry rightGeometry,
        const BooleanOperation booleanOperation)
        :
        left(
            std::move(leftGeometry)),
        right(
            std::move(rightGeometry)),
        operation(
            booleanOperation),
        flowType(
            determineFlowType(
                left.flowType,
                right.flowType,
                operation)),
        bounds(
            determineBounds(
                left,
                right,
                operation)),
        openBox(
            determineOpenBox(
                left,
                right,
                operation))
    {
    }


private:

    static FlowType determineFlowType(
        const FlowType leftFlowType,
        const FlowType rightFlowType,
        const BooleanOperation booleanOperation)
    {
        if(booleanOperation ==
           BooleanOperation::Union)
        {
            if(leftFlowType !=
                   FlowType::Internal ||
               rightFlowType !=
                   FlowType::Internal)
            {
                throw std::invalid_argument(
                    "Geometry union requires two "
                    "internal-flow geometries.");
            }

            return FlowType::Internal;
        }


        if(booleanOperation ==
           BooleanOperation::Difference)
        {
            if(leftFlowType !=
                   FlowType::Internal ||
               rightFlowType !=
                   FlowType::Internal)
            {
                throw std::invalid_argument(
                    "Geometry difference requires two "
                    "internal-flow geometries.");
            }

            return FlowType::Internal;
        }


        if(booleanOperation ==
           BooleanOperation::Intersection)
        {
            if(leftFlowType ==
                   FlowType::External &&
               rightFlowType ==
                   FlowType::External)
            {
                return FlowType::External;
            }

            if(leftFlowType ==
                   FlowType::Internal &&
               rightFlowType ==
                   FlowType::External)
            {
                return FlowType::Internal;
            }

            if(leftFlowType ==
                   FlowType::External &&
               rightFlowType ==
                   FlowType::Internal)
            {
                return FlowType::Internal;
            }

            throw std::invalid_argument(
                "Geometry intersection requires "
                "external-external or "
                "internal-external geometries.");
        }


        throw std::invalid_argument(
            "Unsupported geometry Boolean operation.");
    }


    static BoundingBoxType unionBounds(
        const BoundingBoxType& leftBounds,
        const BoundingBoxType& rightBounds)
    {
        BoundingBoxType result =
            leftBounds;

        for(std::size_t d = 0;
            d < result.min.size();
            ++d)
        {
            result.min[d] =
                std::min(
                    leftBounds.min[d],
                    rightBounds.min[d]);

            result.max[d] =
                std::max(
                    leftBounds.max[d],
                    rightBounds.max[d]);
        }

        return result;
    }


    static bool sameBounds(
        const BoundingBoxType& leftBounds,
        const BoundingBoxType& rightBounds)
    {
        for(std::size_t d = 0;
            d < leftBounds.min.size();
            ++d)
        {
            if(leftBounds.min[d] !=
                   rightBounds.min[d] ||
               leftBounds.max[d] !=
                   rightBounds.max[d])
            {
                return false;
            }
        }

        return true;
    }


    static BoundingBoxType determineBounds(
        const LeftGeometry& leftGeometry,
        const RightGeometry& rightGeometry,
        const BooleanOperation booleanOperation)
    {
        if(booleanOperation ==
           BooleanOperation::Union)
        {
            return unionBounds(
                leftGeometry.bounds,
                rightGeometry.bounds);
        }


        if(booleanOperation ==
           BooleanOperation::Difference)
        {
            return leftGeometry.bounds;
        }


        if(booleanOperation ==
           BooleanOperation::Intersection)
        {
            if(leftGeometry.flowType ==
                   FlowType::External &&
               rightGeometry.flowType ==
                   FlowType::External)
            {
                return unionBounds(
                    leftGeometry.bounds,
                    rightGeometry.bounds);
            }

            if(leftGeometry.flowType ==
                   FlowType::Internal &&
               rightGeometry.flowType ==
                   FlowType::External)
            {
                return leftGeometry.bounds;
            }

            if(leftGeometry.flowType ==
                   FlowType::External &&
               rightGeometry.flowType ==
                   FlowType::Internal)
            {
                return rightGeometry.bounds;
            }
        }


        throw std::invalid_argument(
            "Cannot determine Boolean geometry bounds.");
    }


    static BoundingBoxType determineOpenBox(
        const LeftGeometry& leftGeometry,
        const RightGeometry& rightGeometry,
        const BooleanOperation booleanOperation)
    {
        if(booleanOperation ==
               BooleanOperation::Intersection &&
           leftGeometry.flowType ==
               FlowType::External &&
           rightGeometry.flowType ==
               FlowType::External)
        {
            if(!sameBounds(
                   leftGeometry.openBox,
                   rightGeometry.openBox))
            {
                throw std::invalid_argument(
                    "External geometry intersection "
                    "requires identical open boxes.");
            }

            return leftGeometry.openBox;
        }

        return BoundingBoxType{};
    }
};


template <
    typename LeftGeometry,
    typename RightGeometry>
auto operator|(
    LeftGeometry left,
    RightGeometry right)
{
    return BooleanGeometry<
        LeftGeometry,
        RightGeometry>(
            std::move(left),
            std::move(right),
            BooleanOperation::Union);
}


template <
    typename LeftGeometry,
    typename RightGeometry>
auto operator+(
    LeftGeometry left,
    RightGeometry right)
{
    return BooleanGeometry<
        LeftGeometry,
        RightGeometry>(
            std::move(left),
            std::move(right),
            BooleanOperation::Union);
}


template <
    typename LeftGeometry,
    typename RightGeometry>
auto operator&(
    LeftGeometry left,
    RightGeometry right)
{
    return BooleanGeometry<
        LeftGeometry,
        RightGeometry>(
            std::move(left),
            std::move(right),
            BooleanOperation::Intersection);
}


template <
    typename LeftGeometry,
    typename RightGeometry>
auto operator-(
    LeftGeometry left,
    RightGeometry right)
{
    return BooleanGeometry<
        LeftGeometry,
        RightGeometry>(
            std::move(left),
            std::move(right),
            BooleanOperation::Difference);
}

} // namespace ntic::lbm::geometry