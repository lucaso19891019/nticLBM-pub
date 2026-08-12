#include "unit_expression.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace ntic::lbm::config
{

//==============================================================
// Constructors
//==============================================================

UnitExpression::UnitExpression()
    : dimension_(PhysicalDimension::Dimensionless()),
      scale_(1.0)
{
}

UnitExpression::UnitExpression(
    const Unit& unit)
    : dimension_(unit.dimension()),
      scale_(unit.scale())
{
}

UnitExpression::UnitExpression(
    const PhysicalDimension& dimension,
    double scale)
    : dimension_(dimension),
      scale_(scale)
{
    if (!std::isfinite(scale_) || scale_ <= 0.0)
    {
        throw std::invalid_argument(
            "UnitExpression: scale must be finite "
            "and greater than zero.");
    }
}

//==============================================================
// Access
//==============================================================

const PhysicalDimension&
UnitExpression::dimension() const
{
    return dimension_;
}

double
UnitExpression::scale() const
{
    return scale_;
}

//==============================================================
// Algebra
//==============================================================

UnitExpression
UnitExpression::operator*(
    const UnitExpression& rhs) const
{
    const double resultScale =
        scale_ * rhs.scale_;

    if (!std::isfinite(resultScale)
        || resultScale <= 0.0)
    {
        throw std::overflow_error(
            "UnitExpression: invalid scale produced "
            "during multiplication.");
    }

    return UnitExpression(
        dimension_ * rhs.dimension_,
        resultScale);
}

UnitExpression
UnitExpression::operator/(
    const UnitExpression& rhs) const
{
    if (!std::isfinite(rhs.scale_)
        || rhs.scale_ <= 0.0)
    {
        throw std::domain_error(
            "UnitExpression: division by an expression "
            "with an invalid scale.");
    }

    const double resultScale =
        scale_ / rhs.scale_;

    if (!std::isfinite(resultScale)
        || resultScale <= 0.0)
    {
        throw std::overflow_error(
            "UnitExpression: invalid scale produced "
            "during division.");
    }

    return UnitExpression(
        dimension_ / rhs.dimension_,
        resultScale);
}

UnitExpression
UnitExpression::pow(
    const common::Rational& exponent) const
{
    if (exponent.isZero())
    {
        return UnitExpression();
    }

    const double resultScale =
        std::pow(
            scale_,
            exponent.toDouble());

    if (!std::isfinite(resultScale)
        || resultScale <= 0.0)
    {
        throw std::domain_error(
            "UnitExpression: invalid scale produced "
            "during exponentiation.");
    }

    return UnitExpression(
        dimension_.pow(exponent),
        resultScale);
}

//==============================================================
// Comparison
//==============================================================

bool
UnitExpression::operator==(
    const UnitExpression& rhs) const
{
    return dimension_ == rhs.dimension_
        && scale_ == rhs.scale_;
}

bool
UnitExpression::operator!=(
    const UnitExpression& rhs) const
{
    return !(*this == rhs);
}

bool
UnitExpression::equivalentTo(
    const UnitExpression& rhs,
    double relativeTolerance) const
{
    if (dimension_ != rhs.dimension_)
    {
        return false;
    }

    if (!std::isfinite(relativeTolerance)
        || relativeTolerance < 0.0)
    {
        throw std::invalid_argument(
            "UnitExpression::equivalentTo: relative tolerance "
            "must be finite and non-negative.");
    }

    const double comparisonScale =
        std::max(
            {1.0,
             std::abs(scale_),
             std::abs(rhs.scale_)});

    return std::abs(scale_ - rhs.scale_)
        <= relativeTolerance * comparisonScale;
}

//==============================================================
// Utility
//==============================================================

bool
UnitExpression::isDimensionless() const
{
    return dimension_.isDimensionless();
}

bool
UnitExpression::isStandardSI() const
{
    return scale_ == 1.0;
}

std::string
UnitExpression::toString() const
{
    std::ostringstream oss;

    oss << dimension_.toString()
        << " scale="
        << scale_;

    return oss.str();
}

} // namespace ntic::lbm::config
