#include "quantity.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace ntic::lbm::config
{

//==============================================================
// Constructors
//==============================================================

Quantity::Quantity()
    :
    value_(0.0),
    unit_()
{
}

Quantity::Quantity(
    double value,
    const UnitExpression& unit)
    :
    value_(value),
    unit_(unit)
{
    if (!std::isfinite(value_))
    {
        throw std::invalid_argument(
            "Quantity: value must be finite.");
    }
}

Quantity
Quantity::Dimensionless(
    double value)
{
    return Quantity(
        value,
        UnitExpression());
}

//==============================================================
// Access
//==============================================================

double
Quantity::value() const
{
    return value_;
}

double
Quantity::valueSI() const
{
    return value_ * unit_.scale();
}

double
Quantity::scale() const
{
    return unit_.scale();
}

const UnitExpression&
Quantity::unit() const
{
    return unit_;
}

//==============================================================
// Algebra
//==============================================================

Quantity
Quantity::operator*(
    const Quantity& rhs) const
{
    return Quantity(
        value_ * rhs.value_,
        unit_ * rhs.unit_);
}

Quantity
Quantity::operator/(
    const Quantity& rhs) const
{
    if (rhs.value_ == 0.0)
    {
        throw std::domain_error(
            "Quantity: division by zero.");
    }

    return Quantity(
        value_ / rhs.value_,
        unit_ / rhs.unit_);
}

//==============================================================
// Comparison
//==============================================================

bool
Quantity::operator==(
    const Quantity& rhs) const
{
    return value_ == rhs.value_
        && unit_ == rhs.unit_;
}

bool
Quantity::operator!=(
    const Quantity& rhs) const
{
    return !(*this == rhs);
}

bool
Quantity::equivalentTo(
    const Quantity& rhs,
    double relativeTolerance) const
{
    if (unit_.dimension() != rhs.unit_.dimension())
    {
        return false;
    }

    const double lhsValue =
        valueSI();

    const double rhsValue =
        rhs.valueSI();

    const double scale =
        std::max(
            1.0,
            std::max(
                std::abs(lhsValue),
                std::abs(rhsValue)));

    return std::abs(
               lhsValue-rhsValue)
        <= relativeTolerance*scale;
}

//==============================================================
// Utility
//==============================================================

bool
Quantity::isDimensionless() const
{
    return unit_.isDimensionless();
}

bool
Quantity::isStandardSI() const
{
    return unit_.isStandardSI();
}

std::string
Quantity::toString() const
{
    std::ostringstream oss;

    oss
        << value_
        << " "
        << unit_.toString();

    return oss.str();
}

}
