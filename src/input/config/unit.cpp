#include "unit.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

using ntic::lbm::common::Rational;

namespace
{

constexpr double pi =
    3.141592653589793238462643383279502884;

}

namespace ntic::lbm::config
{

//==============================================================
// Constructors
//==============================================================

Unit::Unit()
    : symbol_(),
      dimension_(PhysicalDimension::Dimensionless()),
      scale_(1.0),
      prefixPolicy_(PrefixPolicy::Forbidden)
{
}

Unit::Unit(
    const std::string& symbol,
    const PhysicalDimension& dimension,
    double scale,
    PrefixPolicy prefixPolicy)
    : symbol_(symbol),
      dimension_(dimension),
      scale_(scale),
      prefixPolicy_(prefixPolicy)
{
    if (symbol_.empty())
    {
        throw std::invalid_argument(
            "Unit: symbol cannot be empty.");
    }

    if (!std::isfinite(scale_) || scale_ <= 0.0)
    {
        throw std::invalid_argument(
            "Unit '" + symbol_
            + "': scale must be finite and greater than zero.");
    }
}

//==============================================================
// Access
//==============================================================

const std::string&
Unit::symbol() const
{
    return symbol_;
}

const PhysicalDimension&
Unit::dimension() const
{
    return dimension_;
}

double
Unit::scale() const
{
    return scale_;
}

PrefixPolicy
Unit::prefixPolicy() const
{
    return prefixPolicy_;
}

bool
Unit::allowsPrefix() const
{
    return prefixPolicy_ == PrefixPolicy::Allowed;
}

//==============================================================
// Comparison
//==============================================================

bool
Unit::operator==(const Unit& rhs) const
{
    return symbol_ == rhs.symbol_
        && dimension_ == rhs.dimension_
        && scale_ == rhs.scale_
        && prefixPolicy_ == rhs.prefixPolicy_;
}

bool
Unit::operator!=(const Unit& rhs) const
{
    return !(*this == rhs);
}

//==============================================================
// Utility
//==============================================================

std::string
Unit::toString() const
{
    std::ostringstream oss;

    oss << symbol_
        << " "
        << dimension_.toString()
        << " scale="
        << scale_
        << " prefix=";

    if (allowsPrefix())
    {
        oss << "allowed";
    }
    else
    {
        oss << "forbidden";
    }

    return oss.str();
}

//==============================================================
// Built-in units
//==============================================================

Unit
Unit::Meter()
{
    return Unit(
        "m",
        PhysicalDimension::Length(),
        1.0,
        PrefixPolicy::Allowed);
}

Unit
Unit::Kilogram()
{
    return Unit(
        "kg",
        PhysicalDimension::Mass(),
        1.0,
        PrefixPolicy::Forbidden);
}

Unit
Unit::Gram()
{
    return Unit(
        "g",
        PhysicalDimension::Mass(),
        1.0e-3,
        PrefixPolicy::Allowed);
}

Unit
Unit::Second()
{
    return Unit(
        "s",
        PhysicalDimension::Time(),
        1.0,
        PrefixPolicy::Allowed);
}

Unit
Unit::Kelvin()
{
    return Unit(
        "K",
        PhysicalDimension::Temperature(),
        1.0,
        PrefixPolicy::Allowed);
}

Unit
Unit::Hertz()
{
    return Unit(
        "Hz",
        PhysicalDimension::Time().pow(
            Rational(-1)),
        1.0,
        PrefixPolicy::Allowed);
}

Unit
Unit::Newton()
{
    return Unit(
        "N",
        PhysicalDimension::Force(),
        1.0,
        PrefixPolicy::Allowed);
}

Unit
Unit::Pascal()
{
    return Unit(
        "Pa",
        PhysicalDimension::Pressure(),
        1.0,
        PrefixPolicy::Allowed);
}

Unit
Unit::Joule()
{
    return Unit(
        "J",
        PhysicalDimension::Energy(),
        1.0,
        PrefixPolicy::Allowed);
}

Unit
Unit::Watt()
{
    return Unit(
        "W",
        PhysicalDimension::Power(),
        1.0,
        PrefixPolicy::Allowed);
}

Unit
Unit::Bar()
{
    return Unit(
        "bar",
        PhysicalDimension::Pressure(),
        1.0e5,
        PrefixPolicy::Forbidden);
}

Unit
Unit::Atmosphere()
{
    return Unit(
        "atm",
        PhysicalDimension::Pressure(),
        101325.0,
        PrefixPolicy::Forbidden);
}

Unit
Unit::Radian()
{
    return Unit(
        "rad",
        PhysicalDimension::Dimensionless(),
        1.0,
        PrefixPolicy::Forbidden);
}

Unit
Unit::Degree()
{
    return Unit(
        "deg",
        PhysicalDimension::Dimensionless(),
        pi / 180.0,
        PrefixPolicy::Forbidden);
}

} // namespace ntic::lbm::config
