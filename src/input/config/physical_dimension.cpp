#include "physical_dimension.hpp"

#include <sstream>

namespace ntic::lbm::config
{

//==============================================================
// Constructors
//==============================================================

PhysicalDimension::PhysicalDimension()
{
    exponent_.fill(common::Rational(0));
}

PhysicalDimension::PhysicalDimension(
    int L,
    int M,
    int T,
    int Theta)
{
    exponent_[0] = common::Rational(L);
    exponent_[1] = common::Rational(M);
    exponent_[2] = common::Rational(T);
    exponent_[3] = common::Rational(Theta);
}

PhysicalDimension::PhysicalDimension(
    common::Rational L,
    common::Rational M,
    common::Rational T,
    common::Rational Theta)
{
    exponent_[0] = L;
    exponent_[1] = M;
    exponent_[2] = T;
    exponent_[3] = Theta;
}

//==============================================================
// Access
//==============================================================

const common::Rational&
PhysicalDimension::operator[](BaseDimension dim) const
{
    return exponent_[static_cast<size_t>(dim)];
}

common::Rational&
PhysicalDimension::operator[](BaseDimension dim)
{
    return exponent_[static_cast<size_t>(dim)];
}

//==============================================================
// Comparison
//==============================================================

bool PhysicalDimension::operator==(const PhysicalDimension& rhs) const
{
    return exponent_ == rhs.exponent_;
}

bool PhysicalDimension::operator!=(const PhysicalDimension& rhs) const
{
    return !(*this == rhs);
}

//==============================================================
// Algebra
//==============================================================

PhysicalDimension
PhysicalDimension::operator*(
    const PhysicalDimension& rhs) const
{
    PhysicalDimension result;

    for (size_t i = 0; i < exponent_.size(); ++i)
    {
        result.exponent_[i] =
            exponent_[i] + rhs.exponent_[i];
    }

    return result;
}

PhysicalDimension
PhysicalDimension::operator/(
    const PhysicalDimension& rhs) const
{
    PhysicalDimension result;

    for (size_t i = 0; i < exponent_.size(); ++i)
    {
        result.exponent_[i] =
            exponent_[i] - rhs.exponent_[i];
    }

    return result;
}

PhysicalDimension
PhysicalDimension::pow(
    const common::Rational& exponent) const
{
    PhysicalDimension result;

    for (size_t i = 0; i < exponent_.size(); ++i)
    {
        result.exponent_[i] =
            exponent_[i] * exponent;
    }

    return result;
}

//==============================================================
// Utility
//==============================================================

bool PhysicalDimension::isDimensionless() const
{
    for (const auto& e : exponent_)
    {
        if (!e.isZero())
            return false;
    }

    return true;
}

bool PhysicalDimension::isIntegral() const
{
    for (const auto& e : exponent_)
    {
        if (!e.isInteger())
            return false;
    }

    return true;
}

std::string
PhysicalDimension::toString() const
{
    static const char* symbol[] =
    {
        "L",
        "M",
        "T",
        "Θ"
    };

    if (isDimensionless())
        return "[1]";

    std::ostringstream oss;

    oss << "[";

    bool first = true;

    for (size_t i = 0; i < exponent_.size(); ++i)
    {
        if (exponent_[i].isZero())
            continue;

        if (!first)
            oss << " ";

        oss << symbol[i];

        if (exponent_[i] != common::Rational(1))
        {
            if (exponent_[i].isInteger())
            {
                oss << "^"
                    << exponent_[i].numerator();
            }
            else
            {
                oss << "^("
                    << exponent_[i].toString()
                    << ")";
            }
        }

        first = false;
    }

    oss << "]";

    return oss.str();
}

//==============================================================
// Base Dimensions
//==============================================================

PhysicalDimension
PhysicalDimension::Dimensionless()
{
    return PhysicalDimension();
}

PhysicalDimension
PhysicalDimension::Length()
{
    return PhysicalDimension(1,0,0,0);
}

PhysicalDimension
PhysicalDimension::Mass()
{
    return PhysicalDimension(0,1,0,0);
}

PhysicalDimension
PhysicalDimension::Time()
{
    return PhysicalDimension(0,0,1,0);
}

PhysicalDimension
PhysicalDimension::Temperature()
{
    return PhysicalDimension(0,0,0,1);
}

//==============================================================
// Derived Dimensions
//==============================================================

PhysicalDimension
PhysicalDimension::Area()
{
    return Length().pow(common::Rational(2));
}

PhysicalDimension
PhysicalDimension::Volume()
{
    return Length().pow(common::Rational(3));
}

PhysicalDimension
PhysicalDimension::Velocity()
{
    return Length() / Time();
}

PhysicalDimension
PhysicalDimension::Acceleration()
{
    return Velocity() / Time();
}

PhysicalDimension
PhysicalDimension::Force()
{
    return Mass() * Acceleration();
}

PhysicalDimension
PhysicalDimension::Pressure()
{
    return Force() / Area();
}

PhysicalDimension
PhysicalDimension::Density()
{
    return Mass() / Volume();
}

PhysicalDimension
PhysicalDimension::DynamicViscosity()
{
    return Pressure() * Time();
}

PhysicalDimension
PhysicalDimension::KinematicViscosity()
{
    return Area() / Time();
}

PhysicalDimension
PhysicalDimension::Energy()
{
    return Force() * Length();
}

PhysicalDimension
PhysicalDimension::Power()
{
    return Energy() / Time();
}

} // namespace ntic::lbm::config
