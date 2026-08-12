#pragma once

#include <array>
#include <string>

#include "rational.hpp"

namespace ntic::lbm::config
{

class PhysicalDimension
{
public:

    enum class BaseDimension : uint8_t
    {
        Length = 0,
        Mass,
        Time,
        Temperature,

        NumDimensions
    };

public:

    //-----------------------------------------
    // Constructors
    //-----------------------------------------

    PhysicalDimension();

    PhysicalDimension(
        int L,
        int M,
        int T,
        int Theta);

    PhysicalDimension(
        common::Rational L,
        common::Rational M,
        common::Rational T,
        common::Rational Theta);

    //-----------------------------------------
    // Access
    //-----------------------------------------

    const common::Rational&
    operator[](BaseDimension dim) const;

    common::Rational&
    operator[](BaseDimension dim);

    //-----------------------------------------
    // Comparison
    //-----------------------------------------

    bool operator==(const PhysicalDimension&) const;

    bool operator!=(const PhysicalDimension&) const;

    //-----------------------------------------
    // Algebra
    //-----------------------------------------

    PhysicalDimension operator*(
        const PhysicalDimension&) const;

    PhysicalDimension operator/(
        const PhysicalDimension&) const;

    PhysicalDimension pow(
        const common::Rational&) const;

    //-----------------------------------------
    // Utility
    //-----------------------------------------

    bool isDimensionless() const;

    bool isIntegral() const;

    std::string toString() const;

    //-----------------------------------------
    // Factory
    //-----------------------------------------

    static PhysicalDimension Dimensionless();

    static PhysicalDimension Length();

    static PhysicalDimension Mass();

    static PhysicalDimension Time();

    static PhysicalDimension Temperature();

    static PhysicalDimension Area();

    static PhysicalDimension Volume();

    static PhysicalDimension Velocity();

    static PhysicalDimension Acceleration();

    static PhysicalDimension Force();

    static PhysicalDimension Pressure();

    static PhysicalDimension Density();

    static PhysicalDimension DynamicViscosity();

    static PhysicalDimension KinematicViscosity();

    static PhysicalDimension Energy();

    static PhysicalDimension Power();

private:

    std::array<common::Rational,4> exponent_;

};

}
