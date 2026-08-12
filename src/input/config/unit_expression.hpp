#pragma once

#include <string>

#include "rational.hpp"
#include "unit.hpp"

namespace ntic::lbm::config
{

class UnitExpression
{
public:
    //----------------------------------------------------------
    // Constructors
    //----------------------------------------------------------

    // Dimensionless expression with scale 1.
    UnitExpression();

    // Constructs an expression from one registered unit.
    explicit UnitExpression(const Unit& unit);

    // Constructs an already evaluated unit expression.
    UnitExpression(
        const PhysicalDimension& dimension,
        double scale);

    //----------------------------------------------------------
    // Access
    //----------------------------------------------------------

    const PhysicalDimension& dimension() const;

    // Multiplicative conversion factor to standard SI units.
    //
    // Examples:
    //
    //   m        -> 1
    //   cm       -> 1e-2
    //   g        -> 1e-3
    //   g/cm^3   -> 1e3
    double scale() const;

    //----------------------------------------------------------
    // Algebra
    //----------------------------------------------------------

    UnitExpression operator*(
        const UnitExpression& rhs) const;

    UnitExpression operator/(
        const UnitExpression& rhs) const;

    UnitExpression pow(
        const common::Rational& exponent) const;

    //----------------------------------------------------------
    // Comparison
    //----------------------------------------------------------

    bool operator==(
        const UnitExpression& rhs) const;

    bool operator!=(
        const UnitExpression& rhs) const;

    // Floating-point tolerant comparison of scale factors.
    bool equivalentTo(
        const UnitExpression& rhs,
        double relativeTolerance = 1.0e-12) const;

    //----------------------------------------------------------
    // Utility
    //----------------------------------------------------------

    bool isDimensionless() const;

    bool isStandardSI() const;

    std::string toString() const;

private:
    PhysicalDimension dimension_;

    double scale_;
};

} // namespace ntic::lbm::config
