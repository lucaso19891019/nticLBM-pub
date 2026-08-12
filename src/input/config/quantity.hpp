#pragma once

#include <string>

#include "unit_expression.hpp"

namespace ntic::lbm::config
{

class Quantity
{
public:

    //----------------------------------------------------------
    // Constructors
    //----------------------------------------------------------

    // Default:
    //      value = 0
    //      unit  = dimensionless
    Quantity();

    Quantity(
        double value,
        const UnitExpression& unit);

    // Helper for dimensionless quantities.
    static Quantity
    Dimensionless(double value);

public:

    //----------------------------------------------------------
    // Access
    //----------------------------------------------------------

    double value() const;

    double valueSI() const;

    double scale() const;

    const UnitExpression&
    unit() const;

public:

    //----------------------------------------------------------
    // Algebra
    //----------------------------------------------------------

    Quantity operator*(
        const Quantity& rhs) const;

    Quantity operator/(
        const Quantity& rhs) const;

public:

    //----------------------------------------------------------
    // Comparison
    //----------------------------------------------------------

    bool operator==(const Quantity& rhs) const;

    bool operator!=(const Quantity& rhs) const;

    bool equivalentTo(
        const Quantity& rhs,
        double relativeTolerance = 1.0e-12) const;

public:

    //----------------------------------------------------------
    // Utility
    //----------------------------------------------------------

    bool isDimensionless() const;

    bool isStandardSI() const;

    std::string toString() const;

private:

    double value_;

    UnitExpression unit_;
};

}
