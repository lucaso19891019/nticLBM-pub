#include "quantity.hpp"
#include "unit.hpp"
#include "unit_expression.hpp"
#include "physical_dimension.hpp"
#include "test_framework.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

using ntic::lbm::common::Rational;

using ntic::lbm::config::PhysicalDimension;
using ntic::lbm::config::Quantity;
using ntic::lbm::config::Unit;
using ntic::lbm::config::UnitExpression;

namespace
{

bool nearlyEqual(
    double lhs,
    double rhs,
    double relativeTolerance = 1.0e-12)
{
    const double scale =
        std::max(
            1.0,
            std::max(
                std::abs(lhs),
                std::abs(rhs)));

    return std::abs(lhs-rhs)
        <= relativeTolerance*scale;
}

}

int main()
{

//----------------------------------------------------------
// Default Constructor
//----------------------------------------------------------

{
    Quantity q;

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            q.value(),
            0.0));

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            q.valueSI(),
            0.0));

    NTIC_EXPECT_TRUE(
        q.isDimensionless());

    NTIC_EXPECT_TRUE(
        q.isStandardSI());
}

//----------------------------------------------------------
// Dimensionless Helper
//----------------------------------------------------------

{
    Quantity q =
        Quantity::Dimensionless(
            123.0);

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            q.value(),
            123.0));

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            q.valueSI(),
            123.0));

    NTIC_EXPECT_TRUE(
        q.isDimensionless());

    NTIC_EXPECT_TRUE(
        q.isStandardSI());
}

//----------------------------------------------------------
// Constructor
//----------------------------------------------------------

{
    Quantity q(

        2.5,

        UnitExpression(
            Unit::Meter()));

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            q.value(),
            2.5));

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            q.valueSI(),
            2.5));

    NTIC_EXPECT_EQ(
        q.unit().dimension(),
        PhysicalDimension::Length());
}

//----------------------------------------------------------
// Invalid Constructor
//----------------------------------------------------------

NTIC_EXPECT_THROW(

    Quantity(

        std::numeric_limits<double>::infinity(),

        UnitExpression()),

    std::invalid_argument);

NTIC_EXPECT_THROW(

    Quantity(

        std::nan(""),

        UnitExpression()),

    std::invalid_argument);

//----------------------------------------------------------
// SI Conversion
//----------------------------------------------------------

{
    Quantity q(

        5.0,

        UnitExpression(
            Unit::Gram()));

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            q.value(),
            5.0));

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            q.valueSI(),
            5.0e-3));
}

//----------------------------------------------------------
// g/cm^3 -> kg/m^3
//----------------------------------------------------------

{
    UnitExpression gram(
        Unit::Gram());

    UnitExpression centimeter(

        PhysicalDimension::Length(),

        1.0e-2);

    UnitExpression density =
        gram
        /
        centimeter.pow(
            Rational(3));

    Quantity rho(

        2.7,

        density);

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            rho.valueSI(),
            2700.0));

    NTIC_EXPECT_EQ(
        rho.unit().dimension(),
        PhysicalDimension::Density());
}

//----------------------------------------------------------
// Multiplication
//----------------------------------------------------------

{
    Quantity a(

        2.0,

        UnitExpression(
            Unit::Meter()));

    Quantity b(

        5.0,

        UnitExpression(
            Unit::Meter()));

    Quantity area =
        a*b;

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            area.value(),
            10.0));

    NTIC_EXPECT_EQ(
        area.unit().dimension(),
        PhysicalDimension::Area());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            area.valueSI(),
            10.0));
}

//----------------------------------------------------------
// Division
//----------------------------------------------------------

{
    Quantity distance(

        20.0,

        UnitExpression(
            Unit::Meter()));

    Quantity time(

        5.0,

        UnitExpression(
            Unit::Second()));

    Quantity velocity =
        distance/time;

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            velocity.value(),
            4.0));

    NTIC_EXPECT_EQ(
        velocity.unit().dimension(),
        PhysicalDimension::Velocity());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            velocity.valueSI(),
            4.0));
}

//----------------------------------------------------------
// Division By Zero
//----------------------------------------------------------

NTIC_EXPECT_THROW(

    Quantity(

        1.0,

        UnitExpression())

    /

    Quantity(

        0.0,

        UnitExpression()),

    std::domain_error);

//----------------------------------------------------------
// Equality
//----------------------------------------------------------

{
    Quantity a(

        5.0,

        UnitExpression(
            Unit::Meter()));

    Quantity b(

        5.0,

        UnitExpression(
            Unit::Meter()));

    Quantity c(

        5.0,

        UnitExpression(
            Unit::Second()));

    NTIC_EXPECT_TRUE(
        a==b);

    NTIC_EXPECT_TRUE(
        a!=c);
}

//----------------------------------------------------------
// equivalentTo()
//----------------------------------------------------------

{
    Quantity a(

        1000.0,

        UnitExpression(
            Unit::Gram()));

    Quantity b(

        1.0,

        UnitExpression(
            Unit::Kilogram()));

    NTIC_EXPECT_TRUE(
        a.equivalentTo(
            b));
}

//----------------------------------------------------------
// Different Dimension
//----------------------------------------------------------

{
    Quantity a(

        1.0,

        UnitExpression(
            Unit::Meter()));

    Quantity b(

        1.0,

        UnitExpression(
            Unit::Second()));

    NTIC_EXPECT_TRUE(
        !a.equivalentTo(
            b));
}

//----------------------------------------------------------
// scale()
//----------------------------------------------------------

{
    Quantity q(

        1.0,

        UnitExpression(
            Unit::Gram()));

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            q.scale(),
            1.0e-3));
}

//----------------------------------------------------------
// isStandardSI()
//----------------------------------------------------------

{
    Quantity a(

        1.0,

        UnitExpression(
            Unit::Kilogram()));

    NTIC_EXPECT_TRUE(
        a.isStandardSI());

    Quantity b(

        1.0,

        UnitExpression(
            Unit::Gram()));

    NTIC_EXPECT_TRUE(
        !b.isStandardSI());
}

//----------------------------------------------------------
// toString()
//----------------------------------------------------------

{
    std::string text =

        Quantity(

            10.0,

            UnitExpression(
                Unit::Meter()))

        .toString();

    NTIC_EXPECT_TRUE(
        text.find("10")
        !=
        std::string::npos);

    NTIC_EXPECT_TRUE(
        text.find("scale=")
        !=
        std::string::npos);
}

//----------------------------------------------------------

ntic::lbm::test::printSummary();

return 0;

}
