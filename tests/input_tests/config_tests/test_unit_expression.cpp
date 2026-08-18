#include "unit_expression.hpp"
#include "unit.hpp"
#include "physical_dimension.hpp"
#include "test_framework.hpp"

#include <algorithm>
#include <cmath>
#include <string>

using ntic::lbm::common::Rational;

using ntic::lbm::config::PhysicalDimension;
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
        {
            1.0,
            std::abs(lhs),
            std::abs(rhs)
        });

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
    UnitExpression expr;

    NTIC_EXPECT_TRUE(
        expr.isDimensionless());

    NTIC_EXPECT_TRUE(
        expr.isStandardSI());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            expr.scale(),
            1.0));
}

//----------------------------------------------------------
// Construct From Unit
//----------------------------------------------------------

{
    UnitExpression expr(
        Unit::Meter());

    NTIC_EXPECT_EQ(
        expr.dimension(),
        PhysicalDimension::Length());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            expr.scale(),
            1.0));
}

{
    UnitExpression expr(
        Unit::Gram());

    NTIC_EXPECT_EQ(
        expr.dimension(),
        PhysicalDimension::Mass());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            expr.scale(),
            1.0e-3));
}

//----------------------------------------------------------
// Construct From Dimension
//----------------------------------------------------------

{
    UnitExpression expr(
        PhysicalDimension::Velocity(),
        123.0);

    NTIC_EXPECT_EQ(
        expr.dimension(),
        PhysicalDimension::Velocity());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            expr.scale(),
            123.0));
}

//----------------------------------------------------------
// Invalid Constructor
//----------------------------------------------------------

NTIC_EXPECT_THROW(

    UnitExpression(
        PhysicalDimension::Length(),
        0.0),

    std::invalid_argument);

NTIC_EXPECT_THROW(

    UnitExpression(
        PhysicalDimension::Length(),
        -1.0),

    std::invalid_argument);

//----------------------------------------------------------
// Multiplication
//----------------------------------------------------------

{
    UnitExpression length(
        Unit::Meter());

    UnitExpression time(
        Unit::Second());

    UnitExpression result =
        length*time;

    NTIC_EXPECT_EQ(
        result.dimension(),
        PhysicalDimension::Length()
        *
        PhysicalDimension::Time());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            result.scale(),
            1.0));
}

//----------------------------------------------------------
// Division
//----------------------------------------------------------

{
    UnitExpression length(
        Unit::Meter());

    UnitExpression time(
        Unit::Second());

    UnitExpression velocity =
        length/time;

    NTIC_EXPECT_EQ(
        velocity.dimension(),
        PhysicalDimension::Velocity());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            velocity.scale(),
            1.0));
}

//----------------------------------------------------------
// Square
//----------------------------------------------------------

{
    UnitExpression length(
        Unit::Meter());

    UnitExpression area =
        length.pow(
            Rational(2));

    NTIC_EXPECT_EQ(
        area.dimension(),
        PhysicalDimension::Area());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            area.scale(),
            1.0));
}

//----------------------------------------------------------
// Cube
//----------------------------------------------------------

{
    UnitExpression length(
        Unit::Meter());

    UnitExpression volume =
        length.pow(
            Rational(3));

    NTIC_EXPECT_EQ(
        volume.dimension(),
        PhysicalDimension::Volume());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            volume.scale(),
            1.0));
}

//----------------------------------------------------------
// Inverse
//----------------------------------------------------------

{
    UnitExpression length(
        Unit::Meter());

    UnitExpression inv =
        length.pow(
            Rational(-1));

    NTIC_EXPECT_EQ(
        inv.dimension(),
        PhysicalDimension::Length().pow(
            Rational(-1)));

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            inv.scale(),
            1.0));
}

//----------------------------------------------------------
// Square Root
//----------------------------------------------------------

{
    UnitExpression area(

        PhysicalDimension::Area(),

        1.0e-4);

    UnitExpression length =
        area.pow(
            Rational(1,2));

    NTIC_EXPECT_EQ(
        length.dimension(),
        PhysicalDimension::Length());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            length.scale(),
            1.0e-2));
}

//----------------------------------------------------------
// Density
//----------------------------------------------------------

{
    UnitExpression kg(
        Unit::Kilogram());

    UnitExpression meter(
        Unit::Meter());

    UnitExpression density =
        kg
        /
        meter.pow(
            Rational(3));

    NTIC_EXPECT_EQ(
        density.dimension(),
        PhysicalDimension::Density());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            density.scale(),
            1.0));
}

//----------------------------------------------------------
// g/cm^3
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

    NTIC_EXPECT_EQ(
        density.dimension(),
        PhysicalDimension::Density());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            density.scale(),
            1000.0));
}

//----------------------------------------------------------
// Dimensionless
//----------------------------------------------------------

{
    UnitExpression length(
        Unit::Meter());

    UnitExpression result =
        length/length;

    NTIC_EXPECT_TRUE(
        result.isDimensionless());

    NTIC_EXPECT_TRUE(
        result.isStandardSI());

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            result.scale(),
            1.0));
}

//----------------------------------------------------------
// Equality
//----------------------------------------------------------

{
    UnitExpression a(
        Unit::Meter());

    UnitExpression b(
        Unit::Meter());

    UnitExpression c(
        Unit::Second());

    NTIC_EXPECT_TRUE(
        a==b);

    NTIC_EXPECT_TRUE(
        a!=c);
}

//----------------------------------------------------------
// equivalentTo()
//----------------------------------------------------------

{
    UnitExpression a(

        PhysicalDimension::Length(),

        1.0);

    UnitExpression b(

        PhysicalDimension::Length(),

        1.0+1.0e-14);

    NTIC_EXPECT_TRUE(
        a.equivalentTo(b));
}

//----------------------------------------------------------
// Different Dimension
//----------------------------------------------------------

{
    UnitExpression a(
        Unit::Meter());

    UnitExpression b(
        Unit::Second());

    NTIC_EXPECT_TRUE(
        !a.equivalentTo(b));
}

//----------------------------------------------------------
// toString()
//----------------------------------------------------------

{
    std::string text =
        UnitExpression(
            Unit::Meter())
        .toString();

    NTIC_EXPECT_TRUE(
        text.find("[L]")
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
