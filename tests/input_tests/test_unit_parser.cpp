#include "unit_parser.hpp"
#include "test_framework.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

using ntic::lbm::common::Rational;

using ntic::lbm::config::PhysicalDimension;
using ntic::lbm::config::PrefixRegistry;
using ntic::lbm::config::UnitExpression;
using ntic::lbm::config::UnitParser;
using ntic::lbm::config::UnitRegistry;

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

    return std::abs(lhs - rhs)
        <= relativeTolerance * scale;
}

void expectDimension(
    const UnitParser& parser,
    const std::string& expression,
    const PhysicalDimension& expected)
{
    const UnitExpression result =
        parser.parse(expression);

    NTIC_EXPECT_EQ(
        result.dimension(),
        expected);
}

void expectScale(
    const UnitParser& parser,
    const std::string& expression,
    double expected,
    double relativeTolerance = 1.0e-12)
{
    const UnitExpression result =
        parser.parse(expression);

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            result.scale(),
            expected,
            relativeTolerance));
}

void expectUnit(
    const UnitParser& parser,
    const std::string& expression,
    const PhysicalDimension& expectedDimension,
    double expectedScale,
    double relativeTolerance = 1.0e-12)
{
    const UnitExpression result =
        parser.parse(expression);

    NTIC_EXPECT_EQ(
        result.dimension(),
        expectedDimension);

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            result.scale(),
            expectedScale,
            relativeTolerance));
}

void expectEquivalent(
    const UnitParser& parser,
    const std::string& lhs,
    const std::string& rhs)
{
    const UnitExpression lhsUnit =
        parser.parse(lhs);

    const UnitExpression rhsUnit =
        parser.parse(rhs);

    NTIC_EXPECT_TRUE(
        lhsUnit.equivalentTo(rhsUnit));
}

} // namespace

int main()
{
    const UnitRegistry unitRegistry;
    const PrefixRegistry prefixRegistry;

    const UnitParser parser(
        unitRegistry,
        prefixRegistry);

    constexpr double pi =
        3.141592653589793238462643383279502884;

    //----------------------------------------------------------
    // Empty Expression
    //----------------------------------------------------------

    {
        const UnitExpression result =
            parser.parse("");

        NTIC_EXPECT_TRUE(
            result.isDimensionless());

        NTIC_EXPECT_TRUE(
            result.isStandardSI());

        NTIC_EXPECT_TRUE(
            nearlyEqual(
                result.scale(),
                1.0));
    }

    {
        const UnitExpression result =
            parser.parse(
                "   \t\r\n   # comment only");

        NTIC_EXPECT_TRUE(
            result.isDimensionless());

        NTIC_EXPECT_TRUE(
            nearlyEqual(
                result.scale(),
                1.0));
    }

    //----------------------------------------------------------
    // Built-in Base Units
    //----------------------------------------------------------

    expectUnit(
        parser,
        "m",
        PhysicalDimension::Length(),
        1.0);

    expectUnit(
        parser,
        "kg",
        PhysicalDimension::Mass(),
        1.0);

    expectUnit(
        parser,
        "g",
        PhysicalDimension::Mass(),
        1.0e-3);

    expectUnit(
        parser,
        "s",
        PhysicalDimension::Time(),
        1.0);

    expectUnit(
        parser,
        "K",
        PhysicalDimension::Temperature(),
        1.0);

    //----------------------------------------------------------
    // Derived SI Units
    //----------------------------------------------------------

    expectUnit(
        parser,
        "Hz",
        PhysicalDimension::Time().pow(
            Rational(-1)),
        1.0);

    expectUnit(
        parser,
        "N",
        PhysicalDimension::Force(),
        1.0);

    expectUnit(
        parser,
        "Pa",
        PhysicalDimension::Pressure(),
        1.0);

    expectUnit(
        parser,
        "J",
        PhysicalDimension::Energy(),
        1.0);

    expectUnit(
        parser,
        "W",
        PhysicalDimension::Power(),
        1.0);

    //----------------------------------------------------------
    // Engineering Pressure Units
    //----------------------------------------------------------

    expectUnit(
        parser,
        "bar",
        PhysicalDimension::Pressure(),
        1.0e5);

    expectUnit(
        parser,
        "atm",
        PhysicalDimension::Pressure(),
        101325.0);

    //----------------------------------------------------------
    // Angle Units Are Dimensionless
    //----------------------------------------------------------

    expectUnit(
        parser,
        "rad",
        PhysicalDimension::Dimensionless(),
        1.0);

    expectUnit(
        parser,
        "deg",
        PhysicalDimension::Dimensionless(),
        pi / 180.0);

    expectDimension(
        parser,
        "rad/s",
        PhysicalDimension::Time().pow(
            Rational(-1)));

    expectScale(
        parser,
        "deg/s",
        pi / 180.0);

    //----------------------------------------------------------
    // Complete Unit Match Has Priority
    //----------------------------------------------------------

    // "kg" must resolve directly as kilogram rather than
    // being interpreted as kilo + gram.
    expectUnit(
        parser,
        "kg",
        PhysicalDimension::Mass(),
        1.0);

    //----------------------------------------------------------
    // Length Prefixes
    //----------------------------------------------------------

    expectUnit(
        parser,
        "Qm",
        PhysicalDimension::Length(),
        1.0e30);

    expectUnit(
        parser,
        "Rm",
        PhysicalDimension::Length(),
        1.0e27);

    expectUnit(
        parser,
        "Ym",
        PhysicalDimension::Length(),
        1.0e24);

    expectUnit(
        parser,
        "Zm",
        PhysicalDimension::Length(),
        1.0e21);

    expectUnit(
        parser,
        "Em",
        PhysicalDimension::Length(),
        1.0e18);

    expectUnit(
        parser,
        "Pm",
        PhysicalDimension::Length(),
        1.0e15);

    expectUnit(
        parser,
        "Tm",
        PhysicalDimension::Length(),
        1.0e12);

    expectUnit(
        parser,
        "Gm",
        PhysicalDimension::Length(),
        1.0e9);

    expectUnit(
        parser,
        "Mm",
        PhysicalDimension::Length(),
        1.0e6);

    expectUnit(
        parser,
        "km",
        PhysicalDimension::Length(),
        1.0e3);

    expectUnit(
        parser,
        "hm",
        PhysicalDimension::Length(),
        1.0e2);

    expectUnit(
        parser,
        "dam",
        PhysicalDimension::Length(),
        1.0e1);

    expectUnit(
        parser,
        "dm",
        PhysicalDimension::Length(),
        1.0e-1);

    expectUnit(
        parser,
        "cm",
        PhysicalDimension::Length(),
        1.0e-2);

    expectUnit(
        parser,
        "mm",
        PhysicalDimension::Length(),
        1.0e-3);

    expectUnit(
        parser,
        "um",
        PhysicalDimension::Length(),
        1.0e-6);

    expectUnit(
        parser,
        "nm",
        PhysicalDimension::Length(),
        1.0e-9);

    expectUnit(
        parser,
        "pm",
        PhysicalDimension::Length(),
        1.0e-12);

    //----------------------------------------------------------
    // Mass Prefixes Are Applied to Gram
    //----------------------------------------------------------

    expectUnit(
        parser,
        "Mg",
        PhysicalDimension::Mass(),
        1.0e3);

    expectUnit(
        parser,
        "kg",
        PhysicalDimension::Mass(),
        1.0);

    expectUnit(
        parser,
        "g",
        PhysicalDimension::Mass(),
        1.0e-3);

    expectUnit(
        parser,
        "mg",
        PhysicalDimension::Mass(),
        1.0e-6);

    expectUnit(
        parser,
        "ug",
        PhysicalDimension::Mass(),
        1.0e-9);

    expectUnit(
        parser,
        "ng",
        PhysicalDimension::Mass(),
        1.0e-12);

    //----------------------------------------------------------
    // Time Prefixes
    //----------------------------------------------------------

    expectUnit(
        parser,
        "ks",
        PhysicalDimension::Time(),
        1.0e3);

    expectUnit(
        parser,
        "ms",
        PhysicalDimension::Time(),
        1.0e-3);

    expectUnit(
        parser,
        "us",
        PhysicalDimension::Time(),
        1.0e-6);

    expectUnit(
        parser,
        "ns",
        PhysicalDimension::Time(),
        1.0e-9);

    //----------------------------------------------------------
    // Derived Unit Prefixes
    //----------------------------------------------------------

    expectUnit(
        parser,
        "kHz",
        PhysicalDimension::Time().pow(
            Rational(-1)),
        1.0e3);

    expectUnit(
        parser,
        "MHz",
        PhysicalDimension::Time().pow(
            Rational(-1)),
        1.0e6);

    expectUnit(
        parser,
        "GHz",
        PhysicalDimension::Time().pow(
            Rational(-1)),
        1.0e9);

    expectUnit(
        parser,
        "mN",
        PhysicalDimension::Force(),
        1.0e-3);

    expectUnit(
        parser,
        "kN",
        PhysicalDimension::Force(),
        1.0e3);

    expectUnit(
        parser,
        "MN",
        PhysicalDimension::Force(),
        1.0e6);

    expectUnit(
        parser,
        "mPa",
        PhysicalDimension::Pressure(),
        1.0e-3);

    expectUnit(
        parser,
        "kPa",
        PhysicalDimension::Pressure(),
        1.0e3);

    expectUnit(
        parser,
        "MPa",
        PhysicalDimension::Pressure(),
        1.0e6);

    expectUnit(
        parser,
        "GPa",
        PhysicalDimension::Pressure(),
        1.0e9);

    expectUnit(
        parser,
        "kJ",
        PhysicalDimension::Energy(),
        1.0e3);

    expectUnit(
        parser,
        "MJ",
        PhysicalDimension::Energy(),
        1.0e6);

    expectUnit(
        parser,
        "mW",
        PhysicalDimension::Power(),
        1.0e-3);

    expectUnit(
        parser,
        "kW",
        PhysicalDimension::Power(),
        1.0e3);

    expectUnit(
        parser,
        "MW",
        PhysicalDimension::Power(),
        1.0e6);

    //----------------------------------------------------------
    // Explicit Multiplication
    //----------------------------------------------------------

    expectDimension(
        parser,
        "kg*m",
        PhysicalDimension::Mass()
        * PhysicalDimension::Length());

    expectDimension(
        parser,
        "N*m",
        PhysicalDimension::Energy());

    expectDimension(
        parser,
        "Pa*s",
        PhysicalDimension::DynamicViscosity());

    expectDimension(
        parser,
        "kg*m*s",
        PhysicalDimension::Mass()
        * PhysicalDimension::Length()
        * PhysicalDimension::Time());

    expectScale(
        parser,
        "g*cm",
        1.0e-5);

    expectScale(
        parser,
        "kN*cm",
        10.0);

    //----------------------------------------------------------
    // Implicit Multiplication
    //----------------------------------------------------------

    expectEquivalent(
        parser,
        "kg m",
        "kg*m");

    expectEquivalent(
        parser,
        "N m",
        "N*m");

    expectEquivalent(
        parser,
        "Pa s",
        "Pa*s");

    expectEquivalent(
        parser,
        "kg(m)",
        "kg*m");

    expectEquivalent(
        parser,
        "(kg)m",
        "kg*m");

    expectEquivalent(
        parser,
        "(kg)(m)",
        "kg*m");

    expectEquivalent(
        parser,
        "kg[m]",
        "kg*m");

    expectEquivalent(
        parser,
        "[kg]m",
        "kg*m");

    expectEquivalent(
        parser,
        "{kg}(m)",
        "kg*m");

    //----------------------------------------------------------
    // Division
    //----------------------------------------------------------

    expectDimension(
        parser,
        "m/s",
        PhysicalDimension::Velocity());

    expectDimension(
        parser,
        "m/s^2",
        PhysicalDimension::Acceleration());

    expectDimension(
        parser,
        "kg/m^3",
        PhysicalDimension::Density());

    expectDimension(
        parser,
        "Pa/s",
        PhysicalDimension::Pressure()
        / PhysicalDimension::Time());

    expectDimension(
        parser,
        "N/m",
        PhysicalDimension::Force()
        / PhysicalDimension::Length());

    expectDimension(
        parser,
        "kg/(m*s)",
        PhysicalDimension::DynamicViscosity());

    expectScale(
        parser,
        "g/cm^3",
        1000.0);

    expectScale(
        parser,
        "mg/mm^3",
        1000.0);

    //----------------------------------------------------------
    // Left Associativity
    //----------------------------------------------------------

    expectEquivalent(
        parser,
        "kg/m/s",
        "(kg/m)/s");

    expectEquivalent(
        parser,
        "m/s/s",
        "(m/s)/s");

    expectEquivalent(
        parser,
        "N/m/m",
        "(N/m)/m");

    //----------------------------------------------------------
    // Integer Exponents
    //----------------------------------------------------------

    expectDimension(
        parser,
        "m^0",
        PhysicalDimension::Dimensionless());

    expectScale(
        parser,
        "cm^0",
        1.0);

    expectDimension(
        parser,
        "m^1",
        PhysicalDimension::Length());

    expectDimension(
        parser,
        "m^2",
        PhysicalDimension::Area());

    expectDimension(
        parser,
        "m^3",
        PhysicalDimension::Volume());

    expectDimension(
        parser,
        "m^-1",
        PhysicalDimension::Length().pow(
            Rational(-1)));

    expectDimension(
        parser,
        "m^-2",
        PhysicalDimension::Length().pow(
            Rational(-2)));

    expectScale(
        parser,
        "cm^2",
        1.0e-4);

    expectScale(
        parser,
        "cm^3",
        1.0e-6);

    expectScale(
        parser,
        "cm^-1",
        1.0e2);

    //----------------------------------------------------------
    // Powers of Bracketed Expressions
    //----------------------------------------------------------

    expectDimension(
        parser,
        "(m/s)^2",
        PhysicalDimension::Velocity().pow(
            Rational(2)));

    expectDimension(
        parser,
        "[m/s]^2",
        PhysicalDimension::Velocity().pow(
            Rational(2)));

    expectDimension(
        parser,
        "{m/s}^2",
        PhysicalDimension::Velocity().pow(
            Rational(2)));

    expectScale(
        parser,
        "(cm/s)^2",
        1.0e-4);

    //----------------------------------------------------------
    // Rational Exponents
    //----------------------------------------------------------

    expectDimension(
        parser,
        "m^(1/2)",
        PhysicalDimension::Length().pow(
            Rational(1, 2)));

    expectDimension(
        parser,
        "m^[1/2]",
        PhysicalDimension::Length().pow(
            Rational(1, 2)));

    expectDimension(
        parser,
        "m^{1/2}",
        PhysicalDimension::Length().pow(
            Rational(1, 2)));

    expectDimension(
        parser,
        "m^(-1/2)",
        PhysicalDimension::Length().pow(
            Rational(-1, 2)));

    expectDimension(
        parser,
        "m^[3/2]",
        PhysicalDimension::Length().pow(
            Rational(3, 2)));

    expectDimension(
        parser,
        "m^{-3/2}",
        PhysicalDimension::Length().pow(
            Rational(-3, 2)));

    expectScale(
        parser,
        "cm^(1/2)",
        1.0e-1);

    expectScale(
        parser,
        "cm^(-1/2)",
        1.0e1);

    //----------------------------------------------------------
    // Same-shaped Bracket Nesting
    //----------------------------------------------------------

    expectEquivalent(
        parser,
        "(m)",
        "m");

    expectEquivalent(
        parser,
        "((m))",
        "m");

    expectEquivalent(
        parser,
        "(((m)))",
        "m");

    expectEquivalent(
        parser,
        "[m]",
        "m");

    expectEquivalent(
        parser,
        "[[m]]",
        "m");

    expectEquivalent(
        parser,
        "[[[m]]]",
        "m");

    expectEquivalent(
        parser,
        "{m}",
        "m");

    expectEquivalent(
        parser,
        "{{m}}",
        "m");

    expectEquivalent(
        parser,
        "{{{m}}}",
        "m");

    //----------------------------------------------------------
    // Valid Mixed Bracket Nesting
    //----------------------------------------------------------

    expectEquivalent(
        parser,
        "[((m))]",
        "m");

    expectEquivalent(
        parser,
        "{[m]}",
        "m");

    expectEquivalent(
        parser,
        "{(m)}",
        "m");

    expectEquivalent(
        parser,
        "{[(m)]}",
        "m");

    expectEquivalent(
        parser,
        "{{[[((m))]]}}",
        "m");

    expectEquivalent(
        parser,
        "{[kg/(m*s)]}",
        "kg/(m*s)");

    //----------------------------------------------------------
    // Equivalent Unit Definitions
    //----------------------------------------------------------

    expectEquivalent(
        parser,
        "N",
        "kg*m/s^2");

    expectEquivalent(
        parser,
        "Pa",
        "N/m^2");

    expectEquivalent(
        parser,
        "Pa",
        "kg/(m*s^2)");

    expectEquivalent(
        parser,
        "J",
        "N*m");

    expectEquivalent(
        parser,
        "J",
        "kg*m^2/s^2");

    expectEquivalent(
        parser,
        "W",
        "J/s");

    expectEquivalent(
        parser,
        "W",
        "kg*m^2/s^3");

    expectEquivalent(
        parser,
        "Hz",
        "s^-1");

    expectEquivalent(
        parser,
        "Pa*s",
        "kg/(m*s)");

    expectEquivalent(
        parser,
        "m^2/s",
        "m*m/s");

    //----------------------------------------------------------
    // Equivalent Scaled Units
    //----------------------------------------------------------

    expectEquivalent(
        parser,
        "MPa",
        "N/mm^2");

    expectEquivalent(
        parser,
        "kPa",
        "kN/m^2");

    expectEquivalent(
        parser,
        "J",
        "kN*mm");

    expectEquivalent(
        parser,
        "W",
        "kJ/ks");

    expectEquivalent(
        parser,
        "kg/m^3",
        "mg/cm^3");

    //----------------------------------------------------------
    // Whitespace and Comments
    //----------------------------------------------------------

    expectEquivalent(
        parser,
        " kg / ( m * s ) ",
        "kg/(m*s)");

    expectEquivalent(
        parser,
        "kg\tm\n/\ts^2",
        "kg*m/s^2");

    expectEquivalent(
        parser,
        "Pa # pressure unit",
        "Pa");

    expectEquivalent(
        parser,
        "kg/(m*s) # dynamic viscosity",
        "Pa*s");

    //----------------------------------------------------------
    // Unknown Units
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse("abc"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("xyz"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("foo^2"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("kg/abc"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("unknown_unit"),
        std::runtime_error);

    //----------------------------------------------------------
    // Forbidden Prefixes
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse("mkg"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("kkg"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("mbar"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("katm"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("mrad"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("mdeg"),
        std::runtime_error);

    //----------------------------------------------------------
    // Invalid Operators
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse("+m"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("-m"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m+"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m-"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m+s"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m-s"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m//s"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m/*s"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m**s"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("*m"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("/m"),
        std::runtime_error);

    //----------------------------------------------------------
    // Invalid Exponents
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse("m^"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m^^2"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m^s"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m^1.5"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m^1e2"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m^(1.5)"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m^(1/0)"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m^(1/)"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m^(/2)"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m^()"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m^[1/2)"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m^{1/2]"),
        std::runtime_error);

    //----------------------------------------------------------
    // Empty Brackets
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse("()"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("[]"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("{}"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("[()]"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("{[]}"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("{[()]}"),
        std::runtime_error);

    //----------------------------------------------------------
    // Missing and Mismatched Brackets
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse("(m"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("[m"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("{m"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m)"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m]"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m}"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("(m]"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("[m)"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("{m]"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("{[m)}"),
        std::runtime_error);

    //----------------------------------------------------------
    // Invalid Mixed Bracket Nesting
    //----------------------------------------------------------

    // A smaller bracket may not contain a larger bracket.
    NTIC_EXPECT_THROW(
        parser.parse("([m])"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("({m})"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("([{m}])"),
        std::runtime_error);

    // A square bracket may not contain braces.
    NTIC_EXPECT_THROW(
        parser.parse("[{m}]"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("[({m})]"),
        std::runtime_error);

    // Same-shaped outer brackets do not relax the rule.
    NTIC_EXPECT_THROW(
        parser.parse("((([m])))"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("[[[{m}]]]"),
        std::runtime_error);

    //----------------------------------------------------------
    // Numeric Coefficients Are Not Unit Expressions
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse("2*m"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("100000*Pa"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("1e5 Pa"),
        std::runtime_error);

    //----------------------------------------------------------
    // Strings Are Not Unit Expressions
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse("\"m\""),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("'Pa'"),
        std::runtime_error);

    //----------------------------------------------------------
    // Unexpected Configuration Symbols
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse("unit=m"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m,"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("m:"),
        std::runtime_error);

    //----------------------------------------------------------

    ntic::lbm::test::printSummary();

    return 0;
}

