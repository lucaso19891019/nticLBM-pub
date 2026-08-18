#include "quantity_parser.hpp"
#include "test_framework.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

using ntic::lbm::config::PrefixRegistry;
using ntic::lbm::config::Quantity;
using ntic::lbm::config::QuantityParser;
using ntic::lbm::config::UnitRegistry;

namespace
{

bool nearlyEqual(double a,double b,double eps=1e-12)
{
    return std::abs(a-b)<=eps*std::max(1.0,std::max(std::abs(a),std::abs(b)));
}

void expectValue(
    const QuantityParser& parser,
    const std::string& text,
    double value)
{
    Quantity q=parser.parse(text);

    NTIC_EXPECT_TRUE(
        nearlyEqual(
            q.value(),
            value));
}

void expectEquivalent(
    const QuantityParser& parser,
    const std::string& a,
    const std::string& b)
{
    NTIC_EXPECT_TRUE(
        parser.parse(a).equivalentTo(
            parser.parse(b)));
}

}

int main()
{
    const UnitRegistry unitRegistry;
    const PrefixRegistry prefixRegistry;
    const QuantityParser parser(
        unitRegistry,
        prefixRegistry);

    //----------------------------------------------------------
    // Dimensionless
    //----------------------------------------------------------

    expectValue(parser,"0",0.0);
    expectValue(parser,"1",1.0);
    expectValue(parser,"-1",-1.0);
    expectValue(parser,"+1",1.0);
    expectValue(parser,"1.25",1.25);
    expectValue(parser,"1e3",1000.0);
    expectValue(parser,"1e-3",1.0e-3);

    NTIC_EXPECT_TRUE(
        parser.parse("1").isDimensionless());

    //----------------------------------------------------------
    // Base Units
    //----------------------------------------------------------

    expectEquivalent(parser,"1 m","100 cm");
    expectEquivalent(parser,"1 kg","1000 g");
    expectEquivalent(parser,"1 s","1000 ms");

    //----------------------------------------------------------
    // Derived Units
    //----------------------------------------------------------

    expectEquivalent(parser,"1 N","1 kg*m/s^2");
    expectEquivalent(parser,"1 Pa","1 N/m^2");
    expectEquivalent(parser,"1 J","1 N*m");
    expectEquivalent(parser,"1 W","1 J/s");
    expectEquivalent(parser,"1 Hz","1 s^-1");

    //----------------------------------------------------------
    // Composite Units
    //----------------------------------------------------------

    expectEquivalent(
        parser,
        "1000 kg/m^3",
        "1 g/cm^3");

    expectEquivalent(
        parser,
        "1 Pa*s",
        "1 kg/(m*s)");

    expectEquivalent(
        parser,
        "1 m^2/s",
        "1 m*m/s");

    //----------------------------------------------------------
    // Angle
    //----------------------------------------------------------

    expectValue(parser,"180 deg",180.0);
    expectValue(parser,"3.141592653589793 rad",3.141592653589793);

    //----------------------------------------------------------
    // White Space
    //----------------------------------------------------------

    expectEquivalent(parser,"1kg","1 kg");
    expectEquivalent(parser,"1   kg","1 kg");
    expectEquivalent(parser,"1\tkg","1 kg");
    expectEquivalent(parser,"1 kg #comment","1 kg");

    //----------------------------------------------------------
    // Invalid
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse(""),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("kg"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("abc"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("1 abc"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("1 m+"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("1 ()"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("1 m^(1/0)"),
        std::runtime_error);

    ntic::lbm::test::printSummary();

    return 0;
}

