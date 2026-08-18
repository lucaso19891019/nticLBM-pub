#include "parameter_parser.hpp"
#include "physical_dimension.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

using namespace ntic::lbm::config;
using ntic::lbm::common::Rational;

namespace
{

const Quantity&
asQuantity(
    const ParsedParameter& parsed)
{
    NTIC_EXPECT_TRUE(
        std::holds_alternative<Quantity>(
            parsed.value));

    return std::get<Quantity>(
        parsed.value);
}

bool
asBool(
    const ParsedParameter& parsed)
{
    NTIC_EXPECT_TRUE(
        std::holds_alternative<bool>(
            parsed.value));

    return std::get<bool>(
        parsed.value);
}

const std::string&
asString(
    const ParsedParameter& parsed)
{
    NTIC_EXPECT_TRUE(
        std::holds_alternative<std::string>(
            parsed.value));

    return std::get<std::string>(
        parsed.value);
}

void
expectQuantity(
    const ParameterParser& parser,
    const std::string& input,
    const std::string& expectedName,
    double expectedValue,
    const PhysicalDimension& expectedDimension)
{
    const ParsedParameter parsed =
        parser.parse(input);

    NTIC_EXPECT_TRUE(
        parsed.name == expectedName);

    const Quantity& quantity =
        asQuantity(parsed);

    NTIC_EXPECT_TRUE(
    std::abs(
        quantity.value()-
        expectedValue)< 1e-12);

    NTIC_EXPECT_EQ(
        quantity.unit().dimension(),
        expectedDimension);
}

void
expectBool(
    const ParameterParser& parser,
    const std::string& input,
    const std::string& expectedName,
    bool expectedValue)
{
    const ParsedParameter parsed =
        parser.parse(input);

    NTIC_EXPECT_TRUE(
        parsed.name == expectedName);

    NTIC_EXPECT_TRUE(
        asBool(parsed) == expectedValue);
}

void
expectString(
    const ParameterParser& parser,
    const std::string& input,
    const std::string& expectedName,
    const std::string& expectedValue)
{
    const ParsedParameter parsed =
        parser.parse(input);

    NTIC_EXPECT_TRUE(
        parsed.name == expectedName);

    NTIC_EXPECT_TRUE(
        asString(parsed) == expectedValue);
}

} // namespace

int main()
{
    UnitRegistry unitRegistry;
    PrefixRegistry prefixRegistry;

    ParameterParser parser(
        unitRegistry,
        prefixRegistry);

    //----------------------------------------------------------
    // Dimensionless quantities
    //----------------------------------------------------------

    expectQuantity(
        parser,
        "reynolds = 1000",
        "reynolds",
        1000.0,
        PhysicalDimension::Dimensionless());

    expectQuantity(
        parser,
        "omega=1.75",
        "omega",
        1.75,
        PhysicalDimension::Dimensionless());

    expectQuantity(
        parser,
        "tolerance = 1.0e-8",
        "tolerance",
        1.0e-8,
        PhysicalDimension::Dimensionless());

    expectQuantity(
        parser,
        "negative_value = -2.5",
        "negative_value",
        -2.5,
        PhysicalDimension::Dimensionless());

    expectQuantity(
        parser,
        "positive_value = +3.5",
        "positive_value",
        3.5,
        PhysicalDimension::Dimensionless());

    expectQuantity(
        parser,
        "zero = 0",
        "zero",
        0.0,
        PhysicalDimension::Dimensionless());

    expectQuantity(
        parser,
        "one = 1",
        "one",
        1.0,
        PhysicalDimension::Dimensionless());

    //----------------------------------------------------------
    // Base SI quantities
    //----------------------------------------------------------

    expectQuantity(
        parser,
        "length = 2 m",
        "length",
        2.0,
        PhysicalDimension::Length());

    expectQuantity(
        parser,
        "mass = 3 kg",
        "mass",
        3.0,
        PhysicalDimension::Mass());

    expectQuantity(
        parser,
        "time = 4 s",
        "time",
        4.0,
        PhysicalDimension::Time());

    expectQuantity(
        parser,
        "temperature = 300 K",
        "temperature",
        300.0,
        PhysicalDimension::Temperature());

    //----------------------------------------------------------
    // Prefixes
    //----------------------------------------------------------

    expectQuantity(
        parser,
        "dx = 5 mm",
        "dx",
        5.0,
        PhysicalDimension::Length());

    expectQuantity(
        parser,
        "radius = 2.5 cm",
        "radius",
        2.5,
        PhysicalDimension::Length());

    expectQuantity(
        parser,
        "distance = 4 km",
        "distance",
        4.0,
        PhysicalDimension::Length());

    expectQuantity(
        parser,
        "dt = 10 ms",
        "dt",
        10.0,
        PhysicalDimension::Time());

    //----------------------------------------------------------
    // Compound quantities
    //----------------------------------------------------------

    expectQuantity(
        parser,
        "velocity = 12 m/s",
        "velocity",
        12.0,
        PhysicalDimension::Velocity());

    expectQuantity(
        parser,
        "acceleration = 9.81 m/s^2",
        "acceleration",
        9.81,
        PhysicalDimension::Acceleration());

    expectQuantity(
        parser,
        "density = 998 kg/m^3",
        "density",
        998.0,
        PhysicalDimension::Density());

    expectQuantity(
        parser,
        "area = 2.5 m^2",
        "area",
        2.5,
        PhysicalDimension::Area());

    expectQuantity(
        parser,
        "volume = 4.0 m^3",
        "volume",
        4.0,
        PhysicalDimension::Volume());

    //----------------------------------------------------------
    // Derived units
    //----------------------------------------------------------

    expectQuantity(
        parser,
        "pressure = 101325 Pa",
        "pressure",
        101325.0,
        PhysicalDimension::Pressure());

    expectQuantity(
        parser,
        "force = 12 N",
        "force",
        12.0,
        PhysicalDimension::Force());

    expectQuantity(
        parser,
        "energy = 15 J",
        "energy",
        15.0,
        PhysicalDimension::Energy());

    expectQuantity(
        parser,
        "power = 20 W",
        "power",
        20.0,
        PhysicalDimension::Power());

    expectQuantity(
        parser,
        "frequency = 50 Hz",
        "frequency",
        50.0,
        PhysicalDimension(
    	    Rational(0),Rational(0),
    	    Rational(-1),Rational(0)));

    //----------------------------------------------------------
    // Separators
    //----------------------------------------------------------

    expectQuantity(
        parser,
        "a=1",
        "a",
        1.0,
        PhysicalDimension::Dimensionless());

    expectQuantity(
        parser,
        "b =1",
        "b",
        1.0,
        PhysicalDimension::Dimensionless());

    expectQuantity(
        parser,
        "c= 1",
        "c",
        1.0,
        PhysicalDimension::Dimensionless());

    expectQuantity(
        parser,
        "d = 1",
        "d",
        1.0,
        PhysicalDimension::Dimensionless());

    expectQuantity(
        parser,
        "e:1",
        "e",
        1.0,
        PhysicalDimension::Dimensionless());

    expectQuantity(
        parser,
        "f : 1",
        "f",
        1.0,
        PhysicalDimension::Dimensionless());

    expectQuantity(
        parser,
        "g 1",
        "g",
        1.0,
        PhysicalDimension::Dimensionless());

    //----------------------------------------------------------
    // Whitespace
    //----------------------------------------------------------

    expectQuantity(
        parser,
        "    density    =    998    kg / m ^ 3    ",
        "density",
        998.0,
        PhysicalDimension::Density());

    expectQuantity(
        parser,
        "\tdt\t=\t1.0e-4\ts\t",
        "dt",
        1.0e-4,
        PhysicalDimension::Time());

    expectString(
        parser,
        "solver\t=\tbicgstab",
        "solver",
        "bicgstab");

    //----------------------------------------------------------
    // Bool: true values
    //----------------------------------------------------------

    const std::vector<std::string> trueValues =
    {
        "true",
        "TRUE",
        "True",
        "tRuE",
        "on",
        "ON",
        "On",
        "oN",
        "yes",
        "YES",
        "Yes",
        "yEs"
    };

    for (const std::string& value : trueValues)
    {
        expectBool(
            parser,
            "enabled = " + value,
            "enabled",
            true);
    }

    //----------------------------------------------------------
    // Bool: false values
    //----------------------------------------------------------

    const std::vector<std::string> falseValues =
    {
        "false",
        "FALSE",
        "False",
        "fAlSe",
        "off",
        "OFF",
        "Off",
        "oFf",
        "no",
        "NO",
        "No",
        "nO"
    };

    for (const std::string& value : falseValues)
    {
        expectBool(
            parser,
            "enabled = " + value,
            "enabled",
            false);
    }

    //----------------------------------------------------------
    // Numeric 0 and 1 are quantities, not bools
    //----------------------------------------------------------

    {
        const ParsedParameter parsed =
            parser.parse("enabled = 1");

        NTIC_EXPECT_TRUE(
            parsed.name == "enabled");

        NTIC_EXPECT_TRUE(
            std::holds_alternative<Quantity>(
                parsed.value));

        NTIC_EXPECT_TRUE(
            !std::holds_alternative<bool>(
                parsed.value));

        NTIC_EXPECT_TRUE(
            std::abs(
                std::get<Quantity>(
                    parsed.value).value()-
            1.0)< 1e-12);
    }

    {
        const ParsedParameter parsed =
            parser.parse("enabled = 0");

        NTIC_EXPECT_TRUE(
            parsed.name == "enabled");

        NTIC_EXPECT_TRUE(
            std::holds_alternative<Quantity>(
                parsed.value));

        NTIC_EXPECT_TRUE(
            !std::holds_alternative<bool>(
                parsed.value));

        NTIC_EXPECT_TRUE(
            std::abs(
                    std::get<Quantity>(
                        parsed.value).value()-
            0.0)< 1e-12);
    }

    //----------------------------------------------------------
    // Unquoted identifier strings
    //----------------------------------------------------------

    expectString(
        parser,
        "solver = bicgstab",
        "solver",
        "bicgstab");

    expectString(
        parser,
        "solver = gmres",
        "solver",
        "gmres");

    expectString(
        parser,
        "backend = cuda",
        "backend",
        "cuda");

    expectString(
        parser,
        "backend = kokkos",
        "backend",
        "kokkos");

    expectString(
        parser,
        "precision = double",
        "precision",
        "double");

    expectString(
        parser,
        "format = vtk",
        "format",
        "vtk");

    expectString(
        parser,
        "collision_model = cumulant",
        "collision_model",
        "cumulant");

    //----------------------------------------------------------
    // Quoted strings
    //----------------------------------------------------------

    expectString(
        parser,
        "solver = \"bicgstab\"",
        "solver",
        "bicgstab");

    expectString(
        parser,
        "solver = bicgstab gmres",
        "solver",
        "bicgstab gmres");

    expectString(
        parser,
        "mesh = \"sphere.stl\"",
        "mesh",
        "sphere.stl");

    expectString(
        parser,
        "title = \"LBM Solver\"",
        "title",
        "LBM Solver");

    expectString(
        parser,
        "empty_string = \"\"",
        "empty_string",
        "");

    expectString(
        parser,
        "path = \"/tmp/case/input.stl\"",
        "path",
        "/tmp/case/input.stl");

    expectString(
        parser,
        "path = \"/tmp/part#1.stl\"",
        "path",
        "/tmp/part#1.stl");

    //----------------------------------------------------------
    // Multi-token unquoted string values
    //----------------------------------------------------------

    expectString(
        parser,
        "collision_model = central moment",
        "collision_model",
        "central moment");

    expectString(
        parser,
        "inlet_boundary = ZH-vel",
        "inlet_boundary",
        "ZH-vel");

    expectString(
        parser,
        "outlet_boundary = Zou-He pressure",
        "outlet_boundary",
        "Zou-He pressure");

    expectString(
        parser,
        "wall_boundary = immersed boundary",
        "wall_boundary",
        "immersed boundary");

    //----------------------------------------------------------
    // Quoted bool-like and number-like strings
    //----------------------------------------------------------

    expectString(
        parser,
        "value = \"true\"",
        "value",
        "true");

    expectString(
        parser,
        "value = \"false\"",
        "value",
        "false");

    expectString(
        parser,
        "value = \"on\"",
        "value",
        "on");

    expectString(
        parser,
        "value = \"1\"",
        "value",
        "1");

    expectString(
        parser,
        "value = \"998 kg/m^3\"",
        "value",
        "998 kg/m^3");

    //----------------------------------------------------------
    // Comments
    //----------------------------------------------------------

    expectQuantity(
        parser,
        "density = 998 kg/m^3 # water density",
        "density",
        998.0,
        PhysicalDimension::Density());

    expectBool(
        parser,
        "enabled = true # enable output",
        "enabled",
        true);

    expectString(
        parser,
        "solver = bicgstab # linear solver",
        "solver",
        "bicgstab");

    expectString(
        parser,
        "path = \"/tmp/a#1.stl\" # mesh path",
        "path",
        "/tmp/a#1.stl");

    //----------------------------------------------------------
    // Parser does not require parameter registration
    //----------------------------------------------------------

    expectQuantity(
        parser,
        "completely_new_parameter = 123",
        "completely_new_parameter",
        123.0,
        PhysicalDimension::Dimensionless());

    expectBool(
        parser,
        "another_unknown_parameter = yes",
        "another_unknown_parameter",
        true);

    expectString(
        parser,
        "custom_backend = experimental",
        "custom_backend",
        "experimental");

    //----------------------------------------------------------
    // Invalid: empty or missing fields
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse(""),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("   "),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("# comment only"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("="),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse(":"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("density ="),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("density :"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("density"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("= 1"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse(": 1"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("123 = 1"),
        std::runtime_error);

    //----------------------------------------------------------
    // Invalid: extra tokens after scalar values
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse("enabled = true false"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("enabled = true extra"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("solver = \"bicgstab\" extra"),
        std::runtime_error);

    //----------------------------------------------------------
    // Invalid: malformed quantity/unit expressions
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        parser.parse("density = 1 unknown_unit"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("density = 1 m+"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("density = 1 m+s"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("density = 1 ()"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("density = 1 m^^2"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("density = 1 m/"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        parser.parse("density = 1 /m"),
        std::runtime_error);

    //----------------------------------------------------------
    // Variant distinction regression tests
    //----------------------------------------------------------

    {
        const ParsedParameter boolParsed =
            parser.parse("value = true");

        const ParsedParameter stringParsed =
            parser.parse("value = \"true\"");

        const ParsedParameter quantityParsed =
            parser.parse("value = 1");

        NTIC_EXPECT_TRUE(
            std::holds_alternative<bool>(
                boolParsed.value));

        NTIC_EXPECT_TRUE(
            std::holds_alternative<std::string>(
                stringParsed.value));

        NTIC_EXPECT_TRUE(
            std::holds_alternative<Quantity>(
                quantityParsed.value));
    }

    //----------------------------------------------------------

    ntic::lbm::test::printSummary();

    return 0;
}

