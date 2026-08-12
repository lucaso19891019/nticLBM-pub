#include "si_unit_formatter.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

using namespace ntic::lbm;
using namespace ntic::lbm::config;

//----------------------------------------------------------
// Test utility
//----------------------------------------------------------

namespace
{

int passed = 0;
int total = 0;

void expectString(
    const PhysicalDimension& dimension,
    const std::string& expected)
{
    ++total;

    const std::string actual =
        SIUnitFormatter::format(
            dimension);

    if (actual != expected)
    {
        throw std::runtime_error(
            "Expected '"
            + expected
            + "', got '"
            + actual
            + "'.");
    }

    ++passed;
}

} // namespace

//----------------------------------------------------------
// Main
//----------------------------------------------------------

int main()
{
    try
    {
        //--------------------------------------------------
        // Dimensionless
        //--------------------------------------------------

        expectString(
            PhysicalDimension::Dimensionless(),
            "");

        //--------------------------------------------------
        // Base dimensions
        //--------------------------------------------------

        expectString(
            PhysicalDimension::Mass(),
            "kg");

        expectString(
            PhysicalDimension::Length(),
            "m");

        expectString(
            PhysicalDimension::Time(),
            "s");

        expectString(
            PhysicalDimension::Temperature(),
            "K");

        //--------------------------------------------------
        // Common derived dimensions
        //--------------------------------------------------

        expectString(
            PhysicalDimension::Area(),
            "m^2");

        expectString(
            PhysicalDimension::Volume(),
            "m^3");

        expectString(
            PhysicalDimension::Velocity(),
            "m/s");

        expectString(
            PhysicalDimension::Acceleration(),
            "m/s^2");

        expectString(
            PhysicalDimension::Force(),
            "(kg*m)/s^2");

        expectString(
            PhysicalDimension::Pressure(),
            "kg/(m*s^2)");

        expectString(
            PhysicalDimension::Density(),
            "kg/m^3");

        expectString(
            PhysicalDimension::DynamicViscosity(),
            "kg/(m*s)");

        expectString(
            PhysicalDimension::KinematicViscosity(),
            "m^2/s");

        expectString(
            PhysicalDimension::Energy(),
            "(kg*m^2)/s^2");

        expectString(
            PhysicalDimension::Power(),
            "(kg*m^2)/s^3");

        //--------------------------------------------------
        // Pure negative dimensions
        //--------------------------------------------------

        expectString(
            PhysicalDimension::Time().pow(
                common::Rational(-1)),
            "s^-1");

        expectString(
            PhysicalDimension::Length().pow(
                common::Rational(-2))
            * PhysicalDimension::Time().pow(
                common::Rational(-1)),
            "m^-2*s^-1");

        //--------------------------------------------------
        // Fractional exponent
        //--------------------------------------------------

        expectString(
            PhysicalDimension::Mass().pow(
                common::Rational(1, 2)),
            "(kg^1/2)");

        expectString(
            PhysicalDimension::Mass().pow(
                common::Rational(-1, 2)),
            "kg^-1/2");

        //--------------------------------------------------
        // Fractional numerator
        //--------------------------------------------------

        expectString(
            PhysicalDimension::Mass().pow(
                common::Rational(1, 2))
            * PhysicalDimension::Length().pow(
                common::Rational(-1)),
            "(kg^1/2)/m");

        //--------------------------------------------------
        // Fractional denominator
        //--------------------------------------------------

        expectString(
            PhysicalDimension::Mass()
            * PhysicalDimension::Length().pow(
                common::Rational(-1, 2)),
            "kg/(m^1/2)");

        //--------------------------------------------------
        // Fractional numerator and denominator
        //--------------------------------------------------

        expectString(
            PhysicalDimension::Mass().pow(
                common::Rational(1, 2))
            * PhysicalDimension::Length().pow(
                common::Rational(-3, 2)),
            "(kg^1/2)/(m^3/2)");

        //--------------------------------------------------
        // Result
        //--------------------------------------------------

        std::cout
            << "\n=================================\n"
            << "All tests passed.\n"
            << "Passed "
            << passed
            << " / "
            << total
            << "\n=================================\n";

        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "\nFAILED\n"
            << error.what()
            << '\n';

        return 1;
    }
}
