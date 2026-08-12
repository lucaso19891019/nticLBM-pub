#include "unit.hpp"
#include "unit_registry.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>


using ntic::lbm::common::Rational;

using ntic::lbm::config::PhysicalDimension;
using ntic::lbm::config::PrefixPolicy;
using ntic::lbm::config::Unit;
using ntic::lbm::config::UnitRegistry;

int main()
{
    //----------------------------------------------------------
    // Default Unit
    //----------------------------------------------------------

    {
        Unit unit;

        NTIC_EXPECT_TRUE(unit.symbol().empty());

        NTIC_EXPECT_EQ(
            unit.dimension(),
            PhysicalDimension::Dimensionless());

        NTIC_EXPECT_TRUE(unit.scale() == 1.0);

        NTIC_EXPECT_TRUE(
            unit.prefixPolicy() == PrefixPolicy::Forbidden);

        NTIC_EXPECT_TRUE(!unit.allowsPrefix());
    }

    //----------------------------------------------------------
    // General Constructor
    //----------------------------------------------------------

    {
        Unit unit(
            "custom_length",
            PhysicalDimension::Length(),
            2.5,
            PrefixPolicy::Allowed);

        NTIC_EXPECT_TRUE(
            unit.symbol() == "custom_length");

        NTIC_EXPECT_EQ(
            unit.dimension(),
            PhysicalDimension::Length());

        NTIC_EXPECT_TRUE(
            unit.scale() == 2.5);

        NTIC_EXPECT_TRUE(
            unit.prefixPolicy() == PrefixPolicy::Allowed);

        NTIC_EXPECT_TRUE(
            unit.allowsPrefix());
    }

    //----------------------------------------------------------
    // Invalid Constructor Arguments
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        Unit(
            "",
            PhysicalDimension::Length(),
            1.0,
            PrefixPolicy::Allowed),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        Unit(
            "invalid_zero",
            PhysicalDimension::Length(),
            0.0,
            PrefixPolicy::Allowed),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        Unit(
            "invalid_negative",
            PhysicalDimension::Length(),
            -1.0,
            PrefixPolicy::Allowed),
        std::invalid_argument);

    //----------------------------------------------------------
    // Meter
    //----------------------------------------------------------

    {
        const Unit meter = Unit::Meter();

        NTIC_EXPECT_TRUE(
            meter.symbol() == "m");

        NTIC_EXPECT_EQ(
            meter.dimension(),
            PhysicalDimension::Length());

        NTIC_EXPECT_TRUE(
            meter.scale() == 1.0);

        NTIC_EXPECT_TRUE(
            meter.prefixPolicy() == PrefixPolicy::Allowed);

        NTIC_EXPECT_TRUE(
            meter.allowsPrefix());
    }

    //----------------------------------------------------------
    // Kilogram
    //----------------------------------------------------------

    {
        const Unit kilogram = Unit::Kilogram();

        NTIC_EXPECT_TRUE(
            kilogram.symbol() == "kg");

        NTIC_EXPECT_EQ(
            kilogram.dimension(),
            PhysicalDimension::Mass());

        NTIC_EXPECT_TRUE(
            kilogram.scale() == 1.0);

        NTIC_EXPECT_TRUE(
            kilogram.prefixPolicy()
            == PrefixPolicy::Forbidden);

        NTIC_EXPECT_TRUE(
            !kilogram.allowsPrefix());
    }

    //----------------------------------------------------------
    // Gram
    //----------------------------------------------------------

    {
        const Unit gram = Unit::Gram();

        NTIC_EXPECT_TRUE(
            gram.symbol() == "g");

        NTIC_EXPECT_EQ(
            gram.dimension(),
            PhysicalDimension::Mass());

        NTIC_EXPECT_TRUE(
            gram.scale() == 1.0e-3);

        NTIC_EXPECT_TRUE(
            gram.prefixPolicy() == PrefixPolicy::Allowed);

        NTIC_EXPECT_TRUE(
            gram.allowsPrefix());
    }

    //----------------------------------------------------------
    // Second
    //----------------------------------------------------------

    {
        const Unit second = Unit::Second();

        NTIC_EXPECT_TRUE(
            second.symbol() == "s");

        NTIC_EXPECT_EQ(
            second.dimension(),
            PhysicalDimension::Time());

        NTIC_EXPECT_TRUE(
            second.scale() == 1.0);

        NTIC_EXPECT_TRUE(
            second.allowsPrefix());
    }

    //----------------------------------------------------------
    // Kelvin
    //----------------------------------------------------------

    {
        const Unit kelvin = Unit::Kelvin();

        NTIC_EXPECT_TRUE(
            kelvin.symbol() == "K");

        NTIC_EXPECT_EQ(
            kelvin.dimension(),
            PhysicalDimension::Temperature());

        NTIC_EXPECT_TRUE(
            kelvin.scale() == 1.0);

        NTIC_EXPECT_TRUE(
            kelvin.allowsPrefix());
    }

    //----------------------------------------------------------
    // Built-in Derived Units
    //----------------------------------------------------------

    {
        const Unit hertz = Unit::Hertz();

        NTIC_EXPECT_TRUE(hertz.symbol() == "Hz");

        NTIC_EXPECT_EQ(
            hertz.dimension(),
            PhysicalDimension::Time().pow(
                Rational(-1)));

        NTIC_EXPECT_TRUE(
            hertz.allowsPrefix());
    }

    {
        const Unit newton = Unit::Newton();

        NTIC_EXPECT_TRUE(newton.symbol() == "N");

        NTIC_EXPECT_EQ(
            newton.dimension(),
            PhysicalDimension::Force());

        NTIC_EXPECT_TRUE(
            newton.allowsPrefix());
    }

    {
        const Unit pascal = Unit::Pascal();

        NTIC_EXPECT_TRUE(pascal.symbol() == "Pa");

        NTIC_EXPECT_EQ(
            pascal.dimension(),
            PhysicalDimension::Pressure());

        NTIC_EXPECT_TRUE(
            pascal.allowsPrefix());
    }

    {
        const Unit joule = Unit::Joule();

        NTIC_EXPECT_TRUE(joule.symbol() == "J");

        NTIC_EXPECT_EQ(
            joule.dimension(),
            PhysicalDimension::Energy());

        NTIC_EXPECT_TRUE(
            joule.allowsPrefix());
    }

    {
        const Unit watt = Unit::Watt();

        NTIC_EXPECT_TRUE(watt.symbol() == "W");

        NTIC_EXPECT_EQ(
            watt.dimension(),
            PhysicalDimension::Power());

        NTIC_EXPECT_TRUE(
            watt.allowsPrefix());
    }

    //----------------------------------------------------------
    // Engineering Units
    //----------------------------------------------------------

    {
        const Unit bar = Unit::Bar();

        NTIC_EXPECT_TRUE(bar.symbol() == "bar");

        NTIC_EXPECT_EQ(
            bar.dimension(),
            PhysicalDimension::Pressure());

        NTIC_EXPECT_TRUE(
            !bar.allowsPrefix());
    }

    {
        const Unit atm = Unit::Atmosphere();

        NTIC_EXPECT_TRUE(atm.symbol() == "atm");

        NTIC_EXPECT_EQ(
            atm.dimension(),
            PhysicalDimension::Pressure());

        NTIC_EXPECT_TRUE(
            !atm.allowsPrefix());
    }

    //----------------------------------------------------------
    // Angle Units
    //----------------------------------------------------------

    {
        const Unit rad = Unit::Radian();

        NTIC_EXPECT_TRUE(rad.symbol() == "rad");

        NTIC_EXPECT_TRUE(
            rad.dimension()
            ==
            PhysicalDimension::Dimensionless());

        NTIC_EXPECT_TRUE(
            !rad.allowsPrefix());
    }

    {
        const Unit deg = Unit::Degree();

        NTIC_EXPECT_TRUE(deg.symbol() == "deg");

        NTIC_EXPECT_TRUE(
            deg.dimension()
            ==
            PhysicalDimension::Dimensionless());

        NTIC_EXPECT_TRUE(
            !deg.allowsPrefix());

        NTIC_EXPECT_TRUE(
            std::abs(
                deg.scale()
                -
                3.141592653589793238462643383279502884
                /180.0)
            <
            1e-12);
    }

    //----------------------------------------------------------
    // Unit Comparison
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(
        Unit::Meter() == Unit::Meter());

    NTIC_EXPECT_TRUE(
        Unit::Meter() != Unit::Second());

    NTIC_EXPECT_TRUE(
        Unit::Kilogram() != Unit::Gram());

    {
        Unit unitA(
            "x",
            PhysicalDimension::Length(),
            1.0,
            PrefixPolicy::Allowed);

        Unit unitB(
            "x",
            PhysicalDimension::Length(),
            1.0,
            PrefixPolicy::Allowed);

        NTIC_EXPECT_TRUE(
            unitA == unitB);
    }

    //----------------------------------------------------------
    // Unit toString
    //----------------------------------------------------------

    {
        const std::string text =
            Unit::Meter().toString();

        NTIC_EXPECT_TRUE(
            text.find("m") != std::string::npos);

        NTIC_EXPECT_TRUE(
            text.find("[L]") != std::string::npos);

        NTIC_EXPECT_TRUE(
            text.find("scale=1") != std::string::npos);

        NTIC_EXPECT_TRUE(
            text.find("prefix=allowed")
            != std::string::npos);
    }

    //----------------------------------------------------------
    // UnitRegistry Built-in Units
    //----------------------------------------------------------

    UnitRegistry registry;

    NTIC_EXPECT_TRUE(
        registry.size() == 14);

    NTIC_EXPECT_TRUE(
        registry.contains("m"));

    NTIC_EXPECT_TRUE(
        registry.contains("kg"));

    NTIC_EXPECT_TRUE(
        registry.contains("g"));

    NTIC_EXPECT_TRUE(
        registry.contains("s"));

    NTIC_EXPECT_TRUE(
        registry.contains("K"));

    NTIC_EXPECT_TRUE(
        registry.contains("Hz"));

    NTIC_EXPECT_TRUE(
        registry.contains("N"));

    NTIC_EXPECT_TRUE(
        registry.contains("Pa"));

    NTIC_EXPECT_TRUE(
        registry.contains("J"));

    NTIC_EXPECT_TRUE(
        registry.contains("W"));

    NTIC_EXPECT_TRUE(
        registry.contains("bar"));

    NTIC_EXPECT_TRUE(
        registry.contains("atm"));

    NTIC_EXPECT_TRUE(
        registry.contains("rad"));

    NTIC_EXPECT_TRUE(
        registry.contains("deg"));

    //----------------------------------------------------------
    // UnitRegistry find()
    //----------------------------------------------------------

    {
        const Unit& meter =
            registry.find("m");

        NTIC_EXPECT_EQ(
            meter,
            Unit::Meter());
    }

    {
        const Unit& kilogram =
            registry.find("kg");

        NTIC_EXPECT_EQ(
            kilogram,
            Unit::Kilogram());
    }

    {
        const Unit& gram =
            registry.find("g");

        NTIC_EXPECT_EQ(
            gram,
            Unit::Gram());
    }

    {
        const Unit& pa =
            registry.find("Pa");

        NTIC_EXPECT_EQ(
            pa,
            Unit::Pascal());
    }

    {
        const Unit& deg =
            registry.find("deg");

        NTIC_EXPECT_EQ(
            deg,
            Unit::Degree());
    }

    //----------------------------------------------------------
    // UnitRegistry Prefix Policy
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(
        registry.canApplyPrefix("m"));

    NTIC_EXPECT_TRUE(
        registry.canApplyPrefix("g"));

    NTIC_EXPECT_TRUE(
        registry.canApplyPrefix("s"));

    NTIC_EXPECT_TRUE(
        registry.canApplyPrefix("K"));

    NTIC_EXPECT_TRUE(
        !registry.canApplyPrefix("kg"));

    NTIC_EXPECT_TRUE(
        registry.canApplyPrefix(Unit::Meter()));

    NTIC_EXPECT_TRUE(
        !registry.canApplyPrefix(Unit::Kilogram()));

    NTIC_EXPECT_TRUE(
        registry.canApplyPrefix("Pa"));

    NTIC_EXPECT_TRUE(
        registry.canApplyPrefix("Hz"));

    NTIC_EXPECT_TRUE(
        !registry.canApplyPrefix("bar"));

    NTIC_EXPECT_TRUE(
        !registry.canApplyPrefix("deg"));

    //----------------------------------------------------------
    // Register Custom Unit
    //----------------------------------------------------------

    {
        Unit minute(
            "min",
            PhysicalDimension::Time(),
            60.0,
            PrefixPolicy::Forbidden);

        const std::size_t oldSize =
            registry.size();

        registry.registerUnit(minute);

        NTIC_EXPECT_TRUE(
            registry.size() == oldSize + 1);

        NTIC_EXPECT_TRUE(
            registry.contains("min"));

        NTIC_EXPECT_EQ(
            registry.find("min"),
            minute);

        NTIC_EXPECT_TRUE(
            !registry.canApplyPrefix("min"));
    }

    //----------------------------------------------------------
    // Duplicate Registration
    //----------------------------------------------------------

    {
        const std::size_t oldSize =
            registry.size();

        NTIC_EXPECT_THROW(
            registry.registerUnit(Unit::Meter()),
            std::runtime_error);

        NTIC_EXPECT_TRUE(
            registry.size() == oldSize);

        NTIC_EXPECT_TRUE(
            registry.contains("m"));

        NTIC_EXPECT_EQ(
            registry.find("m"),
            Unit::Meter());
    }

    //----------------------------------------------------------
    // Unknown Unit
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        registry.find("unknown_unit"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        registry.canApplyPrefix("unknown_unit"),
        std::runtime_error);

    //----------------------------------------------------------
    // clear()
    //----------------------------------------------------------

    registry.clear();

    NTIC_EXPECT_TRUE(
        registry.size() == 0);

    NTIC_EXPECT_TRUE(
        !registry.contains("m"));

    NTIC_EXPECT_TRUE(
        !registry.contains("kg"));

    NTIC_EXPECT_TRUE(
        !registry.contains("g"));

    NTIC_EXPECT_TRUE(
        !registry.contains("s"));

    NTIC_EXPECT_TRUE(
        !registry.contains("K"));

    NTIC_EXPECT_THROW(
        registry.find("m"),
        std::runtime_error);

    //----------------------------------------------------------
    // Register After clear()
    //----------------------------------------------------------

    registry.registerUnit(Unit::Meter());

    NTIC_EXPECT_TRUE(
        registry.size() == 1);

    NTIC_EXPECT_TRUE(
        registry.contains("m"));

    NTIC_EXPECT_EQ(
        registry.find("m"),
        Unit::Meter());

    //----------------------------------------------------------

    ntic::lbm::test::printSummary();

    return 0;
}
