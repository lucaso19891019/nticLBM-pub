#pragma once

#include <string>

#include "physical_dimension.hpp"

namespace ntic::lbm::config
{

// Controls whether an SI prefix may be attached to a unit symbol.
//
// Examples:
//   m  : Allowed   -> mm, cm, km
//   Pa : Allowed   -> kPa, MPa
//   g  : Allowed   -> mg, kg, Mg
//   kg : Forbidden -> mkg, microkg, etc. are invalid
//   min: Forbidden -> mmin is invalid
enum class PrefixPolicy
{
    Forbidden,
    Allowed
};

class Unit
{
public:
    //----------------------------------------------------------
    // Constructors
    //----------------------------------------------------------

    // Creates an empty, dimensionless unit with scale 1.
    Unit();

    Unit(
        const std::string& symbol,
        const PhysicalDimension& dimension,
        double scale,
        PrefixPolicy prefixPolicy);

    //----------------------------------------------------------
    // Access
    //----------------------------------------------------------

    const std::string& symbol() const;

    const PhysicalDimension& dimension() const;

    // Multiplicative conversion factor to the corresponding
    // standard SI unit.
    //
    // Examples:
    //   m  -> 1
    //   g  -> 1e-3, because 1 g = 1e-3 kg
    //   Pa -> 1
    double scale() const;

    PrefixPolicy prefixPolicy() const;

    bool allowsPrefix() const;

    //----------------------------------------------------------
    // Comparison
    //----------------------------------------------------------

    bool operator==(const Unit& rhs) const;

    bool operator!=(const Unit& rhs) const;

    //----------------------------------------------------------
    // Utility
    //----------------------------------------------------------

    std::string toString() const;

    //----------------------------------------------------------
    // Built-in units
    //----------------------------------------------------------

    static Unit Meter();

    static Unit Kilogram();

    // Gram is registered as the prefix root for mass units.
    // Its scale is relative to the SI mass unit kilogram.
    static Unit Gram();

    static Unit Second();

    static Unit Kelvin();

    // Derived SI units

    static Unit Hertz();

    static Unit Newton();

    static Unit Pascal();

    static Unit Joule();

    static Unit Watt();

    // Engineering units

    static Unit Bar();

    static Unit Atmosphere();

    // Angle

    static Unit Radian();

    static Unit Degree();

private:
    std::string symbol_;

    PhysicalDimension dimension_;

    double scale_;

    PrefixPolicy prefixPolicy_;
};

} // namespace ntic::lbm::config
