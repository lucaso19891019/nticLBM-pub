#include "unit_registry.hpp"

#include <stdexcept>

namespace ntic::lbm::config
{

//==============================================================
// Constructor
//==============================================================

UnitRegistry::UnitRegistry()
{
    registerBuiltInUnits();
}

//==============================================================
// Built-in registration
//==============================================================

void
UnitRegistry::registerBuiltInUnits()
{
    registerUnit(Unit::Meter());
    registerUnit(Unit::Kilogram());
    registerUnit(Unit::Gram());
    registerUnit(Unit::Second());
    registerUnit(Unit::Kelvin());

    registerUnit(Unit::Radian());
    registerUnit(Unit::Degree());

    registerUnit(Unit::Pascal());
    registerUnit(Unit::Bar());
    registerUnit(Unit::Atmosphere());

    registerUnit(Unit::Newton());
    registerUnit(Unit::Joule());
    registerUnit(Unit::Watt());
    registerUnit(Unit::Hertz());
}

//==============================================================
// Registration
//==============================================================

void
UnitRegistry::registerUnit(const Unit& unit)
{
    const auto result =
        units_.emplace(unit.symbol(), unit);

    if (!result.second)
    {
        throw std::runtime_error(
            "Unit '" + unit.symbol()
            + "' has already been registered.");
    }
}

//==============================================================
// Query
//==============================================================

bool
UnitRegistry::contains(const std::string& symbol) const
{
    return units_.find(symbol) != units_.end();
}

const Unit&
UnitRegistry::find(const std::string& symbol) const
{
    const auto iter = units_.find(symbol);

    if (iter == units_.end())
    {
        throw std::runtime_error(
            "Unknown unit: '" + symbol + "'.");
    }

    return iter->second;
}

bool
UnitRegistry::canApplyPrefix(
    const std::string& symbol) const
{
    return find(symbol).allowsPrefix();
}

bool
UnitRegistry::canApplyPrefix(
    const Unit& unit) const
{
    return unit.allowsPrefix();
}

std::size_t
UnitRegistry::size() const
{
    return units_.size();
}

void
UnitRegistry::clear()
{
    units_.clear();
}

} // namespace ntic::lbm::config
