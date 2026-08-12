#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>

#include "unit.hpp"

namespace ntic::lbm::config
{

class UnitRegistry
{
public:
    //----------------------------------------------------------
    // Constructor
    //----------------------------------------------------------

    // Constructs a registry containing the built-in units:
    //
    //   m, kg, g, s, K
    //
    // Gram is included because SI mass prefixes are applied
    // to gram rather than kilogram.
    UnitRegistry();

    //----------------------------------------------------------
    // Registration
    //----------------------------------------------------------

    void registerUnit(const Unit& unit);

    //----------------------------------------------------------
    // Query
    //----------------------------------------------------------

    bool contains(const std::string& symbol) const;

    const Unit& find(const std::string& symbol) const;

    bool canApplyPrefix(const std::string& symbol) const;

    bool canApplyPrefix(const Unit& unit) const;

    std::size_t size() const;

    void clear();

private:
    void registerBuiltInUnits();

private:
    std::unordered_map<std::string, Unit> units_;
};

} // namespace ntic::lbm::config
