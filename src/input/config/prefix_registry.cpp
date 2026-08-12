#include "prefix_registry.hpp"

#include <algorithm>
#include <stdexcept>

namespace ntic::lbm::config
{

//==============================================================
// Constructor
//==============================================================

PrefixRegistry::PrefixRegistry()
{
    registerSIPrefixes();
}

//==============================================================
// Built-in SI prefixes
//==============================================================

void
PrefixRegistry::registerSIPrefixes()
{
    registerPrefix(Prefix("Q",  "quetta", 1.0e30));
    registerPrefix(Prefix("R",  "ronna",  1.0e27));
    registerPrefix(Prefix("Y",  "yotta",  1.0e24));
    registerPrefix(Prefix("Z",  "zetta",  1.0e21));
    registerPrefix(Prefix("E",  "exa",    1.0e18));
    registerPrefix(Prefix("P",  "peta",   1.0e15));
    registerPrefix(Prefix("T",  "tera",   1.0e12));
    registerPrefix(Prefix("G",  "giga",   1.0e9));
    registerPrefix(Prefix("M",  "mega",   1.0e6));
    registerPrefix(Prefix("k",  "kilo",   1.0e3));
    registerPrefix(Prefix("h",  "hecto",  1.0e2));
    registerPrefix(Prefix("da", "deca",   1.0e1));

    registerPrefix(Prefix("d",  "deci",   1.0e-1));
    registerPrefix(Prefix("c",  "centi",  1.0e-2));
    registerPrefix(Prefix("m",  "milli",  1.0e-3));
    registerPrefix(Prefix("u",  "micro",  1.0e-6));
    registerPrefix(Prefix("n",  "nano",   1.0e-9));
    registerPrefix(Prefix("p",  "pico",   1.0e-12));
    registerPrefix(Prefix("f",  "femto",  1.0e-15));
    registerPrefix(Prefix("a",  "atto",   1.0e-18));
    registerPrefix(Prefix("z",  "zepto",  1.0e-21));
    registerPrefix(Prefix("y",  "yocto",  1.0e-24));
    registerPrefix(Prefix("r",  "ronto",  1.0e-27));
    registerPrefix(Prefix("q",  "quecto", 1.0e-30));
}

//==============================================================
// Registration
//==============================================================

void
PrefixRegistry::registerPrefix(
    const Prefix& prefix)
{
    const auto result =
        prefixes_.emplace(
            prefix.symbol(),
            prefix);

    if (!result.second)
    {
        throw std::runtime_error(
            "Prefix '" + prefix.symbol()
            + "' has already been registered.");
    }
}

//==============================================================
// Query
//==============================================================

bool
PrefixRegistry::contains(
    const std::string& symbol) const
{
    return prefixes_.find(symbol)
        != prefixes_.end();
}

const Prefix&
PrefixRegistry::find(
    const std::string& symbol) const
{
    const auto iter =
        prefixes_.find(symbol);

    if (iter == prefixes_.end())
    {
        throw std::runtime_error(
            "Unknown prefix: '" + symbol + "'.");
    }

    return iter->second;
}

std::size_t
PrefixRegistry::size() const
{
    return prefixes_.size();
}

//==============================================================
// Prefix matching
//==============================================================

std::vector<std::string>
PrefixRegistry::symbolsByDescendingLength() const
{
    std::vector<std::string> symbols;

    symbols.reserve(prefixes_.size());

    for (const auto& entry : prefixes_)
    {
        symbols.push_back(entry.first);
    }

    std::sort(
        symbols.begin(),
        symbols.end(),
        [](const std::string& lhs,
           const std::string& rhs)
        {
            if (lhs.size() != rhs.size())
            {
                return lhs.size() > rhs.size();
            }

            return lhs < rhs;
        });

    return symbols;
}

//==============================================================
// Maintenance
//==============================================================

void
PrefixRegistry::clear()
{
    prefixes_.clear();
}

} // namespace ntic::lbm::config
