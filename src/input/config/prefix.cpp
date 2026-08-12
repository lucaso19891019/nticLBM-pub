#include "prefix.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace ntic::lbm::config
{

//==============================================================
// Constructors
//==============================================================

Prefix::Prefix()
    : symbol_(),
      name_(),
      factor_(1.0)
{
}

Prefix::Prefix(
    const std::string& symbol,
    const std::string& name,
    double factor)
    : symbol_(symbol),
      name_(name),
      factor_(factor)
{
    if (symbol_.empty())
    {
        throw std::invalid_argument(
            "Prefix: symbol cannot be empty.");
    }

    if (name_.empty())
    {
        throw std::invalid_argument(
            "Prefix '" + symbol_
            + "': name cannot be empty.");
    }

    if (!std::isfinite(factor_) || factor_ <= 0.0)
    {
        throw std::invalid_argument(
            "Prefix '" + symbol_
            + "': factor must be finite and greater than zero.");
    }
}

//==============================================================
// Access
//==============================================================

const std::string&
Prefix::symbol() const
{
    return symbol_;
}

const std::string&
Prefix::name() const
{
    return name_;
}

double
Prefix::factor() const
{
    return factor_;
}

//==============================================================
// Comparison
//==============================================================

bool
Prefix::operator==(const Prefix& rhs) const
{
    return symbol_ == rhs.symbol_
        && name_ == rhs.name_
        && factor_ == rhs.factor_;
}

bool
Prefix::operator!=(const Prefix& rhs) const
{
    return !(*this == rhs);
}

//==============================================================
// Utility
//==============================================================

std::string
Prefix::toString() const
{
    std::ostringstream oss;

    oss << symbol_
        << " ("
        << name_
        << ", factor="
        << factor_
        << ")";

    return oss.str();
}

} // namespace ntic::lbm::config
