#pragma once

#include <string>

namespace ntic::lbm::config
{

class Prefix
{
public:
    //----------------------------------------------------------
    // Constructors
    //----------------------------------------------------------

    Prefix();

    Prefix(
        const std::string& symbol,
        const std::string& name,
        double factor);

    //----------------------------------------------------------
    // Access
    //----------------------------------------------------------

    const std::string& symbol() const;

    const std::string& name() const;

    double factor() const;

    //----------------------------------------------------------
    // Comparison
    //----------------------------------------------------------

    bool operator==(const Prefix& rhs) const;

    bool operator!=(const Prefix& rhs) const;

    //----------------------------------------------------------
    // Utility
    //----------------------------------------------------------

    std::string toString() const;

private:
    std::string symbol_;

    std::string name_;

    double factor_;
};

} // namespace ntic::lbm::config
