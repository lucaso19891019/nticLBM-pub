#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include "prefix.hpp"

namespace ntic::lbm::config
{

class PrefixRegistry
{
public:
    //----------------------------------------------------------
    // Constructor
    //----------------------------------------------------------

    // Registers the supported SI prefixes automatically.
    PrefixRegistry();

    //----------------------------------------------------------
    // Registration
    //----------------------------------------------------------

    void registerPrefix(const Prefix& prefix);

    //----------------------------------------------------------
    // Query
    //----------------------------------------------------------

    bool contains(const std::string& symbol) const;

    const Prefix& find(const std::string& symbol) const;

    std::size_t size() const;

    //----------------------------------------------------------
    // Prefix matching
    //----------------------------------------------------------

    // Returns registered prefix symbols ordered from longest
    // to shortest. This is important because "da" must be
    // attempted before "d".
    std::vector<std::string>
    symbolsByDescendingLength() const;

    //----------------------------------------------------------
    // Maintenance
    //----------------------------------------------------------

    void clear();

private:
    void registerSIPrefixes();

private:
    std::unordered_map<std::string, Prefix> prefixes_;
};

} // namespace ntic::lbm::config
