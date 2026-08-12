#pragma once

#include <string>

#include "lexer.hpp"
#include "parameter_value.hpp"
#include "quantity_parser.hpp"
#include "token_stream.hpp"

namespace ntic::lbm::config
{

//==============================================================
// ParsedParameter
//==============================================================

struct ParsedParameter
{
    std::string name;

    ParameterValue value;
};

//==============================================================
// ParameterParser
//==============================================================

class ParameterParser
{
public:
    //----------------------------------------------------------
    // Constructor
    //----------------------------------------------------------

    ParameterParser(
        const UnitRegistry& unitRegistry,
        const PrefixRegistry& prefixRegistry);

    //----------------------------------------------------------
    // Parsing
    //----------------------------------------------------------

    ParsedParameter
    parse(
        const std::string& text) const;

private:
    //----------------------------------------------------------
    // Value parsing
    //----------------------------------------------------------

    bool
    tryParseBool(
        const Token& token,
        bool& value) const;

    std::string
    rebuildRemainingText(
        TokenStream& stream) const;

    //----------------------------------------------------------
    // Diagnostics
    //----------------------------------------------------------

    [[noreturn]]
    void
    throwError(
        const Token& token,
        const std::string& message) const;

private:
    Lexer lexer_;

    QuantityParser quantityParser_;
};

} // namespace ntic::lbm::config
