
#pragma once

#include <string>

#include "lexer.hpp"
#include "quantity.hpp"
#include "token_stream.hpp"
#include "unit_parser.hpp"

namespace ntic::lbm::config
{

class QuantityParser
{
public:
    QuantityParser(
        const UnitRegistry& unitRegistry,
        const PrefixRegistry& prefixRegistry);

    Quantity
    parse(
        const std::string& text) const;

private:
    double
    parseNumber(
        const Token& token) const;

private:
    Lexer lexer_;
    UnitParser unitParser_;
};

} // namespace ntic::lbm::config

