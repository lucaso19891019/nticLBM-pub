#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "lexer.hpp"
#include "prefix_registry.hpp"
#include "rational.hpp"
#include "token_stream.hpp"
#include "unit_expression.hpp"
#include "unit_registry.hpp"

namespace ntic::lbm::config
{

class UnitParser
{
public:
    //----------------------------------------------------------
    // Constructor
    //----------------------------------------------------------

    UnitParser(
        const UnitRegistry& unitRegistry,
        const PrefixRegistry& prefixRegistry);

    //----------------------------------------------------------
    // Parsing
    //----------------------------------------------------------

    // Parses a complete unit expression.
    //
    // Examples:
    //
    //   m
    //   cm
    //   kg/m^3
    //   kg/(m*s)
    //   kg m^-1 s^-2
    //   {[kg/(m*s)]}
    //
    // An empty expression is interpreted as dimensionless.
    UnitExpression
    parse(const std::string& text) const;

private:
    //----------------------------------------------------------
    // Recursive-descent grammar
    //----------------------------------------------------------

    UnitExpression
    parseExpression(
        TokenStream& stream,
        std::vector<char>& bracketStack) const;

    UnitExpression
    parseFactor(
        TokenStream& stream,
        std::vector<char>& bracketStack) const;

    UnitExpression
    parsePrimary(
        TokenStream& stream,
        std::vector<char>& bracketStack) const;

    //----------------------------------------------------------
    // Unit and prefix resolution
    //----------------------------------------------------------

    UnitExpression
    resolveUnit(const Token& token) const;

    //----------------------------------------------------------
    // Exponents
    //----------------------------------------------------------

    common::Rational
    parseExponent(
        TokenStream& stream) const;

    common::Rational
    parseRational(
        TokenStream& stream) const;

    int
    parseInteger(const Token& token) const;

    //----------------------------------------------------------
    // Grammar helpers
    //----------------------------------------------------------

    bool
    beginsFactor(const Token& token) const;

    bool
    isOpeningBracket(const Token& token) const;

    bool
    isClosingBracket(const Token& token) const;

    char
    matchingClosingBracket(char opening) const;

    int
    bracketLevel(char bracket) const;

    void
    validateBracketNesting(
        char opening,
        const std::vector<char>& bracketStack,
        const Token& token) const;

    //----------------------------------------------------------
    // Diagnostics
    //----------------------------------------------------------

    [[noreturn]]
    void
    throwError(
        const Token& token,
        const std::string& message) const;

private:
    const UnitRegistry& unitRegistry_;

    const PrefixRegistry& prefixRegistry_;

    Lexer lexer_;
};

} // namespace ntic::lbm::config
