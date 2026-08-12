#include "unit_parser.hpp"

#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace ntic::lbm::config
{

//==============================================================
// Constructor
//==============================================================

UnitParser::UnitParser(
    const UnitRegistry& unitRegistry,
    const PrefixRegistry& prefixRegistry)
    : unitRegistry_(unitRegistry),
      prefixRegistry_(prefixRegistry),
      lexer_()
{
}

//==============================================================
// Public parsing interface
//==============================================================

UnitExpression
UnitParser::parse(
    const std::string& text) const
{
    const std::vector<Token> tokens =
        lexer_.tokenize(text);

    TokenStream stream(tokens);

    // Empty unit text means a dimensionless unit expression.
    if (stream.atEnd())
    {
        return UnitExpression();
    }

    std::vector<char> bracketStack;

    UnitExpression result =
        parseExpression(
            stream,
            bracketStack);

    if (!bracketStack.empty())
    {
        throw std::runtime_error(
            "UnitParser: internal bracket stack is not empty.");
    }

    if (!stream.atEnd())
    {
        throwError(
            stream.current(),
            "unexpected token '"
            + stream.current().text()
            + "' after unit expression");
    }

    return result;
}

//==============================================================
// Recursive-descent grammar
//==============================================================

UnitExpression
UnitParser::parseExpression(
    TokenStream& stream,
    std::vector<char>& bracketStack) const
{
    UnitExpression result =
        parseFactor(
            stream,
            bracketStack);

    while (true)
    {
        //------------------------------------------------------
        // Explicit multiplication
        //------------------------------------------------------

        if (stream.match("*"))
        {
            if (!beginsFactor(stream.current()))
            {
                throwError(
                    stream.current(),
                    "expected a unit or bracketed expression "
                    "after '*'");
            }

            result =
                result
                * parseFactor(
                    stream,
                    bracketStack);

            continue;
        }

        //------------------------------------------------------
        // Division
        //------------------------------------------------------

        if (stream.match("/"))
        {
            if (!beginsFactor(stream.current()))
            {
                throwError(
                    stream.current(),
                    "expected a unit or bracketed expression "
                    "after '/'");
            }

            result =
                result
                / parseFactor(
                    stream,
                    bracketStack);

            continue;
        }

        //------------------------------------------------------
        // Implicit multiplication
        //
        // Examples:
        //
        //   kg m
        //   kg(m)
        //   (kg)m
        //   (kg)(m)
        //------------------------------------------------------

        if (beginsFactor(stream.current()))
        {
            result =
                result
                * parseFactor(
                    stream,
                    bracketStack);

            continue;
        }

        break;
    }

    return result;
}

UnitExpression
UnitParser::parseFactor(
    TokenStream& stream,
    std::vector<char>& bracketStack) const
{
    UnitExpression result =
        parsePrimary(
            stream,
            bracketStack);

    if (stream.match("^"))
    {
        const common::Rational exponent =
            parseExponent(stream);

        result =
            result.pow(exponent);
    }

    return result;
}

UnitExpression
UnitParser::parsePrimary(
    TokenStream& stream,
    std::vector<char>& bracketStack) const
{
    //----------------------------------------------------------
    // Unit identifier
    //----------------------------------------------------------

    if (stream.check(TokenType::Identifier))
    {
        const Token token =
            stream.current();

        stream.consume();

        return resolveUnit(token);
    }

    //----------------------------------------------------------
    // Bracketed unit expression
    //----------------------------------------------------------

    if (isOpeningBracket(stream.current()))
    {
        const Token openingToken =
            stream.current();

        const char opening =
            openingToken.text()[0];

        validateBracketNesting(
            opening,
            bracketStack,
            openingToken);

        bracketStack.push_back(opening);

        stream.consume();

        if (stream.atEnd())
        {
            throwError(
                openingToken,
                std::string("missing closing bracket '")
                + matchingClosingBracket(opening)
                + "'");
        }

        if (isClosingBracket(stream.current()))
        {
            throwError(
                stream.current(),
                "empty bracketed unit expression is not allowed");
        }

        UnitExpression result =
            parseExpression(
                stream,
                bracketStack);

        const char expectedClosing =
            matchingClosingBracket(opening);

        if (!stream.check(
                std::string(1, expectedClosing)))
        {
            if (stream.atEnd())
            {
                throwError(
                    openingToken,
                    std::string("missing closing bracket '")
                    + expectedClosing
                    + "'");
            }

            if (isClosingBracket(stream.current()))
            {
                throwError(
                    stream.current(),
                    std::string("mismatched closing bracket '")
                    + stream.current().text()
                    + "'; expected '"
                    + expectedClosing
                    + "'");
            }

            throwError(
                stream.current(),
                std::string("expected closing bracket '")
                + expectedClosing
                + "'");
        }

        stream.consume();

        bracketStack.pop_back();

        return result;
    }

    //----------------------------------------------------------
    // Unexpected closing bracket
    //----------------------------------------------------------

    if (isClosingBracket(stream.current()))
    {
        throwError(
            stream.current(),
            "unexpected closing bracket '"
            + stream.current().text()
            + "'");
    }

    if (stream.atEnd())
    {
        throwError(
            stream.current(),
            "expected a unit expression");
    }

    throwError(
        stream.current(),
        "expected a unit identifier or opening bracket");
}

//==============================================================
// Unit and prefix resolution
//==============================================================

UnitExpression
UnitParser::resolveUnit(
    const Token& token) const
{
    const std::string& symbol =
        token.text();

    //----------------------------------------------------------
    // Complete unit symbols always take precedence.
    //
    // This is essential for "kg": it must be resolved directly,
    // rather than interpreted as kilo + gram.
    //----------------------------------------------------------

    if (unitRegistry_.contains(symbol))
    {
        return UnitExpression(
            unitRegistry_.find(symbol));
    }

    //----------------------------------------------------------
    // Try SI prefix + base unit.
    //
    // Prefixes are attempted from longest to shortest so that
    // "da" is checked before "d".
    //----------------------------------------------------------

    const std::vector<std::string> prefixes =
        prefixRegistry_.symbolsByDescendingLength();

    bool foundForbiddenBase = false;
    std::string forbiddenBase;

    for (const std::string& prefixSymbol : prefixes)
    {
        if (symbol.size() <= prefixSymbol.size())
        {
            continue;
        }

        if (symbol.compare(
                0,
                prefixSymbol.size(),
                prefixSymbol) != 0)
        {
            continue;
        }

        const std::string baseSymbol =
            symbol.substr(prefixSymbol.size());

        if (!unitRegistry_.contains(baseSymbol))
        {
            continue;
        }

        const Unit& baseUnit =
            unitRegistry_.find(baseSymbol);

        if (!baseUnit.allowsPrefix())
        {
            foundForbiddenBase = true;
            forbiddenBase = baseSymbol;
            continue;
        }

        const Prefix& prefix =
            prefixRegistry_.find(prefixSymbol);

        const double combinedScale =
            prefix.factor()
            * baseUnit.scale();

        return UnitExpression(
            baseUnit.dimension(),
            combinedScale);
    }

    if (foundForbiddenBase)
    {
        throwError(
            token,
            "SI prefix cannot be applied to unit '"
            + forbiddenBase
            + "'");
    }

    throwError(
        token,
        "unknown unit '"
        + symbol
        + "'");
}

//==============================================================
// Exponents
//==============================================================

common::Rational
UnitParser::parseExponent(
    TokenStream& stream) const
{
    //----------------------------------------------------------
    // Bracketed exponent
    //
    // Rational exponents must be enclosed in brackets:
    //
    //   m^(1/2)
    //   m^[-1/2]
    //   m^{-3/2}
    //----------------------------------------------------------

    if (isOpeningBracket(stream.current()))
    {
        const Token openingToken =
            stream.current();

        const char opening =
            openingToken.text()[0];

        const char closing =
            matchingClosingBracket(opening);

        stream.consume();

        const common::Rational result =
            parseRational(stream);

        if (!stream.check(
                std::string(1, closing)))
        {
            if (stream.atEnd())
            {
                throwError(
                    openingToken,
                    std::string(
                        "missing closing bracket '")
                    + closing
                    + "' in exponent");
            }

            throwError(
                stream.current(),
                std::string(
                    "expected closing bracket '")
                + closing
                + "' in exponent");
        }

        stream.consume();

        return result;
    }

    //----------------------------------------------------------
    // Unbracketed exponent
    //
    // Only a single integer is accepted here:
    //
    //   m^2
    //   m^-2
    //
    // Requiring brackets around rational exponents avoids
    // ambiguity between:
    //
    //   m^2/s      unit division
    //
    // and:
    //
    //   m^(2/3)    rational exponent
    //----------------------------------------------------------

    if (!stream.check(TokenType::Number))
    {
        throwError(
            stream.current(),
            "expected an integer exponent");
    }

    const int exponent =
        parseInteger(stream.current());

    stream.consume();

    return common::Rational(exponent);
}

common::Rational
UnitParser::parseRational(
    TokenStream& stream) const
{
    if (!stream.check(TokenType::Number))
    {
        throwError(
            stream.current(),
            "expected an integer exponent");
    }

    const int numerator =
        parseInteger(stream.current());

    stream.consume();

    if (!stream.match("/"))
    {
        return common::Rational(numerator);
    }

    if (!stream.check(TokenType::Number))
    {
        throwError(
            stream.current(),
            "expected an integer denominator");
    }

    const int denominator =
        parseInteger(stream.current());

    stream.consume();

    if (denominator == 0)
    {
        throwError(
            stream.current(),
            "unit exponent denominator cannot be zero");
    }

    return common::Rational(
        numerator,
        denominator);
}

int
UnitParser::parseInteger(
    const Token& token) const
{
    if (!token.is(TokenType::Number))
    {
        throwError(
            token,
            "expected an integer");
    }

    const std::string& text =
        token.text();

    //----------------------------------------------------------
    // Rational exponents must use integer numerator and
    // denominator components. Decimal and scientific notation
    // are rejected here.
    //----------------------------------------------------------

    if (text.find('.') != std::string::npos
        || text.find('e') != std::string::npos
        || text.find('E') != std::string::npos)
    {
        throwError(
            token,
            "unit exponent must be an integer or rational number");
    }

    std::size_t processed = 0;
    long long value = 0;

    try
    {
        value =
            std::stoll(
                text,
                &processed,
                10);
    }
    catch (const std::exception&)
    {
        throwError(
            token,
            "invalid integer exponent '"
            + text
            + "'");
    }

    if (processed != text.size())
    {
        throwError(
            token,
            "invalid integer exponent '"
            + text
            + "'");
    }

    if (value
            < static_cast<long long>(
                std::numeric_limits<int>::min())
        || value
            > static_cast<long long>(
                std::numeric_limits<int>::max()))
    {
        throwError(
            token,
            "integer exponent is outside the supported range");
    }

    return static_cast<int>(value);
}

//==============================================================
// Grammar helpers
//==============================================================

bool
UnitParser::beginsFactor(
    const Token& token) const
{
    return token.is(TokenType::Identifier)
        || isOpeningBracket(token);
}

bool
UnitParser::isOpeningBracket(
    const Token& token) const
{
    return token.is("(")
        || token.is("[")
        || token.is("{");
}

bool
UnitParser::isClosingBracket(
    const Token& token) const
{
    return token.is(")")
        || token.is("]")
        || token.is("}");
}

char
UnitParser::matchingClosingBracket(
    char opening) const
{
    switch (opening)
    {
    case '(':
        return ')';

    case '[':
        return ']';

    case '{':
        return '}';

    default:
        throw std::invalid_argument(
            "UnitParser: invalid opening bracket.");
    }
}

int
UnitParser::bracketLevel(
    char bracket) const
{
    switch (bracket)
    {
    case '(':
        return 1;

    case '[':
        return 2;

    case '{':
        return 3;

    default:
        return 0;
    }
}

void
UnitParser::validateBracketNesting(
    char opening,
    const std::vector<char>& bracketStack,
    const Token& token) const
{
    if (bracketStack.empty())
    {
        return;
    }

    const char outer =
        bracketStack.back();

    //----------------------------------------------------------
    // Same-shaped brackets may be nested.
    //
    //   ((m))
    //   [[m]]
    //   {{m}}
    //----------------------------------------------------------

    if (opening == outer)
    {
        return;
    }

    //----------------------------------------------------------
    // Mixed brackets may only become smaller inward:
    //
    //   { [ ( ) ] }
    //
    // Therefore:
    //
    //   [()]  is valid
    //   {[]}  is valid
    //   {()}  is valid
    //
    // but:
    //
    //   ([])  is invalid
    //   [({})] is invalid
    //----------------------------------------------------------

    if (bracketLevel(opening)
        < bracketLevel(outer))
    {
        return;
    }

    throwError(
        token,
        std::string("invalid bracket nesting: '")
        + outer
        + "' cannot contain '"
        + opening
        + "'");
}

//==============================================================
// Diagnostics
//==============================================================

void
UnitParser::throwError(
    const Token& token,
    const std::string& message) const
{
    std::ostringstream oss;

    oss << "UnitParser: "
        << message
        << " at line "
        << token.line()
        << ", column "
        << token.column()
        << ".";

    throw std::runtime_error(
        oss.str());
}

} // namespace ntic::lbm::config
