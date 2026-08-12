
#include "quantity_parser.hpp"

#include <sstream>
#include <stdexcept>

namespace ntic::lbm::config
{

QuantityParser::QuantityParser(
    const UnitRegistry& unitRegistry,
    const PrefixRegistry& prefixRegistry)
    : lexer_(),
      unitParser_(unitRegistry, prefixRegistry)
{
}

Quantity
QuantityParser::parse(
    const std::string& text) const
{
    const std::vector<Token> tokens =
        lexer_.tokenize(text);

    TokenStream stream(tokens);

    if (stream.atEnd())
    {
        throw std::runtime_error(
            "QuantityParser: empty quantity.");
    }

    if (!stream.check(TokenType::Number))
    {
        throw std::runtime_error(
            "QuantityParser: quantity must begin with a number.");
    }

    const double value =
        parseNumber(stream.current());

    stream.consume();

    if (stream.atEnd())
    {
        return Quantity::Dimensionless(value);
    }

    std::string unitText;

    while (!stream.atEnd())
    {
        unitText += stream.current().text();
        stream.consume();
    }

    const UnitExpression unit =
        unitParser_.parse(unitText);

    return Quantity(value, unit);
}

double
QuantityParser::parseNumber(
    const Token& token) const
{
    if (!token.is(TokenType::Number))
    {
        throw std::runtime_error(
            "QuantityParser: expected a numeric literal.");
    }

    std::size_t processed = 0;

    double value = 0.0;

    try
    {
        value =
            std::stod(
                token.text(),
                &processed);
    }
    catch (const std::exception&)
    {
        throw std::runtime_error(
            "QuantityParser: invalid numeric literal \""
            + token.text()
            + "\".");
    }

    if (processed != token.text().size())
    {
        throw std::runtime_error(
            "QuantityParser: invalid numeric literal \""
            + token.text()
            + "\".");
    }

    return value;
}

} // namespace ntic::lbm::config

