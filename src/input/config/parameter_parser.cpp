
#include "parameter_parser.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ntic::lbm::config
{

ParameterParser::ParameterParser(
    const UnitRegistry& unitRegistry,
    const PrefixRegistry& prefixRegistry)
    :
    lexer_(),
    quantityParser_(
        unitRegistry,
        prefixRegistry)
{
}

ParsedParameter
ParameterParser::parse(
    const std::string& text) const
{
    const std::vector<Token> tokens =
        lexer_.tokenize(text);

    TokenStream stream(tokens);

    if (stream.atEnd())
    {
        throwError(
            Token(),
            "empty parameter assignment");
    }

    if (!stream.check(TokenType::Identifier))
    {
        throwError(
            stream.current(),
            "parameter assignment must begin with a parameter name");
    }

    const std::string parameterName =
        stream.current().text();

    stream.consume();

    if (stream.check("=") ||
        stream.check(":"))
    {
        stream.consume();
    }

    if (stream.atEnd())
    {
        throwError(
            Token(
                TokenType::Identifier,
                parameterName,
                1,
                1),
            "parameter '" + parameterName + "' has no value");
    }

    // quoted string
    if (stream.check(TokenType::String))
    {
        const std::string value =
            stream.current().text();

        stream.consume();

        if (!stream.atEnd())
        {
            throwError(
                stream.current(),
                "unexpected content after string");
        }

        return ParsedParameter{
            parameterName,
            ParameterValue(value)};
    }

    // bool
    if (stream.check(TokenType::Identifier))
    {
        bool b=false;

        if (tryParseBool(
                stream.current(),
                b))
        {
            stream.consume();

            if (!stream.atEnd())
            {
                throwError(
                    stream.current(),
                    "unexpected content after bool");
            }

            return ParsedParameter{
                parameterName,
                ParameterValue(b)};
        }
    }

    // quantity
    if (stream.check(TokenType::Number))
    {
        const Quantity quantity =
            quantityParser_.parse(
                rebuildRemainingText(stream));

        return ParsedParameter{
            parameterName,
            ParameterValue(quantity)};
    }

    // identifier expression -> unquoted string
    //
    // Unquoted string values may contain multiple identifiers
    // separated by whitespace and/or hyphens.
    //
    // Examples:
    //
    //     BGK
    //     central moment
    //     ZH-vel
    //     Zou-He velocity
    //
    // Other symbols are rejected. More general strings should
    // still be quoted.

    if (stream.check(TokenType::Identifier))
    {
        std::string value;

        bool previousWasIdentifier = false;

        while (!stream.atEnd())
        {
            if (stream.check(TokenType::Identifier))
            {
                if (previousWasIdentifier)
                {
                    value += ' ';
                }

                value +=
                    stream.current().text();

                stream.consume();

                previousWasIdentifier = true;

                continue;
            }

            if (stream.check("-"))
            {
                value += '-';

                stream.consume();

                previousWasIdentifier = false;

                continue;
            }

            throwError(
                stream.current(),
                "unsupported content in unquoted string value");
        }

        return ParsedParameter{
            parameterName,
            ParameterValue(value)};
    }

    throwError(
        stream.current(),
        "unsupported parameter value");
}

bool
ParameterParser::tryParseBool(
    const Token& token,
    bool& value) const
{
    if (!token.is(TokenType::Identifier))
    {
        return false;
    }

    std::string text =
        token.text();

    std::transform(
        text.begin(),
        text.end(),
        text.begin(),
        [](unsigned char c)
        {
            return static_cast<char>(std::tolower(c));
        });

    if (text=="true" ||
        text=="on" ||
        text=="yes")
    {
        value=true;
        return true;
    }

    if (text=="false" ||
        text=="off" ||
        text=="no")
    {
        value=false;
        return true;
    }

    return false;
}

std::string
ParameterParser::rebuildRemainingText(
    TokenStream& stream) const
{
    std::string text;

    while (!stream.atEnd())
    {
        text +=
            stream.current().text();

        stream.consume();
    }

    return text;
}

void
ParameterParser::throwError(
    const Token& token,
    const std::string& message) const
{
    std::ostringstream oss;

    oss
        << "ParameterParser: "
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

