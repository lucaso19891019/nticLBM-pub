#include "lexer.hpp"

#include <cctype>
#include <sstream>
#include <stdexcept>

namespace ntic::lbm::config
{

std::vector<Token>
Lexer::tokenize(const std::string& text) const
{
    std::vector<Token> tokens;

    std::size_t position = 0;
    std::size_t line = 1;
    std::size_t column = 1;

    while (position < text.size())
    {
        skipWhitespace(text, position, line, column);

        if (position >= text.size())
            break;

        const char current = text[position];

        if (current == '#')
        {
            skipComment(text, position, line, column);
            continue;
        }

        if (current == '"' || current == '\'')
        {
            tokens.push_back(
                readString(text, position, line, column));
            continue;
        }

        if (isIdentifierStart(current))
        {
            tokens.push_back(
                readIdentifier(text, position, line, column));
            continue;
        }

        if (isNumberStart(text, position))
        {
            tokens.push_back(
                readNumber(text, position, line, column));
            continue;
        }

        if (isSymbol(current))
        {
            tokens.push_back(
                readSymbol(text, position, line, column));
            continue;
        }

        std::ostringstream oss;
        oss << "Lexer: unexpected character '"
            << current
            << "' at line "
            << line
            << ", column "
            << column
            << ".";

        throw std::runtime_error(oss.str());
    }

    tokens.emplace_back(
        TokenType::End,
        "",
        line,
        column);

    return tokens;
}

void
Lexer::advance(
    char c,
    std::size_t& position,
    std::size_t& line,
    std::size_t& column) const
{
    ++position;

    if (c == '\n')
    {
        ++line;
        column = 1;
    }
    else
    {
        ++column;
    }
}

void
Lexer::skipWhitespace(
    const std::string& text,
    std::size_t& position,
    std::size_t& line,
    std::size_t& column) const
{
    while (position < text.size())
    {
        const unsigned char c =
            static_cast<unsigned char>(text[position]);

        if (!std::isspace(c))
            break;

        advance(
            text[position],
            position,
            line,
            column);
    }
}

void
Lexer::skipComment(
    const std::string& text,
    std::size_t& position,
    std::size_t& line,
    std::size_t& column) const
{
    while (position < text.size()
           && text[position] != '\n')
    {
        advance(
            text[position],
            position,
            line,
            column);
    }
}

Token
Lexer::readIdentifier(
    const std::string& text,
    std::size_t& position,
    std::size_t line,
    std::size_t& column) const
{
    const std::size_t startPosition = position;
    const std::size_t startColumn = column;

    while (position < text.size()
           && isIdentifierCharacter(text[position]))
    {
        advance(
            text[position],
            position,
            line,
            column);
    }

    return Token(
        TokenType::Identifier,
        text.substr(
            startPosition,
            position - startPosition),
        line,
        startColumn);
}

Token
Lexer::readNumber(
    const std::string& text,
    std::size_t& position,
    std::size_t line,
    std::size_t& column) const
{
    const std::size_t startPosition = position;
    const std::size_t startColumn = column;

    if (position < text.size()
        && (text[position] == '+'
            || text[position] == '-'))
    {
        advance(
            text[position],
            position,
            line,
            column);
    }

    bool hasDigits = false;

    while (position < text.size()
           && std::isdigit(
               static_cast<unsigned char>(
                   text[position])))
    {
        hasDigits = true;
        advance(
            text[position],
            position,
            line,
            column);
    }

    if (position < text.size()
        && text[position] == '.')
    {
        advance(
            text[position],
            position,
            line,
            column);

        while (position < text.size()
               && std::isdigit(
                   static_cast<unsigned char>(
                       text[position])))
        {
            hasDigits = true;
            advance(
                text[position],
                position,
                line,
                column);
        }
    }

    if (!hasDigits)
    {
        std::ostringstream oss;
        oss << "Lexer: invalid number at line "
            << line
            << ", column "
            << startColumn
            << ".";
        throw std::runtime_error(oss.str());
    }

    if (position < text.size()
        && (text[position] == 'e'
            || text[position] == 'E'))
    {
        advance(
            text[position],
            position,
            line,
            column);

        if (position < text.size()
            && (text[position] == '+'
                || text[position] == '-'))
        {
            advance(
                text[position],
                position,
                line,
                column);
        }

        const std::size_t exponentStart =
            position;

        while (position < text.size()
               && std::isdigit(
                   static_cast<unsigned char>(
                       text[position])))
        {
            advance(
                text[position],
                position,
                line,
                column);
        }

        if (position == exponentStart)
        {
            std::ostringstream oss;
            oss << "Lexer: scientific exponent has no digits "
                << "at line "
                << line
                << ", column "
                << startColumn
                << ".";
            throw std::runtime_error(oss.str());
        }
    }

    return Token(
        TokenType::Number,
        text.substr(
            startPosition,
            position - startPosition),
        line,
        startColumn);
}

Token
Lexer::readString(
    const std::string& text,
    std::size_t& position,
    std::size_t& line,
    std::size_t& column) const
{
    const char quote = text[position];
    const std::size_t startLine = line;
    const std::size_t startColumn = column;

    advance(
        text[position],
        position,
        line,
        column);

    std::string value;

    while (position < text.size())
    {
        const char current = text[position];

        if (current == quote)
        {
            advance(
                current,
                position,
                line,
                column);

            return Token(
                TokenType::String,
                value,
                startLine,
                startColumn);
        }

        if (current == '\n'
            || current == '\r')
        {
            std::ostringstream oss;
            oss << "Lexer: unterminated string beginning at line "
                << startLine
                << ", column "
                << startColumn
                << ".";
            throw std::runtime_error(oss.str());
        }

        value.push_back(current);

        advance(
            current,
            position,
            line,
            column);
    }

    std::ostringstream oss;
    oss << "Lexer: unterminated string beginning at line "
        << startLine
        << ", column "
        << startColumn
        << ".";
    throw std::runtime_error(oss.str());
}

Token
Lexer::readSymbol(
    const std::string& text,
    std::size_t& position,
    std::size_t line,
    std::size_t& column) const
{
    const std::size_t startColumn = column;
    const std::string symbol(1, text[position]);

    advance(
        text[position],
        position,
        line,
        column);

    return Token(
        TokenType::Symbol,
        symbol,
        line,
        startColumn);
}

bool
Lexer::isIdentifierStart(char c) const
{
    const unsigned char value =
        static_cast<unsigned char>(c);

    return std::isalpha(value)
        || c == '_';
}

bool
Lexer::isIdentifierCharacter(char c) const
{
    const unsigned char value =
        static_cast<unsigned char>(c);

    return std::isalnum(value)
        || c == '_';
}

bool
Lexer::isNumberStart(
    const std::string& text,
    std::size_t position) const
{
    if (position >= text.size())
        return false;

    const char current = text[position];

    if (std::isdigit(
            static_cast<unsigned char>(current)))
    {
        return true;
    }

    if (current == '.')
    {
        return position + 1 < text.size()
            && std::isdigit(
                static_cast<unsigned char>(
                    text[position + 1]));
    }

    if (current == '+'
        || current == '-')
    {
        if (position + 1 >= text.size())
            return false;

        const char next = text[position + 1];

        if (std::isdigit(
                static_cast<unsigned char>(next)))
        {
            return true;
        }

        return next == '.'
            && position + 2 < text.size()
            && std::isdigit(
                static_cast<unsigned char>(
                    text[position + 2]));
    }

    return false;
}

bool
Lexer::isSymbol(char c) const
{
    switch (c)
    {
    case '*':
    case '/':
    case '^':
    case '=':
    case '+':
    case '-':
    case '(':
    case ')':
    case '[':
    case ']':
    case '{':
    case '}':
    case ',':
    case ':':
        return true;

    default:
        return false;
    }
}

} // namespace ntic::lbm::config

