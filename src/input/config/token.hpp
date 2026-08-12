#pragma once

#include <cstddef>
#include <string>

namespace ntic::lbm::config
{

enum class TokenType
{
    Identifier,
    Number,
    String,
    Symbol,
    End
};

class Token
{
public:
    Token();

    Token(
        TokenType type,
        const std::string& text,
        std::size_t line,
        std::size_t column);

    TokenType type() const;
    const std::string& text() const;
    std::size_t line() const;
    std::size_t column() const;

    bool is(TokenType type) const;
    bool is(const std::string& text) const;

    std::string toString() const;

private:
    TokenType type_;
    std::string text_;
    std::size_t line_;
    std::size_t column_;
};

} // namespace ntic::lbm::config

