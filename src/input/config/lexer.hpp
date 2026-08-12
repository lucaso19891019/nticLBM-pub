#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "token.hpp"

namespace ntic::lbm::config
{

class Lexer
{
public:
    std::vector<Token>
    tokenize(const std::string& text) const;

private:
    void advance(
        char c,
        std::size_t& position,
        std::size_t& line,
        std::size_t& column) const;

    void skipWhitespace(
        const std::string& text,
        std::size_t& position,
        std::size_t& line,
        std::size_t& column) const;

    void skipComment(
        const std::string& text,
        std::size_t& position,
        std::size_t& line,
        std::size_t& column) const;

    Token readIdentifier(
        const std::string& text,
        std::size_t& position,
        std::size_t line,
        std::size_t& column) const;

    Token readNumber(
        const std::string& text,
        std::size_t& position,
        std::size_t line,
        std::size_t& column) const;

    Token readString(
        const std::string& text,
        std::size_t& position,
        std::size_t& line,
        std::size_t& column) const;

    Token readSymbol(
        const std::string& text,
        std::size_t& position,
        std::size_t line,
        std::size_t& column) const;

    bool isIdentifierStart(char c) const;
    bool isIdentifierCharacter(char c) const;

    bool isNumberStart(
        const std::string& text,
        std::size_t position) const;

    bool isSymbol(char c) const;
};

} // namespace ntic::lbm::config

