#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "token.hpp"

namespace ntic::lbm::config
{

class TokenStream
{
public:
    //----------------------------------------------------------
    // Constructor
    //----------------------------------------------------------

    explicit TokenStream(
        const std::vector<Token>& tokens);

    //----------------------------------------------------------
    // Access
    //----------------------------------------------------------

    const Token& current() const;

    const Token& peek(
        std::size_t offset = 1) const;

    //----------------------------------------------------------
    // Navigation
    //----------------------------------------------------------

    void consume();

    bool match(
        TokenType type);

    bool match(
        const std::string& text);

    //----------------------------------------------------------
    // Query
    //----------------------------------------------------------

    bool check(
        TokenType type) const;

    bool check(
        const std::string& text) const;

    bool atEnd() const;

    std::size_t position() const;

    //----------------------------------------------------------
    // Utility
    //----------------------------------------------------------

    void reset();

private:
    const Token& endToken() const;

private:
    std::vector<Token> tokens_;

    std::size_t position_;
};

} // namespace ntic::lbm::config
