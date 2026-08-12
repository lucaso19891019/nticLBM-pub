#include "token_stream.hpp"

#include <stdexcept>

namespace ntic::lbm::config
{

//==============================================================
// Constructor
//==============================================================

TokenStream::TokenStream(
    const std::vector<Token>& tokens)
    : tokens_(tokens),
      position_(0)
{
    if (tokens_.empty())
    {
        throw std::invalid_argument(
            "TokenStream: token list cannot be empty.");
    }

    if (!tokens_.back().is(TokenType::End))
    {
        throw std::invalid_argument(
            "TokenStream: token list must end with an End token.");
    }
}

//==============================================================
// Access
//==============================================================

const Token&
TokenStream::current() const
{
    if (position_ >= tokens_.size())
    {
        return endToken();
    }

    return tokens_[position_];
}

const Token&
TokenStream::peek(
    std::size_t offset) const
{
    if (offset > tokens_.size()
        || position_ >= tokens_.size()
        || offset >= tokens_.size() - position_)
    {
        return endToken();
    }

    return tokens_[position_ + offset];
}

//==============================================================
// Navigation
//==============================================================

void
TokenStream::consume()
{
    if (!atEnd())
    {
        ++position_;
    }
}

bool
TokenStream::match(
    TokenType type)
{
    if (!check(type))
    {
        return false;
    }

    consume();

    return true;
}

bool
TokenStream::match(
    const std::string& text)
{
    if (!check(text))
    {
        return false;
    }

    consume();

    return true;
}

//==============================================================
// Query
//==============================================================

bool
TokenStream::check(
    TokenType type) const
{
    return current().is(type);
}

bool
TokenStream::check(
    const std::string& text) const
{
    return current().is(text);
}

bool
TokenStream::atEnd() const
{
    return current().is(TokenType::End);
}

std::size_t
TokenStream::position() const
{
    return position_;
}

//==============================================================
// Utility
//==============================================================

void
TokenStream::reset()
{
    position_ = 0;
}

//==============================================================
// Private
//==============================================================

const Token&
TokenStream::endToken() const
{
    return tokens_.back();
}

} // namespace ntic::lbm::config
