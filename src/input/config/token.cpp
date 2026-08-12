#include "token.hpp"

#include <sstream>

namespace ntic::lbm::config
{

Token::Token()
    : type_(TokenType::End),
      text_(),
      line_(0),
      column_(0)
{
}

Token::Token(
    TokenType type,
    const std::string& text,
    std::size_t line,
    std::size_t column)
    : type_(type),
      text_(text),
      line_(line),
      column_(column)
{
}

TokenType
Token::type() const
{
    return type_;
}

const std::string&
Token::text() const
{
    return text_;
}

std::size_t
Token::line() const
{
    return line_;
}

std::size_t
Token::column() const
{
    return column_;
}

bool
Token::is(TokenType type) const
{
    return type_ == type;
}

bool
Token::is(const std::string& text) const
{
    return text_ == text;
}

std::string
Token::toString() const
{
    std::ostringstream oss;

    oss << "[";

    switch (type_)
    {
    case TokenType::Identifier:
        oss << "Identifier";
        break;
    case TokenType::Number:
        oss << "Number";
        break;
    case TokenType::String:
        oss << "String";
        break;
    case TokenType::Symbol:
        oss << "Symbol";
        break;
    case TokenType::End:
        oss << "End";
        break;
    }

    oss << "] \""
        << text_
        << "\" ("
        << line_
        << ":"
        << column_
        << ")";

    return oss.str();
}

} // namespace ntic::lbm::config

