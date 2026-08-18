#include "lexer.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>

using ntic::lbm::config::Lexer;
using ntic::lbm::config::Token;
using ntic::lbm::config::TokenType;

namespace
{

void
expectToken(
    const Token& token,
    TokenType type,
    const std::string& text,
    std::size_t line,
    std::size_t column)
{
    NTIC_EXPECT_TRUE(token.type() == type);
    NTIC_EXPECT_TRUE(token.text() == text);
    NTIC_EXPECT_TRUE(token.line() == line);
    NTIC_EXPECT_TRUE(token.column() == column);
}

} // namespace

int main()
{
    const Lexer lexer;

    //----------------------------------------------------------
    // Empty input
    //----------------------------------------------------------

    {
        const auto tokens = lexer.tokenize("");

        NTIC_EXPECT_TRUE(tokens.size() == 1);

        expectToken(
            tokens[0],
            TokenType::End,
            "",
            1,
            1);
    }

    //----------------------------------------------------------
    // Identifiers
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "density kg Pa atm abc123 solver_type _private");

        NTIC_EXPECT_TRUE(tokens.size() == 8);

        expectToken(tokens[0], TokenType::Identifier, "density",     1, 1);
        expectToken(tokens[1], TokenType::Identifier, "kg",          1, 9);
        expectToken(tokens[2], TokenType::Identifier, "Pa",          1, 12);
        expectToken(tokens[3], TokenType::Identifier, "atm",         1, 15);
        expectToken(tokens[4], TokenType::Identifier, "abc123",      1, 19);
        expectToken(tokens[5], TokenType::Identifier, "solver_type", 1, 26);
        expectToken(tokens[6], TokenType::Identifier, "_private",    1, 38);
        expectToken(tokens[7], TokenType::End,        "",            1, 46);
    }

    //----------------------------------------------------------
    // Numbers
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "0 12 -3 +4 3.14 .5 10. -0.25 +.75");

        NTIC_EXPECT_TRUE(tokens.size() == 10);

        expectToken(tokens[0], TokenType::Number, "0",     1, 1);
        expectToken(tokens[1], TokenType::Number, "12",    1, 3);
        expectToken(tokens[2], TokenType::Number, "-3",    1, 6);
        expectToken(tokens[3], TokenType::Number, "+4",    1, 9);
        expectToken(tokens[4], TokenType::Number, "3.14",  1, 12);
        expectToken(tokens[5], TokenType::Number, ".5",    1, 17);
        expectToken(tokens[6], TokenType::Number, "10.",   1, 20);
        expectToken(tokens[7], TokenType::Number, "-0.25", 1, 24);
        expectToken(tokens[8], TokenType::Number, "+.75",  1, 30);
        expectToken(tokens[9], TokenType::End,    "",      1, 34);
    }

    {
        const auto tokens =
            lexer.tokenize(
                "1e-6 2E8 -3.5e+10 +4E-2");

        NTIC_EXPECT_TRUE(tokens.size() == 5);

        expectToken(tokens[0], TokenType::Number, "1e-6",     1, 1);
        expectToken(tokens[1], TokenType::Number, "2E8",      1, 6);
        expectToken(tokens[2], TokenType::Number, "-3.5e+10", 1, 10);
        expectToken(tokens[3], TokenType::Number, "+4E-2",    1, 19);
        expectToken(tokens[4], TokenType::End,    "",         1, 24);
    }

    //----------------------------------------------------------
    // Symbols and brackets
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "*/^=+-()[]{} ,:");

        NTIC_EXPECT_TRUE(tokens.size() == 15);

        expectToken(tokens[0],  TokenType::Symbol, "*", 1, 1);
        expectToken(tokens[1],  TokenType::Symbol, "/", 1, 2);
        expectToken(tokens[2],  TokenType::Symbol, "^", 1, 3);
        expectToken(tokens[3],  TokenType::Symbol, "=", 1, 4);
        expectToken(tokens[4],  TokenType::Symbol, "+", 1, 5);
        expectToken(tokens[5],  TokenType::Symbol, "-", 1, 6);
        expectToken(tokens[6],  TokenType::Symbol, "(", 1, 7);
        expectToken(tokens[7],  TokenType::Symbol, ")", 1, 8);
        expectToken(tokens[8],  TokenType::Symbol, "[", 1, 9);
        expectToken(tokens[9],  TokenType::Symbol, "]", 1, 10);
        expectToken(tokens[10], TokenType::Symbol, "{", 1, 11);
        expectToken(tokens[11], TokenType::Symbol, "}", 1, 12);
        expectToken(tokens[12], TokenType::Symbol, ",", 1, 14);
        expectToken(tokens[13], TokenType::Symbol, ":", 1, 15);
        expectToken(tokens[14], TokenType::End,    "",  1, 16);
    }

    //----------------------------------------------------------
    // Unit expression
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize("kg/(m*s^2)");

        NTIC_EXPECT_TRUE(tokens.size() == 10);

        expectToken(tokens[0], TokenType::Identifier, "kg", 1, 1);
        expectToken(tokens[1], TokenType::Symbol,     "/",  1, 3);
        expectToken(tokens[2], TokenType::Symbol,     "(",  1, 4);
        expectToken(tokens[3], TokenType::Identifier, "m",  1, 5);
        expectToken(tokens[4], TokenType::Symbol,     "*",  1, 6);
        expectToken(tokens[5], TokenType::Identifier, "s",  1, 7);
        expectToken(tokens[6], TokenType::Symbol,     "^",  1, 8);
        expectToken(tokens[7], TokenType::Number,     "2",  1, 9);
        expectToken(tokens[8], TokenType::Symbol,     ")",  1, 10);
        expectToken(tokens[9], TokenType::End,        "",   1, 11);
    }

    {
        const auto tokens =
            lexer.tokenize("{[kg/(m*s)]}");

        NTIC_EXPECT_TRUE(tokens.size() == 12);

        expectToken(tokens[0],  TokenType::Symbol,     "{",  1, 1);
        expectToken(tokens[1],  TokenType::Symbol,     "[",  1, 2);
        expectToken(tokens[2],  TokenType::Identifier, "kg", 1, 3);
        expectToken(tokens[3],  TokenType::Symbol,     "/",  1, 5);
        expectToken(tokens[4],  TokenType::Symbol,     "(",  1, 6);
        expectToken(tokens[5],  TokenType::Identifier, "m",  1, 7);
        expectToken(tokens[6],  TokenType::Symbol,     "*",  1, 8);
        expectToken(tokens[7],  TokenType::Identifier, "s",  1, 9);
        expectToken(tokens[8],  TokenType::Symbol,     ")",  1, 10);
        expectToken(tokens[9],  TokenType::Symbol,     "]",  1, 11);
        expectToken(tokens[10], TokenType::Symbol,     "}",  1, 12);
        expectToken(tokens[11], TokenType::End,        "",   1, 13);
    }

    //----------------------------------------------------------
    // Strings
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "stl_path=\"/home/user/model.stl\"");

        NTIC_EXPECT_TRUE(tokens.size() == 4);

        expectToken(tokens[0], TokenType::Identifier, "stl_path",             1, 1);
        expectToken(tokens[1], TokenType::Symbol,     "=",                    1, 9);
        expectToken(tokens[2], TokenType::String,     "/home/user/model.stl", 1, 10);
        expectToken(tokens[3], TokenType::End,        "",                     1, 32);
    }

    {
        const auto tokens =
            lexer.tokenize(
                "stl_path='part.stl'");

        NTIC_EXPECT_TRUE(tokens.size() == 4);

        expectToken(tokens[2], TokenType::String, "part.stl", 1, 10);
    }

    //----------------------------------------------------------
    // Comments and hashes inside strings
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "stl_path=\"/home/user/#001/model.stl\" # comment");

        NTIC_EXPECT_TRUE(tokens.size() == 4);

        expectToken(
            tokens[2],
            TokenType::String,
            "/home/user/#001/model.stl",
            1,
            10);

        NTIC_EXPECT_TRUE(tokens.back().is(TokenType::End));
    }

    {
        const auto tokens =
            lexer.tokenize(
                "name='part#2.stl'#comment");

        NTIC_EXPECT_TRUE(tokens.size() == 4);

        expectToken(
            tokens[2],
            TokenType::String,
            "part#2.stl",
            1,
            6);
    }

    {
        const auto tokens =
            lexer.tokenize("# comment only");

        NTIC_EXPECT_TRUE(tokens.size() == 1);
        NTIC_EXPECT_TRUE(tokens[0].is(TokenType::End));
    }

    {
        const auto tokens =
            lexer.tokenize("density=998#water");

        NTIC_EXPECT_TRUE(tokens.size() == 4);

        expectToken(tokens[0], TokenType::Identifier, "density", 1, 1);
        expectToken(tokens[1], TokenType::Symbol,     "=",       1, 8);
        expectToken(tokens[2], TokenType::Number,     "998",     1, 9);
        NTIC_EXPECT_TRUE(tokens[3].is(TokenType::End));
    }

    //----------------------------------------------------------
    // Multi-line input and positions
    //----------------------------------------------------------

    {
        const std::string input =
            "# fluid properties\n"
            "density = 998 kg/m^3 # water\n"
            "nu = 1e-6 m^2/s\n"
            "stl_path = '/tmp/part#1.stl'\n";

        const auto tokens =
            lexer.tokenize(input);

        NTIC_EXPECT_TRUE(tokens.size() == 20);

        expectToken(tokens[0],  TokenType::Identifier, "density",         2, 1);
        expectToken(tokens[1],  TokenType::Symbol,     "=",               2, 9);
        expectToken(tokens[2],  TokenType::Number,     "998",             2, 11);
        expectToken(tokens[3],  TokenType::Identifier, "kg",              2, 15);
        expectToken(tokens[4],  TokenType::Symbol,     "/",               2, 17);
        expectToken(tokens[5],  TokenType::Identifier, "m",               2, 18);
        expectToken(tokens[6],  TokenType::Symbol,     "^",               2, 19);
        expectToken(tokens[7],  TokenType::Number,     "3",               2, 20);
        expectToken(tokens[8],  TokenType::Identifier, "nu",              3, 1);
        expectToken(tokens[9],  TokenType::Symbol,     "=",               3, 4);
        expectToken(tokens[10], TokenType::Number,     "1e-6",            3, 6);
        expectToken(tokens[11], TokenType::Identifier, "m",               3, 11);
        expectToken(tokens[12], TokenType::Symbol,     "^",               3, 12);
        expectToken(tokens[13], TokenType::Number,     "2",               3, 13);
        expectToken(tokens[14], TokenType::Symbol,     "/",               3, 14);
        expectToken(tokens[15], TokenType::Identifier, "s",               3, 15);
        expectToken(tokens[16], TokenType::Identifier, "stl_path",        4, 1);
        expectToken(tokens[17], TokenType::Symbol,     "=",               4, 10);
        expectToken(tokens[18], TokenType::String,     "/tmp/part#1.stl", 4, 12);
        expectToken(tokens[19], TokenType::End,        "",                5, 1);
    }

    //----------------------------------------------------------
    // Token helpers
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize("density=1");

        NTIC_EXPECT_TRUE(tokens[0].is(TokenType::Identifier));
        NTIC_EXPECT_TRUE(tokens[0].is("density"));
        NTIC_EXPECT_TRUE(tokens[1].is(TokenType::Symbol));
        NTIC_EXPECT_TRUE(tokens[1].is("="));

        const std::string text =
            tokens[0].toString();

        NTIC_EXPECT_TRUE(
            text.find("Identifier")
            != std::string::npos);

        NTIC_EXPECT_TRUE(
            text.find("density")
            != std::string::npos);

        NTIC_EXPECT_TRUE(
            text.find("1:1")
            != std::string::npos);
    }

    //----------------------------------------------------------
    // Invalid input
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        lexer.tokenize("@"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        lexer.tokenize("density = 1 $"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        lexer.tokenize("\"abc"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        lexer.tokenize("'abc"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        lexer.tokenize("\"abc\nnext"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        lexer.tokenize("1e"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        lexer.tokenize("1e+"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        lexer.tokenize("-2E-"),
        std::runtime_error);

    //----------------------------------------------------------

    ntic::lbm::test::printSummary();

    return 0;
}

