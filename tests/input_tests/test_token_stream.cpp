#include "lexer.hpp"
#include "token_stream.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>
#include <vector>

using ntic::lbm::config::Lexer;
using ntic::lbm::config::Token;
using ntic::lbm::config::TokenStream;
using ntic::lbm::config::TokenType;

int main()
{
    const Lexer lexer;

    //----------------------------------------------------------
    // Constructor
    //----------------------------------------------------------

    {
        const std::vector<Token> tokens =
        {
            Token(
                TokenType::End,
                "",
                1,
                1)
        };

        TokenStream stream(tokens);

        NTIC_EXPECT_TRUE(
            stream.atEnd());

        NTIC_EXPECT_TRUE(
            stream.position() == 0);
    }

    //----------------------------------------------------------
    // Empty Token List
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        TokenStream(
            std::vector<Token>()),
        std::invalid_argument);

    //----------------------------------------------------------
    // Missing End Token
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        TokenStream(
            std::vector<Token>
            {
                Token(
                    TokenType::Identifier,
                    "kg",
                    1,
                    1)
            }),
        std::invalid_argument);

    //----------------------------------------------------------
    // current()
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "kg/m");

        TokenStream stream(tokens);

        NTIC_EXPECT_TRUE(
            stream.current().is(
                TokenType::Identifier));

        NTIC_EXPECT_TRUE(
            stream.current().text()
            == "kg");

        NTIC_EXPECT_TRUE(
            stream.position() == 0);
    }

    //----------------------------------------------------------
    // peek()
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "kg/m");

        TokenStream stream(tokens);

        NTIC_EXPECT_TRUE(
            stream.peek().text()
            == "/");

        NTIC_EXPECT_TRUE(
            stream.peek(2).text()
            == "m");

        NTIC_EXPECT_TRUE(
            stream.peek(3).is(
                TokenType::End));

        NTIC_EXPECT_TRUE(
            stream.peek(100).is(
                TokenType::End));

        NTIC_EXPECT_TRUE(
            stream.position() == 0);
    }

    //----------------------------------------------------------
    // consume()
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "kg/m");

        TokenStream stream(tokens);

        stream.consume();

        NTIC_EXPECT_TRUE(
            stream.current().text()
            == "/");

        NTIC_EXPECT_TRUE(
            stream.position() == 1);

        stream.consume();

        NTIC_EXPECT_TRUE(
            stream.current().text()
            == "m");

        NTIC_EXPECT_TRUE(
            stream.position() == 2);

        stream.consume();

        NTIC_EXPECT_TRUE(
            stream.atEnd());

        NTIC_EXPECT_TRUE(
            stream.position() == 3);
    }

    //----------------------------------------------------------
    // consume() At End
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize("");

        TokenStream stream(tokens);

        NTIC_EXPECT_TRUE(
            stream.atEnd());

        stream.consume();

        NTIC_EXPECT_TRUE(
            stream.atEnd());

        NTIC_EXPECT_TRUE(
            stream.position() == 0);
    }

    //----------------------------------------------------------
    // check(TokenType)
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "kg/m");

        TokenStream stream(tokens);

        NTIC_EXPECT_TRUE(
            stream.check(
                TokenType::Identifier));

        NTIC_EXPECT_TRUE(
            !stream.check(
                TokenType::Number));

        NTIC_EXPECT_TRUE(
            stream.position() == 0);
    }

    //----------------------------------------------------------
    // check(text)
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "kg/m");

        TokenStream stream(tokens);

        NTIC_EXPECT_TRUE(
            stream.check("kg"));

        NTIC_EXPECT_TRUE(
            !stream.check("m"));

        NTIC_EXPECT_TRUE(
            stream.position() == 0);
    }

    //----------------------------------------------------------
    // match(TokenType)
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "kg/m");

        TokenStream stream(tokens);

        NTIC_EXPECT_TRUE(
            stream.match(
                TokenType::Identifier));

        NTIC_EXPECT_TRUE(
            stream.current().text()
            == "/");

        NTIC_EXPECT_TRUE(
            stream.position() == 1);

        NTIC_EXPECT_TRUE(
            !stream.match(
                TokenType::Number));

        NTIC_EXPECT_TRUE(
            stream.position() == 1);
    }

    //----------------------------------------------------------
    // match(text)
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "kg/m");

        TokenStream stream(tokens);

        NTIC_EXPECT_TRUE(
            stream.match("kg"));

        NTIC_EXPECT_TRUE(
            stream.match("/"));

        NTIC_EXPECT_TRUE(
            stream.current().text()
            == "m");

        NTIC_EXPECT_TRUE(
            !stream.match("s"));

        NTIC_EXPECT_TRUE(
            stream.current().text()
            == "m");
    }

    //----------------------------------------------------------
    // Full Consumption
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "kg/(m*s^2)");

        TokenStream stream(tokens);

        NTIC_EXPECT_TRUE(
            stream.match("kg"));

        NTIC_EXPECT_TRUE(
            stream.match("/"));

        NTIC_EXPECT_TRUE(
            stream.match("("));

        NTIC_EXPECT_TRUE(
            stream.match("m"));

        NTIC_EXPECT_TRUE(
            stream.match("*"));

        NTIC_EXPECT_TRUE(
            stream.match("s"));

        NTIC_EXPECT_TRUE(
            stream.match("^"));

        NTIC_EXPECT_TRUE(
            stream.match(
                TokenType::Number));

        NTIC_EXPECT_TRUE(
            stream.match(")"));

        NTIC_EXPECT_TRUE(
            stream.atEnd());
    }

    //----------------------------------------------------------
    // reset()
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "Pa");

        TokenStream stream(tokens);

        NTIC_EXPECT_TRUE(
            stream.match("Pa"));

        NTIC_EXPECT_TRUE(
            stream.atEnd());

        stream.reset();

        NTIC_EXPECT_TRUE(
            !stream.atEnd());

        NTIC_EXPECT_TRUE(
            stream.current().text()
            == "Pa");

        NTIC_EXPECT_TRUE(
            stream.position() == 0);
    }

    //----------------------------------------------------------
    // End Token Inspection
    //----------------------------------------------------------

    {
        const auto tokens =
            lexer.tokenize(
                "kg");

        TokenStream stream(tokens);

        stream.consume();

        NTIC_EXPECT_TRUE(
            stream.atEnd());

        NTIC_EXPECT_TRUE(
            stream.current().is(
                TokenType::End));

        NTIC_EXPECT_TRUE(
            stream.peek().is(
                TokenType::End));

        NTIC_EXPECT_TRUE(
            stream.peek(100).is(
                TokenType::End));
    }

    //----------------------------------------------------------

    ntic::lbm::test::printSummary();

    return 0;
}
