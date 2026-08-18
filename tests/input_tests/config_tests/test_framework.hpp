#pragma once

#include <iostream>
#include <string>
#include <cstdlib>

#define NTIC_EXPECT_THROW(expr, exception_type)        \
{                                                      \
    bool caught = false;                               \
    try                                                \
    {                                                  \
        expr;                                          \
    }                                                  \
    catch(const exception_type&)                       \
    {                                                  \
        caught = true;                                 \
    }                                                  \
    ntic::lbm::test::expectTrue(                       \
        caught,                                        \
        #expr,                                         \
        __FILE__,                                      \
        __LINE__);                                     \
}


namespace ntic::lbm::test
{

inline int totalTests = 0;
inline int passedTests = 0;

inline void expectTrue(
    bool condition,
    const char* expression,
    const char* file,
    int line)
{
    ++totalTests;

    if(condition)
    {
        ++passedTests;
    }
    else
    {
        std::cerr
            << "\nFAILED\n"
            << "File : " << file << '\n'
            << "Line : " << line << '\n'
            << "Expr : " << expression << "\n\n";

        std::exit(EXIT_FAILURE);
    }
}

template<typename T>
inline void expectEqual(
    const T& lhs,
    const T& rhs,
    const char* lhsName,
    const char* rhsName,
    const char* file,
    int line)
{
    ++totalTests;

    if(lhs == rhs)
    {
        ++passedTests;
    }
    else
    {
        std::cerr
            << "\nFAILED\n"
            << "File : " << file << '\n'
            << "Line : " << line << '\n'
            << lhsName << " != " << rhsName
            << "\n\n";

        std::exit(EXIT_FAILURE);
    }
}

inline void printSummary()
{
    std::cout
        << "\n=================================\n"
        << "All tests passed.\n"
        << "Passed "
        << passedTests
        << " / "
        << totalTests
        << "\n=================================\n";
}

}

#define NTIC_EXPECT_TRUE(expr) \
    ntic::lbm::test::expectTrue((expr), #expr, __FILE__, __LINE__)

#define NTIC_EXPECT_EQ(a,b) \
    ntic::lbm::test::expectEqual((a),(b),#a,#b,__FILE__,__LINE__)
