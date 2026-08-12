#include "prefix.hpp"
#include "prefix_registry.hpp"
#include "test_framework.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using ntic::lbm::config::Prefix;
using ntic::lbm::config::PrefixRegistry;

namespace
{

bool nearlyEqual(
    double lhs,
    double rhs,
    double relativeTolerance = 1.0e-12)
{
    const double scale = std::max(
        {1.0, std::abs(lhs), std::abs(rhs)});

    return std::abs(lhs - rhs)
        <= relativeTolerance * scale;
}

} // namespace

int main()
{
    //----------------------------------------------------------
    // Default Prefix
    //----------------------------------------------------------

    {
        Prefix prefix;

        NTIC_EXPECT_TRUE(
            prefix.symbol().empty());

        NTIC_EXPECT_TRUE(
            prefix.name().empty());

        NTIC_EXPECT_TRUE(
            nearlyEqual(prefix.factor(), 1.0));
    }

    //----------------------------------------------------------
    // General Constructor
    //----------------------------------------------------------

    {
        Prefix prefix(
            "k",
            "kilo",
            1.0e3);

        NTIC_EXPECT_TRUE(
            prefix.symbol() == "k");

        NTIC_EXPECT_TRUE(
            prefix.name() == "kilo");

        NTIC_EXPECT_TRUE(
            nearlyEqual(prefix.factor(), 1.0e3));
    }

    //----------------------------------------------------------
    // Invalid Constructor Arguments
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        Prefix(
            "",
            "invalid",
            1.0),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        Prefix(
            "x",
            "",
            1.0),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        Prefix(
            "x",
            "invalid",
            0.0),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        Prefix(
            "x",
            "invalid",
            -1.0),
        std::invalid_argument);

    //----------------------------------------------------------
    // Prefix Comparison
    //----------------------------------------------------------

    {
        const Prefix kiloA(
            "k",
            "kilo",
            1.0e3);

        const Prefix kiloB(
            "k",
            "kilo",
            1.0e3);

        const Prefix milli(
            "m",
            "milli",
            1.0e-3);

        NTIC_EXPECT_TRUE(
            kiloA == kiloB);

        NTIC_EXPECT_TRUE(
            kiloA != milli);
    }

    //----------------------------------------------------------
    // Prefix toString()
    //----------------------------------------------------------

    {
        const std::string text =
            Prefix(
                "M",
                "mega",
                1.0e6).toString();

        NTIC_EXPECT_TRUE(
            text.find("M")
            != std::string::npos);

        NTIC_EXPECT_TRUE(
            text.find("mega")
            != std::string::npos);

        NTIC_EXPECT_TRUE(
            text.find("factor=")
            != std::string::npos);
    }

    //----------------------------------------------------------
    // Built-in Prefix Registry
    //----------------------------------------------------------

    PrefixRegistry registry;

    // 12 positive powers and 12 negative powers.
    NTIC_EXPECT_TRUE(
        registry.size() == 24);

    //----------------------------------------------------------
    // Positive-power SI Prefixes
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(registry.contains("Q"));
    NTIC_EXPECT_TRUE(registry.contains("R"));
    NTIC_EXPECT_TRUE(registry.contains("Y"));
    NTIC_EXPECT_TRUE(registry.contains("Z"));
    NTIC_EXPECT_TRUE(registry.contains("E"));
    NTIC_EXPECT_TRUE(registry.contains("P"));
    NTIC_EXPECT_TRUE(registry.contains("T"));
    NTIC_EXPECT_TRUE(registry.contains("G"));
    NTIC_EXPECT_TRUE(registry.contains("M"));
    NTIC_EXPECT_TRUE(registry.contains("k"));
    NTIC_EXPECT_TRUE(registry.contains("h"));
    NTIC_EXPECT_TRUE(registry.contains("da"));

    //----------------------------------------------------------
    // Negative-power SI Prefixes
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(registry.contains("d"));
    NTIC_EXPECT_TRUE(registry.contains("c"));
    NTIC_EXPECT_TRUE(registry.contains("m"));
    NTIC_EXPECT_TRUE(registry.contains("u"));
    NTIC_EXPECT_TRUE(registry.contains("n"));
    NTIC_EXPECT_TRUE(registry.contains("p"));
    NTIC_EXPECT_TRUE(registry.contains("f"));
    NTIC_EXPECT_TRUE(registry.contains("a"));
    NTIC_EXPECT_TRUE(registry.contains("z"));
    NTIC_EXPECT_TRUE(registry.contains("y"));
    NTIC_EXPECT_TRUE(registry.contains("r"));
    NTIC_EXPECT_TRUE(registry.contains("q"));

    //----------------------------------------------------------
    // Unknown Symbols
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(
        !registry.contains(""));

    NTIC_EXPECT_TRUE(
        !registry.contains("K"));

    NTIC_EXPECT_TRUE(
        !registry.contains("x"));

    NTIC_EXPECT_TRUE(
        !registry.contains("micro"));

    //----------------------------------------------------------
    // Find Positive-power Prefixes
    //----------------------------------------------------------

    {
        const Prefix& quetta =
            registry.find("Q");

        NTIC_EXPECT_TRUE(
            quetta.name() == "quetta");

        NTIC_EXPECT_TRUE(
            nearlyEqual(
                quetta.factor(),
                1.0e30));
    }

    {
        const Prefix& ronna =
            registry.find("R");

        NTIC_EXPECT_TRUE(
            ronna.name() == "ronna");

        NTIC_EXPECT_TRUE(
            nearlyEqual(
                ronna.factor(),
                1.0e27));
    }

    {
        const Prefix& kilo =
            registry.find("k");

        NTIC_EXPECT_TRUE(
            kilo.symbol() == "k");

        NTIC_EXPECT_TRUE(
            kilo.name() == "kilo");

        NTIC_EXPECT_TRUE(
            nearlyEqual(
                kilo.factor(),
                1.0e3));
    }

    {
        const Prefix& deca =
            registry.find("da");

        NTIC_EXPECT_TRUE(
            deca.symbol() == "da");

        NTIC_EXPECT_TRUE(
            deca.name() == "deca");

        NTIC_EXPECT_TRUE(
            nearlyEqual(
                deca.factor(),
                1.0e1));
    }

    //----------------------------------------------------------
    // Find Negative-power Prefixes
    //----------------------------------------------------------

    {
        const Prefix& milli =
            registry.find("m");

        NTIC_EXPECT_TRUE(
            milli.name() == "milli");

        NTIC_EXPECT_TRUE(
            nearlyEqual(
                milli.factor(),
                1.0e-3));
    }

    {
        const Prefix& micro =
            registry.find("u");

        NTIC_EXPECT_TRUE(
            micro.name() == "micro");

        NTIC_EXPECT_TRUE(
            nearlyEqual(
                micro.factor(),
                1.0e-6));
    }

    {
        const Prefix& quecto =
            registry.find("q");

        NTIC_EXPECT_TRUE(
            quecto.name() == "quecto");

        NTIC_EXPECT_TRUE(
            nearlyEqual(
                quecto.factor(),
                1.0e-30));
    }

    //----------------------------------------------------------
    // Verify All Built-in Factors
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("Q").factor(), 1.0e30));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("R").factor(), 1.0e27));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("Y").factor(), 1.0e24));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("Z").factor(), 1.0e21));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("E").factor(), 1.0e18));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("P").factor(), 1.0e15));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("T").factor(), 1.0e12));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("G").factor(), 1.0e9));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("M").factor(), 1.0e6));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("k").factor(), 1.0e3));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("h").factor(), 1.0e2));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("da").factor(), 1.0e1));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("d").factor(), 1.0e-1));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("c").factor(), 1.0e-2));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("m").factor(), 1.0e-3));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("u").factor(), 1.0e-6));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("n").factor(), 1.0e-9));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("p").factor(), 1.0e-12));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("f").factor(), 1.0e-15));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("a").factor(), 1.0e-18));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("z").factor(), 1.0e-21));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("y").factor(), 1.0e-24));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("r").factor(), 1.0e-27));

    NTIC_EXPECT_TRUE(
        nearlyEqual(registry.find("q").factor(), 1.0e-30));

    //----------------------------------------------------------
    // Symbols Ordered by Descending Length
    //----------------------------------------------------------

    {
        const std::vector<std::string> symbols =
            registry.symbolsByDescendingLength();

        NTIC_EXPECT_TRUE(
            symbols.size() == registry.size());

        NTIC_EXPECT_TRUE(
            !symbols.empty());

        // "da" is currently the only two-character prefix,
        // so it must appear first.
        NTIC_EXPECT_TRUE(
            symbols.front() == "da");

        const auto decaIter =
            std::find(
                symbols.begin(),
                symbols.end(),
                "da");

        const auto deciIter =
            std::find(
                symbols.begin(),
                symbols.end(),
                "d");

        NTIC_EXPECT_TRUE(
            decaIter != symbols.end());

        NTIC_EXPECT_TRUE(
            deciIter != symbols.end());

        NTIC_EXPECT_TRUE(
            decaIter < deciIter);

        for (std::size_t i = 1;
             i < symbols.size();
             ++i)
        {
            NTIC_EXPECT_TRUE(
                symbols[i - 1].size()
                >= symbols[i].size());
        }
    }

    //----------------------------------------------------------
    // Register Custom Prefix
    //----------------------------------------------------------

    {
        const Prefix custom(
            "x",
            "custom",
            2.5);

        const std::size_t oldSize =
            registry.size();

        registry.registerPrefix(custom);

        NTIC_EXPECT_TRUE(
            registry.size() == oldSize + 1);

        NTIC_EXPECT_TRUE(
            registry.contains("x"));

        NTIC_EXPECT_EQ(
            registry.find("x"),
            custom);
    }

    //----------------------------------------------------------
    // Duplicate Registration
    //----------------------------------------------------------

    {
        const std::size_t oldSize =
            registry.size();

        NTIC_EXPECT_THROW(
            registry.registerPrefix(
                Prefix(
                    "k",
                    "another-kilo",
                    1.0e3)),
            std::runtime_error);

        NTIC_EXPECT_TRUE(
            registry.size() == oldSize);

        NTIC_EXPECT_TRUE(
            registry.find("k").name()
            == "kilo");

        NTIC_EXPECT_TRUE(
            nearlyEqual(
                registry.find("k").factor(),
                1.0e3));
    }

    //----------------------------------------------------------
    // Unknown Prefix
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        registry.find("unknown"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        registry.find("K"),
        std::runtime_error);

    //----------------------------------------------------------
    // clear()
    //----------------------------------------------------------

    registry.clear();

    NTIC_EXPECT_TRUE(
        registry.size() == 0);

    NTIC_EXPECT_TRUE(
        !registry.contains("k"));

    NTIC_EXPECT_TRUE(
        !registry.contains("m"));

    NTIC_EXPECT_TRUE(
        !registry.contains("da"));

    NTIC_EXPECT_THROW(
        registry.find("k"),
        std::runtime_error);

    NTIC_EXPECT_TRUE(
        registry.symbolsByDescendingLength().empty());

    //----------------------------------------------------------
    // Register After clear()
    //----------------------------------------------------------

    {
        const Prefix kilo(
            "k",
            "kilo",
            1.0e3);

        registry.registerPrefix(kilo);

        NTIC_EXPECT_TRUE(
            registry.size() == 1);

        NTIC_EXPECT_TRUE(
            registry.contains("k"));

        NTIC_EXPECT_EQ(
            registry.find("k"),
            kilo);
    }

    //----------------------------------------------------------

    ntic::lbm::test::printSummary();

    return 0;
}
