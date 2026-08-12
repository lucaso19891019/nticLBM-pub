#include "si_unit_formatter.hpp"

#include <string>
#include <vector>

namespace ntic::lbm::config
{

namespace
{

//----------------------------------------------------------
// Unit factor
//----------------------------------------------------------

struct UnitFactor
{
    std::string symbol;

    common::Rational exponent;
};

//----------------------------------------------------------
// Check fractional exponent
//----------------------------------------------------------

bool
hasFractionalExponent(
    const std::vector<UnitFactor>& factors)
{
    for (const UnitFactor& factor : factors)
    {
        if (!factor.exponent.isInteger())
        {
            return true;
        }
    }

    return false;
}

//----------------------------------------------------------
// Format exponent
//----------------------------------------------------------

std::string
formatExponent(
    const common::Rational& exponent)
{
    if (exponent
        == common::Rational(1))
    {
        return "";
    }

    if (exponent.isInteger())
    {
        return "^"
            + std::to_string(
                exponent.numerator());
    }

    return "^"
        + std::to_string(
            exponent.numerator())
        + "/"
        + std::to_string(
            exponent.denominator());
}

//----------------------------------------------------------
// Format unit factor
//----------------------------------------------------------

std::string
formatFactor(
    const UnitFactor& factor)
{
    return factor.symbol
        + formatExponent(
            factor.exponent);
}

//----------------------------------------------------------
// Join unit factors
//----------------------------------------------------------

std::string
joinFactors(
    const std::vector<UnitFactor>& factors)
{
    std::string result;

    for (std::size_t i = 0;
         i < factors.size();
         ++i)
    {
        if (i != 0)
        {
            result += "*";
        }

        result +=
            formatFactor(
                factors[i]);
    }

    return result;
}

//----------------------------------------------------------
// Format numerator or denominator
//----------------------------------------------------------

std::string
formatGroup(
    const std::vector<UnitFactor>& factors)
{
    const std::string text =
        joinFactors(
            factors);

    //------------------------------------------------------
    // Parentheses are required when:
    //
    // 1. multiple unit factors are present, or
    // 2. a fractional exponent is present.
    //------------------------------------------------------

    if (factors.size() > 1
        || hasFractionalExponent(factors))
    {
        return "("
            + text
            + ")";
    }

    return text;
}

} // namespace

//----------------------------------------------------------
// Format SI base units
//----------------------------------------------------------

std::string
SIUnitFormatter::format(
    const PhysicalDimension& dimension)
{
    using BaseDimension =
        PhysicalDimension::BaseDimension;

    if (dimension.isDimensionless())
    {
        return "";
    }

    std::vector<UnitFactor> positive;
    std::vector<UnitFactor> negative;

    //------------------------------------------------------
    // Add a base SI unit.
    //------------------------------------------------------

    const auto addUnit =
        [&positive, &negative](
            const std::string& symbol,
            const common::Rational& exponent)
        {
            if (exponent.isZero())
            {
                return;
            }

            if (exponent
                > common::Rational(0))
            {
                positive.push_back(
                    {
                        symbol,
                        exponent
                    });
            }
            else
            {
                negative.push_back(
                    {
                        symbol,
                        -exponent
                    });
            }
        };

    //------------------------------------------------------
    // SI base-unit display order:
    //
    // kg, m, s, K
    //------------------------------------------------------

    addUnit(
        "kg",
        dimension[
            BaseDimension::Mass]);

    addUnit(
        "m",
        dimension[
            BaseDimension::Length]);

    addUnit(
        "s",
        dimension[
            BaseDimension::Time]);

    addUnit(
        "K",
        dimension[
            BaseDimension::Temperature]);

    //------------------------------------------------------
    // Pure negative dimension.
    //
    // Preserve negative exponents instead of producing
    // "1 / ...".
    //
    // Examples:
    //
    //     T^-1       -> s^-1
    //     L^-2 T^-1  -> m^-2*s^-1
    //------------------------------------------------------

    if (positive.empty())
    {
        std::vector<UnitFactor> signedFactors;

        signedFactors.reserve(
            negative.size());

        for (const UnitFactor& factor
             : negative)
        {
            signedFactors.push_back(
                {
                    factor.symbol,
                    -factor.exponent
                });
        }

        return joinFactors(
            signedFactors);
    }

    //------------------------------------------------------
    // Positive dimensions only.
    //------------------------------------------------------

    if (negative.empty())
    {
        return formatGroup(
            positive);
    }

    //------------------------------------------------------
    // Mixed positive / negative dimensions.
    //------------------------------------------------------

    return formatGroup(
        positive)
        + "/"
        + formatGroup(
            negative);
}

} // namespace ntic::lbm::config
