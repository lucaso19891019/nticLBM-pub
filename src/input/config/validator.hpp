#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "model_alias_registry.hpp"
#include "parameter_table.hpp"

namespace ntic::lbm::config
{

class Validator
{
public:
    explicit
    Validator(
        ParameterTable& table);

    void
    validate();

    bool
    hasWarnings() const;

    const std::vector<std::string>&
    warnings() const;

    void
    printWarnings(
        std::ostream& output) const;

private:
    void
    validateFluidProperties();

    void
    validateCharacteristicScales();

    void
    validateLatticeScaling();

    void
    validateSimulationControl();

    void
    validateModelSelection();

    void
    validateRecommendations();

    bool
    has(
        const std::string& canonicalName) const;

    double
    quantitySI(
        const std::string& canonicalName) const;

    std::string
    stringValue(
        const std::string& canonicalName) const;

    void
    setQuantitySI(
        const std::string& canonicalName,
        double value);

    void
    setString(
        const std::string& canonicalName,
        const std::string& value);

    void
    warn(
        const std::string& message);

    static bool
    nearlyEqual(
        double lhs,
        double rhs,
        double relativeTolerance = 1.0e-10,
        double absoluteTolerance = 1.0e-14);

    static void
    requirePositive(
        double value,
        const std::string& canonicalName);

    static long long
    checkedPositiveInteger(
        double value,
        const std::string& canonicalName);

private:
    ParameterTable& table_;

    ModelAliasRegistry modelAliasRegistry_;

    std::vector<std::string> warnings_;
};

} // namespace ntic::lbm::config
