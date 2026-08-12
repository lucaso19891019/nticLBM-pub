#pragma once

#include <cstddef>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <vector>

#include "lbm_parameter_registry.hpp"
#include "parameter_parser.hpp"
#include "parameter_value.hpp"
#include "prefix_registry.hpp"
#include "unit_registry.hpp"
#include "parameter_data.hpp"

namespace ntic::lbm::config
{

class Config
{
public:
    Config();

    void parse(const std::string& text);
    void parseFile(const std::string& filePath);

    bool contains(const std::string& identifier) const;
    bool containsBuiltin(const std::string& identifier) const;
    bool containsUserDefined(const std::string& identifier) const;

    std::size_t builtinCount() const;
    std::size_t userDefinedCount() const;
    std::size_t userDefinedBoundCount() const;
    bool hasUnboundUserDefinedParameters() const;

    void bindQuantity(double& variable, const std::string& identifier) const;
    void bindQuantity(float& variable, const std::string& identifier) const;
    void bindQuantity(int& variable, const std::string& identifier) const;
    void bindBool(bool& variable, const std::string& identifier) const;
    void bindString(std::string& variable, const std::string& identifier) const;
    void bindString(char* variable, std::size_t bufferSize,
                    const std::string& identifier) const;

    void bindUserDefinedQuantity(
        double& variable,
        const std::string& identifier,
        const PhysicalDimension& expectedDimension);

    void bindUserDefinedQuantity(
        float& variable,
        const std::string& identifier,
        const PhysicalDimension& expectedDimension);

    void bindUserDefinedQuantity(
        int& variable,
        const std::string& identifier,
        const PhysicalDimension& expectedDimension);

    void bindUserDefinedBool(
        bool& variable,
        const std::string& identifier);

    void bindUserDefinedString(
        std::string& variable,
        const std::string& identifier);

    void bindUserDefinedString(
        char* variable,
        std::size_t bufferSize,
        const std::string& identifier);

    void reportUnboundUserDefinedParameters(std::ostream& output) const;

    UnitRegistry& unitRegistry();
    PrefixRegistry& prefixRegistry();
    const LBMParameterRegistry& lbmParameterRegistry() const;
    ParameterData exportParameterData() const;

private:
    static std::string trim(const std::string& text);
    static std::string prepareLineForParameterParser(const std::string& line);

    void storeParsedParameter(
        const ParsedParameter& parsed,
        std::size_t lineNumber,
        std::unordered_map<std::string, ParameterValue>& builtinValues,
        std::vector<std::string>& userIdentifiers,
        std::vector<ParameterValue>& userValues,
        std::vector<bool>& userBound,
        std::vector<std::size_t>& userLineNumbers) const;

    void validateBuiltinValue(
        const std::string& canonicalName,
        const ParameterValue& value,
        std::size_t lineNumber) const;

    const ParameterValue& builtinValue(const std::string& identifier) const;
    std::size_t findUserDefinedIndex(const std::string& identifier) const;

    static double quantitySIValue(const Quantity& quantity);

    static int checkedIntegerValue(
        double value,
        const std::string& identifier,
        std::size_t lineNumber);

    static void copyStringToBuffer(
        char* destination,
        std::size_t bufferSize,
        const std::string& source,
        const std::string& identifier,
        std::size_t lineNumber);

    static std::string valueTypeName(const ParameterValue& value);

private:
    UnitRegistry unitRegistry_;
    PrefixRegistry prefixRegistry_;
    LBMParameterRegistry lbmParameterRegistry_;
    ParameterParser parameterParser_;

    std::unordered_map<std::string, ParameterValue> builtinValues_;

    std::vector<std::string> userIdentifiers_;
    std::vector<ParameterValue> userValues_;
    std::vector<bool> userBound_;
    std::vector<std::size_t> userLineNumbers_;
};

} // namespace ntic::lbm::config
