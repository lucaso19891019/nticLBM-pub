#include "config.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <variant>

namespace ntic::lbm::config
{

Config::Config()
    : unitRegistry_(),
      prefixRegistry_(),
      lbmParameterRegistry_(),
      parameterParser_(unitRegistry_, prefixRegistry_)
{
}

void Config::parse(const std::string& text)
{
    std::unordered_map<std::string, ParameterValue> builtinValues;
    std::vector<std::string> userIdentifiers;
    std::vector<ParameterValue> userValues;
    std::vector<bool> userBound;
    std::vector<std::size_t> userLineNumbers;

    std::istringstream input(text);
    std::string line;
    std::size_t lineNumber = 0;

    while (std::getline(input, line))
    {
        ++lineNumber;
        const std::string stripped = trim(line);

        if (stripped.empty() || stripped.front() == '#')
            continue;

        try
        {
            const ParsedParameter parsed = parameterParser_.parse(
                prepareLineForParameterParser(line));

            storeParsedParameter(
                parsed,
                lineNumber,
                builtinValues,
                userIdentifiers,
                userValues,
                userBound,
                userLineNumbers);
        }
        catch (const std::exception& error)
        {
            throw std::runtime_error(
                "Config: error at line " + std::to_string(lineNumber)
                + ": " + error.what());
        }
    }

    builtinValues_ = std::move(builtinValues);
    userIdentifiers_ = std::move(userIdentifiers);
    userValues_ = std::move(userValues);
    userBound_ = std::move(userBound);
    userLineNumbers_ = std::move(userLineNumbers);
}

void Config::parseFile(const std::string& filePath)
{
    std::ifstream input(filePath);

    if (!input)
        throw std::runtime_error(
            "Config: unable to open configuration file '" + filePath + "'.");

    std::ostringstream buffer;
    buffer << input.rdbuf();

    if (!input.good() && !input.eof())
        throw std::runtime_error(
            "Config: failed while reading configuration file '"
            + filePath + "'.");

    parse(buffer.str());
}

bool Config::contains(const std::string& identifier) const
{
    return containsBuiltin(identifier) || containsUserDefined(identifier);
}

bool Config::containsBuiltin(const std::string& identifier) const
{
    if (!lbmParameterRegistry_.contains(identifier))
        return false;

    const std::string& canonical = lbmParameterRegistry_.resolve(identifier);
    return builtinValues_.find(canonical) != builtinValues_.end();
}

bool Config::containsUserDefined(const std::string& identifier) const
{
    const std::string normalized = ParameterAliasRegistry::normalize(identifier);

    return std::find(
        userIdentifiers_.begin(),
        userIdentifiers_.end(),
        normalized) != userIdentifiers_.end();
}

std::size_t Config::builtinCount() const
{
    return builtinValues_.size();
}

std::size_t Config::userDefinedCount() const
{
    return userIdentifiers_.size();
}

std::size_t Config::userDefinedBoundCount() const
{
    return static_cast<std::size_t>(
        std::count(userBound_.begin(), userBound_.end(), true));
}

bool Config::hasUnboundUserDefinedParameters() const
{
    return std::find(userBound_.begin(), userBound_.end(), false)
        != userBound_.end();
}

void Config::bindQuantity(double& variable, const std::string& identifier) const
{
    const ParameterValue& value = builtinValue(identifier);

    if (!std::holds_alternative<Quantity>(value))
        throw std::runtime_error(
            "Config: built-in parameter '" + identifier + "' has type "
            + valueTypeName(value) + ", but Quantity was expected.");

    variable = quantitySIValue(std::get<Quantity>(value));
}

void Config::bindQuantity(float& variable, const std::string& identifier) const
{
    double value = 0.0;
    bindQuantity(value, identifier);

    if (!std::isfinite(value)
        || std::abs(value) > static_cast<double>(std::numeric_limits<float>::max()))
        throw std::runtime_error(
            "Config: built-in quantity '" + identifier
            + "' cannot be represented as float.");

    variable = static_cast<float>(value);
}

void Config::bindQuantity(int& variable, const std::string& identifier) const
{
    const ParameterValue& value = builtinValue(identifier);

    if (!std::holds_alternative<Quantity>(value))
        throw std::runtime_error(
            "Config: built-in parameter '" + identifier + "' has type "
            + valueTypeName(value) + ", but Quantity was expected.");

    variable = checkedIntegerValue(
        quantitySIValue(std::get<Quantity>(value)), identifier, 0);
}

void Config::bindBool(bool& variable, const std::string& identifier) const
{
    const ParameterValue& value = builtinValue(identifier);

    if (!std::holds_alternative<bool>(value))
        throw std::runtime_error(
            "Config: built-in parameter '" + identifier + "' has type "
            + valueTypeName(value) + ", but bool was expected.");

    variable = std::get<bool>(value);
}

void Config::bindString(
    std::string& variable,
    const std::string& identifier) const
{
    const ParameterValue& value = builtinValue(identifier);

    if (!std::holds_alternative<std::string>(value))
        throw std::runtime_error(
            "Config: built-in parameter '" + identifier + "' has type "
            + valueTypeName(value) + ", but string was expected.");

    variable = std::get<std::string>(value);
}

void Config::bindString(
    char* variable,
    std::size_t bufferSize,
    const std::string& identifier) const
{
    const ParameterValue& value = builtinValue(identifier);

    if (!std::holds_alternative<std::string>(value))
        throw std::runtime_error(
            "Config: built-in parameter '" + identifier + "' has type "
            + valueTypeName(value) + ", but string was expected.");

    copyStringToBuffer(
        variable, bufferSize, std::get<std::string>(value), identifier, 0);
}

void Config::bindUserDefinedQuantity(
    double& variable,
    const std::string& identifier,
    const PhysicalDimension& expectedDimension)
{
    const std::size_t index = findUserDefinedIndex(identifier);
    const ParameterValue& value = userValues_[index];

    if (!std::holds_alternative<Quantity>(value))
        throw std::runtime_error(
            "Config: user-defined parameter '" + userIdentifiers_[index]
            + "' at line " + std::to_string(userLineNumbers_[index])
            + " has type " + valueTypeName(value)
            + ", but Quantity was expected.");

    const Quantity& quantity = std::get<Quantity>(value);

    if (quantity.unit().dimension() != expectedDimension)
        throw std::runtime_error(
            "Config: user-defined quantity '" + userIdentifiers_[index]
            + "' at line " + std::to_string(userLineNumbers_[index])
            + " has dimension " + quantity.unit().dimension().toString()
            + ", but " + expectedDimension.toString() + " was expected.");

    variable = quantitySIValue(quantity);
    userBound_[index] = true;
}

void Config::bindUserDefinedQuantity(
    float& variable,
    const std::string& identifier,
    const PhysicalDimension& expectedDimension)
{
    const std::size_t index = findUserDefinedIndex(identifier);
    const ParameterValue& stored = userValues_[index];

    if (!std::holds_alternative<Quantity>(stored))
        throw std::runtime_error(
            "Config: user-defined parameter '" + userIdentifiers_[index]
            + "' at line " + std::to_string(userLineNumbers_[index])
            + " has type " + valueTypeName(stored)
            + ", but Quantity was expected.");

    const Quantity& quantity = std::get<Quantity>(stored);

    if (quantity.unit().dimension() != expectedDimension)
        throw std::runtime_error(
            "Config: user-defined quantity '" + userIdentifiers_[index]
            + "' at line " + std::to_string(userLineNumbers_[index])
            + " has an unexpected physical dimension.");

    const double value = quantitySIValue(quantity);

    if (!std::isfinite(value)
        || std::abs(value) > static_cast<double>(std::numeric_limits<float>::max()))
        throw std::runtime_error(
            "Config: user-defined quantity '" + userIdentifiers_[index]
            + "' at line " + std::to_string(userLineNumbers_[index])
            + " cannot be represented as float.");

    variable = static_cast<float>(value);
    userBound_[index] = true;
}

void Config::bindUserDefinedQuantity(
    int& variable,
    const std::string& identifier,
    const PhysicalDimension& expectedDimension)
{
    const std::size_t index = findUserDefinedIndex(identifier);
    const ParameterValue& value = userValues_[index];

    if (!std::holds_alternative<Quantity>(value))
        throw std::runtime_error(
            "Config: user-defined parameter '" + userIdentifiers_[index]
            + "' at line " + std::to_string(userLineNumbers_[index])
            + " has type " + valueTypeName(value)
            + ", but Quantity was expected.");

    const Quantity& quantity = std::get<Quantity>(value);

    if (quantity.unit().dimension() != expectedDimension)
        throw std::runtime_error(
            "Config: user-defined quantity '" + userIdentifiers_[index]
            + "' at line " + std::to_string(userLineNumbers_[index])
            + " has an unexpected physical dimension.");

    variable = checkedIntegerValue(
        quantitySIValue(quantity),
        userIdentifiers_[index],
        userLineNumbers_[index]);

    userBound_[index] = true;
}

void Config::bindUserDefinedBool(
    bool& variable,
    const std::string& identifier)
{
    const std::size_t index = findUserDefinedIndex(identifier);
    const ParameterValue& value = userValues_[index];

    if (!std::holds_alternative<bool>(value))
        throw std::runtime_error(
            "Config: user-defined parameter '" + userIdentifiers_[index]
            + "' at line " + std::to_string(userLineNumbers_[index])
            + " has type " + valueTypeName(value)
            + ", but bool was expected.");

    variable = std::get<bool>(value);
    userBound_[index] = true;
}

void Config::bindUserDefinedString(
    std::string& variable,
    const std::string& identifier)
{
    const std::size_t index = findUserDefinedIndex(identifier);
    const ParameterValue& value = userValues_[index];

    if (!std::holds_alternative<std::string>(value))
        throw std::runtime_error(
            "Config: user-defined parameter '" + userIdentifiers_[index]
            + "' at line " + std::to_string(userLineNumbers_[index])
            + " has type " + valueTypeName(value)
            + ", but string was expected.");

    variable = std::get<std::string>(value);
    userBound_[index] = true;
}

void Config::bindUserDefinedString(
    char* variable,
    std::size_t bufferSize,
    const std::string& identifier)
{
    const std::size_t index = findUserDefinedIndex(identifier);
    const ParameterValue& value = userValues_[index];

    if (!std::holds_alternative<std::string>(value))
        throw std::runtime_error(
            "Config: user-defined parameter '" + userIdentifiers_[index]
            + "' at line " + std::to_string(userLineNumbers_[index])
            + " has type " + valueTypeName(value)
            + ", but string was expected.");

    copyStringToBuffer(
        variable,
        bufferSize,
        std::get<std::string>(value),
        userIdentifiers_[index],
        userLineNumbers_[index]);

    userBound_[index] = true;
}

void Config::reportUnboundUserDefinedParameters(std::ostream& output) const
{
    bool wroteHeader = false;

    for (std::size_t index = 0; index < userIdentifiers_.size(); ++index)
    {
        if (userBound_[index])
            continue;

        if (!wroteHeader)
        {
            output << "Unbound user-defined parameters:\n";
            wroteHeader = true;
        }

        output << "  line " << userLineNumbers_[index] << ": "
               << userIdentifiers_[index] << " ("
               << valueTypeName(userValues_[index]) << ")\n";
    }
}

UnitRegistry& Config::unitRegistry()
{
    return unitRegistry_;
}

PrefixRegistry& Config::prefixRegistry()
{
    return prefixRegistry_;
}

const LBMParameterRegistry& Config::lbmParameterRegistry() const
{
    return lbmParameterRegistry_;
}

std::string Config::trim(const std::string& text)
{
    const std::size_t begin = text.find_first_not_of(" \t\r\n");

    if (begin == std::string::npos)
        return std::string();

    const std::size_t end = text.find_last_not_of(" \t\r\n");
    return text.substr(begin, end - begin + 1);
}

std::string Config::prepareLineForParameterParser(const std::string& line)
{
    bool inSingleQuote = false;
    bool inDoubleQuote = false;
    bool escaped = false;
    std::size_t separator = std::string::npos;

    for (std::size_t index = 0; index < line.size(); ++index)
    {
        const char character = line[index];

        if (escaped)
        {
            escaped = false;
            continue;
        }

        if (character == '\\' && (inSingleQuote || inDoubleQuote))
        {
            escaped = true;
            continue;
        }

        if (character == '\'' && !inDoubleQuote)
        {
            inSingleQuote = !inSingleQuote;
            continue;
        }

        if (character == '"' && !inSingleQuote)
        {
            inDoubleQuote = !inDoubleQuote;
            continue;
        }

        if (!inSingleQuote && !inDoubleQuote
            && (character == '=' || character == ':'))
        {
            separator = index;
            break;
        }
    }

    if (separator == std::string::npos)
        return line;

    const std::string rawIdentifier = trim(line.substr(0, separator));

    if (rawIdentifier.empty())
        return line;

    return ParameterAliasRegistry::normalize(rawIdentifier)
        + " = " + line.substr(separator + 1);
}

void Config::storeParsedParameter(
    const ParsedParameter& parsed,
    std::size_t lineNumber,
    std::unordered_map<std::string, ParameterValue>& builtinValues,
    std::vector<std::string>& userIdentifiers,
    std::vector<ParameterValue>& userValues,
    std::vector<bool>& userBound,
    std::vector<std::size_t>& userLineNumbers) const
{
    if (lbmParameterRegistry_.contains(parsed.name))
    {
        const std::string& canonical = lbmParameterRegistry_.resolve(parsed.name);

        if (builtinValues.find(canonical) != builtinValues.end())
            throw std::runtime_error(
                "duplicate built-in parameter '" + canonical + "'.");

        validateBuiltinValue(canonical, parsed.value, lineNumber);
        builtinValues.emplace(canonical, parsed.value);
        return;
    }

    const std::string normalized =
        ParameterAliasRegistry::normalize(parsed.name);

    for (std::size_t index = 0; index < userIdentifiers.size(); ++index)
    {
        if (userIdentifiers[index] == normalized)
            throw std::runtime_error(
                "duplicate user-defined parameter '" + normalized
                + "'; it was first defined at line "
                + std::to_string(userLineNumbers[index]) + ".");
    }

    userIdentifiers.push_back(normalized);
    userValues.push_back(parsed.value);
    userBound.push_back(false);
    userLineNumbers.push_back(lineNumber);
}

void Config::validateBuiltinValue(
    const std::string& canonicalName,
    const ParameterValue& value,
    std::size_t lineNumber) const
{
    const ParameterValueType expected =
        lbmParameterRegistry_.valueType(canonicalName);

    const bool typeMatches =
        (expected == ParameterValueType::Quantity
         && std::holds_alternative<Quantity>(value))
        || (expected == ParameterValueType::Bool
            && std::holds_alternative<bool>(value))
        || (expected == ParameterValueType::String
            && std::holds_alternative<std::string>(value));

    if (!typeMatches)
        throw std::runtime_error(
            "built-in parameter '" + canonicalName + "' at line "
            + std::to_string(lineNumber) + " has type "
            + valueTypeName(value)
            + ", which does not match its registered type.");

    if (expected != ParameterValueType::Quantity)
        return;

    const Quantity& quantity = std::get<Quantity>(value);
    const PhysicalDimension& required =
        lbmParameterRegistry_.find(canonicalName).dimension();

    if (quantity.unit().dimension() != required)
        throw std::runtime_error(
            "built-in quantity '" + canonicalName + "' at line "
            + std::to_string(lineNumber) + " has dimension "
            + quantity.unit().dimension().toString() + ", but "
            + required.toString() + " is required.");
}

const ParameterValue& Config::builtinValue(
    const std::string& identifier) const
{
    if (!lbmParameterRegistry_.contains(identifier))
        throw std::runtime_error(
            "Config: '" + identifier
            + "' is not a registered LBM parameter.");

    const std::string& canonical = lbmParameterRegistry_.resolve(identifier);
    const auto iterator = builtinValues_.find(canonical);

    if (iterator == builtinValues_.end())
        throw std::runtime_error(
            "Config: built-in parameter '" + canonical
            + "' was not provided.");

    return iterator->second;
}

std::size_t Config::findUserDefinedIndex(
    const std::string& identifier) const
{
    const std::string normalized = ParameterAliasRegistry::normalize(identifier);

    for (std::size_t index = 0; index < userIdentifiers_.size(); ++index)
    {
        if (userIdentifiers_[index] == normalized)
            return index;
    }

    throw std::runtime_error(
        "Config: user-defined parameter '" + identifier + "' was not found.");
}

double Config::quantitySIValue(const Quantity& quantity)
{
    return quantity.value() * quantity.unit().scale();
}

int Config::checkedIntegerValue(
    double value,
    const std::string& identifier,
    std::size_t lineNumber)
{
    if (!std::isfinite(value))
        throw std::runtime_error(
            "Config: quantity '" + identifier
            + "' is not finite and cannot be converted to int.");

    const double rounded = std::round(value);
    const double tolerance = 1.0e-12 * std::max(1.0, std::abs(value));

    if (std::abs(value - rounded) > tolerance)
        throw std::runtime_error(
            "Config: quantity '" + identifier + "'"
            + (lineNumber == 0 ? std::string()
                               : " at line " + std::to_string(lineNumber))
            + " cannot be converted to int without loss.");

    if (rounded < static_cast<double>(std::numeric_limits<int>::min())
        || rounded > static_cast<double>(std::numeric_limits<int>::max()))
        throw std::runtime_error(
            "Config: quantity '" + identifier + "' is outside int range.");

    return static_cast<int>(rounded);
}

void Config::copyStringToBuffer(
    char* destination,
    std::size_t bufferSize,
    const std::string& source,
    const std::string& identifier,
    std::size_t lineNumber)
{
    if (destination == nullptr)
        throw std::invalid_argument(
            "Config: destination buffer for parameter '"
            + identifier + "' is null.");

    if (bufferSize == 0)
        throw std::invalid_argument(
            "Config: destination buffer for parameter '"
            + identifier + "' has zero size.");

    if (source.size() + 1 > bufferSize)
        throw std::runtime_error(
            "Config: destination buffer for parameter '" + identifier + "'"
            + (lineNumber == 0 ? std::string()
                               : " at line " + std::to_string(lineNumber))
            + " is too small.");

    std::memcpy(destination, source.c_str(), source.size() + 1);
}

std::string Config::valueTypeName(const ParameterValue& value)
{
    if (std::holds_alternative<Quantity>(value))
        return "Quantity";

    if (std::holds_alternative<bool>(value))
        return "bool";

    return "string";
}

ParameterData Config::exportParameterData() const
{
    ParameterData data;

    data.builtinValues = builtinValues_;

    data.userDefinedValues.reserve(
        userIdentifiers_.size());

    for (std::size_t i = 0;
         i < userIdentifiers_.size();
         ++i)
    {
        data.userDefinedValues.push_back(
        {
            userIdentifiers_[i],
            userValues_[i],
            userBound_[i],
            userLineNumbers_[i]
        });
    }

    return data;
}

} // namespace ntic::lbm::config
