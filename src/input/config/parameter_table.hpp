#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include "lbm_parameter_registry.hpp"
#include "parameter_data.hpp"
#include "parameter_info.hpp"
#include "parameter_value.hpp"

namespace ntic::lbm::config
{

//==============================================================
// Built-in LBM parameter entry
//==============================================================

struct BuiltinParameter
{
    const ParameterInfo* info;

    bool hasValue;

    ParameterValue value;
};

//==============================================================
// User-defined parameter entry
//==============================================================

struct UserParameter
{
    std::string identifier;

    ParameterValue value;

    bool bound;

    std::size_t lineNumber;
};

//==============================================================
// ParameterTable
//
// This class is only a data container.
//
// It does not parse configuration files.
// It does not validate LBM relationships.
// It does not emit warnings.
// It does not derive or override values.
//==============================================================

class ParameterTable
{
public:
    //----------------------------------------------------------
    // Construction
    //----------------------------------------------------------

    ParameterTable(
        const LBMParameterRegistry& registry,
        const ParameterData& data);

    //----------------------------------------------------------
    // Built-in table access
    //
    // The identifier must be a canonical LBM parameter name.
    //----------------------------------------------------------

    const BuiltinParameter&
    builtin(
        const std::string& canonicalName) const;

    BuiltinParameter&
    builtin(
        const std::string& canonicalName);

    //----------------------------------------------------------
    // User-defined table access
    //----------------------------------------------------------

    const std::vector<UserParameter>&
    userDefined() const;

    std::vector<UserParameter>&
    userDefined();

private:
    //----------------------------------------------------------
    // Construction helpers
    //----------------------------------------------------------

    void
    initializeBuiltins(
        const LBMParameterRegistry& registry);

    void
    importBuiltinValues(
        const LBMParameterRegistry& registry,
        const ParameterData& data);

    void
    importUserDefinedValues(
        const ParameterData& data);

private:
    std::unordered_map<
        std::string,
        BuiltinParameter> builtin_;

    std::vector<UserParameter> userDefined_;
};

} // namespace ntic::lbm::config
