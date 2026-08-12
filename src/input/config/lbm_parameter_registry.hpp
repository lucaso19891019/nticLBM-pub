#pragma once

#include <cstddef>
#include <initializer_list>
#include <string>
#include <unordered_map>

#include "parameter_alias_registry.hpp"
#include "parameter_info.hpp"
#include "parameter_registry.hpp"
#include "parameter_value_type.hpp"

namespace ntic::lbm::config
{

class LBMParameterRegistry
{
public:
    //----------------------------------------------------------
    // Construction
    //----------------------------------------------------------

    LBMParameterRegistry();

    //----------------------------------------------------------
    // Query
    //
    // The name passed to these functions may be either a
    // canonical name or a registered alias.
    //----------------------------------------------------------

    bool
    contains(
        const std::string& name) const;

    const std::string&
    resolve(
        const std::string& name) const;

    const ParameterInfo&
    find(
        const std::string& name) const;

    ParameterValueType
    valueType(
        const std::string& name) const;

    std::size_t
    size() const;

    //----------------------------------------------------------
    // Underlying registries
    //----------------------------------------------------------

    const ParameterRegistry&
    parameterRegistry() const;

    const ParameterAliasRegistry&
    aliasRegistry() const;

private:
    //----------------------------------------------------------
    // Built-in registration
    //----------------------------------------------------------

    void
    registerBuiltins();

    void
    registerParameter(
        const std::string& canonicalName,
        const PhysicalDimension& dimension,
        ParameterCategory category,
        ParameterValueType valueType,
        const std::string& description,
        std::initializer_list<std::string> aliases = {});

private:
    ParameterRegistry parameterRegistry_;

    ParameterAliasRegistry aliasRegistry_;

    std::unordered_map<
        std::string,
        ParameterValueType> valueTypes_;
};

} // namespace ntic::lbm::config
