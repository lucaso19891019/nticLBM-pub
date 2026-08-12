#include "parameter_table.hpp"

#include <array>
#include <stdexcept>
#include <string>
#include <utility>

namespace ntic::lbm::config
{

namespace
{

//==============================================================
// Canonical built-in names
//
// LBMParameterRegistry currently does not expose iteration over
// its registered definitions. ParameterTable therefore keeps
// the frozen canonical-name list used by the current registry.
//
// Every name is checked against LBMParameterRegistry during
// construction, so a mismatch fails immediately.
//==============================================================

constexpr std::array<const char*, 34>
builtinCanonicalNames =
{{
    //----------------------------------------------------------
    // General
    //----------------------------------------------------------

    "case_name",
    "precision",
    "geometry_file",
    "restart_file",

    //----------------------------------------------------------
    // Physical
    //----------------------------------------------------------

    "density",
    "dynamic_viscosity",
    "kinematic_viscosity",
    "characteristic_length",
    "characteristic_velocity",
    "reynolds_number",
    "physical_time",
    "restart_time",
    "time_step",
    "grid_spacing",
    "periodic_x",
    "periodic_y",
    "periodic_z",
    "symmetric_x",
    "symmetric_y",
    "symmetric_z",

    //----------------------------------------------------------
    // Lattice
    //----------------------------------------------------------

    "lattice_model",
    "collision_model",
    "relaxation_time",
    "relaxation_frequency",
    "lattice_characteristic_velocity",
    "mach_number",
    "number_of_levels",
    "wall_boundary",
    "inlet_boundary",
    "outlet_boundary",
    "max_steps",

    //----------------------------------------------------------
    // Output
    //----------------------------------------------------------

    "output_directory",
    "output_interval",
    "checkpoint_interval"
}};

} // namespace

//==============================================================
// Construction
//==============================================================

ParameterTable::ParameterTable(
    const LBMParameterRegistry& registry,
    const ParameterData& data)
    :
    builtin_(),
    userDefined_()
{
    initializeBuiltins(
        registry);

    importBuiltinValues(
        registry,
        data);

    importUserDefinedValues(
        data);
}

//==============================================================
// Built-in table access
//==============================================================

const BuiltinParameter&
ParameterTable::builtin(
    const std::string& canonicalName) const
{
    const auto iter =
        builtin_.find(
            canonicalName);

    if (iter == builtin_.end())
    {
        throw std::runtime_error(
            "ParameterTable: unknown canonical built-in "
            "parameter '"
            + canonicalName
            + "'.");
    }

    return iter->second;
}

BuiltinParameter&
ParameterTable::builtin(
    const std::string& canonicalName)
{
    const auto iter =
        builtin_.find(
            canonicalName);

    if (iter == builtin_.end())
    {
        throw std::runtime_error(
            "ParameterTable: unknown canonical built-in "
            "parameter '"
            + canonicalName
            + "'.");
    }

    return iter->second;
}

//==============================================================
// User-defined table access
//==============================================================

const std::vector<UserParameter>&
ParameterTable::userDefined() const
{
    return userDefined_;
}

std::vector<UserParameter>&
ParameterTable::userDefined()
{
    return userDefined_;
}

//==============================================================
// Construction helpers
//==============================================================

void
ParameterTable::initializeBuiltins(
    const LBMParameterRegistry& registry)
{
    if (registry.size()
        != builtinCanonicalNames.size())
    {
        throw std::runtime_error(
            "ParameterTable: LBMParameterRegistry contains "
            + std::to_string(registry.size())
            + " parameters, but ParameterTable expects "
            + std::to_string(
                builtinCanonicalNames.size())
            + ".");
    }

    builtin_.reserve(
        builtinCanonicalNames.size());

    for (const char* canonicalName
         : builtinCanonicalNames)
    {
        if (!registry.contains(
                canonicalName))
        {
            throw std::runtime_error(
                "ParameterTable: LBMParameterRegistry does not "
                "contain expected canonical parameter '"
                + std::string(canonicalName)
                + "'.");
        }

        const std::string& resolvedName =
            registry.resolve(
                canonicalName);

        if (resolvedName
            != canonicalName)
        {
            throw std::runtime_error(
                "ParameterTable: expected canonical parameter '"
                + std::string(canonicalName)
                + "' resolves to '"
                + resolvedName
                + "'.");
        }

        const ParameterInfo& info =
            registry.find(
                canonicalName);

        const auto result =
            builtin_.emplace(
                canonicalName,
                BuiltinParameter
                {
                    &info,
                    false,
                    ParameterValue()
                });

        if (!result.second)
        {
            throw std::runtime_error(
                "ParameterTable: duplicate canonical parameter '"
                + std::string(canonicalName)
                + "'.");
        }
    }
}

void
ParameterTable::importBuiltinValues(
    const LBMParameterRegistry& registry,
    const ParameterData& data)
{
    for (const auto& item
         : data.builtinValues)
    {
        const std::string& canonicalName =
            item.first;

        const ParameterValue& value =
            item.second;

        //------------------------------------------------------
        // ParameterData must already contain canonical names.
        //------------------------------------------------------

        if (!registry.contains(
                canonicalName))
        {
            throw std::runtime_error(
                "ParameterTable: ParameterData contains unknown "
                "built-in parameter '"
                + canonicalName
                + "'.");
        }

        const std::string& resolvedName =
            registry.resolve(
                canonicalName);

        if (resolvedName
            != canonicalName)
        {
            throw std::runtime_error(
                "ParameterTable: ParameterData built-in key '"
                + canonicalName
                + "' is not canonical; it resolves to '"
                + resolvedName
                + "'.");
        }

        auto iter =
            builtin_.find(
                canonicalName);

        if (iter == builtin_.end())
        {
            throw std::runtime_error(
                "ParameterTable: canonical parameter '"
                + canonicalName
                + "' is registered but absent from the "
                  "static built-in table.");
        }

        iter->second.value =
            value;

        iter->second.hasValue =
            true;
    }
}

void
ParameterTable::importUserDefinedValues(
    const ParameterData& data)
{
    userDefined_.reserve(
        data.userDefinedValues.size());

    for (const UserDefinedParameterData& parameter
         : data.userDefinedValues)
    {
        userDefined_.push_back(
            UserParameter
            {
                parameter.identifier,
                parameter.value,
                parameter.bound,
                parameter.lineNumber
            });
    }
}

} // namespace ntic::lbm::config
