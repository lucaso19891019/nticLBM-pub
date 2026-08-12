#include "config.hpp"
#include "lbm_parameter_registry.hpp"
#include "parameter_table.hpp"
#include "si_unit_formatter.hpp"
#include "validator.hpp"

#include <array>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <variant>

using namespace ntic::lbm::config;

namespace
{

//==============================================================
// Canonical built-in parameter names
//
// Keep the order consistent with LBMParameterRegistry so the
// printed table is easy to read.
//==============================================================

constexpr std::array<const char*, 4> generalParameters =
{{
    "case_name",
    "precision",
    "geometry_file",
    "restart_file"
}};

constexpr std::array<const char*, 16> physicalParameters =
{{
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
    "symmetric_z"
}};

constexpr std::array<const char*, 11> latticeParameters =
{{
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
    "max_steps"
}};

constexpr std::array<const char*, 3> outputParameters =
{{
    "output_directory",
    "output_interval",
    "checkpoint_interval"
}};

//==============================================================
// Value printing
//==============================================================

void
printParameterValue(
    const ParameterValue& value)
{
    if (std::holds_alternative<Quantity>(value))
    {
        const Quantity& quantity =
            std::get<Quantity>(value);

        const double siValue =
            quantity.value()
            * quantity.unit().scale();

        const std::string unit =
            SIUnitFormatter::format(
                quantity.unit().dimension());

        std::cout
            << std::setprecision(16)
            << siValue;

        if (!unit.empty())
        {
            std::cout
                << " "
                << unit;
        }
    }
    else if (std::holds_alternative<bool>(value))
    {
        std::cout
            << (std::get<bool>(value)
                ? "true"
                : "false");
    }
    else
    {
        std::cout
            << '"'
            << std::get<std::string>(value)
            << '"';
    }
}

template<std::size_t N>
void
printBuiltinGroup(
    const std::string& title,
    const std::array<const char*, N>& names,
    const ParameterTable& table)
{
    std::cout
        << title
        << '\n'
        << std::string(title.size(), '-')
        << '\n';

    for (const char* name : names)
    {
        const BuiltinParameter& parameter =
            table.builtin(name);

        std::cout
            << std::left
            << std::setw(34)
            << name
            << " = ";

        if (!parameter.hasValue)
        {
            std::cout
                << "<not set>";
        }
        else
        {
            printParameterValue(
                parameter.value);
        }

        std::cout
            << '\n';
    }

    std::cout
        << '\n';
}

void
printBuiltinParameters(
    const ParameterTable& table)
{
    std::cout
        << "========================================\n"
        << "Built-in Parameters\n"
        << "========================================\n\n";

    printBuiltinGroup(
        "General",
        generalParameters,
        table);

    printBuiltinGroup(
        "Physical",
        physicalParameters,
        table);

    printBuiltinGroup(
        "Lattice",
        latticeParameters,
        table);

    printBuiltinGroup(
        "Output",
        outputParameters,
        table);
}

void
printUserDefinedParameters(
    const ParameterTable& table)
{
    std::cout
        << "========================================\n"
        << "User-defined Parameters\n"
        << "========================================\n\n";

    const std::vector<UserParameter>& parameters =
        table.userDefined();

    if (parameters.empty())
    {
        std::cout
            << "<none>\n\n";

        return;
    }

    for (const UserParameter& parameter : parameters)
    {
        std::cout
            << std::left
            << std::setw(34)
            << parameter.identifier
            << " = ";

        printParameterValue(
            parameter.value);

        std::cout
            << "  [line "
            << parameter.lineNumber
            << ", "
            << (parameter.bound
                ? "bound"
                : "unbound")
            << "]\n";
    }

    std::cout
        << '\n';
}

} // namespace

int
main(
    int argc,
    char** argv)
{
    if (argc != 2)
    {
        std::cerr
            << "Usage:\n"
            << "    "
            << argv[0]
            << " input.nticonf\n";

        return EXIT_FAILURE;
    }

    try
    {
        //------------------------------------------------------
        // Read the configuration file.
        //------------------------------------------------------

        Config config;

        config.parseFile(
            argv[1]);

        //------------------------------------------------------
        // Create the complete parameter table.
        //------------------------------------------------------

        LBMParameterRegistry registry;

        ParameterTable table(
            registry,
            config.exportParameterData());

        //------------------------------------------------------
        // Validate, repair and complete the built-in table.
        //------------------------------------------------------

        Validator validator(
            table);

        validator.validate();

        //------------------------------------------------------
        // Report validation result.
        //------------------------------------------------------

        std::cout
            << "========================================\n"
            << "Configuration validated successfully.\n"
            << "========================================\n";

        if (validator.hasWarnings())
        {
            std::cout
                << "\nWarnings\n"
                << "--------\n";

            validator.printWarnings(
                std::cout);
        }
        else
        {
            std::cout
                << "\nWarnings\n"
                << "--------\n"
                << "<none>\n";
        }

        std::cout
            << '\n';

        //------------------------------------------------------
        // Print the finalized table. This output logic belongs
        // only to this integration test and does not modify the
        // configuration library.
        //------------------------------------------------------

        printBuiltinParameters(
            table);

        printUserDefinedParameters(
            table);

        return EXIT_SUCCESS;
    }
    catch (const std::exception& error)
    {
        std::cerr
            << "========================================\n"
            << "Configuration validation failed.\n"
            << "========================================\n\n"
            << error.what()
            << '\n';

        return EXIT_FAILURE;
    }
}

