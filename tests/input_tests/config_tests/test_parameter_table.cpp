#include "config.hpp"
#include "parameter_data.hpp"
#include "parameter_table.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>
#include <variant>

using namespace ntic::lbm::config;

namespace
{

void
checkMissingBuiltin(
    const ParameterTable& table,
    const std::string& canonicalName)
{
    const BuiltinParameter& entry =
        table.builtin(canonicalName);

    NTIC_EXPECT_TRUE(
        entry.info != nullptr);

    NTIC_EXPECT_EQ(
        entry.info->name(),
        canonicalName);

    NTIC_EXPECT_TRUE(
        !entry.hasValue);
}

void
checkPresentBuiltin(
    const ParameterTable& table,
    const std::string& canonicalName)
{
    const BuiltinParameter& entry =
        table.builtin(canonicalName);

    NTIC_EXPECT_TRUE(
        entry.info != nullptr);

    NTIC_EXPECT_EQ(
        entry.info->name(),
        canonicalName);

    NTIC_EXPECT_TRUE(
        entry.hasValue);
}

void
testEmptyData()
{
    const LBMParameterRegistry registry;

    const ParameterData data;

    ParameterTable table(
        registry,
        data);

    //----------------------------------------------------------
    // A complete static LBM table is created even when no
    // values were read from the configuration file.
    //----------------------------------------------------------

    checkMissingBuiltin(
        table,
        "case_name");

    checkMissingBuiltin(
        table,
        "density");

    checkMissingBuiltin(
        table,
        "kinematic_viscosity");

    checkMissingBuiltin(
        table,
        "mach_number");

    checkMissingBuiltin(
        table,
        "checkpoint_interval");

    NTIC_EXPECT_EQ(
        table.userDefined().size(),
        static_cast<std::size_t>(0));

    //----------------------------------------------------------
    // Unknown and alias names are not valid table keys.
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        table.builtin(
            "unknown_parameter"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        table.builtin(
            "rho"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        table.builtin(
            "Dynamic Viscosity"),
        std::runtime_error);
}

void
testDirectParameterDataImport()
{
    const LBMParameterRegistry registry;

    ParameterData data;

    data.builtinValues.emplace(
        "case_name",
        ParameterValue(
            std::string("direct_case")));

    data.builtinValues.emplace(
        "periodic_x",
        ParameterValue(true));

    data.builtinValues.emplace(
        "reynolds_number",
        ParameterValue(
            Quantity()));

    data.userDefinedValues.push_back(
        UserDefinedParameterData
        {
            "custom_quantity",
            ParameterValue(
                Quantity()),
            false,
            11
        });

    data.userDefinedValues.push_back(
        UserDefinedParameterData
        {
            "custom_bool",
            ParameterValue(false),
            true,
            12
        });

    data.userDefinedValues.push_back(
        UserDefinedParameterData
        {
            "custom_string",
            ParameterValue(
                std::string("abc")),
            false,
            13
        });

    ParameterTable table(
        registry,
        data);

    //----------------------------------------------------------
    // Built-in values.
    //----------------------------------------------------------

    checkPresentBuiltin(
        table,
        "case_name");

    checkPresentBuiltin(
        table,
        "periodic_x");

    checkPresentBuiltin(
        table,
        "reynolds_number");

    checkMissingBuiltin(
        table,
        "density");

    NTIC_EXPECT_TRUE(
        std::holds_alternative<std::string>(
            table.builtin(
                "case_name").value));

    NTIC_EXPECT_EQ(
        std::get<std::string>(
            table.builtin(
                "case_name").value),
        std::string("direct_case"));

    NTIC_EXPECT_TRUE(
        std::holds_alternative<bool>(
            table.builtin(
                "periodic_x").value));

    NTIC_EXPECT_TRUE(
        std::get<bool>(
            table.builtin(
                "periodic_x").value));

    NTIC_EXPECT_TRUE(
        std::holds_alternative<Quantity>(
            table.builtin(
                "reynolds_number").value));

    //----------------------------------------------------------
    // User-defined values preserve order and metadata.
    //----------------------------------------------------------

    const std::vector<UserParameter>& user =
        table.userDefined();

    NTIC_EXPECT_EQ(
        user.size(),
        static_cast<std::size_t>(3));

    NTIC_EXPECT_EQ(
        user[0].identifier,
        std::string("custom_quantity"));

    NTIC_EXPECT_TRUE(
        std::holds_alternative<Quantity>(
            user[0].value));

    NTIC_EXPECT_TRUE(
        !user[0].bound);

    NTIC_EXPECT_EQ(
        user[0].lineNumber,
        static_cast<std::size_t>(11));

    NTIC_EXPECT_EQ(
        user[1].identifier,
        std::string("custom_bool"));

    NTIC_EXPECT_TRUE(
        std::holds_alternative<bool>(
            user[1].value));

    NTIC_EXPECT_TRUE(
        !std::get<bool>(
            user[1].value));

    NTIC_EXPECT_TRUE(
        user[1].bound);

    NTIC_EXPECT_EQ(
        user[1].lineNumber,
        static_cast<std::size_t>(12));

    NTIC_EXPECT_EQ(
        user[2].identifier,
        std::string("custom_string"));

    NTIC_EXPECT_TRUE(
        std::holds_alternative<std::string>(
            user[2].value));

    NTIC_EXPECT_EQ(
        std::get<std::string>(
            user[2].value),
        std::string("abc"));

    NTIC_EXPECT_TRUE(
        !user[2].bound);

    NTIC_EXPECT_EQ(
        user[2].lineNumber,
        static_cast<std::size_t>(13));
}

void
testMutableAccess()
{
    const LBMParameterRegistry registry;

    const ParameterData data;

    ParameterTable table(
        registry,
        data);

    //----------------------------------------------------------
    // Validator-style modification of a built-in entry.
    //----------------------------------------------------------

    BuiltinParameter& mach =
        table.builtin(
            "mach_number");

    NTIC_EXPECT_TRUE(
        !mach.hasValue);

    mach.value =
        ParameterValue(
            Quantity());

    mach.hasValue =
        true;

    checkPresentBuiltin(
        table,
        "mach_number");

    NTIC_EXPECT_TRUE(
        std::holds_alternative<Quantity>(
            table.builtin(
                "mach_number").value));

    //----------------------------------------------------------
    // Validator-style replacement of an existing value.
    //----------------------------------------------------------

    BuiltinParameter& collision =
        table.builtin(
            "collision_model");

    collision.value =
        ParameterValue(
            std::string("cumulant"));

    collision.hasValue =
        true;

    NTIC_EXPECT_EQ(
        std::get<std::string>(
            table.builtin(
                "collision_model").value),
        std::string("cumulant"));

    collision.value =
        ParameterValue(
            std::string("mrt"));

    NTIC_EXPECT_EQ(
        std::get<std::string>(
            table.builtin(
                "collision_model").value),
        std::string("mrt"));

    //----------------------------------------------------------
    // Mutable user-defined access.
    //----------------------------------------------------------

    std::vector<UserParameter>& user =
        table.userDefined();

    user.push_back(
        UserParameter
        {
            "custom",
            ParameterValue(true),
            false,
            42
        });

    NTIC_EXPECT_EQ(
        table.userDefined().size(),
        static_cast<std::size_t>(1));

    user[0].bound =
        true;

    NTIC_EXPECT_TRUE(
        table.userDefined()[0].bound);
}

void
testConfigExportIntegration()
{
    Config config;

    config.parse(
        R"(
case name = table_case
rho = 998 kg/m^3
viscosity = 1.0e-3 Pa*s
Re = 1200
periodic x = true
collision model = cumulant
tau = 0.8
output directory = "results"

custom pressure = 101325 Pa
custom flag = false
custom model = enthalpy
)");

    //----------------------------------------------------------
    // Bind one user-defined value before export. The bound flag
    // must be preserved by ParameterData and ParameterTable.
    //----------------------------------------------------------

    bool customFlag = true;

    config.bindUserDefinedBool(
        customFlag,
        "custom flag");

    NTIC_EXPECT_TRUE(
        !customFlag);

    const ParameterData data =
        config.exportParameterData();

    const LBMParameterRegistry registry;

    ParameterTable table(
        registry,
        data);

    //----------------------------------------------------------
    // Built-in data imported under canonical names.
    //----------------------------------------------------------

    checkPresentBuiltin(
        table,
        "case_name");

    checkPresentBuiltin(
        table,
        "density");

    checkPresentBuiltin(
        table,
        "dynamic_viscosity");

    checkPresentBuiltin(
        table,
        "reynolds_number");

    checkPresentBuiltin(
        table,
        "periodic_x");

    checkPresentBuiltin(
        table,
        "collision_model");

    checkPresentBuiltin(
        table,
        "relaxation_time");

    checkPresentBuiltin(
        table,
        "output_directory");

    checkMissingBuiltin(
        table,
        "kinematic_viscosity");

    checkMissingBuiltin(
        table,
        "characteristic_velocity");

    checkMissingBuiltin(
        table,
        "time_step");

    checkMissingBuiltin(
        table,
        "mach_number");

    NTIC_EXPECT_EQ(
        std::get<std::string>(
            table.builtin(
                "case_name").value),
        std::string("table_case"));

    NTIC_EXPECT_TRUE(
        std::get<bool>(
            table.builtin(
                "periodic_x").value));

    NTIC_EXPECT_EQ(
        std::get<std::string>(
            table.builtin(
                "collision_model").value),
        std::string("cumulant"));

    NTIC_EXPECT_EQ(
        std::get<std::string>(
            table.builtin(
                "output_directory").value),
        std::string("results"));

    //----------------------------------------------------------
    // User-defined data imported in source order.
    //----------------------------------------------------------

    const std::vector<UserParameter>& user =
        table.userDefined();

    NTIC_EXPECT_EQ(
        user.size(),
        static_cast<std::size_t>(3));

    NTIC_EXPECT_EQ(
        user[0].identifier,
        std::string("custom_pressure"));

    NTIC_EXPECT_TRUE(
        std::holds_alternative<Quantity>(
            user[0].value));

    NTIC_EXPECT_TRUE(
        !user[0].bound);

    NTIC_EXPECT_EQ(
        user[1].identifier,
        std::string("custom_flag"));

    NTIC_EXPECT_TRUE(
        std::holds_alternative<bool>(
            user[1].value));

    NTIC_EXPECT_TRUE(
        user[1].bound);

    NTIC_EXPECT_EQ(
        user[2].identifier,
        std::string("custom_model"));

    NTIC_EXPECT_TRUE(
        std::holds_alternative<std::string>(
            user[2].value));

    NTIC_EXPECT_TRUE(
        !user[2].bound);
}

void
testInvalidBuiltinData()
{
    const LBMParameterRegistry registry;

    //----------------------------------------------------------
    // Unknown built-in key.
    //----------------------------------------------------------

    {
        ParameterData data;

        data.builtinValues.emplace(
            "unknown_parameter",
            ParameterValue(true));

        NTIC_EXPECT_THROW(
            ParameterTable(
                registry,
                data),
            std::runtime_error);
    }

    //----------------------------------------------------------
    // Alias key instead of canonical key.
    //----------------------------------------------------------

    {
        ParameterData data;

        data.builtinValues.emplace(
            "rho",
            ParameterValue(
                Quantity()));

        NTIC_EXPECT_THROW(
            ParameterTable(
                registry,
                data),
            std::runtime_error);
    }

    //----------------------------------------------------------
    // Another normalized alias key.
    //----------------------------------------------------------

    {
        ParameterData data;

        data.builtinValues.emplace(
            "viscosity",
            ParameterValue(
                Quantity()));

        NTIC_EXPECT_THROW(
            ParameterTable(
                registry,
                data),
            std::runtime_error);
    }
}

void
testInfoPointersReferenceRegistry()
{
    const LBMParameterRegistry registry;

    const ParameterData data;

    ParameterTable table(
        registry,
        data);

    const BuiltinParameter& density =
        table.builtin(
            "density");

    const ParameterInfo& registryInfo =
        registry.find(
            "density");

    NTIC_EXPECT_TRUE(
        density.info
        == &registryInfo);

    NTIC_EXPECT_EQ(
        density.info->category(),
        ParameterCategory::Physical);

    NTIC_EXPECT_EQ(
        density.info->dimension(),
        PhysicalDimension::Density());

    NTIC_EXPECT_TRUE(
        !density.info->description().empty());
}

} // namespace

int main()
{
    testEmptyData();
    testDirectParameterDataImport();
    testMutableAccess();
    testConfigExportIntegration();
    testInvalidBuiltinData();
    testInfoPointersReferenceRegistry();

    ntic::lbm::test::printSummary();

    return 0;
}

