#include "parameter_alias_registry.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>

using ntic::lbm::config::ParameterAliasRegistry;

int main()
{
    NTIC_EXPECT_EQ(
        ParameterAliasRegistry::normalize("Viscosity"),
        std::string("viscosity"));

    NTIC_EXPECT_EQ(
        ParameterAliasRegistry::normalize("Dynamic Viscosity"),
        std::string("dynamic_viscosity"));

    NTIC_EXPECT_EQ(
        ParameterAliasRegistry::normalize("dynamic-viscosity"),
        std::string("dynamic_viscosity"));

    NTIC_EXPECT_EQ(
        ParameterAliasRegistry::normalize("dynamic__ viscosity"),
        std::string("dynamic_viscosity"));

    NTIC_EXPECT_EQ(
        ParameterAliasRegistry::normalize("  dynamic viscosity  "),
        std::string("dynamic_viscosity"));

    NTIC_EXPECT_EQ(
        ParameterAliasRegistry::normalize("D3Q19"),
        std::string("d3q19"));

    NTIC_EXPECT_THROW(
        ParameterAliasRegistry::normalize(""),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        ParameterAliasRegistry::normalize("123parameter"),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        ParameterAliasRegistry::normalize("dynamic/viscosity"),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        ParameterAliasRegistry::normalize("dynamic.viscosity"),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        ParameterAliasRegistry::normalize("dynamic+viscosity"),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        ParameterAliasRegistry::normalize("dynamic*viscosity"),
        std::invalid_argument);

    ParameterAliasRegistry registry;

    NTIC_EXPECT_EQ(
        registry.size(),
        static_cast<std::size_t>(0));

    registry.registerAlias(
        "viscosity",
        "dynamic_viscosity");

    NTIC_EXPECT_TRUE(registry.contains("Viscosity"));

    NTIC_EXPECT_EQ(
        registry.resolve("VISCOSITY"),
        std::string("dynamic_viscosity"));

    registry.registerAlias(
        "VISCOSITY",
        "Dynamic Viscosity");

    NTIC_EXPECT_EQ(
        registry.size(),
        static_cast<std::size_t>(1));

    NTIC_EXPECT_THROW(
        registry.registerAlias(
            "viscosity",
            "kinematic_viscosity"),
        std::runtime_error);

    registry.registerAliases(
        {
            "dynamic viscosity",
            "dynamic-viscosity",
            "dynamic_viscosity",
            "mu",
            "MU"
        },
        "dynamic_viscosity");

    NTIC_EXPECT_EQ(
        registry.resolve("Dynamic Viscosity"),
        std::string("dynamic_viscosity"));

    NTIC_EXPECT_EQ(
        registry.resolve("mu"),
        std::string("dynamic_viscosity"));

    NTIC_EXPECT_EQ(
        registry.resolve("dynamic_viscosity"),
        std::string("dynamic_viscosity"));

    registry.registerAliases(
        {
            "nu",
            "kinematic viscosity",
            "kinematic-viscosity"
        },
        "kinematic_viscosity");

    NTIC_EXPECT_EQ(
        registry.resolve("NU"),
        std::string("kinematic_viscosity"));

    registry.registerAliases(
        {
            "Re",
            "Reynolds",
            "Reynolds Number"
        },
        "reynolds_number");

    NTIC_EXPECT_EQ(
        registry.resolve("re"),
        std::string("reynolds_number"));

    NTIC_EXPECT_EQ(
        registry.resolve("Reynolds-Number"),
        std::string("reynolds_number"));

    const std::size_t sizeBeforeConflict = registry.size();

    NTIC_EXPECT_THROW(
        registry.registerAliases(
            {
                "new_alias_before_conflict",
                "viscosity",
                "new_alias_after_conflict"
            },
            "another_parameter"),
        std::runtime_error);

    NTIC_EXPECT_EQ(
        registry.size(),
        sizeBeforeConflict);

    NTIC_EXPECT_TRUE(
        !registry.contains("new_alias_before_conflict"));

    NTIC_EXPECT_TRUE(
        !registry.contains("new_alias_after_conflict"));

    NTIC_EXPECT_THROW(
        registry.registerAlias(
            "dynamic/viscosity",
            "dynamic_viscosity"),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        registry.contains("dynamic/viscosity"),
        std::invalid_argument);

    NTIC_EXPECT_TRUE(
        !registry.contains("unknown_parameter"));

    NTIC_EXPECT_THROW(
        registry.resolve("unknown_parameter"),
        std::runtime_error);

    ntic::lbm::test::printSummary();
    return 0;
}
