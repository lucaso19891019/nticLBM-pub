#include "parameter_registry.hpp"
#include "test_framework.hpp"

#include <stdexcept>

using ntic::lbm::config::ParameterCategory;
using ntic::lbm::config::ParameterInfo;
using ntic::lbm::config::ParameterRegistry;
using ntic::lbm::config::PhysicalDimension;

int main()
{
    //----------------------------------------------------------
    // Empty Registry
    //----------------------------------------------------------

    ParameterRegistry registry;

    NTIC_EXPECT_TRUE(
        registry.size() == 0);

    NTIC_EXPECT_TRUE(
        !registry.contains("density"));

    //----------------------------------------------------------
    // Register One Parameter
    //----------------------------------------------------------

    registry.registerParameter(
        ParameterInfo(
            "density",
            PhysicalDimension::Density(),
            ParameterCategory::Physical,
            false,
            false,
            "Fluid density"));

    NTIC_EXPECT_TRUE(
        registry.size() == 1);

    NTIC_EXPECT_TRUE(
        registry.contains("density"));

    //----------------------------------------------------------
    // Find
    //----------------------------------------------------------

    const ParameterInfo& density =
        registry.find("density");

    NTIC_EXPECT_TRUE(
        density.name() == "density");

    NTIC_EXPECT_EQ(
        density.dimension(),
        PhysicalDimension::Density());

    NTIC_EXPECT_TRUE(
        density.category()
        == ParameterCategory::Physical);

    NTIC_EXPECT_TRUE(
        !density.required());

    NTIC_EXPECT_TRUE(
        !density.derived());

    NTIC_EXPECT_TRUE(
        density.description()
        == "Fluid density");

    //----------------------------------------------------------
    // Register Another Parameter
    //----------------------------------------------------------

    registry.registerParameter(
        ParameterInfo(
            "Re",
            PhysicalDimension::Dimensionless(),
            ParameterCategory::Physical,
            true,
            false,
            "Reynolds number"));

    NTIC_EXPECT_TRUE(
        registry.size() == 2);

    NTIC_EXPECT_TRUE(
        registry.contains("Re"));

    //----------------------------------------------------------
    // Find Reynolds Number
    //----------------------------------------------------------

    const ParameterInfo& re =
        registry.find("Re");

    NTIC_EXPECT_EQ(
        re.dimension(),
        PhysicalDimension::Dimensionless());

    NTIC_EXPECT_TRUE(
        re.required());

    NTIC_EXPECT_TRUE(
        !re.derived());

    //----------------------------------------------------------
    // Duplicate Registration
    //----------------------------------------------------------

    std::size_t oldSize = registry.size();

    NTIC_EXPECT_THROW(

        registry.registerParameter(

            ParameterInfo(
                "density",
                PhysicalDimension::Density(),
                ParameterCategory::Physical,
                false,
                false)),

        std::runtime_error);

    // Registry 应保持不变
    NTIC_EXPECT_TRUE(
        registry.size() == oldSize);

    NTIC_EXPECT_TRUE(
        registry.contains("density"));

    NTIC_EXPECT_TRUE(
        registry.contains("Re"));

    //----------------------------------------------------------
    // Unknown Parameter
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(

        registry.find("velocity"),

        std::runtime_error);

    //----------------------------------------------------------
    // contains()
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(
        registry.contains("density"));

    NTIC_EXPECT_TRUE(
        registry.contains("Re"));

    NTIC_EXPECT_TRUE(
        !registry.contains("pressure"));

    //----------------------------------------------------------
    // Registry Integrity
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(
        registry.size() == 2);

    NTIC_EXPECT_TRUE(
        registry.find("density").name() == "density");

    NTIC_EXPECT_TRUE(
        registry.find("Re").name() == "Re");

    //----------------------------------------------------------
    // ParameterInfo::accepts()
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(
        density.accepts(
            PhysicalDimension::Density()));

    NTIC_EXPECT_TRUE(
        !density.accepts(
            PhysicalDimension::Length()));

    NTIC_EXPECT_TRUE(
        !density.accepts(
            PhysicalDimension::Velocity()));

    NTIC_EXPECT_TRUE(
        re.accepts(
            PhysicalDimension::Dimensionless()));

    NTIC_EXPECT_TRUE(
        !re.accepts(
            PhysicalDimension::Density()));

    //----------------------------------------------------------

    ntic::lbm::test::printSummary();

    return 0;
}
