#include "lbm_parameter_registry.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>
#include <vector>

using namespace ntic::lbm::config;

namespace
{

struct ExpectedParameter
{
    const char* name;
    ParameterCategory category;
    ParameterValueType type;
    PhysicalDimension dimension;
};

struct ExpectedAlias
{
    const char* alias;
    const char* canonical;
};

void
checkParameter(
    const LBMParameterRegistry& registry,
    const ExpectedParameter& expected)
{
    NTIC_EXPECT_TRUE(
        registry.contains(expected.name));

    const ParameterInfo& info =
        registry.find(expected.name);

    NTIC_EXPECT_EQ(
        info.name(),
        std::string(expected.name));

    NTIC_EXPECT_EQ(
        info.category(),
        expected.category);

    NTIC_EXPECT_EQ(
        registry.valueType(expected.name),
        expected.type);

    NTIC_EXPECT_EQ(
        info.dimension(),
        expected.dimension);

    NTIC_EXPECT_TRUE(
        !info.required());

    NTIC_EXPECT_TRUE(
        !info.derived());

    NTIC_EXPECT_TRUE(
        !info.description().empty());

    NTIC_EXPECT_EQ(
        registry.resolve(expected.name),
        std::string(expected.name));
}

void
checkAlias(
    const LBMParameterRegistry& registry,
    const ExpectedAlias& expected)
{
    NTIC_EXPECT_TRUE(
        registry.contains(expected.alias));

    NTIC_EXPECT_EQ(
        registry.resolve(expected.alias),
        std::string(expected.canonical));

    NTIC_EXPECT_EQ(
        registry.find(expected.alias).name(),
        std::string(expected.canonical));

    NTIC_EXPECT_EQ(
        registry.valueType(expected.alias),
        registry.valueType(expected.canonical));
}

} // namespace

int main()
{
    const LBMParameterRegistry registry;

    const PhysicalDimension dimensionless =
        PhysicalDimension::Dimensionless();

    const std::vector<ExpectedParameter> parameters =
    {
        {"case_name", ParameterCategory::General,
         ParameterValueType::String, dimensionless},

        {"precision", ParameterCategory::General,
         ParameterValueType::String, dimensionless},

        {"geometry_file", ParameterCategory::General,
         ParameterValueType::String, dimensionless},

        {"restart_file", ParameterCategory::General,
         ParameterValueType::String, dimensionless},

        {"density", ParameterCategory::Physical,
         ParameterValueType::Quantity,
         PhysicalDimension::Density()},

        {"dynamic_viscosity", ParameterCategory::Physical,
         ParameterValueType::Quantity,
         PhysicalDimension::DynamicViscosity()},

        {"kinematic_viscosity", ParameterCategory::Physical,
         ParameterValueType::Quantity,
         PhysicalDimension::KinematicViscosity()},

        {"characteristic_length", ParameterCategory::Physical,
         ParameterValueType::Quantity,
         PhysicalDimension::Length()},

        {"characteristic_velocity", ParameterCategory::Physical,
         ParameterValueType::Quantity,
         PhysicalDimension::Velocity()},

        {"reynolds_number", ParameterCategory::Physical,
         ParameterValueType::Quantity, dimensionless},

        {"physical_time", ParameterCategory::Physical,
         ParameterValueType::Quantity,
         PhysicalDimension::Time()},

        {"restart_time", ParameterCategory::Physical,
         ParameterValueType::Quantity,
         PhysicalDimension::Time()},

        {"time_step", ParameterCategory::Physical,
         ParameterValueType::Quantity,
         PhysicalDimension::Time()},

        {"grid_spacing", ParameterCategory::Physical,
         ParameterValueType::Quantity,
         PhysicalDimension::Length()},

        {"periodic_x", ParameterCategory::Physical,
         ParameterValueType::Bool, dimensionless},

        {"periodic_y", ParameterCategory::Physical,
         ParameterValueType::Bool, dimensionless},

        {"periodic_z", ParameterCategory::Physical,
         ParameterValueType::Bool, dimensionless},

        {"symmetric_x", ParameterCategory::Physical,
         ParameterValueType::Bool, dimensionless},

        {"symmetric_y", ParameterCategory::Physical,
         ParameterValueType::Bool, dimensionless},

        {"symmetric_z", ParameterCategory::Physical,
         ParameterValueType::Bool, dimensionless},

        {"lattice_model", ParameterCategory::Lattice,
         ParameterValueType::String, dimensionless},

        {"collision_model", ParameterCategory::Lattice,
         ParameterValueType::String, dimensionless},

        {"relaxation_time", ParameterCategory::Lattice,
         ParameterValueType::Quantity, dimensionless},

        {"relaxation_frequency", ParameterCategory::Lattice,
         ParameterValueType::Quantity, dimensionless},

	{"lattice_characteristic_velocity", ParameterCategory::Lattice,
    	 ParameterValueType::Quantity, dimensionless},

        {"mach_number", ParameterCategory::Lattice,
         ParameterValueType::Quantity, dimensionless},

        {"number_of_levels", ParameterCategory::Lattice,
         ParameterValueType::Quantity, dimensionless},

        {"wall_boundary", ParameterCategory::Lattice,
         ParameterValueType::String, dimensionless},

        {"inlet_boundary", ParameterCategory::Lattice,
         ParameterValueType::String, dimensionless},

        {"outlet_boundary", ParameterCategory::Lattice,
         ParameterValueType::String, dimensionless},

        {"max_steps", ParameterCategory::Lattice,
         ParameterValueType::Quantity, dimensionless},

        {"output_directory", ParameterCategory::Output,
         ParameterValueType::String, dimensionless},

        {"output_interval", ParameterCategory::Output,
         ParameterValueType::Quantity,
         PhysicalDimension::Time()},

        {"checkpoint_interval", ParameterCategory::Output,
         ParameterValueType::Quantity,
         PhysicalDimension::Time()}
    };

    NTIC_EXPECT_EQ(
        registry.size(),
        parameters.size());

    NTIC_EXPECT_EQ(
        registry.parameterRegistry().size(),
        parameters.size());

    NTIC_EXPECT_TRUE(
        registry.aliasRegistry().size()
        >= parameters.size());

    for (const ExpectedParameter& parameter : parameters)
    {
        checkParameter(
            registry,
            parameter);
    }

    const std::vector<ExpectedAlias> aliases =
    {
        {"Case Name", "case_name"},
        {"CASE", "case_name"},
        {"Floating-Point Precision", "precision"},
        {"mesh file", "geometry_file"},
        {"STL_FILE", "geometry_file"},
        {"checkpoint file", "restart_file"},

        {"rho", "density"},
        {"reference density", "density"},
        {"viscosity", "dynamic_viscosity"},
        {"Dynamic Viscosity", "dynamic_viscosity"},
        {"MU", "dynamic_viscosity"},
        {"nu", "kinematic_viscosity"},
        {"Kinematic-Viscosity", "kinematic_viscosity"},
        {"reference length", "characteristic_length"},
        {"reference velocity", "characteristic_velocity"},
        {"Re", "reynolds_number"},
        {"Reynolds Number", "reynolds_number"},
        {"total physical time", "physical_time"},
        {"end_time", "physical_time"},
        {"checkpoint time", "restart_time"},
        {"dt", "time_step"},
        {"physical time step", "time_step"},
        {"dx", "grid_spacing"},
        {"minimum grid spacing", "grid_spacing"},
        {"grid resolution", "grid_spacing"},
        {"x periodic", "periodic_x"},
        {"PERIODIC Y", "periodic_y"},
        {"z-periodic", "periodic_z"},
        {"symmetry x", "symmetric_x"},
        {"Y Symmetric", "symmetric_y"},
        {"z_symmetry", "symmetric_z"},

        {"velocity set", "lattice_model"},
        {"collision operator", "collision_model"},
        {"tau", "relaxation_time"},
        {"omega", "relaxation_frequency"},
        {"Ma", "mach_number"},
        {"levels", "number_of_levels"},
        {"wall boundary model", "wall_boundary"},
        {"inlet", "inlet_boundary"},
        {"outlet", "outlet_boundary"},
        {"nsteps", "max_steps"},

        {"output dir", "output_directory"},
        {"write interval", "output_interval"},
        {"restart interval", "checkpoint_interval"}
    };

    for (const ExpectedAlias& alias : aliases)
    {
        checkAlias(
            registry,
            alias);
    }

    //----------------------------------------------------------
    // Unknown but syntactically valid names
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(
        !registry.contains("pressure"));

    NTIC_EXPECT_TRUE(
        !registry.contains("temperature"));

    NTIC_EXPECT_TRUE(
        !registry.contains("custom_parameter"));

    NTIC_EXPECT_THROW(
        registry.resolve("custom_parameter"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        registry.find("custom_parameter"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        registry.valueType("custom_parameter"),
        std::runtime_error);

    //----------------------------------------------------------
    // Invalid names are rejected by the alias layer
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
        registry.contains("dynamic/viscosity"),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        registry.resolve("dynamic.viscosity"),
        std::invalid_argument);

    NTIC_EXPECT_THROW(
        registry.find("dynamic+viscosity"),
        std::invalid_argument);

    ntic::lbm::test::printSummary();

    return 0;
}
