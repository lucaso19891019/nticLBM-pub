#include "lbm_parameter_registry.hpp"

#include <stdexcept>
#include <string>

namespace ntic::lbm::config
{

//==============================================================
// Construction
//==============================================================

LBMParameterRegistry::LBMParameterRegistry()
    :
    parameterRegistry_(),
    aliasRegistry_(),
    valueTypes_()
{
    registerBuiltins();
}

//==============================================================
// Query
//==============================================================

bool
LBMParameterRegistry::contains(
    const std::string& name) const
{
    if (!aliasRegistry_.contains(name))
    {
        return false;
    }

    const std::string& canonicalName =
        aliasRegistry_.resolve(name);

    return parameterRegistry_.contains(
        canonicalName);
}

const std::string&
LBMParameterRegistry::resolve(
    const std::string& name) const
{
    return aliasRegistry_.resolve(name);
}

const ParameterInfo&
LBMParameterRegistry::find(
    const std::string& name) const
{
    const std::string& canonicalName =
        resolve(name);

    return parameterRegistry_.find(
        canonicalName);
}

ParameterValueType
LBMParameterRegistry::valueType(
    const std::string& name) const
{
    const std::string& canonicalName =
        resolve(name);

    const auto iter =
        valueTypes_.find(
            canonicalName);

    if (iter == valueTypes_.end())
    {
        throw std::runtime_error(
            "No value type is registered for LBM parameter '"
            + canonicalName
            + "'.");
    }

    return iter->second;
}

std::size_t
LBMParameterRegistry::size() const
{
    return parameterRegistry_.size();
}

//==============================================================
// Underlying registries
//==============================================================

const ParameterRegistry&
LBMParameterRegistry::parameterRegistry() const
{
    return parameterRegistry_;
}

const ParameterAliasRegistry&
LBMParameterRegistry::aliasRegistry() const
{
    return aliasRegistry_;
}

//==============================================================
// Internal registration
//==============================================================

void
LBMParameterRegistry::registerParameter(
    const std::string& canonicalName,
    const PhysicalDimension& dimension,
    ParameterCategory category,
    ParameterValueType valueType,
    const std::string& description,
    std::initializer_list<std::string> aliases)
{
    const std::string normalizedCanonicalName =
        ParameterAliasRegistry::normalize(
            canonicalName);

    if (normalizedCanonicalName
        != canonicalName)
    {
        throw std::invalid_argument(
            "LBM canonical parameter name '"
            + canonicalName
            + "' is not normalized. Expected '"
            + normalizedCanonicalName
            + "'.");
    }

    if (parameterRegistry_.contains(
            canonicalName))
    {
        throw std::runtime_error(
            "LBM parameter '"
            + canonicalName
            + "' has already been registered.");
    }

    if (valueTypes_.find(canonicalName)
        != valueTypes_.end())
    {
        throw std::runtime_error(
            "A value type is already registered for LBM parameter '"
            + canonicalName
            + "'.");
    }

    //----------------------------------------------------------
    // Register aliases first. registerAliases() also registers
    // the canonical name as an alias of itself.
    //----------------------------------------------------------

    aliasRegistry_.registerAliases(
        aliases,
        canonicalName);

    //----------------------------------------------------------
    // All built-in parameters are optional at parsing time.
    // Redundancy, missing data and derivability are handled by
    // the later validation stage.
    //----------------------------------------------------------

    parameterRegistry_.registerParameter(
        ParameterInfo(
            canonicalName,
            dimension,
            category,
            false,
            false,
            description));

    valueTypes_.emplace(
        canonicalName,
        valueType);
}

//==============================================================
// Built-in registration
//==============================================================

void
LBMParameterRegistry::registerBuiltins()
{
    const PhysicalDimension dimensionless =
        PhysicalDimension::Dimensionless();

    //----------------------------------------------------------
    // General
    //----------------------------------------------------------

    registerParameter(
        "case_name",
        dimensionless,
        ParameterCategory::General,
        ParameterValueType::String,
        "Name of the simulation case.",
        {
            "case",
            "case name",
            "casename"
        });

    registerParameter(
        "precision",
        dimensionless,
        ParameterCategory::General,
        ParameterValueType::String,
        "Floating-point precision used by the simulation.",
        {
            "floating point precision",
            "floating_point_precision",
            "floatingpointprecision"
        });

    registerParameter(
        "geometry_file",
        dimensionless,
        ParameterCategory::General,
        ParameterValueType::String,
        "Path to the geometry input file.",
        {
            "geometry",
            "geometry file",
            "geometryfile",
            "mesh",
            "mesh file",
            "mesh_file",
            "meshfile",
            "stl",
            "stl file",
            "stl_file",
            "stlfile"
        });

    registerParameter(
        "restart_file",
        dimensionless,
        ParameterCategory::General,
        ParameterValueType::String,
        "Path to the restart or checkpoint input file.",
        {
            "restart",
            "restart file",
            "restartfile",
            "checkpoint file",
            "checkpoint_file",
            "checkpointfile"
        });

    //----------------------------------------------------------
    // Physical
    //----------------------------------------------------------

    registerParameter(
        "density",
        PhysicalDimension::Density(),
        ParameterCategory::Physical,
        ParameterValueType::Quantity,
        "Reference fluid density.",
        {
            "rho",
            "reference density",
            "reference_density",
            "referencedensity"
        });

    registerParameter(
        "dynamic_viscosity",
        PhysicalDimension::DynamicViscosity(),
        ParameterCategory::Physical,
        ParameterValueType::Quantity,
        "Reference dynamic viscosity. The alias viscosity means dynamic viscosity.",
        {
            "viscosity",
            "dynamic viscosity",
            "dynamicviscosity",
            "mu"
        });

    registerParameter(
        "kinematic_viscosity",
        PhysicalDimension::KinematicViscosity(),
        ParameterCategory::Physical,
        ParameterValueType::Quantity,
        "Reference kinematic viscosity.",
        {
            "kinematic viscosity",
            "kinematicviscosity",
            "nu"
        });

    registerParameter(
        "characteristic_length",
        PhysicalDimension::Length(),
        ParameterCategory::Physical,
        ParameterValueType::Quantity,
        "Characteristic physical length used for nondimensional parameters.",
        {
            "characteristic length",
            "characteristiclength",
            "reference length",
            "reference_length",
            "referencelength"
        });

    registerParameter(
        "characteristic_velocity",
        PhysicalDimension::Velocity(),
        ParameterCategory::Physical,
        ParameterValueType::Quantity,
        "Characteristic physical velocity used for nondimensional parameters.",
        {
            "characteristic velocity",
            "characteristicvelocity",
            "reference velocity",
            "reference_velocity",
            "referencevelocity"
        });

    registerParameter(
        "reynolds_number",
        dimensionless,
        ParameterCategory::Physical,
        ParameterValueType::Quantity,
        "Reynolds number.",
        {
            "re",
            "reynolds",
            "reynolds number",
            "reynoldsnumber"
        });

    registerParameter(
        "physical_time",
        PhysicalDimension::Time(),
        ParameterCategory::Physical,
        ParameterValueType::Quantity,
        "Total physical simulation time.",
        {
            "physical time",
            "physicaltime",
            "total physical time",
            "total_physical_time",
            "totalphysicaltime",
            "end time",
            "end_time",
            "endtime"
        });

    registerParameter(
        "restart_time",
        PhysicalDimension::Time(),
        ParameterCategory::Physical,
        ParameterValueType::Quantity,
        "Physical time represented by the restart state.",
        {
            "restart time",
            "restarttime",
            "checkpoint time",
            "checkpoint_time",
            "checkpointtime"
        });

    registerParameter(
        "time_step",
        PhysicalDimension::Time(),
        ParameterCategory::Physical,
        ParameterValueType::Quantity,
        "Physical time represented by one simulation step.",
        {
            "time step",
            "timestep",
            "physical time step",
            "physical_time_step",
            "physicaltimestep",
            "dt"
        });

    registerParameter(
        "grid_spacing",
        PhysicalDimension::Length(),
        ParameterCategory::Physical,
        ParameterValueType::Quantity,
        "Minimum physical spacing of the finest grid.",
        {
            "grid spacing",
            "gridspacing",
            "minimum grid spacing",
            "minimum_grid_spacing",
            "minimumgridspacing",
            "minimum lattice spacing",
            "minimum_lattice_spacing",
            "minimumlatticespacing",
            "lattice spacing",
            "lattice_spacing",
            "latticespacing",
            "grid resolution",
            "grid_resolution",
            "gridresolution",
            "dx"
        });

    registerParameter(
        "periodic_x",
        dimensionless,
        ParameterCategory::Physical,
        ParameterValueType::Bool,
        "Whether the physical domain is periodic in the x direction.",
        {
            "periodic x",
            "periodicx",
            "x periodic",
            "x_periodic",
            "xperiodic"
        });

    registerParameter(
        "periodic_y",
        dimensionless,
        ParameterCategory::Physical,
        ParameterValueType::Bool,
        "Whether the physical domain is periodic in the y direction.",
        {
            "periodic y",
            "periodicy",
            "y periodic",
            "y_periodic",
            "yperiodic"
        });

    registerParameter(
        "periodic_z",
        dimensionless,
        ParameterCategory::Physical,
        ParameterValueType::Bool,
        "Whether the physical domain is periodic in the z direction.",
        {
            "periodic z",
            "periodicz",
            "z periodic",
            "z_periodic",
            "zperiodic"
        });

    registerParameter(
        "symmetric_x",
        dimensionless,
        ParameterCategory::Physical,
        ParameterValueType::Bool,
        "Whether the physical domain uses symmetry in the x direction.",
        {
            "symmetric x",
            "symmetricx",
            "symmetry x",
            "symmetry_x",
            "symmetryx",
            "x symmetric",
            "x_symmetric",
            "xsymmetric",
            "x symmetry",
            "x_symmetry",
            "xsymmetry"
        });

    registerParameter(
        "symmetric_y",
        dimensionless,
        ParameterCategory::Physical,
        ParameterValueType::Bool,
        "Whether the physical domain uses symmetry in the y direction.",
        {
            "symmetric y",
            "symmetricy",
            "symmetry y",
            "symmetry_y",
            "symmetryy",
            "y symmetric",
            "y_symmetric",
            "ysymmetric",
            "y symmetry",
            "y_symmetry",
            "ysymmetry"
        });

    registerParameter(
        "symmetric_z",
        dimensionless,
        ParameterCategory::Physical,
        ParameterValueType::Bool,
        "Whether the physical domain uses symmetry in the z direction.",
        {
            "symmetric z",
            "symmetricz",
            "symmetry z",
            "symmetry_z",
            "symmetryz",
            "z symmetric",
            "z_symmetric",
            "zsymmetric",
            "z symmetry",
            "z_symmetry",
            "zsymmetry"
        });

    //----------------------------------------------------------
    // Lattice
    //----------------------------------------------------------

    registerParameter(
        "lattice_model",
        dimensionless,
        ParameterCategory::Lattice,
        ParameterValueType::String,
        "Discrete lattice velocity model, such as D2Q9, D3Q19 or D3Q27.",
        {
            "lattice",
            "lattice model",
            "latticemodel",
            "velocity set",
            "velocity_set",
            "velocityset"
        });

    registerParameter(
        "collision_model",
        dimensionless,
        ParameterCategory::Lattice,
        ParameterValueType::String,
        "LBM collision model.",
        {
            "collision",
            "collision model",
            "collisionmodel",
            "collision operator",
            "collision_operator",
            "collisionoperator"
        });

    registerParameter(
        "relaxation_time",
        dimensionless,
        ParameterCategory::Lattice,
        ParameterValueType::Quantity,
        "Dimensionless lattice relaxation time.",
        {
            "relaxation time",
            "relaxationtime",
            "lattice relaxation time",
            "lattice_relaxation_time",
            "latticerelaxationtime",
            "tau"
        });

    registerParameter(
        "relaxation_frequency",
        dimensionless,
        ParameterCategory::Lattice,
        ParameterValueType::Quantity,
        "Dimensionless lattice relaxation frequency.",
        {
            "relaxation frequency",
            "relaxationfrequency",
            "lattice relaxation frequency",
            "lattice_relaxation_frequency",
            "latticerelaxationfrequency",
            "omega"
        });

    registerParameter(
        "lattice_characteristic_velocity",
        dimensionless,
        ParameterCategory::Lattice,
        ParameterValueType::Quantity,
        "Characteristic velocity in lattice units.",
        {
            "lattice characteristic velocity",
            "latticecharacteristicvelocity",
            "lattice velocity",
            "lattice_velocity",
            "latticevelocity",
            "u lat",
            "u_lat",
            "ulat"
        });

    registerParameter(
        "mach_number",
        dimensionless,
        ParameterCategory::Lattice,
        ParameterValueType::Quantity,
        "Lattice Mach number.",
        {
            "mach",
            "mach number",
            "machnumber",
            "ma"
        });

    registerParameter(
        "number_of_levels",
        dimensionless,
        ParameterCategory::Lattice,
        ParameterValueType::Quantity,
        "Number of grid or lattice refinement levels.",
        {
            "number of levels",
            "numberoflevels",
            "num levels",
            "num_levels",
            "numlevels",
            "levels",
            "nlevels"
        });

    registerParameter(
        "wall_boundary",
        dimensionless,
        ParameterCategory::Lattice,
        ParameterValueType::String,
        "Wall boundary treatment used throughout the case.",
        {
            "wall",
            "wall boundary",
            "wallboundary",
            "wall boundary model",
            "wall_boundary_model",
            "wallboundarymodel"
        });

    registerParameter(
        "inlet_boundary",
        dimensionless,
        ParameterCategory::Lattice,
        ParameterValueType::String,
        "Inlet boundary treatment used throughout the case.",
        {
            "inlet",
            "inlet boundary",
            "inletboundary",
            "inlet boundary model",
            "inlet_boundary_model",
            "inletboundarymodel"
        });

    registerParameter(
        "outlet_boundary",
        dimensionless,
        ParameterCategory::Lattice,
        ParameterValueType::String,
        "Outlet boundary treatment used throughout the case.",
        {
            "outlet",
            "outlet boundary",
            "outletboundary",
            "outlet boundary model",
            "outlet_boundary_model",
            "outletboundarymodel"
        });

    registerParameter(
        "max_steps",
        dimensionless,
        ParameterCategory::Lattice,
        ParameterValueType::Quantity,
        "Maximum number of lattice time steps.",
        {
            "max steps",
            "maxsteps",
            "maximum steps",
            "maximum_steps",
            "maximumsteps",
            "number of steps",
            "number_of_steps",
            "numberofsteps",
            "num steps",
            "num_steps",
            "numsteps",
            "nsteps"
        });

    //----------------------------------------------------------
    // Output
    //----------------------------------------------------------

    registerParameter(
        "output_directory",
        dimensionless,
        ParameterCategory::Output,
        ParameterValueType::String,
        "Directory used for simulation output.",
        {
            "output directory",
            "outputdirectory",
            "output dir",
            "output_dir",
            "outputdir",
            "result directory",
            "result_directory",
            "resultdirectory"
        });

    registerParameter(
        "output_interval",
        PhysicalDimension::Time(),
        ParameterCategory::Output,
        ParameterValueType::Quantity,
        "Physical-time interval between regular output operations.",
        {
            "output interval",
            "outputinterval",
            "output time interval",
            "output_time_interval",
            "outputtimeinterval",
            "write interval",
            "write_interval",
            "writeinterval"
        });

    registerParameter(
        "checkpoint_interval",
        PhysicalDimension::Time(),
        ParameterCategory::Output,
        ParameterValueType::Quantity,
        "Physical-time interval between checkpoint operations.",
        {
            "checkpoint interval",
            "checkpointinterval",
            "checkpoint time interval",
            "checkpoint_time_interval",
            "checkpointtimeinterval",
            "restart interval",
            "restart_interval",
            "restartinterval"
        });
}

} // namespace ntic::lbm::config
