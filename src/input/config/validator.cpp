#include "validator.hpp"

#include "physical_dimension.hpp"
#include "quantity.hpp"
#include "unit_expression.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>

namespace ntic::lbm::config
{

namespace
{

constexpr double latticeSoundSpeedSquared =
    1.0 / 3.0;

const double latticeSoundSpeed =
    1.0 / std::sqrt(3.0);

std::string
numberText(
    double value)
{
    std::ostringstream output;

    output.precision(16);

    output << value;

    return output.str();
}

} // namespace

//==============================================================
// Construction
//==============================================================

Validator::Validator(
    ParameterTable& table)
    :
    table_(table),
    modelAliasRegistry_(),
    warnings_()
{
}

//==============================================================
// Public interface
//==============================================================

void
Validator::validate()
{
    warnings_.clear();

    //----------------------------------------------------------
    // The groups are intentionally processed in dependency
    // order. Each group handles completeness, best practice,
    // over-definition and completion for its own parameters.
    //----------------------------------------------------------

    validateFluidProperties();

    validateCharacteristicScales();

    validateLatticeScaling();

    validateSimulationControl();

    validateModelSelection();

    validateRecommendations();
}

bool
Validator::hasWarnings() const
{
    return !warnings_.empty();
}

const std::vector<std::string>&
Validator::warnings() const
{
    return warnings_;
}

void
Validator::printWarnings(
    std::ostream& output) const
{
    for (const std::string& warning : warnings_)
    {
        output
            << "Warning: "
            << warning
            << '\n';
    }
}

//==============================================================
// Group 1: fluid properties
//
//     kinematic_viscosity
//         = dynamic_viscosity / density
//
// Any two quantities determine the third. The preferred input
// pair is density + dynamic_viscosity.
//==============================================================

void
Validator::validateFluidProperties()
{
    const bool hasDensity =
        has("density");

    const bool hasDynamic =
        has("dynamic_viscosity");

    const bool hasKinematic =
        has("kinematic_viscosity");

    const int count =
        static_cast<int>(hasDensity)
        + static_cast<int>(hasDynamic)
        + static_cast<int>(hasKinematic);

    if (count < 2)
    {
        throw std::runtime_error(
            "Fluid properties are incomplete. "
            "At least two of 'density', 'dynamic_viscosity' "
            "and 'kinematic_viscosity' are required.");
    }

    if (hasDensity)
    {
        requirePositive(
            quantitySI("density"),
            "density");
    }

    if (hasDynamic)
    {
        requirePositive(
            quantitySI("dynamic_viscosity"),
            "dynamic_viscosity");
    }

    if (hasKinematic)
    {
        requirePositive(
            quantitySI("kinematic_viscosity"),
            "kinematic_viscosity");
    }

    //----------------------------------------------------------
    // Preferred pair: density + dynamic_viscosity.
    //----------------------------------------------------------

    if (hasDensity && hasDynamic)
    {
        const double density =
            quantitySI("density");

        const double dynamicViscosity =
            quantitySI("dynamic_viscosity");

        const double derivedKinematic =
            dynamicViscosity / density;

        if (hasKinematic)
        {
            const double inputKinematic =
                quantitySI("kinematic_viscosity");

            if (!nearlyEqual(
                    inputKinematic,
                    derivedKinematic))
            {
                warn(
                    "'kinematic_viscosity' is inconsistent "
                    "with 'density' and 'dynamic_viscosity'. "
                    "The input value "
                    + numberText(inputKinematic)
                    + " has been replaced by "
                    + numberText(derivedKinematic)
                    + " in SI units.");
            }
            else
            {
                warn(
                    "'kinematic_viscosity' is redundant because "
                    "it is determined by 'density' and "
                    "'dynamic_viscosity'.");
            }
        }

        setQuantitySI(
            "kinematic_viscosity",
            derivedKinematic);

        return;
    }

    //----------------------------------------------------------
    // Non-preferred but valid pairs.
    //----------------------------------------------------------

    if (hasDensity && hasKinematic)
    {
        warn(
            "'dynamic_viscosity' is normally supplied together "
            "with 'density'; it has been derived from "
            "'density' and 'kinematic_viscosity'.");

        setQuantitySI(
            "dynamic_viscosity",
            quantitySI("density")
            * quantitySI("kinematic_viscosity"));

        return;
    }

    warn(
        "'density' is normally supplied together with "
        "'dynamic_viscosity'; it has been derived from "
        "'dynamic_viscosity' and 'kinematic_viscosity'.");

    setQuantitySI(
        "density",
        quantitySI("dynamic_viscosity")
        / quantitySI("kinematic_viscosity"));
}

//==============================================================
// Group 2: characteristic scales
//
//     reynolds_number
//         = characteristic_velocity
//           * characteristic_length
//           / kinematic_viscosity
//
// Any three quantities determine the fourth. The preferred
// input set is reynolds_number + characteristic_length +
// kinematic_viscosity, with characteristic_velocity derived.
//==============================================================

void
Validator::validateCharacteristicScales()
{
    const bool hasReynolds =
        has("reynolds_number");

    const bool hasLength =
        has("characteristic_length");

    const bool hasVelocity =
        has("characteristic_velocity");

    const bool hasKinematic =
        has("kinematic_viscosity");

    const int count =
        static_cast<int>(hasReynolds)
        + static_cast<int>(hasLength)
        + static_cast<int>(hasVelocity)
        + static_cast<int>(hasKinematic);

    if (count < 3)
    {
        throw std::runtime_error(
            "Characteristic scales are incomplete. "
            "At least three of 'reynolds_number', "
            "'characteristic_length', "
            "'characteristic_velocity' and "
            "'kinematic_viscosity' are required.");
    }

    if (hasReynolds)
    {
        requirePositive(
            quantitySI("reynolds_number"),
            "reynolds_number");
    }

    if (hasLength)
    {
        requirePositive(
            quantitySI("characteristic_length"),
            "characteristic_length");
    }

    if (hasVelocity)
    {
        requirePositive(
            quantitySI("characteristic_velocity"),
            "characteristic_velocity");
    }

    requirePositive(
        quantitySI("kinematic_viscosity"),
        "kinematic_viscosity");

    //----------------------------------------------------------
    // Preferred authoritative set:
    // Re + L + nu -> U
    //----------------------------------------------------------

    if (hasReynolds
        && hasLength
        && hasKinematic)
    {
        const double reynolds =
            quantitySI("reynolds_number");

        const double length =
            quantitySI("characteristic_length");

        const double kinematic =
            quantitySI("kinematic_viscosity");

        const double derivedVelocity =
            reynolds
            * kinematic
            / length;

        if (hasVelocity)
        {
            const double inputVelocity =
                quantitySI("characteristic_velocity");

            if (!nearlyEqual(
                    inputVelocity,
                    derivedVelocity))
            {
                warn(
                    "'characteristic_velocity' is inconsistent "
                    "with 'reynolds_number', "
                    "'characteristic_length' and "
                    "'kinematic_viscosity'. The input value "
                    + numberText(inputVelocity)
                    + " has been replaced by "
                    + numberText(derivedVelocity)
                    + " in SI units.");
            }
            else
            {
                warn(
                    "'characteristic_velocity' is redundant and "
                    "is normally derived from "
                    "'reynolds_number', "
                    "'characteristic_length' and "
                    "'kinematic_viscosity'.");
            }
        }

        setQuantitySI(
            "characteristic_velocity",
            derivedVelocity);

        return;
    }

    //----------------------------------------------------------
    // Other mathematically valid minimal sets.
    //----------------------------------------------------------

    warn(
        "The characteristic-scale input set is valid but is not "
        "the preferred set. The preferred inputs are "
        "'reynolds_number', 'characteristic_length' and "
        "'kinematic_viscosity'.");

    if (!hasReynolds)
    {
        setQuantitySI(
            "reynolds_number",
            quantitySI("characteristic_velocity")
            * quantitySI("characteristic_length")
            / quantitySI("kinematic_viscosity"));

        return;
    }

    if (!hasLength)
    {
        setQuantitySI(
            "characteristic_length",
            quantitySI("reynolds_number")
            * quantitySI("kinematic_viscosity")
            / quantitySI("characteristic_velocity"));

        return;
    }

    //----------------------------------------------------------
    // This branch is only reachable when the original table did
    // not contain kinematic_viscosity. Group 1 normally fills it
    // before this group, so keep this branch as a consistency
    // safeguard.
    //----------------------------------------------------------

    setQuantitySI(
        "kinematic_viscosity",
        quantitySI("characteristic_velocity")
        * quantitySI("characteristic_length")
        / quantitySI("reynolds_number"));
}

//==============================================================
// Group 3: lattice scaling
//
// Lattice units:
//     delta_x_lattice = 1
//     delta_t_lattice = 1
//
// Relations:
//     omega = 1 / tau
//     nu_lattice = (tau - 0.5) / 3
//     L_lattice = characteristic_length / grid_spacing
//     Re = u_lattice * L_lattice / nu_lattice
//     time_step = grid_spacing * u_lattice
//                 / characteristic_velocity
//     Mach = u_lattice / c_s, c_s = 1/sqrt(3)
//
// The preferred seed is relaxation_time. If it is absent,
// relaxation_frequency is preferred next. Directly supplied
// lattice velocity, Mach or physical time step are accepted,
// but generate a best-practice warning.
//==============================================================

void
Validator::validateLatticeScaling()
{
    if (!has("grid_spacing"))
    {
        throw std::runtime_error(
            "Lattice scaling requires "
            "'grid_spacing'.");
    }

    const double gridSpacing =
        quantitySI("grid_spacing");

    requirePositive(
        gridSpacing,
        "grid_spacing");

    const double length =
        quantitySI("characteristic_length");

    const double velocity =
        quantitySI("characteristic_velocity");

    const double reynolds =
        quantitySI("reynolds_number");

    requirePositive(
        length,
        "characteristic_length");

    requirePositive(
        velocity,
        "characteristic_velocity");

    requirePositive(
        reynolds,
        "reynolds_number");

    const double latticeLength =
        length / gridSpacing;

    requirePositive(
        latticeLength,
        "characteristic_length/grid_spacing");

    const bool hasTau =
        has("relaxation_time");

    const bool hasOmega =
        has("relaxation_frequency");

    const bool hasLatticeVelocity =
        has("lattice_characteristic_velocity");

    const bool hasMach =
        has("mach_number");

    const bool hasTimeStep =
        has("time_step");

    if (!hasTau
        && !hasOmega
        && !hasLatticeVelocity
        && !hasMach
        && !hasTimeStep)
    {
        throw std::runtime_error(
            "Lattice scaling is incomplete. Provide "
            "at least one of 'relaxation_time', "
            "'relaxation_frequency', "
            "'lattice_characteristic_velocity', "
            "'mach_number' or 'time_step'.");
    }

    double tau = 0.0;
    double omega = 0.0;
    double latticeVelocity = 0.0;

    //----------------------------------------------------------
    // Preferred authoritative seed: relaxation_time.
    //----------------------------------------------------------

    if (hasTau)
    {
        tau =
            quantitySI("relaxation_time");

        if (!(tau > 0.5))
        {
            throw std::runtime_error(
                "'relaxation_time' must be greater "
                "than 0.5.");
        }

        omega =
            1.0 / tau;

        const double latticeViscosity =
            latticeSoundSpeedSquared
            * (tau - 0.5);

        latticeVelocity =
            reynolds
            * latticeViscosity
            / latticeLength;

        if (hasOmega)
        {
            const double inputOmega =
                quantitySI("relaxation_frequency");

            if (!nearlyEqual(
                    inputOmega,
                    omega))
            {
                warn(
                    "'relaxation_frequency' is inconsistent "
                    "with 'relaxation_time'. The input value "
                    + numberText(inputOmega)
                    + " has been replaced by "
                    + numberText(omega)
                    + ".");
            }
            else
            {
                warn(
                    "'relaxation_frequency' is redundant because "
                    "it is determined by 'relaxation_time'.");
            }
        }
    }
    else if (hasOmega)
    {
        omega =
            quantitySI("relaxation_frequency");

        if (!(omega > 0.0 && omega < 2.0))
        {
            throw std::runtime_error(
                "'relaxation_frequency' must be "
                "strictly between 0 and 2.");
        }

        tau =
            1.0 / omega;

        const double latticeViscosity =
            latticeSoundSpeedSquared
            * (tau - 0.5);

        latticeVelocity =
            reynolds
            * latticeViscosity
            / latticeLength;

        warn(
            "'relaxation_frequency' is accepted, but "
            "'relaxation_time' is the preferred lattice-scaling "
            "input.");
    }
    else
    {
        //------------------------------------------------------
        // Non-preferred seeds: u_lattice, Mach or time_step.
        //------------------------------------------------------

        if (hasLatticeVelocity)
        {
            latticeVelocity =
                quantitySI(
                    "lattice_characteristic_velocity");

            warn(
                "'lattice_characteristic_velocity' is normally "
                "derived from the relaxation parameter and "
                "Reynolds similarity.");
        }
        else if (hasMach)
        {
            latticeVelocity =
                quantitySI("mach_number")
                * latticeSoundSpeed;

            warn(
                "'mach_number' is normally derived from "
                "'lattice_characteristic_velocity'.");
        }
        else
        {
            latticeVelocity =
                velocity
                * quantitySI("time_step")
                / gridSpacing;

            warn(
                "'time_step' is normally derived by the "
                "physical-to-lattice scaling.");
        }

        requirePositive(
            latticeVelocity,
            "lattice_characteristic_velocity");

        const double latticeViscosity =
            latticeVelocity
            * latticeLength
            / reynolds;

        tau =
            latticeViscosity
            / latticeSoundSpeedSquared
            + 0.5;

        if (!(tau > 0.5))
        {
            throw std::runtime_error(
                "The supplied lattice scaling results in "
                "'relaxation_time' not greater than 0.5.");
        }

        omega =
            1.0 / tau;
    }

    requirePositive(
        latticeVelocity,
        "lattice_characteristic_velocity");

    const double derivedMach =
        latticeVelocity
        / latticeSoundSpeed;

    const double derivedTimeStep =
        gridSpacing
        * latticeVelocity
        / velocity;

    //----------------------------------------------------------
    // Reconcile all redundant lattice quantities against the
    // selected authoritative scaling.
    //----------------------------------------------------------

    if (hasLatticeVelocity)
    {
        const double input =
            quantitySI(
                "lattice_characteristic_velocity");

        if (!nearlyEqual(
                input,
                latticeVelocity))
        {
            warn(
                "'lattice_characteristic_velocity' is "
                "inconsistent with the selected lattice scaling. "
                "The input value "
                + numberText(input)
                + " has been replaced by "
                + numberText(latticeVelocity)
                + ".");
        }
    }

    if (hasMach)
    {
        const double input =
            quantitySI("mach_number");

        if (!nearlyEqual(
                input,
                derivedMach))
        {
            warn(
                "'mach_number' is inconsistent with the selected "
                "lattice scaling. The input value "
                + numberText(input)
                + " has been replaced by "
                + numberText(derivedMach)
                + ".");
        }
    }

    if (hasTimeStep)
    {
        const double input =
            quantitySI("time_step");

        if (!nearlyEqual(
                input,
                derivedTimeStep))
        {
            warn(
                "'time_step' is inconsistent with the selected "
                "physical-to-lattice mapping. The input value "
                + numberText(input)
                + " has been replaced by "
                + numberText(derivedTimeStep)
                + " seconds.");
        }
    }

    setQuantitySI(
        "relaxation_time",
        tau);

    setQuantitySI(
        "relaxation_frequency",
        omega);

    setQuantitySI(
        "lattice_characteristic_velocity",
        latticeVelocity);

    setQuantitySI(
        "mach_number",
        derivedMach);

    setQuantitySI(
        "time_step",
        derivedTimeStep);
}

//==============================================================
// Group 4: simulation control
//
//     physical_time = max_steps * time_step
//
// restart_time is intentionally ignored by Validator.
// time_step has already been finalized by lattice scaling.
//==============================================================

void
Validator::validateSimulationControl()
{
    const double timeStep =
        quantitySI("time_step");

    requirePositive(
        timeStep,
        "time_step");

    const bool hasPhysicalTime =
        has("physical_time");

    const bool hasMaxSteps =
        has("max_steps");

    if (!hasPhysicalTime
        && !hasMaxSteps)
    {
        throw std::runtime_error(
            "Simulation control requires either "
            "'physical_time' or 'max_steps'.");
    }

    if (hasMaxSteps)
    {
        checkedPositiveInteger(
            quantitySI("max_steps"),
            "max_steps");
    }

    if (hasPhysicalTime)
    {
        requirePositive(
            quantitySI("physical_time"),
            "physical_time");
    }

    //----------------------------------------------------------
    // max_steps is authoritative when supplied because lattice
    // execution must use an integral number of steps.
    //----------------------------------------------------------

    if (hasMaxSteps)
    {
        const long long maxSteps =
            checkedPositiveInteger(
                quantitySI("max_steps"),
                "max_steps");

        const double derivedPhysicalTime =
            static_cast<double>(maxSteps)
            * timeStep;

        if (hasPhysicalTime)
        {
            const double inputPhysicalTime =
                quantitySI("physical_time");

            if (!nearlyEqual(
                    inputPhysicalTime,
                    derivedPhysicalTime))
            {
                warn(
                    "'physical_time' is inconsistent with "
                    "'max_steps' and the validated 'time_step'. "
                    "The input value "
                    + numberText(inputPhysicalTime)
                    + " has been replaced by "
                    + numberText(derivedPhysicalTime)
                    + " seconds.");
            }
            else
            {
                warn(
                    "'physical_time' is redundant because it is "
                    "determined by 'max_steps' and 'time_step'.");
            }
        }

        setQuantitySI(
            "physical_time",
            derivedPhysicalTime);

        return;
    }

    //----------------------------------------------------------
    // Only physical_time was supplied. Use enough complete
    // lattice steps to reach or exceed it.
    //----------------------------------------------------------

    const double requestedPhysicalTime =
        quantitySI("physical_time");

    const double rawSteps =
        requestedPhysicalTime
        / timeStep;

    const double maximumLongLong =
        static_cast<double>(
            std::numeric_limits<long long>::max());

    if (!std::isfinite(rawSteps)
        || rawSteps > maximumLongLong)
    {
        throw std::runtime_error(
            "The derived 'max_steps' is outside the "
            "supported integer range.");
    }

    const long long maxSteps =
        static_cast<long long>(
            std::ceil(rawSteps));

    if (maxSteps <= 0)
    {
        throw std::runtime_error(
            "The derived 'max_steps' must be positive.");
    }

    const double actualPhysicalTime =
        static_cast<double>(maxSteps)
        * timeStep;

    if (!nearlyEqual(
            requestedPhysicalTime,
            actualPhysicalTime))
    {
        warn(
            "'physical_time' is not an integral multiple of the "
            "validated 'time_step'. It has been increased from "
            + numberText(requestedPhysicalTime)
            + " to "
            + numberText(actualPhysicalTime)
            + " seconds to align with the lattice time step.");
    }

    setQuantitySI(
        "max_steps",
        static_cast<double>(maxSteps));

    setQuantitySI(
        "physical_time",
        actualPhysicalTime);
}

//==============================================================
// Group 5: model selection
//
// Required:
//     lattice_model
//     collision_model
//
// Optional:
//     wall_boundary
//     inlet_boundary
//     outlet_boundary
//
// Existing strings are validated and replaced by their unique
// canonical model spelling.
//==============================================================

void
Validator::validateModelSelection()
{
    if (!has("lattice_model"))
    {
        throw std::runtime_error(
            "'lattice_model' is required.");
    }

    if (!has("collision_model"))
    {
        throw std::runtime_error(
            "'collision_model' is required.");
    }

    setString(
        "lattice_model",
        modelAliasRegistry_.resolveLatticeModel(
            stringValue("lattice_model")));

    setString(
        "collision_model",
        modelAliasRegistry_.resolveCollisionModel(
            stringValue("collision_model")));

    if (has("wall_boundary"))
    {
        setString(
            "wall_boundary",
            modelAliasRegistry_.resolveWallBoundary(
                stringValue("wall_boundary")));
    }

    if (has("inlet_boundary"))
    {
        setString(
            "inlet_boundary",
            modelAliasRegistry_.resolveInletBoundary(
                stringValue("inlet_boundary")));
    }

    if (has("outlet_boundary"))
    {
        setString(
            "outlet_boundary",
            modelAliasRegistry_.resolveOutletBoundary(
                stringValue("outlet_boundary")));
    }
}

//==============================================================
// Final LBM recommendations
//==============================================================

void
Validator::validateRecommendations()
{
    const double tau =
        quantitySI("relaxation_time");

    const double omega =
        quantitySI("relaxation_frequency");

    const double mach =
        quantitySI("mach_number");

    if (mach > 0.1)
    {
        warn(
            "'mach_number' is "
            + numberText(mach)
            + ". A value not greater than 0.1 is generally "
              "recommended for weakly compressible LBM.");
    }

    if (tau < 0.55)
    {
        warn(
            "'relaxation_time' is "
            + numberText(tau)
            + ", which is close to the stability limit 0.5.");
    }
    else if (tau > 2.0)
    {
        warn(
            "'relaxation_time' is "
            + numberText(tau)
            + ". Values above 2.0 are generally not preferred.");
    }

    if (!(omega > 0.0 && omega < 2.0))
    {
        throw std::runtime_error(
            "The finalized 'relaxation_frequency' "
            "must be strictly between 0 and 2.");
    }
}

//==============================================================
// Helpers
//==============================================================

bool
Validator::has(
    const std::string& canonicalName) const
{
    return table_
        .builtin(canonicalName)
        .hasValue;
}

double
Validator::quantitySI(
    const std::string& canonicalName) const
{
    const BuiltinParameter& entry =
        table_.builtin(canonicalName);

    if (!entry.hasValue)
    {
        throw std::runtime_error(
            "Built-in parameter '"
            + canonicalName
            + "' has no value.");
    }

    if (!std::holds_alternative<Quantity>(
            entry.value))
    {
        throw std::runtime_error(
            "Built-in parameter '"
            + canonicalName
            + "' is not a quantity.");
    }

    const Quantity& quantity =
        std::get<Quantity>(
            entry.value);

    const PhysicalDimension& expectedDimension =
        entry.info->dimension();

    if (quantity.unit().dimension()
        != expectedDimension)
    {
        throw std::runtime_error(
            "Built-in quantity '"
            + canonicalName
            + "' has dimension "
            + quantity.unit().dimension().toString()
            + ", but "
            + expectedDimension.toString()
            + " is required.");
    }

    return quantity.value()
        * quantity.unit().scale();
}

std::string
Validator::stringValue(
    const std::string& canonicalName) const
{
    const BuiltinParameter& entry =
        table_.builtin(canonicalName);

    if (!entry.hasValue)
    {
        throw std::runtime_error(
            "Built-in parameter '"
            + canonicalName
            + "' has no value.");
    }

    if (!std::holds_alternative<std::string>(
            entry.value))
    {
        throw std::runtime_error(
            "Built-in parameter '"
            + canonicalName
            + "' is not a string.");
    }

    return std::get<std::string>(
        entry.value);
}

void
Validator::setQuantitySI(
    const std::string& canonicalName,
    double value)
{
    if (!std::isfinite(value))
    {
        throw std::runtime_error(
            "Derived value for '"
            + canonicalName
            + "' is not finite.");
    }

    BuiltinParameter& entry =
        table_.builtin(canonicalName);

    if (entry.info == nullptr)
    {
        throw std::runtime_error(
            "Built-in parameter '"
            + canonicalName
            + "' has no parameter metadata.");
    }

    entry.value =
        ParameterValue(
            Quantity(
                value,
                UnitExpression(
                    entry.info->dimension(),
                    1.0)));

    entry.hasValue =
        true;
}

void
Validator::setString(
    const std::string& canonicalName,
    const std::string& value)
{
    BuiltinParameter& entry =
        table_.builtin(canonicalName);

    entry.value =
        ParameterValue(value);

    entry.hasValue =
        true;
}

void
Validator::warn(
    const std::string& message)
{
    warnings_.push_back(message);
}

bool
Validator::nearlyEqual(
    double lhs,
    double rhs,
    double relativeTolerance,
    double absoluteTolerance)
{
    const double scale =
        std::max(
            std::abs(lhs),
            std::abs(rhs));

    return std::abs(lhs - rhs)
        <= absoluteTolerance
        + relativeTolerance * scale;
}

void
Validator::requirePositive(
    double value,
    const std::string& canonicalName)
{
    if (!std::isfinite(value)
        || value <= 0.0)
    {
        throw std::runtime_error(
            "'"
            + canonicalName
            + "' must be finite and greater than zero.");
    }
}

long long
Validator::checkedPositiveInteger(
    double value,
    const std::string& canonicalName)
{
    requirePositive(
        value,
        canonicalName);

    const double rounded =
        std::round(value);

    if (!nearlyEqual(
            value,
            rounded,
            0.0,
            1.0e-10))
    {
        throw std::runtime_error(
            "'"
            + canonicalName
            + "' must be an integer.");
    }

    if (rounded
        > static_cast<double>(
            std::numeric_limits<long long>::max()))
    {
        throw std::runtime_error(
            "'"
            + canonicalName
            + "' is outside the supported integer range.");
    }

    return static_cast<long long>(
        rounded);
}

} // namespace ntic::lbm::config
