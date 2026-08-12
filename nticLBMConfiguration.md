# nticLBM Configuration

The nticLBM configuration module provides a human-readable input format for defining the physical problem, lattice scaling, LBM models, boundary treatments, simulation duration, and output options used by an nticLBM case.

Configuration files use the `.nticonf` extension.

The configuration system is designed so that users can specify a simulation primarily in **physical units**, while nticLBM checks dimensional consistency and derives related physical and lattice quantities whenever possible.

A typical configuration may look like:

```text
case name = standard

density = 998 kg/m^3
dynamic viscosity = 0.001 Pa*s

re = 200
characteristic length = 0.05 m

grid spacing = 0.0005 m
relaxation time = 0.8

physical time = 0.02 s

lattice model = D3Q19
collision model = Cumulant

wall boundary = IBB
inlet boundary = ZH-vel
outlet boundary = Default
```

## Configuration syntax

Each configuration entry is written as:

```text
parameter = value
```

Whitespace around parameter names, values, and the `=` sign is allowed.

Parameter names are intended to be readable rather than restricted to programming-style identifiers. For example:

```text
characteristic length = 0.05 m
dynamic viscosity = 0.001 Pa*s
relaxation time = 0.8
```

Common alternative forms and aliases are also accepted. For example:

```text
re = 200
rho = 1000 kg/m^3
mu = 0.001 Pa*s
nu = 1e-6 m^2/s
tau = 0.8
omega = 1.25
```

These correspond to the canonical parameters:

```text
reynolds_number
density
dynamic_viscosity
kinematic_viscosity
relaxation_time
relaxation_frequency
```

Parameter aliases are normalized before interpretation, allowing commonly used names with spaces, underscores, or compact spelling.

For example, forms such as:

```text
characteristic length
characteristic_length
characteristiclength
```

refer to the same physical quantity.

## Comments

The `#` character starts a comment outside quoted strings.

For example:

```text
# Water at approximately room conditions
density = 998 kg/m^3

re = 200       # Reynolds number
```

A `#` character contained inside a quoted string is treated as part of the string rather than as the beginning of a comment.

## Strings

Text parameters may be entered as ordinary values or quoted when needed.

Quoted strings are useful for values containing whitespace or characters that would otherwise have special meaning in the configuration file.

Typical string-valued parameters include:

```text
case name
geometry file
restart file
output directory
lattice model
collision model
wall boundary
inlet boundary
outlet boundary
```

## Physical quantities and units

Physical parameters may be specified together with units:

```text
density = 998 kg/m^3
dynamic viscosity = 0.001 Pa*s
characteristic length = 50 mm
grid spacing = 0.5 mm
physical time = 20 ms
```

nticLBM interprets the numerical value together with its unit and converts the result to a consistent SI representation internally.

The configuration system supports dimensional quantities constructed from SI units, prefixes, multiplication, division, and integer powers.

Examples include:

```text
m
mm
s
ms

m/s
m^2/s

kg/m^3
Pa*s
```

This allows users to describe a physical problem using convenient engineering-scale units without manually converting every value before writing the input file.

### Unit expressions

Compound units can be constructed using multiplication and division:

```text
kg/m^3
m/s
m^2/s
Pa*s
```

Integer powers can be used where appropriate:

```text
m^2
m^3
s^-1
```

Grouping is supported in unit expressions when a more complex expression needs to be written explicitly.

## Dimensional validation

Every registered physical parameter has an expected physical dimension.

The configuration system checks that the supplied unit is compatible with the parameter.

For example:

```text
characteristic length = 0.05 m
```

is valid because the parameter requires a length.

A value with an incompatible dimension is rejected rather than silently interpreted.

This prevents errors such as accidentally supplying a time unit to a length parameter or using a dynamic-viscosity unit where kinematic viscosity is expected.

Dimensionless parameters do not require a physical unit.

For example:

```text
re = 200
relaxation time = 0.8
```

## Fluid properties

The basic fluid-property relationship used by the configuration system is

```text
kinematic_viscosity = dynamic_viscosity / density
```

or

$$
\nu = \frac{\mu}{\rho}.
$$

The three associated parameters are:

```text
density
dynamic_viscosity
kinematic_viscosity
```

Any two are sufficient to determine the third.

### Preferred input

The preferred form is:

```text
density = 998 kg/m^3
dynamic viscosity = 0.001 Pa*s
```

nticLBM then derives the kinematic viscosity automatically.

### Alternative input

The following is also accepted:

```text
density = 998 kg/m^3
kinematic viscosity = 1.002e-6 m^2/s
```

In this case, dynamic viscosity can be derived.

Likewise:

```text
dynamic viscosity = 0.001 Pa*s
kinematic viscosity = 1.002e-6 m^2/s
```

allows density to be derived.

### Redundant input

Users may supply all three values:

```text
density = 998 kg/m^3
dynamic viscosity = 0.001 Pa*s
kinematic viscosity = 1.002e-6 m^2/s
```

nticLBM checks whether the supplied values are mutually consistent.

If a redundant value agrees with the value derived from the preferred quantities, a warning may be issued to indicate that the input is unnecessary.

If it is inconsistent, nticLBM reports the inconsistency and uses the value determined from the authoritative physical quantities.

All density and viscosity values must be positive.

## Reynolds-number specification

The characteristic physical scales satisfy

$$
Re = \frac{U L}{\nu},
$$

where

- \(Re\) is the Reynolds number,
- \(U\) is the characteristic physical velocity,
- \(L\) is the characteristic physical length,
- \(\nu\) is the kinematic viscosity.

The corresponding configuration parameters are:

```text
reynolds_number
characteristic_velocity
characteristic_length
kinematic_viscosity
```

Any three determine the fourth.

### Preferred specification

The preferred input is:

```text
re = 200
characteristic length = 0.05 m
```

together with sufficient fluid properties to determine kinematic viscosity.

For example:

```text
density = 998 kg/m^3
dynamic viscosity = 0.001 Pa*s

re = 200
characteristic length = 0.05 m
```

nticLBM then derives the characteristic physical velocity.

This allows users to define a flow using the quantities that are commonly known before the simulation starts: fluid properties, characteristic length, and Reynolds number.

### Explicit characteristic velocity

A characteristic velocity may also be supplied:

```text
characteristic velocity = 0.004 m/s
```

If Reynolds number, length, viscosity, and velocity are all present, nticLBM checks whether they satisfy the Reynolds-number relationship.

A redundant consistent value is accepted with a warning.

An inconsistent value is reported and reconciled according to the preferred parameter relationship.

## Grid spacing

The physical lattice resolution is specified using:

```text
grid spacing = 0.0005 m
```

Aliases such as:

```text
lattice spacing
minimum lattice spacing
grid resolution
```

are accepted.

`grid_spacing` represents the physical spacing of the finest lattice level.

The ratio

$$
L_\mathrm{lattice} = \frac{L_\mathrm{physical}}{\Delta x}
$$

determines the characteristic length measured in lattice cells.

For example:

```text
characteristic length = 0.05 m
grid spacing = 0.0005 m
```

corresponds to a characteristic lattice length of 100 cells.

Both quantities must be positive.

## Physical-to-lattice scaling

nticLBM supports automatic conversion between physical flow scales and lattice quantities.

The principal lattice-scaling quantities are:

```text
relaxation_time
relaxation_frequency
lattice_characteristic_velocity
mach_number
time_step
```

At least one lattice-scaling quantity must be supplied.

The preferred quantity is:

```text
relaxation time
```

### Relaxation time

For example:

```text
relaxation time = 0.8
```

The relaxation time is dimensionless and must satisfy

$$
\tau > 0.5.
$$

From the physical Reynolds number, physical/grid length ratio, and relaxation time, nticLBM can derive the lattice viscosity, lattice characteristic velocity, physical time step, relaxation frequency, and lattice Mach number.

For the standard isothermal lattice relation,

$$
\nu_\mathrm{lat} = \frac{\tau-0.5}{3}.
$$

### Relaxation frequency

The relaxation frequency is related to relaxation time by

$$
\omega = \frac{1}{\tau}.
$$

It may be supplied directly:

```text
relaxation frequency = 1.25
```

or:

```text
omega = 1.25
```

The accepted range is

$$
0 < \omega < 2.
$$

Supplying relaxation frequency is valid, although relaxation time is the preferred lattice-scaling input.

If both `relaxation_time` and `relaxation_frequency` are supplied, their consistency is checked.

### Lattice characteristic velocity

Users may explicitly specify the characteristic velocity in lattice units:

```text
lattice velocity = 0.05
```

This is dimensionless.

When possible, nticLBM normally derives this quantity from the physical Reynolds similarity and the selected relaxation parameter.

### Mach number

The lattice Mach number may also be used as a scaling input:

```text
mach = 0.1
```

For the standard lattice sound speed

$$
c_s = \frac{1}{\sqrt{3}},
$$

the characteristic lattice velocity is related to Mach number by

$$
Ma = \frac{u_\mathrm{lat}}{c_s}.
$$

Mach number is normally derived from the lattice characteristic velocity, but it can also be supplied when the user wants to control the scaling directly.

### Physical time step

A physical time step may be supplied directly:

```text
time step = 1e-5 s
```

or:

```text
dt = 1e-5 s
```

The relationship between physical and lattice velocity is

$$
u_\mathrm{lat} = U_\mathrm{physical} \frac{\Delta t}{\Delta x}.
$$

Normally the time step is derived from the selected lattice scaling. Direct specification is nevertheless supported.

## Automatic consistency checking

Several physical and lattice quantities describe the same underlying scaling relations. nticLBM therefore allows users to provide redundant information, but checks it for consistency.

Examples include:

```text
density
dynamic viscosity
kinematic viscosity
```

```text
reynolds number
characteristic length
characteristic velocity
kinematic viscosity
```

and:

```text
relaxation time
relaxation frequency
lattice characteristic velocity
mach number
time step
```

When redundant information is supplied, the configuration system distinguishes between:

- a valid but redundant value,
- a valid alternative way of defining the same problem,
- an inconsistent value,
- and a physically invalid value.

This permits compact input files while still allowing advanced users to explicitly provide additional scaling quantities for verification.

## Simulation time

The requested physical duration of a simulation is specified using:

```text
physical time = 0.02 s
```

Alternative names such as:

```text
total physical time
end time
```

are accepted.

Once the physical time step is known, nticLBM determines the required number of lattice steps.

Because a simulation executes an integer number of lattice steps, the requested physical duration may not always correspond exactly to an integer step count. The final step count is chosen so that the complete requested interval is represented.

A maximum step count may also be specified directly:

```text
max steps = 10000
```

## Restart time

Restarted simulations may specify the physical time represented by the loaded state:

```text
restart time = 0.5 s
```

Aliases such as:

```text
checkpoint time
```

are also accepted.

A restart/checkpoint file may be specified separately through the corresponding file parameter.

## Lattice models

The discrete lattice velocity model is selected using:

```text
lattice model = D3Q19
```

Currently recognized lattice models are:

```text
D2Q9
D3Q15
D3Q19
D3Q27
```

Common parameter aliases such as:

```text
lattice
velocity set
```

may also be used for the parameter name.

Unsupported lattice model names are rejected during configuration validation.

## Collision models

The LBM collision operator is selected using:

```text
collision model = Cumulant
```

Currently recognized collision models include:

```text
BGK
MRT
TRT
Cascaded
Cumulant
```

The Cascaded model also accepts common aliases associated with the central-moment formulation, including forms such as:

```text
central moment
central-moment
CM
```

These are normalized to the canonical `Cascaded` model name.

Unsupported collision-model names are rejected.

## Wall boundary model

The wall treatment is selected with:

```text
wall boundary = IBB
```

Supported wall treatments currently include:

```text
InterpolatedBB
IBM
```

`InterpolatedBB` accepts aliases such as:

```text
IBB
interpolated bb
interpolated bounce back
```

`IBM` also accepts forms such as:

```text
immersed
immersed boundary
```

The selected value is normalized to the corresponding canonical model name.

## Inlet and outlet boundary models

Inlet and outlet treatments are configured independently:

```text
inlet boundary = ZH-vel
outlet boundary = Default
```

Currently recognized flow-boundary models include:

```text
Default
ZouHeVelocity
ZouHePressure
```

Convenient aliases are supported.

For example, the Zou-He velocity condition may be written as:

```text
ZouHeVelocity
Zou-He velocity
Zou He velocity
ZH velocity
ZH vel
ZH-vel
```

Similarly, the Zou-He pressure condition accepts forms such as:

```text
ZouHePressure
Zou-He pressure
Zou He pressure
ZH pressure
ZH pre
ZH-pre
```

These aliases are normalized to the canonical model names.

## Symmetry flags

The configuration system supports symmetry flags for physical-domain directions.

For example:

```text
symmetric x = true
symmetric y = false
symmetric z = false
```

Equivalent parameter-name aliases using forms such as `symmetry_x`, `x_symmetric`, or `x symmetry` are also recognized.

These parameters are Boolean values.

## Multiple lattice levels

The number of lattice or grid-refinement levels can be specified with:

```text
number of levels = 1
```

Common abbreviated forms such as:

```text
num levels
levels
nlevels
```

are recognized.

The parameter is dimensionless.

## General case information

A configuration may contain general case information such as the case name and paths associated with input, restart, and output data.

For example:

```text
case name = cylinder_Re200
restart file = restart/checkpoint.dat
output directory = results
```

These entries are treated as text rather than physical quantities.

## Parameter aliases

The configuration language intentionally supports readable parameter names and commonly used CFD/LBM notation.

Examples include:

```text
density
rho
```

```text
dynamic viscosity
viscosity
mu
```

```text
kinematic viscosity
nu
```

```text
reynolds number
reynolds
re
```

```text
relaxation time
tau
```

```text
relaxation frequency
omega
```

```text
mach number
mach
ma
```

The canonical parameter name is used internally after parsing, so different aliases do not create different parameters.

Users are encouraged to use descriptive names in production input files and short mathematical aliases when they improve readability.

## Recommended input style

Although several equivalent ways of defining a case are supported, the recommended physical specification is:

```text
density = ...
dynamic viscosity = ...

reynolds number = ...
characteristic length = ...

grid spacing = ...
relaxation time = ...
```

With this set of inputs, nticLBM can derive:

```text
kinematic viscosity
characteristic velocity
relaxation frequency
lattice characteristic velocity
mach number
physical time step
```

This approach keeps the user input tied to physically meaningful quantities while leaving lattice-unit conversion to the configuration system.

A complete example is:

```text
#--------------------------------------------------
# Case
#--------------------------------------------------

case name = cylinder_Re200

#--------------------------------------------------
# Fluid
#--------------------------------------------------

density = 998 kg/m^3
dynamic viscosity = 0.001 Pa*s

#--------------------------------------------------
# Physical scales
#--------------------------------------------------

reynolds number = 200
characteristic length = 0.05 m

#--------------------------------------------------
# Lattice scaling
#--------------------------------------------------

grid spacing = 0.0005 m
relaxation time = 0.8

#--------------------------------------------------
# Simulation time
#--------------------------------------------------

physical time = 0.02 s

#--------------------------------------------------
# LBM
#--------------------------------------------------

lattice model = D3Q19
collision model = Cumulant

#--------------------------------------------------
# Boundaries
#--------------------------------------------------

wall boundary = IBB
inlet boundary = ZH-vel
outlet boundary = Default

#--------------------------------------------------
# Output
#--------------------------------------------------

output directory = results
```

## Errors and warnings

The configuration system distinguishes between conditions that make a case invalid and conditions that are valid but potentially undesirable or redundant.

An error is produced for conditions such as:

- missing physical information required to define the fluid,
- missing quantities required to establish the Reynolds scaling,
- missing lattice-scaling information,
- incompatible physical dimensions,
- non-positive density, viscosity, characteristic length, or grid spacing,
- relaxation time not greater than `0.5`,
- relaxation frequency outside the interval `(0, 2)`,
- unsupported lattice models,
- unsupported collision models,
- unsupported boundary models,
- or values that cannot produce a valid physical-to-lattice scaling.

Warnings are used when the configuration remains usable but contains redundant, non-preferred, or inconsistent information.

Examples include supplying both `tau` and `omega`, explicitly supplying a quantity that would normally be derived, or supplying a redundant viscosity or characteristic velocity.

This behavior allows nticLBM to detect common setup mistakes before a simulation begins while still supporting several practical ways of specifying the same physical problem.

## Minimal information required for a physical LBM case

For the current physical-to-lattice scaling workflow, a typical case requires enough information to determine all of the following:

1. Fluid properties.
2. Reynolds scaling.
3. Characteristic physical length.
4. Grid spacing.
5. At least one lattice-scaling parameter.
6. Lattice and collision models as required by the simulation.
7. Simulation duration or step control as required by the application.

A compact recommended setup is therefore:

```text
density = 998 kg/m^3
dynamic viscosity = 0.001 Pa*s

re = 200
characteristic length = 0.05 m

grid spacing = 0.0005 m
tau = 0.8

physical time = 0.02 s

lattice model = D3Q19
collision model = Cumulant
```

The remaining dependent physical and lattice quantities can then be determined automatically by nticLBM.
