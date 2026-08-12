#!/bin/bash

APP=./build/input/test_full_config

$APP tests/input_tests/inputs/01_minimal.nticonf
$APP tests/input_tests/inputs/02_standard.nticonf
$APP tests/input_tests/inputs/03_redundant_viscosity.nticonf
$APP tests/input_tests/inputs/04_derived_velocity.nticonf
$APP tests/input_tests/inputs/05_derived_time_step.nticonf
$APP tests/input_tests/inputs/06_lattice_velocity.nticonf
$APP tests/input_tests/inputs/07_relaxation_frequency.nticonf
$APP tests/input_tests/inputs/08_high_mach.nticonf
$APP tests/input_tests/inputs/09_large_tau.nticonf
$APP tests/input_tests/inputs/10_invalid_collision.nticonf
$APP tests/input_tests/inputs/11_invalid_lattice.nticonf
$APP tests/input_tests/inputs/12_missing_fluid.nticonf
$APP tests/input_tests/inputs/13_missing_reynolds.nticonf
$APP tests/input_tests/inputs/14_missing_scaling.nticonf
$APP tests/input_tests/inputs/15_negative_density.nticonf
$APP tests/input_tests/inputs/16_negative_viscosity.nticonf
$APP tests/input_tests/inputs/17_inconsistent_viscosity.nticonf
$APP tests/input_tests/inputs/18_inconsistent_velocity.nticonf
$APP tests/input_tests/inputs/19_inconsistent_relaxation.nticonf
$APP tests/input_tests/inputs/20_mach_input.nticonf
$APP tests/input_tests/inputs/21_time_step_input.nticonf
$APP tests/input_tests/inputs/22_redundant_simulation_control.nticonf
$APP tests/input_tests/inputs/23_inconsistent_simulation_control.nticonf
$APP tests/input_tests/inputs/24_model_aliases.nticonf
$APP tests/input_tests/inputs/25_user_defined.nticonf
