#!/bin/bash

APP=./build/input/config/test_full_config
INPUT_DIR=tests/input_tests/config_tests/inputs

$APP $INPUT_DIR/01_minimal.nticonf
$APP $INPUT_DIR/02_standard.nticonf
$APP $INPUT_DIR/03_redundant_viscosity.nticonf
$APP $INPUT_DIR/04_derived_velocity.nticonf
$APP $INPUT_DIR/05_derived_time_step.nticonf
$APP $INPUT_DIR/06_lattice_velocity.nticonf
$APP $INPUT_DIR/07_relaxation_frequency.nticonf
$APP $INPUT_DIR/08_high_mach.nticonf
$APP $INPUT_DIR/09_large_tau.nticonf
$APP $INPUT_DIR/10_invalid_collision.nticonf
$APP $INPUT_DIR/11_invalid_lattice.nticonf
$APP $INPUT_DIR/12_missing_fluid.nticonf
$APP $INPUT_DIR/13_missing_reynolds.nticonf
$APP $INPUT_DIR/14_missing_scaling.nticonf
$APP $INPUT_DIR/15_negative_density.nticonf
$APP $INPUT_DIR/16_negative_viscosity.nticonf
$APP $INPUT_DIR/17_inconsistent_viscosity.nticonf
$APP $INPUT_DIR/18_inconsistent_velocity.nticonf
$APP $INPUT_DIR/19_inconsistent_relaxation.nticonf
$APP $INPUT_DIR/20_mach_input.nticonf
$APP $INPUT_DIR/21_time_step_input.nticonf
$APP $INPUT_DIR/22_redundant_simulation_control.nticonf
$APP $INPUT_DIR/23_inconsistent_simulation_control.nticonf
$APP $INPUT_DIR/24_model_aliases.nticonf
$APP $INPUT_DIR/25_user_defined.nticonf