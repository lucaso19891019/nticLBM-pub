cmake_minimum_required(VERSION 3.18)

project(nticLBM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(
        -Wall
        -Wextra
        -Wpedantic
    )
endif()

#====================================================
# Test output directories
#====================================================

set(INPUT_CONFIG_TEST_OUTPUT_DIR
    "${CMAKE_BINARY_DIR}/input/config"
)

set(INPUT_STL_TEST_OUTPUT_DIR
    "${CMAKE_BINARY_DIR}/input/stl"
)

function(set_input_config_test_output TARGET_NAME)
    set_target_properties(${TARGET_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${INPUT_CONFIG_TEST_OUTPUT_DIR}"
    )
endfunction()

function(set_input_stl_test_output TARGET_NAME)
    set_target_properties(${TARGET_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${INPUT_STL_TEST_OUTPUT_DIR}"
    )
endfunction()

#====================================================
# test_rational
#====================================================

add_executable(test_rational
    tests/input_tests/config_tests/test_rational.cpp
    src/common/rational.cpp
)

target_include_directories(test_rational PRIVATE
    src/common
    tests/input_tests/config_tests
)

set_input_config_test_output(test_rational)

#====================================================
# test_physical_dimension
#====================================================

add_executable(test_physical_dimension
    tests/input_tests/config_tests/test_physical_dimension.cpp
    src/common/rational.cpp
    src/input/config/physical_dimension.cpp
)

target_include_directories(test_physical_dimension PRIVATE
    src/common
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_physical_dimension)

#====================================================
# test_parameter_registry
#====================================================

add_executable(test_parameter_registry
    tests/input_tests/config_tests/test_parameter_registry.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp
    src/input/config/parameter_info.cpp
    src/input/config/parameter_registry.cpp
)

target_include_directories(test_parameter_registry PRIVATE
    src/common
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_parameter_registry)

#====================================================
# test_unit
#====================================================

add_executable(test_unit
    tests/input_tests/config_tests/test_unit.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp
    src/input/config/unit.cpp
    src/input/config/unit_registry.cpp
)

target_include_directories(test_unit PRIVATE
    src/common
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_unit)

#====================================================
# test_prefix_registry
#====================================================

add_executable(test_prefix_registry
    tests/input_tests/config_tests/test_prefix_registry.cpp

    src/input/config/prefix.cpp
    src/input/config/prefix_registry.cpp
)

target_include_directories(test_prefix_registry PRIVATE
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_prefix_registry)

#====================================================
# test_unit_expression
#====================================================

add_executable(test_unit_expression
    tests/input_tests/config_tests/test_unit_expression.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp

    src/input/config/unit.cpp
    src/input/config/unit_expression.cpp
)

target_include_directories(test_unit_expression PRIVATE
    src/common
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_unit_expression)

#====================================================
# test_quantity
#====================================================

add_executable(test_quantity
    tests/input_tests/config_tests/test_quantity.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp
    src/input/config/unit.cpp
    src/input/config/unit_expression.cpp
    src/input/config/quantity.cpp
)

target_include_directories(test_quantity PRIVATE
    src/common
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_quantity)

#====================================================
# test_parameter
#====================================================

add_executable(test_parameter

    tests/input_tests/config_tests/test_parameter.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp
    src/input/config/parameter_info.cpp
    src/input/config/quantity.cpp
    src/input/config/unit.cpp
    src/input/config/unit_expression.cpp
    src/input/config/parameter.cpp
)

target_include_directories(test_parameter PRIVATE
    src/common
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_parameter)

#====================================================
# test_lexer
#====================================================

add_executable(test_lexer
    tests/input_tests/config_tests/test_lexer.cpp

    src/input/config/token.cpp
    src/input/config/lexer.cpp
)

target_include_directories(test_lexer PRIVATE
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_lexer)

#====================================================
# test_token_stream
#====================================================

add_executable(test_token_stream
    tests/input_tests/config_tests/test_token_stream.cpp

    src/input/config/token.cpp
    src/input/config/lexer.cpp
    src/input/config/token_stream.cpp
)

target_include_directories(test_token_stream PRIVATE
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_token_stream)

#====================================================
# test_unit_parser
#====================================================

add_executable(test_unit_parser
    tests/input_tests/config_tests/test_unit_parser.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp

    src/input/config/unit.cpp
    src/input/config/unit_registry.cpp

    src/input/config/prefix.cpp
    src/input/config/prefix_registry.cpp

    src/input/config/unit_expression.cpp

    src/input/config/token.cpp
    src/input/config/lexer.cpp
    src/input/config/token_stream.cpp

    src/input/config/unit_parser.cpp
)

target_include_directories(test_unit_parser PRIVATE
    src/common
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_unit_parser)

#====================================================
# test_quantity_parser
#====================================================

add_executable(test_quantity_parser

    tests/input_tests/config_tests/test_quantity_parser.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp

    src/input/config/unit.cpp
    src/input/config/unit_registry.cpp

    src/input/config/prefix.cpp
    src/input/config/prefix_registry.cpp

    src/input/config/unit_expression.cpp
    src/input/config/quantity.cpp

    src/input/config/token.cpp
    src/input/config/lexer.cpp
    src/input/config/token_stream.cpp

    src/input/config/unit_parser.cpp
    src/input/config/quantity_parser.cpp
)

target_include_directories(test_quantity_parser PRIVATE

    src/common

    src/input/config

    tests/input_tests/config_tests
)

set_input_config_test_output(test_quantity_parser)

#====================================================
# test_parameter_parser
#====================================================

add_executable(test_parameter_parser

    tests/input_tests/config_tests/test_parameter_parser.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp

    src/input/config/unit.cpp
    src/input/config/unit_registry.cpp

    src/input/config/prefix.cpp
    src/input/config/prefix_registry.cpp

    src/input/config/unit_expression.cpp
    src/input/config/quantity.cpp

    src/input/config/token.cpp
    src/input/config/lexer.cpp
    src/input/config/token_stream.cpp

    src/input/config/unit_parser.cpp
    src/input/config/quantity_parser.cpp

    src/input/config/parameter_parser.cpp
)

target_include_directories(test_parameter_parser PRIVATE

    src/common

    src/input/config

    tests/input_tests/config_tests
)

set_input_config_test_output(test_parameter_parser)

#====================================================
# test_parameter_alias_registry
#====================================================

add_executable(test_parameter_alias_registry

    tests/input_tests/config_tests/test_parameter_alias_registry.cpp

    src/input/config/parameter_alias_registry.cpp
)

target_include_directories(
    test_parameter_alias_registry
PRIVATE

    src/input/config

    tests/input_tests/config_tests
)

set_input_config_test_output(test_parameter_alias_registry)

#====================================================
# test_lbm_parameter_registry
#====================================================

add_executable(test_lbm_parameter_registry

    tests/input_tests/config_tests/test_lbm_parameter_registry.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp
    src/input/config/parameter_info.cpp
    src/input/config/parameter_registry.cpp
    src/input/config/parameter_alias_registry.cpp
    src/input/config/lbm_parameter_registry.cpp
)

target_include_directories(
    test_lbm_parameter_registry
PRIVATE

    src/common
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_lbm_parameter_registry)

#====================================================
# test_config
#====================================================

add_executable(test_config

    tests/input_tests/config_tests/test_config.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp

    src/input/config/prefix.cpp
    src/input/config/prefix_registry.cpp

    src/input/config/unit.cpp
    src/input/config/unit_expression.cpp
    src/input/config/unit_registry.cpp

    src/input/config/quantity.cpp

    src/input/config/token.cpp
    src/input/config/lexer.cpp
    src/input/config/token_stream.cpp

    src/input/config/unit_parser.cpp
    src/input/config/quantity_parser.cpp
    src/input/config/parameter_parser.cpp

    src/input/config/parameter_info.cpp
    src/input/config/parameter_registry.cpp
    src/input/config/parameter_alias_registry.cpp
    src/input/config/lbm_parameter_registry.cpp

    src/input/config/config.cpp
)

target_include_directories(
    test_config
PRIVATE

    src/common
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_config)

#====================================================
# test_parameter_table
#====================================================

add_executable(test_parameter_table

    tests/input_tests/config_tests/test_parameter_table.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp

    src/input/config/prefix.cpp
    src/input/config/prefix_registry.cpp

    src/input/config/unit.cpp
    src/input/config/unit_expression.cpp
    src/input/config/unit_registry.cpp

    src/input/config/quantity.cpp

    src/input/config/token.cpp
    src/input/config/lexer.cpp
    src/input/config/token_stream.cpp

    src/input/config/unit_parser.cpp
    src/input/config/quantity_parser.cpp
    src/input/config/parameter_parser.cpp

    src/input/config/parameter_info.cpp
    src/input/config/parameter_registry.cpp
    src/input/config/parameter_alias_registry.cpp
    src/input/config/lbm_parameter_registry.cpp

    src/input/config/config.cpp
    src/input/config/parameter_table.cpp
)

target_include_directories(
    test_parameter_table
PRIVATE

    src/common
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_parameter_table)

#====================================================
# test_model_alias_registry
#====================================================

add_executable(test_model_alias_registry
    tests/input_tests/config_tests/test_model_alias_registry.cpp
    src/input/config/parameter_alias_registry.cpp
    src/input/config/model_alias_registry.cpp
)

target_include_directories(test_model_alias_registry PRIVATE
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_model_alias_registry)

#====================================================
# test_si_unit_formatter
#====================================================

add_executable(test_si_unit_formatter

    tests/input_tests/config_tests/test_si_unit_formatter.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp
    src/input/config/si_unit_formatter.cpp
)

target_include_directories(test_si_unit_formatter PRIVATE

    src/common
    src/input/config
)

set_input_config_test_output(test_si_unit_formatter)

#====================================================
# test_full_config
#====================================================

add_executable(test_full_config

    tests/input_tests/config_tests/test_full_config.cpp

    src/common/rational.cpp

    src/input/config/physical_dimension.cpp

    src/input/config/prefix.cpp
    src/input/config/prefix_registry.cpp

    src/input/config/unit.cpp
    src/input/config/unit_expression.cpp
    src/input/config/unit_registry.cpp

    src/input/config/quantity.cpp

    src/input/config/token.cpp
    src/input/config/lexer.cpp
    src/input/config/token_stream.cpp

    src/input/config/unit_parser.cpp
    src/input/config/quantity_parser.cpp
    src/input/config/parameter_parser.cpp

    src/input/config/parameter_info.cpp
    src/input/config/parameter_registry.cpp
    src/input/config/parameter_alias_registry.cpp
    src/input/config/lbm_parameter_registry.cpp

    src/input/config/config.cpp

    src/input/config/parameter_table.cpp

    # 新增
    src/input/config/model_alias_registry.cpp
    src/input/config/validator.cpp
    src/input/config/si_unit_formatter.cpp
)

target_include_directories(
    test_full_config
PRIVATE

    src/common
    src/input/config
    tests/input_tests/config_tests
)

set_input_config_test_output(test_full_config)

#====================================================
# test_stl_validator
#====================================================

add_executable(test_stl_validator

    tests/input_tests/stl_tests/test_stl_validator.cpp

    src/input/stl/stl_reader.cpp
    src/input/stl/stl_validator.cpp
)

target_include_directories(
    test_stl_validator
PRIVATE

    src/input/stl
)

target_compile_definitions(
    test_stl_validator
PRIVATE

    TEST_STL_DIR="${CMAKE_SOURCE_DIR}/tests/input_tests/stl_tests/stls"
)

set_input_stl_test_output(test_stl_validator)
