#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include "parameter_value.hpp"

namespace ntic::lbm::config
{

struct UserDefinedParameterData
{
    std::string identifier;

    ParameterValue value;

    bool bound;

    std::size_t lineNumber;
};

struct ParameterData
{
    std::unordered_map<
        std::string,
        ParameterValue> builtinValues;

    std::vector<
        UserDefinedParameterData> userDefinedValues;
};

} // namespace ntic::lbm::config
