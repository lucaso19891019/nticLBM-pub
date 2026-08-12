#pragma once

#include <string>
#include <variant>

#include "quantity.hpp"

namespace ntic::lbm::config
{

using ParameterValue =
    std::variant<
        Quantity,
        bool,
        std::string>;

} // namespace ntic::lbm::config
