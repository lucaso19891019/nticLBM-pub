#pragma once

#include <string>

#include "physical_dimension.hpp"

namespace ntic::lbm::config
{

class SIUnitFormatter
{
public:

    static std::string
    format(
        const PhysicalDimension& dimension);
};

} // namespace ntic::lbm::config
