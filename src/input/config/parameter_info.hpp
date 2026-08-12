#pragma once

#include <string>

#include "physical_dimension.hpp"

namespace ntic::lbm::config
{

enum class ParameterCategory
{
    General,
    Physical,
    Lattice,
    Output
};

class ParameterInfo
{
public:

    ParameterInfo();

    ParameterInfo(
        const std::string& name,
        const PhysicalDimension& dimension,
        ParameterCategory category,
        bool required,
        bool derived,
        const std::string& description = "");

public:

    const std::string& name() const;

    const PhysicalDimension& dimension() const;

    ParameterCategory category() const;

    bool required() const;

    bool derived() const;

    const std::string& description() const;

    bool accepts(const PhysicalDimension& dimension) const;

private:

    std::string name_;

    PhysicalDimension dimension_;

    ParameterCategory category_;

    bool required_;

    bool derived_;

    std::string description_;
};

}
