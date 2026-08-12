#include "parameter_info.hpp"

namespace ntic::lbm::config
{

ParameterInfo::ParameterInfo()
    :
    category_(ParameterCategory::Physical),
    required_(false),
    derived_(false)
{
}

ParameterInfo::ParameterInfo(
    const std::string& name,
    const PhysicalDimension& dimension,
    ParameterCategory category,
    bool required,
    bool derived,
    const std::string& description)
    :
    name_(name),
    dimension_(dimension),
    category_(category),
    required_(required),
    derived_(derived),
    description_(description)
{
}

const std::string&
ParameterInfo::name() const
{
    return name_;
}

const PhysicalDimension&
ParameterInfo::dimension() const
{
    return dimension_;
}

ParameterCategory
ParameterInfo::category() const
{
    return category_;
}

bool
ParameterInfo::required() const
{
    return required_;
}

bool
ParameterInfo::derived() const
{
    return derived_;
}

const std::string&
ParameterInfo::description() const
{
    return description_;
}

bool
ParameterInfo::accepts(
    const PhysicalDimension& dimension) const
{
    return dimension == dimension_;
}

}
