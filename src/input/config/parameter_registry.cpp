#include "parameter_registry.hpp"

#include <stdexcept>

namespace ntic::lbm::config
{

ParameterRegistry::ParameterRegistry()
{
}

void
ParameterRegistry::registerParameter(
    const ParameterInfo& info)
{
    auto result =
        parameters_.emplace(
            info.name(),
            info);

    if (!result.second)
    {
        throw std::runtime_error(
            "Parameter '" + info.name() +
            "' has already been registered.");
    }
}

bool
ParameterRegistry::contains(
    const std::string& name) const
{
    return parameters_.find(name)
        != parameters_.end();
}

const ParameterInfo&
ParameterRegistry::find(
    const std::string& name) const
{
    auto iter = parameters_.find(name);

    if (iter == parameters_.end())
    {
        throw std::runtime_error(
            "Unknown parameter: " + name);
    }

    return iter->second;
}

std::size_t
ParameterRegistry::size() const
{
    return parameters_.size();
}

}
