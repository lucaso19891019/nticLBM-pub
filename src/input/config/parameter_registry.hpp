#pragma once

#include <string>
#include <unordered_map>

#include "parameter_info.hpp"

namespace ntic::lbm::config
{

class ParameterRegistry
{
public:

    ParameterRegistry();

    void registerParameter(
        const ParameterInfo& info);

    bool contains(
        const std::string& name) const;

    const ParameterInfo&
    find(
        const std::string& name) const;

    std::size_t size() const;

private:

    std::unordered_map<
        std::string,
        ParameterInfo> parameters_;
};

}
