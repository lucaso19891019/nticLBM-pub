#include "model_alias_registry.hpp"
#include "parameter_alias_registry.hpp"

#include <stdexcept>

namespace ntic::lbm::config
{

ModelAliasRegistry::ModelAliasRegistry()
{
    registerBuiltins();
}

void
ModelAliasRegistry::registerAlias(
    AliasMap& aliases,
    const std::string& alias,
    const std::string& canonicalName)
{
    const std::string key =
        ParameterAliasRegistry::normalize(alias);

    const auto result =
        aliases.emplace(key, canonicalName);

    if (!result.second
        && result.first->second != canonicalName)
    {
        throw std::runtime_error(
            "Conflicting model alias: '" + alias + "'.");
    }
}

bool
ModelAliasRegistry::contains(
    const AliasMap& aliases,
    const std::string& name)
{
    return aliases.find(
        ParameterAliasRegistry::normalize(name))
        != aliases.end();
}

const std::string&
ModelAliasRegistry::resolve(
    const AliasMap& aliases,
    const std::string& name,
    const std::string& kind)
{
    const auto iter =
        aliases.find(
            ParameterAliasRegistry::normalize(name));

    if (iter == aliases.end())
    {
        throw std::runtime_error(
            "Unsupported " + kind + ": '" + name + "'.");
    }

    return iter->second;
}

bool ModelAliasRegistry::containsLatticeModel(const std::string& n) const
{
    return contains(latticeModels_, n);
}

bool ModelAliasRegistry::containsCollisionModel(const std::string& n) const
{
    return contains(collisionModels_, n);
}

bool ModelAliasRegistry::containsWallBoundary(const std::string& n) const
{
    return contains(wallBoundaries_, n);
}

bool ModelAliasRegistry::containsInletBoundary(const std::string& n) const
{
    return contains(inletBoundaries_, n);
}

bool ModelAliasRegistry::containsOutletBoundary(const std::string& n) const
{
    return contains(outletBoundaries_, n);
}

const std::string&
ModelAliasRegistry::resolveLatticeModel(const std::string& n) const
{
    return resolve(latticeModels_, n, "lattice model");
}

const std::string&
ModelAliasRegistry::resolveCollisionModel(const std::string& n) const
{
    return resolve(collisionModels_, n, "collision model");
}

const std::string&
ModelAliasRegistry::resolveWallBoundary(const std::string& n) const
{
    return resolve(wallBoundaries_, n, "wall boundary model");
}

const std::string&
ModelAliasRegistry::resolveInletBoundary(const std::string& n) const
{
    return resolve(inletBoundaries_, n, "inlet boundary model");
}

const std::string&
ModelAliasRegistry::resolveOutletBoundary(const std::string& n) const
{
    return resolve(outletBoundaries_, n, "outlet boundary model");
}

void
ModelAliasRegistry::registerBuiltins()
{
    for (const char* value : {"D2Q9", "D3Q15", "D3Q19", "D3Q27"})
        registerAlias(latticeModels_, value, value);

    for (const char* value : {"BGK", "MRT", "TRT"})
        registerAlias(collisionModels_, value, value);

    registerAlias(collisionModels_, "Cascaded", "Cascaded");
    registerAlias(collisionModels_, "central moment", "Cascaded");
    registerAlias(collisionModels_, "central-moment", "Cascaded");
    registerAlias(collisionModels_, "central_moment", "Cascaded");
    registerAlias(collisionModels_, "CM", "Cascaded");
    registerAlias(collisionModels_, "Cumulant", "Cumulant");

    registerAlias(wallBoundaries_, "InterpolatedBB", "InterpolatedBB");
    registerAlias(wallBoundaries_, "interpolated bb", "InterpolatedBB");
    registerAlias(wallBoundaries_, "interpolated bounce back", "InterpolatedBB");
    registerAlias(wallBoundaries_, "IBB", "InterpolatedBB");
    registerAlias(wallBoundaries_, "IBM", "IBM");
    registerAlias(wallBoundaries_, "immersed", "IBM");
    registerAlias(wallBoundaries_, "immersed boundary", "IBM");

    const auto addFlowBoundary =
        [](AliasMap& map)
        {
            registerAlias(map, "Default", "Default");
            registerAlias(map, "standard", "Default");

            registerAlias(map, "ZouHeVelocity", "ZouHeVelocity");
            registerAlias(map, "Zou-He velocity", "ZouHeVelocity");
            registerAlias(map, "Zou He velocity", "ZouHeVelocity");
            registerAlias(map, "ZH vel", "ZouHeVelocity");
            registerAlias(map, "ZH velocity", "ZouHeVelocity");
            registerAlias(map, "zh-vel", "ZouHeVelocity");

            registerAlias(map, "ZouHePressure", "ZouHePressure");
            registerAlias(map, "Zou-He pressure", "ZouHePressure");
            registerAlias(map, "Zou He pressure", "ZouHePressure");
            registerAlias(map, "ZH pre", "ZouHePressure");
            registerAlias(map, "ZH pressure", "ZouHePressure");
            registerAlias(map, "zh-pre", "ZouHePressure");
        };

    addFlowBoundary(inletBoundaries_);
    addFlowBoundary(outletBoundaries_);
}

} // namespace ntic::lbm::config
