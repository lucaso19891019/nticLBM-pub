#pragma once
#include <string>
#include <unordered_map>

namespace ntic::lbm::config
{

class ModelAliasRegistry
{
public:
    ModelAliasRegistry();

    bool containsLatticeModel(const std::string& name) const;
    bool containsCollisionModel(const std::string& name) const;
    bool containsWallBoundary(const std::string& name) const;
    bool containsInletBoundary(const std::string& name) const;
    bool containsOutletBoundary(const std::string& name) const;

    const std::string& resolveLatticeModel(const std::string& name) const;
    const std::string& resolveCollisionModel(const std::string& name) const;
    const std::string& resolveWallBoundary(const std::string& name) const;
    const std::string& resolveInletBoundary(const std::string& name) const;
    const std::string& resolveOutletBoundary(const std::string& name) const;

private:
    using AliasMap = std::unordered_map<std::string, std::string>;

    static void registerAlias(
        AliasMap& aliases,
        const std::string& alias,
        const std::string& canonicalName);

    static bool contains(
        const AliasMap& aliases,
        const std::string& name);

    static const std::string& resolve(
        const AliasMap& aliases,
        const std::string& name,
        const std::string& kind);

    void registerBuiltins();

    AliasMap latticeModels_;
    AliasMap collisionModels_;
    AliasMap wallBoundaries_;
    AliasMap inletBoundaries_;
    AliasMap outletBoundaries_;
};

} // namespace ntic::lbm::config
