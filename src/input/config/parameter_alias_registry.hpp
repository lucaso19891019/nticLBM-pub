#pragma once

#include <cstddef>
#include <initializer_list>
#include <string>
#include <unordered_map>
#include <vector>

namespace ntic::lbm::config
{

class ParameterAliasRegistry
{
public:
    ParameterAliasRegistry();

    void registerAlias(
        const std::string& alias,
        const std::string& canonicalName);

    void registerAliases(
        const std::vector<std::string>& aliases,
        const std::string& canonicalName);

    void registerAliases(
        std::initializer_list<std::string> aliases,
        const std::string& canonicalName);

    bool contains(const std::string& name) const;

    const std::string& resolve(
        const std::string& name) const;

    std::size_t size() const;

    static std::string normalize(
        const std::string& name);

private:
    void registerNormalizedAlias(
        const std::string& normalizedAlias,
        const std::string& normalizedCanonicalName);

private:
    std::unordered_map<std::string, std::string> aliases_;
};

} // namespace ntic::lbm::config
