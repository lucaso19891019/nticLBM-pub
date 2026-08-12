#include "parameter_alias_registry.hpp"

#include <cctype>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace ntic::lbm::config
{

ParameterAliasRegistry::ParameterAliasRegistry()
    : aliases_()
{
}

void
ParameterAliasRegistry::registerAlias(
    const std::string& alias,
    const std::string& canonicalName)
{
    registerNormalizedAlias(
        normalize(alias),
        normalize(canonicalName));
}

void
ParameterAliasRegistry::registerAliases(
    const std::vector<std::string>& aliases,
    const std::string& canonicalName)
{
    const std::string canonical =
        normalize(canonicalName);

    std::unordered_set<std::string> normalizedAliases;
    normalizedAliases.reserve(aliases.size() + 1);

    // The canonical name always resolves to itself.
    normalizedAliases.insert(canonical);

    for (const std::string& alias : aliases)
    {
        normalizedAliases.insert(normalize(alias));
    }

    // Validate the whole batch first so conflicts do not cause
    // partial registration.
    for (const std::string& alias : normalizedAliases)
    {
        const auto iter = aliases_.find(alias);

        if (iter != aliases_.end()
            && iter->second != canonical)
        {
            throw std::runtime_error(
                "Parameter alias '" + alias
                + "' is already mapped to canonical name '"
                + iter->second
                + "', and cannot be remapped to '"
                + canonical + "'.");
        }
    }

    for (const std::string& alias : normalizedAliases)
    {
        aliases_.emplace(alias, canonical);
    }
}

void
ParameterAliasRegistry::registerAliases(
    std::initializer_list<std::string> aliases,
    const std::string& canonicalName)
{
    registerAliases(
        std::vector<std::string>(
            aliases.begin(),
            aliases.end()),
        canonicalName);
}

bool
ParameterAliasRegistry::contains(
    const std::string& name) const
{
    return aliases_.find(normalize(name))
        != aliases_.end();
}

const std::string&
ParameterAliasRegistry::resolve(
    const std::string& name) const
{
    const auto iter = aliases_.find(normalize(name));

    if (iter == aliases_.end())
    {
        throw std::runtime_error(
            "Unknown parameter alias: '" + name + "'.");
    }

    return iter->second;
}

std::size_t
ParameterAliasRegistry::size() const
{
    return aliases_.size();
}

std::string
ParameterAliasRegistry::normalize(
    const std::string& name)
{
    std::string normalized;
    normalized.reserve(name.size());

    bool separatorPending = false;

    for (char character : name)
    {
        const unsigned char value =
            static_cast<unsigned char>(character);

        if (std::isalpha(value) || std::isdigit(value))
        {
            if (separatorPending && !normalized.empty())
            {
                normalized.push_back('_');
            }

            normalized.push_back(
                std::isalpha(value)
                    ? static_cast<char>(std::tolower(value))
                    : character);

            separatorPending = false;
            continue;
        }

        if (character == '_'
            || character == '-'
            || std::isspace(value))
        {
            if (!normalized.empty())
            {
                separatorPending = true;
            }
            continue;
        }

        throw std::invalid_argument(
            "Invalid character '"
            + std::string(1, character)
            + "' in parameter name '" + name
            + "'. Parameter names may contain only letters, "
              "digits, whitespace, underscores, and hyphens.");
    }

    if (normalized.empty())
    {
        throw std::invalid_argument(
            "Parameter name cannot be empty.");
    }

    if (!std::isalpha(
            static_cast<unsigned char>(normalized.front())))
    {
        throw std::invalid_argument(
            "Parameter name '" + name
            + "' must begin with a letter.");
    }

    return normalized;
}

void
ParameterAliasRegistry::registerNormalizedAlias(
    const std::string& normalizedAlias,
    const std::string& normalizedCanonicalName)
{
    const auto iter = aliases_.find(normalizedAlias);

    if (iter == aliases_.end())
    {
        aliases_.emplace(
            normalizedAlias,
            normalizedCanonicalName);
        return;
    }

    if (iter->second != normalizedCanonicalName)
    {
        throw std::runtime_error(
            "Parameter alias '" + normalizedAlias
            + "' is already mapped to canonical name '"
            + iter->second
            + "', and cannot be remapped to '"
            + normalizedCanonicalName + "'.");
    }
}

} // namespace ntic::lbm::config
