#include "parameter.hpp"

#include <stdexcept>
#include <string>

namespace ntic::lbm::config
{

//==============================================================
// Constructors
//==============================================================

Parameter::Parameter()
    :
    info_(),
    value_(),
    hasValue_(false)
{
}

Parameter::Parameter(
    const ParameterInfo& info)
    :
    info_(info),
    value_(),
    hasValue_(false)
{
}

//==============================================================
// Metadata
//==============================================================

const ParameterInfo&
Parameter::info() const
{
    return info_;
}

const std::string&
Parameter::name() const
{
    return info_.name();
}

const PhysicalDimension&
Parameter::dimension() const
{
    return info_.dimension();
}

bool
Parameter::required() const
{
    return info_.required();
}

//==============================================================
// Value
//==============================================================

bool
Parameter::hasValue() const
{
    return hasValue_;
}

const Quantity&
Parameter::value() const
{
    if (!hasValue_)
    {
        throw std::runtime_error(
            "Parameter '" +
            info_.name() +
            "' has not been assigned a value.");
    }

    return value_;
}

void
Parameter::setValue(
    const Quantity& quantity)
{
    if (!info_.accepts(
            quantity.unit().dimension()))
    {
        throw std::invalid_argument(
            "Parameter '" +
            info_.name() +
            "' expects dimension " +
            info_.dimension().toString() +
            ", but received " +
            quantity.unit().dimension().toString() +
            ".");
    }

    value_ = quantity;

    hasValue_ = true;
}

void
Parameter::clear()
{
    value_ = Quantity();

    hasValue_ = false;
}

//==============================================================
// Utility
//==============================================================

std::string
Parameter::toString() const
{
    if (!hasValue_)
    {
        return info_.name() + " = <unset>";
    }

    return info_.name()
        + " = "
        + value_.toString();
}

} // namespace ntic::lbm::config
