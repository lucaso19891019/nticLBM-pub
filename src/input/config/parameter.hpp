#pragma once

#include <string>

#include "parameter_info.hpp"
#include "quantity.hpp"

namespace ntic::lbm::config
{

class Parameter
{
public:

    //----------------------------------------------------------
    // Constructors
    //----------------------------------------------------------

    Parameter();

    explicit
    Parameter(
        const ParameterInfo& info);

public:

    //----------------------------------------------------------
    // Metadata
    //----------------------------------------------------------

    const ParameterInfo&
    info() const;

    const std::string&
    name() const;

    const PhysicalDimension&
    dimension() const;

    bool
    required() const;

public:

    //----------------------------------------------------------
    // Value
    //----------------------------------------------------------

    bool
    hasValue() const;

    const Quantity&
    value() const;

    void
    setValue(
        const Quantity& quantity);

    void
    clear();

public:

    //----------------------------------------------------------
    // Utility
    //----------------------------------------------------------

    std::string
    toString() const;

private:

    ParameterInfo info_;

    Quantity value_;

    bool hasValue_;
};

}
