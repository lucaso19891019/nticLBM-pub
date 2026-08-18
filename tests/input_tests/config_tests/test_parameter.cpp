#include "parameter.hpp"
#include "parameter_info.hpp"
#include "quantity.hpp"
#include "unit.hpp"
#include "unit_expression.hpp"
#include "physical_dimension.hpp"

#include "test_framework.hpp"

#include <stdexcept>
#include <string>

using ntic::lbm::config::Parameter;
using ntic::lbm::config::ParameterCategory;
using ntic::lbm::config::ParameterInfo;
using ntic::lbm::config::PhysicalDimension;
using ntic::lbm::config::Quantity;
using ntic::lbm::config::Unit;
using ntic::lbm::config::UnitExpression;

int main()
{

//----------------------------------------------------------
// Default Constructor
//----------------------------------------------------------

{
    Parameter parameter;

    NTIC_EXPECT_TRUE(
        !parameter.hasValue());

    NTIC_EXPECT_THROW(
        parameter.value(),
        std::runtime_error);
}

//----------------------------------------------------------
// Constructor
//----------------------------------------------------------

ParameterInfo densityInfo(
    "density",
    PhysicalDimension::Density(),
    ParameterCategory::Physical,
    true,
    false,
    "Fluid density");

Parameter parameter(
    densityInfo);

NTIC_EXPECT_TRUE(
    !parameter.hasValue());

NTIC_EXPECT_TRUE(
    parameter.name()
    ==
    "density");

NTIC_EXPECT_EQ(
    parameter.dimension(),
    PhysicalDimension::Density());

NTIC_EXPECT_TRUE(
    parameter.required());

NTIC_EXPECT_TRUE(
    parameter.info().name()
    ==
    densityInfo.name());

NTIC_EXPECT_EQ(
    parameter.info().dimension(),
    densityInfo.dimension());

NTIC_EXPECT_TRUE(
    parameter.info().category()
    ==
    densityInfo.category());

NTIC_EXPECT_TRUE(
    parameter.info().required()
    ==
    densityInfo.required());

NTIC_EXPECT_TRUE(
    parameter.info().derived()
    ==
    densityInfo.derived());

NTIC_EXPECT_TRUE(
    parameter.info().description()
    ==
    densityInfo.description());

//----------------------------------------------------------
// setValue()
//----------------------------------------------------------

{
    Quantity rho(

        998.0,

        UnitExpression(
            PhysicalDimension::Density(),
            1.0));

    parameter.setValue(
        rho);

    NTIC_EXPECT_TRUE(
        parameter.hasValue());

    NTIC_EXPECT_TRUE(
        parameter.value()
        ==
        rho);
}

//----------------------------------------------------------
// Wrong Dimension
//----------------------------------------------------------

{
    Quantity length(

        1.0,

        UnitExpression(
            Unit::Meter()));

    NTIC_EXPECT_THROW(

        parameter.setValue(
            length),

        std::invalid_argument);
}

//----------------------------------------------------------
// value()
//----------------------------------------------------------

{
    Quantity rho(

        1000.0,

        UnitExpression(
            PhysicalDimension::Density(),
            1.0));

    parameter.setValue(
        rho);

    NTIC_EXPECT_TRUE(

        parameter.value()
        ==
        rho);
}

//----------------------------------------------------------
// clear()
//----------------------------------------------------------

{
    Quantity rho(

        998.0,

        UnitExpression(
            PhysicalDimension::Density(),
            1.0));

    parameter.setValue(
        rho);

    NTIC_EXPECT_TRUE(
        parameter.hasValue());

    parameter.clear();

    NTIC_EXPECT_TRUE(
        !parameter.hasValue());

    NTIC_EXPECT_THROW(

        parameter.value(),

        std::runtime_error);
}

//----------------------------------------------------------
// toString()
//----------------------------------------------------------

{
    std::string text =
        parameter.toString();

    NTIC_EXPECT_TRUE(

        text.find("density")

        !=

        std::string::npos);
}

{
    Quantity rho(

        998.0,

        UnitExpression(
            PhysicalDimension::Density(),
            1.0));

    parameter.setValue(
        rho);

    std::string text =
        parameter.toString();

    NTIC_EXPECT_TRUE(

        text.find("density")

        !=

        std::string::npos);

    NTIC_EXPECT_TRUE(

        text.find("998")

        !=

        std::string::npos);
}

//----------------------------------------------------------
// Replace Existing Value
//----------------------------------------------------------

{
    Quantity rho1(

        900.0,

        UnitExpression(
            PhysicalDimension::Density(),
            1.0));

    Quantity rho2(

        998.0,

        UnitExpression(
            PhysicalDimension::Density(),
            1.0));

    parameter.setValue(
        rho1);

    NTIC_EXPECT_TRUE(

        parameter.value()
        ==
        rho1);

    parameter.setValue(
        rho2);

    NTIC_EXPECT_TRUE(

        parameter.value()
        ==
        rho2);
}

//----------------------------------------------------------
// Multiple Parameters
//----------------------------------------------------------

{
    ParameterInfo viscosityInfo(

        "viscosity",

        PhysicalDimension::KinematicViscosity(),

        ParameterCategory::Physical,

        true,

        false,

        "");

    Parameter viscosity(
        viscosityInfo);

    Quantity nu(

        1.0e-6,

        UnitExpression(
            PhysicalDimension::KinematicViscosity(),
            1.0));

    viscosity.setValue(
        nu);

    NTIC_EXPECT_TRUE(

        viscosity.hasValue());

    NTIC_EXPECT_TRUE(

        viscosity.value()
        ==
        nu);

    NTIC_EXPECT_TRUE(

        viscosity.name()

        ==

        "viscosity");
}

//----------------------------------------------------------

ntic::lbm::test::printSummary();

return 0;

}
