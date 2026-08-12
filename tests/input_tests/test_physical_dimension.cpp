#include "physical_dimension.hpp"
#include "rational.hpp"
#include "test_framework.hpp"

using namespace ntic::lbm::config;
using namespace ntic::lbm::common;

int main()
{
    //----------------------------------------------------------
    // Constructors
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(
        PhysicalDimension().isDimensionless());

    NTIC_EXPECT_TRUE(
        !PhysicalDimension::Length().isDimensionless());

    //----------------------------------------------------------
    // Base Dimensions
    //----------------------------------------------------------

    auto L = PhysicalDimension::Length();
    auto M = PhysicalDimension::Mass();
    auto T = PhysicalDimension::Time();
    auto Theta = PhysicalDimension::Temperature();

    NTIC_EXPECT_EQ(
        L,
        PhysicalDimension(1,0,0,0));

    NTIC_EXPECT_EQ(
        M,
        PhysicalDimension(0,1,0,0));

    NTIC_EXPECT_EQ(
        T,
        PhysicalDimension(0,0,1,0));

    NTIC_EXPECT_EQ(
        Theta,
        PhysicalDimension(0,0,0,1));

    //----------------------------------------------------------
    // Derived Dimensions
    //----------------------------------------------------------

    NTIC_EXPECT_EQ(
        PhysicalDimension::Area(),
        L.pow(Rational(2)));

    NTIC_EXPECT_EQ(
        PhysicalDimension::Volume(),
        L.pow(Rational(3)));

    NTIC_EXPECT_EQ(
        PhysicalDimension::Velocity(),
        L / T);

    NTIC_EXPECT_EQ(
        PhysicalDimension::Acceleration(),
        L / T.pow(Rational(2)));

    NTIC_EXPECT_EQ(
        PhysicalDimension::Force(),
        M * PhysicalDimension::Acceleration());

    NTIC_EXPECT_EQ(
        PhysicalDimension::Pressure(),
        PhysicalDimension::Force()
            / PhysicalDimension::Area());

    NTIC_EXPECT_EQ(
        PhysicalDimension::Density(),
        M / PhysicalDimension::Volume());

    NTIC_EXPECT_EQ(
        PhysicalDimension::DynamicViscosity(),
        PhysicalDimension::Pressure() * T);

    NTIC_EXPECT_EQ(
        PhysicalDimension::KinematicViscosity(),
        PhysicalDimension::Area() / T);

    NTIC_EXPECT_EQ(
        PhysicalDimension::Energy(),
        PhysicalDimension::Force() * L);

    NTIC_EXPECT_EQ(
        PhysicalDimension::Power(),
        PhysicalDimension::Energy() / T);

    //----------------------------------------------------------
    // Algebra
    //----------------------------------------------------------

    NTIC_EXPECT_EQ(
        (L / T) * T,
        L);

    NTIC_EXPECT_EQ(
        (L * T) / T,
        L);

    NTIC_EXPECT_EQ(
        L.pow(Rational(2)),
        L * L);

    NTIC_EXPECT_EQ(
        L.pow(Rational(3)),
        L * L * L);

    //----------------------------------------------------------
    // Rational Exponent
    //----------------------------------------------------------

    auto sqrtL =
        L.pow(Rational(1,2));

    NTIC_EXPECT_TRUE(
        !sqrtL.isIntegral());

    NTIC_EXPECT_EQ(
        sqrtL.pow(Rational(2)),
        L);

    auto invL =
        L.pow(Rational(-1));

    NTIC_EXPECT_EQ(
        invL * L,
        PhysicalDimension::Dimensionless());

    NTIC_EXPECT_EQ(
        PhysicalDimension::Volume().pow(Rational(2,3)),
        PhysicalDimension::Area());

    //----------------------------------------------------------
    // Utility
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(
        PhysicalDimension::Dimensionless().isIntegral());

    NTIC_EXPECT_TRUE(
        PhysicalDimension::Velocity().isIntegral());

    NTIC_EXPECT_TRUE(
        !L.pow(Rational(1,2)).isIntegral());

    NTIC_EXPECT_TRUE(
        PhysicalDimension::Dimensionless().toString()
            == "[1]");

    NTIC_EXPECT_TRUE(
        PhysicalDimension::Length().toString()
            == "[L]");

    NTIC_EXPECT_TRUE(
        PhysicalDimension::Velocity().toString()
            == "[L T^-1]");

    NTIC_EXPECT_TRUE(
        L.pow(Rational(3,2)).toString()
            == "[L^(3/2)]");

    //----------------------------------------------------------

    ntic::lbm::test::printSummary();

    return 0;
}
