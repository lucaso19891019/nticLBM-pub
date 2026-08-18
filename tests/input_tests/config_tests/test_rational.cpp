#include "rational.hpp"
#include "test_framework.hpp"

#include <stdexcept>

using namespace ntic::lbm::common;

int main()
{
    //----------------------------------------------------------
    // Constructors
    //----------------------------------------------------------

    NTIC_EXPECT_EQ(Rational(), Rational(0));
    NTIC_EXPECT_EQ(Rational(5), Rational(5,1));
    NTIC_EXPECT_EQ(Rational(-5), Rational(-5,1));

    //----------------------------------------------------------
    // Normalization
    //----------------------------------------------------------

    NTIC_EXPECT_EQ(Rational(2,4), Rational(1,2));
    NTIC_EXPECT_EQ(Rational(4,2), Rational(2));
    NTIC_EXPECT_EQ(Rational(1,-3), Rational(-1,3));
    NTIC_EXPECT_EQ(Rational(-1,-3), Rational(1,3));
    NTIC_EXPECT_EQ(Rational(0,100), Rational(0));
    NTIC_EXPECT_EQ(Rational(100,-100), Rational(-1));

    //----------------------------------------------------------
    // Comparison
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(Rational(1,2)==Rational(2,4));
    NTIC_EXPECT_TRUE(Rational(1,2)!=Rational(3,4));

    NTIC_EXPECT_TRUE(Rational(1,2)< Rational(2,3));
    NTIC_EXPECT_TRUE(Rational(2,3)> Rational(1,2));

    NTIC_EXPECT_TRUE(Rational(2,3)<=Rational(2,3));
    NTIC_EXPECT_TRUE(Rational(2,3)>=Rational(2,3));

    //----------------------------------------------------------
    // Addition
    //----------------------------------------------------------

    NTIC_EXPECT_EQ(
        Rational(1,2)+Rational(1,3),
        Rational(5,6));

    NTIC_EXPECT_EQ(
        Rational(-1,2)+Rational(1,2),
        Rational(0));

    //----------------------------------------------------------
    // Subtraction
    //----------------------------------------------------------

    NTIC_EXPECT_EQ(
        Rational(5,6)-Rational(1,3),
        Rational(1,2));

    NTIC_EXPECT_EQ(
        Rational(1,2)-Rational(1,2),
        Rational(0));

    //----------------------------------------------------------
    // Multiplication
    //----------------------------------------------------------

    NTIC_EXPECT_EQ(
        Rational(2,3)*Rational(3,4),
        Rational(1,2));

    NTIC_EXPECT_EQ(
        Rational(-2,3)*Rational(3,4),
        Rational(-1,2));

    //----------------------------------------------------------
    // Division
    //----------------------------------------------------------

    NTIC_EXPECT_EQ(
        Rational(2,3)/Rational(4,5),
        Rational(5,6));

    NTIC_EXPECT_EQ(
        Rational(-2,3)/Rational(4,5),
        Rational(-5,6));

    //----------------------------------------------------------
    // Unary minus
    //----------------------------------------------------------

    NTIC_EXPECT_EQ(
        -Rational(1,2),
        Rational(-1,2));

    NTIC_EXPECT_EQ(
        -Rational(-1,2),
        Rational(1,2));

    //----------------------------------------------------------
    // Mixed arithmetic
    //----------------------------------------------------------

    NTIC_EXPECT_EQ(
        Rational(1,2)
        + Rational(1,3)
        - Rational(5,6),
        Rational(0));

    NTIC_EXPECT_EQ(
        Rational(2,3)
        * Rational(3,5)
        / Rational(2,5),
        Rational(1));

    //----------------------------------------------------------
    // Utility
    //----------------------------------------------------------

    NTIC_EXPECT_TRUE(Rational(5).isInteger());
    NTIC_EXPECT_TRUE(!Rational(5,2).isInteger());

    NTIC_EXPECT_TRUE(Rational().isZero());
    NTIC_EXPECT_TRUE(Rational(0,5).isZero());

    NTIC_EXPECT_TRUE(Rational(1,2).toDouble()==0.5);
    NTIC_EXPECT_TRUE(Rational(3,2).toDouble()==1.5);

    NTIC_EXPECT_TRUE(Rational(5).toString()=="5");
    NTIC_EXPECT_TRUE(Rational(-5).toString()=="-5");
    NTIC_EXPECT_TRUE(Rational(3,2).toString()=="3/2");
    NTIC_EXPECT_TRUE(Rational(-3,2).toString()=="-3/2");

    //----------------------------------------------------------
    // Exception
    //----------------------------------------------------------

    NTIC_EXPECT_THROW(
    	Rational(1) / Rational(0),
    	std::runtime_error);

    //----------------------------------------------------------

    ntic::lbm::test::printSummary();

    return 0;
}
