#pragma once

#include <string>

namespace ntic::lbm::common
{

class Rational
{
public:

    //-----------------------------------------
    // Constructors
    //-----------------------------------------

    Rational();

    Rational(int numerator);

    Rational(int numerator,
             int denominator);

    //-----------------------------------------
    // Access
    //-----------------------------------------

    int numerator() const;

    int denominator() const;

    //-----------------------------------------
    // Comparison
    //-----------------------------------------

    bool operator==(const Rational&) const;
    bool operator!=(const Rational&) const;

    bool operator<(const Rational&) const;
    bool operator<=(const Rational&) const;
    bool operator>(const Rational&) const;
    bool operator>=(const Rational&) const;

    //-----------------------------------------
    // Arithmetic
    //-----------------------------------------

    Rational operator+(
        const Rational&) const;

    Rational operator-(
        const Rational&) const;

    Rational operator*(
        const Rational&) const;

    Rational operator/(
        const Rational&) const;

    Rational operator-() const;

    //-----------------------------------------
    // Utility
    //-----------------------------------------

    bool isInteger() const;

    bool isZero() const;

    double toDouble() const;

    std::string toString() const;

private:

    void normalize();

    static int gcd(int a,int b);

private:

    int numerator_;

    int denominator_;

};

}
