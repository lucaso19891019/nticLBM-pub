#include "rational.hpp"

#include <sstream>
#include <stdexcept>
#include <cstdlib>

namespace ntic::lbm::common
{

//==============================================================
// gcd
//==============================================================

int Rational::gcd(int a, int b)
{
    a = std::abs(a);
    b = std::abs(b);

    while (b != 0)
    {
        int t = a % b;
        a = b;
        b = t;
    }

    return (a == 0) ? 1 : a;
}

//==============================================================
// normalize
//==============================================================

void Rational::normalize()
{
    if (denominator_ == 0)
    {
        throw std::runtime_error(
            "Rational: denominator cannot be zero.");
    }

    // zero
    if (numerator_ == 0)
    {
        denominator_ = 1;
        return;
    }

    int g = gcd(numerator_, denominator_);

    numerator_ /= g;
    denominator_ /= g;

    // denominator always positive
    if (denominator_ < 0)
    {
        numerator_   = -numerator_;
        denominator_ = -denominator_;
    }
}

//==============================================================
// Constructors
//==============================================================

Rational::Rational()
    : numerator_(0),
      denominator_(1)
{
}

Rational::Rational(int numerator)
    : numerator_(numerator),
      denominator_(1)
{
}

Rational::Rational(
    int numerator,
    int denominator)
    : numerator_(numerator),
      denominator_(denominator)
{
    normalize();
}

//==============================================================
// Access
//==============================================================

int Rational::numerator() const
{
    return numerator_;
}

int Rational::denominator() const
{
    return denominator_;
}

//==============================================================
// Comparison
//==============================================================

bool Rational::operator==(const Rational& rhs) const
{
    return numerator_   == rhs.numerator_
        && denominator_ == rhs.denominator_;
}

bool Rational::operator!=(const Rational& rhs) const
{
    return !(*this == rhs);
}

bool Rational::operator<(const Rational& rhs) const
{
    return numerator_ * rhs.denominator_
        < rhs.numerator_ * denominator_;
}

bool Rational::operator<=(const Rational& rhs) const
{
    return !(*this > rhs);
}

bool Rational::operator>(const Rational& rhs) const
{
    return rhs < *this;
}

bool Rational::operator>=(const Rational& rhs) const
{
    return !(*this < rhs);
}

//==============================================================
// Arithmetic
//==============================================================

Rational Rational::operator+(const Rational& rhs) const
{
    return Rational(
        numerator_ * rhs.denominator_
            + rhs.numerator_ * denominator_,
        denominator_ * rhs.denominator_);
}

Rational Rational::operator-(const Rational& rhs) const
{
    return Rational(
        numerator_ * rhs.denominator_
            - rhs.numerator_ * denominator_,
        denominator_ * rhs.denominator_);
}

Rational Rational::operator*(const Rational& rhs) const
{
    return Rational(
        numerator_ * rhs.numerator_,
        denominator_ * rhs.denominator_);
}

Rational Rational::operator/(const Rational& rhs) const
{
    if (rhs.isZero())
    {
        throw std::runtime_error(
            "Rational: division by zero.");
    }

    return Rational(
        numerator_ * rhs.denominator_,
        denominator_ * rhs.numerator_);
}

Rational Rational::operator-() const
{
    return Rational(-numerator_, denominator_);
}

//==============================================================
// Utility
//==============================================================

bool Rational::isInteger() const
{
    return denominator_ == 1;
}

bool Rational::isZero() const
{
    return numerator_ == 0;
}

double Rational::toDouble() const
{
    return static_cast<double>(numerator_)
        / static_cast<double>(denominator_);
}

std::string Rational::toString() const
{
    std::ostringstream oss;

    if (denominator_ == 1)
        oss << numerator_;
    else
        oss << numerator_ << "/" << denominator_;

    return oss.str();
}

}

