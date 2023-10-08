#ifndef ASSIGNMENT1_FRACTION_H
#define ASSIGNMENT1_FRACTION_H

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

class Fraction {
private:
    int p, q;

    void simplify() {
        bool sign = (p < 0) ^ (q < 0);
        p = std::abs(p);
        q = std::abs(q);
        int g = std::gcd(p, q);
        p /= g;
        q /= g;
        p *= (sign ? -1 : 1);
    }

    friend std::ostream& operator<<(std::ostream& stream, const Fraction& fraction);
    friend std::istream& operator>>(std::istream& stream, Fraction& fraction);

public:
    Fraction(int p, int q) : p(p), q(q) { simplify(); }
    Fraction(int p) : Fraction(p, 1) { }
    Fraction() : Fraction(0) { }

    Fraction& operator*=(const Fraction& other) {
        p *= other.p;
        q *= other.q;
        simplify();
        return *this;
    }

    Fraction operator*(const Fraction& other) const {
        Fraction f(*this);
        return f *= other;
    }

    Fraction& operator/=(const Fraction& other) {
        p *= other.q;
        q *= other.p;
        simplify();
        return *this;
    }

    Fraction operator/(const Fraction& other) const {
        Fraction f(*this);
        return f /= other;
    }

    Fraction& operator+=(const Fraction& other) {
        int l = std::lcm(q, other.q);
        int mp = l / q;
        p *= mp;
        q *= mp;
        p += (l / other.q) * other.p;
        simplify();
        return *this;
    }

    Fraction operator+(const Fraction& other) const {
        Fraction f(*this);
        return f += other;
    }

    Fraction& operator-=(const Fraction& other) {
        int l = std::lcm(q, other.q);
        int mp = l / q;
        p *= mp;
        q *= mp;
        p -= (l / other.q) * other.p;
        simplify();
        return *this;
    }

    Fraction operator-(const Fraction& other) const {
        Fraction f(*this);
        return f -= other;
    }

    Fraction& operator=(const int& other) {
        return *this = Fraction(other);
    }

    bool operator<(const Fraction& other) const {
        int l = std::lcm(q, other.q);
        int lp = p * (l / q), rp = other.p * (l / other.q);
        return lp < rp;
    }

    bool operator>(const Fraction& other) const {
        return other < *this;
    }

    bool operator==(const Fraction& other) const {
        return !(*this < other || *this > other);
    }

    bool operator<=(const Fraction& other) const {
        return !(*this > other);
    }

    bool operator>=(const Fraction& other) const {
        return !(*this < other);
    }

    bool operator<(const int& other) const {
        return *this < Fraction(other);
    }

    bool operator<=(const int& other) const {
        return *this <= Fraction(other);
    }

    bool operator==(const int& other) const {
        return *this == Fraction(other);
    }

    bool operator>(const int& other) const {
        return *this > Fraction(other);
    }

    bool operator>=(const int& other) const {
        return *this >= Fraction(other);
    }
};

std::ostream& operator<<(std::ostream& stream, const Fraction& fraction) {
    if (fraction.q == 1) {
        return stream << fraction.p;
    }
    return stream << fraction.p << '/' << fraction.q;
}

/**
 * Read fraction as integer number
 * @param stream input
 * @param fraction fraction to input
 * @return updated input
 */
std::istream& operator>>(std::istream& stream, Fraction& fraction) {
    stream >> fraction.p;
    fraction.q = 1;
    fraction.simplify();
    return stream;
}

#endif //ASSIGNMENT1_FRACTION_H
