//
// Created by algys on 4/26/18.
//

#ifndef DLEQ_ECC_HH
#define DLEQ_ECC_HH

#include <assert.h>
#include <utility>
#include "boost/multiprecision/cpp_int.hpp"

using namespace boost::multiprecision;

class Curve {
public:
    uint256_t p;
    uint256_t a;
    uint256_t b;
    uint256_t n;

    struct {
        uint256_t x;
        uint256_t y;
    } G;
};

// Secp256k1 curve
const Curve Secp256k1 = {
    uint256_t("115792089237316195423570985008687907853269984665640564039457584007908834671663"),
    uint256_t("0"),
    uint256_t("7"),
    uint256_t("115792089237316195423570985008687907852837564279074904382605163141518161494337"),
    {
        uint256_t("55066263022277343669578718895168534326250603453777594175500187360389116729240"),
        uint256_t("32670510020758816978083085130507043184471273380659243275938904335757337482424")
    }
};


int1024_t gcdExtended(int1024_t const & a, int1024_t const & b, int1024_t & resX) {
    auto lastR = abs(a), r = abs(b);
    int1024_t x = 0, lastX = 1, quot = 0, t{};
    while (r > 0) {
        quot = lastR / r;
        t = r;
        r = lastR % r;
        lastR = t;

        t = x;
        x = lastX - quot * x;
        lastX = t;
    }

    resX = lastX;
    return lastR;
}

int1024_t modInverse(int1024_t const &a, int1024_t const &m) {
    int1024_t g, x;
    g = gcdExtended(a, m, x);

    assert(g == 1);

    return x % m;
}


class Point {
public:
    Point() :
        _curve{Secp256k1}
    { }

    Point(uint256_t const & x, uint256_t const & y) :
        _x{x}, _y{y}, _curve{Secp256k1}
    { }

    Point(Curve curve, uint256_t const & x, uint256_t const & y) :
        _x{x}, _y{y}, _curve{std::move(curve)}
    { }

    //R = P + Q
    Point add(Point const & point) {
        if (identity())
            return point;
        if (point.identity())
            return *this;

        if (x() == point.x()) {
            if (y() == point.y())
                return dbl();
            else
                return identityPoint();
        }


        auto m = div((int1024_t)y() + curve().p - (int1024_t)point.y(), (int1024_t)x() + curve().p - (int1024_t)point.x());

        return lineIntersect(point, m);
    }

    //R = P+P
    Point dbl() {
        if (identity())
            return *this;

        auto t = div((int1024_t)x() * (int1024_t)x() * 3 + curve().a , (int1024_t)y() * 2);
        return lineIntersect(*this, t);
    }

    Point scalarMult(uint256_t scalar) {
        Point q = *this;
        Point r = identityPoint();
        while(scalar > 0) {
            if (scalar % 2 != 0)
                r = r.add(q);
            scalar = scalar / 2;
            if (scalar != 0)
                q = q.dbl();
        }
        return r;
    }

    bool isOnCurve() {
        auto y2 = ((int1024_t)y() * (int1024_t)y()) % (int1024_t)curve().p;
        auto x3ab = (mul(((int1024_t)x() * (int1024_t)x()) % curve().p + curve().a, (int1024_t)x()) + curve().b) % curve().p;
        return y2 == x3ab;
    }

    friend bool operator== (Point const & lhs, Point const & rhs) {
        return lhs._x == rhs._x && lhs._y == rhs._y;
    }

    friend bool operator!= (Point const & lhs, Point const & rhs) {
        return !(lhs == rhs);
    }

    const Curve & curve() const { return _curve; }
    const uint256_t & x() const { return _x; }
    const uint256_t & y() const { return _y; }

    bool identity() const { return curve().p == x(); }

    Point identityPoint() {
        return {curve().p, 0};
    }

private:
    int1024_t mul(int1024_t const & a, int1024_t const & b) {
        return ((a * b) % curve().p);
    }

    int1024_t div(int1024_t const & num, int1024_t const & den) {
        auto invDen = modInverse(den % curve().p, curve().p);
        return mul(num % curve().p, invDen);
    }

    Point lineIntersect(Point const & point, int1024_t const & m) {
        auto v = ((int1024_t)y() + curve().p - (m * (int1024_t)x()) % curve().p) % curve().p;
        auto x = (m * m + curve().p - (int1024_t)this->x() + curve().p - (int1024_t)point.x()) % curve().p;
        auto y = (curve().p - (m * x) % curve().p + curve().p - v) % curve().p;
        return {curve(), (uint256_t)x, (uint256_t)y};
    }

private:
    Curve _curve{};
    uint256_t _x{};
    uint256_t _y{};
};


#endif //DLEQ_ECC_HH
