//
// Created by algys on 4/26/18.
//

#ifndef DLEQ_ECC_HH
#define DLEQ_ECC_HH

#include <utility>
#include "BigNum.h"

typedef BigNum::Num<256> UInt256;
typedef BigNum::Num<272> UInt272;
typedef BigNum::Num<576> UInt576;

class Curve {
public:
    UInt256 p;
    UInt256 a;
    UInt256 b;
    UInt256 n;

    struct {
        UInt256 x;
        UInt256 y;
    } G;
};

const UInt256 fromString(std::string const & str)
{
    UInt256 ret;
    auto i = 0;
    while (i < str.size() && isdigit(str[i]))
    {
        auto t = ret * 10u;
        ret = t;
        ret = ret + (unsigned)(str[i] - '0');
        i++;
    }
    return ret;
}


// Secp256k1 curve
const Curve Secp256k1 = {
    fromString("115792089237316195423570985008687907853269984665640564039457584007908834671663"), //115792089237316195423570985008687907853269984665640564039457584007908834671663
    UInt256(0),
    UInt256(7),
    fromString("115792089237316195423570985008687907852837564279074904382605163141518161494337"), //115792089237316195423570985008687907852837564279074904382605163141518161494337
    {
        fromString("55066263022277343669578718895168534326250603453777594175500187360389116729240"), //55066263022277343669578718895168534326250603453777594175500187360389116729240
        fromString("32670510020758816978083085130507043184471273380659243275938904335757337482424") //32670510020758816978083085130507043184471273380659243275938904335757337482424
    }
};

class Point {
public:
    Point() :
        _curve{Secp256k1}
    { }

    Point(UInt256 const & x, UInt256 const & y) :
        _x{x}, _y{y}, _curve{Secp256k1}
    { }

    Point(Curve curve, UInt256 const & x, UInt256 const & y) :
        _x{x}, _y{y}, _curve{std::move(curve)}
    { }

    //R = P + Q
    Point add(Point const & point) const{
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

        auto t1 = (UInt272)y() + curve().p - point.y();
        auto t2 = (UInt272)x() + curve().p - point.x();
        auto m = div(t1, t2);

        return lineIntersect(point, m);
    }

    //R = P+P
    Point dbl() const{
        if (identity())
            return *this;

        auto t = div(x() * x() * 3 + curve().a , y() * 2);
        return lineIntersect(*this, t);
    }

    Point scalarMult(UInt256 scalar) const {
        Point q = *this;
        Point r = identityPoint();
        while(scalar > 0) {
            if ((scalar[0] & 1) != 0)
                r = r.add(q);
            scalar = scalar >> 1;
            if (scalar != 0)
                q = q.dbl();
        }
        return r;
    }

    bool isOnCurve() const {
        auto y2 = modExp(y(), UInt256(2), curve().p);
        auto x3ab = mul((UInt272)modExp(x(), UInt256(2), curve().p) + curve().a, x()) + curve().b;
        return y2 == x3ab;
    }

    friend bool operator== (Point const & lhs, Point const & rhs) {
        return lhs._x == rhs._x && lhs._y == rhs._y;
    }

    friend bool operator!= (Point const & lhs, Point const & rhs) {
        return !(lhs == rhs);
    }

    const Curve & curve() const { return _curve; }
    const UInt256 & x() const { return _x; }
    const UInt256 & y() const { return _y; }

    bool identity() const { return curve().p == x(); }

    Point identityPoint() const{
        return {curve().p, 0};
    }

private:
    UInt256 mul(UInt256 const & a, UInt256 const & b) const{
        return ((a * b) % curve().p);
    }

    UInt256 div(UInt576 const & num, UInt576 const & den) const{
        auto invDen = modInv((BigNum::Num<257>)(den % curve().p), (BigNum::Num<257>)curve().p);
        return mul(num % curve().p, invDen);
    }

    Point lineIntersect(Point const & point, UInt256 const & m) const{
        auto v = ((UInt272)y() + curve().p - (m * x()) % curve().p) % curve().p;
        auto x = (m * m + curve().p + curve().p - this->x() - point.x()) % curve().p;
        auto y = ((UInt272)curve().p + curve().p - (m * x) % curve().p - v) % curve().p;
        return {curve(), x, y};
    }

private:
    Curve _curve{};
    UInt256 _x{};
    UInt256 _y{};
};


#endif //DLEQ_ECC_HH
