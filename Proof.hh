//
// Created by algys on 4/19/18.
//

#ifndef DLEQ_PROOF_HH
#define DLEQ_PROOF_HH

#include <array>
#include <cstring>
#include <vector>
#include "Ecc.hh"


struct Proof {
    Proof() = default;

    Point g, m;
    Point h, z;
    Point a, b;

    UInt256 c, r;

    bool verify() {
        if (!h.isOnCurve() || !g.isOnCurve())
            return false;

        if (!z.isOnCurve() || !m.isOnCurve())
            return false;

        auto ch = h.scalarMult(c);
        auto rg = g.scalarMult(r);
        auto aa = rg.add(ch);

        auto cz = z.scalarMult(c);
        auto rm = m.scalarMult(r);
        auto bb = rm.add(cz);

        return (a == aa && b == bb);

    }

    /*
    static Proof generate(Point g, Point m, UInt256 x, std::vector<UInt256> const & cList) {
        // s must be large then Curve.N
        UInt256 s = Secp256k1.n;
        s = s + UInt256(rand());

        std::vector<std::pair<UInt256, UInt256>> cr;
        UInt256 cc, rr;
        for (auto c: cList) {
            cc = c % Secp256k1.n;
            rr = (s - cc * x) % Secp256k1.n;
            cr.push_back(std::make_pair(cc, rr));
        }

        auto a = g.scalarMult(s);
        auto b = m.scalarMult(s);

        return {
            g, m,
            g.scalarMult(x), m.scalarMult(x),
            a, b,
            std::move(cr)
        };
    }
     */
};

#endif //DLEQ_PROOF_HH
