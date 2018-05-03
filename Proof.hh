//
// Created by algys on 4/19/18.
//

#ifndef DLEQ_PROOF_HH
#define DLEQ_PROOF_HH

#include <eosiolib/eosio.hpp>

#include <array>
#include <cstring>
#include <vector>
#include "Ecc.hh"


struct Proof {
    Proof() = default;

    Point g, m;
    Point h, z;
    Point a, b;

    std::vector<std::pair<uint256_t, uint256_t>> cr;

    bool verify() {
        if (!h.isOnCurve() || !g.isOnCurve())
            return false;

        if (!z.isOnCurve() || !m.isOnCurve())
            return false;

        for (auto pair: cr) {
            auto ch = h.scalarMult(pair.first);
            auto rg = g.scalarMult(pair.second);
            auto aa = rg.add(ch);

            auto cz = z.scalarMult(pair.first);
            auto rm = m.scalarMult(pair.second);
            auto bb = rm.add(cz);

            if (a != aa || b != bb)
                return false;
        }

        return true;
    }

    static Proof generate(Point g, Point m, uint256_t x, std::vector<uint256_t> const & cList) {
        // s must be large then Curve.N
        uint256_t s = Secp256k1.n;
        s += 1234;

        std::vector<std::pair<uint256_t, uint256_t>> cr;
        uint256_t cc, rr;
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
};

#endif //DLEQ_PROOF_HH
