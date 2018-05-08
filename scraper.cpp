#include "scraper.hpp"

#include <set>
#include <map>
#include <vector>

#include "Proof.hh"
#include "Ecc.hh"

//count of participants for random generation (min 4)
constexpr static uint8_t participant_cnt = 10;

static const UInt256 seed_for_g = fromString("21343453542133456453212");
static const UInt256 seed_for_m = fromString("56987654354236234435425");

static const Point base_point { Secp256k1, Secp256k1.G.x, Secp256k1.G.y };
static const Point gen_g = base_point.scalarMult(seed_for_g);
static const Point gen_m = base_point.scalarMult(seed_for_m);


struct binuint256_t {
    checksum256 data;
};

struct binpoint {
    binuint256_t x, y;
};

struct binproof {
    binpoint g,m,h,z,a,b;
    binuint256_t c, r;
};

UInt256 fromBin(binuint256_t const & data) {
    UInt256 x;
    memcpy(x.data, data.data.hash, sizeof(data.data));
    return std::move(x);
}

Point fromBin(binpoint const &data) {
    Point point {Secp256k1, fromBin(data.x), fromBin(data.y)};
    eosio_assert(point.isOnCurve(), "Point not in curve");

    return std::move(point);
}

Proof fromBin(binproof const &data) {
    return {
        fromBin(data.g), fromBin(data.m),
        fromBin(data.h), fromBin(data.z),
        fromBin(data.a), fromBin(data.b),
        fromBin(data.c), fromBin(data.r)
    };
}


// @abi table
struct encshare {
    uint64_t id;
    account_name from;
    account_name to;

    binpoint data;
    binpoint commitment;
    binproof proof;
    uint64_t dec_id;

    auto primary_key() const { return id; }

    EOSLIB_SERIALIZE(encshare, (id)(from)(to)(data)(commitment)(proof)(dec_id));
};

// @abi table
struct decshare {
    uint64_t id;
    binpoint s;
    binproof proof;

    auto primary_key() const { return id; }

    EOSLIB_SERIALIZE(decshare, (id)(s)(proof));
};

// @abi table
struct random {
    enum state : uint8_t {
        wait_joins = 0,
        wait_enc_shares,
        wait_dec_shares,
        done,
        error,
    };

    uint64_t id;
    uint8_t state;
    uint8_t joined_cnt;
    std::vector<account_name> participants;

    uint32_t pushed_cnt;
    std::vector<uint64_t> enc_ids;
    std::vector<uint8_t> dec_cnts;

    auto primary_key()const { return id; }

    EOSLIB_SERIALIZE(random, (id)(state)(participants)(enc_ids)(dec_cnts));
};

// @abi table
struct account {
    account_name owner;
    binpoint pubkey;

    auto primary_key()const { return owner; }

    EOSLIB_SERIALIZE(account, (owner)(pubkey));
};

class scraper : public eosio::contract {
public:
    scraper( account_name self ) :
        contract(self),
        _accounts(_self, _self),
        _randoms( _self, _self),
        _encshares(_self, _self),
        _decshares(_self, _self)
    { }

    // @abi action
    void bindkey(account_name sender, binpoint pubkey) {
        require_auth(sender);
        auto acc = _accounts.find(sender);
        if (acc != _accounts.end())
            _accounts.modify(acc, sender, [&] (auto & acc) {
                acc.pubkey = pubkey;
            });
        else
            _accounts.emplace(sender, [&](auto & acc) {
                acc.owner = sender;
                acc.pubkey = pubkey;
            });
    }

    // @abi action
    void initrand(account_name sender) {
        require_auth(sender);
        auto acc = _accounts.find(sender);
        eosio_assert(acc != _accounts.end(), "Need to bind public key");

        _randoms.emplace(sender, [&](auto & rand) {
            rand.id = _randoms.available_primary_key();
            rand.state = random::state::wait_joins;
            rand.participants.reserve(participant_cnt);
            rand.enc_ids.reserve(participant_cnt*participant_cnt);
            rand.dec_cnts.reserve(participant_cnt);
            rand.participants[0] = sender;
            rand.joined_cnt = 1;
            rand.pushed_cnt = 0;
            for(auto k = 0; k < participant_cnt * participant_cnt; ++k)
                rand.enc_ids[k] = -1;
        });
    }

    // @abi action
    void joinrand(uint64_t rand_id, account_name sender) {
        require_auth(sender);
        auto acc = _accounts.find(sender);
        eosio_assert(acc != _accounts.end(), "Need to bind public key");

        auto rand_it = _randoms.find(rand_id);
        eosio_assert(rand_it != _randoms.end(), "Rand not found");

        eosio_assert(participantId(*rand_it, sender) != -1, "Already joinded");
        eosio_assert(rand_it->state == random::state::wait_joins, "Invalid state");

        _randoms.modify(rand_it, sender, [&](auto & rand) {
            rand.participants[rand.joined_cnt++] = sender;
            if (rand.joined_cnt == participant_cnt)
                rand.state = random::state::wait_enc_shares;
        });
    }

    // @abi action
    void pushencshare(uint64_t rand_id, account_name sender, account_name receiver,
                      binpoint data, binpoint commitment, binproof proof)
    {
        require_auth(sender);
        auto senderAccIt = _accounts.find(sender);
        eosio_assert(senderAccIt != _accounts.end(), "Sender not found");

        auto receiverAccIt = _accounts.find(receiver);
        eosio_assert(receiverAccIt != _accounts.end(), "Receiver not found");

        auto rand_it = _randoms.find(rand_id);
        eosio_assert(rand_it != _randoms.end(), "Rand not found");

        auto senderPartId = participantId(*rand_it, sender);
        eosio_assert(senderPartId != -1, "Sender is not participant");

        auto receiverPartId = participantId(*rand_it, receiver);
        eosio_assert(receiverPartId != -1, "Receiver is not participant");

        eosio_assert(rand_it->enc_ids[(uint32_t)senderPartId * participant_cnt + receiverPartId] == -1, "share already pushed");

        auto pk = fromBin(receiverAccIt->pubkey);
        auto pr = fromBin(proof);
        eosio_assert(gen_g == pr.g, "Invalid proof (gen_g != g)");
        eosio_assert(pk == pr.m, "Invalid proof (pk != m)");
        eosio_assert(pr.verify(), "Proof validate failed");

        uint64_t new_share_id;
        _encshares.emplace(sender, [&](auto & share) {
            new_share_id = _encshares.available_primary_key();
            share.id = new_share_id;
            share.from = sender;
            share.to = receiver;
            share.data = data;
            share.commitment = commitment;
            share.proof = proof;
            share.dec_id = (uint64_t)-1;
        });

        _randoms.modify(rand_it, sender, [&](auto & rand) {
            rand.enc_ids[(uint32_t)senderPartId * participant_cnt + receiverPartId] = new_share_id;
            rand.pushed_cnt++;

            if (rand.pushed_cnt == (uint32_t)participant_cnt * participant_cnt)
                rand.state = random::state::wait_dec_shares;
        });
    }

    // @abi action
    void pushdecshare(uint64_t rand_id, account_name sender, account_name from,
                      binpoint s, binproof proof)
    {
        require_auth(sender);
        auto senderAccIt = _accounts.find(sender);
        eosio_assert(senderAccIt != _accounts.end(), "Sender not found");

        auto fromAccIt = _accounts.find(from);
        eosio_assert(fromAccIt != _accounts.end(), "From not found");

        auto rand_it = _randoms.find(rand_id);
        eosio_assert(rand_it != _randoms.end(), "Rand not found");

        eosio_assert(rand_it->state == random::state::wait_dec_shares, "Invalid state");

        auto senderPartId = participantId(*rand_it, from);
        eosio_assert(senderPartId != -1, "Sender is not participant");

        auto receiverPartId = participantId(*rand_it, sender);
        eosio_assert(receiverPartId != -1, "Receiver is not participant");

        auto enc_id = rand_it->enc_ids[(uint32_t)senderPartId * participant_cnt + receiverPartId];
        auto encshare_it = _encshares.find(enc_id);
        eosio_assert(encshare_it != _encshares.end(), "Share not found");
        eosio_assert(encshare_it->dec_id == -1, "Already pushed");

        auto encdata = fromBin(encshare_it->data);
        auto pk = fromBin(senderAccIt->pubkey);
        auto ss = fromBin(s);
        auto pr = fromBin(proof);
        eosio_assert(gen_g == pr.g, "Invalid proof (gen_g != g)");
        eosio_assert(encdata == pr.z, "Invalid proof (encdata != g)");
        eosio_assert(ss == pr.m, "Invalid proof (m != s)");
        eosio_assert(pk == pr.h, "Invalid proof (pk != h)");
        eosio_assert(pr.verify(), "Proof validate failed");

        uint64_t new_dec_id;
        _decshares.emplace(sender, [&](auto & obj) {
            new_dec_id = _decshares.available_primary_key();
            obj.id = new_dec_id;
            obj.s = s;
            obj.proof = proof;
        });

        _encshares.modify(encshare_it, sender, [&](auto & obj) {
            obj.dec_id = new_dec_id;
        });

        _randoms.modify(rand_it, sender, [&](auto & rand) {
            rand.dec_cnts[from]++;
        });
    }

private:
    eosio::multi_index<N(accounts), account> _accounts;
    eosio::multi_index<N(randoms), random> _randoms;
    eosio::multi_index<N(encshare), encshare> _encshares;
    eosio::multi_index<N(decshare), decshare> _decshares;

private:
    uint8_t participantId(random const & obj, account_name acc) {
        for (uint8_t i = 0; i < participant_cnt; ++i)
            if (obj.participants[i] == acc)
                return i;
        return -1;
    }
};

EOSIO_ABI( scraper, (bindkey)(initrand)(joinrand)(pushencshare)(pushdecshare))
