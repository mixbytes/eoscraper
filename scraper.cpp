#include "scraper.hpp"

#include <set>
#include <map>
#include <vector>

#include "Proof.hh"
#include "Ecc.hh"

//count of participants for random generation (min 4)
constexpr static uint8_t participant_cnt = 10;

Point pointFromDatastream(eosio::datastream<unsigned char*> & ds) {
    eosio_assert(ds.remaining() >= 64, "Invalid Point");
    Curve curve;
    uint256_t x;
    uint256_t y;

    import_bits(curve.p, ds.pos(), ds.pos() + 8);
    ds.skip(8);
    import_bits(curve.a, ds.pos(), ds.pos() + 8);
    ds.skip(8);
    import_bits(curve.b, ds.pos(), ds.pos() + 8);
    ds.skip(8);
    import_bits(curve.n, ds.pos(), ds.pos() + 8);
    ds.skip(8);
    import_bits(curve.G.x, ds.pos(), ds.pos() + 8);
    ds.skip(8);
    import_bits(curve.G.x, ds.pos(), ds.pos() + 8);
    ds.skip(8);

    import_bits(x, ds.pos(), ds.pos() + 8);
    ds.skip(8);
    import_bits(y, ds.pos(), ds.pos() + 8);
    ds.skip(8);

    return {curve, x, y};
}

Proof proofFromDatastream(eosio::datastream<unsigned char*> & ds) {
    eosio_assert(ds.remaining() >= 384 + 16 * participant_cnt, "Incorrect Proof");

    Proof proof;
    proof.g = pointFromDatastream(ds);
    proof.m = pointFromDatastream(ds);
    proof.h = pointFromDatastream(ds);
    proof.z = pointFromDatastream(ds);
    proof.a = pointFromDatastream(ds);
    proof.b = pointFromDatastream(ds);

    for (auto k = 0; k < 0; ++k) {
        uint256_t c, r;
        import_bits(c, ds.pos(), ds.pos() + 8);
        ds.skip(8);
        import_bits(r, ds.pos(), ds.pos() + 8);
        ds.skip(8);

        proof.cr.push_back(std::make_pair(c, r));
    }

    return std::move(proof);
}


// @abi table
struct encshare {
    uint64_t id;
    account_name from;
    account_name to;

    std::vector<uint8_t> data;
    checksum256 commitment;
    std::vector<uint8_t> proof_data;
    //Proof proof;
    uint64_t dec_id;

    auto primary_key() const { return id; }

    EOSLIB_SERIALIZE(encshare, (id)(from)(to)(data)(commitment)(proof_data)(dec_id));
};

// @abi table
struct decshare {
    uint64_t id;
    uint64_t x;
    uint64_t y;

    auto primary_key() const { return id; }

    EOSLIB_SERIALIZE(decshare, (id)(x)(y));
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

//@ abi table
struct account {
    account_name owner;
    public_key pubkey;

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
    void bindkey(account_name sender, public_key pubkey) {
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
                      std::vector<uint8_t> data, std::vector<uint8_t> proof_data, checksum256 commitment) {

        auto rand_it = _randoms.find(rand_id);
        eosio_assert(rand_it != _randoms.end(), "Rand not found");

        auto senderPartId = participantId(*rand_it, sender);
        eosio_assert(senderPartId != -1, "Sender is not participant");

        auto receiverPartId = participantId(*rand_it, receiver);
        eosio_assert(receiverPartId != -1, "Receiver is not participant");

        eosio::datastream<unsigned char*> ds(proof_data.data(), proof_data.size());
        Proof proof = proofFromDatastream(ds);
        eosio_assert(proof.verify(), "Proof validate failed");

        eosio_assert(rand_it->enc_ids[(uint32_t)senderPartId * participant_cnt + receiverPartId] == -1, "share already pushed");

        uint64_t new_share_id;
        _encshares.emplace(sender, [&](auto & share) {
            new_share_id = _encshares.available_primary_key();
            share.id = new_share_id;
            share.from = sender;
            share.to = receiver;
            share.data = data;
            share.commitment = commitment;
            share.proof_data = proof_data;
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
    void pushdecshare(uint64_t rand_id, account_name sender, account_name from, uint64_t x, uint64_t y) {
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

        //TODO hash check

        uint64_t new_dec_id;
        _decshares.emplace(sender, [&](auto & obj) {
            new_dec_id = _decshares.available_primary_key();
            obj.id = new_dec_id;
            obj.x = x;
            obj.y = y;
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

    uint8_t participantId(random const & obj, account_name acc) {
        for (uint8_t i = 0; i < participant_cnt; ++i)
            if (obj.participants[i] == acc)
                return i;
        return -1;
    }
};

EOSIO_ABI( scraper, (bindkey)(initrand)(joinrand)(pushencshare)(pushdecshare))
