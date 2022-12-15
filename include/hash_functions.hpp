#include <vector>
#include <cstdint>
#include "kmer_factory.hpp"

#pragma once

class RollingHasher1
{
    public:

        // Hash table size (mod)
        uint64_t q;
        // Prime for rolling hash
        //uint64_t P;
        // Previous hash value for rolling hash
        uint64_t current_hash;
        // k-mer length
        uint64_t m;
        // Bits per character
        uint64_t bpc;
        // Mask for one character
        uint64_t character_mask;
        // Alphabet size
        uint64_t d;
        // Big Power
        uint64_t h;
        // Hashed rollers
        uint64_t hashed_count;

        // Constructor
        RollingHasher1(uint64_t q, uint64_t m);
        // Destructor
        ~RollingHasher1(){}
        // Update rolling hash with only the incoming character
        uint64_t update_rolling_hash_in(uint64_t in);
        // THIS MIGHT HAVE BEEN BUGGY BUT IS FIXED NOW?
        // Update rolling hash with both incoming and outgoing characters
        uint64_t update_rolling_hash_in_and_out(uint64_t in, uint64_t out);
        // Update rolling hash (UNIVERSAL)
        uint64_t update_rolling_hash(uint64_t in, uint64_t out);
        // Return the current hash value
        uint64_t get_current_hash();
        // Reset hasher state
        void reset();
        // Load full contents from k-mer factory
        void load_full_factory(KMerFactory2BC* kmer_factory);
        // Load full contents from canonical k-mer factory
        void load_full_factory_canonical(KMerFactoryCanonical2BC* kmer_factory);

};

class AdderHasher1
{   
    private:

        uint64_t hash_table_slots;
    
    public:

        AdderHasher1(uint64_t slots);

        ~AdderHasher1(){};

        uint64_t calculate_hash(KMerFactoryCanonical2BC * kmer_factory);

};

class ProbeHasher1
{
    public:

        ProbeHasher1(){};

        ~ProbeHasher1(){};

        // Double hahsing probing
        uint64_t probe_1(uint64_t item, uint64_t M);

        // Quadratic probing
        uint64_t probe_2(uint64_t iteration);

};