#include <cmath>
#include <cstdint>

#pragma once

//////////////////////////////////////////////
//
// This file has some useful math functions
//
//////////////////////////////////////////////

namespace mathfunctions
{
    /*
        This function returns the smallest prime 
        that is at least as large as the given integer

        Useful for finding hash table size 
        (maybe a prime at elast 1.33 times the expected amount of stored items)
    */
    uint64_t next_prime(uint64_t at_least);

}