#include "functions_math.hpp"


namespace mathfunctions
{
    
    uint64_t next_prime(uint64_t at_least)
    {
        // The prime candidate is here
        uint64_t prime_candidate = at_least;
        // Check if the given number is less than 2
        if (prime_candidate <= 2)
            return 2;
        // If the given number is even it cannot be a prime so we make it odd
        if (prime_candidate % 2 == 0)
            prime_candidate += 1;
        // Define few useful variables
        uint64_t max_check;
        uint64_t check;
        bool is_prime;
        // Run this loop until a prime is found
        while(true)
        {
            // Prime divisibility checking starts at 3
            check = 3;
            // Prime divisibility checking ends at the square root
            max_check = std::floor(std::sqrt(prime_candidate));
            is_prime=true;
            //std::cout << "Checking for prime "<< prime_candidate << " with max check " << max_check << "\n";
            // Check divisibility for all relevant values
            while(check <= max_check)
            {
                // If divisible, no prime
                if (prime_candidate % check == 0)
                {
                    is_prime = false;
                    break;
                }
                check+=1;
            }
            // When we find a prime, rreturn it
            if (is_prime)
                return prime_candidate;
            // If candidate is not prime, increase it by 2 (skip the even value)
            prime_candidate+=2;
        }
    }
}