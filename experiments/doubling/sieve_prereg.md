# Independent factor-blind exact-residue check

Before running a direct large-index evaluation, take the first two registered 32-bit N in each class. With only those four N as inputs, independently sieve primes up to the largest 2m and multiply prime powers from Legendre's factorial valuation formula modulo each N. Do not supply the factors, call a factoring routine, or use CRT. Require four exact matches including both nonzero residues. This supplements the factor-first exact computation at 32, 48 and 64 bits.
