"""lem85: exact checks of the elementary identities used in the proof of Lemma 8.5' (phase 11).

See notes/lem85_results.md. Module `identities` verifies, by exact integer/rational arithmetic on toy
sizes, the four elementary steps of the proof:
- (E2) lift averaging: averaging E(x,y;ls,.) over the l-1 lifts of a class mod s gives E(x,y;s,.)/(l-1)
  when l is a prime > y and l does not divide s (this is the "chi_l principal" part, with no characters);
- (E3) the Moebius expansion of a character sum over y-smooth n coprime to a modulus (real characters);
- (E4) the unique factorisation n = m k of a y-smooth n > M1 with M1 < m <= P+(m) M1, P+(m) <= P-(k);
- (E5) sum_{j<=X} 1/phi(j) <= 2(1 + log X) on a range of X (log only for the comparison).
These are correctness checks of statements that are also proved in the note. Nothing here is analytic.
"""
