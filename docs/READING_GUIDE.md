# Reading the proof and companion

Start with Theorem `thm:main` and the definitions in the
[manuscript](../paper/Sun_3_7_Supercongruence.pdf). Both sums are truncated at
p−1, with denominators 20 and −16 and common sign (−1)^((p−1)/2).
The target is modulo p² for every prime p≠2,5. The two split classes of
discriminant −20 give different quadratic-form traces; the inert classes
give zero.

## A first pass

1. Read the main theorem and proof outline, then the finite invariant and
   local lifting argument (`lem:invariant`, `lem:lifting`).
2. For the split primes, follow the two periods, the endpoint congruence,
   finite Hecke specialization and exact CM matrices
   (`lem:period`, `lem:endpoint`, `thm:specialization`, `prop:split`).
3. For the inert primes, follow the finite transforms and Hasse polynomial
   through the explicit degree-five isogeny (`lem:transport`, `lem:dixon`,
   `lem:hasse`, `lem:legendre-square`, `thm:isogeny`, `prop:inert`).
4. Read the corollaries after the main proof. The split converse in
   `cor:factor` uses the main theorem; it is not a premise of that theorem.

Each sum receives its own evaluation. Their equality is a conclusion of
those evaluations. Consult the [claim map](CLAIM_MAP.md) and
[sources](SOURCES.md) for the role of each argument and import.

## Points to inspect closely

| Point | What the argument retains |
|---|---|
| Finite invariant | The resonant coefficients at degrees p and 2p are computed separately. |
| Local lifting | The parameter and differential-operator coefficient are units; the second derivative is a unit at the double zero. |
| Full truncation | All terms through p−1 are retained modulo p². |
| Modular identities | Cusp holomorphy and the valence bound precede the finite coefficient certificates. |
| Hecke coefficients | Integral descent and the geometric degree bound hold for arbitrary admissible p. |
| Degree p | A separate endpoint g_p≡g_1 eliminates the degree-p term; comparison only through p−1 would be insufficient. |
| Specialization | A finite polynomial identity is evaluated after its proof, using an exact unit root. |
| CM classes | Both determinant 1 and determinant 2 occur, with the latter's character sign. |
| Transport | Parameters and coefficients are p-integral and 1−4ct is a unit. |
| Explicit map | Both coordinates, including the sign of Y, establish Φ²=[−5]; denominator and resultant checks allow reduction for every p>5. |
| Torsion | The ordinary geometric p-torsion group is cyclic, even if Φ is defined over an extension. |
| Hasse | The identity is polynomial, so it may be evaluated over an extension; supersingularity does not mean zero integer trace over every extension. |
| Exceptional prime | p=3 is treated directly. The source excludes 2 and 5. |

## Computational pass

Run `python -B scripts/verify.py` for the full 19-file replay, using the
declared SymPy dependency, or `python -B scripts/verify.py --core-only`
for the 17-file standard-library subset. Follow [REPRODUCE](../REPRODUCE.md)
for tests, output comparison and inventory checks. The
[evidence guide](EVIDENCE.md) identifies each file's scope.

The modular coefficient and fixed isogeny certificates support specific
exact identities when their written hypotheses are verified. The prime
panels are finite implementation regressions. A successful replay cannot
replace the mathematical reading or establish bibliographic priority.
