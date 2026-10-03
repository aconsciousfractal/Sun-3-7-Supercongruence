# Evidence roles and scope

The all-prime result is a written proof with explicitly attributed imports.
The companion's exact calculations have distinct roles. All 19 files below
are delivered directly in `evidence/`; the verifier reconstructs them before
comparing their bytes with the delivered files.

| Delivered file | Role and finite scope |
|---|---|
| original_sequences.json | The 597 values of f, a and Q from their defining sums, for n=0,…,198. The constructor compares these values with the recurrence coefficients. |
| prime_regression.csv | Prime, branch, representation, sign, two raw target sums, Franel-sum value and signed target for the 301 admissible primes through 1999, using recurrence coefficients. The panel has 70 principal split, 79 nonprincipal split and 152 inert primes. |
| CM_matrices.json | 148 exact matrices for the split primes p>5 through 1999. The prime 3 has its own direct treatment in the proof. |
| finite_invariants.json | 129 check records with coefficient counts and two resonance booleans for the three sequences at 43 primes p>5 through 199. The constructor checks the complete invariant polynomials, comprising 25,431 coefficients and 258 resonant equalities; their coefficient vectors are not stored here. |
| finite_Hasse_bridges.json | Hasse and Legendre coefficient vectors, finite sum values and Vandermonde row counts at the same 43 primes. The constructor checks 2,130 integer Vandermonde rows. These calculations test the finite reductions and their support. |
| hecke_polynomials.json | Six exact reconstructions at p=7,11,29 for the two periods, containing 660 integer coefficients and four exact CM factors. The constructors do not use the desired CM roots to reconstruct the polynomials. |
| modular_certificates.json | Eta and cusp data, series and residual vectors through q²⁴, and six corruption-control outcomes. The constructor compares two eta expansion methods. The proof gives the shorter sufficient bounds. |
| modular_summary.json | Parameters and outcomes of the modular calculations. It summarizes those calculations and is not another experiment. |
| negative_controls.json | Deliberately altered arithmetic, character, normalization and weight data, together with failures of equality at modulus p³, including p=7. |
| arithmetic_summary.json | Parameters and outcomes of the recurrence, invariant, transform, matrix and prime calculations. |
| exact_Legendre_squares.json | Twenty-one degree, coefficient-count and success records for n=0,…,20. The constructor compares 441 coefficients in exact polynomial square identities over the rationals; these vectors are not stored here. |
| Hasse_polynomials.json | Hasse and Legendre coefficient vectors and their evaluations for 43 primes p>5 through 199. The Hasse vectors contain 1,075 coefficients. |
| finite_curves.json | Counts, field parameters and traces for 22 curves at eleven primes p>5 through 43. The constructor checks 16,988 point compositions and 264 sampled homomorphism pairs; individual points are not stored. |
| original_sums.json | The prime, sign and two raw sum values at all 44 admissible primes through 199, recomputed from the defining binomial sums. |
| isogeny_certificate.json | Fixed input and intermediate polynomials over Q[s]/(s²−5), in ascending rational coefficient pairs. These permit reconstruction of the four map identities. |
| isogeny_summary.json | Four zero-residual records with coefficient counts totalling 108 rational pairs, plus denominator, discriminant and resultant data and three corruption controls. The constructor expands all four identities checking the curve equation and both coordinates of Φ²=[−5]. |
| finite_summary.json | Counts and outcomes of the finite Legendre, Hasse, curve and original-sum calculations, including 1,120 successful parameter-shift checks at inert primes through 97 and a transport-domain counterexample. |
| algebra.json | SymPy calculation: 35 named zero-residual records, exact CM and conjugation data, and a nonzero residual for an incorrect covariance formula. |
| symbolic_summary.json | SymPy calculation: seven named zero-residual records, comprising the four map identities and three auxiliary identities, plus the resultant value. |

## Exact certificates and finite panels

The modular residuals are holomorphic modular forms before any coefficient
bound is applied. The four weight-eight identities each have bound 8;
checking q⁰ through q⁸ gives 36 coefficients in total. The weight-two
identities use their own bound 2, and the rational eta relation uses its
pole-degree argument. The longer q²⁴ windows are consistency checks.
Coefficient agreement alone would not establish modularity or cusp holomorphy.
The manuscript uses 48 residual coefficients: 36 for weight eight, six
for weight two, and six for the rational eta relation. The companion also
checks two coefficients of an auxiliary E₄ identity, giving 50 in total;
that identity is not an input to the proof.

The isogeny certificate is a set of fixed identities in characteristic zero.
The written good-reduction argument checks denominators, discriminant and
resultant, so those identities may be reduced for every p>5. The finite point
counts test their implementation; sampling points does not establish that
reduction argument or the map identity for all characteristics.

The generic invariant, resonance, multiplicity and lifting arguments are
in the manuscript. The invariant and local-lift panels test instances.
The Hecke examples test normalization and reconstruction; the integral
descent and degree bound for arbitrary p are written arguments. Likewise,
the finite prime panels cannot prove the theorem for all primes.

## Reproduction and identity

Full mode, `python -B scripts/verify.py`, reconstructs all 19 files and
requires the declared SymPy dependency. `--core-only` reconstructs 17 files
using the standard library and explicitly omits `algebra.json` and
`symbolic_summary.json`. Both modes require the complete 19-file expected
inventory; subset success applies only to the files reconstructed.
See [REPRODUCE](../REPRODUCE.md) for commands and serialization details.

`EVIDENCE_SHA256.txt` identifies the computational core, including the
verification code, tests and delivered evidence. The manuscript prints its
digest. Manuscript and explanatory prose are excluded to avoid a circular
hash reference. `MANIFEST_SHA256.txt` covers the complete distribution.
A digest establishes byte identity, not mathematical truth, novelty or
independent specialist acceptance.
