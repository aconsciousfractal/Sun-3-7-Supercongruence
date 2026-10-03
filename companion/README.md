# Exact companion

Run `python -B scripts/verify.py` from the repository root to reconstruct and
compare all 19 delivered evidence files. Full mode explicitly requires
SymPy 1.14.0. `--core-only` reconstructs 17 files using the standard library
and names the two omitted symbolic files. Verification does not update
delivered evidence. See [reproduction instructions](../REPRODUCE.md).

The constructors use exact integers, rational numbers, modular inverses and
algebraic polynomial domains. Expected evidence is read only after each
constructor has completed; it never supplies input coefficients or target
values to the constructions. All checks use explicit exceptions and remain
active under `python -O`.

| Module | Role |
| --- | --- |
| `core.py` | Integer series, exact eta products, recurrences and modular arithmetic shared by three constructors. |
| `check_arithmetic.py` | Recurrence regression for 301 primes through 1999; direct definitions for 597 sequence values; CM matrices, finite invariants and Hasse bridges. |
| `check_modular.py` | Two eta expansion algorithms through degree 24 and six Hecke polynomials for primes 7, 11 and 29. |
| `check_algebra.py` | 35 exact symbolic identities and a covariance corruption control; requires SymPy. |
| `exact_isogeny.py` | Four full polynomial identities over Q[s]/(s²−5), nondegeneration data and three corruption controls. |
| `finite_audit.py` | A separate implementation of Hasse and Legendre polynomials, curve point checks, and the original binomial sums for 44 primes through 199. |
| `symbolic_audit.py` | A second implementation of the four isogeny identities and three auxiliary identities; requires SymPy. |
| `evidence.py` | Strict file inventory and byte comparison, with no arithmetic input from frozen results. |

The rational-pair isogeny engine and its symbolic counterpart check complete
polynomial identities. The finite-field curve checks, prime regressions and
bounded original-sum evaluations have their stated finite scopes. Shared
mathematical input remains shared even when implementations differ. The
manuscript supplies the uniform arguments and hypotheses needed for the
theorem for all admissible primes.

JSON output is UTF-8 with LF. The CSV has UTF-8 encoding and explicit CRLF
record terminators. Output is written to a temporary directory with automatic
cleanup unless `--out-dir` is supplied. The CLI rejects paths equal to,
above, or below the source or expected-evidence directory before creating
files, and refuses to overwrite any nonempty directory.
