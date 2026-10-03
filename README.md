# Sun's supercongruence for two Apéry-like binomial sums at discriminant −20

This repository accompanies the proof of Zhi-Hong Sun's Conjecture 3.7 in
*Supercongruences involving Apéry-like numbers and binomial coefficients*,
AIMS Mathematics **7** (2022), 2729–2781,
[DOI 10.3934/math.2022153](https://doi.org/10.3934/math.2022153).

Define

$$
f_n=\sum_{j=0}^n\binom nj^3,\qquad
a_n=\sum_{j=0}^n\binom nj^2\binom{2j}j,\qquad
Q_n=\sum_{j=0}^n\binom nj(-8)^{n-j}f_j.
$$

For every prime $p\ne2,5$, put $\varepsilon_p=(-1)^{(p-1)/2$. The manuscript proves

$$
\varepsilon_p\sum_{n=0}^{p-1}\binom{2n}n\frac{a_n}{20^n}
\equiv
\varepsilon_p\sum_{n=0}^{p-1}\binom{2n}n\frac{Q_n}{(-16)^n}
\equiv
\begin{cases}
4x^2-2p,&p=x^2+5y^2,\ p\equiv1,9\pmod{20},\\
2x^2-2p,&2p=x^2+5y^2,\ p\equiv3,7\pmod{20},\\
0,&p\equiv11,13,17,19\pmod{20}
\end{cases}\pmod{p^2}.
$$

The two sums are evaluated separately. For split primes, two modular periods
are specialized at a common complex multiplication point; an independent
coefficient congruence handles degree $p$ in the Hecke-polynomial argument.
For inert primes, an explicit degree-five isogeny with square $[-5]$
forces the Hasse invariant to vanish. Finite binomial identities and a
quadratic differential invariant then give the congruences modulo $p^2$.
The prime $3$, both complex multiplication classes and all denominators
are treated explicitly.

## Read and reproduce

- [Manuscript PDF](paper/Sun_3_7_Supercongruence.pdf) and
  [standalone LaTeX source](paper/Sun_3_7_Supercongruence.tex).
- [Reader's guide](docs/READING_GUIDE.md),
  [claim and evidence map](docs/CLAIM_MAP.md) and [sources](docs/SOURCES.md).
- [Reproduction instructions](REPRODUCE.md),
  [companion](companion/README.md) and [evidence guide](docs/EVIDENCE.md).

From the repository root, with Python 3.10 or later:

~~~text
python -B scripts/verify.py --core-only
python -B -m unittest discover -s tests -v
python -B scripts/check_manifest.py
~~~

Core verification uses only the standard library. The full replay of all
19 evidence files additionally uses SymPy 1.14.0; see REPRODUCE.
Verification runs locally without a network connection.

## Scope and identity

The 301-prime regression uses recurrences; direct evaluation from the
defining sums covers 44 primes. Neither finite panel proves the all-prime
statement. Modular coefficient certificates require the proved holomorphy
and valence bound; the fixed isogeny identities require the stated good
reduction argument. The equality of the two sums modulo $p^3$ already fails
at $p=7$.

The work has not undergone independent specialist human review. No claim of
bibliographic priority is made. See [AI assistance](AI_USE.md).

Version: **0.1.0**. EVIDENCE_SHA256.txt identifies the computational core;
its digest is printed in the manuscript. MANIFEST_SHA256.txt identifies
the complete distribution. Neither inventory is a signature, and a package
version does not denote a journal version of record.

Companion repository address:
[aconsciousfractal/Sun-3-7-Supercongruence](https://github.com/aconsciousfractal/Sun-3-7-Supercongruence).

Code, certificate data and documentation use [MIT](LICENSE); the manuscript
uses [CC BY 4.0](LICENSE_MANUSCRIPT.md). See [licence scope](LICENSE_SCOPE.md)
and [citation metadata](CITATION.cff).
