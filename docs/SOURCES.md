# Sources and mathematical dependencies

The sources below identify the statements used in the manuscript and the
hypotheses retained when applying them. Fixed preprint versions identify
the source text; they do not identify a later version of record. Third-party
PDFs are not redistributed by this repository.

| Source | Locator, hypotheses and use |
|---|---|
| Zhi-Hong Sun, *Supercongruences involving Apéry-like numbers and binomial coefficients*, AIMS Mathematics **7**(2) (2022), 2729–2781, [DOI 10.3934/math.2022153](https://doi.org/10.3934/math.2022153) | Definitions and recurrence parameters for f, a and Q, pp.2729–2730; binomial transforms, Lemma 3.1, p.2736. Theorem 3.4, pp.2738–2739, gives the modulo-p central-binomial transforms used here at 20 and −16. Conjecture 3.7, p.2742 (PDF page 14), is the target: primes p≠2,5, denominators 20 and −16, common sign (−1)^((p−1)/2), and modulus p². |
| Zhi-Hong Sun, *Congruences concerning Legendre polynomials II*, Journal of Number Theory **133** (2013), 1950–1976, [DOI 10.1016/j.jnt.2012.11.004](https://doi.org/10.1016/j.jnt.2012.11.004); text and numbering checked in [arXiv:1012.3898v2](https://arxiv.org/abs/1012.3898v2) | Theorem 2.1(ii), p.6, gives the same cubic character sum; equation (4.1) in the proof of Lemma 4.2, p.16, explicitly extends it to a polynomial parameter. Coefficient extraction gives exactly the Hasse–Legendre identity derived in Section 7. Theorem 4.2, p.17, at m=−1024 gives R_p=P_floor(p/4)(√5/2)² modulo p, and R_p=0 modulo p² when that Legendre value vanishes; Theorem 4.1 is a polynomial congruence supporting specialization over extensions. Theorem 5.3, p.19 of this preprint, gives the principal split value of R_p modulo p. These are antecedents, with the identities needed here proved locally. |
| Zhi-Hong Sun, *Congruences for sums involving Franel numbers*, International Journal of Number Theory **14** (2018), 123–142, [DOI 10.1142/S1793042118500094](https://doi.org/10.1142/S1793042118500094), [author's text](https://maths.hytc.edu.cn/__local/3/45/77/77F7497DDFD531E8688ED13C0D0_262306C8_24EB4.pdf) | Theorem 2.10, PDF p.9, and its proof, p.10: F_(16,p)=4x² modulo p when p=x²+5y² and p≡1,9 modulo 20. This is a direct antecedent for the principal split modulo-p value, not a modulo-p² evaluation of the two target sums. |
| Frits Beukers, *Supercongruences using modular forms*, [arXiv:2403.03301v3](https://arxiv.org/abs/2403.03301v3) | Lemma 3.7, pp.17–18, gives the weight-two Hecke-orbit permutation for an Atkin–Lehner extension of Γ₀(N), p∤N, and a character trivial on Γ₀(N); here N=6. Theorem 6.1, p.24, is the antecedent of the degree-p correction (g_p−g_1)t^p. Its Assumption (II), p.7, has r=1/3 in place of part (4), and its modulus is p³. |
| Armin Straub, *Gessel–Lucas congruences for sporadic sequences*, [arXiv:2301.12248v1](https://arxiv.org/abs/2301.12248v1) | Equation (10) and Theorem 3.2, p.5: A_F(p^r n)≡A_F(p^(r−1)n) modulo p^(2r), for p≥3 and positive n,r. The explicit conversion is A_F(n)=(−1)^n Q_n. It supplies the individual coefficient endpoint, not a value of either target sum. |
| Jeremy Rouse and John J. Webb, *On spaces of modular forms spanned by eta-quotients*, [arXiv:1311.1460v3](https://arxiv.org/abs/1311.1460v3) | Newman modularity conditions and Ligozat's cusp-order formula, p.2. The manuscript checks the eta exponent congruences, square condition and cusp orders before using holomorphy. |
| William Craig, *New Types of Sturm bounds via p-adic transfer methods*, [arXiv:2602.10240v1](https://arxiv.org/abs/2602.10240v1) | Only the classical bound recalled in Theorem 1.1(1), p.2: a holomorphic weight-k form vanishing to order greater than k[SL₂(Z):Γ]/12 is zero. The new transfer results are not used. |
| NIST Digital Library of Mathematical Functions, version 1.2.8, [§23.18, equations 23.18.5–7](https://dlmf.nist.gov/23.18) | Dedekind eta transformation, multiplier and principal square-root convention. These fix exact cusp phases and the eta inversion used at the common CM point. |
| Andrew V. Sutherland, *18.783 Elliptic Curves*, Fall 2023, [Lecture 4](https://math.mit.edu/classes/18.783/2023/LectureNotes4.pdf) | Theorem 4.15, p.4, extends a rational map from a smooth projective curve to a projective curve; Definition 4.20 and its preceding explanation, p.5, identify a nonconstant map of elliptic curves preserving O as an isogeny. Smoothness, nonconstancy and the image of O are checked for the displayed map. |
| Andrew V. Sutherland, *18.783 Elliptic Curves*, Fall 2023, [Lecture 5](https://math.mit.edu/classes/18.783/2023/LectureNotes5.pdf) | §§5.5–5.6, equation (4) and Theorem 5.21, pp.10–11, give the division-polynomial formulas for both coordinates of [5]. The short Weierstrass formulas apply in characteristic different from 2 and 3; the application uses p>5. |
| Andrew V. Sutherland, *18.783 Elliptic Curves*, Fall 2023, [Lecture 13](https://math.mit.edu/classes/18.783/2023/LectureNotes13.pdf) | p.1: for an ordinary elliptic curve in characteristic p, the geometric p-torsion points form a cyclic group of order p. Ordinarity and supersingularity are invariant under field extension. The proof uses geometric points even when the explicit endomorphism is defined over an extension. |
| Christophe Debry, *Beyond two criteria for supersingularity: coefficients of division polynomials*, [arXiv:1303.5002v1](https://arxiv.org/abs/1303.5002v1) | Deuring criterion, p.1: a smooth cubic model in odd characteristic is supersingular exactly when the coefficient of X^(p−1) in f(X)^((p−1)/2) vanishes. The printed statement assumes the cubic splits; the manuscript permits a finite field extension to ensure this. |

## Imported statements and local arguments

The Hecke-orbit lemma is an import. The integral coefficient descent, pole
bound, finite polynomial comparison and specialization are established for
the two periods in the manuscript. The degree-p step follows Beukers's
Theorem 6.1, but the ratio here is r=1/2 and the required modulus is p².
That adaptation is proved locally using separate endpoint congruences;
the modulo-p³ theorem is not invoked under different hypotheses.

The inert argument uses fixed polynomial identities for the explicit
degree-five map and both coordinates of its square. Its reduction for every
p>5 follows from the stated denominators, discriminant and resultant.
The finite Hasse–Legendre identity and the square of a Legendre polynomial
are proved in the text; the corresponding Hasse and quartic reductions
already occur in Sun's Legendre-polynomial paper cited above, including
the polynomial extension of the parameter. No analytic recognition of an elliptic j-value or
reduction theorem for an abstract CM endomorphism is needed.

The transport lemma requires p-integral parameters and coefficients, with
1−4ct a unit. The prime 3 is evaluated directly from the defining sums.
The transport and Dixon reductions have only modulo-p precision.
Consequently Sun's stronger vanishing result for the quartic sum does not
by itself lift either original sum modulo p². Each original sum is lifted
separately using the finite invariant. The cited older results do not
supply that missing precision for the two target sums.
The conclusions of the neighboring Conjectures 3.5 and 3.6 are not inputs.
The invariant and Hecke organization is also used in the author's
[Conjecture 3.5](https://github.com/aconsciousfractal/Sun-3-5-Supercongruence)
and [Conjecture 3.6](https://github.com/aconsciousfractal/Sun-3-6-Supercongruence)
treatments; the present manuscript proves the statements it needs.

## Related identities

The following comparisons use full texts, with a targeted reading of the
listed passages and a search of the complete extracted text. They are
related-work comparisons, not additional imports into the proof.

| Work | Passages compared and relation |
|---|---|
| Li-Li Cui and Zhi-Hong Sun, *Curious identities involving Legendre polynomials and Apéry-like numbers*, [arXiv:2607.12330v1](https://arxiv.org/abs/2607.12330v1), [DOI](https://doi.org/10.48550/arXiv.2607.12330) | Introduction, pp.1–3, Theorem 2.1 and its proof, pp.4–5, and Corollaries 2.1–2.2, p.6: weighted Legendre products and applications to first-kind Apéry-like sequences. These concern sums over polynomial degrees, rather than the individual Legendre square used here. |
| Yu-Tian Li and Zhi-Hong Sun, *Weighted bilinear identities and supercongruences for Apéry-like polynomials*, [arXiv:2609.28098v1](https://arxiv.org/abs/2609.28098v1), [DOI](https://doi.org/10.48550/arXiv.2609.28098) | Introduction and start of Section 2, pp.1–4: weighted products of g_n(x) and v_n(x), with endpoint expansions. The families and weighted square sums differ from the two central-binomial sums evaluated here. |
| Zhi-Hong Sun, *Identities and congruences involving orthogonal polynomials and Apéry-like numbers*, [arXiv:2609.29910v1](https://arxiv.org/abs/2609.29910v1), [DOI](https://doi.org/10.48550/arXiv.2609.29910) | pp.1–5, Theorem 3.1, pp.9–10, and Corollaries 3.3–3.4, p.11: the same f, a and Q definitions and recurrences occur in a general quadratic telescoping framework. The explicit a_n² and Q_n² sums at 81 and 64 differ from the target sums. |

These targeted passages do not identify a proof of the precise target
congruence. This observation does not exclude other equivalent formulations,
other passages or other antecedents. No exhaustive priority review is claimed.
