# Delivered PDF build

The version 0.1.0 manuscript has 18 A4 pages, nine main sections and three
appendices. Theorem 1.1 begins on page 1. Page numbers are the only running
text; there is no table of contents or repeated author/title header.
All 20 used font subsets are embedded.

The delivered build used MiKTeX-pdfTeX 4.23 (MiKTeX 25.12), pdfTeX 1.40.28,
LaTeX2e 2025-11-01 and L3 2025-12-24. Three passes completed with no compiler,
package, reference, duplicate-destination or overfull/underfull warnings.
All 18 rendered pages were visually inspected.

| Artifact | SHA-256 |
| --- | --- |
| TeX source | `54f42bca484d5ca181b419eb4af7c18197605fcba27c36103e8115f8c2198a88` |
| Delivered PDF | `b91bf8716bc9c881836b7f973767e03054733dd41005ee882cd44967ce2df952` |
| EVIDENCE_SHA256.txt | `b5d9762b681919932cbf1c54e0c9296bdd42fe09d44dc8450215b40079277219` |

Run from the repository root with an existing TeX installation:

~~~text
python -B scripts/build_pdf.py
~~~

The script finds `pdflatex` on PATH, or uses the `PDFLATEX` environment
variable. It writes to ignored `paper/build/`, performs three passes and
reports hashes in `build_info.json`. `--update` explicitly replaces the
delivered PDF. It does not install anything and disables shell escape and
automatic MiKTeX package installation.

The build fixes `SOURCE_DATE_EPOCH=1790985600`, `FORCE_SOURCE_DATE=1` and
`TZ=UTC`; the source suppresses volatile dates and identifiers. Byte identity
requires the same TeX engine, packages and fonts. Another installation can
produce an equivalent PDF with different bytes. The package versions below
are the observed build environment, not an automatic dependency installer
or a complete operating-system image.

~~~text
Document Class: article 2025/01/22 v1.4n Standard LaTeX document class
Package: fontenc 2025/07/18 v2.1d Standard LaTeX package
Package: lmodern 2015/05/01 v1.6.1 Latin Modern Fonts
Package: microtype 2026/03/01 v3.2d Micro-typographical refinements (RS)
Package: keyval 2022/05/29 v1.15 key=value parser (DPC)
Package: etoolbox 2025/10/02 v2.5m e-TeX tools for LaTeX (JAW)
Package: geometry 2026/03/07 v6.0 Page Geometry
Package: ifvtex 2019/10/25 v1.7 ifvtex legacy package. Use iftex instead.
Package: iftex 2024/12/12 v1.0g TeX engine tests
Package: amsmath 2025/07/09 v2.17z AMS math features
Package: amstext 2024/11/17 v2.01 AMS text
Package: amsbsy 1999/11/29 v1.2d Bold Symbols
Package: amsopn 2022/04/08 v2.04 operator names
Package: amssymb 2013/01/14 v3.01 AMS font symbols
Package: amsfonts 2013/01/14 v3.01 Basic AMSFonts support
Package: amsthm 2020/05/29 v2.20.6
Package: mathrsfs 1996/01/01 Math RSFS package v1.0 (jk)
Package: booktabs 2020/01/12 v1.61803398 Publication quality tables
Package: array 2025/09/25 v2.6n Tabular extension package (FMi)
Package: hyperref 2025-07-12 v7.01o Hypertext links for LaTeX
Package: kvsetkeys 2022-10-05 v1.19 Key value parser (HO)
Package: kvdefinekeys 2019-12-19 v1.6 Define keys (HO)
Package: pdfescape 2019/12/09 v1.15 Implements pdfTeX's escape features (HO)
Package: ltxcmds 2023-12-04 v1.26 LaTeX kernel commands for general use (HO)
Package: pdftexcmds 2020-06-27 v0.33 Utility functions of pdfTeX for LuaTeX (HO
Package: infwarerr 2019/12/03 v1.5 Providing info/warning/error messages (HO)
Package: hycolor 2020-01-27 v1.10 Color options for hyperref/bookmark (HO)
Package: nameref 2025-06-21 v2.57 Cross-referencing by name of section
Package: refcount 2019/12/15 v3.6 Data extraction from label references (HO)
Package: gettitlestring 2019/12/15 v1.6 Cleanup title references (HO)
Package: kvoptions 2022-06-15 v3.15 Key value format for package options (HO)
Package: stringenc 2019/11/29 v1.12 Convert strings between diff. encodings (HO
Package: intcalc 2019/12/15 v1.3 Expandable calculations with integers (HO)
Package: url 2013/09/16  ver 3.4  Verb mode for urls, etc.
Package: bitset 2019/12/09 v1.3 Handle bit-vector datatype (HO)
Package: bigintcalc 2019/12/15 v1.5 Expandable calculations on big integers (HO
Package: rerunfilecheck 2025-06-21 v1.11 Rerun checks for auxiliary files (HO)
Package: uniquecounter 2019/12/15 v1.4 Provide unlimited unique counter (HO)
~~~

Fonts are Latin Modern, AMS msbm10 and RSFS rsfs10. See
[reproduction instructions](../REPRODUCE.md) for the exact companion,
source archive and inventory checks.
