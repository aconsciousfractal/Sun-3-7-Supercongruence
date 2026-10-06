# Reproduction

Run commands from the repository root. Python 3.10 or later is required.
The delivered PDF is independent of the Python verification environment.

## Exact companion

The standard-library subset reconstructs 17 evidence files:

~~~text
python -B scripts/verify.py --core-only
python -B -O scripts/verify.py --core-only
~~~

To reconstruct all 19 files, use an environment with the optional symbolic
dependency installed:

~~~text
python -m pip install -r requirements-symbolic.txt
python -B scripts/verify.py
python -B -O scripts/verify.py
~~~

The full mode explicitly requires SymPy and fails if it is unavailable.
The core-only mode names the two omitted symbolic files; its success does
not claim they were checked. Constructors calculate evidence independently
before comparison with the delivered files. Expected data never provides
input coefficients.

Temporary output is cleaned automatically. To retain regenerated evidence,
pass --out-dir followed by a new or empty directory outside the repository.
The verifier refuses output equal to, inside, or above the repository or
expected evidence directory before creating files. It does not refresh
delivered evidence. The default expected directory is evidence/.

JSON is serialized deterministically with UTF-8 and LF; the CSV uses UTF-8
and the explicit CRLF record terminator of the standard CSV dialect.
Comparison checks exact bytes and the complete declared inventory, including
missing, extra or corrupted files. Core-only mode still requires a complete
19-file expected inventory, but compares only its 17 reconstructed files.
See [the evidence guide](docs/EVIDENCE.md) for individual roles and scopes.

Run the regression and corruption tests and verify both inventories:

~~~text
python -B -m unittest discover -s tests -v
python -B -O -m unittest discover -s tests -v
python -B scripts/check_manifest.py --evidence
python -B scripts/check_manifest.py
~~~

The arithmetic regression covers every prime from 3 through 1999 except 5.
Direct evaluation from the defining sums covers the 44 admissible primes
through 199. Other checks have their own stated finite scopes. Prime
regression is evidence of implementation agreement, not a proof for all primes.

## Paper

Use an existing pdflatex installation with the packages in the TeX preamble.
The build script does not install TeX or packages; it disables shell escape
and automatic MiKTeX package installation.

~~~text
python -B scripts/build_pdf.py
~~~

The rebuilt PDF and build_info.json go to ignored paper/build/. The delivered
PDF is replaced only by the explicit authoring option --update. Within the
repository, custom build output is restricted to paper/build/ and its
subdirectories. A path outside the repository is also allowed.

The source suppresses volatile dates and IDs. The builder fixes
SOURCE_DATE_EPOCH and performs three passes. Identical bytes require the
recorded engine, package versions and fonts; different installations may
produce an equivalent PDF with different bytes. See [build details](paper/BUILD.md).

## Inventories and source archive

EVIDENCE_SHA256.txt covers the companion, verification code, tests, exact
JSON/CSV evidence and declared verification configuration. Its digest is
printed in the paper. Manuscript and explanatory prose are excluded from
that identity to avoid a circular hash reference. MANIFEST_SHA256.txt
covers the complete distribution except itself. Neither is a signature.

After deliberate changes, the following are explicit authoring commands:

~~~text
python -B scripts/check_manifest.py --evidence --write
python -B scripts/check_manifest.py --write
~~~

When the computational core changes, update its digest in the paper and
rebuild before refreshing the full manifest. Ordinary verification never
rewrites evidence or inventories.

~~~text
python -B scripts/make_archive.py --output <output-directory>/Sun-3-7-Supercongruence-0.1.1.zip
~~~

Replace the placeholder with a writable directory outside the repository
and quote paths containing spaces. The destination must not already exist.
The ZIP contains precisely the verified manifest members, with deterministic
entry metadata. Git state, caches and build products are excluded. An
extracted copy supports the same checks without Git.

Hosted CI is configured for core and full checks. Its configuration is not
evidence that a hosted run has occurred.
