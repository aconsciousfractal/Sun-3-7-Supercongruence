"""Reconstruct exact finite certificates and bounded regressions before byte comparison."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import tempfile
from uuid import uuid4

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT))

from companion.evidence import ALL_FILES, CORE_FILES, SYMBOLIC_FILES, compare_expected, inventory
from companion import check_arithmetic, check_modular, exact_isogeny, finite_audit


def validate_output_dir(directory, expected):
    """Resolve aliases and reject every containment relation before any write."""
    directory = Path(directory).resolve()
    expected = Path(expected).resolve()
    for protected, label in ((ROOT.resolve(), 'repository'), (expected, 'expected evidence')):
        if (directory == protected or directory.is_relative_to(protected)
                or protected.is_relative_to(directory)):
            raise ValueError('output must be disjoint from '+label+' (equal, ancestor or descendant refused)')
    if directory.exists() and (not directory.is_dir() or any(directory.iterdir())):
        raise ValueError('output must be a new or empty directory; existing files are never overwritten')
    return directory


def require_symbolic():
    try:
        import sympy
    except ImportError as exc:
        raise RuntimeError('Full verification requires SymPy 1.14.0 from requirements-symbolic.txt; '
                           '--core-only explicitly omits the two symbolic files') from exc
    if sympy.__version__ != '1.14.0':
        raise RuntimeError('Full verification requires SymPy 1.14.0; found '+sympy.__version__)


def run(directory, core_only=False):
    """Calculate every coefficient from definitions, recurrences or exact algebra."""
    arithmetic = check_arithmetic.run(directory)
    modular = check_modular.run(directory)
    isogeny = exact_isogeny.run(directory)
    finite = finite_audit.run(directory)
    symbolic_counts = None
    if not core_only:
        from companion import check_algebra, symbolic_audit
        algebra = check_algebra.run(directory)
        symbolic = symbolic_audit.run(directory)
        symbolic_counts = {'algebra_identities': algebra['identity_count'],
                           'isogeny_and_auxiliary_identities': len(symbolic['identities'])}
    return {'schema_version': 1, 'status': 'PASS',
            'scope': 'exact finite certificates and bounded regressions',
            'mode': 'core-only' if core_only else 'full',
            'reconstructed_files': len(CORE_FILES) if core_only else len(ALL_FILES),
            'omitted_files': list(SYMBOLIC_FILES) if core_only else [],
            'recurrence_prime_rows': arithmetic['primes'],
            'recurrence_prime_bound': arithmetic['max_prime'],
            'direct_sequence_values': arithmetic['direct_sequence_values'],
            'direct_original_prime_rows': finite['fresh_original_prime_rows'],
            'CM_matrices': arithmetic['CM_matrices'],
            'exact_isogeny_identities': len(isogeny['identities']),
            'exact_isogeny_coefficient_pairs': isogeny['residual_coefficient_pairs'],
            'Hecke_polynomials': modular['hecke_polynomials'],
            'finite_point_compositions': finite['finite_point_compositions'],
            'symbolic_counts': symbolic_counts}


def replay(directory, expected, core_only):
    summary = run(directory, core_only)
    summary['byte_identical_evidence_files'] = compare_expected(directory, expected, core_only)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path,
                        help='retain reconstructed evidence in a new or empty directory outside the repository')
    parser.add_argument('--expected-dir', type=Path, default=ROOT/'evidence',
                        help='complete 19-file expected evidence directory (default: evidence/)')
    parser.add_argument('--core-only', action='store_true',
                        help='explicitly omit algebra.json and symbolic_summary.json; use only the standard library')
    args = parser.parse_args()
    expected = args.expected_dir.resolve()
    if args.out_dir is not None:
        directory = validate_output_dir(args.out_dir, expected)
    else:
        # Validate the whole private container before TemporaryDirectory creates
        # its child. This also covers custom expected directories under the OS
        # temporary root without ever creating a directory inside them.
        container = validate_output_dir(Path(tempfile.gettempdir())/('sun37-replay-'+uuid4().hex), expected)
    inventory(expected, ALL_FILES, 'expected evidence')
    if not args.core_only:
        require_symbolic()
    if args.out_dir is not None:
        directory.mkdir(parents=True, exist_ok=True)
        summary = replay(directory, expected, args.core_only)
    else:
        container.mkdir()
        try:
            with tempfile.TemporaryDirectory(prefix='run-', dir=container) as scratch:
                summary = replay(Path(scratch), expected, args.core_only)
        finally:
            container.rmdir()
    print(json.dumps(summary, sort_keys=True, separators=(',', ':')))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, ArithmeticError, RuntimeError) as exc:
        print('SUN37_REPLAY_FAIL: '+str(exc), file=sys.stderr)
        raise SystemExit(1)
