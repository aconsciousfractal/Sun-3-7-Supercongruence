"""Strict inventory and byte comparison; frozen evidence is never an input oracle."""
from pathlib import Path

ARITHMETIC_FILES = (
    'original_sequences.json', 'prime_regression.csv', 'CM_matrices.json',
    'finite_Hasse_bridges.json', 'finite_invariants.json', 'negative_controls.json',
    'arithmetic_summary.json',
)
MODULAR_FILES = ('modular_certificates.json', 'hecke_polynomials.json',
                 'modular_summary.json')
ISOGENY_FILES = ('isogeny_certificate.json', 'isogeny_summary.json')
FINITE_FILES = ('finite_summary.json', 'Hasse_polynomials.json',
                'finite_curves.json', 'exact_Legendre_squares.json', 'original_sums.json')
CORE_FILES = ARITHMETIC_FILES + MODULAR_FILES + ISOGENY_FILES + FINITE_FILES
SYMBOLIC_FILES = ('algebra.json', 'symbolic_summary.json')
ALL_FILES = CORE_FILES + SYMBOLIC_FILES


def inventory(directory, names, label):
    """Require exactly the named regular files in one flat directory."""
    directory = Path(directory)
    if not directory.is_dir():
        raise ValueError(label+' must be a directory')
    members = list(directory.iterdir())
    if any(not path.is_file() or path.is_symlink() for path in members):
        raise ArithmeticError(label+' must contain only regular evidence files')
    actual = {path.name for path in members}
    expected = set(names)
    if actual != expected:
        missing = sorted(expected-actual)
        extra = sorted(actual-expected)
        raise ArithmeticError(f'{label} inventory mismatch: missing={missing}; extra={extra}')


def compare_expected(actual, expected, core_only=False):
    """Compare reconstructed bytes; core mode explicitly omits two symbolic files."""
    actual = Path(actual).resolve()
    expected = Path(expected).resolve()
    if actual == expected or actual.is_relative_to(expected) or expected.is_relative_to(actual):
        raise ValueError('reconstructed output and expected evidence must be disjoint')
    names = CORE_FILES if core_only else ALL_FILES
    inventory(actual, names, 'reconstructed output')
    inventory(expected, ALL_FILES, 'expected evidence')
    for name in names:
        if (actual/name).read_bytes() != (expected/name).read_bytes():
            raise ArithmeticError('evidence byte mismatch: '+name)
    return len(names)
