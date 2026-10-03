"""Adverse checks of exact arithmetic, mathematical scope and evidence preservation."""
import importlib.util
import io
import json
from fractions import Fraction
from math import comb
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT))
from companion import core, check_arithmetic, check_modular, exact_isogeny, finite_audit
from companion.evidence import ALL_FILES, CORE_FILES, SYMBOLIC_FILES, compare_expected

spec = importlib.util.spec_from_file_location('sun37_verify', ROOT/'scripts/verify.py')
verify = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verify)


class ExactArithmeticTests(unittest.TestCase):
    def test_definitions_recurrences_and_rational_sums_are_distinct_oracles(self):
        top = 24
        f = [sum(comb(n, k)**3 for k in range(n+1)) for n in range(top+1)]
        a = [sum(comb(n, k)**2*comb(2*k, k) for k in range(n+1)) for n in range(top+1)]
        q = [sum(comb(n, k)*(-8)**(n-k)*f[k] for k in range(n+1)) for n in range(top+1)]
        for name, direct in (('f', f), ('a', a), ('Q', q)):
            self.assertEqual(direct, core.sequence(*core.PARAMS[name], top))
        for p in (3, 7, 11, 13, 19):
            for seq, t in ((a, Fraction(1, 20)), (q, Fraction(-1, 16))):
                total = sum((comb(2*n, n)*seq[n]*t**n for n in range(p)), Fraction())
                reduced = total.numerator*pow(total.denominator, -1, p*p) % (p*p)
                self.assertEqual(reduced, core.sum_mod(core.central(seq), p, t))
        with self.assertRaises(ValueError):
            core.residue(Fraction(1, 5), 25)

    def test_exact_arithmetic_rejects_nonintegral_and_float_inputs(self):
        with self.assertRaises(ArithmeticError):
            core.sequence(0, 1, 0, 2)
        with self.assertRaises(ArithmeticError):
            core.inv([2, 1], 3)
        with self.assertRaises(ArithmeticError):
            core.compose([1, 2], [1, 1], 3)
        for value in (True, 3.0, -1):
            with self.subTest(value=value), self.assertRaises(ValueError):
                core.sequence(10, 3, 9, value)
        for value in (True, 0.5, '1/2'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                core.residue(value, 49)
        with self.assertRaises(ValueError):
            core.mpow([1, 1], -1, 2)
        with self.assertRaises(ValueError):
            core.modpow([1, 1], -1, 7)
        with self.assertRaises(ArithmeticError):
            core.need(False, 'deliberate failed check')

    def test_invariant_and_CM_character_corruption_are_detected(self):
        p = 7
        g = core.central(core.sequence(*core.PARAMS['a'], p))
        exact = check_arithmetic.invpoly(g[:p], *core.PARAMS['a'])
        self.assertTrue(all(c % (p*p) == 0 for c in exact))
        self.assertEqual(exact[p], -2*p*p*g[p])
        bad = g[:p]
        bad[1] += 1
        self.assertTrue(any(c % (p*p) for c in check_arithmetic.invpoly(bad, *core.PARAMS['a'])))
        row = check_arithmetic.cm(p, 3, 1, 2)
        values = (row['unit_root_mod_p2'], -row['unit_root_mod_p2'] % 49,
                  2*row['unit_root_mod_p2'] % 49)
        self.assertEqual(len(set(values)), 3)
        with self.assertRaises(ArithmeticError):
            check_arithmetic.cm(p, 1, 1, 2)

    def test_exact_Hecke_factor_rejects_change_invisible_mod_p2(self):
        row = check_modular.hecke(7, 'a')
        coefficients = list(map(Fraction, row['specialized']))
        trace = row['raw_trace']
        self.assertEqual(check_modular.quadratic_remainder(coefficients, trace, 49), [0, 0])
        self.assertNotEqual(check_modular.quadratic_remainder(coefficients, -trace, 49), [0, 0])
        bad = coefficients[:]
        bad[0] += 49
        self.assertEqual([core.residue(x, 49) for x in bad],
                         [core.residue(x, 49) for x in coefficients])
        self.assertNotEqual(check_modular.quadratic_remainder(bad, trace, 49), [0, 0])
        with self.assertRaises(ArithmeticError):
            check_modular.exact_div([2, 4, 7], 2)

    def test_isogeny_curve_identity_rejects_scaling_sign_and_numerator(self):
        e = exact_isogeny
        a, b, k, n, f = e.data()
        A, B = e.neg(n), e.scale(e.power(k, 2), e.r(5))
        def residual(A, B, sign):
            D = e.add(e.mul(e.der(A), B), e.neg(e.mul(A, e.der(B))))
            rhs = e.add(e.mul(e.power(A, 3), B), e.scale(e.mul(A, e.power(B, 3)), a),
                        e.scale(e.power(B, 4), b))
            return e.add(e.mul(f, e.power(D, 2)), e.scale(rhs, e.r(sign)))
        self.assertEqual(residual(A, B, 5), [e.Z])
        for bad_A, bad_B, sign in ((A, e.power(k, 2), 5), (A, B, -5),
                                    (e.add(A, [e.O]), B, 5)):
            with self.subTest(sign=sign):
                self.assertNotEqual(residual(bad_A, bad_B, sign), [e.Z])
        with self.assertRaises(ZeroDivisionError):
            e.ri(e.Z)
        with self.assertRaises(ZeroDivisionError):
            e.divmodp([e.O], [e.Z])
        with self.assertRaises(ValueError):
            e.power([e.O], -1)

    def test_Hasse_square_and_integrality_domain_controls(self):
        f = finite_audit
        for p in (7, 11, 13, 19):
            H = f.hasse(p)
            self.assertEqual(H, f.hasse_translated(p))
            self.assertEqual(H, f.scale(f.legendre(p//4, p), f.leg(6, p), p))
        for n in (0, 1, 7, 20):
            self.assertEqual(f.mul(f.legendre(n), f.legendre(n)), f.square_formula(n))
        t = Fraction(1, 7)
        gate = 1-4*7*t
        difference = comb(12, 6)*t**6-gate**3*comb(12, 6)*(t/gate)**6
        self.assertEqual(difference, Fraction(176, 21609))
        denominator = difference.denominator
        valuation = 0
        while denominator % 7 == 0:
            denominator //= 7
            valuation -= 1
        self.assertEqual(valuation, -4)
        with self.assertRaises(ValueError):
            f.Field(7).inv((0, 0))


class EvidenceAndOutputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='sun37-tests-')
        self.base = Path(self.temp.name)
        self.actual = self.base/'actual'
        self.expected = self.base/'expected'
        self.actual.mkdir()
        self.expected.mkdir()
        for name in ALL_FILES:
            for directory in (self.actual, self.expected):
                (directory/name).write_bytes(b'{"value":1}\n')

    def tearDown(self):
        self.temp.cleanup()

    def cli(self, *arguments, no_site=False):
        flags = ['-B']
        if sys.flags.optimize:
            flags.append('-O')
        if no_site:
            flags.append('-S')
        return subprocess.run([sys.executable, *flags, str(ROOT/'scripts/verify.py'), *map(str, arguments)],
                              cwd=self.base, capture_output=True, text=True)

    def snapshot(self):
        return {p.relative_to(self.expected).as_posix(): p.read_bytes() if p.is_file() else None
                for p in self.expected.rglob('*')}

    def test_UTF8_LF_serialization_is_platform_independent(self):
        path = self.base/'writer.json'
        core.save(path, {'text': 'caffè', 'value': Fraction(1, 2).numerator})
        raw = path.read_bytes()
        self.assertNotIn(b'\r', raw)
        self.assertTrue(raw.endswith(b'\n'))
        self.assertEqual(json.loads(raw), {'text': 'caffè', 'value': 1})

    def test_every_byte_is_checked_and_frozen_input_is_preserved(self):
        self.assertEqual(compare_expected(self.actual, self.expected), 19)
        before = self.snapshot()
        for data in (b'{ "value":1}\n', b'{"value":true}\n', b'{"value":1}\r\n', b'{}', b''):
            (self.actual/ALL_FILES[0]).write_bytes(data)
            with self.subTest(data=data), self.assertRaises(ArithmeticError):
                compare_expected(self.actual, self.expected)
        self.assertEqual(before, self.snapshot())

    def test_missing_extra_directory_and_noncomplete_core_inventory_fail(self):
        (self.actual/ALL_FILES[0]).unlink()
        with self.assertRaises(ArithmeticError):
            compare_expected(self.actual, self.expected)
        (self.actual/ALL_FILES[0]).write_bytes(b'{"value":1}\n')
        (self.actual/'extra.json').write_bytes(b'{}')
        with self.assertRaises(ArithmeticError):
            compare_expected(self.actual, self.expected)
        (self.actual/'extra.json').unlink()
        for name in SYMBOLIC_FILES:
            (self.actual/name).unlink()
        self.assertEqual(compare_expected(self.actual, self.expected, core_only=True), 17)
        (self.expected/SYMBOLIC_FILES[0]).unlink()
        with self.assertRaises(ArithmeticError):
            compare_expected(self.actual, self.expected, core_only=True)
        (self.expected/SYMBOLIC_FILES[0]).mkdir()
        with self.assertRaises(ArithmeticError):
            compare_expected(self.actual, self.expected, core_only=True)

    def test_output_equal_descendant_and_ancestor_are_rejected_before_mkdir(self):
        child = self.expected/'new-child'
        self.assertFalse(child.exists())
        before = self.snapshot()
        paths = (ROOT, ROOT/'guard-regression-new-child', ROOT.parent,
                 self.expected, child, self.base)
        for path in paths:
            result = self.cli('--core-only', '--out-dir', path, '--expected-dir', self.expected)
            with self.subTest(path=path):
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('SUN37_REPLAY_FAIL: output must be disjoint', result.stderr)
        self.assertEqual(before, self.snapshot())
        self.assertFalse(child.exists())
        self.assertFalse((ROOT/'guard-regression-new-child').exists())
        with self.assertRaises(ValueError):
            compare_expected(self.expected, self.expected)

    def test_nonempty_output_is_not_modified_and_bad_inventory_creates_nothing(self):
        before = {p.name: p.read_bytes() for p in self.actual.iterdir()}
        result = self.cli('--core-only', '--out-dir', self.actual, '--expected-dir', self.expected)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('new or empty', result.stderr)
        self.assertEqual(before, {p.name: p.read_bytes() for p in self.actual.iterdir()})
        destination = self.base/'must-not-be-created'
        (self.expected/'extra.json').write_bytes(b'{}')
        result = self.cli('--core-only', '--out-dir', destination, '--expected-dir', self.expected)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('inventory mismatch', result.stderr)
        self.assertFalse(destination.exists())

    def test_missing_explicit_symbolic_dependency_fails_before_output(self):
        destination = self.base/'missing-symbolic-output'
        before = self.snapshot()
        result = self.cli('--out-dir', destination, '--expected-dir', self.expected, no_site=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Full verification requires SymPy 1.14.0', result.stderr)
        self.assertFalse(destination.exists())
        self.assertEqual(before, self.snapshot())

    def test_wrong_symbolic_version_fails_and_path_aliases_resolve_before_guard(self):
        with patch.dict(sys.modules, {'sympy': SimpleNamespace(__version__='1.13.0')}):
            with self.assertRaisesRegex(RuntimeError, 'requires SymPy 1.14.0; found 1.13.0'):
                verify.require_symbolic()
        alias = ROOT/'companion'/'..'/'evidence'/'nonexistent'
        with self.assertRaisesRegex(ValueError, 'disjoint from repository'):
            verify.validate_output_dir(alias, self.expected)
        before = self.snapshot()
        alias = self.expected/'..'/'expected'/'nonexistent'
        with self.assertRaisesRegex(ValueError, 'disjoint from expected evidence'):
            verify.validate_output_dir(alias, self.expected)
        self.assertEqual(before, self.snapshot())

    def test_optional_temporary_output_is_cleaned_and_names_omissions(self):
        before = {p.name for p in self.base.iterdir()}
        observed = []
        def fake_replay(directory, expected, core_only):
            self.assertTrue(directory.is_dir())
            self.assertEqual(expected, self.expected)
            self.assertTrue(core_only)
            observed.append(directory)
            return {'mode': 'core-only', 'omitted_files': list(SYMBOLIC_FILES)}
        output = io.StringIO()
        with patch.object(sys, 'argv', ['verify.py', '--core-only', '--expected-dir', str(self.expected)]), \
                patch.object(verify.tempfile, 'gettempdir', return_value=str(self.base)), \
                patch.object(verify, 'replay', side_effect=fake_replay), redirect_stdout(output):
            self.assertEqual(verify.main(), 0)
        self.assertEqual({p.name for p in self.base.iterdir()}, before)
        self.assertTrue(observed)
        self.assertFalse(observed[0].exists())
        self.assertEqual(json.loads(output.getvalue())['omitted_files'], list(SYMBOLIC_FILES))


if __name__ == '__main__':
    unittest.main()
