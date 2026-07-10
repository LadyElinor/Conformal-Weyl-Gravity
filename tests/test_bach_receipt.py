"""Pytest wrapper for the Bach receipt cases.

Pins not just pass/fail but the truncation residual VALUES, so a silent
convention or normalization drift fails loudly instead of renormalizing.
"""
import sys
import pathlib

import sympy as sp
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from bach_receipt import (  # noqa: E402
    DiagMetric, run_case, frw_weyl_check, CASES, t, r, th, ph, beta, gam, k,
)


def test_mk_exact_is_vacuum_solution():
    assert CASES['mk']()


def test_sds_control_passes():
    assert CASES['sds']()


def test_truncated_form_fails_at_order_beta_gamma():
    assert CASES['trunc']()  # run_case returns True when FAIL was expected
    # Pin the residual values against the cited-convention Bach tensor.
    B = 1 - 2 * beta / r + gam * r - k * r**2
    M = DiagMetric([-B, 1 / B, r**2, r**2 * sp.sin(th)**2], [t, r, th, ph])
    memo = {}
    _, b22 = M.bach(2, 2, memo)
    assert sp.simplify(b22 - 2 * beta * gam / r**2) == 0
    _, b33 = M.bach(3, 3, memo)
    assert sp.simplify(b33 - 2 * beta * gam * sp.sin(th)**2 / r**2) == 0


def test_corrupted_metric_fails_and_bach_is_traceless():
    assert CASES['bad']()


def test_frw_weyl_vanishes_for_generic_scale_factor():
    assert frw_weyl_check()


def test_bach_offdiagonals_vanish_for_mk():
    B = (1 - 3 * beta * gam - beta * (2 - 3 * beta * gam) / r
         + gam * r - k * r**2)
    M = DiagMetric([-B, 1 / B, r**2, r**2 * sp.sin(th)**2], [t, r, th, ph])
    memo = {}
    for a in range(4):
        for b in range(a + 1, 4):
            _, simp = M.bach(a, b, memo)
            assert simp == 0, f"off-diagonal B_{a}{b} nonzero"


def test_divergence_identity_on_nonzero_bach():
    from bach_receipt import divergence_identity_check
    assert divergence_identity_check()


@pytest.mark.slow
def test_family_constraint_is_3uv_eq_w2_minus_1():
    from sympy.core.cache import clear_cache
    clear_cache()  # order-sensitivity mitigation: drop state from prior tests
    from bach_receipt import constraint_discovery
    assert constraint_discovery()
    # Pin the constraint arithmetic directly, independent of the probe:
    w_, u_, v_ = sp.symbols('w u v', real=True)
    mk_w = 1 - 3 * beta * gam
    mk_u = -beta * (2 - 3 * beta * gam)
    mk_v = gam
    assert sp.simplify(3 * mk_u * mk_v - (mk_w**2 - 1)) == 0
