"""Pins for the observable probe, including the covariance cross-check."""
import sys, pathlib
import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from observable_probe import (part_a, part_b, part_c, B_MK, A_PHYS,
                              beta, gam, k, r, S0)


def test_invariant_structure():
    assert part_a()


def test_observable_under_both_profiles():
    assert part_b()


def test_scalar_eom_selects_sds_frame_and_rejects_constant_S():
    exact_ii, const_ok, lam_val = part_c()
    assert exact_ii          # S = S0/(1 - a_phys r) is an exact solution
    assert not const_ok      # constant S fails: R_MK not constant


def test_einstein_failure_is_proportional_to_gamma():
    from bach_receipt import DiagMetric, t, th, ph
    M = DiagMetric([-B_MK, 1/B_MK, r**2, r**2*sp.sin(th)**2], [t, r, th, ph])
    dR = sp.simplify(sp.diff(M.Rs, r))
    assert sp.simplify(dR - 6*gam*(r - 2*beta)/r**3) == 0
    assert sp.simplify(dR.subs(gam, 0)) == 0


def test_quartic_coupling_matches_sds_frame_constant_scalar_relation():
    # 4*lambda must equal 2*s_hat/S0^2 where s_hat is the SdS quadratic
    # coefficient from the frame probe -- conformal covariance cross-check.
    _, _, lam_val = part_c()
    s_hat = -(9*beta**2*gam**2*k - beta*gam**3 - 12*beta*gam*k
              + gam**2 + 4*k) / (3*beta*gam - 2)**2
    assert sp.simplify(4*lam_val - 2*s_hat/S0**2) == 0
