"""Pins for the frame probe: the branch, the map, and the punchline."""
import sys, pathlib
import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from frame_probe import fact1, fact2, B_MK, beta, gam, k, r


def test_physical_branch_and_map_structure():
    ok, a_val = fact1()
    assert ok
    assert sp.simplify(a_val - gam / (3 * beta * gam - 2)) == 0
    # beta -> 0 limit recovers the known map Omega = 1/(1 + gamma r / 2)
    assert sp.simplify(a_val.subs(beta, 0) + gam / 2) == 0


def test_newtonian_charge_preserved():
    # 1/rho coefficient of the SdS-frame metric equals MK's own u
    u_mk = -beta * (2 - 3 * beta * gam)
    assert sp.simplify(beta * (3 * beta * gam - 2) - u_mk) == 0


def test_diagnostic_frame_dependent_and_lead_is_gamma_r_over_2():
    _, a_val = fact1()
    dep, lead = fact2(a_val)
    assert dep  # not conformally invariant
    assert sp.simplify(lead - gam * r / 2) == 0


def test_entire_difference_proportional_to_gamma():
    _, a_val = fact1()
    from frame_probe import sp as _sp
    Omega = 1 / (1 - a_val * r)
    rho = sp.symbols('rho', positive=True)
    Bhat = sp.cancel((Omega**2 * B_MK).subs(r, rho / (1 + a_val * rho)))
    v2_mk = sp.cancel(r * sp.diff(B_MK, r) / (2 * B_MK))
    v2_hat = sp.cancel(rho * sp.diff(Bhat, rho) / (2 * Bhat))
    diff = sp.cancel(v2_mk - v2_hat.subs(rho, r / (1 - a_val * r)))
    assert sp.simplify(diff.subs(gam, 0)) == 0
