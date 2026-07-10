"""Pins for the coherence probes."""
import sys, pathlib
import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from rigid_mass_probe import probe1, probe2, beta, gam, r, S0


def test_anchor_required_and_gamma_driven():
    ok, deficit = probe1()
    assert ok
    # Pin the closed form of the anchoring deficit
    assert sp.simplify(deficit - S0 * gam * (r - beta) / r**2) == 0


def test_curvature_insensitivity_scale():
    ok, worst = probe2()
    assert ok
    assert worst < 1e-70          # rigidity threshold
    assert worst > 1e-85          # sanity: nonzero, right magnitude
