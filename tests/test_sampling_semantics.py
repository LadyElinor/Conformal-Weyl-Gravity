from conformal_probe.core.geometry import component_rng, sample_points
from bach_receipt import r, th, beta, gam, k, Lam


def test_component_rng_is_repeatable_for_same_case_and_component():
    rng1 = component_rng("case-a", (2, 3))
    rng2 = component_rng("case-a", (2, 3))
    seq1 = [rng1.randint(0, 10**6) for _ in range(5)]
    seq2 = [rng2.randint(0, 10**6) for _ in range(5)]
    assert seq1 == seq2


def test_component_rng_is_component_specific():
    rng1 = component_rng("case-a", (2, 3))
    rng2 = component_rng("case-a", (2, 2))
    seq1 = [rng1.randint(0, 10**6) for _ in range(5)]
    seq2 = [rng2.randint(0, 10**6) for _ in range(5)]
    assert seq1 != seq2


def test_sample_points_are_call_order_independent_and_component_keyed():
    symbols = {r, th, beta, gam, k, Lam}
    first = sample_points(symbols, r, th, beta, gam, k, Lam, case_id="mk", component=(0, 0))
    second = sample_points(symbols, r, th, beta, gam, k, Lam, case_id="mk", component=(0, 0))
    other = sample_points(symbols, r, th, beta, gam, k, Lam, case_id="mk", component=(2, 2))
    assert first == second
    assert first != other
    assert len(first) == 3
    assert all(r in point and th in point for point in first)
