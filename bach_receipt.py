"""
RECEIPT v2 (hardened): vacuum Bach equation tests for conformal gravity metrics.

DEFINITION (cited, not from memory):
    B_ab = 2 * D^m D^n C_{mabn} + C_{mabn} R^{mn}
Sources: arXiv:gr-qc/9811086 eq.(2); arXiv:1105.5632 App. C (Wald conventions).
Properties: symmetric, trace-free, conformally invariant (weight -1).

Conventions here: signature (-,+,+,+); Riemann R^a_{bcd} = d_c Gam^a_{bd}
- d_d Gam^a_{bc} + Gam Gam; Ricci R_bd = R^a_{bad}. Weyl decomposition is the
standard 4D one for these conventions.

HARDENING over v1:
  Fix 1: explicit Schwarzschild-de Sitter control case (must PASS).
  Fix 2: definition taken verbatim from cited sources; plus a structural
         identity test -- the trace g^{ab} B_ab must vanish identically even
         for a metric whose Bach tensor is NONZERO (corrupted case). This
         tests the implementation, not the metric.
  Fix 3: all 10 independent components (a <= b) computed, not just diagonal.
  Fix 4: exact-point numerical cross-check. Each raw (pre-simplification)
         component is evaluated at random exact rational parameter points;
         PASS requires symbolic zero AND exact zero at every sample point.
         This removes trust in simplify() as a sole oracle.
  Extra: Weyl tensor implementation check -- conformally flat FRW must give
         identically vanishing Weyl for generic a(t).

Usage: python3 bach_receipt_v2.py [mk|sds|trunc|bad|frw|all]
"""
import sys
import time
import sympy as sp

from conformal_probe.core.geometry import (
    N_DIMENSIONS,
    DiagMetric,
    sample_points,
    simplify_expr as S,
)

t, r, th, ph = sp.symbols('t r theta phi', real=True)
beta, gam, k, Lam = sp.symbols('beta gamma kappa Lambda', real=True)
n = N_DIMENSIONS


def run_case(name, gdiag, coords, expect_pass, check_trace=False):
    t0 = time.time()
    M = DiagMetric(gdiag, coords)
    memo = {}
    results = {}
    all_syms = set().union(*[e.free_symbols for e in gdiag]) | {r, th}
    numeric_ok = True
    for a in range(n):
        for b in range(a, n):
            raw, simp = M.bach(a, b, memo)
            results[(a, b)] = simp
            # Fix 4: exact-point evaluation of the RAW expression
            for pt in sample_points(all_syms, r, th, beta, gam, k, Lam):
                val = sp.cancel(raw.subs(pt))
                sym_zero = (simp == 0)
                num_zero = (val == 0)
                if sym_zero != num_zero:
                    numeric_ok = False
                    print(f"    DISCREPANCY B_{a}{b}: symbolic zero={sym_zero}, "
                          f"exact-point zero={num_zero} at {pt}")
    nonzero = {ab: v for ab, v in results.items() if v != 0}
    passed = (not nonzero) if expect_pass else bool(nonzero)
    verdict = "PASS" if passed else "UNEXPECTED"
    kind = "Bach = 0 (all 10 components)" if not nonzero else \
           f"Bach != 0 ({len(nonzero)} nonzero components)"
    print(f"[{name}] {verdict}: {kind}  "
          f"[exact-point check {'consistent' if numeric_ok else 'INCONSISTENT'}]"
          f"  ({time.time()-t0:.0f}s)")
    for (a, b), v in list(nonzero.items())[:4]:
        print(f"    B_{a}{b} = {v}")
    trace_ok = None
    if check_trace:
        tr = sp.simplify(S(sum(M.gu(a, a) * results[(a, a)] for a in range(n))))
        trace_ok = (tr == 0)
        print(f"    IDENTITY g^ab B_ab = {tr}  "
              f"({'PASS: traceless even with nonzero Bach' if trace_ok else 'FAIL'})")
    return passed and numeric_ok and (trace_ok in (None, True))


def frw_weyl_check():
    """Conformally flat FRW with generic a(t): Weyl must vanish identically."""
    a_t = sp.Function('a')(t)
    x, y, z = sp.symbols('x y z', real=True)
    M = DiagMetric([-sp.Integer(1), a_t**2, a_t**2, a_t**2], [t, x, y, z])
    ok = (len(M.Weyl) == 0)
    print(f"[FRW conformally-flat] {'PASS' if ok else 'FAIL'}: "
          f"Weyl components identically zero = {ok} (generic a(t))")
    return ok


CASES = {
    'mk': lambda: run_case(
        "Mannheim-Kazanas exact",
        [-(1 - 3*beta*gam - beta*(2 - 3*beta*gam)/r + gam*r - k*r**2),
         1/(1 - 3*beta*gam - beta*(2 - 3*beta*gam)/r + gam*r - k*r**2),
         r**2, r**2*sp.sin(th)**2],
        [t, r, th, ph], expect_pass=True),
    'sds': lambda: run_case(
        "Schwarzschild-de Sitter control (Fix 1)",
        [-(1 - 2*beta/r - Lam*r**2/3),
         1/(1 - 2*beta/r - Lam*r**2/3),
         r**2, r**2*sp.sin(th)**2],
        [t, r, th, ph], expect_pass=True),
    'trunc': lambda: run_case(
        "Truncated 1-2b/r+gr-kr^2 (expect O(bg) residual)",
        [-(1 - 2*beta/r + gam*r - k*r**2),
         1/(1 - 2*beta/r + gam*r - k*r**2),
         r**2, r**2*sp.sin(th)**2],
        [t, r, th, ph], expect_pass=False),
    'bad': lambda: run_case(
        "Corrupted gamma*r^3 + trace identity (Fix 2)",
        [-(1 - 2*beta/r + gam*r**3),
         1/(1 - 2*beta/r + gam*r**3),
         r**2, r**2*sp.sin(th)**2],
        [t, r, th, ph], expect_pass=False, check_trace=True),
    'frw': frw_weyl_check,
}

# ---------------------------------------------------------------------------
# Improvement 1: divergence-free identity D^a B_ab = 0
# Tested on the corrupted metric (nonzero Bach), so it probes the
# implementation and the identity, not the metric.
# ---------------------------------------------------------------------------
def divergence_identity_check():
    B = 1 - 2 * beta / r + gam * r**3
    M = DiagMetric([-B, 1 / B, r**2, r**2 * sp.sin(th)**2], [t, r, th, ph])
    memo = {}
    Bab = {}
    for a in range(n):
        for b in range(a, n):
            _, simp = M.bach(a, b, memo)
            Bab[(a, b)] = simp
            Bab[(b, a)] = simp
    ok = True
    for b in range(n):
        div = sum(M.gu(a, a) * (sp.diff(Bab[(a, b)], M.coords[a])
                                - sum(M.G(c, a, a) * Bab[(c, b)] for c in range(n))
                                - sum(M.G(c, a, b) * Bab[(a, c)] for c in range(n)))
                  for a in range(n))
        div = sp.simplify(S(div))
        if div != 0:
            ok = False
            print(f"    (div B)_{b} = {div}")
    print(f"[Divergence identity D^a B_ab = 0] "
          f"{'PASS: holds exactly on a metric with NONZERO Bach' if ok else 'FAIL'}")
    return ok


# ---------------------------------------------------------------------------
# Improvement 2: constraint discovery over the full polynomial family
#     B(r) = w + u/r + v*r + s*r^2
# Impose Bach flatness and let the algebra return the constraint among
# (w, u, v, s). Expected (MK 1989, rediscovered here, not assumed): the
# vacuum condition is exactly  3*u*v = w^2 - 1, with s unconstrained.
# ---------------------------------------------------------------------------
def constraint_discovery():
    w, u, v, s = sp.symbols('w u v s', real=True)
    B = w + u / r + v * r + s * r**2
    M = DiagMetric([-B, 1 / B, r**2, r**2 * sp.sin(th)**2], [t, r, th, ph])
    memo = {}
    conditions = set()
    for a in range(n):
        for b in range(a, n):
            _, simp = M.bach(a, b, memo)
            if simp == 0:
                continue
            num, _ = sp.fraction(sp.together(simp))
            num = sp.expand(num / sp.gcd_terms(num).as_coeff_Mul()[0]) \
                if num.is_Mul else sp.expand(num)
            poly = sp.Poly(sp.expand(sp.fraction(sp.together(simp))[0]),
                           r, sp.sin(th))
            for coeff in poly.coeffs():
                c = sp.factor(coeff)
                if c != 0:
                    conditions.add(c)
    # Reduce the condition set
    sols = sp.solve(list(conditions), [w, u, v, s], dict=True)
    print("[Constraint discovery over B = w + u/r + v*r + s*r^2]")
    print(f"    raw vanishing conditions: {sorted(map(str, conditions))}")
    print(f"    solve() branches: {sols}")
    # Verify the expected constraint branch: u = (w^2 - 1)/(3 v)
    Bc = B.subs(u, (w**2 - 1) / (3 * v))
    Mc = DiagMetric([-Bc, 1 / Bc, r**2, r**2 * sp.sin(th)**2], [t, r, th, ph])
    memo2 = {}
    residuals = [sp.simplify(Mc.bach(a, b, memo2)[1])
                 for a in range(n) for b in range(a, n)]
    ok = all(x == 0 for x in residuals)
    print(f"    VERIFY: substituting 3*u*v = w^2 - 1 kills all 10 components: {ok}")
    # Cross-check against MK parametrization: w = 1-3bg, u = -b(2-3bg), v = g
    lhs = sp.expand(3 * (-beta * (2 - 3 * beta * gam)) * gam)
    rhs = sp.expand((1 - 3 * beta * gam)**2 - 1)
    mk_ok = sp.simplify(lhs - rhs) == 0
    print(f"    MK parametrization satisfies 3uv = w^2-1 identically: {mk_ok}")
    # And the truncation w=1, u=-2b, v=g violates it by exactly -6bg
    viol = sp.simplify(3 * (-2 * beta) * gam - (1**2 - 1))
    print(f"    Truncation violation 3uv - (w^2-1) = {viol}  (source of O(bg) residual)")
    return ok and mk_ok


CASES['div'] = divergence_identity_check
CASES['family'] = constraint_discovery


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    names = list(CASES) if which == 'all' else [which]
    outcomes = {nm: CASES[nm]() for nm in names}
    print()
    print("RECEIPT SUMMARY:", {nm: bool(v) for nm, v in outcomes.items()})
