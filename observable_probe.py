"""
OBSERVABLE PROBE: the conformally invariant rotation observable
after mass generation.

Design: the earlier frame_probe showed the bare diagnostic v^2 = rB'/(2B)
is frame-dependent (diagnostic frame-dependence). This probe adds the
matter sector at the macroscopic classical level and constructs the
invariant that decides which frame's diagnostic is physical -- under an
explicitly stated, DISPUTED premise.

PART A -- the invariant exists.
With particle mass m(x) = h*S(x) and S of conformal weight -1, the
worldline action -h * Integral(S ds) is exactly conformally invariant:
(S/Omega)*(Omega ds) = S ds. Physical massive trajectories are therefore
frame-independent curves: geodesics of the mass-frame metric
ghat = (S/S0)^2 g, which is itself conformally invariant as a composite.
The physical rotation observable is the static-observer diagnostic of
ghat brought to standard static form.

PART B -- the observable under the two contested scalar profiles.
  (i)  S = S0 constant in the MK frame  =>  ghat = g_MK, observable is
       the MK diagnostic, gamma*r/2 support present.   [Mannheim-side]
  (ii) S = S0/(1 - a*r) with a the physical branch of frame_probe
       =>  ghat = Omega^2 g_MK = SdS-type metric, observable is the SdS
       diagnostic, gamma*r/2 support absent.           [Horne/HL-side]

PART C -- what the classical scalar field equation admits and rejects.
The conformally coupled scalar obeys
    Box S - (R/6) S - 4*lambda*S^3 = 0.
The probe computes, on the fixed MK background:
  (1) whether S = S0/(1 - a*r) solves the scalar equation exactly on that
      background (i.e. whether E := Box S - (R/6) S is proportional to
      S^3 with CONSTANT coefficient, which then fixes lambda);
  (2) whether S = S0 constant can solve the equation (it requires
      R_MK to be constant, which the probe tests directly).

VERDICT SHAPE (conditional, per lane discipline):
IF particle masses track a macroscopic classical scalar solving the
conformal EOM in this fixed background, THEN the conformally invariant
rotation observable associated with the admitted SdS-generating profile is
that SdS diagnostic and the linear support cancels -- the Horne
(arXiv:1601.07537) / Hobson-Lasenby (PRD 104, 064014) line, rendered as
computation under that premise.
THIS DOES NOT YET ESTABLISH a coupled gravity-matter solution, uniqueness
of the scalar profile, or frame selection in the stronger physical sense.
It also does not establish that the required quartic coupling is source-
independent; in the current fixed-background construction the derived
lambda depends explicitly on beta.
THE PREMISE IS THE DISPUTE: Mannheim (GRG 54, 99, 2022) holds that the
physical S is built from microscopic vacuum expectation values varying
only within particle interiors, not a macroscopic classical field. That
counter-position is not computable at this level and is NOT adjudicated
here.

Usage: python3 observable_probe.py
"""
import sympy as sp
from bach_receipt import DiagMetric, t, r, th, ph, beta, gam, k

rho = sp.symbols('rho', positive=True)
S0, lam, Om = sp.symbols('S0 lambda Omega', positive=True)

B_MK = 1 - 3*beta*gam - beta*(2 - 3*beta*gam)/r + gam*r - k*r**2
A_PHYS = gam / (3*beta*gam - 2)   # physical branch from frame_probe


def part_a():
    S, ds = sp.symbols('S ds', positive=True)
    action_density = S * ds
    transformed = (S / Om) * (Om * ds)
    inv_action = sp.simplify(transformed - action_density) == 0
    g_sym = sp.symbols('g', positive=True)
    ghat = (S / S0)**2 * g_sym
    ghat_transformed = ((S / Om) / S0)**2 * (Om**2 * g_sym)
    inv_ghat = sp.simplify(ghat_transformed - ghat) == 0
    print("[PART A] Invariant structure:")
    print(f"    worldline action S*ds conformally invariant: {inv_action}")
    print(f"    mass-frame metric (S/S0)^2 g conformally invariant: {inv_ghat}")
    return inv_action and inv_ghat


def part_b():
    def vsq(Bfun, var):
        return sp.cancel(var * sp.diff(Bfun, var) / (2 * Bfun))

    # (i) constant S: ghat = g_MK
    v2_i = vsq(B_MK, r)
    lin_i = sp.simplify(sp.limit(sp.diff(v2_i, r), r, sp.oo))  # asymptotic slope
    # weak-field linear coefficient: expand v2 for small params
    eps = sp.symbols('epsilon', positive=True)
    v2_eps = sp.cancel(v2_i.subs({beta: eps*beta, gam: eps*gam, k: eps*k}))
    lead_i = sp.expand(sp.series(v2_eps, eps, 0, 2).removeO().coeff(eps, 1))
    has_linear_i = sp.simplify(lead_i.coeff(r, 1) - gam/2) == 0

    # (ii) S = S0/(1 - a r): ghat = Omega^2 g_MK -> SdS-type standard form
    a = sp.symbols('a', real=True)
    Omega = 1 / (1 - a * r)
    Bhat = sp.cancel((Omega**2 * B_MK).subs(r, rho / (1 + a * rho)))
    p = sp.Poly(sp.expand(sp.cancel(Bhat) * rho), rho)
    coeffs = {mono - 1: sp.factor(sp.simplify(c.subs(a, A_PHYS))) for mono, c in
              zip(reversed(range(p.degree() + 1)), p.all_coeffs())}
    coeffs = {pw: c for pw, c in coeffs.items() if c != 0}
    no_linear_ii = 1 not in coeffs
    print("[PART B] Observable under contested scalar profiles:")
    print(f"    (i)  S const in MK frame: weak-field linear rotation "
          f"coefficient = gamma/2 present: {has_linear_i}")
    print(f"    (ii) S = S0/(1-a r): mass-frame metric is SdS-type, "
          f"linear term absent: {no_linear_ii}")
    return has_linear_i and no_linear_ii


def part_c():
    M = DiagMetric([-B_MK, 1/B_MK, r**2, r**2*sp.sin(th)**2], [t, r, th, ph])
    R_MK = M.Rs

    def box(Sfun):
        # Box S = (1/sqrt(-g)) d_r ( sqrt(-g) g^rr d_r S ), radial S
        return sp.cancel(sp.diff(r**2 * B_MK * sp.diff(Sfun, r), r) / r**2)

    # Candidate (ii): S = S0/(1 - a r), a = physical branch
    S_ii = S0 / (1 - A_PHYS * r)
    E_ii = sp.cancel(sp.together(box(S_ii) - R_MK * S_ii / 6))
    ratio = sp.simplify(sp.cancel(E_ii / S_ii**3))
    is_exact = sp.simplify(sp.diff(ratio, r)) == 0
    lam_val = sp.simplify(ratio / 4)

    # Candidate (i): S = S0 constant
    E_i = sp.simplify(box(S0) - R_MK * S0 / 6)
    R_is_const = sp.simplify(sp.diff(R_MK, r)) == 0
    const_solves = R_is_const  # constant S solves EOM iff R constant

    print("[PART C] Classical conformal scalar EOM on the fixed MK background:")
    print(f"    R_MK constant (Einstein-space condition): {R_is_const}")
    print(f"    dR_MK/dr = {sp.simplify(sp.diff(R_MK, r))}")
    print(f"    S = S0 constant solves the EOM: {const_solves}")
    print(f"    S = S0/(1 - a_phys r): (Box S - R S/6)/S^3 = {ratio}")
    print(f"    ... constant in r (exact scalar-EOM solution on this background): {is_exact}")
    print(f"    ... required quartic coupling 4*lambda = {sp.simplify(ratio)}")
    print(f"    ... beta-dependence of lambda vanishes: {sp.simplify(sp.diff(lam_val, beta)) == 0}")
    return is_exact, const_solves, lam_val


if __name__ == '__main__':
    a_ok = part_a()
    b_ok = part_b()
    exact_ii, const_ok, lam_val = part_c()
    print()
    print("VERDICT (conditional):")
    print("  Under the macroscopic classical scalar treatment, the scalar EOM")
    print("  on the fixed MK background admits S = S0/(1 - a_phys r) and does")
    print("  NOT admit constant S (since R_MK is not constant). Under that")
    print("  premise, the conformally invariant rotation observable associated")
    print("  with the admitted profile is the SdS diagnostic, so the")
    print("  gamma*r/2 support cancels in that observable.")
    print("  This does NOT yet establish a coupled gravity-matter solution,")
    print("  uniqueness of the scalar profile, or stronger frame selection.")
    print("  It also does NOT show that the derived quartic coupling is")
    print("  source-independent: in this fixed-background construction,")
    print("  lambda depends explicitly on beta.")
    print("  DISPUTED PREMISE (not adjudicated): Mannheim (GRG 54, 99, 2022)")
    print("  denies that the mass-generating S is a macroscopic classical")
    print("  field; microscopic vevs are outside this probe's scope.")
    print()
    print(f"RECEIPT SUMMARY: {{'invariant_structure': {a_ok}, "
          f"'observable_profiles': {b_ok}, "
          f"'scalar_EOM_admits_SdS_profile_on_fixed_background': {exact_ii}, "
          f"'constant_S_fails_in_MK_frame': {not const_ok}}}")
    raise SystemExit(0 if a_ok and b_ok and exact_ii and (not const_ok) else 1)
