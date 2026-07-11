"""
RIGID MASS PROBE: two COHERENCE probes (not adjudication probes) for the
microscopic mass-generation route.

Context: observable_probe established the conditional -- under the
macroscopic classical scalar realization of mass generation, the invariant
rotation observable is the SdS diagnostic and the gamma*r/2 support
cancels. Mannheim's remaining route (GRG 54, 99, 2022) is microscopic /
anomalous mass generation keeping particle masses rigid in the MK frame.
These probes render what that route REQUIRES and what part of it is
already verifiable, without pretending to derive it.

PROBE 1 -- rigid-mass axiom requirement analysis.
Axiom: m(x) = m0 constant in the MK frame at galactic scales.
Consequences rendered as computation:
  (a) the mass-frame metric is then the MK metric and the invariant
      observable retains the gamma*r/2 support;
  (b) but constant S cannot solve the conformal scalar EOM in the MK
      background for ANY constant quartic coupling: the residual is
      J(r) = -(S0/6) * (R_MK(r) - R_bar), whose r-dependent part is
      proportional to gamma -- the anchoring deficit is exactly the
      curvature variation the linear term itself creates;
  (c) therefore the axiom requires an additional frame-anchoring
      principle OUTSIDE the macroscopic scalar sector.
Receipt sentence: "Under the macroscopic scalar EOM, constant mass in
the MK frame is not selected. Treating MK-frame mass rigidity as
physical therefore requires an additional microscopic/anomalous
anchoring principle."

PROBE 2 -- curvature-insensitivity estimate.
Verifies the strongest fair point on the Mannheim side: microscopic
masses cannot respond appreciably to galactic curvature. Gap-equation-
style curvature corrections scale as R/m^2 (m as inverse reduced Compton
wavelength in geometrized units). With the receipted universal parameters
(gamma_0 = 3.06e-30 cm^-1, kappa = 9.54e-54 cm^-2, arXiv:1007.0970) the
ratio at galactic radii is ~1e-80 for the proton and ~1e-73 for the
electron.
CRUCIAL LIMITATION (the whole point): curvature insensitivity shows mass
RIGIDITY. It does not select the conformal FRAME in which the rigid mass
is defined. Rigidity is necessary for Mannheim's route; it is not the
anchoring principle itself.

Usage: python3 rigid_mass_probe.py
"""
import sympy as sp
from bach_receipt import DiagMetric, t, r, th, ph, beta, gam, k

S0, lam = sp.symbols('S0 lambda', positive=True)
B_MK = 1 - 3*beta*gam - beta*(2 - 3*beta*gam)/r + gam*r - k*r**2


def probe1():
    M = DiagMetric([-B_MK, 1/B_MK, r**2, r**2*sp.sin(th)**2], [t, r, th, ph])
    R_MK = sp.simplify(M.Rs)
    # (b) residual of the EOM at S = S0 for arbitrary constant lambda:
    #     E(r) = -(R_MK/6) S0 - 4 lambda S0^3.
    # A constant lambda can absorb only the constant part of R_MK.
    E = -R_MK * S0 / 6 - 4 * lam * S0**3
    # decompose R_MK into constant part and r-dependent part
    R_const = sp.limit(R_MK, r, sp.oo)
    R_var = sp.simplify(R_MK - R_const)
    solvable = sp.simplify(R_var) == 0
    anchor_deficit = sp.factor(-S0 * R_var / 6)
    deficit_gamma_free = sp.simplify(anchor_deficit.subs(gam, 0))
    print("[PROBE 1] Rigid-mass axiom requirement analysis:")
    print(f"    R_MK = {R_MK}")
    print(f"    constant part (absorbable by lambda): {sp.simplify(R_const)}")
    print(f"    r-dependent part (NOT absorbable):    {R_var}")
    print(f"    constant S solvable for some constant lambda: {solvable}")
    print(f"    required anchoring term J(r) = {anchor_deficit}")
    print(f"    J(r) vanishes when gamma -> 0: {deficit_gamma_free == 0}")
    print("    RECEIPT: Under the macroscopic scalar EOM, constant mass in")
    print("    the MK frame is not selected. Treating MK-frame mass rigidity")
    print("    as physical therefore requires an additional microscopic/")
    print("    anomalous anchoring principle.")
    return (not solvable) and (deficit_gamma_free == 0), anchor_deficit


def probe2():
    M = DiagMetric([-B_MK, 1/B_MK, r**2, r**2*sp.sin(th)**2], [t, r, th, ph])
    R_MK = sp.simplify(M.Rs)
    # receipted universal parameters + fiducial 1e11 M_sun galaxy
    vals = {beta: sp.Float(1.4766e16), gam: sp.Float(3.06e-30),
            k: sp.Float(9.54e-54)}
    # particle masses as inverse reduced Compton wavelengths (cm^-1)
    m2 = {'proton':   (1 / sp.Float(2.1031e-14))**2,
          'electron': (1 / sp.Float(3.8616e-11))**2}
    print("[PROBE 2] Curvature-insensitivity estimate "
          "(receipted gamma_0, kappa; arXiv:1007.0970):")
    ok = True
    worst = 0.0
    for kpc in (10, 100):
        rcm = sp.Float(kpc * 3.0857e21)
        Rnum = abs(float(R_MK.subs(vals).subs(r, rcm)))
        line = f"    r = {kpc:>3d} kpc: |R_MK| = {Rnum:.2e} cm^-2"
        for name, msq in m2.items():
            ratio = Rnum / float(msq)
            worst = max(worst, ratio)
            line += f";  R/m^2({name}) = {ratio:.1e}"
            ok = ok and ratio < 1e-70
        print(line)
    print(f"    curvature-induced mass variation negligible (<1e-70): {ok}")
    print("    CRUCIAL LIMITATION: curvature insensitivity shows mass")
    print("    RIGIDITY. It does not select the conformal FRAME in which")
    print("    the rigid mass is defined.")
    return ok, worst


if __name__ == '__main__':
    ok1, deficit = probe1()
    ok2, worst = probe2()
    print()
    print("COHERENCE VERDICT: the microscopic route is coherent in principle")
    print("(rigid masses cannot feel galactic curvature) and demanding in")
    print("structure (it must supply a frame-anchoring term whose profile is")
    print("exactly the gamma-driven curvature variation of the MK metric).")
    print("Neither probe adjudicates whether such an anchoring mechanism")
    print("exists in the quantized theory; that is the conformal-anomaly /")
    print("dimensional-transmutation frontier, outside this lane's scope.")
    print()
    print(f"RECEIPT SUMMARY: {{'anchor_required_and_gamma_driven': {ok1}, "
          f"'curvature_insensitivity_verified': {ok2}}}")
    raise SystemExit(0 if ok1 and ok2 else 1)
