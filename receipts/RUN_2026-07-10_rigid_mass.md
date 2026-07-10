# RUN RECEIPT: rigid_mass_probe.py — 2026-07-10

## Purpose
Two coherence probes (explicitly NOT adjudication probes) for the
microscopic mass-generation route, per external review of same date.
Also includes the stabilization changes: pytest defaults to the quick
suite (addopts = -m "not slow"); the slow family test clears the sympy
global cache on entry to mitigate order-sensitivity; full-suite
single-process runs are documented as unsupported on constrained CI.

## Environment
Python 3.12.3, sympy 1.14.0

## Verbatim output
```
[PROBE 1] Rigid-mass axiom requirement analysis:
    R_MK = 6*beta*gamma/r**2 - 6*gamma/r + 12*kappa
    constant part (absorbable by lambda): 12*kappa
    r-dependent part (NOT absorbable):    6*gamma*(beta - r)/r**2
    constant S solvable for some constant lambda: False
    required anchoring term J(r) = S0*gamma*(-beta + r)/r**2
    J(r) vanishes when gamma -> 0: True
    RECEIPT: Under the macroscopic scalar EOM, constant mass in
    the MK frame is not selected. Treating MK-frame mass rigidity
    as physical therefore requires an additional microscopic/
    anomalous anchoring principle.
[PROBE 2] Curvature-insensitivity estimate (receipted gamma_0, kappa; arXiv:1007.0970):
    r =  10 kpc: |R_MK| = 4.81e-52 cm^-2;  R/m^2(proton) = 2.1e-79;  R/m^2(electron) = 7.2e-73
    r = 100 kpc: |R_MK| = 5.50e-53 cm^-2;  R/m^2(proton) = 2.4e-80;  R/m^2(electron) = 8.2e-74
    curvature-induced mass variation negligible (<1e-70): True
    CRUCIAL LIMITATION: curvature insensitivity shows mass
    RIGIDITY. It does not select the conformal FRAME in which
    the rigid mass is defined.

COHERENCE VERDICT: the microscopic route is coherent in principle
(rigid masses cannot feel galactic curvature) and demanding in
structure (it must supply a frame-anchoring term whose profile is
exactly the gamma-driven curvature variation of the MK metric).
Neither probe adjudicates whether such an anchoring mechanism
exists in the quantized theory; that is the conformal-anomaly /
dimensional-transmutation frontier, outside this lane's scope.

RECEIPT SUMMARY: {'anchor_required_and_gamma_driven': True, 'curvature_insensitivity_verified': True}
```

## Findings ledger
1. R_MK = 6*beta*gamma/r^2 - 6*gamma/r + 12*kappa: the constant part
   (12*kappa) is absorbable by the quartic coupling; the ENTIRE
   r-dependent part is proportional to gamma.
2. Constant S in the MK frame is unsolvable for any constant lambda;
   the required anchoring term has the closed form
   J(r) = S0*gamma*(r - beta)/r^2, vanishing as gamma -> 0. The
   anchoring deficit is exactly the curvature variation the linear term
   itself creates.
3. Curvature insensitivity verified with receipted parameters:
   |R_MK| ~ 5e-52 cm^-2 at 10 kpc; R/m^2 ~ 2e-79 (proton),
   ~ 7e-73 (electron). Rigidity holds; frame selection does not follow.

## Controlling receipt sentences
"Under the macroscopic scalar EOM, constant mass in the MK frame is not
selected. Treating MK-frame mass rigidity as physical therefore requires
an additional microscopic/anomalous anchoring principle."
"Curvature insensitivity shows mass rigidity. It does not select the
conformal frame in which the rigid mass is defined."

Test suite: quick 18 passed (1 deselected); slow lane 1 passed isolated.
