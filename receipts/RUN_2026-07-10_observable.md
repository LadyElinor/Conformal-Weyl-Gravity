# RUN RECEIPT: observable_probe.py — 2026-07-10

## Purpose
Priority-1 module from external review: construct the conformally
invariant rotation observable after mass generation, at the macroscopic
classical scalar level, and let the field equations select the frame.

## Environment
Python 3.12.3, sympy 1.14.0

## Verbatim output
```
[PART A] Invariant structure:
    worldline action S*ds conformally invariant: True
    mass-frame metric (S/S0)^2 g conformally invariant: True
[PART B] Observable under contested scalar profiles:
    (i)  S const in MK frame: weak-field linear rotation coefficient = gamma/2 present: True
    (ii) S = S0/(1-a r): mass-frame metric is SdS-type, linear term absent: True
[PART C] Classical conformal scalar EOM in the exact MK background:
    R_MK constant (Einstein-space condition): False
    dR_MK/dr = 6*gamma*(-2*beta + r)/r**3
    S = S0 constant solves the EOM: False
    S = S0/(1 - a_phys r): (Box S - R S/6)/S^3 = 2*(-9*beta**2*gamma**2*kappa + beta*gamma**3 + 12*beta*gamma*kappa - gamma**2 - 4*kappa)/(S0**2*(9*beta**2*gamma**2 - 12*beta*gamma + 4))
    ... constant in r (EXACT solution): True
    ... required quartic coupling 4*lambda = 2*(-9*beta**2*gamma**2*kappa + beta*gamma**3 + 12*beta*gamma*kappa - gamma**2 - 4*kappa)/(S0**2*(9*beta**2*gamma**2 - 12*beta*gamma + 4))

VERDICT (conditional):
  Under the macroscopic classical scalar treatment, the EOM in the
  MK background is solved exactly by S = S0/(1 - a_phys r) and is
  NOT solved by constant S (since R_MK is not constant). The
  conformally invariant rotation observable is therefore the SdS
  diagnostic: the gamma*r/2 support cancels. This renders the
  Horne / Hobson-Lasenby conclusion as computation.
  DISPUTED PREMISE (not adjudicated): Mannheim (GRG 54, 99, 2022)
  denies that the mass-generating S is a macroscopic classical
  field; microscopic vevs are outside this probe's scope.

RECEIPT SUMMARY: {'invariant_structure': True, 'observable_profiles': True, 'scalar_EOM_selects_SdS_frame': True, 'constant_S_fails_in_MK_frame': True}
```

## Findings ledger
1. The worldline action S*ds and the mass-frame metric (S/S0)^2 g are
   exactly conformally invariant: physical massive trajectories are
   frame-independent curves. The invariant observable exists.
2. Constant S in the MK frame yields the MK diagnostic (gamma*r/2
   present); S = S0/(1 - a_phys r) yields the SdS diagnostic (linear
   term absent). The two contested positions map to two scalar profiles.
3. In the exact MK background the conformal scalar EOM
   Box S - (R/6)S - 4*lambda*S^3 = 0 is solved EXACTLY by
   S = S0/(1 - a_phys r), with constant 4*lambda = 2*s_hat/S0^2.
   Constant S is NOT a solution: dR_MK/dr = 6*gamma*(r - 2*beta)/r^3,
   nonzero precisely in proportion to gamma.
4. Covariance cross-check: the required quartic coupling equals the
   SdS-frame constant-scalar relation -R_hat/(6*S0^2), reproduced
   independently from the MK-frame computation.

## Verdict (conditional -- the lane's controlling statement)
IF particle masses track a macroscopic classical scalar solving the
conformal EOM in this background, THEN the conformally invariant rotation
observable is the SdS diagnostic and the gamma*r/2 support cancels
(Horne arXiv:1601.07537; Hobson-Lasenby PRD 104, 064014, rendered as
computation). The premise is the live dispute: Mannheim (GRG 54, 99,
2022) holds the mass-generating S is microscopic (vevs varying only
within particle interiors), not a macroscopic classical field. That
counter-position is not computable at this level and is NOT adjudicated.

Test suite after addition: 17 passed.
