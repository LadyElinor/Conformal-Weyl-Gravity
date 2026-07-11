# conformal_probe

Receipts lane for conformal (Weyl) gravity claims. First probe: symbolic
verification that the exact Mannheim-Kazanas exterior solution satisfies the
vacuum Bach equations, with positive controls, negative controls, structural
identity tests, and exact-point numerical cross-validation.

## What this lane verifies (and what it does not)

**Verified here.** The exact Mannheim-Kazanas metric function

    B(r) = 1 - 3*beta*gamma - beta*(2 - 3*beta*gamma)/r + gamma*r - kappa*r^2

(Mannheim & Kazanas, ApJ 342, 635-638, 1989; DOI 10.1086/167623) satisfies
B_ab = 0 for all ten independent components of the Bach tensor. The common
truncation `1 - 2*beta/r + gamma*r - kappa*r^2` does **not**: it violates the
vacuum equations with residuals exactly first order in beta*gamma (e.g.
B_22 = 2*beta*gamma/r^2), which quantifies the sense in which the truncation
is a controlled approximation rather than an exact solution.

**Not verified here.** This is an existence receipt for the vacuum exterior
solution only. It does not touch the contested steps of the conformal-gravity
rotation-curve argument: extended-source integration (Flanagan PRD 74, 023002,
2006; Yoon arXiv:1305.0163; Phillips arXiv:1710.05970), conformal-frame
dependence and geodesic choice (Horne arXiv:1601.07537; Hobson & Lasenby
PRD 104, 064014, 2021, arXiv:2103.13451; EPJC 82, 585, 2022; Mannheim
GRG 54, 99, 2022), interior matching, or uniqueness/completeness of the 1989
solution. Passing this receipt licenses exactly one sentence: the published
exterior metric is an exact vacuum solution of the stated field equations
under the stated conventions.

## Definition and conventions (cited, not from memory)

Bach tensor, taken verbatim from the literature:

    B_ab = 2 * D^m D^n C_{mabn} + C_{mabn} R^{mn}

Sources: arXiv:gr-qc/9811086 eq. (2); arXiv:1105.5632 App. C (Wald
conventions). Properties tested as implementation checks: trace-free
(verified exactly on a metric with nonzero Bach). Weyl machinery is
independently checked by annihilation of a conformally flat FRW metric
with generic scale factor a(t). Symmetry and Bach conformal-invariance
should be treated here as expected theoretical properties, not yet as
independently receipted implementation checks.

Conventions: signature (-,+,+,+); Riemann
R^a_{bcd} = d_c Gamma^a_{bd} - d_d Gamma^a_{bc} + Gamma*Gamma;
Ricci R_bd = R^a_{bad}; standard 4D Weyl decomposition for these conventions.
Note: an implementation using C_{acbd} index order with a 1/2 normalization
computes exactly -1/2 times this tensor; zero-tests are unaffected but
residual values differ by -2. (This lane's v1 did so; v2 matches the cited
convention.)

## Units

The script works in geometrized units, c = 1, lengths in arbitrary consistent
units. Dimension ledger for restoring SI:

    [beta]  = length            beta = G M / c^2
    [gamma] = 1/length          Mannheim-O'Brien fitted universal value:
                                gamma_0 = 3.06e-30 cm^-1
    [kappa] = 1/length^2        Mannheim-O'Brien fitted universal value:
                                kappa = 9.54e-54 cm^-2; cosmological, ~Lambda/3

    Both values receipted 2026-07-08 against the primary abstract,
    arXiv:1007.0970 (Mannheim & O'Brien, PRL 106, 121101, 2011;
    PMID 21517292). Same source gives the bound-orbit galaxy size limit
    gamma_0/kappa = 3.21e23 cm.
    [B(r)]  = dimensionless

Bach components with lowered coordinate indices do not all share the same
coordinate-basis length dimension. For example, angular covariant
components carry the expected extra coordinate-basis powers of r. Uniform
1/length^4 scaling is more naturally associated with orthonormal or
appropriately normalized components, not every printed covariant
coordinate component.

## Files

    bach_receipt.py                     the probe (library + CLI)
    tests/test_bach_receipt.py          pytest wrapper, one test per case
    receipts/RUN_2026-07-08_bach.md     frozen run artifact with environment pin
    requirements.txt                    sympy pin

## Running

    pip install -r requirements.txt
    python3 bach_receipt.py all         # ~25 s total on commodity hardware
    python3 bach_receipt.py mk          # single case: mk|sds|trunc|bad|frw
    pytest -q                           # quick suite (slow test deselected by default)
    pytest -q -m slow                   # constraint-discovery receipt, isolated
    # do not run the full suite in one process on constrained CI; the slow
    # test is order-sensitive under sympy global cache growth

A case reports PASS only if the symbolic verdict and the exact-rational-point
evaluations agree for every component. Expected outcomes are pinned in the
test suite, including the *values* of the truncation residuals, so a silent
convention drift will fail CI rather than silently renormalize.

## Critique dispositions (external review, 2026-07-10)

An 11-point external critique was applied to the frame probe. Dispositions:

**Fixed in this version.**
(4, 11) All claims now use "diagnostic frame-dependence": the precise
statement is that the leading gamma*r/2 contribution to the fixed-mass
static-observer rotation diagnostic is frame-dependent under the
ansatz-preserving conformal map. The package nowhere claims that observed
galaxy rotation curves are gauge artifacts; that adjudication requires the
matter sector. (5) frame_probe.py now prints a scaling report with the
receipted Mannheim-O'Brien values: at galactic radii the leading-order
gamma*r/2 statement holds to ~5e-7 relative error, beta*gamma ~ 4.5e-14,
and the Newtonian/linear crossover falls between 10 and 100 kpc.
(10) Branch structure documented below.

**Documented as standing caveats.**
(2) The conformal map is ansatz-preserving, NOT the full conformal freedom
of the theory. Omega + r*Omega' = Omega^2 is forced only by demanding the
transformed metric stay in static spherical form. The probe renders the
algebraic skeleton of Hobson-Lasenby one notch narrower than their full
argument. (9) SymPy simplification plus pytest is reproducible
computational evidence, not formal proof: simplification artifacts, domain
restrictions (r, rho declared positive; beta, gamma, kappa real), branch
choices, and coordinate singularities are risks mitigated by controls and
exact-point sampling, not eliminated. (8) Nothing in this package moves
conformal gravity closer to matching the full cosmological evidence stack
(CMB spectra, BAO, clusters, lensing, GW, pulsars, solar system); the lane
sharpens one dispute only.

**Scoped as future modules (the decisive ones).**
(1, 3) observable_probe: BUILT (observable_probe.py). Constructs the
invariant: with m(x) = h*S(x) and S of weight -1, the worldline action
S*ds is exactly conformally invariant, so physical trajectories are
geodesics of the invariant mass-frame metric (S/S0)^2 g. On the fixed MK
background the classical conformal scalar EOM admits
S = S0/(1 - a_phys r) -- whose mass frame is the SdS-type metric -- and
does NOT admit constant S, because dR_MK/dr = 6*gamma*(r - 2*beta)/r^3
is nonzero in proportion to gamma itself. CONTROLLING CLAIM (canonical phrasing): under the macroscopic classical scalar realization of
mass generation on this fixed background, the invariant rotation
observable associated with the admitted profile is the SdS diagnostic and
the gamma*r/2 support cancels. This does NOT yet establish a coupled
metric-plus-scalar solution, because the probe currently checks the scalar
field equation on a fixed MK background rather than the full coupled
gravity-plus-matter equations. It also does NOT yet establish uniqueness
of the scalar profile or source-independence of the derived quartic
coupling; in the present construction lambda depends explicitly on beta.
The disputed premise -- Mannheim's
position that S is microscopic, varying only within particle interiors --
is stated, cited (GRG 54, 99, 2022), and NOT adjudicated. Clocks, rods,
and spectral-line bookkeeping beyond the worldline level remain future
work. (6) greens_probe covers the extended-source question
(Phillips/Yoon); the frame probe inherits the point-source exterior form.
(7) lensing_probe: null paths are conformally invariant but observed
angles, angular-diameter distances, and mass calibrations are not
automatically so; compare bending-angle constructions across frames.

## Coherence probes (rigid_mass_probe.py)

Probe 1 renders the rigid-mass axiom's requirement: since
R_MK = 6*beta*gamma/r^2 - 6*gamma/r + 12*kappa has ALL of its r-dependence
proportional to gamma, constant mass in the MK frame demands an anchoring
term with the closed form J(r) = S0*gamma*(r - beta)/r^2 from outside the
macroscopic scalar sector. Probe 2 verifies the strongest fair point on
the Mannheim side: with receipted parameters, R/m^2 ~ 2e-79 (proton) at
10 kpc -- microscopic masses cannot feel galactic curvature. Rigidity is
thereby verified; frame selection is not, and the residue is the
conformal-anomaly / dimensional-transmutation question, outside this
lane's scope.

## Branch structure of the conformal map (point 10)

Two values of a eliminate the linear term. The physical branch,
a = gamma/(3*beta*gamma - 2) (limit -gamma/2 as beta -> 0), yields constant
coefficient exactly +1 -- so NO time renormalization is needed -- preserves
the Newtonian charge (1/rho coefficient equals MK's u), and has its
Omega-singular surface 1 - a*r = 0 at negative r, hence nowhere in the
physical domain r > 0. The rejected branch, a = 1/(3*beta), yields constant
coefficient -1 (not a physically normalizable static form) and is singular
at r = 3*beta. The quadratic coefficient shifts between frames, consistent
with kappa's cosmological (frame-embedding) interpretation.

## Known limitations (carried forward deliberately)

1. Existence plus within-family completeness: the `family` probe shows
   Bach flatness over B = w + u/r + v*r + s*r^2 is equivalent to the single
   constraint 3*u*v = w^2 - 1 (s free), so MK is the general Bach-flat
   member of the polynomial family. Full uniqueness over generic B(r)
   (the ODE reduction of MK 1989) remains unverified.
2. Zero-testing rests on sympy simplification plus exact-point sampling at
   three random rational points per component; this is strong evidence, not
   proof against pathological simplifier failure.
3. Divergence-free identity D^a B_ab = 0: now tested (`div` case), verified
   exactly on a metric with nonzero Bach.
4. Vacuum only; no interior solution, no matching conditions, no sources.

## Planned probes for this lane

- `rotation_probe`: v^2(R) decomposition with universal (gamma_0, kappa)
  against SPARC photometry; per-galaxy chi^2 vs NFW and MOND baselines.
- `frame_probe`: BUILT (frame_probe.py). Derives the ansatz-preserving
  conformal map Omega = 1/(1 - a r) taking MK to Schwarzschild-de Sitter
  form, shows the static-observer rotation diagnostic v^2 = r B'/(2B) exhibits
  diagnostic frame-dependence, with the frame-relocatable piece exactly
  gamma*r/2 at leading order (the linear rotation support of the fits).
  Renders (does not adjudicate) the Hobson-Lasenby vs Mannheim dispute;
  mass-generation modeling deliberately out of scope.
- `greens_probe`: fourth-order weak-field Green function over extended
  sources; sign of the emergent linear coefficient (Yoon/Phillips question).
- `pais_uhlenbeck_probe`: PT inner-product construction for the fourth-order
  oscillator; numerical positive-norm verification of the Bender-Mannheim
  no-ghost theorem (PRL 100, 110402, 2008) within its actual scope.
