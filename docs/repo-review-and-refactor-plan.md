# Conformal Weyl Gravity Repo Review and Refactor Plan

Status: first-pass maintainer note
Date: 2026-07-10

## Executive summary

This repository is stronger than most controversial-theory codebases because it already has real receipt discipline, unusually careful claim scoping, and tests that pin symbolic identities rather than merely exercising code paths.

Its clearest strength is the `bach_receipt.py` lane, which verifies a precise mathematical claim with positive controls, negative controls, structural identities, and exact-point cross-checks. Its main weakness is architectural: too much of the repo's intellectual weight currently lives in print-heavy one-file scripts, which makes the boundary between computation, receipt generation, and narrative interpretation harder to defend than it needs to be.

The right direction is not to turn this into a generic library first. The right direction is to remain receipt-first while extracting clean computational kernels and structured result objects so narrative claims cannot blur into symbolic results by accident.

## Current strengths

- Good epistemic hygiene for a contested subject.
- Strongest artifact is `bach_receipt.py`.
- Tests pin exact identities and residual values, not just pass/fail status.
- The repo usually distinguishes:
  - exact mathematical verification,
  - diagnostic consequences,
  - conditional physical interpretation,
  - unresolved or disputed premises.
- Existing receipts preserve work products in a human-auditable form.

## Current weaknesses

- The repository sits awkwardly between executable research note and reusable package.
- Probe logic, CLI behavior, and argumentative printing are too interwoven.
- `tests/` currently rely on `sys.path` insertion rather than package imports.
- There is no machine-readable receipt format alongside markdown receipts.
- The interpretive probes, especially `observable_probe.py`, carry more rhetorical attack surface than necessary because computational result and premise-sensitive framing are still too close together.

## Recommended architectural stance

Stay receipt-first.

Do not optimize first for broad library reuse. Instead:
- preserve the current probe files and workflows,
- extract stable symbolic kernels under a package namespace,
- introduce structured result objects,
- move CLI/report formatting into thin wrappers,
- add JSON receipts only after the computational/core boundaries are clearer.

## Refactor goals

1. Preserve existing scientific scope and receipt workflows.
2. Sharpen the boundary between computation, receipt, and narrative.
3. Make each probe callable as importable code, not just as a print-heavy script.
4. Standardize conditional-claim structure across probes.
5. Prepare for CI and machine-readable receipts without destabilizing the current repository.

## Recommended phased plan

### Phase 0: freeze intent
Write and keep this note in-repo so the architectural direction is explicit before code moves begin.

### Phase 1: add package skeleton without breaking current runners
Create:
- `pyproject.toml`
- `conformal_probe/`
- `conformal_probe/core/`
- `conformal_probe/probes/`
- `conformal_probe/reports/`
- `conformal_probe/claims.py`

This phase should not remove or rename any existing top-level probe scripts.

### Phase 2: first safe extraction
Extract only the reusable symbolic geometry engine (`DiagMetric` and helpers) into a package module such as `conformal_probe/core/geometry.py`.

Do not split too aggressively at first. A single `geometry.py` is safer than prematurely separating `metric.py` and `tensors.py`.

### Phase 3: structured result objects
Each probe should return a structured result object with explicit fields separating:
- verified fact,
- claim class,
- premises,
- disputed premises,
- conditional inference,
- non-claims,
- limitations,
- pinned expressions,
- sources.

### Phase 4: reporting split
Move human-readable markdown receipt formatting and future JSON receipt export into `conformal_probe/reports/`.

Legacy top-level scripts should become thin compatibility wrappers over the package APIs.

### Phase 5: CI and receipt hardening
Add:
- quick test suite on every push/PR,
- isolated slow symbolic suite,
- optional JSON receipt snapshots or pinned-hash comparisons.

## Claim-shape recommendation

Every probe should eventually expose a common structure:

- **Verified fact**
- **Claim class** (`exact_math`, `diagnostic`, `conditional_physics`, `coherence_requirement`)
- **Premises**
- **Disputed premises**
- **Conditional inference**
- **Not adjudicated**
- **Limitations**
- **Sources**

This is especially important for `observable_probe.py`, where the computation is meaningful but premise-sensitive.

## Immediate next implementation step

The first concrete implementation pass should be intentionally modest:
- add package skeleton,
- add `pyproject.toml`,
- add claim/result dataclasses,
- leave current probe scripts and tests working as they are.

That gives the repo a durable target shape without forcing risky code motion in the same change.
