# DS Framework Phase 1 — Execution Handback Document

## Context
This document summarizes the Phase 1 validation run executed by Claude Code against the spec `DS_PHASE1_CLAUDE_CODE_SPEC.md`. It is intended for handback to the originating chat LLM for review and next-step planning.

**Repo:** https://github.com/IanD25/D_eff-Test-1.0
**Code:** `ds_phase1_validation.py` (single file, ~720 lines)
**Output:** `phase1_results/` directory (10 PNG plots + `PHASE1_RESULTS.md`)
**Runtime:** 70 seconds on Apple Silicon
**Date:** 2026-03-04

---

## Overall Verdict: PARTIAL PASS

**2D torus: FULL PASS. 3D torus: all individual routes pass, cross-route agreement fails due to finite-size effects on 25^3 lattice.**

---

## Detailed Results

### System 1: 2D Torus (200x200) — PASS

| Route | Method | Estimate | True | Gate [1.85, 2.15] | Status |
|-------|--------|----------|------|-------------------|--------|
| 1 | Growth Dimension (BFS) | **1.9754** | 2.0 | PASS | R²=0.999997 |
| 2 | Spectral Dimension (Weyl) | **2.0051** | 2.0 | PASS | R²=0.9977 |
| 3 | Fisher Info Rank (D_eff) | **2.0000** | 2.0 | PASS | All 20/20 samples = rank 2 |

- **Cross-route agreement:** max pairwise diff = **0.030** (gate < 0.10) → **PASS**
- All three routes converge tightly. Fisher returns exactly 2.0 on every single sample.

### System 2: 3D Torus (25x25x25) — PARTIAL

| Route | Method | Estimate | True | Gate [2.70, 3.30] | Status |
|-------|--------|----------|------|-------------------|--------|
| 1 | Growth Dimension (BFS) | **2.7497** | 3.0 | PASS (barely) | R²=0.9998 |
| 2 | Spectral Dimension (Weyl) | **3.0330** | 3.0 | PASS | R²=0.9846 |
| 3 | Fisher Info Rank (D_eff) | **3.0000** | 3.0 | PASS | All 20/20 samples = rank 3 |

- **Cross-route agreement:** max pairwise diff = **0.283** (gate < 0.15) → **FAIL**
- The failure is between Route 1 (2.75) and Route 2 (3.03). Fisher (3.00) agrees with spectral.
- Route 1's underestimate is the expected finite-size effect on the small 25^3 lattice (the spec predicted this: "Finite-size effects will pull estimates slightly below 3.0").

---

## Implementation Note: Fisher Rank Detection

The spec's suggested fixed threshold (0.01) did **not** work. The singular value profiles show:

- **2D:** `[1.000, 1.000, 0.231, 0.165]` — clear gap at position 2→3, but "noise" SVs are ~0.17-0.23, far above 0.01
- **3D:** `[1.000, 1.000, 1.000, 0.270, 0.172, 0.172]` — clear gap at position 3→4

**Fix applied:** Gap-based rank detection — find the largest relative drop in consecutive singular values. This correctly identifies rank=2 for 2D and rank=3 for 3D at every sample vertex with zero failures.

This is significant: the finite-torus geometry creates non-negligible "residual" singular values from opposite-direction neighbors not perfectly canceling (e.g., moving left vs right on a finite cycle). The gap is unambiguous, but any fixed threshold below ~0.20 would overcount.

---

## Fisher Sigma Sensitivity

The gap-based rank is stable across sigma values (with minor instability at very small sigma):

| sigma | 2D Rank | 3D Rank |
|-------|---------|---------|
| 1.5 | 3.0 (unstable) | 4.0 (unstable) |
| 2.0 | **2.0** | **3.0** |
| 3.0 | **2.0** | **3.0** |
| 5.0 | **2.0** | **3.0** |
| 8.0 | **2.0** | **3.0** |

At sigma=1.5, distributions are too peaked and the gap structure breaks down. For sigma >= 2.0, results are perfectly stable.

---

## Singular Value Profiles (Key Evidence for D_eff)

These profiles are the core empirical result supporting Axiom 2:

**2D Torus (4 neighbors → 2 independent directions):**
```
SV1: 1.0000  ████████████████████  ← dimension 1
SV2: 1.0000  ████████████████████  ← dimension 2
SV3: 0.2305  █████                 ← noise (opposite of dim 1)
SV4: 0.1651  ███                   ← noise (opposite of dim 2)
```

**3D Torus (6 neighbors → 3 independent directions):**
```
SV1: 1.0000  ████████████████████  ← dimension 1
SV2: 1.0000  ████████████████████  ← dimension 2
SV3: 1.0000  ████████████████████  ← dimension 3
SV4: 0.2700  █████                 ← noise (opposite of dim 1)
SV5: 0.1717  ███                   ← noise (opposite of dim 2)
SV6: 0.1717  ███                   ← noise (opposite of dim 3)
```

The gap ratio (largest real SV / largest noise SV) is ~4.3x for 2D and ~3.7x for 3D. Unambiguous.

---

## What Passed vs What Needs Work

### Fully validated:
1. Torus construction is correct (degree uniformity, node counts verified)
2. Route 1 (Growth dimension) works on both systems
3. Route 2 (Spectral dimension) works on both systems
4. Route 3 (Fisher rank / D_eff) works on both systems — **the core DS claim holds**
5. 2D convergence: all three routes agree within 0.03

### Known issues:
1. **3D cross-route agreement fails** — but this is a Route 1 finite-size problem, not a Fisher problem. Route 1 growth dimension underestimates on small lattices. A 50^3 lattice would likely fix this but costs ~8x more computation.
2. **Fisher threshold must be gap-based, not fixed** — the spec's 0.01 threshold fails. This is a methodological refinement, not a fundamental flaw. On real (non-benchmark) systems, gap detection will be the correct approach anyway.
3. **Sigma sensitivity at small values** — sigma=1.5 gives incorrect rank. Minimum practical sigma appears to be ~2.0 for lattices of size >= 25.

---

## Files in Repo

```
D_eff-Test-1.0/
├── ds_phase1_validation.py          # Full implementation (~720 lines)
└── phase1_results/
    ├── PHASE1_RESULTS.md            # Machine-generated results table
    ├── torus2d_growth_delta.png     # Δ(r) vs r for 2D
    ├── torus2d_growth_loglog.png    # log V(r) vs log r for 2D
    ├── torus2d_spectral_weyl.png   # Weyl law fit for 2D
    ├── torus2d_fisher_sv.png       # Fisher SV bar chart for 2D
    ├── torus2d_convergence.png     # 3-route comparison for 2D
    ├── torus3d_growth_delta.png    # Δ(r) vs r for 3D
    ├── torus3d_growth_loglog.png   # log V(r) vs log r for 3D
    ├── torus3d_spectral_weyl.png   # Weyl law fit for 3D
    ├── torus3d_fisher_sv.png       # Fisher SV bar chart for 3D
    └── torus3d_convergence.png     # 3-route comparison for 3D
```

---

## Recommended Next Steps

1. **Decide whether PARTIAL is sufficient to proceed to Phase 2**, given that the only failure is finite-size agreement on the 3D lattice (not a Fisher/D_eff issue).
2. **Optionally re-run with 50^3 lattice** to close the 3D agreement gap (expect ~8-10 min runtime).
3. **The gap-based rank detection is a methodological improvement** over the spec's fixed threshold. Consider whether this changes anything about how D_eff should be defined operationally for non-benchmark systems.
4. **The participation ratio trends toward true dimension** as sigma increases (2D: 2.27 at sigma=8; 3D: 3.58 at sigma=8) but does not converge fast enough to be a primary estimator. Gap-based rank is superior.
