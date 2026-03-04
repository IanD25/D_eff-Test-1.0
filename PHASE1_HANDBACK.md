# DS Framework Phase 1 — Execution Handback Document (v2)

## Context
This document summarizes the Phase 1 v2 validation run executed by Claude Code. v1 returned PARTIAL PASS due to 3D cross-route agreement failure on the 25^3 lattice. v2 upgrades the 3D torus to 50^3 (125,000 vertices) with tighter gates [2.85, 3.15] to close that gap.

**Repo:** https://github.com/IanD25/D_eff-Test-1.0
**Code:** `ds_phase1_validation.py` (single file, ~770 lines)
**Output:** `phase1_results/` directory (10 PNG plots + `PHASE1_RESULTS.md`)
**Runtime:** 39.2 seconds on Apple Silicon
**Date:** 2026-03-04
**Git tag:** `phase1-v2`

---

## Overall Verdict: FULL PASS

**Both 2D torus (200x200) and 3D torus (50x50x50) pass all individual route gates AND cross-route agreement gates.**

---

## Detailed Results

### System 1: 2D Torus (200x200) — PASS

| Route | Method | Estimate | True | Gate [1.85, 2.15] | Status |
|-------|--------|----------|------|-------------------|--------|
| 1 | Growth Dimension (BFS) | **1.9754** | 2.0 | PASS | R²=0.999997 |
| 2 | Spectral Dimension (Weyl) | **2.0051** | 2.0 | PASS | R²=0.9977 |
| 3 | Fisher Info Rank (D_eff) | **2.0000** | 2.0 | PASS | All 20/20 samples = rank 2 |

- **Cross-route agreement:** max pairwise diff = **0.030** (gate < 0.10) → **PASS**
- Identical to v1. All three routes converge tightly. Fisher returns exactly 2.0 on every sample.

### System 2: 3D Torus (50x50x50) — PASS

| Route | Method | Estimate | True | Gate [2.85, 3.15] | Status |
|-------|--------|----------|------|-------------------|--------|
| 1 | Growth Dimension (BFS) | **2.8898** | 3.0 | PASS | R²=0.999970 |
| 2 | Spectral Dimension (Weyl) | **2.9399** | 3.0 | PASS | R²=0.9843 |
| 3 | Fisher Info Rank (D_eff) | **3.0000** | 3.0 | PASS | All 20/20 samples = rank 3 |

- **Cross-route agreement:** max pairwise diff = **0.110** (gate < 0.15) → **PASS**
- v1 had max diff = 0.283 at 25^3 → now 0.110 at 50^3. Gate cleared with margin.

### v1 → v2 Comparison (3D only)

| Metric | v1 (25^3) | v2 (50^3) | Change |
|--------|-----------|-----------|--------|
| Growth dimension | 2.7497 | 2.8898 | +0.14 (closer to 3.0) |
| Spectral dimension | 3.0330 | 2.9399 | -0.09 (closer to 3.0) |
| Fisher rank | 3.0000 | 3.0000 | unchanged (perfect) |
| Cross-route agreement | 0.283 (FAIL) | 0.110 (PASS) | **Fixed** |
| Gates | [2.70, 3.30] | [2.85, 3.15] | Tightened |

---

## Implementation Changes in v2

### 1. Analytical Eigenvalues (new in v2)
The v1 `eigsh` sparse solver could not handle the 125k-node 3D Laplacian (hung after 50+ minutes). Since torus eigenvalues have a closed-form expression:

```
λ(k1,...,kd) = 2 * Σ (1 - cos(2πkj/n))   for kj = 0,...,n-1
```

v2 computes all eigenvalues analytically in O(n^d) time. This reduced spectral dimension computation from ~20-40 seconds (eigsh) to **0.0 seconds** and eliminated all convergence issues. Results are mathematically exact.

Note: This analytical shortcut is specific to regular torus lattices. Non-benchmark systems in later phases will need the sparse solver or alternative approaches.

### 2. Gap-Based Fisher Rank Detection (carried from v1 fix)
Same as v1 — the spec's fixed 0.01 threshold fails; gap-based detection works perfectly.

---

## Fisher Singular Value Profiles

### v2 (50^3) vs v1 (25^3):

**2D Torus (unchanged):**
```
SV1: 1.0000  ████████████████████  ← dimension 1
SV2: 1.0000  ████████████████████  ← dimension 2
SV3: 0.2305  █████                 ← noise
SV4: 0.1651  ███                   ← noise
```

**3D Torus v2 (50^3):**
```
SV1: 1.0000  ████████████████████  ← dimension 1
SV2: 1.0000  ████████████████████  ← dimension 2
SV3: 1.0000  ████████████████████  ← dimension 3
SV4: 0.2633  █████                 ← noise
SV5: 0.1653  ███                   ← noise
SV6: 0.1653  ███                   ← noise
```

**3D Torus v1 (25^3) for comparison:**
```
SV4: 0.2700 → 0.2633  (slightly cleaner)
SV5: 0.1717 → 0.1653  (slightly cleaner)
SV6: 0.1717 → 0.1653  (slightly cleaner)
```

Gap ratio improved slightly: 3.7x (v1) → 3.8x (v2). As predicted, the larger torus produces slightly cleaner separation.

---

## Fisher Sigma Sensitivity (v2)

| sigma | 2D Rank | 3D Rank |
|-------|---------|---------|
| 1.5 | 3.0 (unstable) | 4.0 (unstable) |
| 2.0 | **2.0** | **3.0** |
| 3.0 | **2.0** | **3.0** |
| 5.0 | **2.0** | **3.0** |
| 8.0 | **2.0** | **3.0** |

Identical stability profile to v1. Sigma >= 2.0 required for correct results.

---

## What Is Now Fully Validated

1. **Torus construction** — correct for both 2D (40k nodes) and 3D (125k nodes)
2. **Route 1 (Growth dimension)** — converges to true dimension on both systems
3. **Route 2 (Spectral dimension)** — converges to true dimension on both systems
4. **Route 3 (Fisher rank / D_eff)** — **returns exact integer dimension on every sample, both systems**
5. **Cross-route agreement** — all three routes agree within gates on both systems
6. **The core DS claim holds:** D_eff = rank(FIM) correctly recovers intrinsic dimension

---

## Remaining Methodological Notes

1. **Gap-based rank detection is essential** — the spec's fixed threshold does not work. Any implementation of D_eff on non-benchmark systems must use gap-based or similar adaptive thresholding.
2. **Sigma must be sufficiently large** — minimum ~2.0 for lattices of side >= 25. For general graphs, sigma should scale with the graph diameter or be tuned via sweep.
3. **Participation ratio is informative but insufficient** — it trends toward true dimension but overshoots (2D: 2.76; 3D: 4.13 at sigma=3). Gap-based rank is the correct primary estimator.
4. **Analytical eigenvalues are a benchmark-only optimization** — real-world systems will need either sparse solvers with better preconditioning, or alternative spectral methods.

---

## Files in Repo

```
D_eff-Test-1.0/
├── ds_phase1_validation.py          # Full implementation (~770 lines)
├── PHASE1_HANDBACK.md               # This document
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

1. **Phase 1 is complete.** All gates pass. Proceed to Phase 2.
2. The gap-based rank detection and analytical eigenvalue optimizations should be documented as methodological refinements for the DS framework.
3. Phase 2 should test D_eff on systems where the dimension is NOT known a priori — this is where the framework either proves its value or reveals limitations.
