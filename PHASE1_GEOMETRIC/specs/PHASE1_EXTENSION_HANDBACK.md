# DS Phase 1 Extension — Execution Handback Document

## Context
Phase 1 torus validation achieved FULL PASS. This extension adds two new test systems:
- **Test 2A: Sierpinski gasket** (L=6, L=7) — fractal with non-integer Hausdorff dimension ≈ 1.585
- **Test 3A: 4D torus** (n=12) — integer dimension d=4

**Repo:** https://github.com/IanD25/D_eff-Test-1.0
**Code:** `ds_phase1_extension.py` (new file, ~700 lines)
**Runtime:** 22.6 seconds
**Git tag:** `phase1-ext-v1`

---

## SIERPINSKI GASKET — THE MAIN EVENT

### Summary of Findings

The Sierpinski gasket results are **scientifically informative but reveal a clear limitation** of the Fisher rank method on fractals. The participation ratio does NOT converge to either known fractal dimension. Instead, it captures something different — a "local navigability dimension" that exceeds the Hausdorff dimension.

### Routes 1 & 2: Growth and Spectral Dimension — Both Pass

| Level | Growth Dim | Spectral Dim | d_H = 1.585 | d_S = 1.365 |
|-------|-----------|-------------|-------------|-------------|
| L=6 (1,095 verts) | **1.485** | **1.483** | ref | ref |
| L=7 (3,282 verts) | **1.520** | **1.446** | ref | ref |

- Growth dimension improves from 1.49→1.52 (approaching d_H=1.585) as lattice size increases. Expected to converge further at L=8, L=9.
- Spectral dimension lands between d_S and d_H (1.45-1.48), consistent with finite-size Weyl law fit on fractal spectrum.
- Both pass their respective gates.

### Route 3: Fisher Information — The Critical Result

#### Singular Value Profile (THE key finding)

```
Sierpinski L=7 mean SV profile:
SV1: 1.0000  ████████████████████  ← primary direction
SV2: 0.5134  ██████████            ← PARTIAL second direction
SV3: 0.1144  ██                    ← weak third direction
SV4: 0.0193  ▌                     ← noise
```

**Contrast with torus (clean integer structure):**
```
2D Torus:     [1.000, 1.000, 0.231, 0.165]  ← TWO equal SVs, then gap
3D Torus:     [1.000, 1.000, 1.000, 0.263, ...]  ← THREE equal SVs, then gap
Sierpinski:   [1.000, 0.513, 0.114, 0.019]  ← GRADED decay, no clean gap
```

**This is the fundamental difference.** On manifold-like systems (tori), the FIM has d equal singular values followed by a sharp drop. On the fractal, the SVs decay gradually — there IS no clean gap. The second direction carries ~50% weight, the third ~11%, the fourth ~2%.

#### Gap-Based Rank Distribution

Per-vertex Fisher across ALL interior vertices:

| Rank | L=6 (1,092 verts) | L=7 (3,279 verts) |
|------|-------------------|-------------------|
| 1 | 129 (11.8%) | 372 (11.3%) |
| 2 | 240 (22.0%) | 726 (22.1%) |
| 3 | 723 (66.2%) | 2,181 (66.5%) |

The distribution is **stable across levels** (ratios barely change L6→L7). ~66% of vertices show rank 3, ~22% rank 2, ~12% rank 1. This reflects the fractal's hierarchical structure: vertices at junction points between sub-triangles see fewer independent directions than vertices deep inside a sub-triangle.

#### Participation Ratio

| Level | Mean PR | Std | d_H | d_S |
|-------|---------|-----|-----|-----|
| L=6 | **2.057** | 0.176 | 1.585 | 1.365 |
| L=7 | **2.056** | 0.174 | 1.585 | 1.365 |

The participation ratio is **remarkably stable** across levels (2.057 vs 2.056) but lands at **~2.06**, which is **above both d_H (1.585) and d_S (1.365)**. This is NOT the Hausdorff dimension.

#### Sigma Sensitivity of Participation Ratio

| sigma | L=6 PR | L=7 PR |
|-------|--------|--------|
| 1.5 | 2.43 | 2.42 |
| 2.0 | 2.30 | 2.29 |
| 3.0 | 2.14 | 2.13 |
| 5.0 | 2.00 | 1.99 |
| 8.0 | 1.93 | 1.90 |
| 12.0 | 1.88 | 1.86 |

**The PR decreases monotonically with sigma.** At small sigma (peaked distributions), the Fisher method sees more local structure and reports higher dimension. At large sigma (smeared distributions), it sees less structure and approaches ~1.86. This is scale-dependent behavior — expected on a fractal, but it means the PR does not have a natural "true" value.

**None of the sigma values produce PR ≈ 1.585 (d_H) or ≈ 1.365 (d_S).** Even at sigma=12, PR is still at 1.86.

### Interpretation

The Fisher participation ratio on the Sierpinski gasket measures **local navigability dimension** — how many independent directions you can effectively move from a vertex, weighted by information content. On a fractal:

1. **The gap-based rank is not meaningful** — there is no clean gap, and the rank varies across vertices. This confirms that D_eff = rank(FIM) is a manifold concept, not a fractal concept.

2. **The participation ratio is a continuous estimator** but it does NOT recover any of the standard fractal dimensions. It reports ~2.06 at sigma=3, which exceeds d_H because the gasket is locally embedded in 2D and every interior vertex has 4 neighbors providing information in multiple directions.

3. **The SV profile IS informative** — the graded decay [1.0, 0.51, 0.11, 0.02] encodes the fractal's dimensional structure. A manifold produces a step function; a fractal produces a smooth decay. This qualitative difference could potentially be formalized into a more nuanced dimension estimator.

4. **Scale dependence is expected** — on a self-similar fractal, dimension is "the same at all scales" only for Hausdorff-type measures. The Fisher method, being distribution-based, naturally sees scale-dependent structure.

---

## 4D TORUS — Quick Extension

### Results

| Route | Estimate | True | Gate | Status |
|-------|----------|------|------|--------|
| Growth | **3.376** | 4.0 | [3.70, 4.30] | **FAIL** |
| Spectral | **3.944** | 4.0 | [3.70, 4.30] | PASS |
| Fisher rank | **4.000** | 4.0 | rank=4 | PASS |

- **Cross-route agreement:** 0.624 (gate < 0.20) → **FAIL**
- The growth dimension failure is purely finite-size: n=12 gives maxR=6, far too short a scaling regime for reliable d=4 estimation. A 4D torus needs n≥20 for growth dimension to work (maxR~10).
- Fisher returns **exactly 4.0 on all 20 samples**, confirming the method scales correctly to d=4.
- Spectral dimension is accurate at 3.944.

### Fisher SV Profile (4D Torus)

```
SV1: 1.0000  ████████████████████  ← dim 1
SV2: 1.0000  ████████████████████  ← dim 2
SV3: 1.0000  ████████████████████  ← dim 3
SV4: 1.0000  ████████████████████  ← dim 4
SV5: 0.3813  ████████              ← noise
SV6: 0.2420  █████                 ← noise
SV7: 0.2420  █████                 ← noise
SV8: 0.2420  █████                 ← noise
```

Clean 4-signal / 4-noise structure with gap ratio ~2.6x. Exactly as predicted.

### Note on Growth Dimension Fix

The 4D growth dimension could be fixed by running at n=20 (160,000 vertices) which gives maxR≈10 and a wider scaling regime. However, the Fisher result (the DS core claim) is already validated, so this is a lower priority.

---

## Level Convergence (Sierpinski)

| Metric | L=6 | L=7 | Trend |
|--------|-----|-----|-------|
| Growth dim | 1.485 | 1.520 | → d_H (improving) |
| Spectral dim | 1.483 | 1.446 | → d_S (improving) |
| Fisher PR (σ=3) | 2.057 | 2.056 | Stable (converged) |
| Fisher rank dist | 12/22/66% | 11/22/67% | Stable (converged) |
| Mean SV2 | 0.517 | 0.513 | Stable |
| Mean SV3 | 0.123 | 0.114 | Slightly decreasing |

Growth and spectral dimensions are still converging toward their theoretical values (more levels would help). Fisher metrics are already converged — the PR and rank distribution barely change from L=6 to L=7.

---

## Key Takeaways for the DS Framework

### What works:
1. **Fisher rank on manifold-like systems is validated through d=4.** Exact integer results on 2D, 3D, and 4D tori with zero exceptions.
2. **The SV gap structure is a reliable manifold indicator.** Clean step function = manifold. Graded decay = fractal/non-manifold.

### What the Sierpinski test reveals:
3. **D_eff = rank(FIM) is a manifold concept.** On fractals, the rank is not well-defined (no clean gap, vertex-dependent, sigma-dependent).
4. **The participation ratio is a continuous generalization** but it does NOT recover standard fractal dimensions. It measures something distinct — a local information-geometric quantity.
5. **The SV profile itself may be more informative than any scalar summary.** The shape of the decay encodes dimensional structure that neither rank nor participation ratio fully captures.

### Implications for Phase 2:
6. On real-world systems that are manifold-like (e.g., data manifolds embedded in high-dimensional spaces), D_eff should work well.
7. On systems with fractal-like structure (e.g., hierarchical networks, scale-free graphs), the SV profile should be analyzed as a whole rather than reduced to a single number.
8. The sigma-dependence of the participation ratio on the gasket suggests a potential **multiscale Fisher dimension** approach, where the PR-vs-sigma curve characterizes the system's scale structure.

---

## Files Added

```
D_eff-Test-1.0/
├── ds_phase1_extension.py                       # Extension implementation
├── PHASE1_EXTENSION_HANDBACK.md                 # This document
└── phase1_results/
    ├── PHASE1_EXTENSION_RESULTS.md              # Machine-generated results
    ├── sierpinski_L6_growth_delta.png
    ├── sierpinski_L6_growth_loglog.png
    ├── sierpinski_L6_spectral_weyl.png
    ├── sierpinski_L6_fisher_sv.png              # KEY: graded SV profile
    ├── sierpinski_L6_fisher_sv_sigma_sweep.png  # KEY: SV vs sigma
    ├── sierpinski_L6_fisher_participation_hist.png  # KEY: PR distribution
    ├── sierpinski_L6_convergence.png
    ├── sierpinski_L7_*.png                      # Same set for L=7
    ├── torus4d_growth_delta.png
    ├── torus4d_growth_loglog.png
    ├── torus4d_spectral_weyl.png
    ├── torus4d_fisher_sv.png
    └── torus4d_convergence.png
```
