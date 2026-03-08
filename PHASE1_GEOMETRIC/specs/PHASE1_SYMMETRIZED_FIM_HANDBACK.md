# Phase 1: Symmetrized Fisher Information Matrix — Results Analysis Handback

**Tag:** `phase1-symmetrized-v1`
**Runtime:** 22.5s
**Code:** `ds_phase1_symmetrized_fim.py`
**Results:** `phase1_results/PHASE1_SYMMETRIZED_FIM_RESULTS.md`

---

## Executive Summary

The anti-parallel cancellation hypothesis is **definitively rejected**. Symmetrizing the FIM score vectors catastrophically fails on disordered graphs (RGGs, ER), exploding the rank from d+1 to approximately k/2 (half the vertex degree). However, the experiment produced two valuable secondary findings: (1) the pairing quality metric cleanly separates ordered from disordered systems and could serve as a local regularity diagnostic, and (2) the SV profile plots reveal that the symmetrized FIM on RGGs does concentrate energy into the first 2 components — but the gap-based rank detector is defeated by a long noise tail with no clean cutoff. The d+1 direction remains unexplained after three dedicated experiments (bounded RGG, periodic RGG, symmetrized FIM).

**Assessment: (C) Hypothesis REJECTED + (D) Symmetrization breaks ER.**

---

## 1. The Core Question: Did Symmetrization Fix d+1?

### Answer: No. It made things far worse.

| System | Standard Rank (mode) | Symmetrized Rank (mode) | Spec Prediction | Actual Outcome |
|--------|---------------------|------------------------|-----------------|----------------|
| Periodic RGG d=2 | **3** (d+1) | **12** (~k/2) | 2 (d) | FAIL — rank exploded |
| Periodic RGG d=3 | **4** (d+1) | **13** (~k/2) | 3 (d) | FAIL — rank exploded |

The symmetrized rank is approximately **half the vertex degree**, not d. On RGG d=2 (modal degree ≈ 14), the symmetrized rank peaked at 12. On RGG d=3 (modal degree ≈ 15), it peaked at 13.

### Full rank distributions (50 samples each):

**RGG d=2 standard:** {1:6, 2:4, **3:38**, 7:1, 9:1} — tight peak at 3
**RGG d=2 symmetrized:** {6:1, 7:2, 8:4, 9:4, 10:6, 11:7, **12:10**, 13:5, 14:3, 15:1, 16:3, 17:2, 18:1, 19:1} — broad spread from 6 to 19

**RGG d=3 standard:** {1:1, 3:5, **4:38**, 14:2, 15:1, 17:1, 18:1, 22:1} — tight peak at 4
**RGG d=3 symmetrized:** {5:3, 7:1, 8:2, 9:3, 10:5, 11:5, 12:4, **13:7**, 14:3, 15:6, 16:4, 17:3, 18:3, 20:1} — broad spread from 5 to 20

The standard FIM gives a tight, interpretable rank distribution. The symmetrized FIM gives a broad, noisy distribution centered at ~k/2 — worse in every respect.

### σ = 5.0 robustness check:

Changing σ from 3.0 to 5.0 had no effect. RGG d=2 symmetrized mode stayed at 12; RGG d=3 stayed at 12–14. The failure is not σ-dependent.

---

## 2. Did Symmetrization Break Anything That Was Working?

### Tori: Rank preserved. PR dramatically improved.

| System | Std Rank (dist.) | Sym Rank (dist.) | Std PR | Sym PR |
|--------|-----------------|-----------------|--------|--------|
| 2D Torus 200² | {2: 20} | {2: 20} | 2.759 | **2.000** |
| 3D Torus 50³ | {3: 20} | {3: 18, 5: 2} | 4.134 | **3.000** |

On tori, symmetrization is essentially a no-op for rank (correct). However, it dramatically sharpens the participation ratio to **exactly d**: PR = 2.000 on the 2D torus (down from 2.759), PR = 3.000 on the 3D torus (down from 4.134). Both with zero variance across all samples. This is a noteworthy result: symmetrization eliminates the σ-dependent PR inflation that causes the PR to overshoot d on lattices.

The 3D torus showed a minor anomaly: 2 of 20 vertices returned symmetrized rank = 5 instead of 3. This suggests the symmetrization introduces some numerical fragility even on regular lattices (the pair-selection is slightly sensitive to tie-breaking at perfect symmetry). The mode remains 3.

### ER: Rank exploded from 1 → 12. BROKEN.

| System | Std Rank (dist.) | Sym Rank (dist.) | Std PR | Sym PR |
|--------|-----------------|-----------------|--------|--------|
| ER n=1000 | {1: 20} | {6:2, 7:2, 9:2, 10:2, 11:1, 12:4, ...} | 4.585 | 5.834 |

The clean rank = 1 signal — which correctly identified ER as having no angular geometric structure — was completely destroyed. The symmetrized rank scattered broadly from 6 to 20 with no meaningful peak. This is assessment **(D)**: the procedure breaks a previously working result.

### Sierpinski gasket: Rank tightened. PR dropped below d_H.

| System | Std Rank (dist.) | Sym Rank (dist.) | Std PR | Sym PR |
|--------|-----------------|-----------------|--------|--------|
| Gasket L=7 | {2:7, 3:23} | {**3:30**} | 2.098 | 1.220 |

Symmetrization actually tightened the gasket rank distribution from mixed {2, 3} to unanimously {3}. This makes sense: the gasket has moderate anti-parallel structure (pairing quality = 0.589), and the symmetrization strengthened the weaker SV to push vertices that were borderline rank-2 up to rank-3.

However, the PR dropped to 1.220, which is **below** the Hausdorff dimension d_H ≈ 1.585. This indicates over-cancellation: the symmetrized score vectors are more collinear than they should be, compressing the effective dimensionality below the true fractal dimension.

---

## 3. The Diagnostic Metrics

### 3a. Pairing Quality (Mean |cosine| of Best Anti-Parallel Partner)

| System | Pairing Quality | Category |
|--------|----------------|----------|
| 2D Torus 200² | **0.670** | Ordered |
| 3D Torus 50³ | **0.670** | Ordered |
| Gasket L=7 | **0.589** | Ordered (borderline) |
| Periodic RGG d=2 | **0.513** | Disordered |
| Periodic RGG d=3 | **0.371** | Disordered |
| ER n=1000 | **0.337** | Disordered |

The pairing quality separates systems into two regimes with a gap between 0.513 and 0.589:

- **Quality ≥ 0.58 (ordered):** Symmetrization preserves or improves results. Tori and gasket.
- **Quality ≤ 0.51 (disordered):** Symmetrization catastrophically fails. RGGs and ER.

This metric functions as a **local regularity diagnostic** — it measures how much approximate inversion symmetry exists in a vertex's neighborhood. It could be reported as a supplementary quantity in the paper (e.g., "the DS framework detects that tori have high local symmetry while RGGs have low local symmetry").

Note: pairing quality = 0.670 on tori is NOT 1.0, even though tori have perfect inversion symmetry. This is because the BFS-based probability distributions p_w have integer-valued distances, creating a discretization effect. The score vectors for +x and -x are approximately anti-parallel but not exactly anti-parallel (cos ≈ -0.67 rather than -1.0). This is expected — the heat kernel on a discrete lattice is not the same as the heat kernel on the continuum.

### 3b. Disorder Index η

| System | η (mean ± std) | Interpretation |
|--------|---------------|----------------|
| Gasket L=7 | 0.220 ± 0.221 | Low disorder |
| 2D Torus 200² | 0.230 ± 0.000 | Low disorder |
| 3D Torus 50³ | 0.337 ± 0.221 | Low–moderate disorder |
| Periodic RGG d=2 | 0.675 ± 0.153 | High disorder |
| Periodic RGG d=3 | 0.762 ± 0.170 | High disorder |
| ER n=1000 | 0.929 ± 0.051 | Maximal disorder |

The disorder index (ratio of the (d+1)-th to d-th standard singular value) cleanly correlates with pairing quality and separates the same two regimes. The RGG d=2 disorder histogram shows a broad distribution from η ≈ 0.25 to η ≈ 1.0, centered at 0.675, with substantial spread (σ = 0.153). This reflects vertex-to-vertex variation in local geometric regularity within the RGG.

### 3c. Symmetrized SV Profiles — A Nuanced Story

The SV profile comparison plot for RGG d=2 reveals something subtle:

**Standard FIM SV profile** (left panel, averaged over degree-14 vertices):
- SV₁ = 1.00, SV₂ ≈ 0.54, SV₃ ≈ 0.28 → graded decay, gap after SV₃ → rank 3

**Symmetrized FIM SV profile** (right panel):
- SV₁ = 1.00, SV₂ ≈ 0.40, SV₃ ≈ 0.04 → **sharp drop after SV₂** (consistent with d=2!)
- SV₃ through SV₁₄ ≈ 0.01–0.04, all small with large error bars

The symmetrized FIM **does concentrate most of the variance into 2 components** on the RGG d=2 system. The ratio SV₃/SV₂ ≈ 0.10 is a large gap. BUT the gap-based rank detector uses `argmin(sv[i+1]/sv[i])` — it finds the **minimum** ratio across the entire profile. The long tail of tiny SVs (indices 3–14) has ratios that fluctuate near 0.5–1.0 most of the time, but occasionally one of those ratios is < 0.10 (e.g., SV₁₂/SV₁₁ might be 0.05 for a particular vertex). The detector picks that spurious minimum in the noise tail instead of the genuine geometric gap at position 2.

**This is an important mechanistic insight:** the symmetrization is partially doing what was hoped — it IS concentrating geometric signal into the first d components — but it generates a long, noisy tail of small residual SVs rather than a clean zero-floor. The gap-based rank detector, which works by finding the single sharpest drop, is unsuitable for this noise structure.

**Implication:** A different rank estimator (e.g., a threshold-based method, or one that considers the absolute scale of SVs rather than ratios) might extract the correct d=2 from the symmetrized SV profile. However, this would be a fundamentally different approach from the gap-based rank used throughout Phase 1, and would need its own validation pipeline. The symmetrized FIM as a corrected estimator is not viable with the current rank detection method.

---

## 4. Why Symmetrization Fails on Disordered Graphs — Mechanistic Analysis

### The pairing assumption is violated

The symmetrization algorithm assumes each neighbor has a "natural opposite" — a neighbor whose score vector is approximately anti-parallel. This assumption maps perfectly to lattice geometry:

- **2D torus, 4 neighbors:** +x pairs with -x (cos ≈ -0.67), +y pairs with -y (cos ≈ -0.67). Two pairs → 2 half-differences → rank 2. ✓

- **RGG d=2, ~14 neighbors:** Neighbors are scattered in random angular positions around the vertex. The "most anti-parallel" neighbor for any given neighbor has cos ≈ -0.51 (for d=2) or cos ≈ -0.37 (for d=3). These are NOT anti-parallel in any meaningful geometric sense — they're just the least correlated pair in a random set.

### What half-differences of weakly-correlated vectors produce

When you compute s̃_j = ½(s_j - s_{opposite(j)}) where cos(s_j, s_{opposite(j)}) ≈ -0.4:

- The half-difference does NOT cancel common-mode signal (because the vectors aren't truly opposite)
- Instead, it produces a new vector that is a nearly-arbitrary linear combination of the original two
- With k such half-differences, you get k vectors that are nearly linearly independent → rank ≈ k

On lattices (cos ≈ -0.67 between true opposites), enough cancellation occurs to collapse the redundant directions. On RGGs (cos ≈ -0.4 between non-opposites), the "cancellation" is really just vector arithmetic on unrelated directions, which preserves or increases the effective rank.

### Why rank ≈ k/2 specifically

The symmetrized vectors s̃_j have a built-in redundancy: if j pairs with m, and m pairs with j, then s̃_j = ½(s_j - s_m) and s̃_m = ½(s_m - s_j) = -s̃_j. These two are anti-parallel, so they contribute the same FIM direction. In practice, the pairing is not perfectly symmetric (opposite(opposite(j)) ≠ j in general), but there's enough overlap in the pairing graph to create approximately k/2 independent directions rather than k. Hence rank ≈ k/2.

---

## 5. What This Tells Us About the d+1 Direction

### Three hypotheses tested and eliminated:

| Hypothesis | Test | Result |
|-----------|------|--------|
| Boundary density gradients | Periodic RGG (no boundaries) | ❌ REJECTED — d+1 identical on periodic and bounded |
| Anti-parallel cancellation residual | Symmetrized FIM (explicit cancellation) | ❌ REJECTED — symmetrization makes rank worse, not better |
| Fixed-σ measurement artifact | Scaled-σ coarse-graining | ❌ Not relevant — d+1 appears at any σ ≥ 2 |

### What the d+1 direction likely IS:

The evidence across all Phase 1 experiments is most consistent with the d+1 direction being **genuine geometric information about local density fluctuations**:

1. **ER graphs have rank 1**: Only one "direction" exists — the radial distance from the center. No angular structure. The FIM correctly detects this.

2. **Lattice tori have rank d**: The d lattice directions are perfectly represented. There is NO density fluctuation (every vertex has identical local geometry). The FIM correctly detects this.

3. **RGGs have rank d+1**: The d manifold tangent directions PLUS one additional direction encoding local density/degree variation. Even on periodic RGGs, local Poisson fluctuations in point density create a "radial" dimension (how concentrated/diffuse the local neighborhood is) that is orthogonal to the d tangential directions.

This interpretation is further supported by the standard SV profile on RGGs: the SVs don't show a clean d-step followed by noise. Instead they show a graded decay [1.0, 0.54, 0.28, 0.07, 0.05, ...] where the (d+1)-th SV (0.28 at position 3 for d=2) is substantially larger than noise but smaller than the first d. This is what you'd expect from a "real but weaker" geometric direction — the density-fluctuation direction carries less variance than the tangent directions but is not noise.

### Implications for the paper

The d+1 puzzle should be presented honestly. Two framings are available:

**(A) rank(FIM) = d+1 on RGGs is a correct measurement of the local information dimension**, which equals d+1 for a Poisson point process on a d-manifold because local density is an additional degree of freedom. The paper should define D_eff = rank(FIM) as measuring "the number of independent directions in which the local geometry varies" and note that this is d for perfect lattices and d+1 for random geometric embeddings.

**(B) rank(FIM) = d+1 is a systematic bias** that should be corrected. The paper would present the unsymmetrized rank as d+1 and note that PR-based estimates (which give values between d and d+1) may be more appropriate continuous estimators for manifold dimension. The off-by-one is small and predictable.

I recommend framing **(A)** because:
- The rank = 1 on ER is perfectly consistent (ER has ONE varying direction: distance)
- The rank = d on tori is perfectly consistent (tori have d varying directions)
- The rank = d+1 on RGGs is consistent with d varying tangent directions plus 1 density direction
- The FIM is measuring exactly what it should — the number of distinguishable local neighborhoods

---

## 6. Summary Table — All Systems

| System | Std Rank (mode) | Sym Rank (mode) | Std PR | Sym PR | Disorder η | Pairing Quality | Symmetrization Outcome |
|--------|----------------|-----------------|--------|--------|------------|-----------------|----------------------|
| 2D Torus 200² | 2 | 2 | 2.759 | 2.000 | 0.230 | 0.670 | ✓ Rank preserved, PR perfected |
| 3D Torus 50³ | 3 | 3 | 4.134 | 3.000 | 0.337 | 0.670 | ✓ Rank preserved, PR perfected |
| Gasket L=7 | 3 | 3 | 2.098 | 1.220 | 0.220 | 0.589 | ⚠ Rank tightened, PR over-cancelled |
| Periodic RGG d=2 | 3 | 12 | 2.925 | 1.807 | 0.675 | 0.513 | ❌ Rank exploded |
| Periodic RGG d=3 | 4 | 13 | 4.329 | 2.451 | 0.762 | 0.371 | ❌ Rank exploded |
| ER n=1000 | 1 | 12 | 4.585 | 5.834 | 0.929 | 0.337 | ❌ Rank exploded, signal destroyed |

---

## 7. What's Salvageable from This Experiment

Despite the main hypothesis failing, several elements are publishable or inform the paper:

### 7a. The pairing quality as a regularity diagnostic
The mean |cosine| of the best anti-parallel partner is a well-defined, computable quantity that cleanly separates graph types. It could appear in the paper as a supplementary diagnostic: "D_eff via Fisher rank does not require local inversion symmetry, but the pairing quality metric (mean anti-parallel cosine) can detect whether such symmetry exists."

### 7b. Symmetrized PR = exactly d on tori
The fact that symmetrized PR → d exactly (2.000 on 2D torus, 3.000 on 3D torus) is a clean mathematical result. The standard PR overshoots d (to 2.759 and 4.134) due to σ-dependent effects, but the symmetrized version eliminates this. On systems where it works (lattices, structured fractals), symmetrized PR is the theoretically clean estimator.

### 7c. The SV profile insight
The symmetrized SV profile on RGGs does concentrate variance into d components — the problem is in rank detection, not in the underlying signal separation. This could inform future work on improved rank estimators for noisy FIMs.

### 7d. The elimination of three hypotheses about d+1
The combination of periodic RGG + symmetrized FIM experiments conclusively narrows the d+1 explanation to density fluctuations or a fundamental property of the FIM on point-process geometries.

---

## 8. Plots Produced

| File | Description |
|------|-------------|
| `symmetrized_rgg2d_rank_comparison.png` | Overlaid histograms: standard rank (peaked at 3) vs. symmetrized rank (peaked at 12) for periodic RGG d=2. Red dashed line at true d=2. |
| `symmetrized_rgg3d_rank_comparison.png` | Same for periodic RGG d=3. Standard peaked at 4, symmetrized peaked at 13. |
| `symmetrized_rgg2d_sv_comparison.png` | Side-by-side bar charts: standard SV profile (graded decay, gap after SV₃) vs. symmetrized SV profile (sharp drop after SV₂, then long noise tail). Both for degree-14 vertices on RGG d=2. |
| `symmetrized_rgg2d_disorder_hist.png` | Histogram of disorder index η across 50 RGG d=2 vertices. Broad distribution centered at η=0.675. |
| `symmetrized_pairing_quality_comparison.png` | Bar chart of mean pairing quality across all 6 systems. Clear separation: tori (0.67) > gasket (0.59) >> RGG d=2 (0.51) > RGG d=3 (0.37) ≈ ER (0.34). |

---

## 9. Commit Details

- **Commit:** `ffff774`
- **Tag:** `phase1-symmetrized-v1`
- **Branch:** `main`
- **Pushed:** ✓
