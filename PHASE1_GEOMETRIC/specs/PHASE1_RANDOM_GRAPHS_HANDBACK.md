# Phase 1 Extension: Random Geometric Graphs + ER Negative Control — Handback

**Tag:** `phase1-random-v1`
**Runtime:** 428s
**Code:** `ds_phase1_random_graphs.py`
**Results:** `phase1_results/PHASE1_RANDOM_GRAPHS_RESULTS.md`

---

## 1. RGG Results — Fisher Rank Recovers Embedding Dimension

### RGG d=2 (5000 vertices, avg degree 13.9, diameter 53)

| Route | Estimate | True d | Gate | Verdict |
|-------|----------|--------|------|---------|
| Growth | 1.573 (R²=0.997) | 2 | [1.7, 2.3] | **FAIL** |
| Spectral | 23.13 (R²=0.943) | 2 | [1.6, 2.4] | **FAIL** |
| Fisher rank | mean=2.90, median=3 | 2 | — | Close |
| Fisher PR | mean=3.21 → 2.86 (σ=5) | 2 | [1.5, 2.5] | **FAIL** (converging) |

**Fisher rank distribution:** {1: 1, 2: 4, **3: 44**, 4: 1} — overwhelmingly rank 3, not 2.

**Sigma sweep PR:** 3.21 → 3.02 → 2.86 (σ = 2, 3, 5). Clear downward trend toward 2 but not converged at σ=5. Higher sigma would likely continue descent but diameter is only 53, so σ=5 is only ~10% of diameter — there may be room to push to σ~15 before finite-size effects.

**PR-degree correlation:** r = 0.325 (moderate). This is higher than ideal (should be ~0), suggesting some local-density dependence in the PR.

### RGG d=3 (5000 vertices, avg degree 13.7, diameter 24)

| Route | Estimate | True d | Gate | Verdict |
|-------|----------|--------|------|---------|
| Growth | 2.309 (R²=0.997) | 3 | [2.7, 3.3] | **FAIL** |
| Spectral | 21.21 (R²=0.941) | 3 | [2.5, 3.5] | **FAIL** |
| Fisher rank | mean=3.60, median=4 | 3 | — | Close |
| Fisher PR | mean=3.98 → 3.92 (σ=5) | 3 | [2.5, 3.5] | **FAIL** (converging) |

**Fisher rank distribution:** {1: 1, 2: 2, 3: 13, **4: 34**} — predominantly rank 4, not 3.

**Sigma sweep PR:** 3.98 → 3.95 → 3.92 (σ = 2, 3, 5). Very slow descent. Diameter is only 24, so σ=5 is ~21% of diameter — less room to push.

### RGG Diagnosis

**Growth dimension underestimates** (1.57 vs 2, 2.31 vs 3). This is a known boundary effect on unit-cube RGGs: vertices near the cube boundary have truncated neighborhoods, reducing ball volumes at large r. The effect is systematic, not noise. Fix would be: (a) use periodic boundary conditions (torus topology), or (b) sample only interior vertices for BFS, or (c) increase N substantially.

**Spectral dimension is meaningless** on RGGs. The Weyl law assumes a smooth manifold Laplacian — the random geometric graph Laplacian with its non-uniform degree distribution violates this badly. The "dimension" of 23 is garbage. Route 2 does not apply to irregular graphs.

**Fisher rank overestimates by 1**: consistently returns d+1 rather than d. The most likely explanation is that the non-uniform density of the RGG introduces an additional informative direction in the FIM — the "density gradient" direction. On a torus every vertex has identical local geometry, so the only information comes from the d spatial directions. On an RGG in a bounded cube, moving toward a denser or sparser region creates a distinguishable distributional shift beyond the d geometric directions. This would show up as a (d+1)-th significant singular value.

**Fisher PR is ~d+1 at small σ, converging toward d** as σ increases. This is consistent with the rank observation: at small σ, the density gradient is informative; at larger σ the smoothing averages out local density variations and only the d global directions remain.

---

## 2. ER Negative Control — Surprising Rank-1 Behavior

### ER n=1000 (avg degree 14.9, diameter 4, clustering 0.015)

| Route | Estimate | Expected | Actual |
|-------|----------|----------|--------|
| Growth | 3.617 (R²=1.0) | Very high / poor fit | Degenerate (only 2 fit points) |
| Spectral | 14.79 (R²=0.969) | Poor Weyl fit | High, no physical meaning |
| Fisher rank | **1.00** (all 30 samples) | High (~degree) | **Rank 1 everywhere** |
| Fisher PR | 4.34 | High (~degree) | Moderate (~4) |

### ER n=2000 (avg degree 14.9, diameter 5, clustering 0.007)

Same pattern: Fisher rank = 1.00 (all 20 samples), PR = 4.76.

### ER Diagnosis — Why Rank=1, Not Rank~Degree?

**The spec predicted** ER Fisher rank ≈ degree (high, unstable, no gap), reasoning that each neighbor direction would be independently informative. **The actual result is rank = 1 everywhere.** This is the opposite extreme.

**Explanation:** On an ER graph with diameter 4, the Boltzmann distribution p(x) ∝ exp(-d(v₀,x)/σ) concentrates almost all mass within BFS distance 2 of v₀ (which already covers most of the graph at degree ~15). Moving from v₀ to ANY neighbor w produces essentially the same distributional shift: w's Boltzmann distribution is almost identical to v₀'s because w reaches the same set of vertices in the same number of hops (±1). All score vectors point in approximately the same direction, making the FIM rank 1.

This is fundamentally different from the RGG, where moving to neighbor w₁ (to the north) produces a different distributional shift than moving to neighbor w₂ (to the east), because spatial locality creates direction-dependent changes.

**In other words:** the ER rank=1 result confirms that the Fisher method detects *directional information structure*, not just connectivity. On the ER graph, there is exactly one informative direction: "closer to v₀ vs farther from v₀" (the radial direction in the BFS tree). There are no angular directions because there is no geometry to distinguish different angular positions.

---

## 3. The Diagnostic Contrast

| Metric | RGG d=2 (geometric) | ER n=1000 (non-geometric) |
|--------|--------------------|--------------------|
| Avg degree | 13.9 | 14.9 |
| Clustering | 0.596 | 0.015 |
| **Fisher rank (mean)** | **2.90** | **1.00** |
| Fisher rank (std) | 0.41 | 0.00 |
| Fisher PR (mean) | 3.21 | 4.34 |
| Fisher PR (std) | 0.32 | 0.18 |
| Growth dim | 1.57 | 3.62 |

**The contrast is clear but differently than predicted.**

The spec predicted: RGG rank ≈ d (low), ER rank ≈ degree (high). What we got: RGG rank ≈ d+1, ER rank = 1. The gap-based rank cleanly separates the two cases, but ER goes *below* d rather than *above* it.

The PR comparison is less clean: RGG PR ≈ 3.2, ER PR ≈ 4.3. Both are in the 3-5 range. The PR does NOT provide a clean diagnostic separation — the gap-based rank does.

**The SV profile comparison** (saved as `rgg_vs_er_sv_comparison.png`) shows the structural difference visually: RGG has 2-3 large SVs then a drop; ER has one dominant SV then gradual decay.

---

## 4. Assessment

**Conclusion: (C) — Partial success with important caveats.**

**What works:**
- Fisher gap-based rank correctly identifies the embedding dimension neighborhood (d or d+1) on random geometric graphs, confirming the method extends beyond idealized lattices.
- The ER negative control produces a qualitatively different Fisher signature (rank=1), confirming the method detects geometric structure specifically — not just graph connectivity.
- The diagnostic contrast between RGG and ER is unambiguous via gap-based rank.

**What doesn't work:**
- Growth dimension underestimates due to boundary effects (unit cube, not torus).
- Spectral dimension is completely inapplicable to irregular graphs — the Weyl law assumption is violated.
- Fisher rank returns d+1 rather than d on RGGs, likely due to the density gradient introduced by the bounded domain. This needs investigation: would an RGG on a torus (periodic boundaries) give rank exactly d?
- Fisher PR is not converged at σ=5 and sits above d. The σ range may be too narrow relative to diameter.
- PR-degree correlation is moderate (r=0.3-0.5), suggesting the measurement is somewhat sensitive to local structure.

**Recommendations for next iteration:**
1. **RGG on periodic domain** — scatter points in [0,1]^d with periodic boundary conditions (torus topology). This eliminates boundary effects on both growth dimension and Fisher rank.
2. **Larger N and wider σ sweep** — N=10000+ with σ up to 15-20 to test convergence of Fisher PR toward d.
3. **The d+1 puzzle** — is the extra Fisher direction really the density gradient? Test by comparing interior-only vertices (where density is more uniform) against boundary vertices.
4. **ER at larger n** — the rank=1 result at diameter 4-5 might change at larger diameter. Test n=5000 or 10000 where diameter grows.

---

## 5. Plots Generated

- `rgg_2d_degree_hist.png`, `rgg_3d_degree_hist.png` — degree distributions
- `rgg_2d_growth.png`, `rgg_3d_growth.png` — growth dimension log-log fits
- `rgg_2d_fisher_rank_hist.png`, `rgg_3d_fisher_rank_hist.png` — rank histograms
- `rgg_2d_fisher_pr_hist.png`, `rgg_3d_fisher_pr_hist.png` — PR histograms
- `rgg_2d_pr_vs_degree.png`, `rgg_3d_pr_vs_degree.png` — PR vs vertex degree
- `rgg_2d_pr_vs_clustering.png`, `rgg_3d_pr_vs_clustering.png` — PR vs clustering
- `rgg_2d_sv_profile.png`, `rgg_3d_sv_profile.png` — SV profiles at modal degree
- `rgg_2d_convergence.png`, `rgg_3d_convergence.png` — three-route convergence
- `er_1000_degree_hist.png`, `er_2000_degree_hist.png` — ER degree distributions
- `er_1000_growth.png`, `er_2000_growth.png` — ER growth dimension
- `er_1000_fisher_rank_hist.png`, `er_2000_fisher_rank_hist.png` — ER rank histograms
- `er_1000_fisher_pr_hist.png`, `er_2000_fisher_pr_hist.png` — ER PR histograms
- `er_1000_sv_profile.png`, `er_2000_sv_profile.png` — ER SV profiles
- `rgg_vs_er_sv_comparison.png` — side-by-side SV profiles (the key plot)
- `rgg_vs_er_rank_comparison.png` — overlaid rank histograms
