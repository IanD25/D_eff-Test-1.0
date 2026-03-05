# Phase 1: Coarse-Graining with Scale-Aware Sigma — Handback

**Tag:** `phase1-coarsegrain-v2`
**Runtime:** 8.3s
**Code:** `ds_phase1_coarsegraining_v2.py`
**Results:** `phase1_results/PHASE1_COARSEGRAIN_V2_RESULTS.md`

---

## 1. The Updated Charts

### 2D Chart: `coarsegrain_2d_deff_vs_level_scaled_sigma.png`
### 3D Chart: `coarsegrain_3d_deff_vs_level_scaled_sigma.png`

Both charts show all routes on a single figure: growth dimension (blue), spectral dimension (red), Fisher PR at fixed σ=3 (green dashed), Fisher PR at scaled σ/side=12.5% Schedule B (green solid), Fisher PR at σ/side=6.25% Schedule C (orange), and Fisher rank at Schedule B (black +). Unreliable levels (σ<2) are shaded in light red.

---

## 2. PR Values at Each Level — Both Protocols

### 2D Torus (128² → 4²)

| Level | Side | σ (fixed) | PR (fixed) | σ (B, 12.5%) | PR (B) | σ (C, 6.25%) | PR (C) |
|-------|------|-----------|-----------|--------------|--------|--------------|--------|
| 0 | 128 | 3.0 | 2.759 | 16.0 | **2.134** | 8.0 | 2.265 |
| 1 | 64 | 3.0 | 2.759 | 8.0 | **2.276** | 4.0 | 2.556 |
| 2 | 32 | 3.0 | 2.767 | 4.0 | **2.579** | 2.0 | 3.157 |
| 3 | 16 | 3.0 | 2.885 | 2.0 | **3.202** | 1.0* | 3.647* |
| 4 | 8 | 3.0 | 3.364 | 1.0* | 3.660* | 0.5* | 2.054* |
| 5 | 4 | 3.0 | 3.991 | 0.5* | 1.991* | 0.25* | 1.114* |

### 3D Torus (32³ → 4³)

| Level | Side | σ (fixed) | PR (fixed) | σ (B, 12.5%) | PR (B) | σ (C, 6.25%) | PR (C) |
|-------|------|-----------|-----------|--------------|--------|--------------|--------|
| 0 | 32 | 3.0 | 4.146 | 4.0 | **3.867** | 2.0 | 4.708 |
| 1 | 16 | 3.0 | 4.323 | 2.0 | **4.776** | 1.0* | 5.050* |
| 2 | 8 | 3.0 | 5.039 | 1.0* | 5.067* | 0.5* | 2.290* |
| 3 | 4 | 3.0 | 5.978 | 0.5* | 2.237* | 0.25* | 1.128* |

---

## 3. Fisher Rank Under Scaled Sigma

### 2D: Rank = 2 at all reliable levels ✓

| Level | Side | σ (B) | Rank (B) |
|-------|------|-------|----------|
| 0 | 128 | 16.0 | **2** |
| 1 | 64 | 8.0 | **2** |
| 2 | 32 | 4.0 | **2** |
| 3 | 16 | 2.0 | **2** |
| 4 | 8 | 1.0* | 3 (unreliable) |

### 3D: Rank = 3 at both reliable levels ✓

| Level | Side | σ (B) | Rank (B) |
|-------|------|-------|----------|
| 0 | 32 | 4.0 | **3** |
| 1 | 16 | 2.0 | **3** |
| 2 | 8 | 1.0* | 4 (unreliable) |

**Rank-based P7 passes perfectly under scaled sigma.** Rank = d at every reliable level.

---

## 4. Is the Scaled-σ PR Flat?

**No.** The scaled-σ PR still increases across levels:

- **2D Schedule B:** 2.134 → 2.276 → 2.579 → 3.202 (change: **+1.068**)
- **3D Schedule B:** 3.867 → 4.776 (change: **+0.909**)
- **P7 monotonicity (PR): FAIL** for both dimensions under scaled sigma.

The scaled-σ protocol **reduces** the PR divergence compared to fixed-σ (the 2D fixed-σ PR goes from 2.759 to 3.991, a change of +1.232, versus +1.068 for Schedule B). But it does NOT eliminate it. The PR is NOT flat.

### Why Scaled-σ Doesn't Fully Fix PR

The root cause is that PR depends on **absolute σ**, not just σ/side. On a 2D torus of any size, PR converges to 2 from above as σ → ∞. At σ=16 (level 0), PR = 2.134 (close to 2). At σ=2 (level 3), PR = 3.202 (further from 2). The σ/side ratio is the same (12.5%) at both levels, but the absolute σ differs by 8×, and the PR at σ=2 is inherently higher than at σ=16 regardless of the lattice size.

In other words: the PR(σ) function on a d-dimensional torus is NOT scale-invariant in σ. The PR at σ/side = 12.5% on a 128² torus is NOT the same as the PR at σ/side = 12.5% on a 16² torus, because the PR also depends on how many lattice spacings σ spans (σ/lattice_spacing), not just σ/diameter.

### The Truly Scale-Invariant Quantity: Rank

The gap-based rank IS scale-invariant: rank = 2 whether σ = 2.0 or σ = 16.0, whether the lattice has 16 or 16,384 vertices. Rank is the correct P7-respecting measure of effective dimension.

---

## 5. At Which Level Does Scaled-σ Break Down?

- **2D Schedule B:** Reliable at levels 0–3 (σ ≥ 2.0). Breaks at level 4 (σ = 1.0).
- **2D Schedule C:** Reliable at levels 0–2 (σ ≥ 2.0). Breaks at level 3 (σ = 1.0).
- **3D Schedule B:** Reliable at levels 0–1 only (σ ≥ 2.0). Breaks at level 2 (σ = 1.0).

The 3D hierarchy is severely limited: only 2 reliable levels with Schedule B. Starting from a larger lattice (e.g., 64³) would give more room.

---

## 6. Assessment

### What this test established:

1. **Fisher rank is perfectly P7-monotone** under scaled sigma. Rank = d at every reliable coarse-graining level. This confirms that the integer D_eff = rank(FIM) respects the Data Processing Inequality.

2. **Fisher PR is NOT P7-monotone** even with scaled sigma. The PR depends on absolute σ, not just σ/side. At small σ (coarse levels), PR inflates above d; at large σ (fine levels), PR converges toward d. This is an inherent property of the participation ratio as a dimension estimator — it's σ-dependent in a way that the rank is not.

3. **The v1 diagnosis was correct but incomplete.** The v1 handback attributed the PR violation to "fixed σ / shrinking lattice" and suggested scaled σ would fix it. Scaled σ reduces the effect but doesn't eliminate it, because PR depends on absolute σ, not the ratio.

4. **Rank is the correct P7-compatible dimension estimator.** PR is useful as a continuous estimator for fractals (where rank ≡ full because there's no gap), but for P7/DPI purposes, rank is the robust quantity.

### Verdict for the paper:

The coarse-graining test confirms:
- **D_eff = rank(FIM) passes P7** — monotonically non-increasing under block-spin coarse-graining, at all reliable measurement scales.
- **PR is a measurement-protocol-dependent continuous extension of rank** — useful for interpolating fractional dimensions but not suitable as a P7 diagnostic without careful σ calibration.

---

## 7. Plots Produced

- `coarsegrain_2d_deff_vs_level_scaled_sigma.png` — the primary 2D chart (replaces Figure 5)
- `coarsegrain_3d_deff_vs_level_scaled_sigma.png` — the 3D chart
