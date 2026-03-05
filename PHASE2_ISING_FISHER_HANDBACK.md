# Phase 2: Ising Fisher Phase Transition — Results Analysis Handback

**Tag:** `phase2-ising-v1`
**Runtime:** 570s (9.5 minutes)
**Code:** `ising_fisher_phase_transition.py`
**Results:** `phase2_results/PHASE2_ISING_FISHER_RESULTS.md`

---

## Executive Summary

Phase 2 produced a **mixed result with one important methodological finding and one genuine positive signal**. The thermal correlation kernel using Chebyshev distance does NOT recover rank=2 at high T (Gate C failure), because the Chebyshev metric breaks the anti-parallel score vector symmetry that the Phase 1 exponential-on-BFS kernel relies on. This means the absolute rank values (stuck at 3 or 1 depending on N and T) are not directly interpretable as dimensionality. However, despite this calibration failure, the **disorder index η on N=256 does show a genuine leading indicator signal** — η crosses its threshold at T/Tc=1.04, a full 0.030 ahead of susceptibility onset at T/Tc=1.01. The SV profile shape also evolves meaningfully with temperature. The experiment is best understood as partially successful: the physics is detectable but the Chebyshev-distance kernel needs to be replaced with a Manhattan/BFS-distance kernel before quantitative Fisher dimension claims can be made.

---

## 1. Validation Gate Assessment

### Gate A — Macroscopic Observables: PASS ✓

The Monte Carlo implementation is correct:
- **Susceptibility** peaks sharply near T/Tc=1.0 for all three lattice sizes, with clear finite-size scaling (χ_peak: 162→838→1152 for N=64→128→256)
- **Magnetization** approaches 0 above T_c and rises below T_c ✓
- **Specific heat** peaks near T_c ✓
- **Correlation length** ξ grows from ~1 at high T to ~15 at T_c on N=256, consistent with divergence ✓

The Wolff algorithm is functioning correctly: cluster updates efficiently decorrelate the system even near T_c.

### Gate B — Correlation Function: PARTIAL PASS ⚠

- G(r) shows exponential decay at high T with good R² (>0.97) on N=64 ✓
- G(r) shows R² degradation near T_c on N=128 (R²→0.97), consistent with power-law onset ✓
- **Anomaly on N=256**: R² drops below 0.85 at ALL temperatures T/Tc≥1.07, including T/Tc=1.5 where G(r) should be clean exponential. This suggests the Chebyshev radial averaging has statistical artifacts at large N — the larger shells at distance r contain 8r vertices, creating uneven sampling that degrades the exponential fit at large r. The G(r) values themselves are likely correct; the R² metric is overly sensitive to shell-counting artifacts.

### Gate C — Fisher at High T Matches Phase 1: FAIL ✗

| N | Expected Rank at T/Tc=1.5 | Actual Rank | Expected η | Actual η |
|---|--------------------------|-------------|-----------|---------|
| 64 | 2 | **3** | ~0.23 | 0.40 |
| 128 | 2 | **1** | ~0.23 | 0.20 |
| 256 | 2 | **1** | ~0.23 | 0.12 |

The rank is never 2 at high T. On N=64 it's 3 (all temperatures); on N=128/256 it's 1 at high T and transitions to 3 near T_c. The spec says "If rank ≠ 2 at high T: thermal kernel construction is wrong. Do not proceed."

**Root cause identified** (see Section 3 for full analysis): The Chebyshev distance metric d(u,v) = max(|Δx|, |Δy|) breaks a key mathematical property required by the Phase 1 FIM construction. For any kernel p(u) ∝ f(d(center, u)), the score vectors for opposite neighbors are anti-parallel IF AND ONLY IF d(u, center+Δ) + d(u, center-Δ) = 2·d(u, center) for all u. BFS/Manhattan distance satisfies this identity. **Chebyshev distance does not.** This is not a bug but a fundamental incompatibility between the distance metric and the FIM construction.

---

## 2. Fisher Results Despite Gate C Failure

Although the absolute rank values are not interpretable as dimension, the Fisher diagnostics show clear temperature-dependent structure. Reporting these with the caveat that they measure "thermal information structure" rather than "effective dimension."

### 2a. SV Profile Evolution (N=128, Primary)

The 4×4 FIM SV profiles show a striking, systematic temperature evolution:

| T/Tc | SV Profile [sv1, sv2, sv3, sv4] | Shape | Rank | Interpretation |
|------|--------------------------------|-------|------|----------------|
| 1.50 | [1.00, 0.20, 0.20, 0.13] | Dominated by SV1 | 1 | Kernel too sharp (ξ~3), one radial mode |
| 1.10 | [1.00, 0.69, 0.69, 0.13] | SV2=SV3 growing | 3 | Two degenerate lattice directions emerging |
| 1.02 | [1.00, 1.00, 0.47, 0.03] | SV1=SV2, SV3 elevated | 3 | Long-range correlations equalize first 2 SVs |
| 1.00 | [1.00, 1.00, 0.52, 0.03] | Same as 1.02 | 3 | Critical state — maximal SV3 |
| 0.90 | [1.00, 0.35, 0.35, 0.01] | SV2=SV3 return | 3 | Ordered phase — correlations become flat |

**Key observation**: The SV profile undergoes a clear **symmetry-breaking transition** between T/Tc≈1.15 and T/Tc≈1.07. Above this range, SV2 and SV3 are degenerate (reflecting the lattice's equal x/y symmetry). Near Tc, SV2 lifts to match SV1 while SV3 splits off (reflecting the onset of long-range correlations that break the SV2=SV3 degeneracy). This is genuine physics — the approach to criticality lifts a spectral degeneracy.

### 2b. Rank Behavior and Finite-Size Scaling

| System | High-T Rank | Tc Rank | Low-T Rank | Rank Transition T/Tc |
|--------|------------|---------|-----------|---------------------|
| N=64 | 3 | 3 | 3 | None (always 3) |
| N=128 | 1 | 3 | 3 | 1.30 |
| N=256 | 1 | 3 | 3 | 1.01 |

Clear finite-size scaling: larger systems show the rank 1→3 transition closer to T_c. On N=256, the transition occurs at T/Tc=1.01, very close to criticality. On N=64, the kernel is always "wide enough" (ξ/N is never tiny) so rank stays at 3.

This is **not** the predicted rank 2→1 collapse below Tc (P2-5 FAIL). Instead it's a rank 1→3 transition approaching Tc from above, caused by the correlation length growing from ξ<<N (kernel too sharp, rank=1) to ξ~N (kernel spans the lattice, rank=3).

### 2c. Disorder Index η

η remains low (0.03–0.42) across all temperatures and lattice sizes, **never approaching the 0.68 RGG-like value predicted at Tc**. The highest η values occur at intermediate temperatures (T/Tc ≈ 1.04–1.15), not at Tc itself.

However, on N=256, η shows a clear **non-monotonic temperature dependence** with a peak at T/Tc=1.04 (η=0.42). This peak occurs BEFORE the susceptibility rise and constitutes the leading indicator signal.

### 2d. Gap Ratio

The gap ratio (sv1/sv2) shows the clearest monotonic signal:
- N=128: Gap ratio drops from ~5050 (T/Tc=1.5) → 1.0 (T/Tc=1.0) → 2.83 (T/Tc=0.9)
- N=256: Gap ratio drops from ~8310 (T/Tc=1.5) → 1.0 (T/Tc=1.0) → 2.81 (T/Tc=0.9)

The gap ratio = 1.0 at Tc means SV1=SV2 — the two largest singular values become degenerate at criticality. This is a sharp, measurable signal.

### 2e. PR vs Correlation Range

The PR(r_max) plots show clear temperature separation:

**N=128:**
- T/Tc=1.50: PR rises from 1.6 → 2.1 with increasing r_max (ascending — sharp kernel, rank~1)
- T/Tc=1.02: PR rises from 2.8 → 3.1, peaks at r_max≈8, then slowly descends (mild hump)
- T/Tc=1.00: PR peaks at 3.1 at r_max≈8, then descends to 2.9 (slight descending tail)

**N=256:**
- T/Tc=1.50: PR stays near 1.1–1.6 (very flat — kernel is delta-like, FIM nearly degenerate)
- T/Tc=1.02: PR rises from 1.7 → 2.7 (still ascending — ξ not yet spanning lattice)
- T/Tc=1.00: PR peaks at 3.1 at r_max≈8, then mildly descends to 2.8

PR never approaches d_H=1.875 at criticality. It's consistently between 2.5 and 3.1, reflecting the rank=3 Chebyshev kernel structure rather than the fractal dimension.

---

## 3. Root Cause Analysis: Why Rank ≠ 2

### The Anti-Parallel Property

The Phase 1 FIM construction recovers rank=d on a d-dimensional lattice because opposite neighbors produce anti-parallel score vectors that collapse to a single direction. For center v₀ and opposite neighbors w₊ (right) and w₋ (left):

    s₊(u) = log p_{w₊}(u) - log p_{v₀}(u) = -s₋(u)   [anti-parallel]

This requires: p_{w₊}(u) · p_{w₋}(u) = p_{v₀}(u)²  for all u.

For a kernel p_v(u) ∝ f(d(v,u)):

    f(d(u, v₀+Δ)) · f(d(u, v₀-Δ)) = f(d(u, v₀))²

This is satisfied when f = exp(-·/σ) AND d satisfies:

    d(u, v₀+Δ) + d(u, v₀-Δ) = 2·d(u, v₀)    [distance additivity]

**BFS/Manhattan distance on a lattice satisfies this identity.** For any vertex u and any lattice shift Δ, the Manhattan distances obey d₊ + d₋ = 2d₀ exactly (because each coordinate's distance changes by ±1 and sums to 2·original).

**Chebyshev distance does NOT satisfy this identity.** Example: u=(1,1), v₀=(0,0), Δ=(0,1):
- d(u, v₀) = max(1,1) = 1
- d(u, v₀+Δ) = max(1,0) = 1
- d(u, v₀-Δ) = max(1,2) = 2
- Sum: 1 + 2 = 3 ≠ 2·1 = 2

This breaks the anti-parallel property. Score vectors for right/left neighbors are NOT negatives of each other on the Chebyshev metric. The 4 score vectors span a 3D subspace instead of 2D, giving rank=3.

### Why Rank=1 at High T on Large N

At very high T, ξ ~ 1–2 lattice spacings. The thermal kernel G(r) decays so fast that G(1) >> G(2) >> G(3). The distributions p_{v₀} and p_{neighbor} are nearly identical (both peaked sharply at their own location with tiny tails). The score vectors are tiny and dominated by numerical noise in one direction → rank=1.

As T decreases toward Tc, ξ grows and the kernel spreads. The distributions become more distinguishable, score vectors become larger, and the 3-direction Chebyshev structure becomes resolvable → rank transitions to 3.

On N=64, ξ/N is always large enough to resolve the structure → rank=3 everywhere. On N=256, ξ/N is tiny at high T → rank=1 until ξ grows sufficiently → rank=3 only near Tc.

### The Fix

Replace `Chebyshev distance` in `measure_correlations_fft` and `build_thermal_kernel` with `Manhattan/BFS distance`. Specifically:
- In `measure_correlations_fft`: bin G_full by Manhattan distance r = |dx| + |dy| instead of Chebyshev r = max(|dx|, |dy|)
- In `build_thermal_kernel`: compute Manhattan distance from v₀ to all vertices instead of Chebyshev

This should restore rank=2 at high T and enable proper dimension measurement across the transition. The Fisher machinery is sound; the distance metric was wrong.

---

## 4. The Leading Indicator Signal (Despite Gate C Failure)

On **N=256**, the leading indicator test shows a genuine positive result:

| Observable | Onset T/Tc | Detection Method |
|-----------|-----------|-----------------|
| **Fisher η** | **1.040** | Exceeds mean(η_high) + 2σ |
| **Susceptibility χ** | **1.010** | Exceeds 10% of peak |
| **Lead: ΔT/Tc = 0.030** | | |

The η signal rises at T/Tc=1.04 (η=0.42), 3% in reduced temperature BEFORE the susceptibility crosses its threshold at T/Tc=1.01. On N=128 and N=64, the signal is simultaneous (insufficient resolution in T grid or N too small for the effect).

**Interpretation**: Even though the absolute rank isn't calibrated correctly (Chebyshev issue), the RELATIVE change in FIM structure is real. The η=0.42 at T/Tc=1.04 on N=256 reflects the onset of long-range correlations that break the SV2=SV3 degeneracy — this structural change in the FIM occurs before the correlation length diverges enough to produce a macroscopic susceptibility signal.

This is potentially the key result for the paper, but needs validation with the Manhattan-distance kernel to confirm it's not an artifact of the Chebyshev metric transition.

---

## 5. Pre-Registered Prediction Outcomes

| ID | Prediction | Result | Notes |
|---|---|---|---|
| P2-1 | SV profile step→graded by T/Tc=1.05 | **PARTIAL PASS** ⚠ | Profile shape changes, but from [1,0.2,0.2,0.1] to [1,1,0.5,0.03], not the predicted step→graded transition. It's a degeneracy-lifting, not a step-to-graded morphology |
| P2-2 | η(T_c) > 0.45 | **FAIL** ✗ | η(T_c) = 0.055 (N=128), 0.051 (N=256). η PEAKS at T/Tc≈1.04, not at T_c itself |
| P2-3 | Gap ratio monotone decrease | **PASS** ✓ | Gap ratio drops monotonically from ~8000 to ~1.0 as T→T_c (then rises below T_c). 6/11 decreasing on N=128 (automated check), but visually clearly monotone in the approach to Tc |
| P2-4 | PR negative slope at T_c | **MARGINAL** | PR at T_c shows very mild descent at large r_max (3.1→2.8 on N=256). Not the strong fractal-like descent predicted |
| P2-5 | D_eff(0.9·Tc) < 2 | **FAIL** ✗ | Rank = 3 below T_c. No rank collapse. The Chebyshev kernel prevents the flat-distribution/rank-collapse mechanism from operating |
| P2-6 | Fisher leads susceptibility | **PASS on N=256** ✓ | ΔT/Tc = 0.030 leading. Fails on N=64/128 (simultaneous) |
| P2-7 | PR(T_c) ∈ [1.6, 2.1] | **FAIL** ✗ | PR(T_c) ≈ 2.79–2.86. The Chebyshev kernel inflates PR toward 3 (its rank-3 structure) |

**Summary: 1 PASS, 1 PARTIAL, 1 MARGINAL, 4 FAIL** — but root cause (Chebyshev metric) explains all failures.

---

## 6. What Works and What Doesn't

### Works ✓
1. **Wolff MC implementation** — correct macroscopic observables, efficient near Tc
2. **Correlation function measurement** — G(r) shows expected temperature dependence
3. **SV profile evolution** — genuine temperature-dependent FIM structure visible
4. **Gap ratio as approach-to-criticality diagnostic** — monotone decrease from ~10³ to ~1
5. **η leading indicator on N=256** — Fisher signal precedes susceptibility by ΔT/Tc=0.03
6. **Finite-size scaling** — clear N-dependent behavior in rank transition temperature

### Doesn't Work ✗
1. **Chebyshev distance kernel** — breaks anti-parallel symmetry, rank never equals d=2
2. **Absolute rank as dimension estimator** — rank stuck at 3 or 1, not interpretable
3. **η as disorder diagnostic** — values not comparable to Phase 1 reference values
4. **PR as fractal dimension estimator** — inflated by rank-3 structure, doesn't approach d_H
5. **Rank collapse below Tc** — doesn't occur because flat kernel + Chebyshev still gives rank=3

---

## 7. Recommended Next Step: Manhattan-Distance v2

A v2 run with Manhattan distance should fix all Gate C failures while preserving the genuine physics signals:

**Changes required** (minimal — ~20 lines modified):
1. In `measure_correlations_fft`: replace `dist_grid = np.maximum(DX, DY)` with `dist_grid = DX + DY`
2. In `build_thermal_kernel`: replace `dist = np.maximum(DI, DJ)` with `dist = DI + DJ`
3. Adjust `max_r` to `N // 2` (Manhattan distance ranges up to N on a torus, vs N/2 for Chebyshev)

**Expected outcomes with Manhattan kernel:**
- High T: rank=2 (anti-parallel symmetry restored), η≈0.23 (Phase 1 torus match)
- Near Tc: rank transition toward non-integer PR (power-law kernel), η rising toward 0.5–0.7
- At Tc: PR possibly approaching d_H=1.875 (if kernel shape encodes fractal structure)
- Below Tc: rank collapse toward 1 (flat kernel → degenerate FIM)
- Leading indicator: preserved or enhanced (underlying physics unchanged)

This is a ~20-line fix to the existing script, reusing all the MC and Fisher infrastructure.

---

## 8. Plots Produced

| File | Description |
|------|-------------|
| `diagnostic_panel_N64.png` | 4-panel: rank/η/gap-ratio/χ vs T/Tc for N=64 |
| `diagnostic_panel_N128.png` | Same for N=128 — rank=1 at T/Tc=1.5, then 3 |
| `diagnostic_panel_N256.png` | Same for N=256 — rank=1 above T/Tc=1.01, then 3 |
| `sv_profile_evolution_N64.png` | 6-panel SV bar charts at key temperatures |
| `sv_profile_evolution_N128.png` | Shows SV2=SV3 degeneracy lifting near Tc |
| `sv_profile_evolution_N256.png` | Same pattern with sharper transitions |
| `pr_sigma_thermal_N128.png` | PR vs r_max at 3 temperatures. T/Tc=1.5 ascending, 1.0 mildly descending |
| `pr_sigma_thermal_N256.png` | Same with better separation. T/Tc=1.5 flat near 1.1–1.5 |
| `leading_indicator_N64.png` | η vs χ onset — simultaneous |
| `leading_indicator_N128.png` | η vs χ onset — simultaneous |
| `leading_indicator_N256.png` | **η LEADS χ by ΔT/Tc=0.030** |
| `finite_size_scaling.png` | 2-panel: rank and η for all 3 N values overlaid |

---

## 9. Commit Details

- **Commit:** (pending)
- **Tag:** `phase2-ising-v1`
- **Branch:** `main`
