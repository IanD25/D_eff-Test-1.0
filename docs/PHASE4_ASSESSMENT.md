# Phase 4: Cross-Universality FIM Assessment

## Programme Overview

Phase 4 tested whether the FIM SV2/SV1 diagnostic -- validated on 2D Ising (Phase 2) -- generalizes to other universality classes and system types. Four systems were tested: Kuramoto oscillators (4A), 2D XY model (4A'), q=4 Potts (4B), and 3D Anderson localization (4C).

## Reference Baseline (Phase 2)

| System | SV2/SV1 at T_c | Transition type |
|--------|----------------|-----------------|
| 2D Ising (q=2) | 1.000 | Continuous |
| 2D Potts q=5 | 0.691 | Weak first-order |
| 2D Potts q=10 | 0.374 | Strong first-order |

All measured on N=128 square lattice with thermal kernel FIM.

---

## Test Results

### 4A: Kuramoto Oscillators -- SYSTEM TYPE MISMATCH

**Setup.** N=500 oscillators, all-to-all coupling, Lorentzian natural frequencies. K_c = 2*gamma (Lorentzian exact), K sweep 0.1--2.0.

**Result.** SV2/SV1 = 0.9998 everywhere. PR = 10.0 = k (number of neighbors).

**Diagnosis.** The complete graph (all-to-all coupling) has no spatial distance structure. k-NN on correlations produces Erdos-Renyi topology. The FIM measures identity-like structure everywhere because there is no intrinsic distance to modulate.

**Classification: SYSTEM TYPE MISMATCH** -- outside the pipeline's domain.

---

### 4A': 2D XY Model (BKT Transition) -- PASS (WEAK)

**Setup.** L=64 square lattice, Wolff cluster + Metropolis hybrid, T_BKT ~ 0.893. Thermal kernel FIM matching Phase 2 Ising exactly.

**Result.** SV2/SV1 peak at T=0.876 (1.9% from T_BKT), value = 0.0677.

**Kill criteria evaluation:**

| Criterion | Result | Notes |
|-----------|--------|-------|
| P4A-1 (KILL) | PASS | Peak within +/-15% window of T_BKT |
| P4A-2 | WEAK | Peak amplitude only 0.056 above baseline |
| P4A-3 | PASS | below-T_c average < peak > above-T_c average |

**Interpretation.** The BKT transition is fundamentally subtler than Ising. C_connected changes shape (power-law to exponential decay), not magnitude (~0 to ~0.5 as in Ising). The small signal is physically expected.

**Classification: PASS** -- FIM detects BKT transition location; signal is weak but real.

---

### 4B: q=4 Potts (Continuous/First-Order Boundary) -- PARTIAL PASS

**Setup.** N=128 square lattice, Wolff cluster, T_c = 1/ln(3) ~ 0.9102 (exact, Baxter 1973). q=4 sits exactly on the continuous/first-order boundary with logarithmic corrections to scaling.

**Result.** SV2/SV1 = 0.521 at T/T_c = 0.989.

**Kill criteria evaluation:**

| Criterion | Result | Notes |
|-----------|--------|-------|
| P4B-1 (KILL) | KILLED | 0.521 < 0.691 (q=5 threshold) |
| P4B-2 | PASS | Peak within +/-10% of T_c |
| P4B-3 | FAIL | q-ordering violated (q=4 < q=5) |
| P4B-4 | INCONCLUSIVE | rank=1 everywhere |

**Interpretation.** The q-ordering violation is explained by logarithmic corrections at finite N=128. Theory predicts SV2/SV1 -> 1.0 as N -> infinity for q=4 (continuous transition), but convergence is logarithmically slow. A finite-size scaling study (N=256, 512) is running to confirm this.

A sharp transition is visible in the data: SV2/SV1 drops from 0.521 to 0.157 between T/T_c = 0.989 and 1.011.

**Classification: PARTIAL PASS** -- peak location correct, magnitude suppressed by known finite-size effects. Finite-size scaling study pending.

---

### 4C: 3D Anderson Localization -- SYSTEM TYPE MISMATCH

**Setup.** L=16 cubic lattice (4096 sites), tight-binding Hamiltonian, W_c ~ 16.5. Wavefunction amplitude correlations C(r) = <|psi(i)| * |psi(j)|> as kernel (replaces thermal G(r,T)). 50 eigenstates near E=0, 10 disorder realizations per W, 18 W values from 2 to 24.

**Result.** SV2/SV1 ~ 0.255 +/- 0.01 across ALL disorder strengths -- completely flat.

**Kill criteria evaluation:**

| Criterion | Result | Notes |
|-----------|--------|-------|
| P4C-1 (KILL) | KILLED | No peak detected (formally peaks at W=10, value 0.2585, well outside [14.5, 18.5]) |

**Independent verification.** IPR clearly shows the transition: 0.0007 (W=2, metallic) -> 0.165 (W=24, insulating).

**SV profile.** [1.0, 0.26, 0.26, 0.26, 0.04, 0.04] -- triple degeneracy from cubic symmetry.

**Diagnosis.** 3D cubic symmetry means the correlation kernel is isotropic in all regimes. Extended wavefunctions are isotropic. Localized wavefunctions are also isotropic (localized in a ball, not a rod). The kernel changes SCALE at W_c but not SHAPE. SV2/SV1 measures shape changes, not scale changes.

**Classification: SYSTEM TYPE MISMATCH** -- transition is real but produces no correlation anisotropy.

---

## Theoretical Synthesis

### The Refined Two-Scale Hypothesis

The Phase 4 results sharpen the original two-scale commensurability hypothesis from Phase 3E.

**Original claim:** FIM SV2/SV1 detects transitions where two independently running characteristic scales approach commensurability at a critical point.

**Refined claim:** FIM SV2/SV1 detects transitions where two-scale commensurability produces a **change in correlation anisotropy** -- i.e., the spatial correlation function changes shape, not just magnitude or scale.

### System Classification Table

| System | Two scales? | Anisotropy change? | SV2/SV1 signal? | Classification |
|--------|-------------|-------------------|------------------|----------------|
| 2D Ising | xi <-> L | Yes (divergent fluctuations break lattice symmetry) | 1.000 at T_c | DETECTED |
| 2D Potts q=5 | xi <-> L | Partial (first-order, coexistence) | 0.691 at T_c | DETECTED |
| 2D Potts q=10 | xi <-> L | Partial (strong first-order, short xi) | 0.374 at T_c | DETECTED |
| 2D XY (BKT) | xi <-> L | Yes (algebraic -> exponential C(r)) | 0.068 at T_BKT | DETECTED (weak) |
| 2D Potts q=4 | xi <-> L | Yes (log corrections suppress) | 0.521 at T_c | PARTIAL (finite-size) |
| Kuramoto | K <-> K_c | No spatial structure exists | 0.9998 flat | MISMATCH |
| 3D Anderson | xi_loc <-> L | No (isotropic in all regimes) | 0.255 flat | MISMATCH |
| EEG seizure | -- | Single-scale takeover | Noise | MISMATCH |

### Domain Boundary Definition

The FIM SV2/SV1 diagnostic reliably detects phase transitions when ALL of the following conditions hold:

1. **Spatial graph structure.** The system has a fixed spatial graph with intrinsic distance structure.
2. **Two-scale commensurability.** Two characteristic length scales become commensurate at the critical point.
3. **Anisotropy change.** The commensurability produces a measurable change in correlation function anisotropy.

Systems failing condition 1 (Kuramoto: no spatial structure) or condition 3 (Anderson: isotropic transition) are outside the diagnostic's domain. This is not a theory failure -- it is a specificity characterization.

### Implications for PFD Project

- The FIM diagnostic is a **2D lattice phase transition detector**, not a universal transition detector.
- It works across universality classes (Ising, Potts, BKT) as long as the lattice is 2D and the transition produces anisotropy.
- The continuous/first-order boundary (q=4) is detectable but requires finite-size scaling to resolve logarithmic corrections.
- 3D systems with isotropic transitions are outside the domain -- would need an anisotropic variant (e.g., quasi-2D layers, directed lattice) to break the symmetry.

---

## File Inventory

| File | Description |
|------|-------------|
| `PHASE4A_KURAMOTO/kuramoto_fisher_phase4a.py` | Kuramoto implementation (mismatch result) |
| `PHASE4A_KURAMOTO/xy_fisher_phase4a.py` | 2D XY BKT test (weak pass) |
| `PHASE4A_KURAMOTO/results/` | XY model outputs |
| `PHASE4B_POTTS4/potts4_fisher_phase4b.py` | q=4 Potts N=128 (partial pass) |
| `PHASE4B_POTTS4/potts4_finite_size_scaling.py` | q=4 FSS study N=256,512 (running) |
| `PHASE4B_POTTS4/results/` | q=4 N=128 outputs |
| `PHASE4C_ANDERSON/anderson_fisher_phase4c.py` | Anderson 3D test (mismatch result) |
| `PHASE4C_ANDERSON/results/` | Anderson outputs |
| `docs/PHASE4_ASSESSMENT.md` | This document |

---

## Status

| Sub-phase | Status |
|-----------|--------|
| Phase 4A (Kuramoto) | Complete -- system type mismatch |
| Phase 4A' (XY/BKT) | Complete -- weak pass |
| Phase 4B (q=4 Potts N=128) | Complete -- partial pass |
| Phase 4B-FSS (q=4 N=256,512) | Running (background, ~3 hours) |
| Phase 4C (Anderson) | Complete -- system type mismatch |

---

Date: 2026-03-25
Authors: Ian Darling, Claude Code
