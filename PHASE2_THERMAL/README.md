# DS Framework — Phase 2: Spectral Phase Transitions in the Fisher Information Matrix

**Status**: ✅ Complete | **Paper**: [DS_Phase2_Thermal_Paper_FINAL_v2.pdf](DS_Phase2_Thermal_Paper_FINAL_v2.pdf) | **Author**: Ian Darling

---

## Overview

Phase 2 replaces the exponential kernel from Phase 1 with the **measured spin-spin correlation function** of physical systems, turning the FIM into a structural probe of phase transitions. The same SVD pipeline — rank, η, PR, SV profile — now reads the **information geometry of correlation structures** rather than spatial geometry.

The key discovery: the FIM singular value spectrum undergoes a **degeneracy swap** at continuous phase transitions that is absent at first-order transitions, enabling parameter-free classification of transition order.

---

## The Three Information-Geometric States

The SV degeneracy pattern defines three states. These are properties of the correlation structure, not the physical mechanism:

| State | SV Profile | Correlation | Physical Analog |
|-------|-----------|-------------|-----------------|
| **1 — Radial-dominated** | [1, a, a, small], a≪1 | Exponential decay, ξ≪L | Paramagnetic / normal regime |
| **2 — Isotropic** | [1, ≈1, ≈1, small] | Power-law / ξ~L | Critical point / crisis onset |
| **3 — Noise-dominated** | [1, ~0.5, ~0.5, ~0.2] | Flat / saturated | Ordered phase / post-crisis lockstep |

**Critical result:** Below T_c, ξ/L → ∞ but SV₂/SV₁ drops. The SV profile is **not** a proxy for ξ/L.

---

## Systems Tested

### Phase 2: 2D Ising Model
T_c = 2/ln(1+√2) ≈ 2.2692. Lattice sizes N = 64, 128, 256. Wolff cluster MC. Manhattan distance.

| T/T_c | SV1 | SV2 | SV3 | SV4 | State |
|-------|-----|-----|-----|-----|-------|
| 1.500 | 1.00 | 0.23 | 0.23 | 0.04 | State 1 |
| 1.020 | 1.00 | 1.00 | 0.34 | 0.09 | State 2 |
| 1.000 | 1.00 | 1.00 | 0.54 | 0.08 | State 2 |
| 0.900 | 1.00 | 0.48 | 0.48 | 0.22 | State 3 |

**Rank = 3.00 ± 0.00 at all temperatures.** η minimum at T_c (0.147). Leading indicator: SV swap precedes χ peak by ΔT/T_c ≈ 0.10–0.13.

### Phase 3D-1: 3D Ising (ν ≈ 0.630)
Degeneracy swap **confirmed**. SV₂/SV₁ = 0.956 at T_c. Rank transition 1→4 sharpens with L.

### Phase 3D-2: Potts q=5 and q=10 (First-Order)
**No swap** at either first-order transition.

| System | Transition | SV₂/SV₁ at T_c | Swap? | Rank Pattern |
|--------|-----------|----------------|-------|-------------|
| 2D Ising | Continuous | 1.000 | YES | 3 always |
| 3D Ising | Continuous | 0.956 | YES | 1→4 at T_c |
| Potts q=5 | 1st-order (weak) | 0.691 | NO | 1→3 at T_c |
| Potts q=10 | 1st-order (strong) | 0.374 | NO | 1 always |

**Classifier:** SV₂/SV₁ > 0.95 → continuous. SV₂/SV₁ < 0.70 → first-order.

---

## Cross-Domain Transfer Results

The identical FIM pipeline (zero domain-specific tuning) was applied to:

**Financial correlations** (~90 S&P 500 stocks, 2005–2024, weekly):
- η extremum at all three crises: **6/6** sub-tests
- Rank→1 during crises: **5/6** sub-tests
- SV₂/SV₁ peaks at 0.65–0.72 (weakly first-order analog)
- VIX correlation r ≈ 0.15 (not a volatility proxy)

**EEG seizure data** (CHB-MIT, 5 patients, 26 seizures):
- Rank→1 during ictal: **5/5**
- η extremum: **4/5**
- SV₂/SV₁ flat (graph too small, diameter ~2–3)

| Diagnostic | Lattice | Financial | EEG |
|-----------|---------|-----------|-----|
| Rank→1 at transition | ✓ | 5/6 | 5/5 |
| η extremum | ✓ | 6/6 | 4/5 |
| SV₂/SV₁ elevation | ✓ | ✓ | ✗ |
| Full degeneracy swap | ✓ | ✗ | ✗ |

**Scope condition:** SV₂/SV₁ amplitude requires graph diameter ≥ 5. Rank and η work at any scale.

---

## Thermal Kernel Construction

Replace the exponential kernel with measured correlations:

$$p_{v_0}(u; T) \propto |G(|v_0 - u|, T)|$$

where $G(r,T) = \langle s_0 s_r \rangle - \langle s_0 \rangle \langle s_r \rangle$.

No free parameters beyond T. Rest of pipeline verbatim from Phase 1: score vectors → FIM → SVD → rank / η / PR / SV profile.

---

## Prediction Outcomes

### Phase 2 (2D Ising): 2 PASS, 2 PARTIAL, 3 FAIL
### Phase 3D (3D Ising + Potts): 4 PASS, 0 PARTIAL, 4 FAIL
### Phase 3B (Financial): 2 PASS, 0 PARTIAL, 2 FAIL, 1 MIXED

**Total across thermal program: 8 PASS, 2 PARTIAL, 9 FAIL.**

Failures cluster in two categories: (a) applying geometric-regime expectations to thermal data, (b) predicting wrong direction of rank change. Two failures (P3B-2 rank direction, P2-2 η direction) produced findings more interesting than the original predictions.

---

## Known Vulnerabilities

1. **q=5 Potts at L=128 ≪ ξ≈2500.** Deep in pseudocritical regime. q=10 result (L/ξ≈25) is clean.
2. **Leading indicator may be FSS artifact.** Conservative estimates: ΔT/T_c ≈ 0.04–0.10.
3. **Only 4 lattice systems tested.** BKT, tricritical, deconfined QCP untested.
4. **SV₂/SV₁ gradient rests on 4 data points.** Functional form not established.

---

## Reproducing Results

### Quick Validation (~10 min)
```bash
cd PHASE2_THERMAL/src/
python ising_fisher_phase2.py          # 2D Ising sweep
python ising_3d_fisher.py              # 3D Ising (Phase 3D-1)
python potts_fisher.py                 # Potts q=5, q=10 (Phase 3D-2)
```

Each script is self-contained. Outputs: CSV results + PNG plots to `results/`.

### Financial Transfer (~6 min, requires QuantConnect LEAN)
```bash
# Deploy fisher_diagnostic_3b_v3_lean.py to QuantConnect
# Download logs.txt after backtest completes
# Results are logged inline — see handbacks/DS_PHASE3B_HANDBACK.md for grep patterns
```

### EEG Transfer (~15 min, downloads data automatically)
```bash
pip install mne wfdb numpy scipy matplotlib
python eeg_fisher_phase3e.py
```

---

## Files

```
PHASE2_THERMAL/
├── README.md                              (this file)
├── DS_Phase2_Thermal_Paper_FINAL_v2.pdf   (paper)
├── DS_Phase2_Thermal_Paper_FINAL_v2.tex   (source)
├── DS_Phase2_Thermal_Paper_FINAL.pdf      (earlier draft)
├── src/
│   ├── ising_fisher_phase2.py             (2D Ising — Phase 2)
│   ├── ising_3d_fisher.py                 (3D Ising — Phase 3D-1)
│   ├── potts_fisher.py                    (Potts q=5,10 — Phase 3D-2)
│   ├── fisher_diagnostic_3b_v3_lean.py    (Financial — Phase 3B, LEAN)
│   └── eeg_fisher_phase3e.py              (EEG seizure — Phase 3E)
├── results/
│   ├── phase2_ising_2d/
│   ├── phase3d1_ising_3d/
│   ├── phase3d2_potts/
│   ├── phase3b_financial/
│   └── phase3e_eeg/
├── figs/                                  (paper figures)
├── specs/
│   ├── DS_PHASE2_ISING_FISHER_SPEC.md
│   ├── DS_PHASE3D1_3D_ISING_SPEC.md
│   ├── DS_PHASE3D2_POTTS5_SPEC.md
│   ├── DS_PHASE3B_SPEC_v3_SIMPLIFIED.md
│   └── DS_PHASE3E_EEG_SPEC.md
└── handbacks/
    ├── PHASE2_ISING_FISHER_HANDBACK.md
    ├── PHASE3D1_3D_ISING_HANDBACK.md
    ├── PHASE3D2_POTTS5_HANDBACK.md
    └── DS_PHASE3B_HANDBACK.md
```

---

## Key Findings

1. **Degeneracy swap is universal for continuous transitions** — confirmed in 2D Ising (ν=1) and 3D Ising (ν≈0.630)
2. **Absent at first-order transitions** — confirmed at two strengths (q=5 weak, q=10 strong)
3. **SV₂/SV₁ classifies transition order** — clean gap between continuous (>0.95) and first-order (<0.70)
4. **Three information-geometric states generalize beyond lattice systems** — rank→1 and η extremum transfer to financial networks (6/6) and neural recordings (5/5)
5. **SV profile ≠ ξ/L proxy** — below T_c they diverge qualitatively
6. **Financial markets behave like weakly first-order systems** — SV₂/SV₁ ≈ 0.65–0.72 at crisis peaks

---

## References

**Paper 2**: Darling, Ian. *Spectral Phase Transitions in the Fisher Information Matrix: Universality of the Degeneracy Swap and Classification of Transition Order.* Preprint (2026).

**Paper 1**: Darling, Ian. *Effective Dimensionality from Fisher Information Rank.* Preprint (2026). See [PHASE1_GEOMETRIC/](../PHASE1_GEOMETRIC/)

---

## AI Disclosure

This research was conducted with AI-assisted development. Claude Opus (Anthropic) contributed to methodology design, analysis, and manuscript development. Claude Code (Anthropic) contributed to implementation. The author is responsible for all scientific claims, research direction, and editorial decisions.
