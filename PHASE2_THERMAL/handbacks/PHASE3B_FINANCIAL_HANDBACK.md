# Phase 3B Handback: Financial Correlation Network Fisher Diagnostics

## 1. THE HEADLINE — Kill Test Result

**Outcome: P3B-1 FAILS — SV2/SV1 does NOT lead market crises.**

The Fisher SV2/SV1 diagnostic on financial correlation networks is a **concurrent/lagging indicator** of market stress, not a leading one. SV2/SV1 rises DURING and AFTER crises as correlations unify, but does not exceed the 2σ threshold before either 2008 or 2020.

| Crisis | SV2/SV1 30d before | SV2/SV1 at crisis | SV2/SV1 peak (post) | Days to 2σ |
|--------|-------------------|------------------|--------------------|----|
| **2008 Lehman** (Sep 15) | 0.759 | 0.797 | 0.890 (Nov 3) | +28d |
| **2020 COVID** (Mar 11) | 0.670 | 0.863 | 0.929 (Mar 23) | +5d |
| **2018 Volmageddon** (Feb 5) | 0.569 | — | 0.817 (Apr 2) | +7d |

**However, the result is scientifically rich.** The SV2/SV1 metric maps directly to the physics SV2/SV1 gradient, and the financial data occupies the expected position between first-order and continuous transitions.

---

## 2. Key Quantitative Results

### 2.1 SV2/SV1 (90d) Statistics

| Metric | Value |
|--------|-------|
| Mean | 0.767 |
| Std | 0.076 |
| Min | 0.536 |
| Max | 0.931 |
| Median | 0.772 |
| VIX correlation | r = 0.591 |

### 2.2 Multi-Scale Comparison

| Window | Mean | Std | Min | Max |
|--------|------|-----|-----|-----|
| 30d | 0.722 | 0.095 | 0.398 | 0.945 |
| 60d | 0.749 | 0.084 | 0.493 | 0.940 |
| 90d | 0.767 | 0.076 | 0.536 | 0.931 |
| 120d | 0.782 | 0.068 | 0.544 | 0.932 |

Longer windows smooth out noise and raise the mean — consistent with longer correlation time-averaging dampening fluctuations.

---

## 3. The Physics Mapping — Where Finance Sits

The financial SV2/SV1 range maps directly onto the physics transition-order gradient:

| System | SV2/SV1 | Interpretation |
|--------|---------|----------------|
| Potts q=10 (strong 1st-order) | 0.374 | Deep disordered |
| **Financial calm (P5)** | **0.623** | — |
| Potts q=5 (weak 1st-order) | 0.691 | — |
| **Financial mean** | **0.767** | — |
| **Financial elevated (P95)** | **0.895** | — |
| **Financial peak crisis** | **0.931** | — |
| 3D Ising (continuous) | 0.956 | Full criticality |
| 2D Ising (continuous) | 1.000 | Full criticality |

**The financial market sits in the Potts q=5 to 3D Ising range** — between weak first-order and continuous transitions. During crises, it approaches (but never reaches) full criticality. This is consistent with the view that:

- **Calm markets** (SV2/SV1 ≈ 0.62): correlation structure is heterogeneous, dominated by one principal factor (market beta), resembling a first-order-like state
- **Crisis markets** (SV2/SV1 ≈ 0.93): all correlations unify ("everything goes to 1"), approaching the critical state where the FIM degeneracies swap
- **The market never fully transitions** (SV2/SV1 < 1.0): interventions (Fed, circuit breakers) or the finite size of the market prevent the full divergence

---

## 4. Crisis-Proximate Behavior

### 4.1 — 2008 Global Financial Crisis

The 2008 crisis developed slowly over 18 months:

| Date | SV2/SV1 | VIX | SPY | Phase |
|------|---------|-----|-----|-------|
| 2007-06-04 | 0.758 | 12.8 | 108.6 | Pre-crisis baseline |
| 2007-10-01 | 0.820 | 18.0 | 108.5 | Subprime concerns emerge |
| 2008-03-03 | 0.824 | 26.5 | 95.4 | Bear Stearns |
| 2008-06-09 | 0.811 | 23.6 | — | Summer calm |
| 2008-08-25 | 0.724 | — | — | **False quiet** |
| 2008-09-15 | 0.797 | 25.7 | 90.9 | Lehman day |
| 2008-11-03 | 0.890 | 59.9 | 70.6 | **Peak crisis** |
| 2009-03-02 | 0.896 | 46.4 | 54.2 | Bottom |

**Key finding**: SV2/SV1 was elevated (0.82) throughout early 2008, but DROPPED to 0.724 in August before Lehman — a **false quiet** that preceded the acute crash. The metric then spiked to 0.890 AFTER Lehman as all correlations unified.

### 4.2 — 2020 COVID Crash

COVID was a sudden shock:

| Date | SV2/SV1 | VIX | Phase |
|------|---------|-----|-------|
| 2020-01-13 | 0.671 | 12.6 | Pre-COVID calm |
| 2020-02-10 | 0.670 | — | Still calm |
| 2020-02-24 | 0.621 | — | **Lowest point** |
| 2020-03-02 | 0.781 | — | First jump |
| 2020-03-09 | 0.863 | — | Z=1.76σ (2 days before WHO pandemic) |
| 2020-03-16 | 0.923 | 57.8 | Peak fear |
| 2020-03-23 | 0.929 | — | **Peak SV2/SV1** |

**For COVID, there's a near-miss**: SV2/SV1 reached 1.76σ on March 9, 2 days before the WHO declared a pandemic — tantalizingly close to the 2σ threshold but not meeting it.

### 4.3 — 2018 Volmageddon

| Date | SV2/SV1 | Context |
|------|---------|---------|
| 2018-01-29 | 0.569 | Pre-event extreme calm |
| 2018-02-12 | 0.762 | Post-Volmageddon |
| 2018-04-02 | 0.817 | Peak |

Volmageddon followed an extreme calm period (SV2/SV1 = 0.569, near-minimum). P3B-4 confirms that the 2008 peak (0.890) > 2018 peak (0.817), distinguishing systemic from technical crises.

---

## 5. Prediction Outcomes

| ID | Prediction | Result | Notes |
|---|---|---|---|
| **P3B-1** | **SV2/SV1 rises before crashes** | **FAIL (KILL)** | Concurrent, not leading. 2σ exceeded +28d (2008), +5d (2020) |
| P3B-2 | Rank transitions before crashes | FAIL | Rank ≈ 1.0 throughout. Financial graph lacks lattice symmetry |
| P3B-3 | η feature near crashes | PASS | η (= SV2/SV1 here) has local extrema near all three crises |
| P3B-4 | 2008 peak > 2018 peak | PASS | 0.890 > 0.817. Systemic > technical |
| P3B-5 | Short-window leads long-window | PASS | 30d exceeds 1σ 14 days before 90d (2008), 7 days (2020) |

### Kill Test Assessment

P3B-1 FAILS under the strict pre-registered criterion (2σ within 30 days before crash). The SV2/SV1 metric is a **concurrent indicator** — it tracks the correlation regime in real-time but does not anticipate it. This makes physical sense: the Fisher kernel measures the CURRENT correlation structure, and the correlation structure changes DURING the crisis, not before.

However, two mitigating observations:
1. **2020 near-miss**: 1.76σ reached 2 days before the WHO pandemic declaration
2. **2008 slow build**: SV2/SV1 was elevated (0.82) for months before Lehman, just not meeting the 2σ threshold because the rolling baseline had already adjusted upward

---

## 6. Signal Characterization

### 6.1 Regime Classification

| Regime | SV2/SV1 Range | Weeks | Fraction | Interpretation |
|--------|--------------|-------|----------|----------------|
| Calm | < 0.728 (P25) | 235 | 25% | Heterogeneous correlations |
| Normal | 0.728 – 0.812 | 471 | 50% | Typical market |
| Elevated | > 0.812 (P75) | 235 | 25% | Stress building |
| Crisis | > 0.861 (P90) | 94 | 10% | Correlations near-unified |

### 6.2 Short-Window Lead Effect

P3B-5 confirms that the 30d window responds faster than the 90d window:

| Crisis | 30d first >1σ | 90d first >1σ | Lead |
|--------|--------------|--------------|------|
| 2008 | Sep 22 | Oct 6 | 14 days |
| 2020 | Mar 2 | Mar 9 | 7 days |

This mirrors the Phase 2/3D physics result: the Fisher diagnostic η responds before susceptibility χ. In the financial context, shorter windows capture the onset of correlation change before longer windows fully absorb it.

---

## 7. Why the Kill Test Fails — Physical Interpretation

The failure has a clear physical explanation:

1. **The Fisher kernel measures the CURRENT correlation function C(r)**, not a predictive one. C(r) changes because market participants become more correlated during stress — this is the crisis itself, not a precursor.

2. **In the Ising model, ξ diverges AT T_c, not before.** The SV degeneracy swap occurs AT the critical point. Similarly, financial correlations unify AT the crisis, not before. The Fisher diagnostic faithfully tracks the underlying structure, but the structure changes concurrently with the crisis.

3. **A leading signal would require a precursor in the correlation structure** — some signature that the network is "approaching criticality" before it arrives. The SV2/SV1 metric doesn't capture this because it measures the instantaneous correlation regime, not the dynamics of approach.

4. **The rank stays at 1** throughout: unlike physics lattices where rank transitions from 1→d+1 at T_c (providing a sharp transition detector), the financial k-NN graph has no lattice symmetry to break. The FIM has one dominant direction (market factor) at all times.

---

## 8. Comparison Across All Phases

| System | Phase | SV2/SV1 Range | Swap? | Rank Pattern | Leading? |
|--------|-------|--------------|-------|-------------|---------|
| 2D Ising | Phase 2 | 1.000 at T_c | YES | 3 always | Fisher leads χ |
| 3D Ising | Phase 3D-1 | 0.956 at T_c | YES | 1→4 at T_c | Fisher leads χ |
| Potts q=5 | Phase 3D-2 | 0.691 at T_c | NO | 1→3 at T_c | Fisher leads χ |
| Potts q=10 | Phase 3D-2 | 0.374 at T_c | NO | 1 always | Fisher leads χ |
| **Financial** | **Phase 3B** | **0.54–0.93** | **NO** | **1 always** | **Concurrent** |

The financial market behaves like a **Potts-like system with variable q**: calm periods resemble q=5–10 (first-order), crisis periods approach continuous (q→4), but the transition never completes.

---

## 9. What DOES Work

Despite the kill test failure, Phase 3B demonstrates:

1. **The SV2/SV1 metric IS responsive to market stress** (r=0.59 with VIX)
2. **It discriminates crisis severity**: 2008 (0.890) > 2018 (0.817) > baseline (0.767)
3. **Short windows lead long windows** (P3B-5), consistent with physics
4. **The physics SV2/SV1 gradient maps meaningfully to finance**: markets sit between first-order and continuous, moving toward criticality during crises
5. **The metric provides regime classification**: calm (25%), normal (50%), elevated (25%), crisis (10%)

---

## 10. Recommendation

**The Fisher SV2/SV1 diagnostic successfully transfers from physics to finance as a regime classifier, but NOT as a leading indicator.** The kill test result is definitive: the metric tracks, not predicts, crises.

**Next steps (if desired):**

1. **Test rate-of-change**: dSV2/SV1/dt might lead the level — the acceleration of correlation change could precede the crisis peak
2. **Test at shorter horizons**: 10d or 5d windows might capture intra-week correlation shifts that precede weekly-scale stress
3. **Alternative graph construction**: minimum spanning tree or planar maximally filtered graph might create more informative distance structure than k-NN
4. **Paper scope**: The physics results (Phases 2, 3D-1, 3D-2) stand on their own. Phase 3B can be an appendix showing the finance application as a regime classifier — which is itself a novel result.

---

## 11. Runtime

| Component | Time |
|-----------|------|
| LEAN backtest (QuantConnect cloud) | 618.8s (~10.3 min) |
| Post-processing | ~5s |
| **Total** | **~10.5 min** |

Data: 85–89 assets, 941 weekly computation dates, 2005-01-03 to 2024-12-30.
