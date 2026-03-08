# Phase 3B-2: Sector Fisher Decomposition — Results

## 1. Kill Test Result

**P3B2-1: PASS**

| Parameter | Value |
|-----------|-------|
| Financials 2σ crossing | 2006-08-21 |
| Market 2σ crossing | 2008-10-20 |
| Lead time (Fin → Mkt) | 791 days |
| Financials peak z-score | 2.23σ |
| Market peak z-score | 2.90σ |
| Financials z at Lehman | 0.22σ |
| Market z at Lehman | -0.29σ |

**PASS: Financials led Market by ≥ 30 days.**

---

## 2. All Prediction Results

| ID | Prediction | Result | Notes |
|----|------------|--------|-------|
| **P3B2-1** | **Financials leads Market 2007-2008** | **PASS** | Lead = 791d. Fin 2σ: 2006-08-21, Mkt 2σ: 2008-10-20 |
| P3B2-2 | Energy leads during 2014-2016 oil crash | PASS | Lead = 350d. Energy peak |z| = 5.38σ |
| P3B2-3 | Tech leads during 2022 rate selloff | PASS | Lead = 63d. Tech peak |z| = 2.70σ |
| P3B2-4 | Spread turns positive before SV2/SV1 peak | PASS | 3/3 crises show positive spread before peak |
| P3B2-5 | High Fisher-temp assets in leading sector | FAIL | 0/3 crises: hottest sector = expected leading sector |
| P3B2-6 | Cross-sector correlation rises during crises | FAIL | Calm r=0.291, Crisis r=0.781 |

---

## 3. Sector SV2/SV1 Statistics (2005–2024)

| Sector | Mean | Std | Min | Max |
|--------|------|-----|-----|-----|
| Market | 0.7742 | 0.0757 | 0.5500 | 0.9352 |
| Consumer | 0.9249 | 0.0373 | 0.7792 | 0.9846 |
| Energy | 0.9756 | 0.0184 | 0.8280 | 0.9967 |
| Financials | 0.9670 | 0.0148 | 0.8943 | 0.9923 |
| HealthCare | 0.9241 | 0.0407 | 0.7169 | 0.9883 |
| Industrials | 0.9432 | 0.0309 | 0.8201 | 0.9892 |
| Technology | 0.9462 | 0.0258 | 0.8418 | 0.9919 |
| Utilities | 0.9831 | 0.0108 | 0.9382 | 0.9979 |

**Key observation:** All sectors have SV₂/SV₁ mean > 0.90 vs Market mean = 0.767. Sectors are near-critical by default due to small N (7-15 assets) and high intra-sector correlations.

---

## 4. Key Scientific Finding: Inverted Behavior

Sector SV₂/SV₁ behaves OPPOSITE to market-level SV₂/SV₁ during crises:

- **Market**: SV₂/SV₁ RISES during crises (all assets co-move, correlations unify → near-critical)
- **Sectors**: SV₂/SV₁ can RISE before crisis onset (intra-sector stress concentrates) then may DROP during acute crisis (intra-sector differentiation increases as some names are hit harder than others)

In 2008, Financials peaked at ~0.986 in Sep 2007 (14 months before Lehman) as subprime contagion unified the sector, then slightly declined as some banks failed (Bear Stearns, Lehman) while others survived. The sector peak preceded the market peak by ~12+ months.

---

## 5. P3B2-4 Detail — Window Spread

| Crisis | Result | Notes |
|--------|--------|-------|
| 2008 GFC | PASS | Spread positive 2007-03-05, SV peak 2008-07-14, lead=497d |
| 2014 Oil | PASS | Spread positive 2014-01-13, SV peak 2015-11-23, lead=679d |
| 2022 Tech | PASS | Spread positive 2021-06-07, SV peak 2022-12-19, lead=560d |

---

## 6. P3B2-5 Detail — Fisher Temperature by Sector

Note: Per-asset Fisher temperatures were not exported to QC charts (per-asset data stored in ObjectStore JSON, not plotted). Assessment uses sector-average SV₂/SV₁ z-score as proxy.

| Crisis | Expected Leading Sector | Observed Hottest Sector | Result |
|--------|------------------------|------------------------|--------|
| 2008 GFC | Financials | Consumer | FAIL |
| 2014 Oil | Energy | Utilities | FAIL |
| 2022 Tech | Technology | Utilities | FAIL |

---

## 7. P3B2-6: Cross-Sector Correlation

Calm 2012-2014: r=0.291, COVID 2020: r=0.781

Criterion: Crisis > 0.8 AND Calm < 0.8. Result: FAIL.

---

## 8. Runtime

| Component | Time |
|-----------|------|
| LEAN backtest (QuantConnect cloud) | 618.8s (~10.3 min) |
| Post-processing | ~5s |
| **Total** | **~10.5 min** |

Data: 7 sectors + market, 941 weekly dates, 2005-01-03 to 2024-12-30.
