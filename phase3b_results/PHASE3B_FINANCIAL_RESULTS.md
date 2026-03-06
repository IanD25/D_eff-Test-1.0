# Phase 3B Results: Financial Correlation Network Fisher Diagnostics

Generated: 2026-03-06 01:40:37

## 1. Data Summary

- Date range: 2005-01-03 to 2024-12-30
- Total computation dates: 941
- Valid 90d results: 941 (100.0%)
- Windows: [30, 60, 90, 120]
- k_nn: 10
- Source: phase3b_results/raw_data/fisher_financial_results.json

## 2. SV2/SV1 (90d) Summary Statistics

| Metric | Value |
|--------|-------|
| Mean | 0.7671 |
| Std | 0.0759 |
| Min | 0.5364 |
| Max | 0.9311 |
| Median | 0.7720 |

## 3. Rank Summary (90d)

| Metric | Value |
|--------|-------|
| Mean | 1.02 |
| Std | 0.11 |
| Min | 1.0 |
| Max | 2.6 |

## 4. Prediction Assessment

| ID | Prediction | Result | Notes |
|---|---|---|---|
| **P3B-1** | **SV2/SV1 rises before crashes** | **FAIL (KILL)** | Kill test |
| P3B-2 | Rank transitions before crashes | FAIL | |
| P3B-3 | eta feature near crashes | PASS | |
| P3B-4 | 2008 peak > 2018 peak | PASS | |
| P3B-5 | Short-window leads long-window | PASS | |

### P3B-1 Detail (Kill Test)

**P3B-1_2020_covid_60d:** pass=True, date=2020-03-09, zscore=1.7605444091006681, sv2_sv1=0.8629166, days_before=2, 

### P3B-4 Detail
- 2008 peak SV2/SV1: 0.8898372
- 2018 peak SV2/SV1: 0.8165295

## 5. Crisis-Proximate Behavior

### Near 2008 Lehman (Sep 15, 2008)
| 2008-07-21 | SV2/SV1=0.765 | rank=1.0 | eta=0.765 |
| 2008-07-28 | SV2/SV1=0.752 | rank=1.0 | eta=0.752 |
| 2008-08-04 | SV2/SV1=0.754 | rank=1.0 | eta=0.754 |
| 2008-08-11 | SV2/SV1=0.755 | rank=1.0 | eta=0.755 |
| 2008-08-18 | SV2/SV1=0.759 | rank=1.0 | eta=0.759 |
| 2008-08-25 | SV2/SV1=0.724 | rank=1.0 | eta=0.724 |
| 2008-09-08 | SV2/SV1=0.770 | rank=1.0 | eta=0.770 |
| 2008-09-15 | SV2/SV1=0.797 | rank=1.0 | eta=0.797 |
| 2008-09-22 | SV2/SV1=0.801 | rank=1.0 | eta=0.801 |
| 2008-09-29 | SV2/SV1=0.812 | rank=1.0 | eta=0.812 |
| 2008-10-06 | SV2/SV1=0.841 | rank=1.0 | eta=0.841 |
| 2008-10-13 | SV2/SV1=0.860 | rank=1.0 | eta=0.860 |

### Near COVID (Mar 11, 2020)
| 2020-01-13 | SV2/SV1=0.671 | rank=1.0 | eta=0.671 |
| 2020-01-27 | SV2/SV1=0.697 | rank=1.0 | eta=0.697 |
| 2020-02-03 | SV2/SV1=0.686 | rank=1.0 | eta=0.686 |
| 2020-02-10 | SV2/SV1=0.670 | rank=1.0 | eta=0.670 |
| 2020-02-24 | SV2/SV1=0.621 | rank=1.0 | eta=0.621 |
| 2020-03-02 | SV2/SV1=0.781 | rank=1.0 | eta=0.781 |
| 2020-03-09 | SV2/SV1=0.863 | rank=1.0 | eta=0.863 |
| 2020-03-16 | SV2/SV1=0.923 | rank=1.0 | eta=0.923 |
| 2020-03-23 | SV2/SV1=0.929 | rank=1.0 | eta=0.929 |
| 2020-03-30 | SV2/SV1=0.923 | rank=1.0 | eta=0.923 |
| 2020-04-06 | SV2/SV1=0.927 | rank=1.0 | eta=0.927 |

## 6. Figures

| # | File | Description |
|---|------|-------------|
| 1 | sv2_sv1_vs_vix.png | HEADLINE: SV2/SV1 vs VIX with crisis markers |
| 2 | diagnostic_panel_timeseries.png | 4-panel diagnostic time series |
| 3 | multiscale_sv2_sv1.png | Multi-scale SV2/SV1 |
| 4 | pr_window_profile.png | PR vs window at key dates |
| 5 | crisis_zoom_2008.png | 2008 crisis zoom |
| 6 | crisis_zoom_2020.png | COVID crisis zoom |
| 7 | crisis_zoom_2018.png | 2018 Volmageddon zoom |
| 8 | sv_profile_snapshots.png | SV bar charts at key dates |
| 9 | correlation_function_key_dates.png | C(r) at key dates |

## 7. Comparison to Physics Results

| System | Context | SV2/SV1 at transition | Swap? |
|--------|---------|----------------------|-------|
| 2D Ising | T_c (continuous) | 1.000 | YES |
| 3D Ising | T_c (continuous) | 0.956 | YES |
| Potts q=5 | T_c (1st-order, weak) | 0.691 | NO |
| Potts q=10 | T_c (1st-order, strong) | 0.374 | NO |
| Financial | Crisis (if applicable) | See above | TBD |

