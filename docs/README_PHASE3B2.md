# DS Phase 3B-2: Sector-Level Fisher Decomposition

## Implementation Complete

### Files Created:

1. **`FisherSectorAlgorithm.py`** - LEAN backtest algorithm
   - Computes market-level Fisher (k=10, 90-day window)
   - Computes 7 sector-level Fisher decompositions (k=5/4 based on sector size)
   - Tracks 30d and 90d windows for spread signal
   - Extracts per-asset Fisher temperatures
   - Saves all results to ObjectStore

2. **Analysis Scripts** - Post-processing and visualization
   - `fisher_sector_analysis_part1.py` - Figures 1-3 + core functions
   - `fisher_sector_analysis_part2.py` - Figures 4-7
   - `fisher_sector_analysis_part3.py` - Figures 8-9 + predictions + report

### Key Features Implemented:

**Dual-Level Computation:**
- Market level: Full ~85 asset universe, k=10, 90-day window (Phase 3B baseline)
- Sector level: 7 GICS sectors independently, k=5 (or k=4 for small sectors)

**Three Testable Signals:**
1. **Sector divergence** - Which sector's SV2/SV1 rises first?
2. **Per-asset Fisher temperature** - Individual asset criticality ranking
3. **Window spread** - 30d-90d spread for active transition detection

**Kill Test (P3B2-1):**
- Financials must exceed 2σ ≥30 days before Market in 2007-2008
- If fails → sector decomposition hypothesis rejected

### How to Use:

1. **Upload to QuantConnect:**
   - Copy `FisherSectorAlgorithm.py` to your QuantConnect project
   - Ensure required packages: numpy, scipy, pandas

2. **Run Backtest:**
   - Period: 2005-01-01 to 2024-12-31
   - Runtime: ~28 minutes (estimated)
   - Results saved to ObjectStore as JSON

3. **Download Results:**
   - Export results from ObjectStore as `fisher_sector_results.json`

4. **Run Analysis:**
   ```bash
   # Combine analysis scripts
   cat fisher_sector_analysis_part1.py fisher_sector_analysis_part2.py fisher_sector_analysis_part3.py > fisher_sector_analysis.py
   
   # Update json_path in main() function
   # Run analysis
   python fisher_sector_analysis.py
   ```

5. **Output:**
   - 9 figures in `phase3b2_results/` directory
   - Results report: `PHASE3B2_SECTOR_RESULTS.md`
   - Prediction test results with PASS/FAIL

### Validation Gates:

**Gate A:** ≥5 of 7 sectors produce valid Fisher diagnostics at ≥90% of dates
**Gate B:** Each sector's SV2/SV1 std > 0.05 and not perfectly correlated with Market
**Gate C:** Data integrity (PSD matrices, aligned dates, no NaN cascades)

### Pre-Registered Predictions:

| ID | Prediction | PASS Criterion |
|----|------------|----------------|
| **P3B2-1** | **Financials leads Market in 2007-2008** | ≥30 day lead (KILL TEST) |
| P3B2-2 | Energy leads during 2014-2016 oil crash | ≥30 day lead |
| P3B2-3 | Tech leads during 2022 rate selloff | ≥30 day lead |
| P3B2-4 | Spread > 0 precedes SV2/SV1 peak | ≥2 of 3 crises |
| P3B2-5 | High temp assets cluster in leading sector | >50% from leading sector |
| P3B2-6 | Cross-sector correlation increases in crises | Crisis >0.8, Calm <0.5 |

### Figures Generated:

1. `sector_sv2_sv1_timeseries.png` - All sectors + market
2. `financials_vs_market_2007_2008.png` - **KILL TEST** visualization
3. `sector_lead_lag_heatmap.png` - Crisis lead/lag heatmap
4. `sector_spread_timeseries.png` - 30d-90d spread per sector
5. `energy_vs_market_2014_2016.png` - Energy during oil crash
6. `tech_vs_market_2022.png` - Tech during rate hike selloff
7. `top_fisher_temperature_assets.png` - Top-10 assets per week
8. `sector_rank_comparison.png` - Sector vs market rank
9. `cross_sector_correlation_of_sv2sv1.png` - Correlation matrix

### Known Risks & Mitigations:

- **Small sector subgraphs** (N=7-15): Flag sectors with diameter <3 as unreliable
- **Noisy sector SV2/SV1**: Use 4-week rolling average for plotting
- **Missing tickers pre-2010**: Filter to active tickers, minimum 7 per sector
- **Survivorship bias**: Use tickers that survived entire period

### Scientific Objective:

Test if **sector-level stress precedes market-wide stress**. If Financials showed elevated SV2/SV1 in early 2007 before the broad market responded in 2008, that's a cross-sectional leading signal — identifying *where* stress is building, not just *when*.

---

**Status:** Ready for QuantConnect deployment
**Target Runtime:** ≤45 minutes on M4 MacBook Air
**Kill Condition:** P3B2-1 failure triggers hypothesis rejection