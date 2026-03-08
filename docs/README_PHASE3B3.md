# DS Phase 3B-3: Fisher Regime Attractor — Multi-Scale Gradient System

## Overview

The Fisher Regime Attractor (FRA) is a sophisticated multi-scale diagnostic system that outputs:

1. **FRA(t) ∈ [-1, +1]** — A single signed scalar representing the market's position on the criticality spectrum
2. **Cluster heat vector** — Continuous gradient scores for dynamically-defined asset communities
3. **Asset temperature ranking** — Every asset scored by its local Fisher state relative to its community

All outputs are **continuous gradients**, not binary signals. No hardcoded thresholds.

## Scientific Objective

Build a system that tracks market regime transitions through:
- **Macro level**: Market-wide correlation structure (FRA)
- **Meso level**: Dynamic community stress (cluster heat)
- **Micro level**: Individual asset positioning (asset temperature)

The FRA trajectory reveals regime cycles and crisis onset patterns.

## Implementation Files

### Core Algorithm
- **`FisherRegimeAttractor.py`** — Full implementation with real data support
- **`FisherRegimeAttractor_Fast.py`** — Optimized version for faster computation

### Visualization & Analysis
- **`fra_visualization.py`** — Generates 7 required figures and summary report

## Key Concepts

### Fisher Regime Attractor (FRA)

```
FRA(t) = tanh(-z(t) / 2)

where:
  z(t) = (SV2/SV1(t) - mean_252(t)) / std_252(t)
  
Interpretation:
  -1 = Maximum unification (panic, single-factor)
   0 = Equilibrium
  +1 = Maximum dispersion (stock-picking)
```

### Multi-Window FRA

Compute at 30d, 60d, 90d windows:
- **FRA momentum** = FRA_30d - FRA_90d
- Negative momentum = short-term tightening = crisis approaching
- Positive momentum = short-term loosening = recovery

### Dynamic Clustering

Each week, cluster assets by correlation structure:
1. Convert correlation to distance: D_ij = sqrt(2 * (1 - corr_ij))
2. Build affinity matrix with Gaussian kernel
3. Compute Laplacian eigenvectors
4. Determine k from eigengap (no hardcoding)
5. Run k-means on top-k eigenvectors
6. Merge small clusters (N < 12) into nearest neighbors

### Cluster Heat Vector

Per-cluster Fisher temperature:
```
Cluster_heat(c, t) = z-scored SV2/SV1 within cluster c
```

Positive = unification stress, Negative = dispersion opportunity

### Asset Temperature Ranking

Per-asset composite score:
```
Composite(i, t) = 0.4 * FRA(t) 
                + 0.4 * Cluster_heat(c_i, t)
                + 0.2 * Asset_temp(i, t)

where:
  Asset_temp(i, t) = SV2/SV1(i) - mean(SV2/SV1 in cluster c_i)
```

All assets ranked by |Composite| descending.

## Usage

### With Real Data

```python
from FisherRegimeAttractor import FisherRegimeAttractor

# Initialize
fra = FisherRegimeAttractor(
    data_path='your_price_data.csv',  # CSV: date, ticker, close
    start_date='2006-01-01',
    end_date='2024-12-31'
)

# Load data
price_data = fra.load_data()

# Compute regime
fra.compute_regime(price_data)

# Save results
fra.save_results()

# Generate report
fra.generate_report()
```

### With Synthetic Data (Demo)

```python
from FisherRegimeAttractor_Fast import FisherRegimeAttractor

fra = FisherRegimeAttractor()
price_data = fra.create_synthetic_data()
fra.compute_regime(price_data)
fra.save_results()
fra.generate_report()
```

### Generate Visualizations

```bash
python fra_visualization.py
```

This creates 7 figures in `phase3b3_results/`:
1. FRA trajectory
2. Phase portrait (FRA vs velocity)
3. Cluster evolution
4. Top assets heatmap
5. Cluster heat distribution
6. Composite score distribution
7. FRA component decomposition

## Output Files

### JSON Results
- **`fra_history.json`** — FRA values and derivatives
  ```json
  {
    "date": "2024-01-08",
    "fra_30d": 0.123,
    "fra_60d": 0.045,
    "fra_90d": -0.087
  }
  ```

- **`cluster_history.json`** — Cluster assignments and heats
  ```json
  {
    "date": "2024-01-08",
    "fra": -0.087,
    "n_clusters": 8,
    "cluster_heats": {"0": 0.45, "1": -0.23, ...},
    "top_10_assets": [
      {
        "name": "AAPL",
        "composite": 0.234,
        "macro_fra": -0.087,
        "meso_heat": 0.45,
        "micro_temp": 0.12,
        "cluster": 0
      },
      ...
    ]
  }
  ```

### Figures
- `figure1_fra_trajectory.png` — FRA over time
- `figure2_phase_portrait.png` — Regime cycles
- `figure3_cluster_evolution.png` — Cluster dynamics
- `figure4_top_assets_heatmap.png` — Asset rankings
- `figure5_cluster_heat_distribution.png` — Heat statistics
- `figure6_composite_distribution.png` — Score statistics
- `figure7_fra_components.png` — Macro/meso/micro breakdown

### Reports
- `analysis_report.md` — Summary statistics and interpretation
- `PHASE3B3_ANALYSIS_REPORT.md` — Detailed analysis report

## Pre-Registered Predictions

| ID | Prediction | PASS Criterion |
|----|------------|----------------|
| P3B3-1 | FRA drifts toward -1 before 2008 and 2020 | FRA < -0.3 within 60 days before crisis |
| P3B3-2 | FRA velocity goes negative before extreme | dFRA/dt < -0.05 at least 2 weeks before FRA < -0.5 |
| P3B3-3 | FRA momentum goes negative before crisis | Momentum < -0.1 at least 2 weeks before FRA < -0.5 |
| P3B3-4 | Financials-heavy cluster has highest heat in 2007 | Cluster with >40% Financials has max \|heat\| by 2007-Q3 |
| P3B3-5 | Cluster count decreases before crises | N_clusters drops by ≥2 in 60 days before crisis |
| P3B3-6 | Phase portrait shows clockwise loops | Visual: trajectory orbits crisis attractor |

## Computation Budget

| Component | Per date | Total (~950 weekly) |
|-----------|----------|---------------------|
| History request | ~0.1s | ~95s |
| Correlation matrix | ~0.02s | ~19s |
| Market Fisher | ~0.3s | ~285s |
| Spectral clustering | ~0.1s | ~95s |
| Cluster Fisher | ~0.5s | ~475s |
| Per-asset temps | ~0.3s | ~285s |
| **Per window** | **~1.3s** | **~21 min** |
| **All 3 windows** | | **~25 min** |
| Post-processing | | ~3 min |
| **Total** | | **~28 min** |

Within 30-minute budget. If tight, drop 60d window (keep 30d and 90d).

## Key Parameters

```python
COMPUTE_FREQUENCY = 5              # every 5 trading days (weekly)
WINDOWS = [30, 60, 90]            # multi-scale windows
LOOKBACK_ZSCORE = 252             # 1-year rolling z-score
MIN_CLUSTER_SIZE = 12             # minimum for valid Fisher
MAX_CLUSTERS = 20                 # cap on spectral k
N_FISHER_SAMPLES = 20             # per cluster
UNIVERSE_SIZE = 200               # top 200 most liquid
```

## Interpretation Guide

### FRA Trajectory Patterns

| FRA | Velocity | Acceleration | Regime |
|-----|----------|-------------|--------|
| Near 0 | ~0 | ~0 | Stable equilibrium |
| Drifting toward -1 | Negative | ~0 | Gradual stress buildup |
| Snapping toward -1 | Strongly negative | Negative | Acute shock onset |
| Near -1 | ~0 | Positive | Crisis plateau, recovery starting |
| Rising from -1 | Positive | Positive | Recovery accelerating |
| Drifting toward +1 | Positive | ~0 | Dispersion building |

### Phase Portrait Interpretation

- **Loops around origin**: Normal regime cycling
- **Spiral toward -1**: Crisis onset
- **Plateau at -1**: Crisis plateau
- **Spiral away from -1**: Recovery phase
- **Clockwise motion**: Typical market dynamics

### Cluster Heat Interpretation

- **Positive heat**: Cluster showing unification (stress)
- **Negative heat**: Cluster showing dispersion (opportunity)
- **Large magnitude**: Anomalous state, actionable signal
- **Cluster reorganization**: Sudden membership changes = regime boundary

## Advanced Usage

### Custom Weights

Adjust macro/meso/micro weights:
```python
composite = (w_macro * fra 
           + w_meso * cluster_heat 
           + w_micro * asset_temp)
```

Default: 0.4, 0.4, 0.2 (macro and meso equally weighted)

### Sector Overlay

After clustering, compute GICS composition:
```
Cluster 3: 60% Financials, 25% Real Estate, 15% Consumer
```

Provides interpretability without constraining clustering.

### Real-Time Monitoring

Stream new data weekly:
1. Update price history
2. Recompute FRA and clusters
3. Generate watchlist
4. Track regime transitions

## Known Limitations

1. **Small clusters**: Spectral clustering on <50 assets may be unstable
2. **Regime persistence**: FRA can stay in extreme states for extended periods
3. **Survivorship bias**: Only includes currently-traded assets
4. **Lookback period**: 252-day z-score requires 1 year of history

## References

- Phase 3B: Market-level Fisher diagnostics
- Phase 3B-2: Sector-level decomposition
- Phase 3B-3: Multi-scale gradient system (this implementation)

---

**Status:** Ready for deployment  
**Target Runtime:** ≤30 minutes weekly  
**Universe:** ~200 most liquid S&P 500 names  
**Output:** Continuous gradients, no binary signals