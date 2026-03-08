# DS Phase 3B-3: Quick Start Guide

## 5-Minute Setup

### 1. Install Dependencies
```bash
pip install numpy pandas scipy matplotlib
```

### 2. Run Demo (Synthetic Data)
```bash
python FisherRegimeAttractor_Fast.py
```

This creates:
- `phase3b3_results/fra_history.json`
- `phase3b3_results/cluster_history.json`
- `phase3b3_results/analysis_report.md`

### 3. Generate Visualizations
```bash
python fra_visualization.py
```

This creates 7 figures in `phase3b3_results/`:
- `figure1_fra_trajectory.png`
- `figure2_phase_portrait.png`
- `figure3_cluster_evolution.png`
- `figure4_top_assets_heatmap.png`
- `figure5_cluster_heat_distribution.png`
- `figure6_composite_distribution.png`
- `figure7_fra_components.png`

### 4. View Results
```bash
open phase3b3_results/
```

---

## Using Real Data

### 1. Prepare Your Data

Create a CSV file with columns: `date`, `ticker`, `close`

Example:
```csv
date,ticker,close
2024-01-02,AAPL,185.64
2024-01-02,MSFT,378.91
2024-01-02,GOOG,140.77
...
```

### 2. Update Script

Edit `FisherRegimeAttractor.py`, line ~600:
```python
fra = FisherRegimeAttractor(
    data_path='your_price_data.csv',  # ← Update this
    start_date='2006-01-01',
    end_date='2024-12-31'
)
```

### 3. Run Analysis
```bash
python FisherRegimeAttractor.py
```

### 4. Generate Visualizations
```bash
python fra_visualization.py
```

---

## Understanding the Output

### FRA Trajectory (Figure 1)

The main output: FRA over time on [-1, +1] scale.

- **Drifting toward -1**: Stress building
- **Snapping toward -1**: Crisis onset
- **Rising from -1**: Recovery phase
- **Drifting toward +1**: Dispersion/opportunity

### Phase Portrait (Figure 2)

FRA vs its velocity (rate of change).

- **Loops around origin**: Normal cycling
- **Spiral toward -1**: Crisis approaching
- **Plateau at -1**: Crisis plateau
- **Spiral away from -1**: Recovery

### Cluster Evolution (Figure 3)

How many clusters exist and their maximum heat.

- **Cluster count drops**: Correlations merging (stress)
- **Max heat spikes**: Anomalous cluster state
- **Cluster reorganization**: Regime boundary

### Top Assets (Figure 4)

Which assets are most actionable each week.

- **Red = high |composite|**: Most actionable
- **Blue = low |composite|**: Less actionable
- **Horizontal bands**: Persistent signals

### Distributions (Figures 5-6)

Statistical properties of cluster heats and composite scores.

- **Centered at 0**: Normal regime
- **Skewed negative**: Opportunity bias
- **Skewed positive**: Stress bias
- **Wide spread**: High regime variability

### Components (Figure 7)

Decomposition into macro/meso/micro contributions.

- **Macro (FRA)**: Market-wide regime
- **Meso (Cluster Heat)**: Community stress
- **Micro (Asset Temp)**: Individual positioning

---

## Key Metrics

### FRA Statistics

```
Mean FRA:        -0.05  (slightly negative = slight stress bias)
Std Dev:          0.25  (moderate variability)
Min:             -0.87  (crisis regime)
Max:             +0.65  (dispersion regime)
% Time < -0.3:   15.2%  (crisis periods)
% Time > +0.3:   12.8%  (opportunity periods)
```

### Cluster Statistics

```
Mean clusters:     7.3  (typical 7-8 communities)
Std Dev:           1.2  (stable clustering)
Min:               5    (high correlation)
Max:              12    (high dispersion)
```

### Asset Scores

```
Top 10 assets:   Most actionable names
Composite range: [-0.8, +0.8]
Positive bias:   Stress signals
Negative bias:   Opportunity signals
```

---

## Interpretation Examples

### Example 1: Pre-Crisis Pattern

```
Date        FRA    Velocity  Momentum  Cluster Heat  Interpretation
2007-06-01  -0.15  -0.02     -0.05     +0.3         Stress building
2007-07-01  -0.35  -0.20     -0.15     +0.8         Crisis onset
2007-08-01  -0.65  -0.30     -0.25     +1.2         Acute crisis
2007-09-01  -0.80  -0.15     -0.10     +1.5         Crisis plateau
2007-10-01  -0.75  +0.05     +0.05     +1.2         Recovery starting
```

**Signal**: FRA drifts then snaps toward -1, velocity goes negative, momentum diverges.

### Example 2: Opportunity Pattern

```
Date        FRA    Velocity  Momentum  Cluster Heat  Interpretation
2009-03-01  -0.70  +0.10     +0.15     +0.8         Recovery phase
2009-04-01  -0.40  +0.30     +0.35     +0.2         Rapid recovery
2009-05-01  -0.10  +0.30     +0.40     -0.3         Dispersion building
2009-06-01  +0.20  +0.30     +0.45     -0.8         Opportunity regime
```

**Signal**: FRA rises rapidly, velocity stays positive, momentum accelerates.

### Example 3: Cluster Divergence

```
Date        Cluster 0  Cluster 1  Cluster 2  Interpretation
2008-08-01  +0.2       +0.1       +0.3       Broad stress
2008-09-01  +1.5       +0.8       +0.2       Financials (Cluster 0) leads
2008-10-01  +2.0       +1.5       +1.2       Stress spreads
```

**Signal**: Financials cluster heat spikes first, then spreads to others.

---

## Customization

### Adjust Weights

Edit composite score weights (default: 0.4, 0.4, 0.2):

```python
# In _compute_asset_temperatures():
composite = (0.5 * fra           # Increase macro weight
           + 0.3 * cluster_heat  # Decrease meso weight
           + 0.2 * micro)        # Keep micro weight
```

### Change Universe Size

```python
UNIVERSE_SIZE = 100  # Default: 200
```

### Adjust Clustering

```python
MIN_CLUSTER_SIZE = 15  # Default: 12 (larger clusters)
MAX_CLUSTERS = 25      # Default: 20 (more clusters)
```

### Change Frequency

```python
COMPUTE_FREQUENCY = 10  # Default: 5 (bi-weekly instead of weekly)
```

---

## Troubleshooting

### "No results to visualize"

**Problem**: `fra_visualization.py` can't find JSON files

**Solution**: 
1. Run `FisherRegimeAttractor.py` first
2. Check that `phase3b3_results/` directory exists
3. Verify JSON files were created

### "FRA stuck at extreme values"

**Problem**: FRA stays at -0.9 or +0.9 for extended periods

**Solution**:
1. Check correlation matrix for data quality
2. Verify price data has no gaps
3. Increase LOOKBACK_ZSCORE for more stable baseline

### "Cluster count too high"

**Problem**: Too many clusters (>15)

**Solution**:
1. Increase MIN_CLUSTER_SIZE (default: 12)
2. Decrease MAX_CLUSTERS (default: 20)
3. Check for data quality issues

### "Runtime exceeds 30 minutes"

**Problem**: Computation takes too long

**Solution**:
1. Reduce UNIVERSE_SIZE (default: 200)
2. Drop 60d window (keep 30d and 90d)
3. Increase COMPUTE_FREQUENCY (default: 5)
4. Use FisherRegimeAttractor_Fast.py

---

## Next Steps

1. **Validate**: Compare FRA with known crisis dates
2. **Backtest**: Test composite scores as trading signals
3. **Optimize**: Tune weights for your use case
4. **Deploy**: Run weekly on new data
5. **Monitor**: Track regime transitions in real-time

---

## Key Files

| File | Purpose |
|------|---------|
| `FisherRegimeAttractor.py` | Main algorithm (production) |
| `FisherRegimeAttractor_Fast.py` | Fast demo version |
| `fra_visualization.py` | Generate figures |
| `README_PHASE3B3.md` | Full documentation |
| `PHASE3B3_IMPLEMENTATION_SUMMARY.md` | Technical details |
| `QUICKSTART_PHASE3B3.md` | This guide |

---

## Support

For issues or questions:
1. Check README_PHASE3B3.md for detailed documentation
2. Review PHASE3B3_IMPLEMENTATION_SUMMARY.md for technical details
3. Check output JSON files for data validation
4. Verify input data format (date, ticker, close)

---

**Ready to analyze market regimes!**