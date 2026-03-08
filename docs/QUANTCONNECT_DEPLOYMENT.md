# DS Phase 3B-3: QuantConnect Deployment Guide

## ✅ Fixed: QuantConnect-Compatible Version

The error you encountered was because the previous version wasn't a proper QuantConnect algorithm. I've created **`FisherRegimeAttractor_QC.py`** which properly inherits from `QCAlgorithm`.

---

## 🚀 How to Deploy to QuantConnect

### Step 1: Copy the Code

1. Open your QuantConnect project
2. Create a new file or replace `main.py` with the contents of `FisherRegimeAttractor_QC.py`
3. Make sure the class name is `FisherRegimeAttractor` (it inherits from `QCAlgorithm`)

### Step 2: Verify the Class Structure

The algorithm must have:
```python
class FisherRegimeAttractor(QCAlgorithm):
    def initialize(self):
        # Setup code
    
    def compute_regime(self):
        # Main computation
```

✅ This is now correct in `FisherRegimeAttractor_QC.py`

### Step 3: Run the Backtest

1. Click "Backtest" in QuantConnect
2. Select date range: 2006-01-01 to 2024-12-31
3. Click "Run"

### Step 4: Monitor Progress

The algorithm will:
- Warm up for 400 days
- Compute weekly (every Monday)
- Plot FRA and cluster dynamics
- Save results to ObjectStore

### Step 5: Download Results

After backtest completes:
1. Go to "Object Store"
2. Download `fisher_regime_attractor_results`
3. This contains all FRA and cluster data

---

## 📊 What the Algorithm Does

### Weekly Computation

Every Monday after market open:

1. **Get active securities** — Top 200 by dollar volume
2. **Compute correlation matrix** — 90-day returns
3. **Market Fisher** — SV2/SV1 on full universe
4. **FRA** — Convert to [-1, +1] scale
5. **Clustering** — Spectral clustering with adaptive k
6. **Cluster heat** — Per-community Fisher temperature
7. **Asset scores** — Composite ranking

### Output

**Charts (visible in QuantConnect):**
- Fisher Regime Attractor (FRA and velocity)
- Cluster Dynamics (count and max heat)

**ObjectStore (download after backtest):**
- `fisher_regime_attractor_results` — JSON with all data

---

## 🔧 Configuration

### Key Parameters (in code)

```python
WINDOWS = [90]              # 90-day window (can add 30, 60)
LOOKBACK_ZSCORE = 252       # 1-year baseline
MIN_CLUSTER_SIZE = 12       # Minimum cluster size
MAX_CLUSTERS = 20           # Maximum clusters
N_FISHER_SAMPLES = 20       # Samples per cluster
```

### Adjust Universe Size

In `coarse_selection()`:
```python
return [x.symbol for x in sorted_by_volume[:200]]  # Change 200 to desired size
```

### Change Computation Frequency

In `initialize()`:
```python
self.schedule.on(
    self.date_rules.every(DayOfWeek.MONDAY),  # Change to other days
    self.time_rules.after_market_open("SPY", 30),
    self.compute_regime
)
```

---

## 📈 Expected Output

### Charts in QuantConnect

**Fisher Regime Attractor:**
- FRA (90d) — Main regime signal
- FRA Velocity — Rate of change

**Cluster Dynamics:**
- N Clusters — Number of communities
- Max |Heat| — Strongest cluster signal

### ObjectStore JSON

```json
{
  "fra_history": [
    -0.087, -0.045, 0.023, ...
  ],
  "cluster_history": [
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
    },
    ...
  ]
}
```

---

## ✅ Troubleshooting

### Error: "Please ensure that one class inherits from QCAlgorithm"

**Solution:** Use `FisherRegimeAttractor_QC.py` (the new version)

### Error: "Unable to import python module"

**Solution:** 
1. Make sure the file is named `main.py` or update the reference
2. Check that the class name matches the file name
3. Verify no syntax errors in the code

### Algorithm runs but no data

**Solution:**
1. Check that universe has enough securities (need ≥50)
2. Verify date range has sufficient data
3. Check logs for error messages

### Runtime too long

**Solution:**
1. Reduce universe size (change 200 to 100)
2. Increase computation frequency (change MONDAY to every 2 weeks)
3. Drop 30d/60d windows (keep only 90d)

---

## 📊 Post-Processing

After downloading results from ObjectStore:

### Option 1: Use Visualization Script

```bash
# Save ObjectStore JSON to phase3b3_results/cluster_history.json
python fra_visualization.py
```

This generates 7 figures.

### Option 2: Manual Analysis

```python
import json

with open('fisher_regime_attractor_results.json', 'r') as f:
    data = json.load(f)

# Access FRA history
fra_values = data['fra_history']

# Access cluster data
for point in data['cluster_history']:
    date = point['date']
    fra = point['fra']
    n_clusters = point['n_clusters']
    top_10 = point['top_10_assets']
    
    print(f"{date}: FRA={fra:.3f}, Clusters={n_clusters}")
    for asset in top_10:
        print(f"  {asset['name']}: {asset['composite']:.3f}")
```

---

## 🎯 Key Metrics to Monitor

### FRA Values

- **< -0.5**: Crisis regime
- **-0.3 to -0.5**: Elevated stress
- **-0.3 to +0.3**: Normal regime
- **+0.3 to +0.5**: Elevated opportunity
- **> +0.5**: Dispersion regime

### Cluster Count

- **Decreasing**: Correlations merging (stress)
- **Stable**: Normal regime
- **Increasing**: Correlations fragmenting (opportunity)

### Max Heat

- **> 2.0**: Anomalous cluster state
- **1.0 to 2.0**: Elevated stress
- **-1.0 to 1.0**: Normal range
- **< -2.0**: Strong opportunity

---

## 📝 Example Interpretation

### Scenario 1: Pre-Crisis Pattern

```
Date        FRA    N_Clusters  Max_Heat  Interpretation
2024-01-08  -0.15  8           +0.3      Stress building
2024-01-15  -0.35  7           +0.8      Crisis onset
2024-01-22  -0.65  6           +1.2      Acute crisis
2024-01-29  -0.80  5           +1.5      Crisis plateau
```

**Signal**: FRA drifting then snapping toward -1, clusters merging, heat rising

### Scenario 2: Recovery Pattern

```
Date        FRA    N_Clusters  Max_Heat  Interpretation
2024-02-05  -0.75  5           +1.2      Crisis plateau
2024-02-12  -0.40  6           +0.8      Recovery starting
2024-02-19  -0.10  7           +0.2      Rapid recovery
2024-02-26  +0.20  8           -0.3      Opportunity regime
```

**Signal**: FRA rising, clusters expanding, heat decreasing

---

## 🔗 Integration with Trading

### Use FRA for Position Sizing

```python
# In your trading algorithm
fra = current_fra_value

if fra < -0.5:
    # Crisis: reduce position size
    position_size = base_size * 0.5
elif fra > 0.3:
    # Opportunity: increase position size
    position_size = base_size * 1.5
else:
    # Normal: standard size
    position_size = base_size
```

### Use Cluster Heat for Sector Rotation

```python
# Identify stressed sectors
for cluster_id, heat in cluster_heats.items():
    if heat > 1.0:
        # This cluster is stressed
        # Reduce exposure to this sector
        pass
    elif heat < -1.0:
        # This cluster has opportunity
        # Increase exposure to this sector
        pass
```

### Use Asset Scores for Stock Selection

```python
# Get top 10 most actionable assets
top_10_assets = cluster_history[-1]['top_10_assets']

for asset in top_10_assets:
    name = asset['name']
    composite = asset['composite']
    
    if composite > 0.3:
        # Stress signal: consider reducing
        pass
    elif composite < -0.3:
        # Opportunity signal: consider increasing
        pass
```

---

## 📚 Files Reference

| File | Purpose |
|------|---------|
| `FisherRegimeAttractor_QC.py` | ✅ QuantConnect algorithm (USE THIS) |
| `FisherRegimeAttractor.py` | Standalone Python version |
| `FisherRegimeAttractor_Fast.py` | Fast demo version |
| `fra_visualization.py` | Post-processing visualization |
| `README_PHASE3B3.md` | Full documentation |

---

## ✅ Deployment Checklist

- [ ] Copy `FisherRegimeAttractor_QC.py` to QuantConnect
- [ ] Verify class inherits from `QCAlgorithm`
- [ ] Set date range: 2006-01-01 to 2024-12-31
- [ ] Run backtest
- [ ] Monitor charts during execution
- [ ] Download results from ObjectStore
- [ ] Run post-processing visualization
- [ ] Analyze results
- [ ] Integrate with trading strategy

---

## 🚀 Next Steps

1. **Deploy**: Copy code to QuantConnect and run
2. **Monitor**: Watch charts during backtest
3. **Download**: Get results from ObjectStore
4. **Analyze**: Run visualization script
5. **Integrate**: Use FRA in your trading strategy
6. **Optimize**: Tune parameters for your use case

---

**Status: ✅ Ready for QuantConnect Deployment**

The `FisherRegimeAttractor_QC.py` file is now properly formatted for QuantConnect and should run without errors.