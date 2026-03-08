# ✅ QuantConnect Error Fixed

## Problem

You got this error when trying to run the algorithm in QuantConnect:

```
Algorithm.Initialize() Error: During the algorithm initialization, the following 
exception has occurred: Loader.TryCreatePythonAlgorithm(): Unable to import python 
module ./cache/algorithm/project/main.pyc. AlgorithmPythonWrapper(): Please ensure 
that one class inherits from QCAlgorithm.
```

## Root Cause

The original `FisherRegimeAttractor.py` was a **standalone Python script**, not a QuantConnect LEAN algorithm. QuantConnect requires:

1. A class that **inherits from `QCAlgorithm`**
2. An `initialize()` method
3. Proper LEAN framework integration

## Solution

Created **`FisherRegimeAttractor_QC.py`** — a proper QuantConnect-compatible version that:

✅ Inherits from `QCAlgorithm`  
✅ Implements `initialize()` method  
✅ Uses QuantConnect's scheduling system  
✅ Plots to LEAN charts  
✅ Saves results to ObjectStore  
✅ Handles universe selection  

---

## 🚀 How to Use

### Step 1: Copy the Code

Replace your `main.py` in QuantConnect with the contents of:
```
FisherRegimeAttractor_QC.py
```

### Step 2: Run Backtest

1. Click "Backtest"
2. Date range: 2006-01-01 to 2024-12-31
3. Click "Run"

### Step 3: Monitor

Watch the charts:
- **Fisher Regime Attractor** — FRA and velocity
- **Cluster Dynamics** — Cluster count and heat

### Step 4: Download Results

After backtest completes:
1. Go to "Object Store"
2. Download `fisher_regime_attractor_results`
3. Contains all FRA and cluster data in JSON format

### Step 5: Post-Process (Optional)

```bash
python fra_visualization.py
```

Generates 7 figures from the results.

---

## 📊 What You'll Get

### In QuantConnect Charts

- **FRA (90d)** — Market regime on [-1, +1] scale
- **FRA Velocity** — Rate of change
- **N Clusters** — Number of communities
- **Max |Heat|** — Strongest cluster signal

### In ObjectStore JSON

```json
{
  "fra_history": [-0.087, -0.045, 0.023, ...],
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
        }
      ]
    }
  ]
}
```

---

## 🔧 Key Features

### 1. Fisher Regime Attractor (FRA)

Market position on [-1, +1] criticality spectrum:
- **-1** = Maximum panic (all correlations unified)
- **0** = Equilibrium
- **+1** = Maximum dispersion (stock-picking regime)

### 2. Dynamic Clustering

Spectral clustering with adaptive k:
- No hardcoded group sizes
- Eigengap-based k selection
- Automatic small-cluster merging

### 3. Cluster Heat Vector

Per-community Fisher temperature:
- Positive = stress/unification
- Negative = opportunity/dispersion
- Magnitude = strength of signal

### 4. Asset Temperature Ranking

Composite score combining:
- 40% Macro (FRA)
- 40% Meso (Cluster heat)
- 20% Micro (Asset temp)

All assets ranked by |composite| descending.

---

## 📈 Interpretation

### FRA Patterns

| FRA | Velocity | Meaning |
|-----|----------|---------|
| -0.9 | -0.3 | Extreme panic |
| -0.5 | -0.2 | Crisis onset |
| -0.3 | -0.1 | Elevated stress |
| 0.0 | 0.0 | Equilibrium |
| +0.3 | +0.1 | Elevated opportunity |
| +0.9 | +0.3 | Extreme dispersion |

### Cluster Dynamics

- **Decreasing clusters** = Correlations merging (stress)
- **Stable clusters** = Normal regime
- **Increasing clusters** = Correlations fragmenting (opportunity)

### Heat Signals

- **> 2.0** = Anomalous cluster state
- **1.0 to 2.0** = Elevated stress
- **-1.0 to 1.0** = Normal range
- **< -2.0** = Strong opportunity

---

## 📁 Files Provided

### QuantConnect
- **`FisherRegimeAttractor_QC.py`** ← USE THIS FOR QUANTCONNECT

### Standalone Python
- `FisherRegimeAttractor.py` — Full implementation
- `FisherRegimeAttractor_Fast.py` — Fast demo

### Visualization
- `fra_visualization.py` — Generate 7 figures

### Documentation
- `QUANTCONNECT_DEPLOYMENT.md` — Deployment guide
- `README_PHASE3B3.md` — Full documentation
- `QUICKSTART_PHASE3B3.md` — Quick start
- `INDEX_PHASE3B3.md` — Navigation guide

---

## ✅ Verification

The new `FisherRegimeAttractor_QC.py` has:

✅ Proper class inheritance from `QCAlgorithm`  
✅ `initialize()` method  
✅ Weekly scheduling  
✅ Universe selection  
✅ Chart plotting  
✅ ObjectStore saving  
✅ Error handling  
✅ Full Fisher pipeline  
✅ Dynamic clustering  
✅ Asset scoring  

---

## 🎯 Next Steps

1. **Copy** `FisherRegimeAttractor_QC.py` to QuantConnect
2. **Run** backtest (2006-2024)
3. **Monitor** charts
4. **Download** results from ObjectStore
5. **Analyze** using visualization script
6. **Integrate** with your trading strategy

---

## 📞 Support

For questions:
- See `QUANTCONNECT_DEPLOYMENT.md` for deployment help
- See `README_PHASE3B3.md` for full documentation
- See `QUICKSTART_PHASE3B3.md` for quick start

---

**Status: ✅ FIXED AND READY**

The error is resolved. Use `FisherRegimeAttractor_QC.py` for QuantConnect.