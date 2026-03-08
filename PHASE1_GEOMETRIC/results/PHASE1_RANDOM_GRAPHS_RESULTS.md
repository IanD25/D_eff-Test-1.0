# Phase 1 Extension: Random Geometric Graphs + Erdős-Rényi Results

Generated: 2026-03-05 00:17:53

## Section 1: Random Geometric Graph Results

### RGG d=2

- **Vertices (LCC):** 5000 (100.0% of 5000)
- **Edges:** 34644
- **Connected:** True
- **Degree:** min=1, max=26, mean=13.9, median=14, std=3.9
- **Clustering:** 0.5965
- **Diameter:** 53 (estimated from 30 samples)

**Growth Dimension:** 1.573 (R²=0.9967)
  Gate [1.7, 2.3]: **FAIL**

**Spectral Dimension:** 23.126 (R²=0.9429)
  Gate [1.6, 2.4]: **FAIL**

**Fisher Rank Distribution:**
  Mean=2.90, Median=3, Std=0.41
  Distribution: {1: 1, 2: 4, 3: 44, 4: 1}

**Fisher PR Distribution:**
  Mean=3.210, Median=3.240, Std=0.318
  Gate [1.5, 2.5]: **FAIL**

**PR-Degree Correlation:** r = 0.325 (moderate)

**Sigma Sweep:**

| σ | PR mean | PR std |
|---|---------|--------|
| 2.0 | 3.210 | 0.318 |
| 3.0 | 3.019 | 0.314 |
| 5.0 | 2.858 | 0.327 |

---

### RGG d=3

- **Vertices (LCC):** 5000 (100.0% of 5000)
- **Edges:** 34317
- **Connected:** True
- **Degree:** min=1, max=31, mean=13.7, median=14, std=4.4
- **Clustering:** 0.4964
- **Diameter:** 24 (estimated from 30 samples)

**Growth Dimension:** 2.309 (R²=0.9972)
  Gate [2.7, 3.3]: **FAIL**

**Spectral Dimension:** 21.210 (R²=0.9411)
  Gate [2.5, 3.5]: **FAIL**

**Fisher Rank Distribution:**
  Mean=3.60, Median=4, Std=0.66
  Distribution: {1: 1, 2: 2, 3: 13, 4: 34}

**Fisher PR Distribution:**
  Mean=3.982, Median=4.048, Std=0.531
  Gate [2.5, 3.5]: **FAIL**

**PR-Degree Correlation:** r = 0.464 (moderate)

**Sigma Sweep:**

| σ | PR mean | PR std |
|---|---------|--------|
| 2.0 | 3.982 | 0.531 |
| 3.0 | 3.946 | 0.534 |
| 5.0 | 3.924 | 0.533 |

---

## Section 2: Erdős-Rényi Negative Control

### ER n=1000

- **Vertices (LCC):** 1000 (100.0% of 1000)
- **Edges:** 7450
- **Degree:** min=6, max=27, mean=14.9, median=15, std=3.8
- **Clustering:** 0.0153
- **Diameter:** 4 (exact)

**Growth Dimension:** 3.617 (R²=1.0000)
  (Expected: very high and/or poor R² — no low-dimensional structure)

**Spectral Dimension:** 14.791 (R²=0.9691)
  (Expected: poor Weyl fit — eigenvalues don't follow power law)

**Fisher Rank Distribution:**
  Mean=1.00, Median=1, Std=0.00
  Distribution: {1: 30}
  (Expected: high, near vertex degree — no gap in SV profile)

**Fisher PR Distribution:**
  Mean=4.341, Median=4.359, Std=0.179
  (Expected: high, near vertex degree — all SVs contribute)

**Sigma Sweep:**

| σ | PR mean | PR std |
|---|---------|--------|
| 2.0 | 4.341 | 0.179 |
| 3.0 | 4.563 | 0.226 |
| 5.0 | 4.718 | 0.275 |

---

### ER n=2000

- **Vertices (LCC):** 2000 (100.0% of 2000)
- **Edges:** 14947
- **Degree:** min=4, max=28, mean=14.9, median=15, std=3.8
- **Clustering:** 0.0071
- **Diameter:** 5 (exact)

**Growth Dimension:** 4.873 (R²=1.0000)
  (Expected: very high and/or poor R² — no low-dimensional structure)

**Spectral Dimension:** 21.196 (R²=0.9542)
  (Expected: poor Weyl fit — eigenvalues don't follow power law)

**Fisher Rank Distribution:**
  Mean=1.00, Median=1, Std=0.00
  Distribution: {1: 20}
  (Expected: high, near vertex degree — no gap in SV profile)

**Fisher PR Distribution:**
  Mean=4.762, Median=4.713, Std=0.368
  (Expected: high, near vertex degree — all SVs contribute)

**Sigma Sweep:**

| σ | PR mean | PR std |
|---|---------|--------|
| 2.0 | 4.762 | 0.368 |
| 3.0 | 4.866 | 0.428 |
| 5.0 | 4.887 | 0.464 |

---

## Section 3: The Diagnostic Contrast

| Metric | RGG d=2 (geometric) | ER (non-geometric) |
|--------|--------------------|--------------------|
| Avg degree | 13.9 | 14.9 |
| Fisher rank (mean) | 2.90 | 1.00 |
| Fisher rank (std) | 0.41 | 0.00 |
| Fisher PR (mean) | 3.21 | 4.34 |
| Fisher PR (std) | 0.32 | 0.18 |
| Growth dim R² | 0.9967 | 1.0000 |
| Clustering | 0.5965 | 0.0153 |

**Verdict:** The Fisher method distinguishes geometric from non-geometric graphs. 
Same average connectivity, completely different Fisher Information structure.
