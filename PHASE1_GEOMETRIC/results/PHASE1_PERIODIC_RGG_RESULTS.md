# Phase 1: Periodic RGG + Boundary Analysis — Results

Generated: 2026-03-05 00:50:03

## Test A: Periodic Random Geometric Graphs

### Periodic RGG d=2

- **Vertices (LCC):** 5000 (100.0%)
- **Edges:** 35415
- **Degree:** min=4, max=26, mean=14.2, median=14, std=3.7
- **Clustering:** 0.5876
- **Diameter:** 30 (est.)

**Growth Dimension:** 2.011 (R²=1.0000)
  Gate [1.7, 2.3]: **PASS**

**Spectral Dimension:** 23.985 (R²=0.9420)

**Fisher Rank Distribution:**
  Mean=2.90, Median=3, Std=0.41
  Distribution: {1: 1, 2: 4, 3: 44, 4: 1}

**Fisher PR Distribution:**
  Mean=3.094, Median=3.144, Std=0.394

**PR-Degree Correlation:** r = 0.441

**Sigma Sweep:**

| σ | PR mean | PR std | Rank mean | Rank mode |
|---|---------|--------|-----------|-----------|
| 2.0 | 3.094 | 0.394 | 2.90 | 3 |
| 3.0 | 2.917 | 0.374 | 3.42 | 3 |
| 5.0 | 2.791 | 0.358 | 4.64 | 3 |
| 8.0 | 2.745 | 0.354 | 5.76 | 3 |

---

### Periodic RGG d=3

- **Vertices (LCC):** 5000 (100.0%)
- **Edges:** 38315
- **Degree:** min=5, max=34, mean=15.3, median=15, std=4.0
- **Clustering:** 0.4706
- **Diameter:** 13 (est.)

**Growth Dimension:** 2.774 (R²=0.9992)
  Gate [2.7, 3.3]: **PASS**

**Spectral Dimension:** 22.667 (R²=0.9522)

**Fisher Rank Distribution:**
  Mean=4.70, Median=4, Std=2.75
  Distribution: {2: 1, 3: 8, 4: 35, 5: 1, 12: 2, 13: 2, 14: 1}

**Fisher PR Distribution:**
  Mean=4.189, Median=4.125, Std=0.488

**PR-Degree Correlation:** r = 0.628

**Sigma Sweep:**

| σ | PR mean | PR std | Rank mean | Rank mode |
|---|---------|--------|-----------|-----------|
| 2.0 | 4.189 | 0.488 | 4.70 | 4 |
| 3.0 | 4.227 | 0.489 | 5.08 | 4 |
| 5.0 | 4.286 | 0.487 | 5.22 | 4 |
| 8.0 | 4.328 | 0.484 | 5.24 | 4 |

---

### Comparison: Periodic vs Bounded RGG

| Metric | Bounded d=2 | Periodic d=2 | Bounded d=3 | Periodic d=3 |
|--------|-------------|--------------|-------------|--------------|
| Growth dim | 1.573 | 2.011 | 2.309 | 2.774 |
| Fisher rank (mean) | 2.90 | 2.90 | 3.60 | 4.70 |
| Fisher rank (mode) | 3 | 3 | 4 | 4 |
| Fisher PR (σ=2) | 3.210 | 3.094 | 3.982 | 4.189 |
| PR-deg corr | 0.325 | 0.441 | 0.464 | 0.628 |
| Degree std | 3.9 | 3.7 | 4.4 | 4.0 |

## Test B: Interior vs Boundary Analysis (Bounded RGG)

### Bounded RGG d=2

- Interior: 2473 vertices (49.5%), margin=0.15
- Boundary: 2526 vertices (50.5%)

| Metric | Interior | Boundary | Full Sample |
|--------|----------|----------|-------------|
| N vertices | 50 | 50 | 100 |
| Mean degree | 14.4 | 12.9 | 13.7 |
| Fisher rank (mean) | 2.74 | 2.88 | 2.81 |
| Fisher rank (mode) | 3 | 3 | 3 |
| Fisher PR (mean) | 3.039 | 2.998 | 3.018 |
| Fisher PR (std) | 0.508 | 0.500 | 0.504 |

---

### Bounded RGG d=3

- Interior: 1679 vertices (33.6%), margin=0.15
- Boundary: 3320 vertices (66.4%)

| Metric | Interior | Boundary | Full Sample |
|--------|----------|----------|-------------|
| N vertices | 50 | 50 | 100 |
| Mean degree | 15.4 | 13.2 | 14.3 |
| Fisher rank (mean) | 5.68 | 4.74 | 5.21 |
| Fisher rank (mode) | 4 | 4 | 4 |
| Fisher PR (mean) | 4.184 | 3.845 | 4.015 |
| Fisher PR (std) | 0.491 | 0.740 | 0.650 |

---
