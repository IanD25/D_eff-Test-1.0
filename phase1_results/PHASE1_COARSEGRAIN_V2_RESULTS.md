# Phase 1: Coarse-Graining with Scale-Aware Sigma — Results

Generated: 2026-03-05 08:54:10

## 2D Torus Hierarchy

### All Routes — Comparison Table

| Level | Side | Growth | Spectral | σ (fixed) | Rank (fixed) | PR (fixed) | σ (B) | Rank (B) | PR (B) | σ (C) | Rank (C) | PR (C) |
|-------|------|--------|----------|---------|------------|----------|---------|------------|----------|---------|------------|----------|
| 0 | 128 | 1.953 | 2.789 | 3.000 | 2.0 | 2.759 | 16.000 | 2.0 | 2.134 | 8.000 | 2.0 | 2.265 |
| 1 | 64 | 1.906 | 2.790 | 3.000 | 2.0 | 2.759 | 8.000 | 2.0 | 2.276 | 4.000 | 2.0 | 2.556 |
| 2 | 32 | 1.805 | 2.791 | 3.000 | 2.0 | 2.767 | 4.000 | 2.0 | 2.579 | 2.000 | 2.0 | 3.157 |
| 3 | 16 | N/A | 2.777 | 3.000 | 2.0 | 2.885 | 2.000 | 2.0 | 3.202 | 1.000 | 3.0* | 3.647* |
| 4 | 8 | N/A | 2.713 | 3.000 | 2.0 | 3.364 | 1.000 | 3.0* | 3.660* | 0.500 | 1.0* | 2.054* |
| 5 | 4 | N/A | 1.865 | 3.000 | 1.0 | 3.991 | 0.500 | 1.0* | 1.991* | 0.250 | 1.0* | 1.114* |
| 6 | 2 | N/A | N/A | 3.000 | N/A | N/A | 0.250 | N/A | N/A | 0.125 | N/A | N/A |

*Asterisk (*) indicates σ < 2.0 — results may be unreliable.*

### P7 Monotonicity (Schedule B, reliable levels)
- PR at first reliable level: 2.134
- PR at last reliable level: 3.202
- Change: +1.068
- Monotonically non-increasing (ε=0.05): **FAIL**

---

## 3D Torus Hierarchy

### All Routes — Comparison Table

| Level | Side | Growth | Spectral | σ (fixed) | Rank (fixed) | PR (fixed) | σ (B) | Rank (B) | PR (B) | σ (C) | Rank (C) | PR (C) |
|-------|------|--------|----------|---------|------------|----------|---------|------------|----------|---------|------------|----------|
| 0 | 32 | 2.769 | 3.910 | 3.000 | 3.0 | 4.146 | 4.000 | 3.0 | 3.867 | 2.000 | 3.0 | 4.708 |
| 1 | 16 | 2.508 | 3.901 | 3.000 | 3.0 | 4.323 | 2.000 | 3.0 | 4.776 | 1.000 | 4.0* | 5.050* |
| 2 | 8 | N/A | 3.816 | 3.000 | 3.0 | 5.039 | 1.000 | 4.0* | 5.067* | 0.500 | 1.0* | 2.290* |
| 3 | 4 | N/A | 3.709 | 3.000 | 1.0 | 5.978 | 0.500 | 1.0* | 2.237* | 0.250 | 1.0* | 1.128* |
| 4 | 2 | N/A | N/A | 3.000 | 1.0 | 2.996 | 0.250 | 1.0* | 1.100* | 0.125 | 1.0* | 1.002* |

*Asterisk (*) indicates σ < 2.0 — results may be unreliable.*

### P7 Monotonicity (Schedule B, reliable levels)
- PR at first reliable level: 3.867
- PR at last reliable level: 4.776
- Change: +0.909
- Monotonically non-increasing (ε=0.05): **FAIL**

---
