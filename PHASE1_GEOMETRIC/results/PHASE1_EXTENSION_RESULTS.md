# DS Phase 1 Extension Results
## Date: 2026-03-04 18:53:30
## Runtime: 22.6 seconds

## Test 2A: Sierpinski Gasket

Reference dimensions: d_H = 1.5850, d_S = 1.3652

### Level 6 (1095 vertices, 1092 interior)

| Route | Method | Estimate | Gate | Status |
|-------|--------|----------|------|--------|
| 1 | Growth Dimension | 1.4853 (R²=0.9993) | [1.40, 1.75] | PASS |
| 2 | Spectral Dimension | 1.4829 (R²=0.9522) | [1.20, 1.55] | PASS |
| 3 | Fisher PR (mean) | 2.1235 ± 0.1709 | — | — |
| 3 | Fisher Rank (mode) | 3 | — | — |

Fisher SV profile: [1.     0.5174 0.1228 0.026 ]
Fisher SV std:     [0.     0.1525 0.0293 0.0238]

Per-vertex Fisher (all 1092 interior): mean PR = 2.0569 ± 0.1756
Per-vertex rank distribution: {1: 129, 2: 240, 3: 723}

**Sigma sensitivity:**

| sigma | Mean Rank | Mean PR | PR Std |
|-------|-----------|---------|--------|
| 1.5 | 2.73 | 2.4252 | 0.1833 |
| 2.0 | 2.60 | 2.2972 | 0.2046 |
| 3.0 | 2.73 | 2.1385 | 0.1901 |
| 5.0 | 2.53 | 2.0019 | 0.1402 |
| 8.0 | 2.00 | 1.9255 | 0.1144 |
| 12.0 | 1.87 | 1.8777 | 0.1426 |

### Level 7 (3282 vertices, 3279 interior)

| Route | Method | Estimate | Gate | Status |
|-------|--------|----------|------|--------|
| 1 | Growth Dimension | 1.5198 (R²=0.9992) | [1.40, 1.75] | PASS |
| 2 | Spectral Dimension | 1.4464 (R²=0.9525) | [1.20, 1.55] | PASS |
| 3 | Fisher PR (mean) | 2.1050 ± 0.1214 | — | — |
| 3 | Fisher Rank (mode) | 3 | — | — |

Fisher SV profile: [1.     0.5134 0.1144 0.0193]
Fisher SV std:     [0.     0.1057 0.0262 0.0172]

Per-vertex Fisher (all 3279 interior): mean PR = 2.0563 ± 0.1739
Per-vertex rank distribution: {1: 372, 2: 726, 3: 2181}

**Sigma sensitivity:**

| sigma | Mean Rank | Mean PR | PR Std |
|-------|-----------|---------|--------|
| 1.5 | 2.87 | 2.4227 | 0.1226 |
| 2.0 | 2.73 | 2.2892 | 0.1391 |
| 3.0 | 2.87 | 2.1287 | 0.1318 |
| 5.0 | 2.60 | 1.9855 | 0.1023 |
| 8.0 | 2.07 | 1.9010 | 0.1055 |
| 12.0 | 1.93 | 1.8607 | 0.1383 |

## Test 3A: 4D Torus

### 4D Torus (n=12, 20736 vertices)

| Route | Method | Estimate | Gate | Status |
|-------|--------|----------|------|--------|
| 1 | Growth Dimension | 3.3758 (R²=0.9996) | [3.70, 4.30] | FAIL |
| 2 | Spectral Dimension | 3.9439 (R²=0.9485) | [3.70, 4.30] | PASS |
| 3 | Fisher Rank | 4.0000 | rank=4 | PASS |

Cross-route agreement: 0.6242 (gate < 0.20) → FAIL

Fisher SV profile: [1.     1.     1.     1.     0.3813 0.242  0.242  0.242 ]

**4D Torus sigma sweep:**

| sigma | Mean Rank | PR |
|-------|-----------|-----|
| 1.5 | 5.00 | 6.8086 |
| 2.0 | 5.00 | 6.4929 |
| 3.0 | 4.00 | 6.0369 |
| 5.0 | 4.00 | 5.7306 |
| 8.0 | 4.00 | 5.6151 |
