# DS Framework Phase 1 Validation Results
## Date: 2026-03-04 16:24:42
## Runtime: 70.0 seconds

### 2D Torus (200x200)

| Route | Method | Estimate | Gate [1.85, 2.15] | Status |
|-------|--------|----------|-------------------|--------|
| 1 | Growth Dimension | 1.9754 (R²=1.0000) | [1.85, 2.15] | PASS |
| 2 | Spectral Dimension | 2.0051 (R²=0.9977) | [1.85, 2.15] | PASS |
| 3 | Fisher Info Rank | 2.0000 (PR=2.759) | [1.85, 2.15] | PASS |

Cross-route agreement: 0.0296 (gate: < 0.1) -> PASS

### 3D Torus (25x25x25)

| Route | Method | Estimate | Gate [2.70, 3.30] | Status |
|-------|--------|----------|-------------------|--------|
| 1 | Growth Dimension | 2.7497 (R²=0.9998) | [2.7, 3.3] | PASS |
| 2 | Spectral Dimension | 3.0330 (R²=0.9846) | [2.7, 3.3] | PASS |
| 3 | Fisher Info Rank | 3.0000 (PR=4.169) | [2.7, 3.3] | PASS |

Cross-route agreement: 0.2833 (gate: < 0.15) -> FAIL

### Phase 1 Overall: **PARTIAL**

- **PASS**: All individual gates pass AND both convergence gates pass
- **PARTIAL**: Some routes pass, others marginal
- **FAIL**: Core route(s) fail

### Diagnostics

**2D Torus (200x200) — Fisher sigma sweep:**

| sigma | Mean Rank | Participation Ratio |
|-------|-----------|---------------------|
| 1.5 | 3.00 | 3.482 |
| 2.0 | 2.00 | 3.156 |
| 3.0 | 2.00 | 2.759 |
| 5.0 | 2.00 | 2.437 |
| 8.0 | 2.00 | 2.265 |

**2D Torus (200x200) — Fisher SV profile:** [1.     1.     0.2305 0.1651]

**3D Torus (25x25x25) — Fisher sigma sweep:**

| sigma | Mean Rank | Participation Ratio |
|-------|-----------|---------------------|
| 1.5 | 4.00 | 5.129 |
| 2.0 | 3.00 | 4.714 |
| 3.0 | 3.00 | 4.169 |
| 5.0 | 3.00 | 3.763 |
| 8.0 | 3.00 | 3.583 |

**3D Torus (25x25x25) — Fisher SV profile:** [1.     1.     1.     0.27   0.1717 0.1717]
