# DS Framework Phase 1 Validation Results
## Date: 2026-03-04 18:13:40
## Runtime: 39.2 seconds

### 2D Torus (200x200)

| Route | Method | Estimate | Gate [1.85, 2.15] | Status |
|-------|--------|----------|-------------------|--------|
| 1 | Growth Dimension | 1.9754 (R²=1.0000) | [1.85, 2.15] | PASS |
| 2 | Spectral Dimension | 2.0051 (R²=0.9977) | [1.85, 2.15] | PASS |
| 3 | Fisher Info Rank | 2.0000 (PR=2.759) | [1.85, 2.15] | PASS |

Cross-route agreement: 0.0296 (gate: < 0.1) -> PASS

### 3D Torus (50x50x50)

| Route | Method | Estimate | Gate [2.85, 3.15] | Status |
|-------|--------|----------|-------------------|--------|
| 1 | Growth Dimension | 2.8898 (R²=1.0000) | [2.85, 3.15] | PASS |
| 2 | Spectral Dimension | 2.9399 (R²=0.9843) | [2.85, 3.15] | PASS |
| 3 | Fisher Info Rank | 3.0000 (PR=4.134) | [2.85, 3.15] | PASS |

Cross-route agreement: 0.1102 (gate: < 0.15) -> PASS

### Phase 1 Overall: **PASS**

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

**3D Torus (50x50x50) — Fisher sigma sweep:**

| sigma | Mean Rank | Participation Ratio |
|-------|-----------|---------------------|
| 1.5 | 4.00 | 5.127 |
| 2.0 | 3.00 | 4.707 |
| 3.0 | 3.00 | 4.134 |
| 5.0 | 3.00 | 3.665 |
| 8.0 | 3.00 | 3.436 |

**3D Torus (50x50x50) — Fisher SV profile:** [1.     1.     1.     0.2633 0.1653 0.1653]
