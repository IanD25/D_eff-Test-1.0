# DS Phase 1: Coarse-Graining Test (P7 / DPI) Results
## Date: 2026-03-04 20:55:21
## Runtime: 11.3 seconds

### 2D Torus, 2×2 Blocking (128² → 4²)

| Level | Size | Growth | Spectral | Fisher Rank | Fisher PR (σ=3) | ΔPR | Status |
|-------|------|--------|----------|-------------|-----------------|-----|--------|
| 0 | 128² | 1.953 | 2.789 | 2.00 | 2.7586 | — | — |
| 1 | 64² | 1.906 | 2.790 | 2.00 | 2.7586 | +0.0000 | MARGINAL |
| 2 | 32² | 1.805 | 2.791 | 2.00 | 2.7668 | +0.0081 | MARGINAL |
| 3 | 16² | N/A | 2.777 | 2.00 | 2.8852 | +0.1185 | FAIL |
| 4 | 8² | N/A | 2.713 | 2.00 | 3.3636 | +0.4784 | FAIL |
| 5 | 4² | N/A | 1.865 | 1.00 | 3.9911 | +0.6275 | FAIL |

**P7 Verdict: FAIL**

### 2D Torus, 3×3 Blocking (81² → 9²)

| Level | Size | Growth | Spectral | Fisher Rank | Fisher PR (σ=3) | ΔPR | Status |
|-------|------|--------|----------|-------------|-----------------|-----|--------|
| 0 | 81² | 1.927 | 2.789 | 2.00 | 2.7586 | — | — |
| 1 | 27² | 1.729 | 2.792 | 2.00 | 2.7754 | +0.0168 | MARGINAL |
| 2 | 9² | N/A | 2.929 | 2.00 | 3.2160 | +0.4406 | FAIL |

**P7 Verdict: FAIL**

### 3D Torus, 2×2×2 Blocking (32³ → 4³)

| Level | Size | Growth | Spectral | Fisher Rank | Fisher PR (σ=3) | ΔPR | Status |
|-------|------|--------|----------|-------------|-----------------|-----|--------|
| 0 | 32³ | 2.769 | 3.910 | 3.00 | 4.1459 | — | — |
| 1 | 16³ | 2.508 | 3.901 | 3.00 | 4.3232 | +0.1773 | FAIL |
| 2 | 8³ | N/A | 3.816 | 3.00 | 5.0391 | +0.7159 | FAIL |
| 3 | 4³ | N/A | 3.709 | 1.00 | 5.9778 | +0.9387 | FAIL |

**P7 Verdict: FAIL**

## Overall P7 Assessment

**FAILURES in: 2D/2×2, 2D/3×3, 3D/2×2×2**

## Block-Spin Sanity Checks
- 2D 128→64 (2×2): ✓ Matches direct torus construction
- 2D 81→27 (3×3): ✓ Matches direct torus construction
- 3D 32→16 (2×2×2): ✓ Matches direct torus construction
