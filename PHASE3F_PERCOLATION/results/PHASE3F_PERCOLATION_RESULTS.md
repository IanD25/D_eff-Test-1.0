# Phase 3F: Percolation FDS Results
Generated: 2026-03-21 20:29:27
Sigma: 3.0 | Realizations: 5 | Centers/realization: 20

## Gate Verification

**Gate A (2D):** Giant fraction at p≈0.5000: ['0.712', '0.700', '0.647', '0.559'] → PASS

**Gate B (2D):** PR at threshold (d_f=1.896):
  - L=32: PR=1.750 (7.7% from d_f)
  - L=64: PR=1.749 (7.7% from d_f)
  - L=128: PR=1.797 (5.2% from d_f)
  - L=256: PR=1.779 (6.1% from d_f)

**Gate A (3D):** Giant fraction at p≈0.2500: ['0.199', '0.241', '0.177'] → FAIL

**Gate B (3D):** PR at threshold (d_f=2.523):
  - L=16: PR=1.566 (37.9% from d_f)
  - L=32: PR=1.666 (34.0% from d_f)
  - L=64: PR=1.598 (36.6% from d_f)

## Primary Result: Kill Condition Assessment

### 2D (p_c = 0.5, measured at p = 0.5000)

| L | SV₂/SV₁ | PR | Rank | η | Verdict |
|---|---------|-----|------|---|--------|
| 32 | 0.3463 | 1.750 | 1.3 | 0.2782 | ❌ KILL |
| 64 | 0.3344 | 1.749 | 1.4 | 0.2478 | ❌ KILL |
| 128 | 0.3580 | 1.797 | 1.4 | 0.2626 | ❌ KILL |
| 256 | 0.3458 | 1.779 | 1.4 | 0.2399 | ❌ KILL |

### 3D (p_c = 0.24881, measured at p = 0.2500)

| L | SV₂/SV₁ | PR | Rank | η | Verdict |
|---|---------|-----|------|---|--------|
| 16 | 0.2545 | 1.566 | 1.2 | 0.2086 | ❌ KILL |
| 32 | 0.2904 | 1.666 | 1.4 | 0.2146 | ❌ KILL |
| 64 | 0.2752 | 1.598 | 1.3 | 0.2024 | ❌ KILL |

## Overall Verdict

**🔴 KILL SHOT TRIGGERED:** SV₂/SV₁ < 0.70 across all lattice sizes in both dimensions.

The two-dimensional ontology is **ABANDONED**. The SV swap is a thermal artifact, not a universal signature of dimensional coupling transitions.

## Pre-Registered Prediction Outcomes

**P3F-1 (PRIMARY):** Mean SV₂/SV₁ at p_c = 0.3150 → FAIL

**P3F-3 (2D):** SV₂/SV₁ peak at p=0.5750 (p_c=0.5000) → FAIL

**P3F-3 (3D):** SV₂/SV₁ peak at p=0.2080 (p_c=0.2488) → FAIL

**P3F-4 (2D):** SV₂/SV₁ at p_c vs L: [(32, '0.3463'), (64, '0.3344'), (128, '0.3580'), (256, '0.3458')] → FAIL

**P3F-4 (3D):** SV₂/SV₁ at p_c vs L: [(16, '0.2545'), (32, '0.2904'), (64, '0.2752')] → FAIL

## Next Step Recommendation

- FDS thermal regime is domain-specific to correlation-based transitions.
- Paper 2 stands. Broader ontological claim does not.
- Next: submit papers as-is.
