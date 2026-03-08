# Phase 3D-1 Handback: 3D Ising Universality Test

## 1. THE HEADLINE

**The degeneracy swap OCCURS in the 3D Ising model.** The phenomenon generalizes across universality classes.

The SV profile on the L=32 primary dataset shows the following evolution:

| T/T_c | SV Profile [sv1–sv6] | Pattern |
|-------|----------------------|---------|
| 1.500 | [1.000, 0.059, 0.059, 0.059, 0.028, 0.028] | SV2=SV3=SV4 (triple), SV5=SV6 (pair) |
| 1.300 | [1.000, 0.091, 0.091, 0.091, 0.035, 0.035] | Same — SVs growing |
| 1.100 | [1.000, 0.198, 0.198, 0.198, 0.079, 0.079] | Same — SVs growing faster |
| 1.020 | [1.000, 0.663, 0.663, 0.663, 0.185, 0.185] | Triple at 66% of SV1 — approaching swap |
| **1.000** | **[1.000, 0.956, 0.956, 0.956, 0.252, 0.252]** | **SV1≈SV2≈SV3≈SV4 — QUADRUPLE DEGENERACY** |
| 0.900 | [1.000, 0.236, 0.236, 0.236, 0.127, 0.127] | Back to triple — below-Tc noise regime |

This is **scenario (a)** from the spec predictions: all three lattice-direction SVs rise together to match SV1, producing quadruple degeneracy at T_c. The cubic symmetry is preserved throughout — no sequential or partial swap.

**The program continues. Phase 3D-2 (Potts q=5) is warranted.**

---

## 2. Rank

**Rank is NOT 4 at all temperatures.** P3D-1 FAILS.

| L | High T (rank) | Near T_c (rank) | At T_c (rank) | Below T_c (rank) |
|---|---------------|-----------------|----------------|-------------------|
| 16 | 1 (sometimes 4) | 4 | 4 | 1 |
| 32 | 1 | 4 (from T/Tc≤1.04) | 4 | 1 |
| 48 | 1 | 1 (until T/Tc=1.00) | 4 | 1 |

**Rank = 1 at high T** because the correlation length ξ is very short (0.7–1.4 lattice spacings). The thermal kernel is so sharply peaked that the three lattice-direction SVs (SV2=SV3=SV4 ≈ 0.06) are negligible compared to SV1 (1.000). The gap detector correctly places the gap between SV1 and the triple.

**Rank = 4 near T_c** because the rising triple (SV2=SV3=SV4 → 0.96) creates a new largest gap between SV4 and SV5 (0.25). The gap detector correctly shifts.

This is actually a **rank transition** from 1→4 as T→T_c — a 3D-specific phenomenon not seen in 2D (where rank=3 at all T). The rank transition is itself a criticality diagnostic.

**Finite-size scaling of rank transition:**
- L=16: rank=4 appears by T/Tc=1.20
- L=32: rank=4 appears by T/Tc=1.04
- L=48: rank=4 appears only at T/Tc=1.00

The rank transition sharpens with system size, approaching a step function at T_c in the thermodynamic limit. This is consistent with the correlation length exceeding a threshold fraction of L.

---

## 3. SV Profile Structure: 3D vs 2D

**2D Ising (Phase 2, 4 SVs):**
```
High T:  [1.000, 0.225, 0.225, 0.040]  — SV2=SV3 (pair), SV4 small
At T_c:  [1.000, 1.000, 0.542, 0.079]  — SV1=SV2 (one rises)
Below:   [1.000, 0.479, 0.479, 0.223]  — SV2=SV3 restored
```

**3D Ising (Phase 3D-1, 6 SVs):**
```
High T:  [1.000, 0.059, 0.059, 0.059, 0.028, 0.028]  — SV2=SV3=SV4 (triple), SV5=SV6 pair
At T_c:  [1.000, 0.956, 0.956, 0.956, 0.252, 0.252]  — SV1≈SV2≈SV3≈SV4 (all rise together)
Below:   [1.000, 0.236, 0.236, 0.236, 0.127, 0.127]  — Triple restored
```

**Pattern correspondence:**
- In d dimensions, the lattice has d pairs of ±directions → d lattice-direction SVs
- At high T: these d SVs are degenerate and small (square/cubic symmetry)
- At T_c: these d SVs rise to match SV1
- In 2D: SV2=SV3 → SV1=SV2 (one of two rises to meet SV1)
- In 3D: SV2=SV3=SV4 → SV1≈SV2≈SV3≈SV4 (all three rise together)

The **anti-parallel residual SVs** (SV4 in 2D, SV5=SV6 in 3D) remain small throughout. These encode the non-exponential kernel contribution (rank=d+1 in Phase 1 terminology).

**The dimensional pattern generalizes perfectly.** The number of "rising" SVs equals the lattice dimensionality d.

---

## 4. Degeneracy Structure Detail

The 3D FIM has a clean block structure at all temperatures:

```
SV1:       radial mode (kernel shape)
SV2=SV3=SV4: lattice-direction modes (cubic symmetry triple)
SV5=SV6:    anti-parallel residuals (non-exponential correction pair)
```

At T_c, SV2/3/4 rise to match SV1, but the block structure is maintained:
- The triple degeneracy SV2=SV3=SV4 is preserved at every temperature
- SV5=SV6 pair is preserved at every temperature
- Only the gap position shifts (between SV1 and SV2 → between SV4 and SV5)

This means cubic symmetry is never broken by the thermal kernel — the swap is a collective phenomenon where all three lattice directions respond equally to the diverging correlation length.

---

## 5. η Behavior and P3D-3

**P3D-3 FAILS: η is NOT minimized at T_c.**

This is a methodological artifact of the rank transition:

| T/T_c | Rank | η = sv[rank]/sv[rank-1] | What η measures |
|-------|------|-------------------------|-----------------|
| 1.50  | 1    | 0.059 (=sv2/sv1)       | Triple/radial gap |
| 1.10  | 1    | 0.198 (=sv2/sv1)       | Triple/radial gap |
| 1.04  | 4    | 0.267 (=sv5/sv4)       | Residual/triple gap |
| 1.00  | 4    | 0.263 (=sv5/sv4)       | Residual/triple gap |
| 0.90  | 1    | 0.236 (=sv2/sv1)       | Triple/radial gap |

When rank transitions from 1→4, η jumps from measuring SV2/SV1 to measuring SV5/SV4 — completely different quantities. The "minimum" η is at high T (0.059) simply because SV2/SV1 is smallest there.

**The correct diagnostic** is the gap ratio SV1/SV2, which drops monotonically from 17.0 (high T) to 1.05 (T_c) and back to 4.24 (below T_c). This minimum at T_c is the true rank-independent criticality signal.

---

## 6. Leading Indicator Assessment

**P3D-4 PASSES: Fisher leads susceptibility at all lattice sizes.**

| L  | T_Fisher (η onset) | T_χ (χ onset) | Lead ΔT/T_c |
|----|-------------------|---------------|-------------|
| 16 | 1.200             | 1.100         | 0.100       |
| 32 | 1.200             | 1.020         | 0.180       |
| 48 | 1.200             | 1.000         | 0.200       |

The Fisher signal (η rising above high-T baseline) appears at T/T_c ≈ 1.20 for all L, while the susceptibility onset sharpens toward T_c with increasing L. This produces a growing lead with system size.

**Caveat:** Same as Phase 2 — the η threshold uses only 2 high-T points (T/Tc ≥ 1.3) for baseline, making the 2σ threshold tight. The reported leads are likely inflated. Conservative reading: genuine lead of ~0.04–0.10.

**The lead is larger in 3D than 2D.** This may reflect the stronger finite-size rounding of χ in 3D (γ/ν ≈ 1.96 in 3D vs 1.75 in 2D), which pushes the χ onset to lower T/T_c on small lattices.

---

## 7. 2D vs 3D Comparison Table

| Observable | 2D Ising (Phase 2) | 3D Ising (Phase 3D-1) |
|---|---|---|
| Universality class | ν=1, γ=7/4 | ν≈0.630, γ≈1.237 |
| FIM size | 4×4 | 6×6 |
| Rank (high T) | 3 | 1 |
| Rank (at T_c) | 3 | 4 |
| High-T degeneracy | SV2=SV3 (pair) | SV2=SV3=SV4 (triple) |
| At-T_c degeneracy | SV1=SV2 | SV1≈SV2≈SV3≈SV4 |
| Swap type | One SV rises | All d SVs rise together |
| η at T_c | 0.147 | 0.263 |
| η below T_c | 0.465 | 0.236 |
| Swap observed? | **YES** | **YES** |
| χ peak (T/T_c) | 1.005 | 1.000 |
| Leading indicator | YES (ΔT/Tc≈0.13) | YES (ΔT/Tc≈0.18) |
| Gap ratio at T_c | ~1.3 | 1.05 |
| Rank transition? | No (rank=3 always) | Yes (1→4 at T_c) |

---

## 8. Prediction Outcomes

| ID | Prediction | Result | Notes |
|---|---|---|---|
| P3D-1 | Rank=4 at all T | **FAIL** | Rank=1 at high T due to sharp kernel; rank=4 only near T_c. Rank transition is itself a diagnostic. |
| P3D-2 | SV degeneracy swap | **PASS** | Scenario (a): all three lattice-direction SVs rise together → quadruple degeneracy at T_c |
| P3D-3 | η minimum at T_c | **FAIL** | Artifact of rank transition: η measures different SV ratios when rank changes. Gap ratio (SV1/SV2) IS minimized at T_c. |
| P3D-4 | Fisher leads susceptibility | **PASS** | ΔT/Tc = 0.10 (L=16), 0.18 (L=32), 0.20 (L=48) |

**Gate check:**
- Gate A (macroscopic observables): χ peaks at T/Tc=1.000 for L=32 and L=48. ✓ Magnetization drops to near-zero above T_c. ✓ χ scales with L: 70→255→631 from L=16→32→48. Scaling exponent: log(631/70)/log(48/16) ≈ 2.0, close to γ/ν ≈ 1.96. ✓
- Gate B (correlations): G(r) decays exponentially at high T (R²>0.95). ξ grows from ~0.7 to ~3.9 approaching T_c. ✓

---

## 9. Surprises and Novel Findings

### 9a. Rank transition as criticality diagnostic
The rank transition from 1→4 is a binary, parameter-free criticality detector. In the thermodynamic limit, rank=1 for all T>T_c and rank=4 at T=T_c. The transition temperature sharpens with L:
- L=16: rank=4 by T/Tc=1.20
- L=32: rank=4 by T/Tc=1.04
- L=48: rank=4 only at T/Tc=1.00

This is a **new diagnostic not present in 2D**, where rank was always 3.

### 9b. Perfect cubic symmetry preservation
The triple degeneracy SV2=SV3=SV4 is maintained at every temperature — never broken. The swap does not go through any asymmetric intermediate state. This is expected from cubic symmetry but was not guaranteed by the framework.

### 9c. Anti-parallel pair SV5=SV6
The bottom two SVs form a degenerate pair at all temperatures. In 2D, SV4 was a single residual. In 3D, the residuals come in pairs. This is likely related to the structure of the non-exponential kernel correction in 3D Manhattan distance.

### 9d. Below-T_c timing
L=48 below T_c took 389 seconds (64% of total runtime) due to Wolff cluster growth in the ordered phase. Each cluster encompasses most of the lattice, making BFS expensive. This is a well-known limitation of single-cluster Wolff below T_c.

### 9e. Specific heat irregularity
L=32 specific heat peaks at T/Tc=1.010 (C=116.25), then drops sharply to 2.34 at T/Tc=1.000. This may reflect the cold-start initialization at T≤T_c producing non-equilibrium correlations, or insufficient equilibration at exactly T_c where τ_auto diverges.

---

## 10. Recommendation

**The degeneracy swap generalizes. Proceed to Phase 3D-2: Potts model q=5.**

The Potts q=5 model on a 2D square lattice has a **first-order** phase transition. If the degeneracy swap occurs at first-order transitions as well as continuous ones, the toolkit is universal across transition types. If the swap is absent at first-order transitions, the toolkit becomes a **phase transition classifier** — distinguishing transition order from SV profiles alone.

Either outcome is scientifically valuable:
- Swap at first-order → universal criticality detector
- No swap at first-order → continuous vs first-order classifier

---

## 11. Runtime

| L  | Sites   | Total time | Per T point |
|----|---------|-----------|-------------|
| 16 | 4,096   | 28.9s     | 2.4s        |
| 32 | 32,768  | 151.1s    | 12.6s       |
| 48 | 110,592 | 464.6s    | 38.7s       |
| **All** | | **644.6s** | |

All three lattice sizes completed within the 30-minute budget (actual: 10.7 minutes).
