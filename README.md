# Fisher Diagnostic Suite

A multi-observable toolkit based on the Fisher Information Matrix that reads the information geometry of discrete systems. One construction, two regimes:

**Geometric regime** (exponential kernel): exact dimension estimation — integers on manifolds, Hausdorff dimension on fractals, DPI-monotone under coarse-graining.

**Thermal regime** (correlation kernel): structural probe — parameter-free detection and classification of phase transitions via spectral signatures in the FIM.

**Author**: Ian Darling | **License**: MIT

---

## Quick Start

```bash
# Phase 1: Dimension estimation on a 2D torus (< 2 min)
cd PHASE1_GEOMETRIC/src/
python fisher_phase1.py

# Phase 2: Ising model phase transition detection (< 10 min)
cd PHASE2_THERMAL/src/
python ising_fisher_phase2.py
```

---

## Results at a Glance

### Geometric Regime (Phase 1)
| System | True d | Fisher Rank | Error |
|--------|--------|------------|-------|
| 2D Torus | 2 | **2.000** | 0% |
| 3D Torus | 3 | **3.000** | 0% |
| 4D Torus | 4 | **4.000** | 0% |
| Sierpinski Carpet | 1.893 | PR → **1.905** | 0.6% |
| ER Random Graph | — | **1** (correct negative control) | — |

### Thermal Regime (Phase 2)
| System | Transition | SV₂/SV₁ at T_c | Swap? |
|--------|-----------|----------------|-------|
| 2D Ising | Continuous | 1.000 | YES |
| 3D Ising | Continuous | 0.956 | YES |
| Potts q=5 | 1st-order | 0.691 | NO |
| Potts q=10 | 1st-order | 0.374 | NO |

### Cross-Domain Transfer (Phase 3)
| Domain | Rank→1 at transition | η extremum | Systems |
|--------|---------------------|------------|---------|
| Lattice models | ✓ | ✓ | 4 models |
| Financial networks | 5/6 | 6/6 | ~90 S&P 500 stocks, 2005–2024 |
| EEG seizures | 5/5 | 4/5 | 5 patients, 26 seizures |

---

## Repository Structure

```
D_eff-Test-1.0/
├── README.md                    ← You are here
├── PHASE1_GEOMETRIC/            ← Dimension estimation (Paper 1)
│   ├── README.md
│   ├── DS_Phase1_MANUSCRIPT_FINAL.pdf
│   ├── src/
│   └── results/
├── PHASE2_THERMAL/              ← Phase transitions + cross-domain transfer (Paper 2)
│   ├── README.md
│   ├── DS_Phase2_Thermal_Paper_FINAL_v2.pdf
│   ├── src/
│   └── results/
└── docs/                        ← Project documentation
```

Each phase is self-contained with its own README, source code, and results. All experiments reproduce on consumer hardware (M4 MacBook Air, 24GB RAM) in under 20 minutes per phase.

---

## Papers

1. **Geometric regime**: *Effective Dimensionality from Fisher Information Rank: Operational Validation on Tori, Fractals, Random Graphs, and Coarse-Grained Lattices.* [PDF](PHASE1_GEOMETRIC/DS_Phase1_MANUSCRIPT_FINAL.pdf)

2. **Thermal regime**: *Spectral Phase Transitions in the Fisher Information Matrix: Universality of the Degeneracy Swap and Classification of Transition Order.* [PDF](PHASE2_THERMAL/DS_Phase2_Thermal_Paper_FINAL_v2.pdf)

---

## AI Disclosure

This research was conducted with AI-assisted development. Claude Opus (Anthropic) contributed to methodology design, analysis, and manuscript development. Claude Code (Anthropic) contributed to implementation. The author is responsible for all scientific claims, research direction, and editorial decisions.
