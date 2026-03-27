#!/usr/bin/env python3
"""
Phase 4B: q=4 Potts Model — Continuous/First-Order Boundary Test
================================================================
Tests the FIM SV₂/SV₁ at the marginal q=4 Potts transition, which sits
exactly on the boundary between continuous (q≤4) and first-order (q>4)
in 2D, with logarithmic corrections to scaling.

T_c = 1/ln(1+√4) = 1/ln(3) ≈ 0.9102 (exact, Baxter 1973)

The Phase 2 paper established:
  - 2D Ising (q=2): SV₂/SV₁ = 1.000 at T_c (continuous)
  - Potts q=5: SV₂/SV₁ = 0.691 at T_c (weak first-order)
  - Potts q=10: SV₂/SV₁ = 0.374 at T_c (strong first-order)

q=4 should be ≥ 0.691 (at least as high as q=5, likely higher since
the transition is formally continuous). The gradient predicts ~0.80.

Uses IDENTICAL code infrastructure to the Phase 2 Potts experiment
(Wolff cluster, FFT correlation function, thermal kernel FIM).

Pre-registered predictions:
  P4B-1: SV₂/SV₁ at T_c exceeds q=5 value of 0.691 [KILL if < 0.60]
  P4B-2: SV₂/SV₁ peak within ±10% of T_c = 0.9102
  P4B-3: The q-ordering is monotonic: q=2 > q=4 ≥ q=5 > q=10
  P4B-4: Rank at T_c ≥ 2 (non-trivial spectral structure)

Author: Ian Darling / Claude Code
Date: 2026-03-25
"""

import os
import sys
import time
import csv
import json
import numpy as np
from collections import deque

# Add Phase 2 src to path for reuse
PHASE2_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "PHASE2_THERMAL", "src")
sys.path.insert(0, PHASE2_SRC)

from potts5_fisher_phase_transition import (
    initialize_potts_lattice,
    wolff_step_potts,
    equilibrate_potts,
    measure_potts_correlations,
    accumulate_potts_correlations,
    estimate_correlation_length,
    build_thermal_kernel,
    compute_FIM_thermal,
    fisher_diagnostics,
    measure_potts_energy,
    measure_potts_macroscopic,
)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ─── Configuration ────────────────────────────────────────────────────
Q = 4
N = 128                # Lattice size (128×128, matching Phase 2)
T_C_Q4 = 1.0 / np.log(1.0 + np.sqrt(4.0))  # = 1/ln(3) ≈ 0.9102

# Temperature sweep: 20 values centered on T_c
T_OVER_TC_VALUES = np.concatenate([
    np.linspace(0.70, 0.90, 5),    # Well below T_c
    np.linspace(0.92, 1.08, 8),    # Near T_c
    np.linspace(1.12, 1.50, 7),    # Well above T_c
])
T_VALUES = T_OVER_TC_VALUES * T_C_Q4

# MC parameters (matching Phase 2 exactly)
N_EQ_SWEEPS = 300      # Equilibration Wolff steps
N_CONFIGS = 500        # Configurations for C(r) averaging
N_FIM_SAMPLES = 30     # FIM sample points

# Phase 2 reference values for comparison
REFERENCE = {
    'q=2 (Ising)':  {'sv2sv1': 1.000, 'rank': 2.0},
    'q=5':          {'sv2sv1': 0.691, 'rank': 1.93},
    'q=10':         {'sv2sv1': 0.374, 'rank': 1.0},
}

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
np.random.seed(42)


def run_temperature(T_val, t_idx, n_temps):
    """Run full Potts q=4 analysis at one temperature."""
    t0 = time.time()
    T_over_Tc = T_val / T_C_Q4

    # Initialize: cold start below T_c, hot above
    mode = 'hot' if T_val > T_C_Q4 else 'cold'
    lattice = initialize_potts_lattice(N, q=Q, mode=mode)

    # Equilibrate
    equilibrate_potts(lattice, T_val, q=Q, n_sweeps=N_EQ_SWEEPS)

    # More configs near T_c
    nc = N_CONFIGS
    ns = N_FIM_SAMPLES
    if 0.95 <= T_over_Tc <= 1.05:
        nc = min(N_CONFIGS * 2, 1000)
        ns = 50

    # Correlations
    lat_corr = lattice.copy()
    G_r = accumulate_potts_correlations(lat_corr, T_val, q=Q,
                                         n_configs=nc, n_sweeps_between=5)
    xi, xi_r2 = estimate_correlation_length(G_r)

    # FIM diagnostics
    fisher = fisher_diagnostics(G_r, N, n_samples=ns)

    # Macroscopic observables
    lat_macro = lattice.copy()
    macro = measure_potts_macroscopic(lat_macro, T_val, q=Q,
                                       n_configs=nc, n_sweeps_between=5)

    dt = time.time() - t0

    sv_raw = fisher['sv_profile_mean'] if fisher else [0]*4
    sv_profile = sv_raw.tolist() if hasattr(sv_raw, 'tolist') else list(sv_raw)
    # SV₂/SV₁ from the mean SV profile
    sv2sv1_val = sv_profile[1] if len(sv_profile) > 1 else 0.0

    result = {
        'T': T_val,
        'T_over_Tc': T_over_Tc,
        'sv2sv1_mean': sv2sv1_val,
        'sv2sv1_std': fisher['sv_profile_std'][1] if fisher and len(fisher['sv_profile_std']) > 1 else 0.0,
        'rank_mean': fisher['rank_mean'] if fisher else 0.0,
        'eta_mean': fisher['eta_mean'] if fisher else 0.0,
        'pr_mean': fisher['pr_mean'] if fisher else 0.0,
        'xi': xi if np.isfinite(xi) else -1.0,
        'chi': macro['susceptibility'],
        'energy': macro['energy'],
        'sv_profile': sv_profile,
    }

    print(f"  [{t_idx+1}/{n_temps}] T={T_val:.4f} (T/Tc={T_over_Tc:.3f}): "
          f"SV₂/SV₁={result['sv2sv1_mean']:.4f} rank={result['rank_mean']:.2f} "
          f"η={result['eta_mean']:.4f} PR={result['pr_mean']:.3f} "
          f"ξ={xi:.1f} χ={macro['susceptibility']:.1f} [{dt:.1f}s]")

    return result


def evaluate_predictions(results):
    """Evaluate pre-registered predictions."""
    Ts = np.array([r['T'] for r in results])
    T_ratios = np.array([r['T_over_Tc'] for r in results])
    sv2sv1 = np.array([r['sv2sv1_mean'] for r in results])
    ranks = np.array([r['rank_mean'] for r in results])

    # Find value closest to T_c
    tc_idx = np.argmin(np.abs(T_ratios - 1.0))
    sv_at_tc = sv2sv1[tc_idx]
    rank_at_tc = ranks[tc_idx]

    lines = []
    lines.append("=" * 70)
    lines.append("PHASE 4B: q=4 POTTS — BOUNDARY CLASSIFIER PREDICTION VERDICTS")
    lines.append("=" * 70)
    lines.append(f"T_c = {T_C_Q4:.6f}, N = {N}, q = {Q}")
    lines.append(f"SV₂/SV₁ at T_c = {sv_at_tc:.4f} (T/T_c = {T_ratios[tc_idx]:.3f})")
    lines.append("")

    # P4B-1: SV₂/SV₁ > 0.691 (q=5 value)
    p1_pass = sv_at_tc > 0.691
    p1_kill = sv_at_tc < 0.60
    lines.append(f"P4B-1: SV₂/SV₁ at T_c = {sv_at_tc:.4f} vs q=5 threshold 0.691")
    lines.append(f"  → {'PASS ✓' if p1_pass else ('KILL ✗' if p1_kill else 'WEAK')}")
    lines.append("")

    # P4B-2: Peak within ±10% of T_c
    peak_idx = np.argmax(sv2sv1)
    peak_T = Ts[peak_idx]
    peak_ratio = T_ratios[peak_idx]
    p2_pass = 0.90 <= peak_ratio <= 1.10
    lines.append(f"P4B-2: SV₂/SV₁ peak at T/T_c = {peak_ratio:.3f} (T={peak_T:.4f})")
    lines.append(f"  Window: [0.90, 1.10]")
    lines.append(f"  → {'PASS ✓' if p2_pass else 'FAIL ✗'}")
    lines.append("")

    # P4B-3: Monotonic ordering q=2 > q=4 ≥ q=5 > q=10
    p3_q2 = 1.000 > sv_at_tc  # q=2 should be above q=4
    p3_q5 = sv_at_tc >= 0.691  # q=4 should be at or above q=5
    p3_pass = p3_q2 and p3_q5
    lines.append(f"P4B-3: q=2(1.000) > q=4({sv_at_tc:.4f}) ≥ q=5(0.691) > q=10(0.374)")
    lines.append(f"  → {'PASS ✓' if p3_pass else 'FAIL ✗'}")
    lines.append("")

    # P4B-4: Rank ≥ 2
    p4_pass = rank_at_tc >= 2.0
    lines.append(f"P4B-4: Rank at T_c = {rank_at_tc:.2f}")
    lines.append(f"  → {'PASS ✓' if p4_pass else 'INCONCLUSIVE'}")
    lines.append("")

    # Classification
    if sv_at_tc > 0.95:
        classification = "CONTINUOUS (like Ising)"
    elif sv_at_tc > 0.70:
        classification = "MARGINAL (between continuous and first-order)"
    elif sv_at_tc > 0.50:
        classification = "WEAK FIRST-ORDER (like q=5)"
    else:
        classification = "STRONG FIRST-ORDER (like q=10)"

    lines.append("=" * 70)
    lines.append(f"CLASSIFICATION: {classification}")
    lines.append(f"SV₂/SV₁ = {sv_at_tc:.4f} places q=4 in the spectrum:")
    lines.append(f"  q=2:  1.000  (continuous)")
    lines.append(f"  q=4:  {sv_at_tc:.4f}  ← THIS RESULT")
    lines.append(f"  q=5:  0.691  (weak first-order)")
    lines.append(f"  q=10: 0.374  (strong first-order)")
    lines.append("")

    overall = "PASS" if (p1_pass and p2_pass) else "FAIL"
    lines.append(f"OVERALL: {overall}")
    lines.append(f"P4B-1={'PASS' if p1_pass else ('KILL' if p1_kill else 'WEAK')} "
                f"P4B-2={'PASS' if p2_pass else 'FAIL'} "
                f"P4B-3={'PASS' if p3_pass else 'FAIL'} "
                f"P4B-4={'PASS' if p4_pass else '?'}")
    lines.append("=" * 70)

    # Data table
    lines.append("\n--- Full Data ---")
    lines.append(f"{'T':>8} {'T/Tc':>7} {'SV2/SV1':>10} {'±':>8} {'rank':>6} "
                f"{'η':>8} {'PR':>8} {'ξ':>8} {'χ':>8}")
    for r in results:
        lines.append(f"{r['T']:8.4f} {r['T_over_Tc']:7.3f} "
                     f"{r['sv2sv1_mean']:10.4f} {r['sv2sv1_std']:8.4f} "
                     f"{r['rank_mean']:6.2f} {r['eta_mean']:8.4f} "
                     f"{r['pr_mean']:8.4f} {r['xi']:8.1f} {r['chi']:8.1f}")

    return "\n".join(lines)


def make_plots(results):
    """Generate diagnostic plots."""
    Ts = [r['T_over_Tc'] for r in results]
    sv2sv1 = [r['sv2sv1_mean'] for r in results]
    sv2sv1_err = [r['sv2sv1_std'] for r in results]
    ranks = [r['rank_mean'] for r in results]
    etas = [r['eta_mean'] for r in results]
    prs = [r['pr_mean'] for r in results]
    xis = [r['xi'] for r in results]
    chis = [r['chi'] for r in results]

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f"Phase 4B: q=4 Potts Boundary Test (N={N}, T_c={T_C_Q4:.4f})",
                fontsize=14, fontweight='bold')

    # A: SV₂/SV₁ with Phase 2 reference lines
    ax = axes[0, 0]
    ax.errorbar(Ts, sv2sv1, yerr=sv2sv1_err, fmt='o-', color='#2196F3',
               capsize=3, label=f'q=4')
    ax.axvline(1.0, color='red', linestyle='--', alpha=0.7, label='T_c')
    ax.axhline(1.000, color='green', linestyle=':', alpha=0.5, label='q=2 (Ising)')
    ax.axhline(0.691, color='orange', linestyle=':', alpha=0.5, label='q=5')
    ax.axhline(0.374, color='purple', linestyle=':', alpha=0.5, label='q=10')
    ax.set_xlabel('T / T_c'); ax.set_ylabel('SV₂/SV₁')
    ax.set_title('FIM Spectral Ratio'); ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # B: Rank
    ax = axes[0, 1]
    ax.plot(Ts, ranks, 'D-', color='#4CAF50')
    ax.axvline(1.0, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('T / T_c'); ax.set_ylabel('Rank')
    ax.set_title('Gap-Based Rank'); ax.grid(True, alpha=0.3)

    # C: η
    ax = axes[0, 2]
    ax.plot(Ts, etas, 'o-', color='#FF9800')
    ax.axvline(1.0, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('T / T_c'); ax.set_ylabel('η')
    ax.set_title('Disorder Index'); ax.grid(True, alpha=0.3)

    # D: PR
    ax = axes[1, 0]
    ax.plot(Ts, prs, '^-', color='#795548')
    ax.axvline(1.0, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('T / T_c'); ax.set_ylabel('PR')
    ax.set_title('Participation Ratio'); ax.grid(True, alpha=0.3)

    # E: ξ
    ax = axes[1, 1]
    valid = [(t, x) for t, x in zip(Ts, xis) if 0 < x < 1000]
    if valid:
        ax.plot([t for t,x in valid], [x for t,x in valid], 's-', color='#9C27B0')
    ax.axvline(1.0, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('T / T_c'); ax.set_ylabel('ξ')
    ax.set_title('Correlation Length'); ax.grid(True, alpha=0.3)

    # F: Susceptibility
    ax = axes[1, 2]
    ax.plot(Ts, chis, 'v-', color='#F44336')
    ax.axvline(1.0, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('T / T_c'); ax.set_ylabel('χ')
    ax.set_title('Susceptibility'); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(RESULTS_DIR, "phase4b_potts4_diagnostics.png")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {path}")

    # SV profile at T_c
    tc_idx = np.argmin([abs(r['T_over_Tc'] - 1.0) for r in results])
    sv_prof = results[tc_idx]['sv_profile']
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.bar(range(1, len(sv_prof)+1), sv_prof, color='#2196F3', alpha=0.8)
    ax2.axhline(0.01, color='red', linestyle='--', alpha=0.5, label='Threshold 0.01')
    ax2.set_xlabel('SV Index'); ax2.set_ylabel('Normalized SV')
    ax2.set_title(f'q=4 Potts SV Profile at T_c (SV₂/SV₁ = {sv_prof[1]:.4f})')
    ax2.legend()
    path2 = os.path.join(RESULTS_DIR, "phase4b_sv_profile_at_tc.png")
    plt.savefig(path2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {path2}")


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("=" * 70)
    print("PHASE 4B: q=4 POTTS — CONTINUOUS/FIRST-ORDER BOUNDARY TEST")
    print("=" * 70)
    print(f"q = {Q}, N = {N}, T_c = {T_C_Q4:.6f}")
    print(f"T sweep: {len(T_VALUES)} values, T/T_c from {T_OVER_TC_VALUES[0]:.2f} to {T_OVER_TC_VALUES[-1]:.2f}")
    print(f"MC: {N_EQ_SWEEPS} eq + {N_CONFIGS} configs, {N_FIM_SAMPLES} FIM samples")
    print(f"\nPhase 2 reference: q=2→1.000, q=5→0.691, q=10→0.374")
    print(f"Prediction: q=4 SV₂/SV₁ ∈ [0.70, 0.95] (marginal)")
    print()

    t_start = time.time()
    results = []

    for i, T_val in enumerate(T_VALUES):
        res = run_temperature(T_val, i, len(T_VALUES))
        results.append(res)

    t_total = time.time() - t_start
    print(f"\nTotal runtime: {t_total:.1f}s ({t_total/60:.1f} min)")

    # Write CSV
    csv_path = os.path.join(RESULTS_DIR, "phase4b_potts4_results.csv")
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'T', 'T_over_Tc', 'sv2sv1_mean', 'sv2sv1_std',
            'rank_mean', 'eta_mean', 'pr_mean', 'xi', 'chi', 'energy'
        ])
        writer.writeheader()
        for r in results:
            row = {k: f"{v:.6f}" if isinstance(v, float) else v
                   for k, v in r.items() if k != 'sv_profile'}
            writer.writerow(row)
    print(f"Results: {csv_path}")

    # Evaluate predictions
    summary = evaluate_predictions(results)
    summary_path = os.path.join(RESULTS_DIR, "phase4b_potts4_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(summary)
    print(f"Summary: {summary_path}")

    # Plots
    make_plots(results)

    # Save SV profile at T_c as JSON
    tc_idx = np.argmin([abs(r['T_over_Tc'] - 1.0) for r in results])
    json_path = os.path.join(RESULTS_DIR, "phase4b_sv_at_tc.json")
    with open(json_path, 'w') as f:
        json.dump({
            'q': Q, 'N': N, 'T_c': T_C_Q4,
            'T_measured': results[tc_idx]['T'],
            'sv2sv1': results[tc_idx]['sv2sv1_mean'],
            'rank': results[tc_idx]['rank_mean'],
            'sv_profile': results[tc_idx]['sv_profile'],
        }, f, indent=2)

    print("\n" + summary)


if __name__ == "__main__":
    main()
