#!/usr/bin/env python3
"""
DS Phase 3D-2: Potts q=5 First-Order Transition Test
=====================================================================
Tests whether the SV degeneracy swap occurs at a first-order phase transition.

q=5 Potts on 2D square lattice: first-order transition at T_c = 1/ln(1+√5).
q=10 Potts control: strong first-order at T_c = 1/ln(1+√10), ξ ≈ 5.

If no swap → toolkit is a phase transition CLASSIFIER (continuous vs first-order).
If swap → toolkit is a UNIVERSAL transition detector.

Adapted from ising_fisher_phase_transition.py (Phase 2 v2, Manhattan distance).
"""

import os
import sys
import json
import time
import datetime
from collections import deque, Counter

import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = "phase3d2_results"
RAW_DIR = os.path.join(OUTPUT_DIR, "raw_data")
os.makedirs(RAW_DIR, exist_ok=True)

# Critical temperatures (exact, Baxter 1973)
T_C_POTTS5 = 1.0 / np.log(1.0 + np.sqrt(5.0))     # ≈ 0.85153
T_C_POTTS10 = 1.0 / np.log(1.0 + np.sqrt(10.0))    # ≈ 0.70131
T_C_ISING = 2.0 / np.log(1 + np.sqrt(2))            # ≈ 2.26919 (reference)

# Phase 2 reference data (2D Ising, N=128, v2 Manhattan)
ISING_SV_PROFILES = {
    1.50: [1.000, 0.225, 0.225, 0.040],
    1.10: [1.000, 0.674, 0.674, 0.170],
    1.02: [1.000, 1.000, 0.344, 0.090],
    1.00: [1.000, 1.000, 0.542, 0.079],
    0.90: [1.000, 0.479, 0.479, 0.223],
}

np.random.seed(42)


# ============================================================
# MODULE 1: Potts Monte Carlo (Wolff cluster algorithm)
# ============================================================

def initialize_potts_lattice(N, q=5, mode='random'):
    """Returns N×N array of spins in {0, 1, ..., q-1}."""
    if mode == 'cold':
        return np.zeros((N, N), dtype=np.int8)
    elif mode in ('hot', 'random'):
        return np.random.randint(0, q, size=(N, N)).astype(np.int8)
    else:
        raise ValueError(f"Unknown mode: {mode}")


def wolff_step_potts(lattice, T, q=5, J=1.0):
    """
    Execute one Wolff cluster flip for q-state Potts model.

    P_add = 1 - exp(-J/T)  (NOT 2J/T as in Ising).
    Cluster: all same-state neighbors bonded with probability P_add.
    Flip: reassign entire cluster to a uniformly random OTHER state.
    """
    N = lattice.shape[0]
    p_add = 1.0 - np.exp(-J / T)

    # Choose random seed
    si = np.random.randint(N)
    sj = np.random.randint(N)
    seed_state = lattice[si, sj]

    # BFS cluster growth
    visited = np.zeros((N, N), dtype=np.bool_)
    visited[si, sj] = True
    queue = deque()
    queue.append((si, sj))
    cluster = [(si, sj)]

    while queue:
        ci, cj = queue.popleft()
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            ni, nj = (ci + di) % N, (cj + dj) % N
            if not visited[ni, nj] and lattice[ni, nj] == seed_state:
                if np.random.random() < p_add:
                    visited[ni, nj] = True
                    queue.append((ni, nj))
                    cluster.append((ni, nj))

    # Choose new state uniformly from {0,...,q-1} \ {seed_state}
    new_state = np.random.randint(0, q - 1)
    if new_state >= seed_state:
        new_state += 1
    new_state = np.int8(new_state)

    # Flip cluster
    for ci, cj in cluster:
        lattice[ci, cj] = new_state

    return lattice


def wolff_sweep_potts(lattice, T, q=5, n_steps=None):
    """Execute n_steps Wolff cluster flips. Default: N*N."""
    N = lattice.shape[0]
    if n_steps is None:
        n_steps = N * N
    for _ in range(n_steps):
        wolff_step_potts(lattice, T, q)
    return lattice


def equilibrate_potts(lattice, T, q=5, n_sweeps=300):
    """
    Equilibrate Potts lattice at temperature T.
    Higher n_sweeps than Ising (300 vs 200) for metastability near first-order T_c.
    """
    for _ in range(n_sweeps):
        wolff_step_potts(lattice, T, q)
    return lattice


# ============================================================
# MODULE 2: Potts Correlation Function Measurement
# ============================================================

def measure_potts_correlations(lattice, q=5, max_r=None):
    """
    Measure Potts correlation function G(r) = <δ(s₀, s_r)> - 1/q.

    FFT approach: for each state s, compute autocorrelation of indicator(lattice==s),
    then sum over all q states to get <δ(s₀, s_r)>.
    """
    N = lattice.shape[0]
    if max_r is None:
        max_r = N // 4

    N2 = N * N
    autocorr_total = np.zeros((N, N), dtype=np.float64)

    for s in range(q):
        indicator = (lattice == s).astype(np.float64)
        F = np.fft.fft2(indicator)
        power = F * np.conj(F)
        autocorr_s = np.real(np.fft.ifft2(power)) / N2
        autocorr_total += autocorr_s

    # autocorr_total[dx, dy] = <δ(s₀, s_{(dx,dy)})>
    # Radial average using Manhattan distance
    ix = np.arange(N)
    dx = np.minimum(ix, N - ix)
    DX, DY = np.meshgrid(dx, dx)
    dist_grid = DX + DY

    G_r = np.zeros(max_r + 1)
    for r in range(max_r + 1):
        mask = (dist_grid == r)
        if np.any(mask):
            G_r[r] = np.mean(autocorr_total[mask]) - 1.0 / q

    return G_r


def accumulate_potts_correlations(lattice, T, q=5, n_configs=500,
                                   n_sweeps_between=5):
    """
    Measure G(r) averaged over n_configs configurations.
    Between configs: n_sweeps_between Wolff sweeps.
    """
    N = lattice.shape[0]
    max_r = N // 4
    G_sum = np.zeros(max_r + 1)

    for i in range(n_configs):
        for _ in range(n_sweeps_between):
            wolff_step_potts(lattice, T, q)
        G_r = measure_potts_correlations(lattice, q=q, max_r=max_r)
        G_sum += G_r[:max_r + 1]

    return G_sum / n_configs


def estimate_correlation_length(G_r):
    """
    Fit G(r) ~ A * exp(-r/xi) for r in [1, len(G_r)//2].
    Returns (xi, R^2).
    """
    max_fit = len(G_r) // 2
    r_vals = np.arange(1, max_fit + 1)
    g_vals = G_r[1:max_fit + 1]

    pos = g_vals > 0
    if np.sum(pos) < 3:
        return float('inf'), 0.0

    r_pos = r_vals[pos]
    log_g = np.log(g_vals[pos])

    try:
        def exp_decay(r, log_A, inv_xi):
            return log_A - r * inv_xi

        popt, _ = curve_fit(exp_decay, r_pos, log_g, p0=[0.0, 0.1], maxfev=5000)
        inv_xi = popt[1]
        if inv_xi <= 0:
            return float('inf'), 0.0
        xi = 1.0 / inv_xi

        predicted = exp_decay(r_pos, *popt)
        ss_res = np.sum((log_g - predicted) ** 2)
        ss_tot = np.sum((log_g - np.mean(log_g)) ** 2)
        r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

        return xi, r_sq
    except Exception:
        return float('inf'), 0.0


# ============================================================
# MODULE 3: Thermal Fisher Kernel and FIM (verbatim Phase 2)
# ============================================================

def build_thermal_kernel(G_r, N, v0):
    """
    Build thermal probability distribution p_{v0}(u; T) using |G(r)|.
    v0 = (i0, j0) tuple. Returns flat array of shape (N*N,).
    Manhattan distance, periodic BC.
    """
    i0, j0 = v0
    max_r = len(G_r) - 1

    ix = np.arange(N)
    di = np.minimum(np.abs(ix - i0), N - np.abs(ix - i0))
    dj = np.minimum(np.abs(ix - j0), N - np.abs(ix - j0))
    DI, DJ = np.meshgrid(di, dj, indexing='ij')
    dist = DI + DJ

    weights = np.zeros((N, N), dtype=np.float64)
    for r in range(max_r + 1):
        mask = (dist == r)
        weights[mask] = np.abs(G_r[r])

    beyond = dist > max_r
    if np.any(beyond):
        weights[beyond] = np.abs(G_r[max_r])

    p = weights.flatten()
    total = np.sum(p)
    if total < 1e-30:
        p = np.ones(N * N) / (N * N)
    else:
        p = p / total

    return p


def compute_FIM_thermal(G_r, N, v0):
    """
    Compute 4×4 Fisher Information Matrix at vertex v0.
    v0 = (i0, j0) tuple. Manhattan distance, periodic BC.
    """
    i0, j0 = v0

    p_v0 = build_thermal_kernel(G_r, N, v0)
    log_p_v0 = np.log(p_v0 + 1e-30)

    neighbors = [
        ((i0 - 1) % N, j0),
        ((i0 + 1) % N, j0),
        (i0, (j0 - 1) % N),
        (i0, (j0 + 1) % N),
    ]
    k = len(neighbors)

    score_vectors = np.zeros((k, N * N))
    for j, w in enumerate(neighbors):
        p_w = build_thermal_kernel(G_r, N, w)
        log_p_w = np.log(p_w + 1e-30)
        score_vectors[j, :] = log_p_w - log_p_v0

    weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
    FIM = weighted_scores @ weighted_scores.T

    return FIM


def gap_based_rank(sv_norm):
    """Gap-based rank from normalized singular values."""
    if len(sv_norm) <= 1:
        return 1
    ratios = sv_norm[1:] / np.maximum(sv_norm[:-1], 1e-15)
    return int(np.argmin(ratios)) + 1


def participation_ratio(svs):
    """PR = (sum sv_i)^2 / sum(sv_i^2)"""
    s = np.sum(svs)
    s2 = np.sum(svs ** 2)
    if s2 < 1e-30:
        return 0.0
    return s * s / s2


def disorder_index(sv_norm, rank):
    """eta = sv[rank] / sv[rank-1]"""
    if rank >= len(sv_norm) or rank < 1:
        return 0.0
    return sv_norm[rank] / sv_norm[rank - 1] if sv_norm[rank - 1] > 1e-15 else 0.0


def sv_profile(FIM):
    """Returns normalized SV profile."""
    svs = np.linalg.svd(FIM, compute_uv=False)
    if svs[0] < 1e-30:
        return svs
    return svs / svs[0]


def fisher_diagnostics(G_r, N, n_samples=30):
    """
    Sample n_samples random vertices. Compute FIM, rank, PR, eta, sv_profile.
    Returns diagnostics dict.
    """
    ranks = []
    prs = []
    etas = []
    sv_profiles_list = []
    gap_ratios = []

    margin = max(2, N // 8)
    sampled = 0
    attempts = 0
    while sampled < n_samples and attempts < n_samples * 5:
        i = np.random.randint(margin, N - margin)
        j = np.random.randint(margin, N - margin)
        attempts += 1

        FIM = compute_FIM_thermal(G_r, N, (i, j))
        svs = np.linalg.svd(FIM, compute_uv=False)

        if svs[0] < 1e-30:
            continue

        sv_n = svs / svs[0]
        r = gap_based_rank(sv_n)
        pr = participation_ratio(svs)
        eta = disorder_index(sv_n, r)
        gr = svs[0] / svs[1] if len(svs) > 1 and svs[1] > 1e-30 else float('inf')

        ranks.append(r)
        prs.append(pr)
        etas.append(eta)
        sv_profiles_list.append(sv_n)
        gap_ratios.append(gr)
        sampled += 1

    if not ranks:
        return None

    sv_arr = np.array(sv_profiles_list)
    gr_arr = np.array([g for g in gap_ratios if g < 1e10])

    return {
        'rank_mean': float(np.mean(ranks)),
        'rank_std': float(np.std(ranks)),
        'rank_distribution': dict(Counter(ranks)),
        'pr_mean': float(np.mean(prs)),
        'pr_std': float(np.std(prs)),
        'eta_mean': float(np.mean(etas)),
        'eta_std': float(np.std(etas)),
        'sv_profile_mean': sv_arr.mean(axis=0).tolist(),
        'sv_profile_std': sv_arr.std(axis=0).tolist(),
        'gap_ratio_mean': float(np.mean(gr_arr)) if len(gr_arr) > 0 else float('inf'),
    }


# ============================================================
# MODULE 4: Temperature Sweep and Diagnostics
# ============================================================

def measure_potts_energy(lattice, J=1.0):
    """Energy per spin for Potts: -J * (same-state neighbor pairs) / N^2."""
    N = lattice.shape[0]
    # Count same-state neighbor pairs (right + down)
    same_right = np.sum(lattice == np.roll(lattice, -1, axis=1))
    same_down = np.sum(lattice == np.roll(lattice, -1, axis=0))
    return -J * (same_right + same_down) / (N * N)


def measure_potts_macroscopic(lattice, T, q=5, n_configs=500,
                               n_sweeps_between=5):
    """
    Measure macroscopic Potts observables over n_configs configurations.
    Order parameter: m = (q * max_fraction - 1) / (q - 1).
    """
    N = lattice.shape[0]
    N2 = N * N

    ms = []
    m2s = []
    Es = []
    E2s = []

    for _ in range(n_configs):
        for _ in range(n_sweeps_between):
            wolff_step_potts(lattice, T, q)

        # Order parameter
        counts = np.bincount(lattice.flatten(), minlength=q)
        max_frac = np.max(counts) / N2
        m_val = (q * max_frac - 1.0) / (q - 1.0)
        ms.append(m_val)
        m2s.append(m_val ** 2)

        # Energy
        e_val = measure_potts_energy(lattice)
        Es.append(e_val)
        E2s.append(e_val ** 2)

    m_mean = np.mean(ms)
    m2_mean = np.mean(m2s)
    E_mean = np.mean(Es)
    E2_mean = np.mean(E2s)

    chi = N2 * (m2_mean - m_mean ** 2)
    C = N2 * (E2_mean - E_mean ** 2) / (T * T)

    return {
        'magnetization': float(m_mean),
        'susceptibility': float(chi),
        'energy': float(E_mean),
        'specific_heat': float(C),
        'energy_values': [float(e) for e in Es],  # For bimodality check
    }


def temperature_sweep_potts(N, q, T_c, temperatures,
                              n_eq_sweeps=300, n_configs=500):
    """
    Full temperature sweep for Potts model.
    Returns list of result dicts, one per temperature.
    """
    results = []

    for t_idx, T in enumerate(temperatures):
        t0 = time.time()
        T_over_Tc = T / T_c

        mode = 'hot' if T > T_c else 'cold'
        lattice = initialize_potts_lattice(N, q=q, mode=mode)

        equilibrate_potts(lattice, T, q=q, n_sweeps=n_eq_sweeps)

        nc = n_configs
        ns = 30
        if 0.98 <= T_over_Tc <= 1.05:
            nc = min(n_configs * 2, 1000)
            ns = 50

        lat_corr = lattice.copy()
        G_r = accumulate_potts_correlations(lat_corr, T, q=q,
                                             n_configs=nc, n_sweeps_between=5)

        xi, xi_r2 = estimate_correlation_length(G_r)

        fisher = fisher_diagnostics(G_r, N, n_samples=ns)

        lat_macro = lattice.copy()
        macro = measure_potts_macroscopic(lat_macro, T, q=q,
                                           n_configs=nc, n_sweeps_between=5)

        dt = time.time() - t0

        result = {
            'T': float(T),
            'T_over_Tc': float(T_over_Tc),
            'N': N,
            'q': q,
            'G_r': G_r.tolist(),
            'xi': float(xi) if np.isfinite(xi) else -1.0,
            'xi_fit_R2': float(xi_r2),
            'time_sec': float(dt),
        }
        if fisher is not None:
            result.update({f'fisher_{k}': v for k, v in fisher.items()})
        result.update({f'macro_{k}': v for k, v in macro.items()})

        # Check energy bimodality near T_c
        if 0.99 <= T_over_Tc <= 1.01:
            e_vals = macro['energy_values']
            e_hist, e_bins = np.histogram(e_vals, bins=30)
            # Simple bimodality check: is there a valley between two peaks?
            peaks = []
            for i in range(1, len(e_hist) - 1):
                if e_hist[i] > e_hist[i-1] and e_hist[i] > e_hist[i+1]:
                    peaks.append(i)
            result['energy_bimodal'] = len(peaks) >= 2
        else:
            result['energy_bimodal'] = False

        # Don't store raw energy values in JSON (too large)
        del result['macro_energy_values']

        rank_str = f"{fisher['rank_mean']:.2f}" if fisher else "N/A"
        eta_str = f"{fisher['eta_mean']:.3f}" if fisher else "N/A"
        pr_str = f"{fisher['pr_mean']:.3f}" if fisher else "N/A"
        chi_str = f"{macro['susceptibility']:.1f}"

        print(f"  [{t_idx+1}/{len(temperatures)}] T={T:.4f} (T/Tc={T_over_Tc:.3f}): "
              f"rank={rank_str} eta={eta_str} PR={pr_str} chi={chi_str} "
              f"xi={xi:.2f} R2={xi_r2:.3f}  [{dt:.1f}s]")

        results.append(result)

    return results


# ============================================================
# MODULE 5: Plots and Results Report
# ============================================================

def plot_diagnostic_panel(results, N, q, tag=""):
    """Figure 1: 4-panel diagnostic plot."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Fisher Diagnostic Panel — Potts q={q}, N={N}", fontsize=14, fontweight='bold')

    T_ratios = [r['T_over_Tc'] for r in results]

    ax = axes[0, 0]
    ranks = [r.get('fisher_rank_mean', np.nan) for r in results]
    rank_stds = [r.get('fisher_rank_std', 0) for r in results]
    ax.errorbar(T_ratios, ranks, yerr=rank_stds, fmt='o-', color='tab:blue', capsize=3)
    ax.axhline(3, color='gray', ls='--', alpha=0.5, label='D_eff=3 (Ising)')
    ax.axhline(2, color='gray', ls=':', alpha=0.5, label='D_eff=2')
    ax.axvline(1.0, color='red', ls='--', alpha=0.5, label='T_c')
    ax.set_xlabel('T / T_c')
    ax.set_ylabel('D_eff (gap-based rank)')
    ax.set_title('(a) Effective Dimension')
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    etas = [r.get('fisher_eta_mean', np.nan) for r in results]
    eta_stds = [r.get('fisher_eta_std', 0) for r in results]
    ax.errorbar(T_ratios, etas, yerr=eta_stds, fmt='s-', color='tab:orange', capsize=3)
    ax.axhline(0.23, color='blue', ls='--', alpha=0.4, label='Torus η=0.23')
    ax.axhline(0.68, color='green', ls='--', alpha=0.4, label='RGG η=0.68')
    ax.axvline(1.0, color='red', ls='--', alpha=0.5)
    ax.set_xlabel('T / T_c')
    ax.set_ylabel('Disorder Index η')
    ax.set_title('(b) Disorder Index')
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    grs = [r.get('fisher_gap_ratio_mean', np.nan) for r in results]
    ax.plot(T_ratios, grs, 'D-', color='tab:green')
    ax.axhline(1.0, color='gray', ls='--', alpha=0.5, label='No gap')
    ax.axvline(1.0, color='red', ls='--', alpha=0.5)
    ax.set_xlabel('T / T_c')
    ax.set_ylabel('Gap Ratio (sv₁/sv₂)')
    ax.set_title('(c) Gap Ratio')
    ax.set_yscale('log')
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    chis = [r.get('macro_susceptibility', np.nan) for r in results]
    ax.plot(T_ratios, chis, '^-', color='tab:red')
    ax.axvline(1.0, color='red', ls='--', alpha=0.5, label='T_c')
    ax.set_xlabel('T / T_c')
    ax.set_ylabel('Susceptibility χ')
    ax.set_title('(d) Susceptibility')
    ax.legend(fontsize=8)

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f"diagnostic_panel_q{q}_N{N}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_sv_evolution(results, N, q):
    """Figure 2: SV profile evolution at representative temperatures."""
    target_ratios = [1.50, 1.10, 1.02, 1.00, 0.90]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle(f"SV Profile Evolution — Potts q={q}, N={N}", fontsize=14, fontweight='bold')
    axes_flat = axes.flatten()

    for idx, t_target in enumerate(target_ratios):
        ax = axes_flat[idx]
        best = min(results, key=lambda r: abs(r['T_over_Tc'] - t_target))
        sv_mean = np.array(best.get('fisher_sv_profile_mean', [0, 0, 0, 0]))
        sv_std = np.array(best.get('fisher_sv_profile_std', [0, 0, 0, 0]))
        rank = best.get('fisher_rank_mean', 0)
        eta = best.get('fisher_eta_mean', 0)

        if t_target >= 1.10:
            colors = ['tab:blue'] * len(sv_mean)
        elif t_target >= 1.00:
            colors = ['tab:purple'] * len(sv_mean)
        else:
            colors = ['tab:red'] * len(sv_mean)

        ax.bar(range(1, len(sv_mean) + 1), sv_mean, yerr=sv_std, capsize=4,
               color=colors, alpha=0.7, edgecolor='black')
        ax.set_xlabel('SV Index')
        ax.set_ylabel('Normalized SV')
        ax.set_title(f"T/Tc={best['T_over_Tc']:.3f}  rank={rank:.1f}  η={eta:.2f}")
        ax.set_ylim(0, 1.15)

    # Last panel: detail at T_c
    ax = axes_flat[5]
    tc_result = min(results, key=lambda r: abs(r['T_over_Tc'] - 1.0))
    sv_mean = np.array(tc_result.get('fisher_sv_profile_mean', [0, 0, 0, 0]))
    sv_std = np.array(tc_result.get('fisher_sv_profile_std', [0, 0, 0, 0]))
    ax.bar(range(1, len(sv_mean) + 1), sv_mean, yerr=sv_std, capsize=4,
           color='tab:purple', alpha=0.7, edgecolor='black')
    for i, (m, s) in enumerate(zip(sv_mean, sv_std)):
        ax.text(i + 1, m + s + 0.03, f"{m:.3f}", ha='center', va='bottom', fontsize=9)
    ax.set_xlabel('SV Index')
    ax.set_ylabel('Normalized SV')
    ax.set_title(f"T=T_c (detail): q={q}")
    ax.set_ylim(0, 1.15)

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f"sv_profile_evolution_q{q}_N{N}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_sv_comparison_three_way(results_q5, results_q10):
    """
    Figure 3: THE COMPARISON — Ising vs Potts q=5 vs Potts q=10 SV profiles.
    3 rows × 5 columns.
    """
    comparison_ratios = [1.50, 1.10, 1.02, 1.00, 0.90]
    fig, axes = plt.subplots(3, 5, figsize=(22, 12))
    fig.suptitle("SV Profile Comparison: Ising vs Potts q=5 vs Potts q=10 (all N=128)",
                 fontsize=14, fontweight='bold')

    row_labels = ['2D Ising\n(continuous)', 'Potts q=5\n(weak 1st-order)', 'Potts q=10\n(strong 1st-order)']
    row_colors = ['tab:blue', 'tab:orange', 'tab:red']

    for col, t_target in enumerate(comparison_ratios):
        # Row 0: Ising (hardcoded)
        ax = axes[0, col]
        sv_ising = ISING_SV_PROFILES.get(t_target, [0, 0, 0, 0])
        ax.bar(range(1, 5), sv_ising, color='tab:blue', alpha=0.7, edgecolor='black')
        ax.set_ylim(0, 1.15)
        ax.set_title(f"T/Tc={t_target:.2f}")
        if col == 0:
            ax.set_ylabel(row_labels[0])

        # Check degeneracy
        if abs(sv_ising[0] - sv_ising[1]) < 0.05:
            ax.annotate('SV1≈SV2', xy=(1.5, 1.05), ha='center', fontsize=7, color='red')
        elif abs(sv_ising[1] - sv_ising[2]) < 0.05:
            ax.annotate('SV2≈SV3', xy=(2.5, sv_ising[1] + 0.05), ha='center', fontsize=7, color='blue')

        # Row 1: Potts q=5
        ax = axes[1, col]
        if results_q5:
            best = min(results_q5, key=lambda r: abs(r['T_over_Tc'] - t_target))
            sv_q5 = np.array(best.get('fisher_sv_profile_mean', [0, 0, 0, 0]))
            sv_q5_std = np.array(best.get('fisher_sv_profile_std', [0, 0, 0, 0]))
            ax.bar(range(1, 5), sv_q5, yerr=sv_q5_std, capsize=3,
                   color='tab:orange', alpha=0.7, edgecolor='black')
            if abs(sv_q5[0] - sv_q5[1]) < 0.05:
                ax.annotate('SV1≈SV2', xy=(1.5, 1.05), ha='center', fontsize=7, color='red')
            elif len(sv_q5) >= 3 and abs(sv_q5[1] - sv_q5[2]) < 0.05:
                ax.annotate('SV2≈SV3', xy=(2.5, sv_q5[1] + 0.05), ha='center', fontsize=7, color='blue')
        ax.set_ylim(0, 1.15)
        if col == 0:
            ax.set_ylabel(row_labels[1])

        # Row 2: Potts q=10
        ax = axes[2, col]
        if results_q10:
            best10 = min(results_q10, key=lambda r: abs(r['T_over_Tc'] - t_target))
            sv_q10 = np.array(best10.get('fisher_sv_profile_mean', [0, 0, 0, 0]))
            sv_q10_std = np.array(best10.get('fisher_sv_profile_std', [0, 0, 0, 0]))
            ax.bar(range(1, 5), sv_q10, yerr=sv_q10_std, capsize=3,
                   color='tab:red', alpha=0.7, edgecolor='black')
            if abs(sv_q10[0] - sv_q10[1]) < 0.05:
                ax.annotate('SV1≈SV2', xy=(1.5, 1.05), ha='center', fontsize=7, color='red')
            elif len(sv_q10) >= 3 and abs(sv_q10[1] - sv_q10[2]) < 0.05:
                ax.annotate('SV2≈SV3', xy=(2.5, sv_q10[1] + 0.05), ha='center', fontsize=7, color='blue')
        ax.set_ylim(0, 1.15)
        ax.set_xlabel('SV Index')
        if col == 0:
            ax.set_ylabel(row_labels[2])

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, "sv_comparison_ising_vs_potts.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_leading_indicator(results, N, q, T_c):
    """Figure 4: Leading indicator comparison — eta vs chi."""
    T_ratios = [r['T_over_Tc'] for r in results]
    etas = [r.get('fisher_eta_mean', np.nan) for r in results]
    chis = [r.get('macro_susceptibility', np.nan) for r in results]

    high_T_etas = [e for e, t in zip(etas, T_ratios) if t >= 1.3 and not np.isnan(e)]
    if high_T_etas:
        eta_threshold = np.mean(high_T_etas) + 2 * np.std(high_T_etas)
    else:
        eta_threshold = 0.35

    valid_chis = [c for c in chis if not np.isnan(c)]
    chi_threshold = 0.1 * max(valid_chis) if valid_chis else 0

    T_fisher = None
    T_chi = None

    for t, e in sorted(zip(T_ratios, etas), reverse=True):
        if not np.isnan(e) and e > eta_threshold and T_fisher is None:
            T_fisher = t

    for t, c in sorted(zip(T_ratios, chis), reverse=True):
        if not np.isnan(c) and c > chi_threshold and T_chi is None:
            T_chi = t

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()

    ax1.plot(T_ratios, etas, 's-', color='tab:blue', label='Fisher η', markersize=6)
    ax2.plot(T_ratios, chis, '^-', color='tab:red', label='Susceptibility χ', markersize=6)

    ax1.axhline(eta_threshold, color='tab:blue', ls=':', alpha=0.4)
    ax2.axhline(chi_threshold, color='tab:red', ls=':', alpha=0.4)

    if T_fisher is not None:
        ax1.axvline(T_fisher, color='tab:blue', ls='--', alpha=0.6, label=f'η onset: T/Tc={T_fisher:.3f}')
    if T_chi is not None:
        ax1.axvline(T_chi, color='tab:red', ls='--', alpha=0.6, label=f'χ onset: T/Tc={T_chi:.3f}')

    ax1.axvline(1.0, color='black', ls='--', alpha=0.3)

    ax1.set_xlabel('T / T_c')
    ax1.set_ylabel('Disorder Index η', color='tab:blue')
    ax2.set_ylabel('Susceptibility χ', color='tab:red')

    delta = (T_fisher - T_chi) if T_fisher is not None and T_chi is not None else 0
    if delta > 0.01:
        verdict = f"Fisher LEADING by ΔT/Tc = {delta:.3f}"
    elif delta < -0.01:
        verdict = f"Fisher LAGGING by ΔT/Tc = {-delta:.3f}"
    else:
        verdict = "Simultaneous onset"

    ax1.set_title(f"Leading Indicator — Potts q={q}, N={N}: {verdict}")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=8)

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f"leading_indicator_q{q}_N{N}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")

    return T_fisher, T_chi, verdict


def plot_finite_size_scaling(all_results_q5):
    """Figure 5: Finite-size scaling — q=5 all N overlaid."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Finite-Size Scaling — Potts q=5", fontsize=14, fontweight='bold')

    colors = {64: 'tab:blue', 128: 'tab:orange', 256: 'tab:red'}

    for N_val, results in sorted(all_results_q5.items()):
        T_ratios = [r['T_over_Tc'] for r in results]
        ranks = [r.get('fisher_rank_mean', np.nan) for r in results]
        etas = [r.get('fisher_eta_mean', np.nan) for r in results]

        ax1.plot(T_ratios, ranks, 'o-', color=colors.get(N_val, 'black'),
                 label=f'N={N_val}', markersize=5)
        ax2.plot(T_ratios, etas, 's-', color=colors.get(N_val, 'black'),
                 label=f'N={N_val}', markersize=5)

    ax1.axvline(1.0, color='black', ls='--', alpha=0.3)
    ax1.axhline(3, color='gray', ls='--', alpha=0.3)
    ax1.set_xlabel('T / T_c')
    ax1.set_ylabel('D_eff (rank mean)')
    ax1.set_title('Effective Dimension')
    ax1.legend()

    ax2.axvline(1.0, color='black', ls='--', alpha=0.3)
    ax2.set_xlabel('T / T_c')
    ax2.set_ylabel('Disorder Index η')
    ax2.set_title('Disorder Index')
    ax2.legend()

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, "finite_size_scaling_q5.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_correlation_length_comparison(results_q5_128, results_q10_128):
    """Figure 6: Correlation length vs T/T_c for q=5 and q=10."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Correlation Length Comparison", fontsize=14, fontweight='bold')

    # ξ vs T/Tc
    for results, label, color in [
        (results_q5_128, 'Potts q=5 (N=128)', 'tab:orange'),
        (results_q10_128, 'Potts q=10 (N=128)', 'tab:red'),
    ]:
        if results:
            T_ratios = [r['T_over_Tc'] for r in results]
            xis = [r['xi'] if r['xi'] > 0 else np.nan for r in results]
            ax1.plot(T_ratios, xis, 'o-', color=color, label=label, markersize=5)

    ax1.axvline(1.0, color='black', ls='--', alpha=0.3)
    ax1.set_xlabel('T / T_c')
    ax1.set_ylabel('Correlation Length ξ')
    ax1.set_title('ξ vs Temperature')
    ax1.legend()

    # R² of exponential fit vs T/Tc
    for results, label, color in [
        (results_q5_128, 'Potts q=5', 'tab:orange'),
        (results_q10_128, 'Potts q=10', 'tab:red'),
    ]:
        if results:
            T_ratios = [r['T_over_Tc'] for r in results]
            r2s = [r['xi_fit_R2'] for r in results]
            ax2.plot(T_ratios, r2s, 's-', color=color, label=label, markersize=5)

    ax2.axvline(1.0, color='black', ls='--', alpha=0.3)
    ax2.axhline(0.90, color='gray', ls='--', alpha=0.5, label='R²=0.90 threshold')
    ax2.set_xlabel('T / T_c')
    ax2.set_ylabel('R² (exponential fit)')
    ax2.set_title('Exponential Fit Quality')
    ax2.legend()
    ax2.set_ylim(0.5, 1.05)

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, "correlation_length_comparison.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def classify_swap(sv_profile_mean, threshold=0.20):
    """
    Classify whether a degeneracy swap occurred.
    Swap = SV2 within threshold of SV1.
    """
    if len(sv_profile_mean) < 2:
        return "INCOMPLETE"
    sv1, sv2 = sv_profile_mean[0], sv_profile_mean[1]
    if abs(sv1 - sv2) / max(sv1, 1e-10) <= threshold:
        return "SWAP"
    else:
        return "NO SWAP"


def write_results_report(all_q5, results_q10, verdicts_q5, verdict_q10, total_time):
    """Write the full results markdown report for Phase 3D-2."""
    lines = []
    lines.append("# Phase 3D-2 Results: Potts q=5 First-Order Transition Test\n")
    lines.append(f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    # Section 1
    lines.append("## 1. System Parameters\n")
    lines.append(f"- T_c (q=5) = {T_C_POTTS5:.5f} (exact: 1/ln(1+√5))")
    lines.append(f"- T_c (q=10) = {T_C_POTTS10:.5f} (exact: 1/ln(1+√10))")
    lines.append(f"- T_c (Ising) = {T_C_ISING:.6f} (reference)")
    lines.append(f"- q=5 lattice sizes: {sorted(all_q5.keys())}")
    lines.append(f"- q=10 lattice size: 128")
    lines.append(f"- FIM dimensionality: 4×4 (same lattice as Ising)")
    lines.append(f"- Distance metric: 2D Manhattan")
    lines.append(f"- Total runtime: {total_time:.1f}s\n")

    # Section 2: Correlation Function
    lines.append("## 2. Correlation Function Behavior\n")
    for q_val, results_dict in [("q=5", all_q5), ("q=10", {128: results_q10} if results_q10 else {})]:
        for N_val in sorted(results_dict.keys()):
            results = results_dict[N_val]
            lines.append(f"### {q_val}, N = {N_val}\n")
            lines.append("| T/Tc | ξ | R² (exp fit) | Notes |")
            lines.append("|------|---|-------------|-------|")
            for r in results:
                xi_str = f"{r['xi']:.2f}" if r['xi'] > 0 else "∞"
                r2_str = f"{r['xi_fit_R2']:.3f}"
                note = ""
                if r['xi_fit_R2'] < 0.85:
                    note = "Power law onset"
                if r['xi'] < 0:
                    xi_str = "divergent"
                    note = "Power law / divergent"
                bimodal = "BIMODAL" if r.get('energy_bimodal', False) else ""
                if bimodal:
                    note = f"{note} {bimodal}".strip()
                lines.append(f"| {r['T_over_Tc']:.3f} | {xi_str} | {r2_str} | {note} |")
            lines.append("")

    # Section 3: Fisher Diagnostics
    lines.append("## 3. Fisher Diagnostics vs Temperature (q=5)\n")
    for N_val in sorted(all_q5.keys()):
        lines.append(f"### N = {N_val}\n")
        lines.append("| T/Tc | Rank (mean±std) | η (mean±std) | PR (mean) | Gap Ratio | χ | C |")
        lines.append("|------|----------------|-------------|-----------|-----------|---|---|")
        for r in all_q5[N_val]:
            rank_m = r.get('fisher_rank_mean', 0)
            rank_s = r.get('fisher_rank_std', 0)
            eta_m = r.get('fisher_eta_mean', 0)
            eta_s = r.get('fisher_eta_std', 0)
            pr_m = r.get('fisher_pr_mean', 0)
            gr = r.get('fisher_gap_ratio_mean', 0)
            chi = r.get('macro_susceptibility', 0)
            C_val = r.get('macro_specific_heat', 0)
            lines.append(f"| {r['T_over_Tc']:.3f} | {rank_m:.2f}±{rank_s:.2f} | "
                        f"{eta_m:.3f}±{eta_s:.3f} | {pr_m:.3f} | {gr:.2f} | "
                        f"{chi:.1f} | {C_val:.2f} |")
        lines.append("")

    # Section 4: SV Profile Evolution — PRIMARY RESULT
    lines.append("## 4. SV Profile Evolution — THE PRIMARY RESULT\n")

    # q=5 at N=128
    if 128 in all_q5:
        lines.append("### q=5 Potts (N=128)\n")
        for target in [1.50, 1.10, 1.02, 1.00, 0.90]:
            best = min(all_q5[128], key=lambda r: abs(r['T_over_Tc'] - target))
            sv = best.get('fisher_sv_profile_mean', [])
            sv_str = ", ".join([f"{v:.3f}" for v in sv])
            swap_status = classify_swap(sv)
            lines.append(f"- T/Tc={best['T_over_Tc']:.3f}: [{sv_str}]  → {swap_status}")
        lines.append("")

        tc_q5 = min(all_q5[128], key=lambda r: abs(r['T_over_Tc'] - 1.0))
        swap_q5 = classify_swap(tc_q5.get('fisher_sv_profile_mean', []))
        lines.append(f"**q=5 SWAP STATUS AT T_c: {swap_q5}**\n")

    # q=10 at N=128
    if results_q10:
        lines.append("### q=10 Potts (N=128)\n")
        for target in [1.50, 1.10, 1.02, 1.00, 0.90]:
            best = min(results_q10, key=lambda r: abs(r['T_over_Tc'] - target))
            sv = best.get('fisher_sv_profile_mean', [])
            sv_str = ", ".join([f"{v:.3f}" for v in sv])
            swap_status = classify_swap(sv)
            lines.append(f"- T/Tc={best['T_over_Tc']:.3f}: [{sv_str}]  → {swap_status}")
        lines.append("")

        tc_q10 = min(results_q10, key=lambda r: abs(r['T_over_Tc'] - 1.0))
        swap_q10 = classify_swap(tc_q10.get('fisher_sv_profile_mean', []))
        lines.append(f"**q=10 SWAP STATUS AT T_c: {swap_q10}**\n")

    # Section 5: Three-Way Comparison
    lines.append("## 5. Three-Way Comparison: Ising vs q=5 vs q=10\n")
    lines.append("SV profiles at T/T_c = 1.00:\n")
    lines.append(f"- Ising:  {ISING_SV_PROFILES[1.00]}")
    if 128 in all_q5:
        tc_q5 = min(all_q5[128], key=lambda r: abs(r['T_over_Tc'] - 1.0))
        sv_q5_tc = tc_q5.get('fisher_sv_profile_mean', [])
        lines.append(f"- q=5:   [{', '.join([f'{v:.3f}' for v in sv_q5_tc])}]")
    if results_q10:
        tc_q10 = min(results_q10, key=lambda r: abs(r['T_over_Tc'] - 1.0))
        sv_q10_tc = tc_q10.get('fisher_sv_profile_mean', [])
        lines.append(f"- q=10:  [{', '.join([f'{v:.3f}' for v in sv_q10_tc])}]")
    lines.append("")

    # Classification
    swap_ising = "SWAP"
    swap_q5_val = classify_swap(tc_q5.get('fisher_sv_profile_mean', [])) if 128 in all_q5 else "N/A"
    swap_q10_val = classify_swap(tc_q10.get('fisher_sv_profile_mean', [])) if results_q10 else "N/A"

    lines.append(f"| Model | Transition Order | Swap at T_c? |")
    lines.append(f"|-------|-----------------|-------------|")
    lines.append(f"| Ising | Continuous | {swap_ising} |")
    lines.append(f"| Potts q=5 | First-order (weak) | {swap_q5_val} |")
    lines.append(f"| Potts q=10 | First-order (strong) | {swap_q10_val} |")
    lines.append("")

    # Section 6: Leading Indicator
    lines.append("## 6. Leading Indicator Assessment\n")
    for N_val, (T_f, T_c_onset, verdict) in sorted(verdicts_q5.items()):
        t_f_str = f"{T_f:.3f}" if T_f is not None else "N/A"
        t_c_str = f"{T_c_onset:.3f}" if T_c_onset is not None else "N/A"
        lines.append(f"- q=5 N={N_val}: T_Fisher={t_f_str}, T_chi={t_c_str} → **{verdict}**")
    if verdict_q10:
        T_f, T_c_onset, verdict = verdict_q10
        t_f_str = f"{T_f:.3f}" if T_f is not None else "N/A"
        t_c_str = f"{T_c_onset:.3f}" if T_c_onset is not None else "N/A"
        lines.append(f"- q=10 N=128: T_Fisher={t_f_str}, T_chi={t_c_str} → **{verdict}**")
    lines.append("")

    # Section 7: Falsifiability
    lines.append("## 7. Falsifiability Assessment\n")
    lines.append("| ID | Prediction | Result |")
    lines.append("|---|---|---|")
    if 128 in all_q5:
        lines.append(f"| P3D2-1 | No swap at q=5 T_c | {swap_q5_val} |")
    if results_q10:
        lines.append(f"| P3D2-2 | No swap at q=10 T_c | {swap_q10_val} |")

    if 128 in all_q5:
        tc_r2 = min(all_q5[128], key=lambda r: abs(r['T_over_Tc'] - 1.0))['xi_fit_R2']
        p3 = "PASS" if tc_r2 > 0.90 else f"FAIL (R²={tc_r2:.3f})"
        lines.append(f"| P3D2-3 | G(r) exponential at T_c (q=5) | {p3} |")

        ranks = [r.get('fisher_rank_mean', 0) for r in all_q5[128]]
        all_3 = all(abs(r - 3.0) < 0.5 for r in ranks)
        p4 = "PASS" if all_3 else f"FAIL"
        lines.append(f"| P3D2-4 | Rank=3 at all T (q=5) | {p4} |")

        if 128 in verdicts_q5:
            _, _, v = verdicts_q5[128]
            p5 = "PASS" if "LEADING" in v else f"FAIL ({v})"
            lines.append(f"| P3D2-5 | Fisher leads susceptibility (q=5) | {p5} |")
    lines.append("")

    # Section 8: Conclusion
    lines.append("## 8. Conclusion: Classifier or Universal Detector?\n")
    if swap_q5_val == "SWAP" and swap_q10_val == "SWAP":
        lines.append("**(C) Swap at both q=5 and q=10 → UNIVERSAL DETECTOR**")
        lines.append("The toolkit responds to any structural change, regardless of transition order.")
    elif swap_q5_val == "SWAP" and swap_q10_val == "NO SWAP":
        lines.append("**(B) Swap at q=5 but not q=10 → SIZE-DEPENDENT**")
        lines.append("ξ/L is the control parameter. q=5 has ξ≈2500 >> L, so pseudocritical.")
    elif swap_q5_val == "NO SWAP" and swap_q10_val == "NO SWAP":
        lines.append("**(A) No swap at either → CLASSIFIER CONFIRMED**")
        lines.append("The toolkit distinguishes continuous from first-order transitions.")
    elif swap_q5_val == "NO SWAP" and swap_q10_val == "SWAP":
        lines.append("**(D) No swap at q=5 but swap at q=10 → UNEXPECTED**")
        lines.append("Investigate further.")
    else:
        lines.append("Classification inconclusive.")
    lines.append("")

    # Section 9: Runtime
    lines.append("## 9. Runtime Log\n")
    for N_val in sorted(all_q5.keys()):
        times = [r['time_sec'] for r in all_q5[N_val]]
        lines.append(f"- q=5 N={N_val}: {sum(times):.1f}s total ({np.mean(times):.1f}s per T point)")
    if results_q10:
        times = [r['time_sec'] for r in results_q10]
        lines.append(f"- q=10 N=128: {sum(times):.1f}s total ({np.mean(times):.1f}s per T point)")
    lines.append(f"- **Grand total: {total_time:.1f}s**")

    report_path = os.path.join(OUTPUT_DIR, "PHASE3D2_POTTS5_RESULTS.md")
    with open(report_path, 'w') as f:
        f.write('\n'.join(lines))
    print(f"  Saved {report_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("DS Phase 3D-2: Potts q=5 First-Order Transition Test")
    print("=" * 70)
    print(f"T_c (q=5) = {T_C_POTTS5:.5f}")
    print(f"T_c (q=10) = {T_C_POTTS10:.5f}")
    print(f"Output: {OUTPUT_DIR}/")
    print()

    # q=5 temperature grid (12 points)
    temperatures_q5 = T_C_POTTS5 * np.array([
        1.50, 1.30, 1.20, 1.15, 1.10, 1.07,
        1.04, 1.02, 1.01, 1.005, 1.000, 0.900,
    ])

    # q=10 temperature grid (9 points)
    temperatures_q10 = T_C_POTTS10 * np.array([
        1.50, 1.20, 1.10, 1.04, 1.02, 1.01, 1.005, 1.000, 0.900,
    ])

    lattice_sizes_q5 = [64, 128, 256]
    all_q5 = {}
    verdicts_q5 = {}
    t_total_start = time.time()

    # ---- q=5 sweeps ----
    for N_val in lattice_sizes_q5:
        print(f"\n{'=' * 70}")
        print(f"  Potts q=5, N={N_val}")
        print(f"{'=' * 70}")

        results = temperature_sweep_potts(N_val, q=5, T_c=T_C_POTTS5,
                                           temperatures=temperatures_q5,
                                           n_eq_sweeps=300, n_configs=500)
        all_q5[N_val] = results

        raw_path = os.path.join(RAW_DIR, f"results_q5_N{N_val}.json")
        with open(raw_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"  Raw data saved to {raw_path}")

    # ---- q=10 sweep (N=128 only) ----
    print(f"\n{'=' * 70}")
    print(f"  Potts q=10, N=128 (strong first-order control)")
    print(f"{'=' * 70}")

    results_q10 = temperature_sweep_potts(128, q=10, T_c=T_C_POTTS10,
                                            temperatures=temperatures_q10,
                                            n_eq_sweeps=300, n_configs=500)

    raw_path = os.path.join(RAW_DIR, "results_q10_N128.json")
    with open(raw_path, 'w') as f:
        json.dump(results_q10, f, indent=2, default=str)
    print(f"  Raw data saved to {raw_path}")

    total_time = time.time() - t_total_start

    # ---- PLOTS ----
    print(f"\n{'=' * 70}")
    print("  GENERATING PLOTS")
    print(f"{'=' * 70}")

    for N_val in lattice_sizes_q5:
        plot_diagnostic_panel(all_q5[N_val], N_val, q=5)
        plot_sv_evolution(all_q5[N_val], N_val, q=5)
        T_f, T_c_onset, verdict = plot_leading_indicator(all_q5[N_val], N_val, q=5, T_c=T_C_POTTS5)
        verdicts_q5[N_val] = (T_f, T_c_onset, verdict)

    # q=10 plots
    plot_diagnostic_panel(results_q10, 128, q=10)
    plot_sv_evolution(results_q10, 128, q=10)
    T_f10, T_c10, verdict10 = plot_leading_indicator(results_q10, 128, q=10, T_c=T_C_POTTS10)
    verdict_q10 = (T_f10, T_c10, verdict10)

    # Three-way comparison (KEY FIGURE)
    results_q5_128 = all_q5.get(128, [])
    plot_sv_comparison_three_way(results_q5_128, results_q10)

    # Finite-size scaling
    plot_finite_size_scaling(all_q5)

    # Correlation length comparison
    plot_correlation_length_comparison(results_q5_128, results_q10)

    # Results report
    write_results_report(all_q5, results_q10, verdicts_q5, verdict_q10, total_time)

    print(f"\n{'=' * 70}")
    print(f"  TOTAL TIME: {total_time:.1f}s")
    print(f"{'=' * 70}")

    # Print headline
    print(f"\n{'*' * 70}")
    if 128 in all_q5:
        tc_q5 = min(all_q5[128], key=lambda r: abs(r['T_over_Tc'] - 1.0))
        sv_q5 = tc_q5.get('fisher_sv_profile_mean', [])
        swap_q5 = classify_swap(sv_q5)
        print(f"  q=5  at T_c: SV=[{', '.join([f'{v:.3f}' for v in sv_q5])}] → {swap_q5}")

    if results_q10:
        tc_q10 = min(results_q10, key=lambda r: abs(r['T_over_Tc'] - 1.0))
        sv_q10 = tc_q10.get('fisher_sv_profile_mean', [])
        swap_q10 = classify_swap(sv_q10)
        print(f"  q=10 at T_c: SV=[{', '.join([f'{v:.3f}' for v in sv_q10])}] → {swap_q10}")

    print(f"\n  Ising at T_c: SV={ISING_SV_PROFILES[1.00]} → SWAP")

    # Classification
    if swap_q5 == "SWAP" and swap_q10 == "SWAP":
        print(f"\n  CONCLUSION: (C) UNIVERSAL DETECTOR — swap at all transition orders")
    elif swap_q5 == "SWAP" and swap_q10 == "NO SWAP":
        print(f"\n  CONCLUSION: (B) SIZE-DEPENDENT — ξ/L controls swap")
    elif swap_q5 == "NO SWAP" and swap_q10 == "NO SWAP":
        print(f"\n  CONCLUSION: (A) CLASSIFIER — distinguishes continuous vs first-order")
    elif swap_q5 == "NO SWAP" and swap_q10 == "SWAP":
        print(f"\n  CONCLUSION: (D) UNEXPECTED — investigate")
    print(f"{'*' * 70}")


if __name__ == '__main__':
    main()
