#!/usr/bin/env python3
"""
Q1: Scale-Dependent Information Sum Rule
=========================================
D_eff(r_max, T) on 2D Ising with Windowed Thermal Kernel

Tests whether scale-integrated FIM participation ratio I_PR(T)
has thermodynamic meaning by comparing to Onsager exact entropy
S(T) and specific heat C(T).

Builds on Phase 2 thermal FDS pipeline (ising_fisher_phase_transition.py).

Key modification: windowed kernel restricts observation radius to r_max,
allowing measurement of information accumulation across scales.
"""

import os
import sys
import json
import time
import datetime
from collections import deque, Counter
from multiprocessing import Pool, cpu_count

import numpy as np
from scipy.special import ellipk
from scipy.integrate import quad
from scipy.optimize import curve_fit
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# CONSTANTS AND PARAMETERS
# ============================================================

T_C = 2.0 / np.log(1 + np.sqrt(2))   # 2.26919...

T_OVER_TC = [
    0.70, 0.80, 0.85, 0.90, 0.93, 0.95, 0.97,
    0.99, 1.00, 1.01, 1.03, 1.05, 1.10,
    1.20, 1.30, 1.40, 1.60, 2.00
]

TEMPERATURES = [t * T_C for t in T_OVER_TC]

R_MAX_VALUES = [2, 3, 4, 5, 7, 10, 15, 20, 25, 30, 40, 50, 55, 64]

N_LATTICE = 128
N_EQ_SWEEPS = 500
N_CONFIGS = 500
N_SWEEPS_BETWEEN = 5
N_FIM_CENTERS = 20

OUTPUT_DIR = "q1_results"
RAW_DIR = os.path.join(OUTPUT_DIR, "raw_data")
os.makedirs(RAW_DIR, exist_ok=True)

np.random.seed(42)

# ============================================================
# MODULE 1: Ising Monte Carlo (Wolff cluster algorithm)
# Copied from Phase 2 ising_fisher_phase_transition.py
# ============================================================

def initialize_lattice(N, mode='random'):
    """Returns NxN array of spins in {+1, -1}."""
    if mode == 'cold':
        return np.ones((N, N), dtype=np.int8)
    elif mode == 'hot' or mode == 'random':
        return np.where(np.random.random((N, N)) < 0.5, 1, -1).astype(np.int8)
    else:
        raise ValueError(f"Unknown mode: {mode}")


def wolff_step(lattice, T, J=1.0):
    """Execute one Wolff cluster flip."""
    N = lattice.shape[0]
    p_add = 1.0 - np.exp(-2.0 * J / T)

    si = np.random.randint(N)
    sj = np.random.randint(N)
    seed_spin = lattice[si, sj]

    visited = np.zeros((N, N), dtype=np.bool_)
    visited[si, sj] = True
    queue = deque()
    queue.append((si, sj))
    cluster = [(si, sj)]

    while queue:
        ci, cj = queue.popleft()
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            ni, nj = (ci + di) % N, (cj + dj) % N
            if not visited[ni, nj] and lattice[ni, nj] == seed_spin:
                if np.random.random() < p_add:
                    visited[ni, nj] = True
                    queue.append((ni, nj))
                    cluster.append((ni, nj))

    for ci, cj in cluster:
        lattice[ci, cj] = -lattice[ci, cj]

    return lattice


def wolff_sweep(lattice, T, n_steps):
    """Execute n_steps Wolff cluster flips."""
    for _ in range(n_steps):
        wolff_step(lattice, T)
    return lattice


# ============================================================
# MODULE 2: Correlation Function Measurement
# Copied from Phase 2
# ============================================================

def measure_correlations_fft(lattice):
    """
    Measure connected G(r) using FFT-based autocorrelation.
    Returns G[r] for r = 0, 1, ..., N//2 using Manhattan distance.
    """
    N = lattice.shape[0]
    lat_float = lattice.astype(np.float64)
    m = np.mean(lat_float)
    m_sq = m * m

    F = np.fft.fft2(lat_float)
    power = F * np.conj(F)
    G_full = np.real(np.fft.ifft2(power)) / (N * N)

    max_r = N // 2
    G_r = np.zeros(max_r + 1)

    ix = np.arange(N)
    dx = np.minimum(ix, N - ix)
    DX, DY = np.meshgrid(dx, dx)
    dist_grid = DX + DY  # Manhattan distance

    for r in range(max_r + 1):
        mask = (dist_grid == r)
        if np.any(mask):
            G_r[r] = np.mean(G_full[mask]) - m_sq

    return G_r


def accumulate_correlations(lattice, T, n_configs=500, n_sweeps_between=5):
    """
    Measure G(r) averaged over n_configs independent configurations.
    Returns averaged G(r) array.
    """
    N = lattice.shape[0]
    max_r = N // 2
    G_sum = np.zeros(max_r + 1)

    for i in range(n_configs):
        for _ in range(n_sweeps_between):
            wolff_step(lattice, T)
        G_r = measure_correlations_fft(lattice)
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
# MODULE 3: Windowed Thermal Kernel and FIM Construction
# Modified from Phase 2 to add r_max windowing
# ============================================================

def build_windowed_kernel(G_r, N, v0, r_max):
    """
    Build probability distribution p_{v0}(u) using thermal kernel
    windowed at radius r_max.

    For each vertex u on NxN periodic lattice:
        r = Manhattan distance from v0 to u
        if 0 < r <= r_max:
            weight(u) = |G_r[r]|
        else:
            weight(u) = 0

    G_r[0] (self-correlation) is excluded (r > 0 condition).
    Normalise so sum(weights) = 1.

    Returns flat array of shape (N*N,) or None if degenerate.
    """
    i0, j0 = v0
    max_r_available = len(G_r) - 1

    ix = np.arange(N)
    di = np.minimum(np.abs(ix - i0), N - np.abs(ix - i0))
    dj = np.minimum(np.abs(ix - j0), N - np.abs(ix - j0))
    DI, DJ = np.meshgrid(di, dj, indexing='ij')
    dist = DI + DJ  # Manhattan distance

    # Build weights: |G(r)| for 0 < r <= r_max
    weights = np.zeros((N, N), dtype=np.float64)
    effective_r_max = min(r_max, max_r_available)
    for r in range(1, effective_r_max + 1):
        mask = (dist == r)
        weights[mask] = np.abs(G_r[r])

    p = weights.flatten()
    total = np.sum(p)
    if total < 1e-30:
        return None
    p = p / total

    return p


def compute_windowed_fim(G_r, N, v0, r_max):
    """
    Build 4x4 FIM at centre v0 using windowed kernel.

    1. Build p_{v0} using build_windowed_kernel
    2. Build p_{wj} for each of 4 Manhattan neighbours of v0
       (using same G_r and r_max, but centred on wj)
    3. Score vectors: s_j(u) = log p_{wj}(u) - log p_{v0}(u)
       Only for vertices u where BOTH p_{v0}(u) > 0 AND p_{wj}(u) > 0
    4. FIM: F_ij = sum_u s_i(u) * s_j(u) * p_{v0}(u)
    5. SVD -> singular values

    Returns: dict with sv_profile, rank, pr, eta, sv2_sv1
    Or None if kernel construction failed.
    """
    i0, j0 = v0

    p_v0 = build_windowed_kernel(G_r, N, v0, r_max)
    if p_v0 is None:
        return None

    neighbors = [
        ((i0 - 1) % N, j0),
        ((i0 + 1) % N, j0),
        (i0, (j0 - 1) % N),
        (i0, (j0 + 1) % N),
    ]
    k = len(neighbors)

    score_vectors = np.zeros((k, N * N))
    valid_mask = p_v0 > 1e-30

    for j, w in enumerate(neighbors):
        p_w = build_windowed_kernel(G_r, N, w, r_max)
        if p_w is None:
            return None

        # Only compute scores where both distributions have support
        both_valid = valid_mask & (p_w > 1e-30)
        if np.sum(both_valid) < 4:
            return None

        log_p_v0 = np.zeros(N * N)
        log_p_w = np.zeros(N * N)
        log_p_v0[both_valid] = np.log(p_v0[both_valid])
        log_p_w[both_valid] = np.log(p_w[both_valid])
        score_vectors[j, :] = log_p_w - log_p_v0

    # FIM: weighted by sqrt(p_v0) for numerical stability
    weighted_scores = score_vectors * np.sqrt(p_v0)[np.newaxis, :]
    FIM = weighted_scores @ weighted_scores.T

    # SVD
    svs = np.linalg.svd(FIM, compute_uv=False)
    if svs[0] < 1e-30:
        return None

    sv_norm = svs / svs[0]

    # Diagnostics
    rank = gap_based_rank(sv_norm)
    pr = participation_ratio(svs)
    eta = disorder_index(sv_norm, rank)
    sv2_sv1 = float(sv_norm[1]) if len(sv_norm) > 1 else 0.0

    return {
        'sv_profile': sv_norm.tolist(),
        'rank': rank,
        'pr': pr,
        'eta': eta,
        'sv2_sv1': sv2_sv1,
    }


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
    """eta = sv[rank] / sv[rank-1] -- ratio of (d+1)-th to d-th SV."""
    if rank >= len(sv_norm) or rank < 1:
        return 0.0
    return sv_norm[rank] / sv_norm[rank - 1] if sv_norm[rank - 1] > 1e-15 else 0.0


def run_windowed_sweep(G_r, N, r_max, n_centers=20):
    """
    Average FIM diagnostics over n_centers random vertices.
    Returns: dict with mean and std of each diagnostic, or None if all failed.
    """
    ranks = []
    prs = []
    etas = []
    sv2sv1s = []
    sv_profiles = []

    margin = max(2, N // 8)
    attempts = 0
    while len(ranks) < n_centers and attempts < n_centers * 5:
        i = np.random.randint(margin, N - margin)
        j = np.random.randint(margin, N - margin)
        attempts += 1

        result = compute_windowed_fim(G_r, N, (i, j), r_max)
        if result is None:
            continue

        ranks.append(result['rank'])
        prs.append(result['pr'])
        etas.append(result['eta'])
        sv2sv1s.append(result['sv2_sv1'])
        sv_profiles.append(result['sv_profile'])

    if not ranks:
        return None

    return {
        'rank_mean': float(np.mean(ranks)),
        'rank_std': float(np.std(ranks)),
        'pr_mean': float(np.mean(prs)),
        'pr_std': float(np.std(prs)),
        'eta_mean': float(np.mean(etas)),
        'eta_std': float(np.std(etas)),
        'sv2sv1_mean': float(np.mean(sv2sv1s)),
        'sv2sv1_std': float(np.std(sv2sv1s)),
        'sv_profile_mean': np.mean(sv_profiles, axis=0).tolist(),
        'sv_profile_std': np.std(sv_profiles, axis=0).tolist(),
        'n_valid': len(ranks),
    }


# ============================================================
# MODULE 4: Onsager Exact Thermodynamics
# ============================================================

def onsager_energy(T, J=1.0):
    """Exact internal energy per site, 2D Ising square lattice."""
    K = J / T
    kappa = 2 * np.sinh(2 * K) / np.cosh(2 * K) ** 2
    kappa_sq = min(kappa ** 2, 1.0 - 1e-12)
    K1 = ellipk(kappa_sq)
    return -J * (1 / np.tanh(2 * K)) * (1 + (2 / np.pi) * (2 * np.tanh(2 * K) ** 2 - 1) * K1)


def onsager_free_energy(T, J=1.0):
    """Exact free energy per site, 2D Ising square lattice."""
    K = J / T
    kappa = 2 * np.sinh(2 * K) / np.cosh(2 * K) ** 2

    def integrand(phi):
        val = 1 - kappa ** 2 * np.sin(phi) ** 2
        val = max(val, 1e-30)
        return np.log(0.5 * (1 + np.sqrt(val)))

    integral, _ = quad(integrand, 0, np.pi)
    return -T * (np.log(2) + integral / (2 * np.pi) + 0.5 * np.log(np.cosh(2 * K)))


def onsager_entropy(T, J=1.0):
    """Exact entropy per site."""
    return (onsager_energy(T, J) - onsager_free_energy(T, J)) / T


def onsager_specific_heat(T, J=1.0, dT=0.001):
    """Specific heat from numerical derivative of energy."""
    return (onsager_energy(T + dT, J) - onsager_energy(T - dT, J)) / (2 * dT)


# ============================================================
# MODULE 5: Sum Rule Computation
# ============================================================

def compute_sum_rule(all_fim_results, r_max_values):
    """
    For each T, integrate PR(T, r_max) and D_eff(T, r_max) over r_max
    using trapezoidal rule.

    all_fim_results: dict keyed by (T_idx, r_max) -> sweep result dict
    r_max_values: sorted list of r_max values

    Returns: I_PR[T_idx], I_rank[T_idx]
    """
    # Collect unique T indices
    t_indices = sorted(set(k[0] for k in all_fim_results.keys()))

    I_PR = {}
    I_rank = {}

    for t_idx in t_indices:
        pr_values = []
        rank_values = []
        r_values = []

        for r_max in r_max_values:
            key = (t_idx, r_max)
            if key in all_fim_results and all_fim_results[key] is not None:
                pr_values.append(all_fim_results[key]['pr_mean'])
                rank_values.append(all_fim_results[key]['rank_mean'])
                r_values.append(r_max)

        if len(r_values) >= 2:
            _trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
            I_PR[t_idx] = float(_trapz(pr_values, r_values))
            I_rank[t_idx] = float(_trapz(rank_values, r_values))
        else:
            I_PR[t_idx] = 0.0
            I_rank[t_idx] = 0.0

    return I_PR, I_rank


# ============================================================
# MODULE 6: Phase A — Monte Carlo
# ============================================================

def run_phase_a(temperatures, N=128, n_eq=500, n_configs=500):
    """
    Phase A: MC at each temperature -> G(r, T).
    Returns dict of T_idx -> G_r array.
    """
    print("\n" + "=" * 70)
    print("  PHASE A: Monte Carlo — Correlation Functions")
    print("=" * 70)

    correlation_data = {}

    for t_idx, T in enumerate(temperatures):
        t0 = time.time()
        T_over_Tc = T / T_C

        # Hot start above Tc, cold start below
        mode = 'hot' if T > T_C else 'cold'
        lattice = initialize_lattice(N, mode=mode)

        # Equilibrate
        for _ in range(n_eq):
            wolff_step(lattice, T)

        # Accumulate correlations
        G_r = accumulate_correlations(lattice, T, n_configs=n_configs,
                                      n_sweeps_between=N_SWEEPS_BETWEEN)

        # Correlation length
        xi, xi_r2 = estimate_correlation_length(G_r)

        dt = time.time() - t0
        xi_str = f"{xi:.2f}" if np.isfinite(xi) and xi > 0 else "div"
        print(f"  [{t_idx+1}/{len(temperatures)}] T/Tc={T_over_Tc:.3f}: "
              f"xi={xi_str}  G(1)={G_r[1]:.4f}  [{dt:.1f}s]")

        correlation_data[t_idx] = {
            'T': float(T),
            'T_over_Tc': float(T_over_Tc),
            'G_r': G_r.tolist(),
            'xi': float(xi) if np.isfinite(xi) else -1.0,
            'xi_fit_R2': float(xi_r2),
        }

    return correlation_data


# ============================================================
# MODULE 7: Phase B — Windowed FIM Sweep
# ============================================================

def run_phase_b(correlation_data, r_max_values, N=128, n_centers=20):
    """
    Phase B: Windowed FIM at each (T, r_max) pair.
    Returns dict of (T_idx, r_max) -> sweep result dict.
    """
    print("\n" + "=" * 70)
    print("  PHASE B: Windowed FIM Sweep")
    print("=" * 70)

    all_fim_results = {}
    total_pairs = len(correlation_data) * len(r_max_values)
    pair_count = 0

    for t_idx in sorted(correlation_data.keys()):
        T_over_Tc = correlation_data[t_idx]['T_over_Tc']
        G_r = np.array(correlation_data[t_idx]['G_r'])

        for r_max in r_max_values:
            pair_count += 1
            t0 = time.time()

            result = run_windowed_sweep(G_r, N, r_max, n_centers=n_centers)

            dt = time.time() - t0

            if result is not None:
                pr_str = f"{result['pr_mean']:.3f}"
                rank_str = f"{result['rank_mean']:.2f}"
                sv2_str = f"{result['sv2sv1_mean']:.3f}"
                nv = result['n_valid']
            else:
                pr_str = "FAIL"
                rank_str = "FAIL"
                sv2_str = "FAIL"
                nv = 0

            # Print progress every few pairs or at interesting points
            if pair_count % 14 == 0 or r_max == r_max_values[-1]:
                print(f"  [{pair_count}/{total_pairs}] T/Tc={T_over_Tc:.3f} r_max={r_max:3d}: "
                      f"PR={pr_str} rank={rank_str} SV2/SV1={sv2_str} "
                      f"(n={nv}) [{dt:.1f}s]")

            all_fim_results[(t_idx, r_max)] = result

    return all_fim_results


# ============================================================
# MODULE 8: Phase C — Onsager Thermodynamics
# ============================================================

def run_phase_c(temperatures):
    """
    Phase C: Compute exact Onsager thermodynamics at each T.
    Returns dict of T_idx -> thermodynamic values.
    """
    print("\n" + "=" * 70)
    print("  PHASE C: Onsager Exact Thermodynamics")
    print("=" * 70)

    onsager_data = {}

    for t_idx, T in enumerate(temperatures):
        T_over_Tc = T / T_C

        # Avoid exact T_c singularity
        T_calc = T
        if abs(T_over_Tc - 1.0) < 0.001:
            T_calc = T_C * 1.001

        u = onsager_energy(T_calc)
        f = onsager_free_energy(T_calc)
        s = onsager_entropy(T_calc)
        c = onsager_specific_heat(T_calc)

        onsager_data[t_idx] = {
            'T': float(T),
            'T_over_Tc': float(T_over_Tc),
            'energy': float(u),
            'free_energy': float(f),
            'entropy': float(s),
            'specific_heat': float(c),
        }

        print(f"  T/Tc={T_over_Tc:.3f}: S={s:.4f}  C={c:.4f}  u={u:.4f}")

    return onsager_data


# ============================================================
# MODULE 9: Plots
# ============================================================

def plot_pr_accumulation_curves(all_fim_results, correlation_data,
                                 r_max_values, temperatures):
    """
    KEY PLOT: PR(r_max) accumulation curves at 6 selected temperatures.
    """
    target_ratios = [0.80, 0.95, 0.99, 1.00, 1.05, 1.40]
    colors = ['tab:blue', 'tab:cyan', 'tab:green', 'tab:red',
              'tab:orange', 'tab:purple']

    fig, ax = plt.subplots(figsize=(12, 7))

    for t_target, color in zip(target_ratios, colors):
        # Find closest T index
        t_idx = min(range(len(temperatures)),
                    key=lambda i: abs(temperatures[i] / T_C - t_target))
        T_over_Tc = temperatures[t_idx] / T_C

        pr_vals = []
        pr_errs = []
        r_vals = []

        for r_max in r_max_values:
            key = (t_idx, r_max)
            if key in all_fim_results and all_fim_results[key] is not None:
                pr_vals.append(all_fim_results[key]['pr_mean'])
                pr_errs.append(all_fim_results[key]['pr_std'])
                r_vals.append(r_max)

        if r_vals:
            ax.errorbar(r_vals, pr_vals, yerr=pr_errs, fmt='o-',
                       color=color, label=f'T/Tc={T_over_Tc:.2f}',
                       capsize=3, markersize=4)

    ax.set_xlabel('Window radius r_max', fontsize=12)
    ax.set_ylabel('Participation Ratio PR', fontsize=12)
    ax.set_title('Information Accumulation Curves: PR(r_max) at Fixed T', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, "pr_accumulation_curves.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_sum_rule_vs_entropy(I_PR, onsager_data, temperatures):
    """
    KEY PLOT: I_PR(T) vs S(T) scatter plot.
    """
    fig, ax = plt.subplots(figsize=(10, 7))

    # Separate T > Tc and T <= Tc
    t_above = []
    s_above = []
    ipr_above = []
    t_below = []
    s_below = []
    ipr_below = []

    for t_idx in sorted(I_PR.keys()):
        T_over_Tc = temperatures[t_idx] / T_C
        s_val = onsager_data[t_idx]['entropy']
        ipr_val = I_PR[t_idx]

        if T_over_Tc > 1.0:
            t_above.append(T_over_Tc)
            s_above.append(s_val)
            ipr_above.append(ipr_val)
        else:
            t_below.append(T_over_Tc)
            s_below.append(s_val)
            ipr_below.append(ipr_val)

    ax.scatter(s_above, ipr_above, c='tab:blue', s=60, label='T > Tc',
              edgecolors='black', zorder=5)
    ax.scatter(s_below, ipr_below, c='tab:red', s=60, marker='s',
              label='T < Tc', edgecolors='black', zorder=5)

    # Annotate with T/Tc
    for s, ipr, t_ratio in zip(s_above, ipr_above, t_above):
        ax.annotate(f'{t_ratio:.2f}', (s, ipr), fontsize=7,
                   textcoords='offset points', xytext=(5, 5))
    for s, ipr, t_ratio in zip(s_below, ipr_below, t_below):
        ax.annotate(f'{t_ratio:.2f}', (s, ipr), fontsize=7,
                   textcoords='offset points', xytext=(5, 5))

    # Spearman correlation for T > Tc
    if len(s_above) >= 3:
        rho, pval = spearmanr(s_above, ipr_above)
        ax.set_title(f'Scale-Integrated PR vs Onsager Entropy\n'
                    f'T > Tc: Spearman rho = {rho:.3f} (p = {pval:.4f})',
                    fontsize=13)
    else:
        ax.set_title('Scale-Integrated PR vs Onsager Entropy', fontsize=13)

    ax.set_xlabel('Onsager Entropy S(T)', fontsize=12)
    ax.set_ylabel('I_PR(T) = integral PR(r_max) dr_max', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, "sum_rule_vs_entropy.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_sum_rule_vs_specific_heat(I_PR, onsager_data, temperatures):
    """
    KEY PLOT: I_PR(T) vs C(T) scatter plot.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    t_ratios = []
    c_vals = []
    ipr_vals = []

    for t_idx in sorted(I_PR.keys()):
        T_over_Tc = temperatures[t_idx] / T_C
        t_ratios.append(T_over_Tc)
        c_vals.append(onsager_data[t_idx]['specific_heat'])
        ipr_vals.append(I_PR[t_idx])

    # Left: I_PR(T) and C(T) vs T/Tc
    ax1.plot(t_ratios, ipr_vals, 'o-', color='tab:blue', label='I_PR(T)',
            markersize=5)
    ax1_twin = ax1.twinx()
    ax1_twin.plot(t_ratios, c_vals, 's-', color='tab:red', label='C(T)',
                 markersize=5)
    ax1.axvline(1.0, color='black', ls='--', alpha=0.3, label='T_c')
    ax1.set_xlabel('T / T_c', fontsize=12)
    ax1.set_ylabel('I_PR(T)', color='tab:blue', fontsize=12)
    ax1_twin.set_ylabel('Specific Heat C(T)', color='tab:red', fontsize=12)
    ax1.set_title('I_PR and C(T) vs Temperature', fontsize=13)
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1_twin.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=9)

    # Right: scatter I_PR vs C
    ax2.scatter(c_vals, ipr_vals, c=[abs(t - 1.0) for t in t_ratios],
               cmap='coolwarm_r', s=60, edgecolors='black', zorder=5)
    for c, ipr, t_ratio in zip(c_vals, ipr_vals, t_ratios):
        ax2.annotate(f'{t_ratio:.2f}', (c, ipr), fontsize=7,
                    textcoords='offset points', xytext=(5, 5))

    if len(c_vals) >= 3:
        rho, pval = spearmanr(c_vals, ipr_vals)
        ax2.set_title(f'I_PR vs C(T)\nSpearman rho = {rho:.3f} (p = {pval:.4f})',
                     fontsize=13)
    else:
        ax2.set_title('I_PR vs C(T)', fontsize=13)

    ax2.set_xlabel('Specific Heat C(T)', fontsize=12)
    ax2.set_ylabel('I_PR(T)', fontsize=12)
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, "sum_rule_vs_specific_heat.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_heatmaps(all_fim_results, temperatures, r_max_values):
    """
    Heatmaps: SV2/SV1, rank, PR as functions of (T/Tc, r_max).
    """
    n_T = len(temperatures)
    n_R = len(r_max_values)
    t_ratios = [T / T_C for T in temperatures]

    # Build 2D arrays
    sv2sv1_grid = np.full((n_T, n_R), np.nan)
    rank_grid = np.full((n_T, n_R), np.nan)
    pr_grid = np.full((n_T, n_R), np.nan)

    for t_idx in range(n_T):
        for r_idx, r_max in enumerate(r_max_values):
            key = (t_idx, r_max)
            if key in all_fim_results and all_fim_results[key] is not None:
                sv2sv1_grid[t_idx, r_idx] = all_fim_results[key]['sv2sv1_mean']
                rank_grid[t_idx, r_idx] = all_fim_results[key]['rank_mean']
                pr_grid[t_idx, r_idx] = all_fim_results[key]['pr_mean']

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # SV2/SV1 heatmap
    im0 = axes[0].imshow(sv2sv1_grid, aspect='auto', origin='lower',
                          cmap='RdYlBu_r', vmin=0, vmax=1)
    axes[0].set_title('SV2/SV1(T, r_max)', fontsize=12)
    axes[0].set_xlabel('r_max index')
    axes[0].set_ylabel('T/Tc index')
    axes[0].set_xticks(range(n_R))
    axes[0].set_xticklabels([str(r) for r in r_max_values], rotation=45, fontsize=7)
    axes[0].set_yticks(range(n_T))
    axes[0].set_yticklabels([f'{t:.2f}' for t in t_ratios], fontsize=7)
    plt.colorbar(im0, ax=axes[0])

    # Rank heatmap
    im1 = axes[1].imshow(rank_grid, aspect='auto', origin='lower',
                          cmap='viridis', vmin=1, vmax=4)
    axes[1].set_title('D_eff(T, r_max)', fontsize=12)
    axes[1].set_xlabel('r_max index')
    axes[1].set_xticks(range(n_R))
    axes[1].set_xticklabels([str(r) for r in r_max_values], rotation=45, fontsize=7)
    axes[1].set_yticks(range(n_T))
    axes[1].set_yticklabels([f'{t:.2f}' for t in t_ratios], fontsize=7)
    plt.colorbar(im1, ax=axes[1])

    # PR heatmap
    im2 = axes[2].imshow(pr_grid, aspect='auto', origin='lower',
                          cmap='magma')
    axes[2].set_title('PR(T, r_max)', fontsize=12)
    axes[2].set_xlabel('r_max index')
    axes[2].set_xticks(range(n_R))
    axes[2].set_xticklabels([str(r) for r in r_max_values], rotation=45, fontsize=7)
    axes[2].set_yticks(range(n_T))
    axes[2].set_yticklabels([f'{t:.2f}' for t in t_ratios], fontsize=7)
    plt.colorbar(im2, ax=axes[2])

    plt.suptitle('Windowed FIM Diagnostics: (T/Tc, r_max) Plane', fontsize=14)
    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, "sv2sv1_heatmap.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")

    # Individual heatmaps for rank and PR
    for grid, name, cmap, vlims in [
        (rank_grid, "rank_heatmap", "viridis", (1, 4)),
        (pr_grid, "pr_heatmap", "magma", (None, None)),
    ]:
        fig, ax = plt.subplots(figsize=(10, 8))
        vmin, vmax = vlims
        im = ax.imshow(grid, aspect='auto', origin='lower', cmap=cmap,
                       vmin=vmin, vmax=vmax)
        ax.set_xlabel('r_max', fontsize=12)
        ax.set_ylabel('T / Tc', fontsize=12)
        ax.set_xticks(range(n_R))
        ax.set_xticklabels([str(r) for r in r_max_values], rotation=45)
        ax.set_yticks(range(n_T))
        ax.set_yticklabels([f'{t:.2f}' for t in t_ratios])
        ax.set_title(f'{name.replace("_", " ").title()}', fontsize=13)
        plt.colorbar(im, ax=ax)
        plt.tight_layout()
        fname = os.path.join(OUTPUT_DIR, f"{name}.png")
        plt.savefig(fname, dpi=150)
        plt.close()
        print(f"  Saved {fname}")


def plot_onsager_thermodynamics(onsager_data, temperatures):
    """Onsager S(T), u(T), C(T) reference curves."""
    t_ratios = [temperatures[i] / T_C for i in sorted(onsager_data.keys())]
    S_vals = [onsager_data[i]['entropy'] for i in sorted(onsager_data.keys())]
    u_vals = [onsager_data[i]['energy'] for i in sorted(onsager_data.keys())]
    C_vals = [onsager_data[i]['specific_heat'] for i in sorted(onsager_data.keys())]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    axes[0].plot(t_ratios, S_vals, 'o-', color='tab:blue')
    axes[0].axvline(1.0, color='red', ls='--', alpha=0.5)
    axes[0].set_xlabel('T / Tc')
    axes[0].set_ylabel('Entropy S(T)')
    axes[0].set_title('Onsager Exact Entropy')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t_ratios, u_vals, 's-', color='tab:green')
    axes[1].axvline(1.0, color='red', ls='--', alpha=0.5)
    axes[1].set_xlabel('T / Tc')
    axes[1].set_ylabel('Energy u(T)')
    axes[1].set_title('Onsager Exact Energy')
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t_ratios, C_vals, '^-', color='tab:red')
    axes[2].axvline(1.0, color='red', ls='--', alpha=0.5)
    axes[2].set_xlabel('T / Tc')
    axes[2].set_ylabel('Specific Heat C(T)')
    axes[2].set_title('Onsager Exact Specific Heat')
    axes[2].grid(True, alpha=0.3)

    plt.suptitle('Onsager Exact Thermodynamics — 2D Ising', fontsize=14)
    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, "onsager_thermodynamics.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_phase2_sanity_check(all_fim_results, temperatures, r_max_values):
    """Phase 2 sanity check: r_max=64 diagnostics vs temperature."""
    r_max_full = max(r_max_values)  # Should be 64

    t_ratios = []
    ranks = []
    rank_errs = []
    sv2sv1s = []
    sv2sv1_errs = []
    etas = []
    eta_errs = []
    prs = []
    pr_errs = []

    for t_idx in range(len(temperatures)):
        key = (t_idx, r_max_full)
        if key in all_fim_results and all_fim_results[key] is not None:
            t_ratios.append(temperatures[t_idx] / T_C)
            ranks.append(all_fim_results[key]['rank_mean'])
            rank_errs.append(all_fim_results[key]['rank_std'])
            sv2sv1s.append(all_fim_results[key]['sv2sv1_mean'])
            sv2sv1_errs.append(all_fim_results[key]['sv2sv1_std'])
            etas.append(all_fim_results[key]['eta_mean'])
            eta_errs.append(all_fim_results[key]['eta_std'])
            prs.append(all_fim_results[key]['pr_mean'])
            pr_errs.append(all_fim_results[key]['pr_std'])

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Phase 2 Sanity Check: r_max = {r_max_full} (full lattice)',
                fontsize=14, fontweight='bold')

    # Rank
    axes[0, 0].errorbar(t_ratios, ranks, yerr=rank_errs, fmt='o-',
                        color='tab:blue', capsize=3)
    axes[0, 0].axvline(1.0, color='red', ls='--', alpha=0.5)
    axes[0, 0].axhline(3, color='gray', ls=':', alpha=0.5, label='rank=3')
    axes[0, 0].set_ylabel('D_eff (rank)')
    axes[0, 0].set_title('(a) Effective Dimension')
    axes[0, 0].legend()

    # SV2/SV1
    axes[0, 1].errorbar(t_ratios, sv2sv1s, yerr=sv2sv1_errs, fmt='s-',
                        color='tab:orange', capsize=3)
    axes[0, 1].axvline(1.0, color='red', ls='--', alpha=0.5)
    axes[0, 1].axhline(1.0, color='gray', ls=':', alpha=0.5, label='SV2/SV1=1')
    axes[0, 1].set_ylabel('SV2/SV1')
    axes[0, 1].set_title('(b) SV Ratio')
    axes[0, 1].legend()

    # eta
    axes[1, 0].errorbar(t_ratios, etas, yerr=eta_errs, fmt='^-',
                        color='tab:green', capsize=3)
    axes[1, 0].axvline(1.0, color='red', ls='--', alpha=0.5)
    axes[1, 0].set_xlabel('T / Tc')
    axes[1, 0].set_ylabel('Disorder Index eta')
    axes[1, 0].set_title('(c) Disorder Index')

    # PR
    axes[1, 1].errorbar(t_ratios, prs, yerr=pr_errs, fmt='D-',
                        color='tab:purple', capsize=3)
    axes[1, 1].axvline(1.0, color='red', ls='--', alpha=0.5)
    axes[1, 1].set_xlabel('T / Tc')
    axes[1, 1].set_ylabel('Participation Ratio')
    axes[1, 1].set_title('(d) PR')

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, "phase2_sanity_check.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


def plot_accumulation_onset_vs_xi(all_fim_results, correlation_data,
                                    temperatures, r_max_values):
    """r_onset vs xi(T) for T > Tc."""
    t_above_tc = []
    xi_vals = []
    r_onset_vals = []

    for t_idx in range(len(temperatures)):
        T_over_Tc = temperatures[t_idx] / T_C
        if T_over_Tc <= 1.05:
            continue

        xi = correlation_data[t_idx]['xi']
        if xi <= 0 or not np.isfinite(xi):
            continue

        # Compute r_onset: r_max where PR first reaches 90% of plateau
        pr_vals = []
        r_vals = []
        for r_max in r_max_values:
            key = (t_idx, r_max)
            if key in all_fim_results and all_fim_results[key] is not None:
                pr_vals.append(all_fim_results[key]['pr_mean'])
                r_vals.append(r_max)

        if len(pr_vals) < 3:
            continue

        plateau = pr_vals[-1]  # PR at largest r_max
        threshold = 0.9 * plateau

        r_onset = r_vals[-1]  # Default to max if never reaches
        for r, pr in zip(r_vals, pr_vals):
            if pr >= threshold:
                r_onset = r
                break

        t_above_tc.append(T_over_Tc)
        xi_vals.append(xi)
        r_onset_vals.append(r_onset)

    fig, ax = plt.subplots(figsize=(8, 6))

    if len(xi_vals) >= 2:
        ax.scatter(xi_vals, r_onset_vals, c='tab:blue', s=60,
                  edgecolors='black', zorder=5)
        for xi, r_on, t_ratio in zip(xi_vals, r_onset_vals, t_above_tc):
            ax.annotate(f'{t_ratio:.2f}', (xi, r_on), fontsize=8,
                       textcoords='offset points', xytext=(5, 5))

        # Diagonal reference
        max_val = max(max(xi_vals), max(r_onset_vals))
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='r_onset = xi')

        if len(xi_vals) >= 3:
            rho, pval = spearmanr(xi_vals, r_onset_vals)
            ax.set_title(f'Accumulation Onset vs Correlation Length\n'
                        f'Spearman rho = {rho:.3f} (p = {pval:.4f})',
                        fontsize=13)
        else:
            ax.set_title('Accumulation Onset vs Correlation Length', fontsize=13)

    ax.set_xlabel('Correlation length xi(T)', fontsize=12)
    ax.set_ylabel('r_onset (90% of plateau)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, "accumulation_onset_vs_xi.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Saved {fname}")


# ============================================================
# MODULE 10: Results Report
# ============================================================

def evaluate_predictions(all_fim_results, I_PR, I_rank, onsager_data,
                          correlation_data, temperatures, r_max_values):
    """Evaluate all pre-registered predictions. Returns verdict strings."""
    verdicts = {}

    # --- P_Q1_5: Phase 2 sanity check (GATE) ---
    r_full = max(r_max_values)
    tc_idx = min(range(len(temperatures)),
                 key=lambda i: abs(temperatures[i] / T_C - 1.0))
    key_tc = (tc_idx, r_full)
    p5_pass = True
    p5_notes = []

    if key_tc in all_fim_results and all_fim_results[key_tc] is not None:
        res = all_fim_results[key_tc]
        rank_tc = res['rank_mean']
        sv2sv1_tc = res['sv2sv1_mean']
        eta_tc = res['eta_mean']

        # Check rank ~ 3 at T_c
        if abs(rank_tc - 3.0) > 1.0:
            p5_pass = False
            p5_notes.append(f"rank={rank_tc:.2f} (expected ~3)")
        else:
            p5_notes.append(f"rank={rank_tc:.2f} OK")

        # Check SV2/SV1 ~ 1 at T_c
        if sv2sv1_tc < 0.8:
            p5_pass = False
            p5_notes.append(f"SV2/SV1={sv2sv1_tc:.3f} (expected ~1)")
        else:
            p5_notes.append(f"SV2/SV1={sv2sv1_tc:.3f} OK")

        # Check eta minimum near T_c
        p5_notes.append(f"eta={eta_tc:.3f}")
    else:
        p5_pass = False
        p5_notes.append("No data at (T_c, r_max=64)")

    verdicts['P_Q1_5'] = ('PASS' if p5_pass else 'FAIL', '; '.join(p5_notes))

    # --- P_Q1_1: Three-regime accumulation structure ---
    p1_results = {}

    # High T (T/Tc > 1.3): plateau before r_max=20
    high_t_pass = True
    high_t_notes = []
    for t_idx in range(len(temperatures)):
        T_over_Tc = temperatures[t_idx] / T_C
        if T_over_Tc < 1.3:
            continue

        pr_20 = None
        pr_64 = None
        for r_max in r_max_values:
            key = (t_idx, r_max)
            if key in all_fim_results and all_fim_results[key] is not None:
                if r_max == 20:
                    pr_20 = all_fim_results[key]['pr_mean']
                if r_max == max(r_max_values):
                    pr_64 = all_fim_results[key]['pr_mean']

        if pr_20 is not None and pr_64 is not None and pr_64 > 0:
            change = abs(pr_64 - pr_20) / pr_64
            if change > 0.05:
                high_t_pass = False
            high_t_notes.append(f"T/Tc={T_over_Tc:.2f}: delta={change:.3f}")

    p1_results['high_T'] = ('PASS' if high_t_pass else 'FAIL',
                            '; '.join(high_t_notes) if high_t_notes else 'no data')

    # T ~ Tc: continuous rise
    tc_pass = False
    tc_notes = []
    for t_idx in range(len(temperatures)):
        T_over_Tc = temperatures[t_idx] / T_C
        if not (0.97 < T_over_Tc < 1.03):
            continue

        pr_10 = None
        pr_64 = None
        for r_max in r_max_values:
            key = (t_idx, r_max)
            if key in all_fim_results and all_fim_results[key] is not None:
                if r_max == 10:
                    pr_10 = all_fim_results[key]['pr_mean']
                if r_max == max(r_max_values):
                    pr_64 = all_fim_results[key]['pr_mean']

        if pr_10 is not None and pr_64 is not None and pr_10 > 0:
            ratio = pr_64 / pr_10
            if ratio > 1.3:
                tc_pass = True
            tc_notes.append(f"T/Tc={T_over_Tc:.2f}: PR(64)/PR(10)={ratio:.3f}")

    p1_results['near_Tc'] = ('PASS' if tc_pass else 'FAIL',
                             '; '.join(tc_notes) if tc_notes else 'no data')

    # Low T (T/Tc < 0.90): high noise
    low_t_pass = False
    low_t_notes = []
    for t_idx in range(len(temperatures)):
        T_over_Tc = temperatures[t_idx] / T_C
        if T_over_Tc > 0.90:
            continue

        noisy_count = 0
        total_count = 0
        for r_max in r_max_values:
            key = (t_idx, r_max)
            if key in all_fim_results and all_fim_results[key] is not None:
                total_count += 1
                if all_fim_results[key]['pr_std'] > 0.3:
                    noisy_count += 1

        if total_count > 0:
            frac = noisy_count / total_count
            if frac > 0.5:
                low_t_pass = True
            low_t_notes.append(f"T/Tc={T_over_Tc:.2f}: {noisy_count}/{total_count} noisy")

    p1_results['low_T'] = ('PASS' if low_t_pass else 'FAIL',
                           '; '.join(low_t_notes) if low_t_notes else 'no data')

    all_p1_pass = all(v[0] == 'PASS' for v in p1_results.values())
    verdicts['P_Q1_1'] = ('PASS' if all_p1_pass else 'MIXED',
                          str(p1_results))

    # --- P_Q1_2: I_PR monotonic with S for T > Tc ---
    s_above = []
    ipr_above = []
    for t_idx in sorted(I_PR.keys()):
        T_over_Tc = temperatures[t_idx] / T_C
        if T_over_Tc > 1.0:
            s_above.append(onsager_data[t_idx]['entropy'])
            ipr_above.append(I_PR[t_idx])

    if len(s_above) >= 3:
        rho_s, pval_s = spearmanr(s_above, ipr_above)
        if rho_s > 0.85:
            p2_verdict = 'PASS'
        elif rho_s < 0.5:
            p2_verdict = 'FAIL'
        else:
            p2_verdict = 'INCONCLUSIVE'
        verdicts['P_Q1_2'] = (p2_verdict, f"rho={rho_s:.3f}, p={pval_s:.4f}, "
                              f"n={len(s_above)} points above Tc")
    else:
        verdicts['P_Q1_2'] = ('INCONCLUSIVE', 'Insufficient data above Tc')

    # --- P_Q1_3: I_PR peaks near T_c ---
    all_t_ratios = [temperatures[t_idx] / T_C for t_idx in sorted(I_PR.keys())]
    all_ipr = [I_PR[t_idx] for t_idx in sorted(I_PR.keys())]

    if len(all_ipr) >= 3:
        max_idx = int(np.argmax(all_ipr))
        t_at_max = all_t_ratios[max_idx]
        if abs(t_at_max - 1.0) < 0.05:
            p3_verdict = 'PASS'
        else:
            # Check if monotonic (no extremum)
            diffs = [all_ipr[i+1] - all_ipr[i] for i in range(len(all_ipr)-1)]
            if all(d >= 0 for d in diffs) or all(d <= 0 for d in diffs):
                p3_verdict = 'FAIL'
            else:
                p3_verdict = 'FAIL'
        verdicts['P_Q1_3'] = (p3_verdict,
                              f"I_PR max at T/Tc={t_at_max:.3f}, "
                              f"I_PR_max={all_ipr[max_idx]:.2f}")
    else:
        verdicts['P_Q1_3'] = ('INCONCLUSIVE', 'Insufficient data')

    # --- P_Q1_4: r_onset tracks xi ---
    xi_vals = []
    r_onset_vals = []
    for t_idx in range(len(temperatures)):
        T_over_Tc = temperatures[t_idx] / T_C
        if T_over_Tc <= 1.05:
            continue
        xi = correlation_data[t_idx]['xi']
        if xi <= 0 or not np.isfinite(xi):
            continue

        pr_vals = []
        r_vals = []
        for r_max in r_max_values:
            key = (t_idx, r_max)
            if key in all_fim_results and all_fim_results[key] is not None:
                pr_vals.append(all_fim_results[key]['pr_mean'])
                r_vals.append(r_max)
        if len(pr_vals) < 3:
            continue

        plateau = pr_vals[-1]
        threshold = 0.9 * plateau
        r_onset = r_vals[-1]
        for r, pr in zip(r_vals, pr_vals):
            if pr >= threshold:
                r_onset = r
                break

        xi_vals.append(xi)
        r_onset_vals.append(r_onset)

    if len(xi_vals) >= 3:
        rho_xi, pval_xi = spearmanr(xi_vals, r_onset_vals)
        if rho_xi > 0.8:
            p4_verdict = 'PASS'
        elif rho_xi < 0.4:
            p4_verdict = 'FAIL'
        else:
            p4_verdict = 'INCONCLUSIVE'
        verdicts['P_Q1_4'] = (p4_verdict, f"rho={rho_xi:.3f}, p={pval_xi:.4f}")
    else:
        verdicts['P_Q1_4'] = ('INCONCLUSIVE', f'Only {len(xi_vals)} points')

    return verdicts


def write_results_report(verdicts, I_PR, I_rank, onsager_data,
                          all_fim_results, correlation_data,
                          temperatures, r_max_values, total_time):
    """Write Q1_SUM_RULE_RESULTS.md"""
    lines = []

    # First line: headline verdict
    p2_v = verdicts.get('P_Q1_2', ('?', ''))[0]
    p3_v = verdicts.get('P_Q1_3', ('?', ''))[0]
    p2_info = verdicts.get('P_Q1_2', ('?', ''))[1]
    p3_info = verdicts.get('P_Q1_3', ('?', ''))[1]

    if p3_v == 'PASS':
        headline = (f"EXTREMUM POSITIVE: I_PR peaks at T_c — "
                   f"scale-integrated information ~ entropy production capacity")
    elif p2_v == 'PASS':
        headline = (f"SUM RULE POSITIVE: I_PR correlates with S — "
                   f"{p2_info}")
    elif p2_v == 'FAIL' and p3_v == 'FAIL':
        headline = (f"SUM RULE NEGATIVE: I_PR uncorrelated with S ({p2_info}) "
                   f"and no extremum at T_c ({p3_info}) — programme reassesses")
    else:
        headline = f"INCONCLUSIVE: Mixed results — see detailed assessment"

    lines.append(f'"{headline}"\n')
    lines.append(f"# Q1: Scale-Dependent Information Sum Rule — Results\n")
    lines.append(f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Total runtime: {total_time:.1f}s ({total_time/60:.1f} min)\n")

    # Prediction verdicts
    lines.append("## Pre-Registered Predictions\n")
    lines.append("| ID | Prediction | Verdict | Details |")
    lines.append("|---|---|---|---|")
    for pid in ['P_Q1_5', 'P_Q1_1', 'P_Q1_2', 'P_Q1_3', 'P_Q1_4']:
        v, detail = verdicts.get(pid, ('N/A', 'Not evaluated'))
        desc = {
            'P_Q1_5': 'Phase 2 sanity check (GATE)',
            'P_Q1_1': 'Three-regime accumulation structure',
            'P_Q1_2': 'I_PR monotonic with S(T) for T > Tc',
            'P_Q1_3': 'I_PR peaks at T_c',
            'P_Q1_4': 'r_onset tracks xi(T)',
        }
        lines.append(f"| {pid} | {desc.get(pid, pid)} | **{v}** | {detail} |")
    lines.append("")

    # Kill condition assessment
    lines.append("## Kill Condition Assessment\n")
    if p2_v == 'FAIL' and p3_v == 'FAIL':
        lines.append("**KILL TRIGGERED**: Both P_Q1_2 and P_Q1_3 failed. "
                     "D_eff(r_max) has no thermodynamic content discoverable "
                     "by this construction. Q2 is not motivated.\n")
    elif verdicts.get('P_Q1_5', ('?',))[0] == 'FAIL':
        lines.append("**GATE FAILED**: P_Q1_5 (Phase 2 sanity) failed. "
                     "Fix windowed construction before evaluating other predictions.\n")
    else:
        lines.append("Kill condition NOT triggered.\n")

    # Sum rule values
    lines.append("## Scale-Integrated PR Values\n")
    lines.append("| T/Tc | I_PR | I_rank | S(T) | C(T) |")
    lines.append("|------|------|--------|------|------|")
    for t_idx in sorted(I_PR.keys()):
        T_over_Tc = temperatures[t_idx] / T_C
        ipr = I_PR[t_idx]
        irank = I_rank.get(t_idx, 0)
        s = onsager_data[t_idx]['entropy']
        c = onsager_data[t_idx]['specific_heat']
        lines.append(f"| {T_over_Tc:.3f} | {ipr:.2f} | {irank:.2f} | "
                    f"{s:.4f} | {c:.4f} |")
    lines.append("")

    # Accumulation curve shapes
    lines.append("## PR Accumulation Curve Shapes\n")
    target_ratios = [0.80, 0.95, 0.99, 1.00, 1.05, 1.40]
    for t_target in target_ratios:
        t_idx = min(range(len(temperatures)),
                    key=lambda i: abs(temperatures[i] / T_C - t_target))
        T_over_Tc = temperatures[t_idx] / T_C
        lines.append(f"### T/Tc = {T_over_Tc:.3f}\n")
        lines.append("| r_max | PR (mean +/- std) | Rank | SV2/SV1 |")
        lines.append("|-------|-------------------|------|---------|")
        for r_max in r_max_values:
            key = (t_idx, r_max)
            if key in all_fim_results and all_fim_results[key] is not None:
                r = all_fim_results[key]
                lines.append(f"| {r_max} | {r['pr_mean']:.3f} +/- {r['pr_std']:.3f} | "
                           f"{r['rank_mean']:.2f} | {r['sv2sv1_mean']:.3f} |")
            else:
                lines.append(f"| {r_max} | FAILED | - | - |")
        lines.append("")

    # Next steps
    lines.append("## Next Step Recommendation\n")
    if p3_v == 'PASS':
        lines.append("**Proceed to Q2**: I_PR peaks at T_c, suggesting correlation "
                    "with specific heat C(T). Q2 should test the specific functional "
                    "form of this relationship using existing Phase 2 data.\n")
    elif p2_v == 'PASS':
        lines.append("**Proceed to Q2**: I_PR correlates monotonically with S(T) "
                    "above T_c. Q2 should examine whether this is a linear, power-law, "
                    "or other functional relationship.\n")
    elif p2_v == 'FAIL' and p3_v == 'FAIL':
        lines.append("**Programme reassesses**: The windowed kernel approach does not "
                    "reveal thermodynamic content. Consider:\n"
                    "- Q1b: Coarse-graining (block-spin) construction\n"
                    "- Alternative observables beyond PR\n"
                    "- Accept that D_eff is geometric, not thermodynamic\n")
    else:
        lines.append("**Collect more data**: Inconclusive results suggest more "
                    "temperatures or larger lattice may resolve the picture.\n")

    lines.append(f"\n---\n*End of Q1 Results. Runtime: {total_time:.1f}s*\n")

    report_path = os.path.join(OUTPUT_DIR, "Q1_SUM_RULE_RESULTS.md")
    with open(report_path, 'w') as f:
        f.write('\n'.join(lines))
    print(f"  Saved {report_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("  Q1: Scale-Dependent Information Sum Rule")
    print("  D_eff(r_max, T) on 2D Ising with Windowed Thermal Kernel")
    print("=" * 70)
    print(f"T_c = {T_C:.6f}")
    print(f"Lattice: N = {N_LATTICE}")
    print(f"Temperatures: {len(TEMPERATURES)} points, T/Tc in "
          f"[{min(T_OVER_TC):.2f}, {max(T_OVER_TC):.2f}]")
    print(f"Window radii: {R_MAX_VALUES}")
    print(f"FIM centres per (T, r_max): {N_FIM_CENTERS}")
    print(f"Output: {OUTPUT_DIR}/")
    print()

    t_total_start = time.time()

    # Phase A: Monte Carlo (cache-aware)
    corr_path = os.path.join(RAW_DIR, "correlation_functions.json")
    if os.path.exists(corr_path):
        print("\n  [CACHE] Loading correlation data from previous run...")
        with open(corr_path, 'r') as f:
            correlation_data_raw = json.load(f)
        # Re-key with int indices
        correlation_data = {int(k): v for k, v in correlation_data_raw.items()}
        print(f"  Loaded {len(correlation_data)} temperature points from {corr_path}")
    else:
        correlation_data = run_phase_a(
            TEMPERATURES, N=N_LATTICE, n_eq=N_EQ_SWEEPS, n_configs=N_CONFIGS
        )
        with open(corr_path, 'w') as f:
            json.dump(correlation_data, f, indent=2, default=str)
        print(f"  Saved {corr_path}")

    # Phase B: Windowed FIM sweep (cache-aware)
    fim_path = os.path.join(RAW_DIR, "windowed_fim_results.json")
    if os.path.exists(fim_path):
        print("\n  [CACHE] Loading FIM data from previous run...")
        with open(fim_path, 'r') as f:
            fim_raw = json.load(f)
        all_fim_results = {}
        for key_str, result in fim_raw.items():
            t_idx, r_max = key_str.split('_')
            all_fim_results[(int(t_idx), int(r_max))] = result
        print(f"  Loaded {len(all_fim_results)} (T, r_max) pairs from {fim_path}")
    else:
        all_fim_results = run_phase_b(
            correlation_data, R_MAX_VALUES, N=N_LATTICE, n_centers=N_FIM_CENTERS
        )
        fim_save = {}
        for (t_idx, r_max), result in all_fim_results.items():
            key_str = f"{t_idx}_{r_max}"
            fim_save[key_str] = result
        with open(fim_path, 'w') as f:
            json.dump(fim_save, f, indent=2, default=str)
        print(f"  Saved {fim_path}")

    # Phase C: Onsager thermodynamics
    onsager_data = run_phase_c(TEMPERATURES)

    onsager_path = os.path.join(RAW_DIR, "onsager_values.json")
    with open(onsager_path, 'w') as f:
        json.dump(onsager_data, f, indent=2, default=str)
    print(f"  Saved {onsager_path}")

    # Phase D: Sum rule computation
    print("\n" + "=" * 70)
    print("  PHASE D: Sum Rule Computation")
    print("=" * 70)

    I_PR, I_rank = compute_sum_rule(all_fim_results, R_MAX_VALUES)

    sum_rule_path = os.path.join(RAW_DIR, "sum_rule_values.json")
    with open(sum_rule_path, 'w') as f:
        json.dump({'I_PR': I_PR, 'I_rank': I_rank}, f, indent=2, default=str)
    print(f"  Saved {sum_rule_path}")

    for t_idx in sorted(I_PR.keys()):
        T_over_Tc = TEMPERATURES[t_idx] / T_C
        s = onsager_data[t_idx]['entropy']
        c = onsager_data[t_idx]['specific_heat']
        print(f"  T/Tc={T_over_Tc:.3f}: I_PR={I_PR[t_idx]:.2f}  "
              f"I_rank={I_rank[t_idx]:.2f}  S={s:.4f}  C={c:.4f}")

    # Evaluate predictions
    verdicts = evaluate_predictions(
        all_fim_results, I_PR, I_rank, onsager_data,
        correlation_data, TEMPERATURES, R_MAX_VALUES
    )

    print("\n" + "=" * 70)
    print("  PREDICTION VERDICTS")
    print("=" * 70)
    for pid in ['P_Q1_5', 'P_Q1_1', 'P_Q1_2', 'P_Q1_3', 'P_Q1_4']:
        v, detail = verdicts.get(pid, ('N/A', ''))
        print(f"  {pid}: {v}")
        print(f"    {detail}")

    # Plots
    print("\n" + "=" * 70)
    print("  GENERATING PLOTS")
    print("=" * 70)

    plot_pr_accumulation_curves(all_fim_results, correlation_data,
                                R_MAX_VALUES, TEMPERATURES)
    plot_sum_rule_vs_entropy(I_PR, onsager_data, TEMPERATURES)
    plot_sum_rule_vs_specific_heat(I_PR, onsager_data, TEMPERATURES)
    plot_heatmaps(all_fim_results, TEMPERATURES, R_MAX_VALUES)
    plot_onsager_thermodynamics(onsager_data, TEMPERATURES)
    plot_phase2_sanity_check(all_fim_results, TEMPERATURES, R_MAX_VALUES)
    plot_accumulation_onset_vs_xi(all_fim_results, correlation_data,
                                   TEMPERATURES, R_MAX_VALUES)

    total_time = time.time() - t_total_start

    # Write results report
    write_results_report(verdicts, I_PR, I_rank, onsager_data,
                         all_fim_results, correlation_data,
                         TEMPERATURES, R_MAX_VALUES, total_time)

    print(f"\n{'=' * 70}")
    print(f"  TOTAL TIME: {total_time:.1f}s ({total_time/60:.1f} min)")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
